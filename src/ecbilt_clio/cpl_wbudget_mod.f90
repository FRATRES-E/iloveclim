!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2026 Didier M. Roche (a.k.a. dmr) | iLOVECLIM / FRATRES coding group

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!   Style sheet: v1.0.0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [cpl_wbudget_mod]
!
!>     @brief Water inventory of the coupled system ECBilt + LBM + CLIO (ROADMAP 3b, wiso:3b:0000).
!
!      DESCRIPTION:
!
!>     Pure diagnostics: nothing here changes the model state (results bit-identical with and without them). Water
!!     only, independent of the isotopes (WISO = 0 or 1). Once a day, after the CLIO step, the reservoirs
!!         A    atmosphere, column water rmoisg                           (ECBilt grid, darea)
!!         L    land, soil water bmoisg + snow dsnow, on the land fraction (ECBilt grid, darea*fractl)
!!         O    ocean, free surface eta (CLIO grid, aire*tms, rho0)
!!         I    sea ice, (1-albq)*hgbq (rhog)
!!         S    snow on sea ice, (1-albq)*hnbq (rhon)
!!     and three fluxes of the day are taken:
!!         Fa   water given by the atmosphere and the land to the ocean coupler, on the ECBilt grid:
!!              -sumofwf*fracto*darea*86400*rowat (sumofwf = E - P - runoff over the ocean part, daily mean)
!!         Fo   water received by CLIO on its grid: fwat*aire (fwat, kg/m2/day, from Fa by at2oc)
!!         Fw   water put into the ocean volume by thersf: -phiss(:,:,0)*aire*rho0 (phiss(:,:,0), m per step, is what
!!              conti3d turns into the surface w that drives eta); Fo - Fw is what thersf keeps for the snow on sea ice
!!     (CLIO's aire is in 10^12 m2, converted here, wiso:3b:0001)
!!     CLIO keeps its volume (wiso:3b:0002): what eta does not take of the surface volume flux leaves each column through
!!     its bottom (flucor: the bottom w goes into phifs as a flux of bottom water, then w = 0 there), the global drift
!!     of each scalar is redistributed (deriv), and with vcor > 0 last year's global mean salt flux is removed from the
!!     surface (thersf, zfluxms). The freshwater the ocean receives therefore shows in its salinity, not its volume:
!!         Fb   water leaving through the bottoms of the columns, -w_bottom*dts*aire*rho0, w_bottom = phifs(2)/S_bottom
!!         Fs   freshwater implied by the change of the mean salinity, -dS/S*V*rho0 (V the ocean volume)
!!         Fc   freshwater equivalent of the vcor salt correction, zfluxms*vcor*dz(ks2)*area/S*rho0 (salt added = water
!!              removed)
!!     with the residual Rb = dO - Fw + Fb (eta + bottom against the surface volume flux).
!!     One line per model month (30 days) goes to outputdata/globals/cpl_wbudget<fini>.txt: the reservoirs at the end
!!     of the month, their changes over the month, Fa and Fo summed, and the three residuals
!!         Ral  = dA + dL + Fa          what the atmosphere and the land lose that does not reach the coupler
!!                                      (atmospheric clipping, see wbudget<fini>.txt, and land bookkeeping)
!!         Rcp  = Fa - Fo               lost or gained in the regridding atmosphere -> ocean
!!         Roi  = dO + dI + dS - Fo     ocean + sea ice not explained by what CLIO received (levitating sea ice,
!!                                      volume correction, snow on ice, ...), split into
!!         Rw   = dO - Fw               the free surface against the volume flux it was given (barotropic scheme)
!!         Rs   = Fw + dS - Fo          snow on sea ice against what thersf kept for it
!!                                      (so that Roi = Rw + Rs + dI)
!!     and dT = dA + dL + dO + dI + dS = Ral + Roi - Rcp, the drift of the total. All in global means, mm of water
!!     (kg/m2 over the area of the ECBilt globe). Sources outside these reservoirs (ice sheets, icebergs, CONSEAU,
!!     uncorrected freshwater fluxes) are not counted: their flags are written in the header.
!!
!!     Contract (public entry points), called from emic:
!!       cpl_wb_start()   : once, before the time loop (reservoirs at the start).
!!       cpl_wb_day()     : once a day, after clio (skipped when the ocean is prescribed).
!!     Active only with BUDGET_DIAG = 1 (choixcomposantes.h, wiso:3b:0009): otherwise both are empty.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module cpl_wbudget_mod

        use global_constants_mod, only: ip, dblp=>dp, days_month360d_i

        implicit none

        private

        public :: cpl_wb_start, cpl_wb_day


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Constants and module state
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip), parameter :: NR   = 5                     !< reservoirs
        integer(ip), parameter :: R_A  = 1, R_L = 2, R_O = 3, R_I = 4, R_S = 5
        integer(ip), parameter :: NF   = 6                     !< fluxes
        integer(ip), parameter :: F_A  = 1, F_O = 2, F_W = 3, F_B = 4, F_S = 5, F_C = 6
        real(dblp),  parameter :: SEC_DAY = 86400.0_dblp
        real(dblp),  parameter :: AIRE_M2 = 1.0e12_dblp        !< CLIO's aire is in 10^12 m2 (defgrid_kernel_mod)

        logical,     save :: ready     = .false.
        real(dblp),  save :: globe     = 0.0_dblp              !< area of the ECBilt globe (m2)
        real(dblp),  save :: res_last(NR)                      !< reservoirs at the end of the last day (kg)
        real(dblp),  save :: sal_last  = 0.0_dblp              !< mean ocean salinity at the end of the last day (psu)
        real(dblp),  save :: d_month(NR)                       !< changes summed over the month (kg)
        real(dblp),  save :: f_month(NF)                       !< fluxes summed over the month (kg)
        integer(ip), save :: nday      = 0
        integer(ip), save :: cur_month = -1, cur_year = -1
        integer(ip), save :: unit      = -1
        logical,     save :: opened    = .false.


      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Entry points
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine cpl_wb_start()

          use comatm, only: nlon, darea

          real(dblp) :: vol

#if ( BUDGET_DIAG == 1 )
          globe = real(nlon, dblp)*sum(darea)
          call reservoirs(res_last)
          call salinity(sal_last, vol)
          d_month = 0.0_dblp
          f_month = 0.0_dblp
          nday    = 0
          ready   = .true.

#endif
        end subroutine cpl_wb_start

        subroutine cpl_wb_day()

          use comemic_mod, only: imonth, iyear

          real(dblp) :: r(NR), f(NF)

#if ( BUDGET_DIAG == 1 )
          if (.not. ready) return

          ! the calendar month changed with days pending (incomplete month, e.g. a restart in the middle of a month)
          if (imonth /= cur_month .or. iyear /= cur_year) then
            if (nday > 0) call write_month()
            cur_month = imonth
            cur_year  = iyear
            d_month   = 0.0_dblp
            f_month   = 0.0_dblp
            nday      = 0
          endif

          call reservoirs(r)
          call fluxes(f)
          d_month  = d_month + (r - res_last)
          f_month  = f_month + f
          res_last = r
          nday     = nday + 1

          if (nday == days_month360d_i) then
            call write_month()
            d_month = 0.0_dblp
            f_month = 0.0_dblp
            nday    = 0
          endif

#endif
        end subroutine cpl_wb_day

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Reservoirs and fluxes (kg)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine reservoirs(r)

          use comatm,      only: nlat, nlon, darea, iwater
          use comphys,     only: rmoisg, rowat
          use comland_mod, only: bmoisg, dsnow, fractl, epsl
          use bloc0_mod,   only: eta, tms, is1, is2, js1, js2, ks2
          use bloc_mod,    only: aire
          use const_mod,   only: rho0
          use ice_mod,     only: hgbq, hnbq, albq, rhog, rhon

          real(dblp), intent(out) :: r(NR)

          integer(ip) :: i, j

          r = 0.0_dblp
          do j = 1, nlon
            do i = 1, nlat
              r(R_A) = r(R_A) + rmoisg(i,j,iwater)*darea(i)
              if (fractl(i,j) > epsl) then
                r(R_L) = r(R_L) + (bmoisg(i,j,iwater) + dsnow(i,j,iwater))*fractl(i,j)*darea(i)
              endif
            enddo
          enddo
          r(R_A) = r(R_A)*rowat
          r(R_L) = r(R_L)*rowat

          do j = js1, js2
            do i = is1(j), is2(j)
              r(R_O) = r(R_O) + eta(i,j)*aire(i,j)*tms(i,j,ks2)
              r(R_I) = r(R_I) + (1.0_dblp - albq(i,j))*hgbq(i,j)*aire(i,j)*tms(i,j,ks2)
              r(R_S) = r(R_S) + (1.0_dblp - albq(i,j))*hnbq(i,j)*aire(i,j)*tms(i,j,ks2)
            enddo
          enddo
          r(R_O) = r(R_O)*rho0*AIRE_M2
          r(R_I) = r(R_I)*rhog*AIRE_M2
          r(R_S) = r(R_S)*rhon*AIRE_M2

        end subroutine reservoirs

        subroutine fluxes(f)

          use comatm,      only: nlat, nlon, darea, iwater
          use comphys,     only: rowat
          use comemic_mod, only: fracto
          use comcoup_mod, only: sumofwf
          use bloc0_mod,   only: tms, is1, is2, js1, js2, ks2, phiss, phifs, scal, kfs, dts, dz, zfluxms, vcor
          use bloc_mod,    only: aire
          use const_mod,   only: rho0
          use ice_mod,     only: fwat

          real(dblp), intent(out) :: f(NF)

          integer(ip) :: i, j
          real(dblp)  :: sal, vol, surf, sb

          f    = 0.0_dblp
          surf = 0.0_dblp
          do j = 1, nlon
            do i = 1, nlat
              f(F_A) = f(F_A) - sumofwf(i,j,iwater)*fracto(i,j)*darea(i)
            enddo
          enddo
          f(F_A) = f(F_A)*SEC_DAY*rowat

          do j = js1, js2
            do i = is1(j), is2(j)
#if ( ISOOCN >= 1 )
              f(F_O) = f(F_O) + fwat(i,j,1)*aire(i,j)*tms(i,j,ks2)
#else
              f(F_O) = f(F_O) + fwat(i,j)*aire(i,j)*tms(i,j,ks2)
#endif
              f(F_W) = f(F_W) - phiss(i,j,0)*aire(i,j)*tms(i,j,ks2)
              surf   = surf + aire(i,j)*tms(i,j,ks2)
              sb = scal(i,j,kfs(i,j),2)
              if (tms(i,j,ks2) > 0.0_dblp .and. sb > 0.0_dblp) f(F_B) = f(F_B) - phifs(i,j,2)/sb*aire(i,j)
            enddo
          enddo
          f(F_O) = f(F_O)*AIRE_M2
          f(F_W) = f(F_W)*rho0*AIRE_M2
          f(F_B) = f(F_B)*dts(ks2)*rho0*AIRE_M2
          surf   = surf*AIRE_M2

          call salinity(sal, vol)
          f(F_S) = -(sal - sal_last)/sal_last*vol*rho0
#if ( ISOOCN >= 1 )
          f(F_C) = zfluxms(2)*vcor*dz(ks2)*surf/sal*rho0
#else
          f(F_C) = zfluxms*vcor*dz(ks2)*surf/sal*rho0
#endif
          sal_last = sal

        end subroutine fluxes

        subroutine salinity(sal, vol)                ! volume-mean salinity (psu) and ocean volume (m3)

          use para0_mod,   only: kmax
          use bloc0_mod,   only: tms, is1, is2, js1, js2, scal, dz
          use bloc_mod,    only: aire

          real(dblp), intent(out) :: sal, vol

          integer(ip) :: i, j, k
          real(dblp)  :: v

          sal = 0.0_dblp
          vol = 0.0_dblp
          do k = 1, kmax
            do j = js1, js2
              do i = is1(j), is2(j)
                v   = aire(i,j)*dz(k)*tms(i,j,k)
                vol = vol + v
                sal = sal + v*scal(i,j,k,2)
              enddo
            enddo
          enddo
          sal = sal/vol
          vol = vol*AIRE_M2

        end subroutine salinity

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Output
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine write_month()

          use comemic_mod, only: fini, irunlabel, iobclin

          real(dblp) :: r(NR), d(NR), f(NF), ral, rcp, roi, rw, rs, rb, dt

          if (.not. opened) then
            open(newunit=unit, file='outputdata/globals/cpl_wbudget'//fini//'.txt', form='formatted',                 &
                 status='replace', action='write')
            write(unit, '(a,i0)') '# Water inventory of the coupled system (cpl_wbudget_mod), global means in mm, '//  &
                                  'irunlabel = ', irunlabel
            write(unit, '(a)') '# A atmosphere, L land, O ocean (eta), I sea ice, S snow on sea ice: end of month;'//   &
                               ' d*: change over the month'
            write(unit, '(a)') '# Fa = -sumofwf (ECBilt grid), Fo = fwat (CLIO grid), Fw = -phiss(0) (volume flux to eta);'// &
                               ' Ral = dA + dL + Fa, Rcp = Fa - Fo, Roi = dO + dI + dS - Fo = Rw + Rs + dI,'//         &
                               ' Rw = dO - Fw, Rs = Fw + dS - Fo, dT = Ral + Roi - Rcp'
            write(unit, '(a)') '# Fb = water out through the column bottoms, Fs = freshwater from the mean salinity,'//   &
                               ' Fc = freshwater equivalent of the vcor salt correction, Rb = dO - Fw + Fb'
            write(unit, '(a,7(a,i0))') '# not counted, flags:', ' ICEBERG=', ICEBERG, ' ISM=', ISM, ' CONSEAU=',     &
                                       CONSEAU, ' UNCORFLUX=', UNCORFLUX, ' APPLY_UNCORFWF=', APPLY_UNCORFWF,           &
                                       ' iobclin=', iobclin, ' ISOOCN=', ISOOCN
            write(unit, '(a)') '# year month nday     A             L             O             I             S'//     &
                               '             dA            dL            dO            dI            dS'//            &
                               '            Fa            Fo            Fw            Fb            Fs'//            &
                               '            Fc            Ral           Rcp           Roi           Rw'//            &
                               '            Rs            Rb            dT'
            opened = .true.
          endif

          r   = res_last/globe
          d   = d_month/globe
          f   = f_month/globe
          ral = d(R_A) + d(R_L) + f(F_A)
          rcp = f(F_A) - f(F_O)
          roi = d(R_O) + d(R_I) + d(R_S) - f(F_O)
          rw  = d(R_O) - f(F_W)
          rs  = f(F_W) + d(R_S) - f(F_O)
          rb  = d(R_O) - f(F_W) + f(F_B)
          dt  = sum(d)
          write(unit, '(i6,i4,i6,25es14.6)') cur_year, cur_month, nday, r, d, f, ral, rcp, roi, rw, rs, rb, dt
          flush(unit)

        end subroutine write_month

      end module cpl_wbudget_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
