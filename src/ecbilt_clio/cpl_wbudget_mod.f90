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
!!     and two fluxes of the day are taken:
!!         Fa   water given by the atmosphere and the land to the ocean coupler, on the ECBilt grid:
!!              -sumofwf*fracto*darea*86400*rowat (sumofwf = E - P - runoff over the ocean part, daily mean)
!!         Fo   water received by CLIO on its grid: fwat*aire (fwat, kg/m2/day, from Fa by at2oc)
!!     One line per model month (30 days) goes to outputdata/globals/cpl_wbudget<fini>.txt: the reservoirs at the end
!!     of the month, their changes over the month, Fa and Fo summed, and the three residuals
!!         Ral  = dA + dL + Fa          what the atmosphere and the land lose that does not reach the coupler
!!                                      (atmospheric clipping, see wbudget<fini>.txt, and land bookkeeping)
!!         Rcp  = Fa - Fo               lost or gained in the regridding atmosphere -> ocean
!!         Roi  = dO + dI + dS - Fo     ocean + sea ice not explained by what CLIO received (levitating sea ice,
!!                                      volume correction, snow on ice, ...)
!!     and dT = dA + dL + dO + dI + dS = Ral + Roi - Rcp, the drift of the total. All in global means, mm of water
!!     (kg/m2 over the area of the ECBilt globe). Sources outside these reservoirs (ice sheets, icebergs, CONSEAU,
!!     uncorrected freshwater fluxes) are not counted: their flags are written in the header.
!!
!!     Contract (public entry points), called from emic:
!!       cpl_wb_start()   : once, before the time loop (reservoirs at the start).
!!       cpl_wb_day()     : once a day, after clio (skipped when the ocean is prescribed).
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
        integer(ip), parameter :: NF   = 2                     !< fluxes
        integer(ip), parameter :: F_A  = 1, F_O = 2
        real(dblp),  parameter :: SEC_DAY = 86400.0_dblp

        logical,     save :: ready     = .false.
        real(dblp),  save :: globe     = 0.0_dblp              !< area of the ECBilt globe (m2)
        real(dblp),  save :: res_last(NR)                      !< reservoirs at the end of the last day (kg)
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

          globe = real(nlon, dblp)*sum(darea)
          call reservoirs(res_last)
          d_month = 0.0_dblp
          f_month = 0.0_dblp
          nday    = 0
          ready   = .true.

        end subroutine cpl_wb_start

        subroutine cpl_wb_day()

          use comemic_mod, only: imonth, iyear

          real(dblp) :: r(NR), f(NF)

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
          r(R_O) = r(R_O)*rho0
          r(R_I) = r(R_I)*rhog
          r(R_S) = r(R_S)*rhon

        end subroutine reservoirs

        subroutine fluxes(f)

          use comatm,      only: nlat, nlon, darea, iwater
          use comphys,     only: rowat
          use comemic_mod, only: fracto
          use comcoup_mod, only: sumofwf
          use bloc0_mod,   only: tms, is1, is2, js1, js2, ks2
          use bloc_mod,    only: aire
          use ice_mod,     only: fwat

          real(dblp), intent(out) :: f(NF)

          integer(ip) :: i, j

          f = 0.0_dblp
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
            enddo
          enddo

        end subroutine fluxes

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Output
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine write_month()

          use comemic_mod, only: fini, irunlabel, iobclin

          real(dblp) :: r(NR), d(NR), f(NF), ral, rcp, roi, dt

          if (.not. opened) then
            open(newunit=unit, file='outputdata/globals/cpl_wbudget'//fini//'.txt', form='formatted',                 &
                 status='replace', action='write')
            write(unit, '(a,i0)') '# Water inventory of the coupled system (cpl_wbudget_mod), global means in mm, '//  &
                                  'irunlabel = ', irunlabel
            write(unit, '(a)') '# A atmosphere, L land, O ocean (eta), I sea ice, S snow on sea ice: end of month;'//   &
                               ' d*: change over the month'
            write(unit, '(a)') '# Fa = -sumofwf (ECBilt grid), Fo = fwat (CLIO grid); Ral = dA + dL + Fa,'//           &
                               ' Rcp = Fa - Fo, Roi = dO + dI + dS - Fo, dT = Ral + Roi - Rcp'
            write(unit, '(a,7(a,i0))') '# not counted, flags:', ' ICEBERG=', ICEBERG, ' ISM=', ISM, ' CONSEAU=',     &
                                       CONSEAU, ' UNCORFLUX=', UNCORFLUX, ' APPLY_UNCORFWF=', APPLY_UNCORFWF,           &
                                       ' iobclin=', iobclin, ' ISOOCN=', ISOOCN
            write(unit, '(a)') '# year month nday     A             L             O             I             S'//     &
                               '             dA            dL            dO            dI            dS'//            &
                               '            Fa            Fo            Ral           Rcp           Roi'//            &
                               '           dT'
            opened = .true.
          endif

          r   = res_last/globe
          d   = d_month/globe
          f   = f_month/globe
          ral = d(R_A) + d(R_L) + f(F_A)
          rcp = f(F_A) - f(F_O)
          roi = d(R_O) + d(R_I) + d(R_S) - f(F_O)
          dt  = sum(d)
          write(unit, '(i6,i4,i6,18es14.6)') cur_year, cur_month, nday, r, d, f, ral, rcp, roi, dt
          flush(unit)

        end subroutine write_month

      end module cpl_wbudget_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
