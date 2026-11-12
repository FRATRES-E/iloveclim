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
!      MODULE: [wiso_budget_mod]
!
!>     @brief Isotopologue budget and delta-bound diagnostics of the ECBilt atmosphere (ROADMAP step 2c, wiso:2c:0006).
!
!      DESCRIPTION:
!
!>     Pure diagnostics: nothing here changes the model state. One pair of files per species s = w16, w17, w18, w2h
!!     (iwat16 .. iwat2h), one line per model month, same month logic as the water budget (atm_wbudget_mod):
!!
!!     outputdata/globals/wisobudget_<s><fini>.txt, global means in mm (rsmow = 1: mm of water at standard ratio):
!!         X    = sum(a q R)/sum(a), the isotopologue inventory (W for the water)
!!         E, Pls, Pcv   evaporation, precipitation (dyrain + dysnow, corain + cosnow), as in the water budget
!!         Rtr  change over the transport (wiso_vap_transport): last-resort clip of the ratio limiter, dry cells
!!         Rmo  change over ec_moisture not explained by its E and precipitation (clip, ec_moisbalance at the ratio)
!!         Rcv  change over ec_convec_two not explained by its precipitation (clip at 0)
!!         Rot  change between ec_convec_two and the next step (nothing expected)
!!         R    = dX - (E - Pls - Pcv) = Rtr + Rmo + Rcv + Rot
!!         eps  = R/X - Rw/W, Rw the same residual of the water at the same points (non-conservation beyond the one
!!                of the water, in relative terms: the drift of the global mean delta it causes)
!!     the water columns W and Rw are repeated for reference.
!!
!!     outputdata/globals/wisobounds_<s><fini>.txt, deltas in permil, over the month, cells with water above
!!     1e-10 of the field maximum:
!!         min / max delta of the vapour, of the rain and of the snow (torain, tosnow)
!!         w2h: d-excess (dD - 8 d18O) of vapour and precipitation, min / max; w17: 17O-excess (per meg), min / max
!!         nbad   cell-steps with a non-finite ratio, a ratio <= 0, or a delta outside [delta_check_min, max]
!!                (wisoatm.nml; w16: |delta| > 1e-6, i.e. iwat16 departing from the water)
!!         no fractionation only (iso_passive: alpha = 1 everywhere, all ratios are mixtures of the sources): strict check against
!!         the source range [smin, smax] = initial vapour (at the first step of the run, restart included) and every
!!         evaporation ratio met since: nout cell-steps beyond it by more than 1e-9 (relative), worst excess, smin,
!!         smax.
!!
!!     Contract (public entry points):
!!       wiso_bud_step_begin()  : ecbilt0, with atm_wb_step_begin (before ec_moisture).
!!       wiso_bud_mark(stage)   : wiso_vap_mod, after the transport (WB_AT_TR), after the sources of ec_moisture
!!                                (WB_AT_MO, precipitation of ec_moisture known) and after convection (WB_AT_CV).
!!       wiso_bud_step_end()    : ecbilt0, with atm_wb_step_end.
!
!      REFERENCES: ROADMAP 2c, D7, D40 (deltas as fractions inside the code); atm_wbudget_mod (step 1e).
!
!      AUTHORS: dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module wiso_budget_mod

#if ( WISOATM == 1 )

        use global_constants_mod, only: ip, dblp=>dp, days_month360d_i
        use comatm,               only: nlat, nlon, iwater, iwat16, iwat17, iwat18, iwat2h

        implicit none

        private

        public :: WB_AT_TR, WB_AT_MO, WB_AT_CV
        public :: wiso_bud_step_begin, wiso_bud_mark, wiso_bud_step_end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Constants and state
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip), parameter :: WB_AT_TR = 1, WB_AT_MO = 2, WB_AT_CV = 3     !< marks inside the step
        integer(ip), parameter :: NS = iwat2h - iwat16 + 1                     !< species
        integer(ip), parameter :: NB = NS + 1                                  !< budgets: water (1) and species (2:)
        character(len=3), parameter :: SNAME(NS) = ['w16', 'w17', 'w18', 'w2h']

        ! budget terms
        integer(ip), parameter :: T_E = 1, T_PLS = 2, T_PCV = 3, T_RTR = 4, T_RMO = 5, T_RCV = 6, T_ROT = 7,      &
                                  T_R = 8, T_DX = 9, NT = 9
        ! bound statistics: min/max of vapour, rain, snow deltas, excess of vapour and precipitation
        integer(ip), parameter :: B_VMIN = 1, B_VMAX = 2, B_RMIN = 3, B_RMAX = 4, B_SMIN = 5, B_SMAX = 6,         &
                                  B_XVMIN = 7, B_XVMAX = 8, B_XPMIN = 9, B_XPMAX = 10, B_WORST = 11, NBS = 11
        real(dblp),  parameter :: WB_MM    = 1000.0_dblp          !< m -> mm
        real(dblp),  parameter :: WB_TINY  = 1.0e-10_dblp         !< water below this fraction of its max: no ratio
        real(dblp),  parameter :: WB_TOL   = 1.0e-9_dblp          !< strict check tolerance (relative to the ratio)
        real(dblp),  parameter :: BIG      = huge(1.0_dblp)

        logical,     save :: ready = .false.
        real(dblp),  save :: wgt(nlat)                             !< latitude weights, sum over cells = 1
        real(dblp),  save :: x_last(NB), x_mark(NB)                !< inventory at the last mark (m), end of step
        real(dblp),  save :: p_mo(NB)                              !< precipitation of ec_moisture this step (m)
        real(dblp),  save :: step(NT, NB), month(NT, NB)           !< terms (m)
        real(dblp),  save :: x_end(NB)                             !< inventory at the end of the last step (m)
        real(dblp),  save :: bst(NBS, NS)                          !< bound statistics of the month (fractions)
        integer(ip), save :: nbad(NS), nout(NS)
        real(dblp),  save :: smin(NS), smax(NS)                    !< source range (strict check)
        logical,     save :: have_src = .false.
        integer(ip), save :: nstep = 0, cur_month = -1, cur_year = -1
        integer(ip), save :: ubud(NS) = -1, ubnd(NS) = -1
        logical,     save :: opened = .false.

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Step bracket and marks
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wiso_bud_step_begin()

          use comdyn,      only: pw
          use comphys,     only: rmoisg
          use comemic_mod, only: imonth, iyear

          if (.not. ready) then
            wgt(:) = 2.0_dblp*sqrt(real(nlon, dblp))*pw(:,1)/(2.0_dblp*real(nlon, dblp))   ! as atm_wbudget_mod
            ready = .true.
          endif
          if (imonth /= cur_month .or. iyear /= cur_year) then
            if (nstep > 0) call write_month()
            cur_month = imonth
            cur_year  = iyear
            call reset_month()
          endif
          if (.not. have_src) call init_sources(rmoisg)
          step = 0.0_dblp
          p_mo = 0.0_dblp
          call inventory(rmoisg, x_mark)
          x_last = x_mark

        end subroutine wiso_bud_step_begin

        ! after a stage of the step: change of the inventory since the last mark, minus the fluxes of the stage
        subroutine wiso_bud_mark(stage)

          use comphys,     only: rmoisg, dyrain, dysnow, ivavm, imsink
          use comsurf_mod, only: evap
          use comatm,      only: dtime

          integer(ip), intent(in) :: stage

          real(dblp)  :: x(NB), e(NB)
          integer(ip) :: b, k

          if (.not. ready) return
          call inventory(rmoisg, x)
          select case (stage)
          case (WB_AT_TR)
            step(T_RTR, :) = step(T_RTR, :) + (x - x_last)
          case (WB_AT_MO)
            do b = 1, NB
              k = spec(b)
              e(b)    = real(imsink, dblp)*dtime*gmean(evap(:,:,k))
              p_mo(b) = real(ivavm, dblp)*dtime*gmean(dyrain(:,:,k) + dysnow(:,:,k))
            enddo
            step(T_RMO, :) = step(T_RMO, :) + (x - x_last) - (e - p_mo)
          case (WB_AT_CV)
            ! precipitation of convection: total at the end of the step minus that of ec_moisture
            do b = 1, NB
              e(b) = ptot(spec(b)) - p_mo(b)
            enddo
            step(T_RCV, :) = step(T_RCV, :) + (x - x_last) + e
          end select
          x_last = x

        end subroutine wiso_bud_mark

        subroutine wiso_bud_step_end()

          use comphys,     only: rmoisg, dyrain, dysnow, corain, cosnow, ivavm, imsink
          use comsurf_mod, only: evap
          use comatm,      only: dtime
          use comemic_mod, only: iatm

          real(dblp)  :: x(NB)
          integer(ip) :: b, k

          if (.not. ready) return
          call inventory(rmoisg, x)
          step(T_ROT, :) = x - x_last
          do b = 1, NB
            k = spec(b)
            step(T_E,   b) = real(imsink, dblp)*dtime*gmean(evap(:,:,k))
            step(T_PLS, b) = real(ivavm, dblp)*dtime*gmean(dyrain(:,:,k) + dysnow(:,:,k))
            step(T_PCV, b) = real(ivavm, dblp)*dtime*gmean(corain(:,:,k) + cosnow(:,:,k))
          enddo
          step(T_DX, :) = x - x_mark
          step(T_R,  :) = step(T_DX, :) - (step(T_E, :) - step(T_PLS, :) - step(T_PCV, :))
          month = month + step
          x_end = x
          nstep = nstep + 1
          call bounds_step()

          if (nstep == days_month360d_i*iatm) then
            call write_month()
            call reset_month()
          endif

        end subroutine wiso_bud_step_end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Bounds
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! source range at the first step of the run (restart included): initial vapour of the cells with water
        subroutine init_sources(rmois)

          real(dblp), intent(in) :: rmois(:,:,:)

          integer(ip) :: s
          logical     :: m(nlat, nlon)

          m = rmois(:,:,iwater) > WB_TINY*maxval(rmois(:,:,iwater)) .and. rmois(:,:,iwater) > 0.0_dblp
          do s = 1, NS
            smin(s) = minval(rmois(:,:,iwat16+s-1)/rmois(:,:,iwater), mask = m)
            smax(s) = maxval(rmois(:,:,iwat16+s-1)/rmois(:,:,iwater), mask = m)
          enddo
          have_src = .true.

        end subroutine init_sources

        subroutine bounds_step()

          use comphys,     only: rmoisg, torain, tosnow
          use comsurf_mod, only: evap
          use iso_nml_mod, only: iso_passive, delta_check_min, delta_check_max

          real(dblp)  :: rv(nlat, nlon, NS), rr(nlat, nlon, NS), rs(nlat, nlon, NS), re(nlat, nlon, NS)
          real(dblp)  :: pw_(nlat, nlon), rp18(nlat, nlon), rp(nlat, nlon), x
          logical     :: mv(nlat, nlon), mr(nlat, nlon), ms(nlat, nlon), me(nlat, nlon), mp(nlat, nlon)
          integer(ip) :: s, k

          mv = wet(rmoisg(:,:,iwater))
          mr = wet(torain(:,:,iwater))
          ms = wet(tosnow(:,:,iwater))
          me = evap(:,:,iwater) > 0.0_dblp
          pw_ = torain(:,:,iwater) + tosnow(:,:,iwater)
          mp = wet(pw_)
          do s = 1, NS
            k = iwat16 + s - 1
            rv(:,:,s) = ratio(rmoisg(:,:,k), rmoisg(:,:,iwater), mv)
            rr(:,:,s) = ratio(torain(:,:,k), torain(:,:,iwater), mr)
            rs(:,:,s) = ratio(tosnow(:,:,k), tosnow(:,:,iwater), ms)
            re(:,:,s) = ratio(evap(:,:,k), evap(:,:,iwater), me)

            call minmax(rv(:,:,s) - 1.0_dblp, mv, bst(B_VMIN, s), bst(B_VMAX, s))
            call minmax(rr(:,:,s) - 1.0_dblp, mr, bst(B_RMIN, s), bst(B_RMAX, s))
            call minmax(rs(:,:,s) - 1.0_dblp, ms, bst(B_SMIN, s), bst(B_SMAX, s))

            ! plausibility: non-finite, <= 0, outside the window
            nbad(s) = nbad(s) + count_bad(rv(:,:,s), mv, delta_check_min(k), delta_check_max(k))                   &
                              + count_bad(rr(:,:,s), mr, delta_check_min(k), delta_check_max(k))                   &
                              + count_bad(rs(:,:,s), ms, delta_check_min(k), delta_check_max(k))

            ! strict check (alpha = 1): sources met so far, then the fields against them
            if (iso_passive()) then
              if (any(me)) then
                smin(s) = min(smin(s), minval(re(:,:,s), mask = me))
                smax(s) = max(smax(s), maxval(re(:,:,s), mask = me))
              endif
              call strict(rv(:,:,s), mv, s)
              call strict(rr(:,:,s), mr, s)
              call strict(rs(:,:,s), ms, s)
            endif
          enddo

          ! excesses: d-excess (w2h file) and 17O-excess (w17 file), vapour and total precipitation
          s = iwat2h - iwat16 + 1
          call minmax((rv(:,:,s) - 1.0_dblp) - 8.0_dblp*(rv(:,:,iwat18-iwat16+1) - 1.0_dblp), mv,                    &
                      bst(B_XVMIN, s), bst(B_XVMAX, s))
          rp18 = ratio(torain(:,:,iwat18) + tosnow(:,:,iwat18), pw_, mp)
          rp   = ratio(torain(:,:,iwat2h) + tosnow(:,:,iwat2h), pw_, mp)
          call minmax((rp - 1.0_dblp) - 8.0_dblp*(rp18 - 1.0_dblp), mp, bst(B_XPMIN, s), bst(B_XPMAX, s))
          s = iwat17 - iwat16 + 1
          call minmax(log(max(rv(:,:,s), tiny(x))) - 0.528_dblp*log(max(rv(:,:,iwat18-iwat16+1), tiny(x))), mv,       &
                      bst(B_XVMIN, s), bst(B_XVMAX, s))
          rp = ratio(torain(:,:,iwat17) + tosnow(:,:,iwat17), pw_, mp)
          call minmax(log(max(rp, tiny(x))) - 0.528_dblp*log(max(rp18, tiny(x))), mp, bst(B_XPMIN, s), bst(B_XPMAX, s))

        end subroutine bounds_step

        subroutine strict(r, m, s)

          real(dblp),  intent(in) :: r(:,:)
          logical,     intent(in) :: m(:,:)
          integer(ip), intent(in) :: s

          real(dblp) :: ex(nlat, nlon)

          ex = 0.0_dblp
          where (m) ex = max(smin(s) - r, r - smax(s), 0.0_dblp)/max(abs(r), tiny(1.0_dblp))
          nout(s) = nout(s) + count(m .and. ex > WB_TOL)
          bst(B_WORST, s) = max(bst(B_WORST, s), maxval(ex, mask = m .and. ex > WB_TOL))

        end subroutine strict

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Helpers
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        pure function spec(b) result(k)              ! budget index -> species index (1 = water)
          integer(ip), intent(in) :: b
          integer(ip)             :: k
          k = iwater
          if (b > 1) k = iwat16 + b - 2
        end function spec

        subroutine inventory(rmois, x)
          real(dblp), intent(in)  :: rmois(:,:,:)
          real(dblp), intent(out) :: x(NB)
          integer(ip) :: b
          do b = 1, NB
            x(b) = gmean(rmois(:,:,spec(b)))
          enddo
        end subroutine inventory

        function ptot(k) result(p)                   ! precipitation of the whole step so far (m)
          use comphys, only: dyrain, dysnow, corain, cosnow, ivavm
          use comatm,  only: dtime
          integer(ip), intent(in) :: k
          real(dblp)              :: p
          p = real(ivavm, dblp)*dtime*gmean(dyrain(:,:,k) + dysnow(:,:,k) + corain(:,:,k) + cosnow(:,:,k))
        end function ptot

        function gmean(f) result(gm)
          real(dblp), intent(in) :: f(:,:)
          real(dblp)             :: gm
          integer(ip) :: i
          gm = 0.0_dblp
          do i = 1, nlat
            gm = gm + wgt(i)*sum(f(i,:))
          enddo
        end function gmean

        function wet(q) result(m)
          real(dblp), intent(in) :: q(:,:)
          logical                :: m(size(q, 1), size(q, 2))
          m = q > WB_TINY*maxval(q) .and. q > 0.0_dblp
        end function wet

        function ratio(qi, q, m) result(r)
          real(dblp), intent(in) :: qi(:,:), q(:,:)
          logical,    intent(in) :: m(:,:)
          real(dblp)             :: r(size(q, 1), size(q, 2))
          r = 1.0_dblp
          where (m) r = qi/q
        end function ratio

        subroutine minmax(f, m, fmin, fmax)
          real(dblp), intent(in)    :: f(:,:)
          logical,    intent(in)    :: m(:,:)
          real(dblp), intent(inout) :: fmin, fmax
          if (.not. any(m)) return
          fmin = min(fmin, minval(f, mask = m))
          fmax = max(fmax, maxval(f, mask = m))
        end subroutine minmax

        function count_bad(r, m, dmin, dmax) result(n)
          real(dblp), intent(in) :: r(:,:)
          logical,    intent(in) :: m(:,:)
          real(dblp), intent(in) :: dmin, dmax
          integer(ip)            :: n
          ! r /= r: NaN; abs(r) > BIG: infinite
          n = count(m .and. (r /= r .or. abs(r) > BIG .or. r <= 0.0_dblp .or. r - 1.0_dblp < dmin                   &
                             .or. r - 1.0_dblp > dmax))
        end function count_bad

        subroutine reset_month()
          month = 0.0_dblp
          nstep = 0
          nbad  = 0
          nout  = 0
          bst(B_VMIN:B_SMAX:2, :)  =  BIG
          bst(B_VMAX:B_SMAX:2, :)  = -BIG
          bst(B_XVMIN, :) =  BIG
          bst(B_XVMAX, :) = -BIG
          bst(B_XPMIN, :) =  BIG
          bst(B_XPMAX, :) = -BIG
          bst(B_WORST, :) = 0.0_dblp
        end subroutine reset_month

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Output
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine write_month()

          use comemic_mod, only: fini, irunlabel
          use iso_nml_mod, only: iso_passive

          real(dblp)  :: m(NT, NB), eps, b(NBS), xs
          integer(ip) :: s

          if (.not. opened) call open_files()
          m = month*WB_MM
          do s = 1, NS
            eps = 0.0_dblp
            if (x_end(s+1) > 0.0_dblp .and. x_end(1) > 0.0_dblp) eps = m(T_R, s+1)/(x_end(s+1)*WB_MM)                &
                                                                     - m(T_R, 1)/(x_end(1)*WB_MM)
            write(ubud(s), '(i6,i4,i6,13es13.4)') cur_year, cur_month, nstep, x_end(s+1)*WB_MM, m(T_DX, s+1),        &
                  m(T_E, s+1), m(T_PLS, s+1), m(T_PCV, s+1), m(T_RTR, s+1), m(T_RMO, s+1), m(T_RCV, s+1),             &
                  m(T_ROT, s+1), m(T_R, s+1), eps, x_end(1)*WB_MM, m(T_R, 1)
            flush(ubud(s))

            b = bst(:, s)*1000.0_dblp                         ! deltas, d-excess: permil
            where (abs(bst(:, s)) >= BIG) b = 0.0_dblp        ! no cell met in the month
            xs = 1000.0_dblp
            if (s == iwat17 - iwat16 + 1) xs = 1.0e6_dblp     ! 17O-excess: per meg
            b(B_XVMIN:B_XPMAX) = b(B_XVMIN:B_XPMAX)/1000.0_dblp*xs
            if (iso_passive()) then
              write(ubnd(s), '(i6,i4,i6,10es12.4,i10,i10,es12.3,2es12.4)') cur_year, cur_month, nstep, b(1:10),        &
                    nbad(s), nout(s), bst(B_WORST, s), (smin(s) - 1.0_dblp)*1000.0_dblp,                             &
                    (smax(s) - 1.0_dblp)*1000.0_dblp
            else
              write(ubnd(s), '(i6,i4,i6,10es12.4,i10)') cur_year, cur_month, nstep, b(1:10), nbad(s)
            endif
            flush(ubnd(s))
          enddo

        contains

          subroutine open_files()
            integer(ip) :: t
            do t = 1, NS
              open(newunit=ubud(t), file='outputdata/globals/wisobudget_'//SNAME(t)//fini//'.txt',                  &
                   form='formatted', status='replace', action='write')
              write(ubud(t), '(a,a,a,i0)') '# ECBilt isotopologue budget (wiso_budget_mod), species ', SNAME(t),      &
                   ', global means in mm (rsmow = 1), irunlabel = ', irunlabel
              write(ubud(t), '(a)') '# R = dX - (E - Pls - Pcv) = Rtr + Rmo + Rcv + Rot; eps = R/X - Rw/W '//         &
                   '(relative non-conservation beyond the water)'
              write(ubud(t), '(a)') '# year month nstep    X_end        dX           E            Pls          Pcv'//  &
                   '          Rtr          Rmo          Rcv          Rot          R            eps          W_end'//   &
                   '        Rw'
              open(newunit=ubnd(t), file='outputdata/globals/wisobounds_'//SNAME(t)//fini//'.txt',                  &
                   form='formatted', status='replace', action='write')
              write(ubnd(t), '(a,a,a,i0)') '# ECBilt isotopologue delta bounds (wiso_budget_mod), species ',          &
                   SNAME(t), ', deltas in permil over the month, irunlabel = ', irunlabel
              write(ubnd(t), '(a)') '# xs: w2h d-excess (permil), w17 17O-excess (per meg), 0 otherwise; '//        &
                   'nbad: non-finite, <= 0 or outside the check window; passive: nout cells beyond the '//     &
                   'source range [smin, smax] (rel. 1e-9), worst relative excess'
              write(ubnd(t), '(a)') '# year month nstep   vap_min     vap_max     rain_min    rain_max    snow_min'//   &
                   '    snow_max    xs_vap_min  xs_vap_max  xs_prc_min  xs_prc_max      nbad      nout     worst'//   &
                   '        smin        smax'
            enddo
            opened = .true.
          end subroutine open_files

        end subroutine write_month

#endif

      end module wiso_budget_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
