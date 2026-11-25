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
!      MODULE: [wisolbm_budget_mod]
!
!>     @brief Water isotope budget and delta-bound diagnostics of the LBM land stores (ROADMAP 3a, wiso:3a:0005).
!
!      DESCRIPTION:
!
!>     Pure diagnostics: nothing here changes the model state. Active with land_iso_mode = 'bucket' only (lbm_wiso).
!!     Called by the LBM driver around each land step (lbm_landsurface), on the land points (fractl > epsl), weights
!!     dareas*fractl over the area of the globe: global means in mm, as the atmospheric budget (wiso_budget_mod), so that
!!     the two can be added. One pair of files per species s = w16, w17, w18, w2h, one line per model month:
!!
!!     outputdata/globals/wisolbm_budget_<s><fini>.txt (rsmow = 1: mm of water at the standard ratio):
!!         L    = B + S, the land inventory (soil water bmoisg + snow dsnow) at the end of the month; B and S also
!!         dL   change over the month
!!         P    rain and snow fall on land (rainf, snowf); E land evaporation (evapl); Q runoff of the points (runofl)
!!         R    = dL - (P - E - Q): what the land processes did not conserve (the clip of negative soil water, and
!!              round-off); the water columns (W, Rw) are repeated: R of an isotope should follow the water's
!!         eps  = R/L - Rw/W (non-conservation beyond the water's, relative)
!!
!!     outputdata/globals/wisolbm_bounds_<s><fini>.txt, deltas in permil over the month, land points with water above
!!     1e-10 of the field maximum: min / max delta of the soil water, the snow, the runoff and the land evaporation;
!!     w2h: d-excess, w17: 17O-excess (per meg) of the soil water and of the snow, min / max; nbad: point-steps with a
!!     non-finite ratio, a ratio <= 0, or a delta outside [delta_check_min, max] of wisoatm.nml (w16: |delta| > 1e-6).
!!
!!     Contract (public entry points):
!!       wlb_step_begin() : lbm_landsurface, before the loop over the points.
!!       wlb_step_end()   : lbm_landsurface, after it (the runoff of the points is known: lbm_landprecip).
!
!      REFERENCES: ROADMAP 3a, D49 (atmospheric budget), D64, D67, D68, O58.
!
!      AUTHORS: dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module wisolbm_budget_mod

#if ( WISO == 1 )

        use global_constants_mod, only: ip, dblp=>dp, days_month360d_i
        use comatm,               only: nlat, nlon, iwater, iwat16, iwat17, iwat18, iwat2h

        implicit none

        private

        public :: wlb_step_begin, wlb_step_end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Constants and state
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip), parameter :: NS = iwat2h - iwat16 + 1                     !< species
        integer(ip), parameter :: NB = NS + 1                                  !< budgets: water (1) and species (2:)
        character(len=3), parameter :: SNAME(NS) = ['w16', 'w17', 'w18', 'w2h']

        ! budget terms
        integer(ip), parameter :: T_P = 1, T_E = 2, T_Q = 3, T_R = 4, T_DL = 5, NT = 5
        ! bound statistics: min / max of soil water, snow, runoff, evaporation deltas; excesses of soil water and snow
        integer(ip), parameter :: B_BMIN = 1, B_BMAX = 2, B_SMIN = 3, B_SMAX = 4, B_QMIN = 5, B_QMAX = 6,         &
                                  B_EMIN = 7, B_EMAX = 8, B_XBMIN = 9, B_XBMAX = 10, B_XSMIN = 11, B_XSMAX = 12, NBS = 12
        real(dblp),  parameter :: WLB_MM   = 1000.0_dblp          !< m -> mm
        real(dblp),  parameter :: WLB_TINY = 1.0e-10_dblp         !< water below this fraction of its max: no ratio
        real(dblp),  parameter :: BIG      = huge(1.0_dblp)

        logical,     save :: ready = .false.
        real(dblp),  save :: wgt(nlat, nlon)                       !< dareas*fractl/area of the globe on land points
        logical,     save :: lnd(nlat, nlon)                       !< land points (fractl > epsl)
        real(dblp),  save :: l_begin(NB)                           !< inventory at the beginning of the step (m)
        real(dblp),  save :: month(NT, NB)                         !< terms of the month (m)
        real(dblp),  save :: b_end(NB), s_end(NB)                  !< soil water and snow at the end of the last step (m)
        real(dblp),  save :: bst(NBS, NS)                          !< bound statistics of the month (fractions)
        integer(ip), save :: nbad(NS)
        integer(ip), save :: nstep = 0, cur_month = -1, cur_year = -1
        integer(ip), save :: ubud(NS) = -1, ubnd(NS) = -1
        logical,     save :: opened = .false.

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Step bracket
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wlb_step_begin()

          use comland_mod, only: bmoisg, dsnow, lbm_wiso
          use comemic_mod, only: imonth, iyear

          real(dblp) :: b(NB), s(NB)

          if (.not. lbm_wiso) return
          call setup()                                   ! every step: fractl may change (ec_landcoverupdate)
          if (imonth /= cur_month .or. iyear /= cur_year) then
            if (nstep > 0) call write_month()
            cur_month = imonth
            cur_year  = iyear
            call reset_month()
          endif
          call inventory(bmoisg, b)
          call inventory(dsnow, s)
          l_begin = b + s

        end subroutine wlb_step_begin

        subroutine wlb_step_end()

          use comland_mod, only: bmoisg, dsnow, rainf, snowf, evapl, runofl, dtland, lbm_wiso
          use comemic_mod, only: iatm, ilan

          real(dblp)  :: b(NB), s(NB), st(NT, NB)

          if (.not. lbm_wiso) return
          if (.not. ready) return
          call inventory(bmoisg, b)
          call inventory(dsnow, s)
          call inventory(rainf + snowf, st(T_P, :))
          call inventory(evapl, st(T_E, :))
          call inventory(runofl, st(T_Q, :))
          st(T_P:T_Q, :) = dtland*st(T_P:T_Q, :)
          st(T_DL, :) = (b + s) - l_begin
          st(T_R,  :) = st(T_DL, :) - (st(T_P, :) - st(T_E, :) - st(T_Q, :))
          month = month + st
          b_end = b
          s_end = s
          nstep = nstep + 1
          call bounds_step()

          if (nstep == days_month360d_i*iatm*ilan) then
            call write_month()
            call reset_month()
          endif

        end subroutine wlb_step_end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Bounds
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine bounds_step()

          use comland_mod, only: bmoisg, dsnow, runofl, evapl
          use iso_nml_mod, only: delta_check_min, delta_check_max

          real(dblp)  :: rb(nlat, nlon, NS), rs(nlat, nlon, NS), rq(nlat, nlon, NS), re(nlat, nlon, NS), x
          logical     :: mb(nlat, nlon), ms(nlat, nlon), mq(nlat, nlon), me(nlat, nlon)
          integer(ip) :: s, k, s18, s17, s2h

          mb = wet(bmoisg(:,:,iwater))
          ms = wet(dsnow(:,:,iwater))
          mq = wet(runofl(:,:,iwater))
          me = wet(evapl(:,:,iwater))
          do s = 1, NS
            k = iwat16 + s - 1
            rb(:,:,s) = ratio(bmoisg(:,:,k), bmoisg(:,:,iwater), mb)
            rs(:,:,s) = ratio(dsnow(:,:,k), dsnow(:,:,iwater), ms)
            rq(:,:,s) = ratio(runofl(:,:,k), runofl(:,:,iwater), mq)
            re(:,:,s) = ratio(evapl(:,:,k), evapl(:,:,iwater), me)
            call minmax(rb(:,:,s) - 1.0_dblp, mb, bst(B_BMIN, s), bst(B_BMAX, s))
            call minmax(rs(:,:,s) - 1.0_dblp, ms, bst(B_SMIN, s), bst(B_SMAX, s))
            call minmax(rq(:,:,s) - 1.0_dblp, mq, bst(B_QMIN, s), bst(B_QMAX, s))
            call minmax(re(:,:,s) - 1.0_dblp, me, bst(B_EMIN, s), bst(B_EMAX, s))
            nbad(s) = nbad(s) + count_bad(rb(:,:,s), mb, delta_check_min(k), delta_check_max(k))                   &
                              + count_bad(rs(:,:,s), ms, delta_check_min(k), delta_check_max(k))                   &
                              + count_bad(rq(:,:,s), mq, delta_check_min(k), delta_check_max(k))                   &
                              + count_bad(re(:,:,s), me, delta_check_min(k), delta_check_max(k))
          enddo

          ! excesses of the soil water and the snow: d-excess (w2h file), 17O-excess (w17 file)
          s18 = iwat18 - iwat16 + 1
          s17 = iwat17 - iwat16 + 1
          s2h = iwat2h - iwat16 + 1
          call minmax((rb(:,:,s2h) - 1.0_dblp) - 8.0_dblp*(rb(:,:,s18) - 1.0_dblp), mb, bst(B_XBMIN, s2h), bst(B_XBMAX, s2h))
          call minmax((rs(:,:,s2h) - 1.0_dblp) - 8.0_dblp*(rs(:,:,s18) - 1.0_dblp), ms, bst(B_XSMIN, s2h), bst(B_XSMAX, s2h))
          call minmax(log(max(rb(:,:,s17), tiny(x))) - 0.528_dblp*log(max(rb(:,:,s18), tiny(x))), mb,                &
                      bst(B_XBMIN, s17), bst(B_XBMAX, s17))
          call minmax(log(max(rs(:,:,s17), tiny(x))) - 0.528_dblp*log(max(rs(:,:,s18), tiny(x))), ms,                &
                      bst(B_XSMIN, s17), bst(B_XSMAX, s17))

        end subroutine bounds_step

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Helpers
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine setup()

          use comland_mod, only: dareas, fractl, epsl

          real(dblp)  :: area
          integer(ip) :: i

          area = real(nlon, dblp)*sum(dareas(1:nlat))
          lnd = fractl > epsl
          do i = 1, nlat
            wgt(i,:) = merge(dareas(i)*fractl(i,:)/area, 0.0_dblp, lnd(i,:))
          enddo
          ready = .true.

        end subroutine setup

        pure function spec(b) result(k)              ! budget index -> species index (1 = water)
          integer(ip), intent(in) :: b
          integer(ip)             :: k
          k = iwater
          if (b > 1) k = iwat16 + b - 2
        end function spec

        subroutine inventory(f, x)                     ! weighted sum over the land points, per budget
          real(dblp), intent(in)  :: f(:,:,:)
          real(dblp), intent(out) :: x(NB)
          integer(ip) :: b
          do b = 1, NB
            x(b) = sum(wgt*f(:,:,spec(b)))
          enddo
        end subroutine inventory

        function wet(q) result(m)
          real(dblp), intent(in) :: q(:,:)
          logical                :: m(size(q, 1), size(q, 2))
          m = lnd .and. q > WLB_TINY*maxval(q, mask = lnd) .and. q > 0.0_dblp
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
          bst(1:NBS:2, :) =  BIG
          bst(2:NBS:2, :) = -BIG
        end subroutine reset_month

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Output
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine write_month()

          use comemic_mod, only: fini, irunlabel

          real(dblp)  :: m(NT, NB), eps, b(NBS), xs, l_end(NB)
          integer(ip) :: s

          if (.not. opened) call open_files()
          m = month*WLB_MM
          l_end = (b_end + s_end)*WLB_MM
          do s = 1, NS
            eps = 0.0_dblp
            if (l_end(s+1) > 0.0_dblp .and. l_end(1) > 0.0_dblp) eps = m(T_R, s+1)/l_end(s+1) - m(T_R, 1)/l_end(1)
            write(ubud(s), '(i6,i4,i6,12es13.4)') cur_year, cur_month, nstep, l_end(s+1), b_end(s+1)*WLB_MM,          &
                  s_end(s+1)*WLB_MM, m(T_DL, s+1), m(T_P, s+1), m(T_E, s+1), m(T_Q, s+1), m(T_R, s+1), eps,           &
                  l_end(1), m(T_R, 1), m(T_P, 1)
            flush(ubud(s))

            b = bst(:, s)*1000.0_dblp                         ! deltas, d-excess: permil
            where (abs(bst(:, s)) >= BIG) b = 0.0_dblp        ! no point met in the month
            xs = 1000.0_dblp
            if (s == iwat17 - iwat16 + 1) xs = 1.0e6_dblp     ! 17O-excess: per meg
            b(B_XBMIN:B_XSMAX) = b(B_XBMIN:B_XSMAX)/1000.0_dblp*xs
            write(ubnd(s), '(i6,i4,i6,12es12.4,i10)') cur_year, cur_month, nstep, b(1:NBS), nbad(s)
            flush(ubnd(s))
          enddo

        contains

          subroutine open_files()
            integer(ip) :: t
            do t = 1, NS
              open(newunit=ubud(t), file='outputdata/globals/wisolbm_budget_'//SNAME(t)//fini//'.txt',               &
                   form='formatted', status='replace', action='write')
              write(ubud(t), '(a,a,a,i0)') '# LBM land isotope budget (wisolbm_budget_mod), species ', SNAME(t),      &
                   ', global means in mm (rsmow = 1), land points, irunlabel = ', irunlabel
              write(ubud(t), '(a)') '# L = B + S (soil water, snow); R = dL - (P - E - Q); eps = R/L - Rw/W '//      &
                   '(relative non-conservation beyond the water); W, Rw, Pw: the water'
              write(ubud(t), '(a)') '# year month nstep    L_end        B_end        S_end        dL           P'//    &
                   '            E            Q            R            eps          W_end        Rw           Pw'
              open(newunit=ubnd(t), file='outputdata/globals/wisolbm_bounds_'//SNAME(t)//fini//'.txt',               &
                   form='formatted', status='replace', action='write')
              write(ubnd(t), '(a,a,a,i0)') '# LBM land isotope delta bounds (wisolbm_budget_mod), species ',          &
                   SNAME(t), ', deltas in permil over the month, irunlabel = ', irunlabel
              write(ubnd(t), '(a)') '# soil water (bm), snow (sn), runoff (ro), land evaporation (ev); xs: w2h '//   &
                   'd-excess (permil), w17 17O-excess (per meg), 0 otherwise; nbad: non-finite, <= 0 or outside '//  &
                   'the check window'
              write(ubnd(t), '(a)') '# year month nstep   bm_min      bm_max      sn_min      sn_max      ro_min'//    &
                   '      ro_max      ev_min      ev_max      xs_bm_min   xs_bm_max   xs_sn_min   xs_sn_max     nbad'
            enddo
            opened = .true.
          end subroutine open_files

        end subroutine write_month

#endif

      end module wisolbm_budget_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
