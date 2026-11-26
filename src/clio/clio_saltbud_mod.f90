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
!      MODULE: [clio_saltbud_mod]
!
!>     @brief Global salt budget of CLIO split by routine and by term (ROADMAP 3b, wiso:3b:0002).
!
!      DESCRIPTION:
!
!>     Pure diagnostics: nothing here changes the model state. The global salt content sum(S*aire*dz*tms) is taken
!!     at marks placed in the CLIO step (clio.f), so that its change is split by stage:
!!         1 icdyna, icdadv, thersf    2 brines    3 diftur ... isoadvec    4 scale    5 scali
!!         6 etat (convective adjustment)    7 the rest (barotropic mode, icebergs, outputs, and the forcing of the
!!         next step, up to its icdyna)
!!     Inside scale (salt only), the explicit terms are summed with the formulas of scale itself, on the points that
!!     scale updates (js1..js2, is1..is2, times tms):
!!         srf  surface flux phiss(:,:,2)          rap  surface restoring rappes*(S - scalr)
!!         bot  bottom outflow cczdt*phifs(:,:,2)   der  global drift correction -deriv*dts (upper ocean)
!!         phv  volume sources -phivs               int  internal restoring rapint*(scalr - S)
!!     with the surface flux further split (wiso:3b:0003): brine of sea-ice growth (fons > 0) and melt (fons < 0),
!!     from fluxbrines, and the ice-shelf melt compensation ficess (thersf, south iceberg zone), which should be the
!!     opposite of phv (the ice-shelf melt at depth); srf - brine - ics is what is left of the surface flux;
!!     and the totals that scale applies: Xs = sum of scaldt*tms (all explicit terms; Xs - the six = what the
!!     horizontal schemes put into scaldt) and As = sum of (scalat - S)*tms (advection and diffusion, should be 0).
!!     All in freshwater equivalent, global means in mm over the CLIO earth (4 pi rterre^2): -dSalt/S*rho0,
!!     positive = salt lost = freshening. One line per 30 steps in outputdata/globals/clio_saltbud.txt.
!!     Sea ice and snow on sea ice (wiso:3b:0004, 0005): mass at the marks
!!         0 start of the step    1 after icdyna    2 after icdadv    3 in thersf, after its section 4 (particular
!!         cases)    4 after thersf
!!     so that the changes split into icdyna, icdadv (transport), particular cases, thermodynamics and the rest, with
!!     the thermodynamic mass changes that thersf reports (dmgbq: ice, the basis of the brine flux fons; dmgwi: snow
!!     ice; dmnbq: snow). Masses are weighted by the sea-ice model's own cell area (dynami_mod area, m2), which the
!!     transport conserves; dmgbq is also given weighted by the ocean area aire (the weight of the brine salt flux in
!!     scale): a difference between the two is a salt imbalance from the two metrics. The ratio area/aire (min, max
!!     over the ocean cells) is written in the header. Mm of water over the CLIO earth, one line per 30 steps in
!!     outputdata/globals/clio_icebud.txt.
!!     Inside icdadv (wiso:3b:0006): ice and snow volume (times rhog, rhon) per hemisphere (south j < jeq, north
!!     j >= jeq) on the ocean cells of the update (js1..js2, is1..is2, times tms) at the sub-stages of the routine:
!!         E   extensive fields of 2.2 against the state at entry (should be 0)
!!         A   advection (adv or advx/advy), Ac the part of it from the clip at zero inside adv (iameth = 1 only)
!!         D   diffusion (diffus, before the clip)    Dc  the clip at zero after diffus
!!         L   section 2.5 (clip of the means, concentration cap, removal of thin ice and snow)
!!     and the mean ice mass on the land cells (tms = 0) of the same rows after the advection and at the end.
!!     One line per 30 steps in outputdata/globals/clio_icdadv.txt; E + A + D + Dc + L = I_adv (S_adv) of clio_icebud.
!!     The diffusion of the ice volume (wiso:3b:0007): diffus solves fld1 = fld0 + dt*(alfa*div(fld1) + (1-alfa)*div(fld0))
!!     (Crank-Nicolson, alfa = 0.5) by relaxed Jacobi iterations; from fld0 and the returned fld1 the face fluxes are
!!     rebuilt as diffus builds them (grxh, gryh times dxs2, dxs1), and the mass change D on the updated cells splits into
!!         T    sum of area*dt*(alfa*div(fld1) + (1-alfa)*div(fld0)): what the discrete equation moves; it is the flux
!!              through faces to cells that diffus does not update (not ocean, outside js1..js2, is1..is2; the cyclic
!!              halo counted as its partner), 0 for a closed domain
!!         Def  D - T: the defect of the unconverged iteration
!!     per hemisphere (T of one hemisphere includes the flux across the equator; T_s + T_n is the leak), with the
!!     number of iterations (mean, max) and the last change per iteration (max), the mass change
!!     in the Bering raccord cells (iberp/iberpm, jberp and ibera/iberam, jbera, where they are updated), and the five
!!     cells with the largest cumulative |Def| and |flux to non-updated cells| since the start of the run.
!!     outputdata/globals/clio_icdiff.txt (one line per 30 steps) and clio_icdiff_top.txt (cumulative top cells).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module clio_saltbud_mod

        use global_constants_mod, only: ip, dblp=>dp

        implicit none

        private

        public :: sb_mark, sb_scale_bc, sb_scale_k, sb_end_step, sb_icesh, im_mark
        public :: ia_state, ia_ext, ia_addclip, ia_cur, ia_diffus, ia_dfinfo

        integer(ip), parameter :: NST = 7                      !< stages
        integer(ip), parameter :: NT  = 11                     !< scale terms: srf rap bot der phv int Xs As brf brm ics
        integer(ip), parameter :: T_SRF = 1, T_RAP = 2, T_BOT = 3, T_DER = 4, T_PHV = 5, T_INT = 6, T_XS = 7, T_AS = 8
        integer(ip), parameter :: T_BRF = 9, T_BRM = 10, T_ICS = 11
        integer(ip), parameter :: NSTEP_OUT = 30

        logical,     save :: have_last = .false.
        real(dblp),  save :: last      = 0.0_dblp              !< salt content at the last mark (psu * 1e12 m3)
        real(dblp),  save :: sal_mean  = 35.0_dblp             !< volume-mean salinity at the last mark
        real(dblp),  save :: d(NST)    = 0.0_dblp
        real(dblp),  save :: t(NT)     = 0.0_dblp
        integer(ip), save :: nstep     = 0
        integer(ip), save :: unit      = -1
        logical,     save :: opened    = .false.
        real(dblp),  save :: ficess_l  = 0.0_dblp              !< last ice-shelf salt compensation of thersf (psu/step)
        integer(ip), parameter :: NIS = 5                      !< ice stages: icdyna, icdadv, particular, thermo, rest
        integer(ip), parameter :: NIM = 4                      !< thermodynamic masses: dmgbq, dmgwi, dmnbq, dmgbq on aire
        logical,     save :: ice_have  = .false.
        real(dblp),  save :: ice_last(2) = 0.0_dblp           !< ice and snow mass at the last ice mark (kg / 1e12)
        real(dblp),  save :: di(2,NIS) = 0.0_dblp
        real(dblp),  save :: dm(NIM)   = 0.0_dblp
        integer(ip), save :: unit_i    = -1
        logical,     save :: opened_i  = .false.
        integer(ip), parameter :: NIA = 6                      !< icdadv marks: 0 entry, 1 2.2, 2 adv, 3 dif, 4 clip, 5 end
        integer(ip), save :: ia_cur    = 0                     !< species whose adv clip is counted (1 ice, 2 snow, 0 none)
        real(dblp),  save :: ia_m(2,2,0:NIA-1) = 0.0_dblp     !< last mass at each mark (species, hemisphere)
        real(dblp),  save :: ia_dd(2,2,NIA)    = 0.0_dblp     !< accumulated: E A D Dc L, Ac (species, hemisphere)
        real(dblp),  save :: ia_land(2,2)      = 0.0_dblp     !< ice mass on land cells (hemisphere; after adv, end)
        integer(ip), save :: ia_n      = 0                     !< icdadv calls in the interval
        integer(ip), save :: unit_a    = -1
        logical,     save :: opened_a  = .false.
        integer(ip), parameter :: NTOP = 5
        real(dblp),  save :: id_d(2,3)  = 0.0_dblp            !< D, T, Def (hemisphere)
        real(dblp),  save :: id_ber     = 0.0_dblp            !< mass change of the Bering raccord cells
        real(dblp),  save :: id_amx     = 0.0_dblp            !< max last change per iteration
        integer(ip), save :: id_kit     = 0, id_kmax = 0, id_n = 0, id_kcur = 0
        real(dblp),  save :: id_amxcur  = 0.0_dblp
        real(dblp),  allocatable, save :: id_cdef(:,:), id_clk(:,:)   !< cumulative per-cell defect and leak (kg)
        integer(ip), save :: unit_d = -1, unit_t = -1
        logical,     save :: opened_d = .false.


      contains

        subroutine sb_mark(stage)                    ! stage 1..6 after a routine; 0 at the start of a step (closes 7)

          integer(ip), intent(in) :: stage

          real(dblp) :: s

          s = content()
          if (have_last) then
            if (stage == 0) then
              d(NST) = d(NST) + (s - last)
            else
              d(stage) = d(stage) + (s - last)
            endif
          endif
          last      = s
          have_last = .true.

        end subroutine sb_mark

        subroutine im_mark(stage)    ! sea ice, snow: 0 start, 1 after icdyna, 2 after icdadv, 3 thersf section 4, 4 after thersf

          use bloc0_mod,   only: tms, is1, is2, js1, js2, ks2
          use bloc_mod,    only: aire
          use dynami_mod,  only: area
          use ice_mod,     only: hgbq, hnbq, albq, rhog, rhon, dmgbq, dmgwi, dmnbq

          integer(ip), intent(in) :: stage

          real(dblp)  :: m(2), a
          integer(ip) :: i, j, k

          m = 0.0_dblp
          do j = js1, js2
            do i = is1(j), is2(j)
              a    = (1.0_dblp - albq(i,j))*area(i,j)*tms(i,j,ks2)
              m(1) = m(1) + a*hgbq(i,j)*rhog
              m(2) = m(2) + a*hnbq(i,j)*rhon
            enddo
          enddo
          if (ice_have) then
            k = stage
            if (stage == 0) k = NIS
            di(:,k) = di(:,k) + (m - ice_last)
          endif
          if (stage == 4) then
            do j = js1, js2
              do i = is1(j), is2(j)
                a     = area(i,j)*tms(i,j,ks2)
                dm(1) = dm(1) + dmgbq(i,j)*a
                dm(2) = dm(2) + dmgwi(i,j)*a
                dm(3) = dm(3) + dmnbq(i,j)*a
                dm(4) = dm(4) + dmgbq(i,j)*aire(i,j)*1.0e12_dblp*tms(i,j,ks2)
              enddo
            enddo
          endif
          ice_last = m
          ice_have = .true.

        end subroutine im_mark

        ! icdadv: masses from the model state (hgbq, hnbq, albq) at entry (k = 0) or at the end (k = 5)
        subroutine ia_state(k)

          use bloc0_mod,   only: tms, is1, is2, js1, js2, ks2, jeq
          use dynami_mod,  only: area
          use ice_mod,     only: hgbq, hnbq, albq, rhog, rhon

          integer(ip), intent(in) :: k

          integer(ip) :: i, j, h
          real(dblp)  :: a, m(2,2), ml(2)

          m  = 0.0_dblp
          ml = 0.0_dblp
          do j = js1, js2
            h = 1
            if (j >= jeq) h = 2
            do i = is1(j), is2(j)
              a      = (1.0_dblp - albq(i,j))*area(i,j)
              m(1,h) = m(1,h) + a*hgbq(i,j)*rhog*tms(i,j,ks2)
              m(2,h) = m(2,h) + a*hnbq(i,j)*rhon*tms(i,j,ks2)
              ml(h)  = ml(h)  + a*hgbq(i,j)*rhog*(1.0_dblp - tms(i,j,ks2))
            enddo
          enddo
          ia_m(:,:,k) = m
          if (k == 0) then
            ia_n = ia_n + 1
          else
            ia_dd(:,:,5)  = ia_dd(:,:,5) + (m - ia_m(:,:,4))      ! L
            ia_land(:,2)  = ia_land(:,2) + ml
          endif

        end subroutine ia_state

        ! icdadv: one species (isp 1 ice, 2 snow) of the transported fields at mark k (1 2.2, 2 adv, 3 dif, 4 clip);
        ! fld extensive (volume per cell, ext) or per unit area (times area)
        subroutine ia_ext(k, isp, fld, ext)

          use para0_mod,   only: imax, jmax
          use bloc0_mod,   only: tms, is1, is2, js1, js2, ks2, jeq
          use dynami_mod,  only: area
          use ice_mod,     only: rhog, rhon

          integer(ip), intent(in) :: k, isp
          real(dblp),  intent(in) :: fld(imax,jmax)
          logical,     intent(in) :: ext

          integer(ip) :: i, j, h
          real(dblp)  :: rho, v, m(2), ml(2)

          rho = rhog
          if (isp == 2) rho = rhon
          m  = 0.0_dblp
          ml = 0.0_dblp
          do j = js1, js2
            h = 1
            if (j >= jeq) h = 2
            do i = is1(j), is2(j)
              v = fld(i,j)
              if (.not. ext) v = v*area(i,j)
              m(h)  = m(h)  + v*rho*tms(i,j,ks2)
              ml(h) = ml(h) + v*rho*(1.0_dblp - tms(i,j,ks2))
            enddo
          enddo
          ia_m(isp,:,k)   = m
          ia_dd(isp,:,k)  = ia_dd(isp,:,k) + (m - ia_m(isp,:,k-1))   ! k = 1 E, 2 A, 3 D, 4 Dc
          if (k == 2 .and. isp == 1) ia_land(:,1) = ia_land(:,1) + ml

        end subroutine ia_ext

        ! adv (iameth = 1): the clip at zero of the corrector call, for the species set in ia_cur
        subroutine ia_addclip(i, j, s)

          use bloc0_mod,   only: tms, ks2, jeq
          use ice_mod,     only: rhog, rhon

          integer(ip), intent(in) :: i, j
          real(dblp),  intent(in) :: s                        !< value before the clip (volume per cell)

          integer(ip) :: h
          real(dblp)  :: rho

          if (ia_cur == 0 .or. s >= 0.0_dblp) return
          rho = rhog
          if (ia_cur == 2) rho = rhon
          h = 1
          if (j >= jeq) h = 2
          ia_dd(ia_cur,h,6) = ia_dd(ia_cur,h,6) - s*rho*tms(i,j,ks2)

        end subroutine ia_addclip

        ! diffus: iterations used and last change per iteration
        subroutine ia_dfinfo(k, amx)

          integer(ip), intent(in) :: k
          real(dblp),  intent(in) :: amx

          id_kcur   = k
          id_amxcur = amx

        end subroutine ia_dfinfo

        ! icdadv: diffusion of the ice volume per unit area, fld0 -> fld1 (diffus, before the clip), split into the flux
        ! through faces to non-updated cells (T) and the defect of the iteration (Def)
        subroutine ia_diffus(dt, difhx, difhy, fld0, fld1)

          use para0_mod,   only: imax, jmax
          use bloc0_mod,   only: tms, is1, is2, js1, js2, ks2, jeq, jcl1, jcl2, ims1, ims2, iberp, jberp, ibera, jbera,  &
                                 iberpm, iberam
          use dynami_mod,  only: area, dxs1, dxs2
          use ice_mod,     only: rhog

          real(dblp), intent(in) :: dt
          real(dblp), intent(in) :: difhx(imax,jmax), difhy(imax,jmax), fld0(imax,jmax), fld1(imax,jmax)

          real(dblp), parameter :: ALFA = 0.5_dblp
          real(dblp)  :: fx(imax,jmax), fy(imax,jmax), tc, dc, f
          logical     :: inu(imax,jmax)
          integer(ip) :: i, j, h, ip1, jp1, ii

          if (.not. allocated(id_cdef)) then
            allocate(id_cdef(imax,jmax), id_clk(imax,jmax))
            id_cdef = 0.0_dblp
            id_clk  = 0.0_dblp
          endif
          inu = .false.
          do j = js1, js2
            do i = is1(j), is2(j)
              inu(i,j) = tms(i,j,ks2) > 0.0_dblp
            enddo
          enddo
          ! face fluxes (volume per step) as diffus: grxh(i,j) on the face (i,j)|(i+1,j), gryh on (i,j)|(i,j+1)
          fx = 0.0_dblp
          fy = 0.0_dblp
          do j = js1-1, js2
            jp1 = j + 1
            do i = is1(j)-1, is2(j)
              ip1 = i + 1
              fx(i,j) = dt*dxs2(i,j)*difhx(i,j)*(ALFA*(fld1(ip1,j) - fld1(i,j)) + (1.0_dblp-ALFA)*(fld0(ip1,j) - fld0(i,j)))
              fy(i,j) = dt*dxs1(i,j)*difhy(i,j)*(ALFA*(fld1(i,jp1) - fld1(i,j)) + (1.0_dblp-ALFA)*(fld0(i,jp1) - fld0(i,j)))
            enddo
          enddo
          do j = js1, js2
            h = 1
            if (j >= jeq) h = 2
            do i = is1(j), is2(j)
              if (.not. inu(i,j)) cycle
              tc = (fx(i,j) - fx(i-1,j) + fy(i,j) - fy(i,j-1))*rhog                 ! as bkappa*area
              dc = area(i,j)*(fld1(i,j) - fld0(i,j))*rhog - tc
              id_d(h,1) = id_d(h,1) + area(i,j)*(fld1(i,j) - fld0(i,j))*rhog
              id_d(h,2) = id_d(h,2) + tc
              id_d(h,3) = id_d(h,3) + dc
              id_cdef(i,j) = id_cdef(i,j) + dc
              if ((i == iberp .or. i == iberpm) .and. j == jberp) id_ber = id_ber + tc + dc
              if ((i == ibera .or. i == iberam) .and. j == jbera) id_ber = id_ber + tc + dc
            enddo
          enddo
          ! unmatched faces: one side updated, the other not (the cyclic halo stands for its partner)
          do j = js1-1, js2
            jp1 = j + 1
            do i = is1(j)-1, is2(j)
              ip1 = i + 1
              f  = fx(i,j)*rhog
              ii = ip1
              if (j >= jcl1 .and. j <= jcl2 .and. ip1 == ims2+1) ii = ims1
              if (j >= jcl1 .and. j <= jcl2 .and. i == ims1-1) then
                ! face ims1-1 | ims1 is the face ims2 | ims2+1 seen from the other side
                if (inu(ims1,j) .and. .not. inu(ims2,j)) id_clk(ims1,j) = id_clk(ims1,j) - f
              else if (inu(i,j) .and. .not. inu(ii,j)) then
                id_clk(i,j) = id_clk(i,j) + f
              else if (inu(ii,j) .and. .not. inu(i,j)) then
                id_clk(ii,j) = id_clk(ii,j) - f
              endif
              if (jp1 <= jmax) then
                f = fy(i,j)*rhog
                if (inu(i,j) .and. .not. inu(i,jp1)) then
                  id_clk(i,j) = id_clk(i,j) + f
                else if (inu(i,jp1) .and. .not. inu(i,j)) then
                  id_clk(i,jp1) = id_clk(i,jp1) - f
                endif
              endif
            enddo
          enddo
          id_n    = id_n + 1
          id_kit  = id_kit + id_kcur
          id_kmax = max(id_kmax, id_kcur)
          id_amx  = max(id_amx, id_amxcur)

        end subroutine ia_diffus

        subroutine write_icdiff()

          use const_mod,   only: pi, rterre
          use bloc_mod,    only: numit
          use bloc0_mod,   only: iberp, jberp, ibera, jbera
          use para0_mod,   only: imax, jmax

          real(dblp)  :: f, w(imax,jmax)
          integer(ip) :: n, loc(2), td(2,NTOP), tl(2,NTOP)
          real(dblp)  :: vd(NTOP), vl(NTOP)

          if (.not. allocated(id_cdef)) return
          if (.not. opened_d) then
            open(newunit=unit_d, file='outputdata/globals/clio_icdiff.txt', form='formatted', status='replace',       &
                 action='write')
            write(unit_d, '(a)') '# CLIO ice-volume diffusion in icdadv (clio_saltbud_mod), mm of water over 30 steps on'//&
                                 ' the updated ocean cells, south (s) and north (n): D mass change (= I_D of'//            &
                                 ' clio_icdadv), T what the equation moves (T_s + T_n: flux through faces to'//          &
                                 ' non-updated cells; each alone includes the flux across the equator), Def = D - T'//    &
                                 ' defect of the'//                                                                         &
                                 ' unconverged iteration; Ber: mass change of the Bering raccord cells; iterations'//     &
                                 ' (mean, max) and the largest last change per iteration (m)'
            write(unit_d, '(a,4i5)') '# Bering raccord cells iberp, jberp, ibera, jbera:', iberp, jberp, ibera, jbera
            write(unit_d, '(a)') '#  numit  D_s        D_n        T_s        T_n        Def_s      Def_n      Ber'//         &
                                 '        it_mean   it_max  amx_max'
            open(newunit=unit_t, file='outputdata/globals/clio_icdiff_top.txt', form='formatted', status='replace',   &
                 action='write')
            write(unit_t, '(a)') '# five cells (i, j, kg/1e12 cumulative since the start) with the largest |Def| and'//    &
                                 ' the largest |flux to non-updated cells| of the ice-volume diffusion'
            opened_d = .true.
          endif
          f = 1.0_dblp/(4.0_dblp*pi*rterre**2)
          write(unit_d, '(i8,7es11.3,f9.2,i7,es11.3)') numit, id_d(:,1)*f, id_d(:,2)*f, id_d(:,3)*f, id_ber*f,           &
                real(id_kit, dblp)/real(max(id_n,1), dblp), id_kmax, id_amx
          flush(unit_d)
          w = abs(id_cdef)
          do n = 1, NTOP
            loc = maxloc(w)
            td(:,n) = loc
            vd(n)   = id_cdef(loc(1),loc(2))*1.0e-12_dblp
            w(loc(1),loc(2)) = -1.0_dblp
          enddo
          w = abs(id_clk)
          do n = 1, NTOP
            loc = maxloc(w)
            tl(:,n) = loc
            vl(n)   = id_clk(loc(1),loc(2))*1.0e-12_dblp
            w(loc(1),loc(2)) = -1.0_dblp
          enddo
          write(unit_t, '(i8,a,5(2i4,es11.3),a,5(2i4,es11.3))') numit, ' Def:', (td(:,n), vd(n), n = 1, NTOP),          &
                ' Leak:', (tl(:,n), vl(n), n = 1, NTOP)
          flush(unit_t)
          id_d    = 0.0_dblp
          id_ber  = 0.0_dblp
          id_amx  = 0.0_dblp
          id_kit  = 0
          id_kmax = 0
          id_n    = 0

        end subroutine write_icdiff

        function content() result(s)

          use para0_mod, only: kmax
          use bloc0_mod, only: tms, is1, is2, js1, js2, scal, dz
          use bloc_mod,  only: aire

          real(dblp)  :: s, v, vol
          integer(ip) :: i, j, k

          s   = 0.0_dblp
          vol = 0.0_dblp
          do k = 1, kmax
            do j = js1, js2
              do i = is1(j), is2(j)
                v   = aire(i,j)*dz(k)*tms(i,j,k)
                vol = vol + v
                s   = s + v*scal(i,j,k,2)
              enddo
            enddo
          enddo
          sal_mean = s/vol

        end function content

        subroutine sb_icesh(ficess)                  ! from thersf: the ice-shelf salt compensation of this step

          real(dblp), intent(in) :: ficess

          ficess_l = ficess

        end subroutine sb_icesh

        ! the explicit boundary terms of scale for the salt, as scale builds them (section 3), on the updated points
        subroutine sb_scale_bc(cczdt, k1000, k2000, j50S)

          use para0_mod, only: kmax
          use bloc0_mod, only: tms, is1, is2, js1, js2, ks2, scal, scalr, dz, dts, phiss, phivs, phifs, rappes,     &
                               kfs, scpme, spvr
          use bloc_mod,  only: aire, deriv
          use para0_mod, only: imax, jmax
          use ice_mod,   only: iicebers1, iicebers2, jicebers1, jicebers2

          real(dblp),  intent(in) :: cczdt(kmax)
          integer(ip), intent(in) :: k1000, k2000, j50S

          integer(ip) :: i, j, k, kfd
          real(dblp)  :: v, sscal

          real*8 fluxbrines(imax,jmax)                 ! as in thersf.f, brines_mod.f90 (fons*unsdz/rho0, psu/step)
          common / common_brines / fluxbrines

          do j = js1, js2
            do k = 1, kmax
              if (scpme(2) == spvr) then
                sscal = -deriv(2)*dts(k)
              else if ((j <= j50S) .and. (k >= k2000)) then
                sscal = -deriv(2)*dts(k)
              else if (k >= k1000) then
                sscal = -deriv(2)*dts(k)
              else
                sscal = 0.0_dblp
              endif
              do i = is1(j), is2(j)
                v = aire(i,j)*dz(k)*tms(i,j,k)
                t(T_DER) = t(T_DER) + sscal*v
                t(T_PHV) = t(T_PHV) - phivs(i,j,k,2)*v
              enddo
            enddo
            do i = is1(j), is2(j)
              v = aire(i,j)*dz(ks2)*tms(i,j,ks2)
              t(T_SRF) = t(T_SRF) - phiss(i,j,2)*v
              t(T_RAP) = t(T_RAP) - rappes(i,j,2)*(scal(i,j,ks2,2) - scalr(i,j,ks2,2))*v
              if (fluxbrines(i,j) > 0.0_dblp) then
                t(T_BRF) = t(T_BRF) + fluxbrines(i,j)*v          ! salt gained = fons*unsdz/rho0
              else
                t(T_BRM) = t(T_BRM) + fluxbrines(i,j)*v
              endif
              if (i >= iicebers1 .and. i <= iicebers2 .and. j >= jicebers1 .and. j <= jicebers2)                     &
                t(T_ICS) = t(T_ICS) + ficess_l*v
              kfd = kfs(i,j)
              t(T_BOT) = t(T_BOT) + cczdt(kfd)*phifs(i,j,2)*aire(i,j)*dz(kfd)*tms(i,j,kfd)
            enddo
          enddo

        end subroutine sb_scale_bc

        ! level k of the final update of scale for the salt (section 7): internal restoring, totals applied
        subroutine sb_scale_k(k, scalat, scaldt)

          use para0_mod, only: imax, jmax, kmax
          use bloc0_mod, only: tms, is1, is2, js1, js2, scal, scalr, dz, rapint
          use bloc_mod,  only: aire, nrap, ijrap

          integer(ip), intent(in) :: k
          real(dblp),  intent(in) :: scalat(imax,jmax,kmax), scaldt(imax,jmax,kmax)

          integer(ip) :: i, j, n
          real(dblp)  :: v

          do n = 1, nrap(k,2)
            i = 1 + mod(ijrap(n,k,2), imax)
            j = 1 + ijrap(n,k,2)/imax
            t(T_INT) = t(T_INT) + rapint(n,k,2)*(scalr(i,j,k,2) - scal(i,j,k,2))*aire(i,j)*dz(k)*tms(i,j,k)
          enddo
          do j = js1, js2
            do i = is1(j), is2(j)
              v = aire(i,j)*dz(k)*tms(i,j,k)
              t(T_XS) = t(T_XS) + scaldt(i,j,k)*tms(i,j,k)*v
              t(T_AS) = t(T_AS) + (scalat(i,j,k) - scal(i,j,k,2))*v
            enddo
          enddo

        end subroutine sb_scale_k

        subroutine sb_end_step()                     ! after the last mark of the step (etat)

          nstep = nstep + 1
          if (nstep == NSTEP_OUT) then
            call write_out()
            call write_ice()
            call write_icdadv()
            call write_icdiff()
            d     = 0.0_dblp
            t     = 0.0_dblp
            di    = 0.0_dblp
            dm    = 0.0_dblp
            ia_dd   = 0.0_dblp
            ia_land = 0.0_dblp
            ia_n    = 0
            nstep = 0
          endif

        end subroutine sb_end_step

        subroutine write_out()

          use const_mod, only: rho0, pi, rterre
          use bloc_mod,  only: numit

          real(dblp) :: f, dd(NST), tt(NT), dtot, tsum, rest

          if (.not. opened) then
            open(newunit=unit, file='outputdata/globals/clio_saltbud.txt', form='formatted', status='replace',       &
                 action='write')
            write(unit, '(a)') '# CLIO global salt budget (clio_saltbud_mod), freshwater equivalent, global means in'//  &
                               ' mm over 30 steps (-dSalt/S*rho0 over 4 pi rterre^2; positive = salt lost)'
            write(unit, '(a)') '# stages: thersf brines pre-scale(diftur..isoadvec) scale scali etat rest; total;'//      &
                               ' scale terms: srf rap bot der phv int, Xs (scaldt applied), Xs-sum, As (scalat-S),'//   &
                               ' scale - Xs - As; surface split: brine freeze, brine melt, ice-shelf compensation,'//   &
                               ' srf - brine - ics'
            write(unit, '(a)') '#  numit  S_mean     thersf     brines     prescl     scale      scali      etat'//    &
                               '       rest       total      srf        rap        bot        der        phv'//       &
                               '        int        Xs         Xs-sum     As         sc-Xs-As   brf        brm'//      &
                               '        ics        srf-rest'
            opened = .true.
          endif

          f    = -1.0e12_dblp/sal_mean*rho0/(4.0_dblp*pi*rterre**2)
          dd   = d*f
          tt   = t*f
          dtot = sum(dd)
          tsum = sum(tt(T_SRF:T_INT))
          rest = tt(T_SRF) - tt(T_BRF) - tt(T_BRM) - tt(T_ICS)
          write(unit, '(i8,f9.4,22es11.3)') numit, sal_mean, dd, dtot, tt(T_SRF:T_INT), tt(T_XS), tt(T_XS) - tsum,      &
                                          tt(T_AS), dd(4) - tt(T_XS) - tt(T_AS), tt(T_BRF), tt(T_BRM), tt(T_ICS), rest
          flush(unit)

        end subroutine write_out

        subroutine write_ice()

          use const_mod,   only: pi, rterre
          use bloc_mod,    only: numit, aire
          use bloc0_mod,   only: tms, is1, is2, js1, js2, ks2
          use dynami_mod,  only: area

          real(dblp)  :: f, rmin, rmax, r
          integer(ip) :: i, j

          if (.not. opened_i) then
            rmin = huge(1.0_dblp)
            rmax = -huge(1.0_dblp)
            do j = js1, js2
              do i = is1(j), is2(j)
                if (tms(i,j,ks2) > 0.0_dblp .and. aire(i,j) > 0.0_dblp) then
                  r    = area(i,j)/(aire(i,j)*1.0e12_dblp)
                  rmin = min(rmin, r)
                  rmax = max(rmax, r)
                endif
              enddo
            enddo
            open(newunit=unit_i, file='outputdata/globals/clio_icebud.txt', form='formatted', status='replace',      &
                 action='write')
            write(unit_i, '(a)') '# CLIO sea-ice and snow mass budget (clio_saltbud_mod), global means in mm of water'// &
                                 ' over 30 steps, weighted by the sea-ice area (dynami_mod); stages: icdyna, icdadv,'//   &
                                 ' part (thersf section 4), thermo (rest of thersf), rest; thersf masses dmgbq (ice),'//  &
                                 ' dmgwi (snow ice), dmnbq (snow); dmgbq_aire: dmgbq weighted by the ocean area aire'
            write(unit_i, '(a,2es14.6)') '# area/aire over the ocean cells, min and max:', rmin, rmax
            write(unit_i, '(a)') '#  numit  I_dyna     I_adv      I_part     I_thermo   I_rest     I_total    dmgbq'//    &
                                 '      I_th-dmgbq dmgbq_aire dmgwi      S_dyna     S_adv      S_part     S_thermo'//   &
                                 '   S_rest     S_total    dmnbq      S_th-dmnbq'
            opened_i = .true.
          endif
          f = 1.0_dblp/(4.0_dblp*pi*rterre**2)
          write(unit_i, '(i8,18es11.3)') numit, di(1,:)*f, sum(di(1,:))*f, dm(1)*f, (di(1,4) - dm(1))*f, dm(4)*f,       &
                                        dm(2)*f, di(2,:)*f, sum(di(2,:))*f, dm(3)*f, (di(2,4) - dm(3))*f
          flush(unit_i)

        end subroutine write_ice

        subroutine write_icdadv()

          use const_mod,   only: pi, rterre
          use bloc_mod,    only: numit

          real(dblp)  :: f, fl
          integer(ip) :: isp, h

          if (.not. opened_a) then
            open(newunit=unit_a, file='outputdata/globals/clio_icdadv.txt', form='formatted', status='replace',      &
                 action='write')
            write(unit_a, '(a)') '# CLIO icdadv sub-stages (clio_saltbud_mod), ice (I) and snow (S) on the ocean cells,'//  &
                                 ' south (s) and north (n), global means in mm of water over 30 steps: E entry, A'//      &
                                 ' advection, Ac its clip at zero (part of A), D diffusion, Dc clip after diffusion,'//    &
                                 ' L section 2.5; E+A+D+Dc+L = I_adv; Lnd_a, Lnd_e: mean ice mass on land cells (mm)'//    &
                                 ' after advection and at the end'
            write(unit_a, '(a)') '#  numit  I_E_s      I_A_s      I_Ac_s     I_D_s      I_Dc_s     I_L_s      I_E_n'//    &
                                 '      I_A_n      I_Ac_n     I_D_n      I_Dc_n     I_L_n      S_E_s      S_A_s'//       &
                                 '      S_Ac_s     S_D_s      S_Dc_s     S_L_s      S_E_n      S_A_n      S_Ac_n'//      &
                                 '     S_D_n      S_Dc_n     S_L_n      Lnd_a_s    Lnd_a_n    Lnd_e_s    Lnd_e_n'
            opened_a = .true.
          endif
          f  = 1.0_dblp/(4.0_dblp*pi*rterre**2)
          fl = f/real(max(ia_n, 1), dblp)
          write(unit_a, '(i8,28es11.3)') numit,                                                                        &
                ((ia_dd(isp,h,1)*f, ia_dd(isp,h,2)*f, ia_dd(isp,h,6)*f, ia_dd(isp,h,3)*f, ia_dd(isp,h,4)*f,           &
                  ia_dd(isp,h,5)*f, h = 1, 2), isp = 1, 2), ia_land(:,1)*fl, ia_land(:,2)*fl
          flush(unit_a)

        end subroutine write_icdadv

      end module clio_saltbud_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
