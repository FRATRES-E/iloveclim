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
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module clio_saltbud_mod

        use global_constants_mod, only: ip, dblp=>dp

        implicit none

        private

        public :: sb_mark, sb_scale_bc, sb_scale_k, sb_end_step, sb_icesh, im_mark

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
            d     = 0.0_dblp
            t     = 0.0_dblp
            di    = 0.0_dblp
            dm    = 0.0_dblp
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

      end module clio_saltbud_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
