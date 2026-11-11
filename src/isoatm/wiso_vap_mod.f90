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
!      MODULE: [wiso_vap_mod]
!
!>     @brief Isotopologues of the atmospheric vapour: ratio memory and transport (ROADMAP step 2c, wiso:2c:0001).
!
!      DESCRIPTION:
!
!>     The isotopologue amounts rmoisg(:,:,iwat16:iwat2h) are kept as R_i * rmoisg(:,:,iwater) (D43): the ratios R_i
!!     (rsmow = 1: isotopologue mass per water mass) are the state, held here in wiso_rvap. Where the water is below
!!     WISO_RATIO_TINY of its field maximum, an amount ratio is noise over noise: the cell keeps its last ratio (the
!!     memory), so that vapour coming back into a dry cell finds a sensible composition.
!!
!!     Transport (D30, D32), in place of the water-only call of ec_moisture:
!!       1. ratios from the amounts (memory where the water is ~0);
!!       2. water with atm_fvmois_step, exactly the call made with WISOATM = 0 (D4 by construction);
!!       3. iwat17, iwat18, iwat2h in ratio mode on a copy of the water before the step, with the same fluxes
!!          (atm_fvmois_ratio: bounded, monotone ratios);
!!       4. iwat16 in density mode, independent of the water (D32; check of the transport);
!!       5. amounts rebuilt as R_i * water; the ratio of iwat16 is updated from its amount.
!!
!!     Sources and sinks of ec_moisture (wiso:2c:0004): evaporation E and large-scale condensation C act together over
!!     the step in the water update (forward Euler, both constant). The isotopologues follow the same process in
!!     continuous time: with q = q0 + (E - C) t, an evaporation flux E_i and a condensate at alpha R,
!!         q dR/dt = E_i - (E + (alpha - 1) C) R,   i.e.   q (E - C) dR/dq = E_i - (E + (alpha - 1) C) R
!!     (from d(Rq)/dt = E_i - alpha R C and dq/dt = E - C; condensate at alpha R removed at once, open system)
!!     hence, with b = E + (alpha - 1) C and R* = E_i / b,
!!         R1 = R* + (R0 - R*) (q1/q0)^(-b/(E - C))
!!     (amounts over the step: E = dtime*imsink*evap, C = dtime*ivavm*vemoisg). Limits: C = 0 gives the mixing
!!     (R0 q0 + E_i)/(q0 + E); E = 0 gives Rayleigh, R0 (q1/q0)^(alpha - 1); E = C the exponential limit (no
!!     branch); q1 -> 0 or q0 ~ 0 give R*. R1 lies between R0 and R*: bounded, whatever the water content. The
!!     condensate follows by mass balance, R0 q0 + E_i - R1 q1, and goes to rain and snow in proportion to the water
!!     ones. The other changes of the water in ec_moisture (negative clip, ec_moisbalance) act at R1 (D7).
!!     alpha = 1 until step 2d (R1 = R0 where E = 0, condensate at the vapour ratio).
!!
!!     Convection (wiso:2c:0005): each iteration of the adjustment loop of ec_convec_two is a pure condensation of
!!     the cell, C = q_before - q_after (E = 0 in the solution above: Rayleigh, R1 = R0 (q1/q0)^(alpha - 1)). The
!!     condensate, R0 q0 - R1 q1, goes to dyrain, dysnow, corain, cosnow in proportion to their water increments in
!!     that iteration (caps and rain/snow choice inherited; DOWNSCALING = 2 cells at the cell level; DOWNSTS levels
!!     without isotopologues). A clip of the vapour at 0 removes all its isotopologues, the precipitation keeps the
!!     ratio of the condensate (D7). wiso_vap_conv_begin / _end bracket each iteration, wiso_vap_conv_total sums
!!     torain, tosnow and rebuilds the vapour amounts.
!!
!!     wiso_vap_follow_water rebuilds the amounts from the ratio memory after a change of the water alone.
!!
!!     Contract
!!       - wiso_vap_transport(rmois): FV mode, once per step, where ec_moisture transports the water; rmois is
!!         rmoisg(nlat,nlon,nwisos), all species updated.
!!       - wiso_vap_follow_water(rmois): rmois(:,:,iwat16:iwat2h) = wiso_rvap * rmois(:,:,iwater).
!!       - The memory is set at the first call from the amounts, r_vap_init (wisoatm.nml) where the water is ~0, or
!!         from the isotope restart (wiso_vap_set_memory, wiso:2c:0002); wiso_vap_get_memory gives it for writing.
!
!      REFERENCES: ROADMAP D4, D30, D32, D43; fvt_ffsl_mod (ratio mode, Zalesak 1979 limiter).
!
!      AUTHORS: dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module wiso_vap_mod

#if ( WISOATM == 1 )

        use global_constants_mod, only: ip, dblp=>dp
        use comatm,               only: nlat, nlon, nwisos, iwater, iwat16, iwat17, iwat2h

        implicit none

        private

        public :: wiso_vap_transport, wiso_vap_follow_water, wiso_vap_get_memory, wiso_vap_set_memory
        public :: wiso_vap_keep_water, wiso_vap_sources, wiso_vap_cont_ratio
        public :: wiso_vap_conv_begin, wiso_vap_conv_end, wiso_vap_conv_total

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Module constants and state
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp),  parameter :: WISO_RATIO_TINY = 1.0e-10_dblp           !< water below this fraction of its maximum:
                                                                           !< the cell keeps its ratio (as fvt_ffsl_mod)
        integer(ip), parameter :: WISO_NR = iwat2h - iwat17 + 1            !< species in ratio mode (17O, 18O, 2H)
        real(dblp),  parameter :: WISO_ALPHA = 1.0_dblp                    !< condensate / vapour: 1 until step 2d

        real(dblp), allocatable, save :: wiso_rvap(:,:,:)                  !< (nlat,nlon,iwat16:iwat2h) vapour ratios
        real(dblp), allocatable, save :: wiso_qw0(:,:)                     !< water before the sources of the step
        real(dblp), allocatable, save :: wiso_dr0(:,:), wiso_ds0(:,:)      !< water rain, snow before ec_moisture adds
        real(dblp),              save :: wiso_cq0, wiso_cp0(4)             !< convection: cell water, dyrain, dysnow,
                                                                           !< corain, cosnow at the start of an iteration

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Transport of all species, in place of the water-only transport
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wiso_vap_transport(rmois)

          use atm_fvmois_mod, only: atm_fvmois_step, atm_fvmois_ratio, atm_fvmois_density

          real(dblp), intent(inout) :: rmois(:,:,:)   !< (nlat,nlon,nwisos) column amounts

          real(dblp)  :: qc(nlat, nlon), rr(nlat, nlon, WISO_NR)

          call wiso_vap_ratios(rmois)
          qc = rmois(:,:,iwater)
          rr = wiso_rvap(:,:,iwat17:iwat2h)

          call atm_fvmois_step(rmois(:,:,iwater))
          call atm_fvmois_ratio(qc, rr)
          call atm_fvmois_density(rmois(:,:,iwat16))

          wiso_rvap(:,:,iwat17:iwat2h) = rr
          call wiso_vap_ratio16(rmois)
          call wiso_vap_follow_water(rmois, only_ratio_mode = .true.)

        end subroutine wiso_vap_transport

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Sources and sinks of ec_moisture: evaporation and large-scale condensation (wiso:2c:0003, 0004)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! water vapour, rain and snow before the sources of the step (ec_moisture, before its condensation loop)
        subroutine wiso_vap_keep_water(rmois, dyr, dys)

          real(dblp), intent(in) :: rmois(:,:,:)   !< (nlat,nlon,nwisos)
          real(dblp), intent(in) :: dyr(:,:,:)     !< (nlat,nlon,nwisos) dyrain
          real(dblp), intent(in) :: dys(:,:,:)     !< (nlat,nlon,nwisos) dysnow

          if (.not. allocated(wiso_qw0)) then
            allocate(wiso_qw0(nlat, nlon), wiso_dr0(nlat, nlon), wiso_ds0(nlat, nlon))
          endif
          wiso_qw0 = rmois(:,:,iwater)
          wiso_dr0 = dyr(:,:,iwater)
          wiso_ds0 = dys(:,:,iwater)

        end subroutine wiso_vap_keep_water

        ! after ec_moisture: evaporation and large-scale condensation together, then the rest of the water changes
        subroutine wiso_vap_sources(rmois, ev, fac, cond, dyr, dys)

          real(dblp), intent(inout) :: rmois(:,:,:)   !< (nlat,nlon,nwisos) column amounts; water already updated
          real(dblp), intent(in)    :: ev(:,:,:)      !< (nlat,nlon,nwisos) evaporation (m/s)
          real(dblp), intent(in)    :: fac            !< dtime*imsink, as in the water update
          real(dblp), intent(in)    :: cond(:,:)      !< (nlat,nlon) water condensed over the step, dtime*ivavm*vemoisg
          real(dblp), intent(inout) :: dyr(:,:,:)     !< (nlat,nlon,nwisos) dyrain: isotopologues added
          real(dblp), intent(inout) :: dys(:,:,:)     !< (nlat,nlon,nwisos) dysnow: isotopologues added

          real(dblp)  :: qtiny, q0, q1, ea, ca, eia, r0, r1, ci, rc, dr, ds
          integer(ip) :: i, j, k

          if (.not. allocated(wiso_qw0)) error stop 'wiso_vap_sources: wiso_vap_keep_water not called'
          if (.not. allocated(wiso_rvap)) call wiso_vap_ratios(rmois)

          qtiny = WISO_RATIO_TINY*max(maxval(wiso_qw0), tiny(1.0_dblp))
          do j = 1, nlon
            do i = 1, nlat
              q0 = wiso_qw0(i, j)
              ea = fac*ev(i, j, iwater)
              ca = cond(i, j)
              q1 = q0 + ea - ca
              dr = dyr(i, j, iwater) - wiso_dr0(i, j)
              ds = dys(i, j, iwater) - wiso_ds0(i, j)
              do k = iwat16, iwat2h
                eia = fac*ev(i, j, k)
                r0  = wiso_rvap(i, j, k)
                if (q0 > qtiny) r0 = rmois(i, j, k)/q0
                r1  = wiso_vap_cont_ratio(q0, r0, ea, eia, ca, WISO_ALPHA, qtiny)
                ! condensate: mass balance; to rain and snow as the water ones
                ci = rmois(i, j, k) + eia - max(q1, 0.0_dblp)*r1
                rc = r1
                if (ca > 0.0_dblp) rc = ci/ca
                dyr(i, j, k) = dyr(i, j, k) + rc*dr
                dys(i, j, k) = dys(i, j, k) + rc*ds
                wiso_rvap(i, j, k) = r1
              enddo
            enddo
          enddo
          call wiso_vap_follow_water(rmois)

        end subroutine wiso_vap_sources

        ! Ratio at the end of a step with constant evaporation (ea, of ratio eia/ea) and condensation (ca, at alpha R)
        ! from q0 at ratio r0 (amounts over the step); see the module description.
        pure function wiso_vap_cont_ratio(q0, r0, ea, eia, ca, alpha, qtiny) result(r1)

          real(dblp), intent(in) :: q0, r0, ea, eia, ca, alpha, qtiny
          real(dblp)             :: r1

          real(dblp)            :: b, rstar, u, g, q1, x

          b  = ea + (alpha - 1.0_dblp)*ca
          q1 = q0 + ea - ca
          if (b == 0.0_dblp) then                      ! nothing changes the ratio (alpha = 1, no evaporation)
            r1 = r0
            return
          endif
          rstar = eia/b
          if (q0 <= qtiny .or. q1 <= 0.0_dblp) then    ! dry start, or everything condensed: the steady state
            r1 = rstar
            return
          endif
          ! (q1/q0)^(-b/(ea-ca)) = exp(-b g/q0), g = ln(1+z)/z, z = (ea-ca)/q0 > -1; g = ln(u)/(u-1) with u = 1+z
          ! as rounded (accurate for any z, Goldberg 1991), g = 1 at u = 1: smooth through ea = ca
          u = 1.0_dblp + (ea - ca)/q0
          if (u == 1.0_dblp) then
            g = 1.0_dblp
          else
            g = log(u)/(u - 1.0_dblp)
          endif
          x  = exp(-b*g/q0)
          r1 = rstar + (r0 - rstar)*x

        end function wiso_vap_cont_ratio

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Convection (wiso:2c:0005): one call pair per iteration of the adjustment loop of ec_convec_two
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wiso_vap_conv_begin(qw, dyr, dys, cor, cos)

          real(dblp), intent(in) :: qw                   !< water vapour of the cell
          real(dblp), intent(in) :: dyr, dys, cor, cos   !< water dyrain, dysnow, corain, cosnow of the cell

          wiso_cq0 = qw
          wiso_cp0 = [dyr, dys, cor, cos]

        end subroutine wiso_vap_conv_begin

        subroutine wiso_vap_conv_end(i, j, qw, dyr, dys, cor, cos)

          integer(ip), intent(in)    :: i, j                    !< cell
          real(dblp),  intent(in)    :: qw                      !< water vapour of the cell after the iteration
          real(dblp),  intent(inout) :: dyr(:), dys(:), cor(:), cos(:)   !< (nwisos) of the cell: isotopologues added

          real(dblp)  :: c, dp(4), qtiny, r0, r1, ci, rc
          integer(ip) :: k

          if (.not. allocated(wiso_rvap)) error stop 'wiso_vap_conv_end: no ratio memory (ec_moisture not called)'
          c = wiso_cq0 - qw
          if (c <= 0.0_dblp) return                     ! nothing condensed in this iteration: ratio unchanged
          dp = [dyr(iwater), dys(iwater), cor(iwater), cos(iwater)] - wiso_cp0
          qtiny = WISO_RATIO_TINY*max(maxval(wiso_qw0), tiny(1.0_dblp))
          do k = iwat16, iwat2h
            r0 = wiso_rvap(i, j, k)
            r1 = wiso_vap_cont_ratio(wiso_cq0, r0, 0.0_dblp, 0.0_dblp, c, WISO_ALPHA, qtiny)
            ci = r0*wiso_cq0 - r1*max(qw, 0.0_dblp)
            rc = ci/c
            dyr(k) = dyr(k) + rc*dp(1)
            dys(k) = dys(k) + rc*dp(2)
            cor(k) = cor(k) + rc*dp(3)
            cos(k) = cos(k) + rc*dp(4)
            wiso_rvap(i, j, k) = r1
          enddo

        end subroutine wiso_vap_conv_end

        ! end of ec_convec_two: total precipitation of the isotopologues, vapour amounts from the ratios
        subroutine wiso_vap_conv_total(rmois, tor, tos, dyr, dys, cor, cos)

          real(dblp), intent(inout) :: rmois(:,:,:)              !< (nlat,nlon,nwisos)
          real(dblp), intent(inout) :: tor(:,:,:), tos(:,:,:)    !< torain, tosnow: isotopologues set
          real(dblp), intent(in)    :: dyr(:,:,:), dys(:,:,:), cor(:,:,:), cos(:,:,:)

          integer(ip) :: k

          do k = iwat16, iwat2h
            tor(:,:,k) = dyr(:,:,k) + cor(:,:,k)
            tos(:,:,k) = dys(:,:,k) + cos(:,:,k)
          enddo
          call wiso_vap_follow_water(rmois)

        end subroutine wiso_vap_conv_total

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Amounts from the ratio memory
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wiso_vap_follow_water(rmois, only_ratio_mode)

          real(dblp),        intent(inout) :: rmois(:,:,:)      !< (nlat,nlon,nwisos) column amounts
          logical, optional, intent(in)    :: only_ratio_mode   !< leave iwat16 (transported as a density) alone

          integer(ip) :: k, k0

          if (.not. allocated(wiso_rvap)) call wiso_vap_ratios(rmois)
          k0 = iwat16
          if (present(only_ratio_mode)) then
            if (only_ratio_mode) k0 = iwat17
          endif
          do k = k0, iwat2h
            rmois(:,:,k) = wiso_rvap(:,:,k)*rmois(:,:,iwater)
          enddo

        end subroutine wiso_vap_follow_water

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Ratio memory for the restart (wiso:2c:0002)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wiso_vap_get_memory(rmois, r)

          real(dblp), intent(in)  :: rmois(:,:,:)   !< (nlat,nlon,nwisos), to set the memory if not done yet
          real(dblp), intent(out) :: r(:,:,:)       !< (nlat,nlon,iwat2h-iwat16+1)

          if (.not. allocated(wiso_rvap)) call wiso_vap_ratios(rmois)
          r = wiso_rvap

        end subroutine wiso_vap_get_memory

        subroutine wiso_vap_set_memory(r)

          real(dblp), intent(in) :: r(:,:,:)        !< (nlat,nlon,iwat2h-iwat16+1)

          if (.not. allocated(wiso_rvap)) allocate(wiso_rvap(nlat, nlon, iwat16:iwat2h))
          wiso_rvap = r

        end subroutine wiso_vap_set_memory

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Ratio memory from the amounts
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! All species; allocates and initialises the memory at the first call (r_vap_init where the water is ~0)
        subroutine wiso_vap_ratios(rmois)

          use iso_nml_mod, only: r_vap_init

          real(dblp), intent(in) :: rmois(:,:,:)

          integer(ip) :: k

          if (.not. allocated(wiso_rvap)) then
            allocate(wiso_rvap(nlat, nlon, iwat16:iwat2h))
            do k = iwat16, iwat2h
              wiso_rvap(:,:,k) = r_vap_init(k)
            enddo
          endif
          do k = iwat16, iwat2h
            call ratio_where_wet(rmois(:,:,iwater), rmois(:,:,k), wiso_rvap(:,:,k))
          enddo

        end subroutine wiso_vap_ratios

        ! iwat16 only (after its density transport)
        subroutine wiso_vap_ratio16(rmois)

          real(dblp), intent(in) :: rmois(:,:,:)

          call ratio_where_wet(rmois(:,:,iwater), rmois(:,:,iwat16), wiso_rvap(:,:,iwat16))

        end subroutine wiso_vap_ratio16

        ! r = qi / q where q > WISO_RATIO_TINY * max(q); r unchanged elsewhere
        subroutine ratio_where_wet(q, qi, r)

          real(dblp), intent(in)    :: q(:,:)
          real(dblp), intent(in)    :: qi(:,:)
          real(dblp), intent(inout) :: r(:,:)

          real(dblp) :: qtiny

          qtiny = WISO_RATIO_TINY*maxval(q)
          where (q > qtiny .and. q > 0.0_dblp) r = qi/q

        end subroutine ratio_where_wet

#endif

      end module wiso_vap_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
