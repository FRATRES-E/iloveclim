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
!!     wiso_vap_follow_water rebuilds the amounts from the ratio memory after a change of the water alone. Until the
!!     isotopologue sources and sinks exist (wiso:2c:0002 to 0004), every water source and sink is applied this way
!!     at the local vapour ratio (alpha = 1, no source composition): a uniform ratio must then stay uniform in a full
!!     model run, a check of the transport and of the bookkeeping. Each later step replaces one of these calls.
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

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Module constants and state
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp),  parameter :: WISO_RATIO_TINY = 1.0e-10_dblp           !< water below this fraction of its maximum:
                                                                           !< the cell keeps its ratio (as fvt_ffsl_mod)
        integer(ip), parameter :: WISO_NR = iwat2h - iwat17 + 1            !< species in ratio mode (17O, 18O, 2H)

        real(dblp), allocatable, save :: wiso_rvap(:,:,:)                  !< (nlat,nlon,iwat16:iwat2h) vapour ratios

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
