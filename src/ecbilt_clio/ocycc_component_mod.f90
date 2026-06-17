!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2026 iLOVECLIM / LUDUS coding group

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#include "choixcomposantes.h"

#if ( OCYCC == 1 )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [ocycc_component_mod]
!
!>     @brief Orchestration adapter exposing OCYCC (ocean carbon cycle) as a coupled component.
!
!      DESCRIPTION:
!>     OCYCC is already self-contained in module ocycc_main: ocycc_ini() and ocycc_step() each handle their OWN grid sync with
!!     CLIO internally (sync_lcm_ocycc(1) on entry to read the CLIO-advected tracers from scal(:,:,:,3:nsmax), sync_lcm_ocycc(2)
!!     on exit to deposit them back). The biogeochemical tracers are advected transparently BY CLIO inside scal; OCYCC never
!!     advects. There is therefore nothing to disentangle: this adapter only routes the existing entry points through the
!!     generic contract, so the OCYCC #if branches leave init_coupled_components / emic.f and live only in the registry.
!!
!!       init                  : ocycc_ini(day,month)   -- emic.f init_coupled_components, last call AFTER init_clio.
!!       step(PHASE_BEFORE_OCEAN): ocycc_step(iday,imonth) -- emic.f L599, once per day, BEFORE ec_co2oc/clio.
!!       finalize              : no-op (see note below on out_cycc / CYCC).
!!
!>     @note OCYCC has NO sub-cadence: it runs every day. wants_phase(PHASE_BEFORE_OCEAN) is always .true.
!>     @note CYCC is a SEPARATE model: the TERRESTRIAL carbon cycle, unrelated to OCYCC (oceanic). They are not variants of
!!           one another; out_cycc (emic.f L959) merely happens to be a poorly-factored output routine mixing global and
!!           OCYCC concerns. It is deliberately left in emic.f, not pulled into this oceanic component.
!>     @note ocycc_ini must run after CLIO init (it reads scal via sync_lcm_ocycc(1)); placing it in init_all AFTER
!!           init_coupled_components preserves the original ordering exactly.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module ocycc_component_mod

       use global_constants_mod,  only: ip
       use coupled_component_mod, only: coupled_component_t, PHASE_BEFORE_OCEAN

       implicit none

       private

       public :: ocycc_component_t, ocycc_component

       type, extends(coupled_component_t) :: ocycc_component_t
       contains
         procedure :: init        => ocycc_init
         procedure :: step        => ocycc_step_phase
         procedure :: finalize    => ocycc_finalize
         procedure :: wants_phase => ocycc_wants_phase
       end type ocycc_component_t

! dmr --- module-level target instance, referenced by the static registry
       type(ocycc_component_t), save, target :: ocycc_component

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   init: delegate to ocycc_ini, which internally does fait_pointer / sync_lcm_ocycc(1) / INIT_OCC / sync_lcm_ocycc(2).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine ocycc_init(self, day, month)
         use ocycc_main, only: ocycc_ini
         class(ocycc_component_t), intent(inout) :: self
         integer(ip),              intent(in)    :: day, month
         logical :: ok

         self%name = "OCYCC"
         ok = ocycc_ini(day, month)

       end subroutine ocycc_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   wants_phase: OCYCC acts every day, before the ocean integration.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical function ocycc_wants_phase(self, phase, iday) result(wants)
         class(ocycc_component_t), intent(in) :: self
         integer(ip),              intent(in) :: phase, iday
         wants = (phase == PHASE_BEFORE_OCEAN)
       end function ocycc_wants_phase

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   step: the daily biogeochemical integration. ocycc_step handles its own CLIO<->OCYCC sync sandwich internally.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine ocycc_step_phase(self, phase, iday, jstep, kstep)
         use comemic_mod, only: cal_iday=>iday, imonth
         use ocycc_main,  only: ocycc_step
         class(ocycc_component_t), intent(inout) :: self
         integer(ip),              intent(in)    :: phase, iday, jstep, kstep
         logical :: ok

         ! NB: ocycc_step / sync_timer expect the CALENDAR day-of-month (1..30) and month, held in comemic_mod by the ECBilt
         ! core -- NOT the absolute run day (the dispatch argument `iday`, 1..ntotday). emic.f L599 passed comemic's iday.
         ok = ocycc_step(cal_iday, imonth)

       end subroutine ocycc_step_phase

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   finalize: no-op for now. The carbon-cycle output routine out_cycc (emic.f L959, guarded by CYCC>=2) is intentionally
! dmr   NOT pulled in here: CYCC is the TERRESTRIAL carbon cycle (distinct from OCYCC, the oceanic one), and out_cycc mixes
! dmr   global and OCYCC-specific output. Moving it into this component would give OCYCC a false responsibility over
! dmr   terrestrial-carbon output. It stays in emic.f until that routine is untangled.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine ocycc_finalize(self)
         class(ocycc_component_t), intent(inout) :: self
         ! intentionally empty: see note above on out_cycc / CYCC
       end subroutine ocycc_finalize

      end module ocycc_component_mod

#endif /* OCYCC == 1 */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
