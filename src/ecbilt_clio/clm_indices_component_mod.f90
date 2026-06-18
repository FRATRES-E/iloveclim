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

#if ( CLM_INDICES >= 1 )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [clm_indices_component_mod]
!
!>     @brief Adapter exposing the climate-indices diagnostics (CLM_INDICES) as a coupled component.
!
!      DESCRIPTION:
!>     The component now has a daily step (originally hidden in atmdiag0.f), plus setup and teardown.
!!       init     : GLOBAL_RE_INIT()              -- emic.f L488 (after init_coupled_components / topo).
!!       step     : DAILYSTEP_FOR_CLIM_INDICES()  -- atmdiag0.f L167, once/day. Rapatriated here from ec_selectout.
!!       finalize : GLOBAL_FINALIZE()             -- emic.f L992.
!!
!>     @note The daily step was buried in ec_selectout (atmdiag0.f), called via ec_ecbilt -> ec_atmout -> ec_selectout and
!!           gated by mod(istep,iatm)==0 (once per day). It reads end-of-day atmospheric state from modules, takes no args, so
!!           it relocates cleanly to step(PHASE_BEFORE_OCEAN). Registered FIRST so its step precedes FROG/OCYCC (faithful to
!!           the original order, where DAILYSTEP ran inside ec_ecbilt, before the j==iatm FROG/OCYCC block).
!!
!>     @note Name clash to be aware of: GLOBAL_RE_INIT is exported BOTH by CLIMATE_INDICES_MOD and by GRID_IO_NC (the latter
!!           is aliased grid_io_reinit in emic.f L146). Here we import ONLY the CLIMATE_INDICES_MOD one, so no collision.
!>     @note Registration order: this component must be registered LAST so its finalize runs after all others (emic.f L992
!!           is the last finalize). finalize_all iterates the registry in registration order.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module clm_indices_component_mod

       use global_constants_mod,  only: ip
       use coupled_component_mod, only: coupled_component_t, PHASE_BEFORE_OCEAN

       implicit none

       private

       public :: clm_indices_component_t, clm_indices_component

       type, extends(coupled_component_t) :: clm_indices_component_t
       contains
         procedure :: init        => clm_init
         procedure :: step        => clm_step
         procedure :: finalize    => clm_finalize
         procedure :: wants_phase => clm_wants_phase
       end type clm_indices_component_t

! dmr --- module-level target instance, referenced by the static registry
       type(clm_indices_component_t), save, target :: clm_indices_component

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   init: GLOBAL_RE_INIT (set up the climate-index accumulators). emic.f L488.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine clm_init(self, day, month)
         use CLIMATE_INDICES_MOD, only: GLOBAL_RE_INIT
         class(clm_indices_component_t), intent(inout) :: self
         integer(ip),                    intent(in)    :: day, month

         self%name = "CLM_INDICES"
         call GLOBAL_RE_INIT()

       end subroutine clm_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   wants_phase: act once per day, end of day, before the ocean. Mirrors the original call site in ec_selectout, which
! dmr   ec_atmout gates with mod(istep,iatm)==0 (i.e. after the last atmospheric step of the day).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical function clm_wants_phase(self, phase, iday) result(wants)
         class(clm_indices_component_t), intent(in) :: self
         integer(ip),                    intent(in) :: phase, iday
         wants = (phase == PHASE_BEFORE_OCEAN)
       end function clm_wants_phase

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   step: DAILYSTEP_FOR_CLIM_INDICES -- daily accumulation of climate indices. Formerly buried in ec_selectout
! dmr   (atmdiag0.f L167), reached via ec_ecbilt -> ec_atmout -> ec_selectout once per day. Reads the end-of-day atmospheric
! dmr   state from modules (comdyn/comphys/comatm), takes no arguments. Placed in BEFORE_OCEAN per non-regression check.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine clm_step(self, phase, iday, jstep, kstep)
         use CLIMATE_INDICES_MOD, only: DAILYSTEP_FOR_CLIM_INDICES
         class(clm_indices_component_t), intent(inout) :: self
         integer(ip),                    intent(in)    :: phase, iday, jstep, kstep
         call DAILYSTEP_FOR_CLIM_INDICES()
       end subroutine clm_step

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   finalize: GLOBAL_FINALIZE (write out the indices). emic.f L992 -- last finalize.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine clm_finalize(self)
         use CLIMATE_INDICES_MOD, only: GLOBAL_FINALIZE
         class(clm_indices_component_t), intent(inout) :: self
         call GLOBAL_FINALIZE()
       end subroutine clm_finalize

      end module clm_indices_component_mod

#endif /* CLM_INDICES >= 1 */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
