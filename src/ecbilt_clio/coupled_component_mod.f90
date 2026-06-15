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

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [coupled_component_mod]
!
!>     @brief Abstract contract that every optional coupled component of iLOVECLIM must fulfil.
!
!      DESCRIPTION:
!>     Defines an abstract derived type [coupled_component_t] exposing the lifecycle of a component (init / step / finalize)
!!     and a set of OPTIONAL coupling hooks (import_state / export_state) that default to no-op. A concrete component (FROG,
!!     CARAIB, OCYCC, MEDUSA ...) extends this type in its own module and overrides only what it needs. The orchestrator
!!     (coupler_core_mod) drives components solely through this contract: it never references a specific component.
!
!>     @note  The "phase" argument passed to step() lets a component react to the generic coupling phases (see cartography
!!            document, section 2) without the orchestrator knowing anything component-specific.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module coupled_component_mod

       use global_constants_mod, only: str_len, ip

       implicit none

       private

       public :: coupled_component_t
       public :: PHASE_DAY_BEGIN, PHASE_ATM_STEP, PHASE_LAND_STEP, PHASE_DAY_END, PHASE_YEAR_END

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Generic coupling phase identifiers. A component is queried with the current phase and decides whether it acts.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer(ip), parameter :: PHASE_DAY_BEGIN = 1   ! once, start of day i      (e.g. ocean_to_coupler)
       integer(ip), parameter :: PHASE_ATM_STEP  = 2   ! every atmospheric step j
       integer(ip), parameter :: PHASE_LAND_STEP = 3   ! every land step k
       integer(ip), parameter :: PHASE_DAY_END   = 4   ! once, end of day (j==iatm) (e.g. integrate_ocean, MEDUSA daily)
       integer(ip), parameter :: PHASE_YEAR_END  = 5   ! once per model year       (e.g. FROG, CARAIB, isotopes)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   The abstract contract.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       type, abstract :: coupled_component_t
         character(len=str_len) :: name = "unnamed"      !< human-readable component name (logging)
         logical                :: active = .false.      !< set by the registry when the component is compiled in
       contains
         procedure(comp_init_i),     deferred :: init     !< allocate / read restart / open files
         procedure(comp_step_i),     deferred :: step     !< advance for one coupling occurrence at a given phase
         procedure(comp_finalize_i), deferred :: finalize !< close files / write restart / free memory
         ! --- optional coupling hooks: default no-op implementations provided below ---
         procedure :: import_state => comp_import_noop    !< coupler -> component (bring shared fields in)
         procedure :: export_state => comp_export_noop    !< component -> coupler (push fields out)
         procedure :: wants_phase  => comp_wants_default  !< does this component act at this (phase, day)? default .false.
       end type coupled_component_t

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Abstract interfaces for the deferred procedures.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       abstract interface

         subroutine comp_init_i(self, day, month)
           import :: coupled_component_t, ip
           class(coupled_component_t), intent(inout) :: self
           integer(ip),                intent(in)    :: day, month
         end subroutine comp_init_i

         subroutine comp_step_i(self, phase, iday, jstep, kstep)
           import :: coupled_component_t, ip
           class(coupled_component_t), intent(inout) :: self
           integer(ip),                intent(in)    :: phase          !< one of PHASE_*
           integer(ip),                intent(in)    :: iday           !< current day in the run [1..ntotday]
           integer(ip),                intent(in)    :: jstep, kstep   !< atmospheric / land sub-step indices
         end subroutine comp_step_i

         subroutine comp_finalize_i(self)
           import :: coupled_component_t
           class(coupled_component_t), intent(inout) :: self
         end subroutine comp_finalize_i

       end interface

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Default (no-op) implementations of the optional hooks. A concrete component overrides only those it needs.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine comp_import_noop(self, iday)
         class(coupled_component_t), intent(inout) :: self
         integer(ip),                intent(in)    :: iday
         ! intentionally empty: components with no inbound coupling keep this default
       end subroutine comp_import_noop

       subroutine comp_export_noop(self, iday)
         class(coupled_component_t), intent(inout) :: self
         integer(ip),                intent(in)    :: iday
         ! intentionally empty
       end subroutine comp_export_noop

       logical function comp_wants_default(self, phase, iday) result(wants)
         class(coupled_component_t), intent(in) :: self
         integer(ip),                intent(in) :: phase, iday
         wants = .false.   ! by default a component does nothing; concrete types override with their timer logic
       end function comp_wants_default

      end module coupled_component_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
