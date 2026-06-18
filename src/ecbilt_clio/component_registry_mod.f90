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
!      MODULE: [component_registry_mod]
!
!>     @brief Static registry of the optional coupled components that are compiled in.
!
!      DESCRIPTION:
!>     This is the ONLY module in the new coupler core that contains component-specific #if branches. It builds, at startup,
!!     a fixed-size array of pointers to the concrete component instances that are active for this build. The orchestrator
!!     then iterates over this array blindly. Decision actee: STATIC registry (no runtime/dynamic registration). Adding a
!!     component = adding a USE + a pointer association here, nothing in the orchestrator or the main program.
!
!>     @note component_ptr_t is PUBLIC so the orchestrator can declare iterators over registry_components without needing the
!!           concrete component types in scope.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module component_registry_mod

       use global_constants_mod,  only: ip
       use coupled_component_mod, only: coupled_component_t

! -----------------------------------------------------------------------------------------------------------------------------
! dmr  Pull in the concrete component modules, guarded by their activation flags. Each provides a `target` instance of a type
! dmr  extending coupled_component_t. FROG is the pilot; uncomment when frog_component_mod is wired in.
! -----------------------------------------------------------------------------------------------------------------------------

#if ( CLM_INDICES >= 1 )
       use clm_indices_component_mod, only: clm_indices_component   ! registered FIRST (step precedes FROG/OCYCC)
#endif
#if ( FROG_EXP > 0 )
       use frog_component_mod,   only: frog_component   ! type(frog_component_t), target, save
#endif
#if ( MEDUSA == 1 )
!      use medusa_component_mod, only: medusa_component
#endif
#if ( ISM == 2 || ISM == 3 )
!      use grisli_component_mod, only: grisli_component
#endif
#if ( CARAIB > 0 )
!      use caraib_component_mod, only: caraib_component
#endif
#if ( OCYCC == 1 )
       use ocycc_component_mod,  only: ocycc_component
#endif
#if ( CYCC >= 2 )
       use cycc_component_mod,   only: cycc_component   ! TEMPORARY port (terrestrial carbon), to unblock OCYCC non-reg
#endif

       implicit none

       private

       public :: registry_build, registry_components, registry_count, component_ptr_t

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Compile-time count of active optional components. Each active flag adds 1.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer(ip), parameter :: N_COMPONENTS = 0                                                                              &
#if ( FROG_EXP > 0 )
     &                          + 1                                                                                            &
#endif
#if ( MEDUSA == 1 )
     &                          + 1                                                                                            &
#endif
#if ( ISM == 2 || ISM == 3 )
     &                          + 1                                                                                            &
#endif
#if ( CARAIB > 0 )
     &                          + 1                                                                                            &
#endif
#if ( OCYCC == 1 )
     &                          + 1                                                                                            &
#endif
#if ( CYCC >= 2 )
     &                          + 1                                                                                            &
#endif
#if ( CLM_INDICES >= 1 )
     &                          + 1                                                                                            &
#endif
     &                          + 0

! dmr --- A polymorphic pointer wrapper holds heterogeneous concrete components without dynamic allocation.
       type :: component_ptr_t
         class(coupled_component_t), pointer :: p => null()
       end type component_ptr_t

       type(component_ptr_t), save :: registry_components(max(N_COMPONENTS,1))
       integer(ip),          save :: registry_count = 0

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Associate the registry slots with the active component instances. Called once, before init_all.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine registry_build()

! dmr   NOTE on ordering: registration order drives BOTH step dispatch order AND finalize order. For CLM_INDICES these two
! dmr   wants conflict: its step must precede FROG/OCYCC (faithful to ec_ecbilt), but its finalize was LAST in emic.f (L992).
! dmr   Decision: prioritise step order (register FIRST); finalize order is assumed irrelevant (GLOBAL_FINALIZE just closes
! dmr   files, disjoint from other finalizes). Revisit if a finalize-order dependency surfaces.

         registry_count = 0

#if ( CLM_INDICES >= 1 )
         call add_component(clm_indices_component)   ! FIRST: step (DAILYSTEP) precedes FROG/OCYCC, faithful to ec_ecbilt
#endif
#if ( FROG_EXP > 0 )
         call add_component(frog_component)
#endif
#if ( MEDUSA == 1 )
!        call add_component(medusa_component)
#endif
#if ( ISM == 2 || ISM == 3 )
!        call add_component(grisli_component)
#endif
#if ( CARAIB > 0 )
!        call add_component(caraib_component)
#endif
#if ( OCYCC == 1 )
         call add_component(ocycc_component)
#endif
#if ( CYCC >= 2 )
         call add_component(cycc_component)   ! after OCYCC: ECO2(0) follows ocycc_ini, faithful to emic.f L250<L476
#endif

       end subroutine registry_build

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Helper: append one concrete component to the registry array. `comp` must have the `target` attribute at its module.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine add_component(comp)
         class(coupled_component_t), target, intent(inout) :: comp
         registry_count = registry_count + 1
         registry_components(registry_count)%p => comp
         registry_components(registry_count)%p%active = .true.
       end subroutine add_component

      end module component_registry_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
