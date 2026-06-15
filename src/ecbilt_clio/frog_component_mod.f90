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

#if ( FROG_EXP > 0 )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [frog_component_mod]   *** DRAFT - signatures to be confirmed against the real FROG source ***
!
!>     @brief Concrete coupled component wrapping FROG behind the generic contract. PILOT for the restructuring.
!
!      DESCRIPTION:
!>     Collects the three FROG blocks that are scattered across emic.f today into ONE place, each as a method:
!!       - init     : INITIALIZE_FROG, INIT_CPL2FROG, GET_COUPLING_STEP, INITIALIZE_FROGVARS
!!       - step     : DAILY_UPDATE_FROGVARS / RESET_FROGVARS_TIMER (sub-daily accumulation), then STEPFWD_FROG, GET_FROGVARS,
!!                    SET_FROG_FEEDBACK, FEEDBACK_FROG at the coupling cadence
!!       - finalize : WRITE_FROGRESTART
!!     The orchestrator never sees any of these calls; it only calls init/step/finalize/wants_phase.
!
!>     @warning Every USE-imported name, argument list and return type below is a BEST-GUESS reconstructed from the call
!!              surface seen in emic.f. They WILL need to be checked against the actual FROG module interfaces (subroutine
!!              vs function, argument order/types, which module exports what). Marked with TODO(frog) where most uncertain.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module frog_component_mod

       use global_constants_mod,  only: ip, dblp=>dp
       use coupled_component_mod, only: coupled_component_t, PHASE_ATM_STEP, PHASE_YEAR_END
       use coupling_timer_mod,    only: coupling_timer_t

       implicit none

       private

       public :: frog_component_t, frog_component

       type, extends(coupled_component_t) :: frog_component_t
         type(coupling_timer_t) :: timer
         integer(ip)            :: coupling_step = 0   !< value returned by GET_COUPLING_STEP at init
       contains
         procedure :: init        => frog_init
         procedure :: step        => frog_step
         procedure :: finalize    => frog_finalize
         procedure :: wants_phase => frog_wants_phase
       end type frog_component_t

! dmr --- module-level target instance, referenced by the registry
       type(frog_component_t), save, target :: frog_component

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   init: open FROG, set up the coupler->FROG bridge, read the coupling cadence, zero the accumulation buffers.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine frog_init(self, day, month)
         use comemic_mod, only: ntotday
         ! TODO(frog): confirm the owning module name(s) and exact interfaces of the routines below.
         ! use frog_mod, only: INITIALIZE_FROG, INIT_CPL2FROG, GET_COUPLING_STEP, INITIALIZE_FROGVARS
         class(frog_component_t), intent(inout) :: self
         integer(ip),             intent(in)    :: day, month

         self%name = "FROG"

!        call INITIALIZE_FROG()
!        call INIT_CPL2FROG()
!        self%coupling_step = GET_COUPLING_STEP()      ! TODO(frog): function or subroutine-with-out-arg?
!        call INITIALIZE_FROGVARS()

         self%timer%every_n_days = max(self%coupling_step, 1_ip)   ! TODO(frog): is coupling_step expressed in days?
         self%timer%ntotday      = ntotday
         self%timer%force_last   = .true.

       end subroutine frog_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   wants_phase: FROG accumulates every atmospheric step, and integrates at its coupling cadence (year-end here).
! dmr   We return .true. for BOTH phases and let frog_step branch on `phase`. (Daily accumulation vs cadence integration.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical function frog_wants_phase(self, phase, iday) result(wants)
         class(frog_component_t), intent(in) :: self
         integer(ip),             intent(in) :: phase, iday
         select case (phase)
           case (PHASE_ATM_STEP)
             wants = .true.                              ! sub-daily accumulation always runs
           case (PHASE_YEAR_END)
             wants = self%timer%is_due(iday)             ! integrate at cadence (or on the last day)
           case default
             wants = .false.
         end select
       end function frog_wants_phase

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   step: branch on phase. Accumulate on PHASE_ATM_STEP; integrate + feed back on the coupling cadence.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine frog_step(self, phase, iday, jstep, kstep)
         ! TODO(frog): confirm owning modules / interfaces.
         ! use frog_mod, only: DAILY_UPDATE_FROGVARS, RESET_FROGVARS_TIMER, STEPFWD_FROG, GET_FROGVARS,
         !                     SET_FROG_FEEDBACK, FEEDBACK_FROG
         class(frog_component_t), intent(inout) :: self
         integer(ip),             intent(in)    :: phase, iday, jstep, kstep

         select case (phase)

           case (PHASE_ATM_STEP)
!            call DAILY_UPDATE_FROGVARS()                ! accumulate atmospheric forcing into FROG buffers

           case (PHASE_YEAR_END)
!            call STEPFWD_FROG()                         ! advance FROG over the coupling interval
!            call GET_FROGVARS()                         ! pull FROG outputs
!            call SET_FROG_FEEDBACK()                    ! prepare feedback onto the climate state
!            call FEEDBACK_FROG()                        ! apply feedback
!            call RESET_FROGVARS_TIMER()                 ! zero the accumulation buffers for the next interval

         end select

       end subroutine frog_step

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   finalize: write the FROG restart.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine frog_finalize(self)
         ! use frog_mod, only: WRITE_FROGRESTART
         class(frog_component_t), intent(inout) :: self
!        call WRITE_FROGRESTART()
       end subroutine frog_finalize

      end module frog_component_mod

#endif /* FROG_EXP > 0 */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
