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
!      MODULE: [frog_component_mod]
!
!>     @brief Orchestration adapter exposing FROG (permafrost) as a coupled component. PILOT for the restructuring.
!
!      DESCRIPTION:
!>     This module is the "WHEN" layer for FROG. It does NOT translate fields (that is CPL2FROG_mod's job, the "WHAT" layer)
!!     and it does NOT implement permafrost physics (that is main_lib_FROG, the library). It only groups, behind the generic
!!     coupled_component_t contract, the FROG calls that emic.f scatters across init / daily loop / yearly block / shutdown:
!!
!!       init                 : INITIALIZE_FROG (lib), INIT_CPL2FROG(GET_COUPLING_STEP()) (cpl+lib), then the priming
!!                              sequence DAILY_UPDATE_FROGVARS / INITIALIZE_FROGVARS(GET_FROGVARS()) / RESET_FROGVARS_TIMER
!!                              -- reproduces emic.f L494-506 exactly.
!!       step(PHASE_BEFORE_OCEAN): DAILY_UPDATE_FROGVARS() -- once per day, before ec_co2oc/clio (emic.f L589).
!!       step(PHASE_YEAR_END) : STEPFWD_FROG(GET_FROGVARS()) / SET_FROG_FEEDBACK(FEEDBACK_FROG()) / RESET_FROGVARS_TIMER
!!                              -- once per model year, emic.f L701-706 (after the j loop).
!!       finalize             : WRITE_FROGRESTART(-1) -- emic.f L967.
!!
!>     @note Cadence: integration is ANNUAL (mod(i, days_year360d_i)==0). GET_COUPLING_STEP() returns the number of days in
!!           the coupling window; it only sizes the accumulation buffers at init (it is NOT the integration period here).
!>     @note The init-time DAILY_UPDATE_FROGVARS (L502) primes storing_time_step=1 before the first INITIALIZE_FROGVARS; it is
!!           intentionally kept and is distinct from the per-day accumulation in step().
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module frog_component_mod

       use global_constants_mod,  only: ip, days_year360d_i
       use coupled_component_mod, only: coupled_component_t, PHASE_BEFORE_OCEAN, PHASE_YEAR_END
       use coupling_timer_mod,    only: coupling_timer_t

       implicit none

       private

       public :: frog_component_t, frog_component

       type, extends(coupled_component_t) :: frog_component_t
         type(coupling_timer_t) :: timer            !< annual cadence (every days_year360d_i days)
       contains
         procedure :: init        => frog_init
         procedure :: step        => frog_step
         procedure :: finalize    => frog_finalize
         procedure :: wants_phase => frog_wants_phase
       end type frog_component_t

! dmr --- module-level target instance, referenced by the static registry
       type(frog_component_t), save, target :: frog_component

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   init: emic.f L494-506. Library init, then coupler-bridge init sized by GET_COUPLING_STEP(), then prime the buffers.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine frog_init(self, day, month)
         use comemic_mod,    only: ntotday
         use infodisplay_mod, only: write_im
         use main_lib_FROG,  only: INITIALIZE_FROG, GET_COUPLING_STEP, INITIALIZE_FROGVARS
         use CPL2FROG_mod,   only: INIT_CPL2FROG, GET_FROGVARS, DAILY_UPDATE_FROGVARS, RESET_FROGVARS_TIMER
         class(frog_component_t), intent(inout) :: self
         integer(ip),             intent(in)    :: day, month
         logical :: well_done

         self%name = "FROG"

         well_done = INITIALIZE_FROG()
         if (well_done) call write_im("FROG INITIALIZATION COMPLETE", "FROG COMPONENT")

         well_done = INIT_CPL2FROG(GET_COUPLING_STEP())

         call DAILY_UPDATE_FROGVARS()                         ! prime storing_time_step = 1 (emic.f L502)
         well_done = INITIALIZE_FROGVARS(GET_FROGVARS())      ! seed the library buffers (emic.f L504)
         call RESET_FROGVARS_TIMER()                          ! emic.f L506

         self%timer%every_n_days = days_year360d_i            ! annual integration cadence
         self%timer%ntotday      = ntotday
         self%timer%force_last   = .false.                    ! emic.f does NOT force a final-day FROG step; stay faithful

       end subroutine frog_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   wants_phase: act every day (accumulation) and at the annual cadence (integration).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical function frog_wants_phase(self, phase, iday) result(wants)
         class(frog_component_t), intent(in) :: self
         integer(ip),             intent(in) :: phase, iday
         select case (phase)
           case (PHASE_BEFORE_OCEAN)
             wants = .true.                                   ! DAILY_UPDATE_FROGVARS every day (emic.f L589)
           case (PHASE_YEAR_END)
             wants = self%timer%is_due(iday)                  ! STEPFWD only when mod(i, days_year360d_i)==0
           case default
             wants = .false.
         end select
       end function frog_wants_phase

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   step: branch on phase. Daily accumulation vs annual integration+feedback.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine frog_step(self, phase, iday, jstep, kstep)
         use main_lib_FROG, only: STEPFWD_FROG, FEEDBACK_FROG
         use CPL2FROG_mod,  only: DAILY_UPDATE_FROGVARS, RESET_FROGVARS_TIMER, GET_FROGVARS, SET_FROG_FEEDBACK
         class(frog_component_t), intent(inout) :: self
         integer(ip),             intent(in)    :: phase, iday, jstep, kstep
         logical :: well_done

         select case (phase)

           case (PHASE_BEFORE_OCEAN)                          ! emic.f L589
             call DAILY_UPDATE_FROGVARS()

           case (PHASE_YEAR_END)                              ! emic.f L701-706
             well_done = STEPFWD_FROG(GET_FROGVARS())
             well_done = SET_FROG_FEEDBACK(FEEDBACK_FROG())
             call RESET_FROGVARS_TIMER()

         end select

       end subroutine frog_step

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   finalize: write the FROG restart (emic.f L967).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine frog_finalize(self)
         use main_lib_FROG, only: WRITE_FROGRESTART
         class(frog_component_t), intent(inout) :: self
         logical :: well_done
         well_done = WRITE_FROGRESTART(-1_ip)
       end subroutine frog_finalize

      end module frog_component_mod

#endif /* FROG_EXP > 0 */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
