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
!      MODULE: [coupling_timer_mod]
!
!>     @brief Small reusable helper encapsulating the cadence at which a component couples.
!
!      DESCRIPTION:
!>     Today the coupling cadence is hard-wired in the main loop as mod(i, ndays_medusastep), mod(i, timCplGRISday), or a
!!     yearly mod(i, days_year360d_i). This type moves that decision into a component property so the orchestrator can ask
!!     "is it due today?" without knowing the component. It also flags the last day of the run (several components force a
!!     final call on i==ntotday).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module coupling_timer_mod

       use global_constants_mod, only: ip

       implicit none

       private

       public :: coupling_timer_t

       type :: coupling_timer_t
         integer(ip) :: every_n_days = 1        !< coupling period in days (1 = every day, 360 = yearly)
         integer(ip) :: ntotday      = huge(1)  !< total days in the run, set at init; enables the last-day override
         logical     :: force_last   = .true.   !< if .true., is_due() also returns .true. on the final day
       contains
         procedure :: is_due => timer_is_due
       end type coupling_timer_t

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Returns .true. when the component should be stepped on day iday.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical function timer_is_due(self, iday) result(due)
         class(coupling_timer_t), intent(in) :: self
         integer(ip),             intent(in) :: iday

         due = ( mod(iday, self%every_n_days) == 0 )
         if (self%force_last .and. (iday == self%ntotday)) due = .true.

       end function timer_is_due

      end module coupling_timer_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
