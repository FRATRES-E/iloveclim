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
! dmr   MAIN PROGRAM FOR THE iLOVECLIM COUPLED EARTH SYSTEM MODEL  (restructured)
!
!>      The whole run now reads as a lifecycle. The orchestrator (coupler_core_mod) owns the coupling phases; the static
!!      registry (component_registry_mod) owns which optional components are compiled in. This file contains NO component
!!      #if branches and no physics: it is the table of contents of an iLOVECLIM run. Formerly "program emic".
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      program main

       use global_constants_mod, only: ip
       use comemic_mod,          only: iday, imonth
       use coupler_core_mod,     only: coupler_t
       use infodisplay_mod,      only: write_im
       use error0_mod,           only: ec_error

       implicit none

       type(coupler_t) :: coupler
       integer(ip)     :: i

       call write_im("Execution started", "ILOVECLIM MAIN")

       call coupler%setup()                 ! ex-ec_initemic + registry_build
       call coupler%init_all(iday, imonth)  ! core + optional components init

       call write_im("Finalized init", "ILOVECLIM MAIN")

       do i = 1, coupler%ntotday
         call coupler%run_day(i)            ! the generic phases + component dispatch
       end do

       call coupler%finalize_all()          ! optional components close / write restart

       call write_im("Execution ended", "ILOVECLIM MAIN")
       call ec_error(999)

      end program main

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
