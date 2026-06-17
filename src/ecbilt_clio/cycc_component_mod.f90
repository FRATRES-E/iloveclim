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

#if ( CYCC >= 2 )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [cycc_component_mod]   *** TEMPORARY PORT -- faithful to emic.f, NOT the final design ***
!
!>     @brief Temporary adapter exposing the CYCC (terrestrial carbon cycle) calls behind the contract, to unblock the OCYCC
!!            non-regression. CYCC and OCYCC interact through the shared atmospheric CO2 state, so OCYCC cannot be validated
!!            with CYCC absent. This component reproduces the emic.f CYCC calls 1:1; a proper redesign is deferred.
!
!      DESCRIPTION:
!>     Ported emic.f calls (all under #if CYCC >= 2):
!!       init                  : ECO2(0,fractn(1,1,nld),darea)  -- emic.f L476 (after change of mask/topo; topo inactive in
!!                               the target run, so end-of-init placement is faithful here).
!!       step(PHASE_AFTER_OCEAN): ECO2(1,fractn(1,1,nld),darea) -- emic.f L722, ONCE PER DAY after the j loop (NOT annual,
!!                               NOT cadenced -- it sits under #if CYCC only, past enddo of the j loop at L655).
!!       finalize              : out_cycc(-2,fractn(1,1,nld),darea) -- emic.f L959 (closes C_reservoirs.txt).
!!
!>     @warning TEMPORARY. Known issues deferred to the real CYCC port:
!!       - out_cycc is poorly factored (mixes global and OCYCC-specific output); reintroduced here ONLY for bit-faithful
!!         restart/output, to be untangled later.
!!       - send_caraib2lbm (emic.f L672, guard CARAIB>1 & CYCC>1) is NOT here: it belongs to the future CARAIB component.
!!       - ECO2(0) ordering vs topo/mask services holds only because those #if blocks are inactive in this run.
!!       - ECO2 / out_cycc / C14ATM_DP are legacy external routines (no owning module); declared external below.
!!
!>     @note KC14 is a SUB-OPTION of CYCC (not a component): under #if KC14, C14ATM_DP is called once per year (KENDY==1),
!!           BEFORE ECO2(1), in step(PHASE_AFTER_OCEAN) -- emic.f L716-722. KENDY is set ONLY by sync_timer, which is called
!!           ONLY from OCYCC (ocycc_step). So KC14 implicitly depends on OCYCC being active and stepped BEFORE CYCC in the
!!           same day; the phase order (OCYCC in BEFORE_OCEAN, CYCC in AFTER_OCEAN) guarantees KENDY is fresh. This implicit
!!           dependency existed unguarded in emic.f and is preserved as-is.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module cycc_component_mod

       use global_constants_mod,  only: ip
       use coupled_component_mod, only: coupled_component_t, PHASE_AFTER_OCEAN

       implicit none

       private

       public :: cycc_component_t, cycc_component

       type, extends(coupled_component_t) :: cycc_component_t
       contains
         procedure :: init        => cycc_init
         procedure :: step        => cycc_step
         procedure :: finalize    => cycc_finalize
         procedure :: wants_phase => cycc_wants_phase
       end type cycc_component_t

! dmr --- module-level target instance, referenced by the static registry
       type(cycc_component_t), save, target :: cycc_component

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   init: ECO2(0) -- initialise atmospheric CO2 in coupled-carbon mode. emic.f L476.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine cycc_init(self, day, month)
         use comsurf_mod,  only: nld, fractn
         use comatm,       only: darea
         use carbone_co2,  only: PA0_C, PA_C
         class(cycc_component_t), intent(inout) :: self
         integer(ip),             intent(in)    :: day, month
         external :: ECO2

         self%name = "CYCC"
         call ECO2(0_ip, fractn(1,1,nld), darea)
         write(*,*) "cycc_component init : PA0_C, PA_C", PA0_C, PA_C   ! faithful to emic.f L477 debug write

       end subroutine cycc_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   wants_phase: CYCC acts once per day, after the ocean (post j-loop in emic.f). No sub-cadence.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical function cycc_wants_phase(self, phase, iday) result(wants)
         class(cycc_component_t), intent(in) :: self
         integer(ip),             intent(in) :: phase, iday
         wants = (phase == PHASE_AFTER_OCEAN)
       end function cycc_wants_phase

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   step: ECO2(1) -- daily update of atmospheric CO2. emic.f L722 (after the j loop, every day).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine cycc_step(self, phase, iday, jstep, kstep)
         use comsurf_mod, only: nld, fractn
         use comatm,      only: darea
#if ( KC14 == 1 )
         use mod_sync_time, only: KENDY      ! end-of-year flag, set by sync_timer (via OCYCC's ocycc_step)
#endif
         class(cycc_component_t), intent(inout) :: self
         integer(ip),             intent(in)    :: phase, iday, jstep, kstep
         external :: ECO2
#if ( KC14 == 1 )
         external :: C14ATM_DP
#endif

#if ( KC14 == 1 )
         if (KENDY == 1) then                ! once per year, emic.f L716-718. Runs BEFORE ECO2(1) (emic.f order).
           call C14ATM_DP
         end if
#endif

         call ECO2(1_ip, fractn(1,1,nld), darea)

       end subroutine cycc_step

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   finalize: out_cycc -- close C_reservoirs.txt. emic.f L959. (Reintroduced for faithful output despite the known
! dmr   factoring debt; will move when out_cycc is untangled.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine cycc_finalize(self)
         use comsurf_mod, only: nld, fractn
         use comatm,      only: darea
         class(cycc_component_t), intent(inout) :: self
         external :: out_cycc

         call out_cycc(-2_ip, fractn(1,1,nld), darea)

       end subroutine cycc_finalize

      end module cycc_component_mod

#endif /* CYCC >= 2 */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
