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
!      MODULE: [coupler_core_mod]
!
!>     @brief Generic orchestrator of the iLOVECLIM coupled run. Knows the coupling PHASES, not the components.
!
!      DESCRIPTION:
!>     coupler_core_mod owns the run lifecycle (setup / init_all / run_day / finalize_all) and drives the generic coupling
!!     phases described in the cartography document. The physics-core exchanges (ECBilt / CLIO / LBM) stay as today's
!!     ec_* / clio calls -- decision actee: do NOT abstract field exchange yet. Optional components are reached only through
!!     the abstract contract, via the static registry. The aim is that this module never grows a component-specific #if.
!
!>     @note Legacy core routines (ec_co2at, ec_at2co, ec_co2la, ec_sumfluxland, ec_sumfluxocean, clio, veget) are NOT in a
!!           module: they are external procedures. They MUST be declared with an explicit external interface below, otherwise
!!           gfortran infers a phantom interface from the call site that then clashes with the module's own contained
!!           procedures (error: "has an explicit interface from a previous declaration"). This is the fix for the build error.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module coupler_core_mod

       use global_constants_mod,   only: ip
       use coupled_component_mod,  only: coupled_component_t                                                                   &
     &                                , PHASE_DAY_BEGIN, PHASE_ATM_STEP, PHASE_LAND_STEP, PHASE_DAY_END, PHASE_YEAR_END
       use component_registry_mod, only: registry_build, registry_components, registry_count

       implicit none

       private

       public :: coupler_t

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Explicit interfaces for the legacy external (non-module) core routines. This is what fixes the compilation error.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       interface
         subroutine ec_co2at()
         end subroutine ec_co2at
         subroutine ec_at2co()
         end subroutine ec_at2co
         subroutine ec_co2la()
         end subroutine ec_co2la
         subroutine ec_sumfluxland(k)
           use global_constants_mod, only: ip
           integer(ip), intent(in) :: k
         end subroutine ec_sumfluxland
         subroutine ec_sumfluxocean(i, j)
           use global_constants_mod, only: ip
           integer(ip), intent(in) :: i, j
         end subroutine ec_sumfluxocean
         subroutine clio(i, iyearlabel, ntotday)
           use global_constants_mod, only: ip
           integer(ip), intent(in) :: i, iyearlabel, ntotday
         end subroutine clio
         subroutine veget(i, j, dtime, epss, pco2veg, frac, darea, tsurf)
           use global_constants_mod, only: ip, dblp=>dp
           integer(ip),     intent(in) :: i, j
           real(kind=dblp), intent(in) :: dtime, epss, pco2veg
           real(kind=dblp), intent(in) :: frac(*), darea(*), tsurf(*)
         end subroutine veget
       end interface

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   The orchestrator type.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       type :: coupler_t
         integer(ip) :: ntotday = 0
         integer(ip) :: iatm    = 1   ! atmospheric steps per day
         integer(ip) :: ilan    = 1   ! land steps per day
       contains
         procedure :: setup        => coupler_setup
         procedure :: init_all     => coupler_init_all
         procedure :: run_day      => coupler_run_day
         procedure :: finalize_all => coupler_finalize_all
         procedure :: run_ocean_to_coupler
         procedure :: run_atmos_step
         procedure :: run_land_step
         procedure :: accumulate_ocean_fluxes
         procedure :: run_coupler_to_ocean
         procedure :: run_vegetation
         procedure :: run_diagnostics
         procedure :: end_of_day
         procedure :: dispatch_components
       end type coupler_t

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   setup: coupler-level init (ex-ec_initemic) + build the static registry. NO component physics here.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine coupler_setup(self)
         use comemic_mod,          only: iatm, ilan, ntotday
         use initcoupledmodel_mod, only: ec_initemic
         class(coupler_t), intent(inout) :: self
         logical :: ok

         ok = ec_initemic()              ! reads emic.param/namelist/fractoc/darea, sets time counters (unchanged)
         self%ntotday = ntotday
         self%iatm    = iatm
         self%ilan    = ilan

         call registry_build()           ! associate the active optional components (the only component-aware step)

       end subroutine coupler_setup

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   init_all: core component init (ECBilt/CLIO/LBM via existing routine) then loop over registered optional components.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine coupler_init_all(self, day, month)
         use initcoupledmodel_mod, only: init_coupled_components
         class(coupler_t), intent(inout) :: self
         integer(ip),      intent(in)    :: day, month
         integer(ip) :: c
         logical :: ok

         ok = init_coupled_components(day, month)   ! ECBilt + CLIO + LBM + ec_initcoup (+ OCYCC/ICEBERG today; to be moved)

         do c = 1, registry_count
           call registry_components(c)%p%init(day, month)
         end do

       end subroutine coupler_init_all

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   run_day: one full day of the integration. Mirrors emic.f's i-loop body; calls the generic phases in order.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine coupler_run_day(self, iday)
         class(coupler_t), intent(inout) :: self
         integer(ip),      intent(in)    :: iday
         integer(ip) :: j, k

         call self%run_ocean_to_coupler(iday)                 ! ec_oc2co(i)             [emic.f L531]
         call self%dispatch_components(PHASE_DAY_BEGIN, iday, 0_ip, 0_ip)

         do j = 1, self%iatm
           call self%run_atmos_step(iday, j)                  ! atmosphere             [emic.f L543-555]
           do k = 1, self%ilan
             call self%run_land_step(iday, j, k)              ! land                   [emic.f L560-570]
           end do
           call self%accumulate_ocean_fluxes(iday, j)         ! ec_sumfluxocean(i,j)   [emic.f L574, EVERY j]

           if (j == self%iatm) then                           ! end of the day         [emic.f L586]
             call self%dispatch_components(PHASE_DAY_END, iday, j, 0_ip)   ! FROG daily, OCYCC step (pre-co2oc)
             call self%run_coupler_to_ocean(iday)             ! ec_co2oc/clio/carbon   [emic.f L601-618]
           end if

           call self%run_vegetation(iday, j)                  ! veget under flgveg     [emic.f L627-632]
           call self%run_diagnostics(iday, j)                 ! IPCC_output (+ downscaling) [emic.f L636+]
           call self%dispatch_components(PHASE_ATM_STEP, iday, j, 0_ip)   ! per-j component hooks (downscaling acc.)
         end do

         call self%dispatch_components(PHASE_YEAR_END, iday, 0_ip, 0_ip)   ! CARAIB/FROG/isotopes self-gate on yearly timer

         call self%end_of_day(iday)                           ! palaeo/bathy/daily_io + ec_writestate [emic.f L902-944]

       end subroutine coupler_run_day

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   finalize_all: optional components first (write restart / close), then core closures handled by caller.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine coupler_finalize_all(self)
         class(coupler_t), intent(inout) :: self
         integer(ip) :: c
         do c = 1, registry_count
           call registry_components(c)%p%finalize()
         end do
       end subroutine coupler_finalize_all

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   dispatch_components: ask each registered component whether it acts at this phase/day, then step it. Agnostic.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine dispatch_components(self, phase, iday, jstep, kstep)
         class(coupler_t), intent(inout) :: self
         integer(ip),      intent(in)    :: phase, iday, jstep, kstep
         integer(ip) :: c
         do c = 1, registry_count
           if (registry_components(c)%p%wants_phase(phase, iday)) then
             call registry_components(c)%p%import_state(iday)
             call registry_components(c)%p%step(phase, iday, jstep, kstep)
             call registry_components(c)%p%export_state(iday)
           end if
         end do
       end subroutine dispatch_components

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   The generic phases. Bodies hold the EXISTING core ec_* / clio calls (unchanged), to be filled during migration.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine run_ocean_to_coupler(self, iday)
         use ocean2coupl_com, only: ec_oc2co
         class(coupler_t), intent(inout) :: self
         integer(ip),      intent(in)    :: iday
         call ec_oc2co(iday)
       end subroutine run_ocean_to_coupler

       subroutine run_atmos_step(self, iday, j)
         use ecbilt0_mod, only: ec_update, ec_ecbilt
         use atmphys_mod, only: ec_fluxes
         use comsurf_mod, only: noc, nse
         class(coupler_t), intent(inout) :: self
         integer(ip),      intent(in)    :: iday, j
         call ec_update(iday, j)
         call ec_co2at
         call ec_at2co
         call ec_fluxes(noc)
         call ec_fluxes(nse)
         call ec_ecbilt(iday, j)
       end subroutine run_atmos_step

       subroutine run_land_step(self, iday, j, k)
         use landmodel_mod, only: ec_la2co, ec_co2la, ec_lbm, ec_lae2co
         use atmphys_mod,   only: ec_fluxes
         use comsurf_mod,   only: nld
         class(coupler_t), intent(inout) :: self
         integer(ip),      intent(in)    :: iday, j, k
         call ec_la2co
         call ec_fluxes(nld)
         call ec_co2la
         call ec_lbm(iday, j, k)
         call ec_lae2co
         call ec_sumfluxland(k)
       end subroutine run_land_step

       subroutine run_coupler_to_ocean(self, iday)
         use coupl2ocean_com,       only: ec_co2oc
         use atmos_composition_mod, only: atmos_carbon_update
         use comemic_mod,           only: irunlabel, iyear
         class(coupler_t), intent(inout) :: self
         integer(ip),      intent(in)    :: iday
         logical :: ok
         call ec_co2oc(iday)
         call clio(iday, irunlabel + iyear, self%ntotday)
         ok = atmos_carbon_update()
       end subroutine run_coupler_to_ocean

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   accumulate_ocean_fluxes: emic.f L574 -- called EVERY atmospheric step j (after the land loop, BEFORE the j==iatm test).
! dmr   It must run iatm times per day, not once. Misplacing it inside the j==iatm block broke non-regression.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine accumulate_ocean_fluxes(self, iday, j)
         class(coupler_t), intent(inout) :: self
         integer(ip),      intent(in)    :: iday, j
         call ec_sumfluxocean(iday, j)
       end subroutine accumulate_ocean_fluxes

       subroutine run_vegetation(self, iday, j)
         use comemic_mod,           only: flgveg
         use comsurf_mod,           only: nld, epss, fractn, tempsgn
         use comatm,                only: darea, dtime
         use atmos_composition_mod, only: get_PCO2_VEG
         class(coupler_t), intent(inout) :: self
         integer(ip),      intent(in)    :: iday, j
         if (flgveg) then
           call veget(iday, j, dtime, epss, get_PCO2_VEG(), fractn(1,1,nld), darea, tempsgn(1,1,nld))
         end if
       end subroutine run_vegetation

       subroutine run_diagnostics(self, iday, j)
         use ipcc_output_mod, only: ipcc_output
         class(coupler_t), intent(inout) :: self
         integer(ip),      intent(in)    :: iday, j
         logical :: ok
         ok = ipcc_output(iday, j)
         ! downscaling / CLIO_OUT_NEWGEN per-j accumulation handled via dispatch_components(PHASE_ATM_STEP)
       end subroutine run_diagnostics

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   end_of_day: emic.f L902-944 -- once per day, AFTER the j loop. palaeo/bathy timers, daily_io_nc, and crucially
! dmr   ec_writestate(i,ntotday) (restart write). Component-specific palaeo/bathy/IO go through dispatch (PHASE_YEAR_END or
! dmr   a dedicated service later); the restart write is core and stays here.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine end_of_day(self, iday)
         use emic_write_state_mod, only: ec_writestate
         class(coupler_t), intent(inout) :: self
         integer(ip),      intent(in)    :: iday
         logical :: ok
         ok = ec_writestate(iday, self%ntotday)
       end subroutine end_of_day

      end module coupler_core_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
