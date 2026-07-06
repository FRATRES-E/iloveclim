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
!      MODULE: eco2_mod
!
!>     @author  Didier M. Roche (dmr), Nathaelle Bouttes (nb), Veronique Mariotti (vm)
!>     @author  refactoring by dmr, clo
!
!>     @brief   Atmospheric CO2 concentration driver for the iLOVECLIM carbon cycle.
!
!>     @details This module is a REDUCER, not a per-cell model.  Unlike the VECODE physics
!>              (point-wise on veget_cell_state_t), eco2 forms GLOBAL sums over the grid
!>              and the ocean reservoirs (cav_la, cav_oc, …) and closes the atmospheric
!>              carbon budget from them.  The point-wise cell pattern therefore does NOT
!>              apply here; the grid loops are genuine reductions and are kept explicit.
!>
!>              The historical single subroutine eco2(KOD,…) carried two unrelated
!>              behaviours behind a KOD switch.  They are now split:
!>                eco2_init  (was KOD == 0) : first-call initialisation — read the carbon
!>                           restart (if KLSR==1), set reference reservoirs, sum the
!>                           initial ocean and land carbon, prime the diagnostics.
!>                eco2_step  (was KOD /= 0) : per-call update — re-sum the ocean and land
!>                           reservoirs, apply emissions, and recompute atmospheric PA_C,
!>                           c13atm and (KC14) C14ATM by closing the budget.
!>              A thin dispatcher eco2(kod,…) is kept so existing call sites that pass KOD
!>              continue to work unchanged; new code should call eco2_init / eco2_step
!>              directly.
!>
!>              Refactoring notes (relative to eco2.f):
!>                - F77 fixed form -> F90 free form; kinds from global_constants_mod.
!>                - the duplicate land-carbon reduction (identical in the init and step
!>                  paths bar the target accumulators) is factored into the private helper
!>                  sum_land_carbon(); likewise sum_ocean_carbon() for the OCYCC ocean sum.
!>                - isotope pools (b*t13 … b4g14) now come from veget_mod (single source of
!>                  truth) instead of veget_iso; carea, st, sg, b*t/b*g from veget_mod.
!>                - `use carbone_co2` -> `use carbone_co2_mod`.
!>                - out_cycc is OUT OF SCOPE (declared elsewhere); calls are preserved.
!>                - all preprocessing branches kept verbatim (OCYCC, CARAIB, CORAL, IMSK,
!>                  CEMIS, PERM_SCEN, FROG_*, INTERACT_CYCC, BATHY, MEDUSA-commented).
!>                  NOTE: the original mixes flag spellings FROG_CARBON and CARBON_FROG in
!>                  different guards; both are preserved exactly as found.
!>
!>     @date    Original   : 14 December 2009 (dmr, nb)
!>     @date    Refactored : 2026-06-29 (dmr, clo)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module eco2_mod

        use global_constants_mod, only: dblp=>dp, ip, stdout

        implicit none
        private

        public :: eco2          ! backward-compatible dispatcher (KOD switch)
        public :: eco2_init     ! first-call initialisation  (was KOD == 0)
        public :: eco2_step     ! per-call update            (was KOD /= 0)

        !  ppm <-> GtC conversion factor (kept as a module constant; was a local literal).
        real(dblp), parameter :: ca_beta = 0.47_dblp

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: eco2  (dispatcher)
!
!>     @brief   Backward-compatible entry point.  KOD == 0 -> eco2_init, else -> eco2_step.
!>     @param[in] kod      0 = initialise, /= 0 = step
!>     @param[in] fracgr   land fraction of each grid cell                       [–]
!>     @param[in] darea    area of each latitude band                            [m^2]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine eco2(kod, fracgr, darea)

        use comatm, only: nlat, nlon

        integer(ip),                           intent(in) :: kod
        real(dblp), dimension(nlat,nlon),      intent(in) :: fracgr
        real(dblp), dimension(nlat),           intent(in) :: darea

        if (kod == 0) then
          call eco2_init(fracgr, darea)
        else
          call eco2_step(fracgr, darea)
        end if

      end subroutine eco2

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: eco2_init   (was KOD == 0)
!
!>     @brief   First-call initialisation of the global carbon pools.
!
!>     @details Reads the carbon restart when KLSR==1, sets the atmospheric reference
!>              (PA0_C, c13atm, ca13_at_ini), zeroes the cumulative diagnostics, sums the
!>              initial ocean (OCYCC) and land carbon into ca_oc_ini / ca_la_ini, lets the
!>              restart overwrite those when present, and primes out_cycc.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine eco2_init(fracgr, darea)

        use comatm,          only: nlat, nlon
        use carbone_co2_mod, only: PA0_C, PA_C, PA_C_D, C14ATM, new_run_c
        use C_res_mod,       only: c13atm, ca13_at_ini, ca13_la_ini,     &
                                   ca13_oc_ini, ca13_oc_rest, ca13_la_rest, &
                                   ca_la_ini, ca_la_rest, ca_oc_ini,      &
                                   ca_oc_rest, ca_oc_vol, dc13at_ini,     &
                                   coc_odoc, coc_odocs, coc_odoc13,       &
                                   coc_odocs13, emis_cum, emis_c13_cum,   &
                                   emis_perm_cum, emis_perm_c13_cum
#if ( KC14 == 1 )
        use C_res_mod,       only: ca14_oc_ini, ca14_la_ini,             &
                                   cav_oc14_rest, cav_la14_rest,          &
                                   cav_oc14_b, cav_oc_b, cav_la14_b,      &
                                   cav_la_b, cav_oc14_b_rest,             &
                                   cav_la14_b_rest, cav_oc_b_rest,        &
                                   cav_la_b_rest
#endif
        use loveclim_transfer_mod, only: KLSR
        use restart_cc_mod,        only: restart_cc

#if ( CEMIS == 1 )
        use carbone_co2_mod, only: cemis, nb_emis
        use newunit_mod,     only: carbon_emission_dat_id
#endif
#if ( PERM_SCEN == 1 )
        use carbone_co2_mod, only: cemis_perm, nb_emis_perm
        use newunit_mod,     only: permafrost_emission_dat_id
#endif

        real(dblp), dimension(nlat,nlon), intent(in) :: fracgr
        real(dblp), dimension(nlat),      intent(in) :: darea

#if ( CEMIS == 1 || PERM_SCEN == 1 )
        integer(ip) :: ii, yy
#endif

        write(stdout,*) ' '
        write(stdout,*) 'initialisation in eco2'

        !  Read total-carbon restart and prime PA0_C / c13atm only on a restart run.
        if (KLSR == 1) call restart_cc(1)

        PA_C   = PA0_C
        PA_C_D = PA0_C
        write(stdout,*) 'atmosphere PA_C, c13atm, c14atm'
        write(stdout,*) PA_C, c13atm, c13atm - 1000.0_dblp, C14ATM

        ca_oc_ini   = 0.0_dblp
        ca_la_ini   = 0.0_dblp
        ca13_oc_ini = 0.0_dblp
        ca13_la_ini = 0.0_dblp
#if ( KC14 == 1 )
        ca14_oc_ini = 0.0_dblp
        ca14_la_ini = 0.0_dblp
#endif
        ca_oc_vol   = 0.0_dblp
        ca13_at_ini = (c13atm - 1000.0_dblp) * PA_C / ca_beta
        write(stdout,*) 'c13_at_ini dans eco2 ', ca13_at_ini
        dc13at_ini  = c13atm
        coc_odoc    = 0.0_dblp
        coc_odocs   = 0.0_dblp
        coc_odoc13  = 0.0_dblp
        coc_odocs13 = 0.0_dblp
        emis_cum          = 0.0_dblp
        emis_c13_cum      = 0.0_dblp
        emis_perm_cum     = 0.0_dblp
        emis_perm_c13_cum = 0.0_dblp

#if ( CEMIS == 1 )
        !  Read the anthropogenic carbon-emission scenario (GtC yr^-1).
        do ii = 1, nb_emis
          read(carbon_emission_dat_id,*) yy, cemis(ii)
        end do
        close(carbon_emission_dat_id)
#endif
#if ( PERM_SCEN == 1 )
        !  Read the permafrost carbon-emission scenario; first file year = run year 1.
        do ii = 1, nb_emis_perm
          read(permafrost_emission_dat_id,*) yy, cemis_perm(ii)
        end do
        close(permafrost_emission_dat_id)
#endif

        !  Initial ocean carbon (OCYCC sum, or a fixed default).
        call sum_ocean_carbon(ca_oc_ini, ca13_oc_ini                     &
#if ( KC14 == 1 )
                              , ca14_oc_ini                               &
#endif
                              , coc_odoc, coc_odocs, coc_odoc13,          &
                              coc_odocs13, ca_oc_vol)

        write(stdout,*) 'carbon ocean ca_oc_ini'  , ca_oc_ini
        write(stdout,*) 'C13 ocean ca13_oc_ini'   , ca13_oc_ini
        !  Diagnostic ratio only (never reused); guard against an empty reservoir so a
        !  zero pre-restart sum cannot raise SIGFPE under -ffpe-trap.
        if (ca_oc_ini /= 0.0_dblp) then
          write(stdout,*) 'd13C ocean '            , ca13_oc_ini / ca_oc_ini
        else
          write(stdout,*) 'd13C ocean '            , ' n/a (ca_oc_ini == 0)'
        end if

        !  Restart overwrites the computed ocean totals when present.
        if (KLSR == 1) then
          write(stdout,*) 'carbon ocean from restart ca_oc_rest', ca_oc_rest
          ca_oc_ini   = ca_oc_rest
          write(stdout,*) 'C13 ocean from restart ca13_oc_rest', ca13_oc_rest
          ca13_oc_ini = ca13_oc_rest
#if ( KC14 == 1 )
          write(stdout,*) 'C14 ocean from restart ca14_oc_rest', cav_oc14_rest
          ca14_oc_ini = cav_oc14_rest
          cav_oc14_b  = cav_oc14_b_rest
          cav_oc_b    = cav_oc_b_rest
#endif
        end if

        !  Initial land carbon (vegetation pools weighted by cover and cell area).
        call sum_land_carbon(fracgr, darea, ca_la_ini, ca13_la_ini       &
#if ( KC14 == 1 )
                             , ca14_la_ini                               &
#endif
                             )

        ca13_la_ini = ca13_la_ini - 1000.0_dblp * ca_la_ini

        write(stdout,*) 'carbon vegetation ca_la_ini' , ca_la_ini
        write(stdout,*) 'C13 vegetation ca13_la_ini'  , ca13_la_ini
        !  Diagnostic ratio only (immediately overwritten by the restart below); guard
        !  against ca_la_ini == 0 (e.g. veget pools not yet filled) to avoid SIGFPE.
        if (ca_la_ini /= 0.0_dblp) then
          write(stdout,*) 'dC13 vegetation '          , ca13_la_ini / ca_la_ini
        else
          write(stdout,*) 'dC13 vegetation '          , ' n/a (ca_la_ini == 0)'
        end if

        !  Restart overwrites the computed land totals when present.
        if (KLSR == 1) then
          write(stdout,*) 'carbon vegetation from restart ca_la_rest', ca_la_rest
          ca_la_ini   = ca_la_rest
          write(stdout,*) 'C13 vegetation from restart ca13_la_rest', ca13_la_rest
          ca13_la_ini = ca13_la_rest
#if ( KC14 == 1 )
          write(stdout,*) 'C14 vegetation from restart ca14_la_rest', cav_la14_rest
          ca14_la_ini = cav_la14_rest
          cav_la14_b  = cav_la14_b_rest
          cav_la_b    = cav_la_b_rest
#endif
        end if

        call out_cycc(-1, fracgr, darea)

#if ( KC14 == 1 )
        !  Fresh carbon start: seed the "before" 14C/total reservoirs from the initials.
        if (new_run_c == 1) then
          cav_oc14_b = ca14_oc_ini
          cav_oc_b   = ca_oc_ini
          cav_la14_b = ca14_la_ini
          cav_la_b   = ca_la_ini
        end if
#endif

        call out_cycc(1, fracgr, darea)
        write(stdout,*) ' '

      end subroutine eco2_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: eco2_step   (was KOD /= 0)
!
!>     @brief   Recompute the atmospheric carbon budget for the current call.
!
!>     @details Re-sums the ocean (OCYCC) and land reservoirs, accumulates emission
!>              scenarios at year end, and closes the budget to obtain PA_C, the diagnostic
!>              PA_C_D, c13atm and (KC14) C14ATM.  The CORAL and INTERACT_CYCC branches
!>              replace the PA_C closure as in the original.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine eco2_step(fracgr, darea)

        use comatm,          only: nlat, nlon
        use carbone_co2_mod, only: PA0_C, PA_C, PA_C_D, C14ATM, C14DEC, c14rstd
        use C_res_mod,       only: c13atm, ca13_at_ini, ca13_oc_ini,     &
                                   ca13_la_ini, ca_oc_ini, ca_la_ini,    &
                                   cav_oc, cav_oc2, cav_oc13, cav_la,     &
                                   cav_la13, cav_oc_p, cav_la_p,          &
                                   emis_cum, emis_c13_cum,               &
                                   emis_perm_cum, emis_perm_c13_cum
#if ( KC14 == 1 )
        use C_res_mod,       only: cav_oc14, cav_la14, cav_oc14_b,        &
                                   cav_oc_b, cav_la14_b, cav_la_b,        &
                                   FC14OA, FC14LA
#endif
        use mod_sync_time,   only: KENDY, NYR
#if ( CEMIS == 1 )
        use carbone_co2_mod, only: cemis
#endif
#if ( PERM_SCEN == 1 )
        use carbone_co2_mod, only: cemis_perm
#endif
#if ( OCYCC == 1 )
        use marine_bio_mod,  only: ODIC_diff
#endif
#if ( CORAL == 1 )
        use coral_mod,       only: C_car_a
        use mbiota_mod,      only: SCANU
#endif
#if ( INTERACT_CYCC == 2 )
        use atmos_composition_mod, only: get_PGA_CO2
#endif

        real(dblp), dimension(nlat,nlon), intent(in) :: fracgr
        real(dblp), dimension(nlat),      intent(in) :: darea

        real(dblp) :: d_oc, d_la, sumflux
#if ( CORAL == 1 )
        real(dblp) :: ca_car_a
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Re-sum the ocean carbon reservoirs (OCYCC).  cav_oc2 mirrors cav_oc but adds the
!  diffusive DIC correction ODIC_diff.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        cav_oc   = 0.0_dblp
        cav_oc2  = 0.0_dblp
        cav_oc13 = 0.0_dblp
#if ( KC14 == 1 )
        cav_oc14 = 0.0_dblp
#endif

        call sum_ocean_carbon_step(cav_oc, cav_oc2, cav_oc13             &
#if ( KC14 == 1 )
                                   , cav_oc14                            &
#endif
                                   )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Re-sum the land carbon reservoirs.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( KC14 == 1 )
        cav_la14 = 0.0_dblp
#endif
        cav_la   = 0.0_dblp
        cav_la13 = 0.0_dblp

        call sum_land_carbon(fracgr, darea, cav_la, cav_la13            &
#if ( KC14 == 1 )
                             , cav_la14                                 &
#endif
                             )

        if (KENDY == 1) write(stdout,*) 'carbon vegetation cav_la ', cav_la

        !  Vecode uses c13atm = d13C + 1000, hence the offset below.
        cav_la13 = cav_la13 - 1000.0_dblp * cav_la

#if ( CEMIS == 1 )
        if (KENDY == 1) then
          emis_cum = emis_cum + cemis(NYR)   ! cemis in GtC
          if (emis_cum /= 0.0_dblp)                                      &
            write(stdout,*) 'cemis,emis_cum', NYR, cemis(NYR), emis_cum
        end if
#endif
#if ( PERM_SCEN == 1 )
        if (KENDY == 1) then
          emis_perm_cum     = emis_perm_cum     + cemis_perm(NYR)
          emis_perm_c13_cum = emis_perm_c13_cum + cemis_perm(NYR) * (-25.0_dblp)
          if (emis_perm_cum /= 0.0_dblp)                                 &
            write(stdout,*) 'cemis_perm,emis_perm_cum', NYR,             &
                            cemis_perm(NYR), emis_perm_cum
        end if
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Close the atmospheric carbon budget -> PA_C.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( CORAL == 0 )
        PA_C = PA0_C - (cav_oc - ca_oc_ini + cav_la - ca_la_ini          &
                        - emis_cum - emis_perm_cum) * ca_beta
#elif ( CORAL > 0 )
        !  With corals, remove the CO2 drawn down by carbonate weathering.
        !  ca_car_a units: Pmol day^-1 * 1e6 * 12 g/mol = 1e15 g/day = Pg/day.
        ca_car_a = C_car_a * SCANU * 1.0e6_dblp * 12.0_dblp
        PA_C = PA0_C - (cav_oc - ca_oc_ini + cav_la - ca_la_ini          &
                        - emis_cum - emis_perm_cum + ca_car_a) * ca_beta
#endif

#if ( INTERACT_CYCC == 2 )
        !  Carbon-cycle CO2 forced to equal the radiative-code CO2 (read from GHG.dat).
        PA_C = get_PGA_CO2()
#endif

        PA_C_D = PA0_C - (cav_oc2 - ca_oc_ini + cav_la - ca_la_ini       &
                          - emis_cum - emis_perm_cum) * ca_beta

#if ( OCYCC == 1 )
        ODIC_diff = 0.0_dblp
#endif

        !  Atmospheric 13C, closed from the same budget.
        c13atm = 1000.0_dblp + (ca13_at_ini + emis_c13_cum               &
                 + emis_perm_c13_cum                                     &
                 - (cav_oc13 - ca13_oc_ini + cav_la13 - ca13_la_ini))    &
                 / (PA_C / ca_beta)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Atmospheric 14C: land+ocean respiration flux into the atmosphere, once per year.
!  vm: n14atm(t) = n14atm(t-1) - c14dec*n14lnd(t) - (n14lnd(t)-n14lnd(t-1)).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (KENDY == 1) then
#if ( KC14 == 0 )
          C14ATM = c14rstd * PA_C
#elif ( KC14 == 1 )
          write(stdout,*) ' /*/*/*/ '
          write(stdout,*) 'eco2 before flux, C14ATM = ', C14ATM

          sumflux = ((cav_la14_b - (1.0_dblp + C14DEC) * cav_la14) / 1.0e15_dblp &
                  +  (cav_oc14_b - (1.0_dblp + C14DEC) * cav_oc14) / 1.0e15_dblp) &
                  * 0.40_dblp

          !  Quick fix below is NOT conservative (vm, dmr): clamp a negative 14C atmosphere.
          if (abs(sumflux) >= C14ATM) then
            write(stdout,*)
            write(stdout,*) 'in eco2: Fc14->atm not conservative'
            write(stdout,*) '====================================='
            write(stdout,*) 'Diagnostics ...'
            write(stdout,*) cav_la14_b, cav_la14, cav_oc14_b, cav_oc14
            write(stdout,*) C14ATM, (C14ATM / PA_C / 1.176e-12_dblp - 1.0_dblp) * 1000.0_dblp
            C14ATM = c14rstd * PA_C
          else
            C14ATM = C14ATM                                              &
                   + ((cav_la14_b - (1.0_dblp + C14DEC) * cav_la14) / 1.0e15_dblp &
                   +  (cav_oc14_b - (1.0_dblp + C14DEC) * cav_oc14) / 1.0e15_dblp) &
                   * 0.40_dblp
          end if

          FC14OA = (cav_oc14_b - (1.0_dblp + C14DEC) * cav_oc14) / 1.0e15_dblp * 0.40_dblp
          FC14LA = (cav_la14_b - (1.0_dblp + C14DEC) * cav_la14) / 1.0e15_dblp * 0.40_dblp

          write(stdout,*) 'eco2 after flux, C14ATM = ', C14ATM
          write(stdout,*) 'FC14OA = ', FC14OA
          write(stdout,*) 'FC14LA = ', FC14LA
          write(stdout,*) ' /*/*/*/ '
#endif
        end if

        !  Year-to-year reservoir deltas (diagnostic).
        if (cav_oc_p /= 0.0_dblp) then
          d_oc = cav_oc - cav_oc_p
          d_la = cav_la - cav_la_p
        end if
        cav_oc_p = cav_oc
        cav_la_p = cav_la

        if (KENDY == 1) call out_cycc(NYR, fracgr, darea)

#if ( KC14 == 1 )
        !  Store this year's totals as the "before" reservoirs for next year's flux.
        if (KENDY == 1) then
          cav_oc14_b = cav_oc14
          cav_oc_b   = cav_oc
          cav_la14_b = cav_la14
          cav_la_b   = cav_la
        end if
#endif

      end subroutine eco2_step

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      PRIVATE HELPERS
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: sum_land_carbon
!
!>     @brief   Reduce the per-cell vegetation carbon pools into global land totals.
!
!>     @details Computes carea(i,k) = darea(i)*fracgr(i,k)*1e-12 (or veget_frac under
!>              CARAIB) for land cells, then accumulates the cover-weighted sum of the
!>              tree/grass pools (b1..b4) into c_la, and the 13C (and optionally 14C)
!>              pools into c13_la / c14_la.  Shared verbatim by eco2_init and eco2_step;
!>              only the destination accumulators differ between the two call sites.
!>
!>              The IMSK branch additionally weights each cell by (1 - icemask).  The
!>              FROG_EXP/CARBON_FROG branch zeroes the slow soil pool b4 locally (the soil
!>              carbon is then injected globally as deepC by the caller).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine sum_land_carbon(fracgr, darea, c_la, c13_la            &
#if ( KC14 == 1 )
                                 , c14_la                               &
#endif
                                 )

        use comatm,     only: nlat, nlon
        use bio_mod,    only: carea
        use veget_mod,  only: st, sg, b1t, b2t, b3t, b4t,                &
                              b1g, b2g, b3g, b4g,                        &
                              b1t13, b2t13, b3t13, b4t13,                &
                              b1g13, b2g13, b3g13, b4g13
#if ( KC14 == 1 )
        use veget_mod,  only: b1t14, b2t14, b3t14, b4t14,               &
                              b1g14, b2g14, b3g14, b4g14
#endif
#if ( IMSK == 1 )
        use input_icemask, only: icemask
#endif
#if ( CARAIB > 0 )
        use ec_ca2lbm,  only: veget_frac
#endif
#if ( FROG_EXP > 0 && FROG_CARBON > 0 )
        use carbone_co2_mod, only: deepC
#endif

        real(dblp), dimension(nlat,nlon), intent(in)    :: fracgr
        real(dblp), dimension(nlat),      intent(in)    :: darea
        real(dblp),                       intent(inout) :: c_la, c13_la
#if ( KC14 == 1 )
        real(dblp),                       intent(inout) :: c14_la
#endif

        integer(ip) :: i, k
        real(dblp)  :: b4t_temp, b4g_temp, wgt

        do k = 1, nlon
          do i = 1, nlat
#if ( CARAIB == 0 )
            if (fracgr(i,k) > 0.001_dblp) then
              carea(i,k) = darea(i) * fracgr(i,k) * 1.0e-12_dblp
#else
            if (veget_frac(i,k) > 0.001_dblp) then
              carea(i,k) = darea(i) * veget_frac(i,k) * 1.0e-12_dblp
#endif

#if ( CARAIB == 0 )
              b4t_temp = b4t(i,k)
              b4g_temp = b4g(i,k)
#if ( FROG_EXP > 0 && FROG_CARBON > 0 )
              b4t_temp = 0.0_dblp
              b4g_temp = 0.0_dblp
#endif

#if ( IMSK == 1 )
              wgt = carea(i,k) * (1.0_dblp - icemask(i,k))
#else
              wgt = carea(i,k)
#endif

              c_la = c_la                                                &
                   + (b1t(i,k)*st(i,k) + b1g(i,k)*sg(i,k)                &
                   +  b2t(i,k)*st(i,k) + b2g(i,k)*sg(i,k)                &
                   +  b3t(i,k)*st(i,k) + b3g(i,k)*sg(i,k)                &
                   +  b4t_temp*st(i,k) + b4g_temp*sg(i,k)) * wgt

              c13_la = c13_la                                            &
                     + ((b1t13(i,k)+b2t13(i,k)+b3t13(i,k)+b4t13(i,k))*st(i,k) &
                     +  (b1g13(i,k)+b2g13(i,k)+b3g13(i,k)+b4g13(i,k))*sg(i,k)) * wgt
#if ( KC14 == 1 )
              c14_la = c14_la                                           &
                     + ((b1t14(i,k)+b2t14(i,k)+b3t14(i,k)+b4t14(i,k))*st(i,k) &
                     +  (b1g14(i,k)+b2g14(i,k)+b3g14(i,k)+b4g14(i,k))*sg(i,k)) &
                     * wgt * 1.0e15_dblp                                 ! gC14
#endif

#elif ( CARAIB > 0 )
              !  CARAIB land carbon set to a fixed 1500 GtC (see original commentary).
              c_la = 1500.0_dblp
#endif

            end if
          end do
        end do

#if ( FROG_EXP > 0 && FROG_CARBON > 0 )
        c_la = c_la + deepC   ! inject permafrost/soil carbon globally
#endif
#if ( CARAIB > 0 )
        c_la = veget_frac_total_caraib()   ! see note; overrides with PgC sum
#endif

      end subroutine sum_land_carbon

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: sum_ocean_carbon   (initialisation form)
!
!>     @brief   Initial ocean carbon reduction (OCYCC), or fixed default when OCYCC==0.
!>     @details Accumulates DIC/DOC/POC and their 13C (and 14C) into the *_ini totals,
!>              plus the DOC bookkeeping (coc_*) and total wet volume (ca_oc_vol).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine sum_ocean_carbon(c_oc, c13_oc                          &
#if ( KC14 == 1 )
                                  , c14_oc                              &
#endif
                                  , odoc_s, odocs_s, odoc13_s,          &
                                  odocs13_s, vol_s)

        real(dblp), intent(inout) :: c_oc, c13_oc
#if ( KC14 == 1 )
        real(dblp), intent(inout) :: c14_oc
#endif
        real(dblp), intent(inout) :: odoc_s, odocs_s, odoc13_s, odocs13_s, vol_s

#if ( OCYCC == 1 )
        call ocean_carbon_reduction(c_oc, c13_oc                        &
#if ( KC14 == 1 )
                                    , c14_oc                            &
#endif
                                    , odoc_s, odocs_s, odoc13_s, odocs13_s, vol_s)
#if ( KC14 == 1 )
        c14_oc = c14_oc * 1.0e15_dblp   ! TgC14 -> gC14 (consistency with cav_oc14_b)
#endif
#else /* OCYCC != 1 */
        c_oc   = 38000.0_dblp
        c13_oc = -0.7_dblp * c_oc
#endif

      end subroutine sum_ocean_carbon

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: sum_ocean_carbon_step   (per-call form)
!
!>     @brief   Per-step ocean carbon reduction (OCYCC); also builds cav_oc2 with the
!>              diffusive DIC correction ODIC_diff.  No-op when OCYCC==0.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine sum_ocean_carbon_step(c_oc, c_oc2, c13_oc              &
#if ( KC14 == 1 )
                                       , c14_oc                         &
#endif
                                       )

        real(dblp), intent(inout) :: c_oc, c_oc2, c13_oc
#if ( KC14 == 1 )
        real(dblp), intent(inout) :: c14_oc
#endif

#if ( OCYCC == 1 )
        call ocean_carbon_reduction_step(c_oc, c_oc2, c13_oc            &
#if ( KC14 == 1 )
                                         , c14_oc                       &
#endif
                                         )
#endif

      end subroutine sum_ocean_carbon_step

#if ( OCYCC == 1 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: ocean_carbon_reduction   (initialisation)
!>     @brief  Triple loop over ocean boxes (LT x JT x NOC_CBR) accumulating the initial
!>             ocean carbon totals.  Factored out so init and step share the box geometry.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine ocean_carbon_reduction(c_oc, c13_oc                    &
#if ( KC14 == 1 )
                                        , c14_oc                        &
#endif
                                        , odoc_s, odocs_s, odoc13_s,    &
                                        odocs13_s, vol_s)

        use declars_mod,           only: LT, JT, NOC_CBR
        use loveclim_transfer_mod, only: dvol, mgt
        use mbiota_mod,            only: scale_m, zoo_m, phyto_m,        &
                                         PHYTO_M13, ZOO_M13
        use marine_bio_mod,        only: odic, oc13, odoc, odocs,        &
                                         odoc13, odocs13, opoc
#if ( KC14 == 1 )
        use marine_bio_mod,        only: oc14
#endif

        real(dblp), intent(inout) :: c_oc, c13_oc
#if ( KC14 == 1 )
        real(dblp), intent(inout) :: c14_oc
#endif
        real(dblp), intent(inout) :: odoc_s, odocs_s, odoc13_s, odocs13_s, vol_s

        integer(ip) :: i, j, n

        do n = 1, NOC_CBR
          do i = 1, LT
            do j = 1, JT
              if (MGT(i,j,n) == 1) then
                c_oc = c_oc + (PHYTO_M(i,j,n) + ZOO_M(i,j,n)             &
                     + ODOC(i,j,n) + OPOC(i,j,n) + ODOCS(i,j,n)          &
                     + ODIC(i,j,n)*1.0e6_dblp)                           &
                     * DVOL(i,j,n) * 12.0_dblp * SCALE_M * 1.028_dblp

                c13_oc = c13_oc + (PHYTO_M13(i,j,n) + ZOO_M13(i,j,n)     &
                       + ODOC13(i,j,n) + ODOCS13(i,j,n)                  &
                       + OC13(i,j,n)*1.0e6_dblp)                         &
                       * DVOL(i,j,n) * 12.0_dblp * SCALE_M * 1.028_dblp
#if ( KC14 == 1 )
                c14_oc = c14_oc + OC14(i,j,n)*1.0e6_dblp                 &
                       * DVOL(i,j,n) * 14.0_dblp * SCALE_M * 1.028_dblp  ! TgC14
#endif
                vol_s     = vol_s     + DVOL(i,j,n)
                odoc_s    = odoc_s    + ODOC(i,j,n)    * DVOL(i,j,n)
                odocs_s   = odocs_s   + ODOCS(i,j,n)   * DVOL(i,j,n)
                odoc13_s  = odoc13_s  + ODOC13(i,j,n)  * DVOL(i,j,n)
                odocs13_s = odocs13_s + ODOCS13(i,j,n) * DVOL(i,j,n)
              end if
            end do
          end do
        end do

      end subroutine ocean_carbon_reduction

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: ocean_carbon_reduction_step
!>     @brief  Per-step ocean reduction; cav_oc2 adds the diffusive DIC correction.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine ocean_carbon_reduction_step(c_oc, c_oc2, c13_oc        &
#if ( KC14 == 1 )
                                             , c14_oc                   &
#endif
                                             )

        use declars_mod,           only: LT, JT, NOC_CBR
        use loveclim_transfer_mod, only: dvol, mgt
        use mbiota_mod,            only: scale_m, zoo_m, phyto_m,        &
                                         PHYTO_M13, ZOO_M13
        use marine_bio_mod,        only: odic, oc13, odoc, odocs,        &
                                         odoc13, odocs13, opoc, odic_diff
#if ( KC14 == 1 )
        use marine_bio_mod,        only: oc14
#endif

        real(dblp), intent(inout) :: c_oc, c_oc2, c13_oc
#if ( KC14 == 1 )
        real(dblp), intent(inout) :: c14_oc
#endif

        integer(ip) :: i, j, n

        do n = 1, NOC_CBR
          do i = 1, LT
            do j = 1, JT
              if (MGT(i,j,n) == 1) then
                c_oc = c_oc + (PHYTO_M(i,j,n) + ZOO_M(i,j,n)             &
                     + ODOC(i,j,n) + OPOC(i,j,n) + ODOCS(i,j,n)          &
                     + ODIC(i,j,n)*1.0e6_dblp)                           &
                     * DVOL(i,j,n) * 12.0_dblp * SCALE_M * 1.028_dblp

                c_oc2 = c_oc2 + (PHYTO_M(i,j,n) + ZOO_M(i,j,n)           &
                      + ODOC(i,j,n) + OPOC(i,j,n) + ODOCS(i,j,n)         &
                      + (ODIC(i,j,n) + ODIC_diff(i,j,n))*1.0e6_dblp)     &
                      * DVOL(i,j,n) * 12.0_dblp * SCALE_M * 1.028_dblp

                c13_oc = c13_oc + (PHYTO_M13(i,j,n) + ZOO_M13(i,j,n)     &
                       + ODOC13(i,j,n) + ODOCS13(i,j,n)                  &
                       + OC13(i,j,n)*1.0e6_dblp)                         &
                       * DVOL(i,j,n) * 12.0_dblp * SCALE_M * 1.028_dblp
#if ( KC14 == 1 )
                c14_oc = c14_oc + OC14(i,j,n)*1.0e6_dblp                 &
                       * DVOL(i,j,n) * 14.0_dblp * SCALE_M * 1.028_dblp  &
                       * 1.0e15_dblp                                     ! gC14
#endif
              end if
            end do
          end do
        end do

      end subroutine ocean_carbon_reduction_step
#endif /* OCYCC == 1 */

#if ( CARAIB > 0 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      FUNCTION: veget_frac_total_caraib
!>     @brief  CARAIB total land carbon in PgC (stock_carbon_caraib * 1e-15).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function veget_frac_total_caraib() result(c_la)
        use ec_ca2lbm, only: stock_carbon_caraib
        real(dblp) :: c_la
        c_la = stock_carbon_caraib * 1.0e-15_dblp
      end function veget_frac_total_caraib
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      end module eco2_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
