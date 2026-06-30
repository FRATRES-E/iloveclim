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
!      MODULE: C_res_mod
!
!>     @author  N. Bouttes (nb), D. Roche (dmr), refactoring by dmr, clo
!
!>     @brief   Carbon-reservoir state for the iLOVECLIM carbon cycle.
!
!>     @details Pure state module (declarations only, no procedures).  Holds the initial,
!>              dynamic (running) and restart values of the carbon reservoirs for the ocean
!>              and the land, together with their 13C and 14C counterparts and the
!>              atmospheric 13C.  Replaces the historical include file C_RESERVOIRS.inc.
!>
!>              Refactoring notes (relative to the original C_res_mod.f90):
!>                - REAL*8 / REAL(KIND=8) replaced by real(dblp) from global_constants_mod
!>                  (dblp is 8 bytes, so this is binary-compatible with rest_cc.dat).
!>                - integer unit numbers and type(fileDescriptor) I/O handles kept AS-IS;
!>                  the I/O machinery is out of scope for this pass.
!>                - the unused CO2-cycle arrays NCO2, NCO2max, TPSCO2, PCO2M were removed
!>                  (no reference anywhere in the source tree; see MIGRATION.md).
!>                - French inline comments translated to English, intent preserved.
!>                - all preprocessing guards kept verbatim (KC14, CORAL, COASTAL).
!>
!>     @date    Original   : 04 December 2007 (nb); verified for C2.5 04 March 2008 (nb, dmr)
!>     @date    Refactored : 2026-06-30 (dmr, clo)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module C_res_mod

        use global_constants_mod, only: dblp=>dp, ip
        use file_libs,            only: fileDescriptor

        implicit none
        public

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Output-file unit numbers (legacy integer units; I/O machinery unchanged)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        !  Unit for C_reservoir.txt
        integer(ip) :: c_res_fich
#if ( KC14 == 1 )
        integer(ip) :: c14_res_fich, c_flux_fich, dc14_fich, c14_bilan_fich
#endif
#if ( CORAL == 1 )
        integer(ip) :: coral_res_fich, c_nino_fich
#endif
#if ( COASTAL == 1 )
        integer(ip) :: coastal_res_fich
#endif

        !  fileDescriptor-based handles (kept as-is; see header note)
        type(fileDescriptor), public, save :: c13_res_fich, c_ocean_fich,    &
                                              fPOC_fich, fCAL_fich

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Total carbon reservoirs (ocean / land): initial, running, restart      [GtC]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: ca_oc_ini, ca_la_ini
        real(dblp) :: ca_oc_rest, ca_la_rest
        real(dblp) :: alk_oc_rest
        real(dblp) :: cav_oc = 0.0_dblp, cav_la = 0.0_dblp, cav_oc2

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  13C reservoirs and atmospheric 13C
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: ca13_oc_ini, ca13_la_ini
        real(dblp) :: ca13_oc_rest, ca13_la_rest
        real(dblp) :: ca13_at_ini, dc13at_ini

        !  DOC bookkeeping and ocean wet volume
        real(dblp) :: ca_oc_vol, coc_odoc, coc_odocs, coc_odoc13, coc_odocs13
        real(dblp) :: cav_oc13, cav_la13

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Cumulative carbon emissions (anthropogenic and permafrost)              [GtC]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: emis_cum
        real(dblp) :: emis_c13_cum
        real(dblp) :: emis_perm_cum
        real(dblp) :: emis_perm_c13_cum

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Previous-year reservoirs (diagnostic deltas)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: cav_oc_p = 0.0_dblp
        real(dblp) :: cav_la_p = 0.0_dblp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Atmospheric 13C restart / init scalars
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: c13atm_rest
        real(dblp) :: c13atm, c13init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  14C reservoirs and fluxes (ocean / land)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( KC14 == 1 )
        real(dblp) :: ca14_oc_ini, ca14_la_ini
        real(dblp) :: cav_oc14_b, cav_oc14 = 0.0_dblp, cav_oc_b
        real(dblp) :: cav_la14_b, cav_la14 = 0.0_dblp, cav_la_b
        real(dblp) :: cav_oc14_rest, cav_la14_rest
        real(dblp) :: cav_la14_b_rest, cav_oc14_b_rest
        real(dblp) :: cav_oc_b_rest, cav_la_b_rest

        !  Air-sea / air-land carbon fluxes (12C and 14C)
        real(dblp) :: FC12OA, FC12LA, FC14OA, FC14LA
#endif

      end module C_res_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
