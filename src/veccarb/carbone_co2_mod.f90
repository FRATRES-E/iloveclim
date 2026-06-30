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
!      MODULE: carbone_co2_mod
!
!>     @author  Didier M. Roche (dmr), Nathaelle Bouttes (nb), Veronique Mariotti (vm)
!>     @author  refactoring by dmr, clo
!
!>     @brief   Atmospheric carbon state for the iLOVECLIM / VECODE carbon cycle.
!
!>     @details This module declares the atmospheric reservoirs and run-control scalars
!>              of the carbon cycle: the total atmospheric carbon (PA_C and friends),
!>              the 14C atmospheric activity (C14ATM) with its radioactive-decay and
!>              standard-ratio constants, and the optional production-file, emission,
!>              oxygen and permafrost fields gated by preprocessing flags.
!>
!>              Refactoring notes (relative to the original carbone_co2.f90):
!>                - kinds now come from global_constants_mod (dblp, ip) instead of the
!>                  default real/integer; literals carry an explicit _dblp suffix.
!>                - the COMATM==0 branch (#include "veget.h" / "comsurf.h") has been
!>                  removed: COMATM==1 is now assumed unconditionally, so nlat/nlon are
!>                  always taken from comatm.
!>                - the historical name was the bare module `carbone_co2`; it is renamed
!>                  carbone_co2_mod for consistency with the *_mod convention.  Callers
!>                  must update `use carbone_co2` -> `use carbone_co2_mod` (see MIGRATION.md).
!>                - isotope grid state that used to live in veget_iso (b*t1[34] …) is NOT
!>                  here: it is part of the per-cell / grid state in veget_mod, the single
!>                  source of truth after the VECODE refactoring.
!>
!>     @date    Original   : 14 December 2009 (dmr, nb)
!>     @date    Last mod.  : 02 December 2026 (EA, production file in ec14)
!>     @date    Refactored : 2026-06-29 (dmr, clo)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module carbone_co2_mod

        use global_constants_mod, only: dblp=>dp, ip

        ! COMATM == 1 is now assumed unconditionally (see header note).
        use comatm,               only: nlat, nlon

        implicit none
        public

        private :: nlat, nlon   ! grid dimensions visible internally only

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Run-control flag
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        !  Fresh-start switch for the carbon reservoirs:
        !    1 = no carbon restart available, initialise pools from equilibrium
        !    0 = isotope / reservoir restart present, read from rest_cc.dat
        integer(ip), parameter :: new_run_c = 1

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Atmospheric carbon — total CO2 reservoir
!    PA_C   : current atmospheric carbon                          [GtC]
!    PA0_C  : reference (initial / restart) atmospheric carbon    [GtC]
!    PA_C_D : diagnostic atmospheric carbon (alt. ocean term)     [GtC]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: PA_C, PA0_C, PA_C_D

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Atmospheric carbon-14
!    C14DEC  : 14C radioactive decay constant (1/8240 yr^-1)      [yr^-1]
!    c14rstd : standard modern 14C activity ratio                 [–]
!    C14ATM* : atmospheric 14C activity (current / initial / rest)[–]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp), parameter :: C14DEC  = 1.0_dblp / 8240.0_dblp
        real(dblp), parameter :: c14rstd = 1.176e-12_dblp
        real(dblp) :: C14ATM, C14ATM0, C14ATM_rest

#if ( KC14 == 1 )
        !  14C mass-budget bookkeeping (cosmogenic production path)
        real(dblp) :: PRODC14      ! cosmogenic 14C production        [ppm 14C yr^-1]
        real(dblp) :: MPRODC14     ! same, expressed as mass          [gC14 yr^-1]
        real(dblp) :: LIMC14MASSE  ! steady-state 14C mass limit      [gC14]
        real(dblp) :: N14LIM       ! steady-state 14C mole limit      [mol C14]
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  14C cosmogenic production read from file (KC14P == 1)
!  Modified 02/12/2026 (EA): production file consumed in ec14 / c14_prod_mod.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( KC14P == 1 )
        integer(ip), parameter :: NC14max = 2
        real(dblp), dimension(NC14max) :: PC14M     ! 14C production samples     [–]
        real(dblp), dimension(NC14max) :: TPSC14    ! sample time axis (yr, <0 = BP)
        integer(ip) :: KTIME       ! 1 = real time, else fixed
        integer(ip) :: n           ! loop / record index
        integer(ip) :: NC14        ! number of records actually read
        integer(ip) :: NYR         ! current model year
        integer(ip) :: NYR0        ! resolved calendar year
        integer(ip) :: NYR01       ! -NYR0
        integer(ip) :: NYRSR       ! start-of-run reference year
        real(dblp)  :: PC14VAR     ! interpolated 14C production         [–]
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Atmospheric oxygen reservoir (O2ATM == 1)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( O2ATM == 1 )
        real(dblp) :: PA_O, PA0_O  ! current / reference atmospheric O2  [arbitrary]
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Anthropogenic carbon emission scenario (CEMIS == 1)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( CEMIS == 1 )
        integer(ip), parameter :: nb_emis = 450    ! number of data lines in emission file
        real(dblp), dimension(nb_emis) :: cemis    ! annual carbon emission           [GtC yr^-1]
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Permafrost carbon emission scenario (PERM_SCEN == 1)
!  NOTE (refactoring): the original carbone_co2.f90 did NOT declare these, yet eco2.f
!  imported them from `carbone_co2` under PERM_SCEN.  The original therefore could not
!  build with PERM_SCEN==1 against the archived module.  Declared here, mirroring CEMIS,
!  so the path is consistent.  Adjust nb_emis_perm to the real file length if needed.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( PERM_SCEN == 1 )
        integer(ip), parameter :: nb_emis_perm = 450    ! number of data lines in permafrost file
        real(dblp), dimension(nb_emis_perm) :: cemis_perm  ! annual permafrost emission   [GtC yr^-1]
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Permafrost carbon coupling (FROG_EXP > 0 && FROG_CARBON > 0)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( FROG_EXP > 0 && FROG_CARBON > 0 )
        real(dblp) :: deepC        ! soil/permafrost carbon replacing b4 [GtC]
#endif

      end module carbone_co2_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
