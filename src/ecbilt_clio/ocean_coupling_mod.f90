!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Copyright 2026 Didier M. Roche (a.k.a. dmr)

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

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#include "choixcomposantes.h"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      MODULE OCEAN_COUPLING_MOD

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  How the atmosphere gets its ocean / sea-ice boundary conditions (ROADMAP step 0a).
!
!  ocean_mode (namelist group oceanctl in the global "namelist" file):
!    "coupled" : ocn_bc gathered from CLIO every day (default, historical behaviour)
!    "record"  : as coupled, and ocn_bc is written to <ocean_bc_path>ocean_bc_<year>.nc
!    "replay"  : ocn_bc read from <ocean_bc_path>ocean_bc_<year>.nc; CLIO is initialised (restart read)
!                but never stepped: ec_co2oc, clio, CLIO restart files and CLIO outputs are skipped.
!
!  Acceptance: replay of a record run (same restart, same executable) gives bit-identical atmosphere and land.
!  For now record/replay is only allowed for the standard choixcomposantes.h (and included switch headers).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       USE global_constants_mod, ONLY: dblp=>dp, ip, str_len

       IMPLICIT NONE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   History
! dmr           0.1.0: created (coupled | record | replay)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      CHARACTER(LEN=5), PARAMETER :: version_mod ="0.1.0"

      PRIVATE
      PUBLIC :: ocean_coupling_read_nml, ocean_bc_update, ocean_is_replay

      CHARACTER(LEN=8)      , SAVE :: ocean_mode    = "coupled"
      CHARACTER(LEN=str_len), SAVE :: ocean_bc_path = "outputdata/coupler/"

      LOGICAL, SAVE :: is_record = .false., is_replay = .false.

      REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: ocn_bc_init   ! record mode: init gather, checked against day 1

      NAMELIST /oceanctl/ ocean_mode, ocean_bc_path

      CONTAINS

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ocean_coupling_read_nml(nml_unit)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Read the optional group oceanctl from the (already open, positioned) global namelist file.
!  Absent group => coupled.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        use, intrinsic :: iso_fortran_env, only: iostat_end

        INTEGER, INTENT(IN) :: nml_unit
        INTEGER             :: ios
        LOGICAL             :: std_config

        read(nml_unit, NML=oceanctl, iostat=ios)
        if (ios == iostat_end) then
           ocean_mode = "coupled"
        else if (ios /= 0) then
           write(*,*) "ocean_coupling: error reading namelist group oceanctl, iostat = ", ios
           STOP 1
        endif

        select case (trim(ocean_mode))
        case ("coupled")
        case ("record")
           is_record = .true.
        case ("replay")
           is_replay = .true.
        case default
           write(*,*) "ocean_coupling: unknown ocean_mode = ", trim(ocean_mode), " (coupled | record | replay)"
           STOP 1
        end select

        std_config = .true.
#if ( VEGGIE != 0 || FAST_OUTPUT != 0 || CLM_INDICES != 0 || BIOM_GEN != 0 || FROG_EXP != 0 || CARAIB != 0 || \
      CARAIB_FORC_W != 0 || HOURLY_RAD != 0 || OXYISO != 0 || D17ISO != 0 || WAXISO != 0 || VEG_LUH != 0 || \
      IMSK != 1 || COMATM != 1 || ROUTEAU != 1 || EVAPTRS != 1 || EVAPSI != 1 || CLAQUIN != 0 || F_PALAEO != 0 || \
      ICEBERG != 0 || ISOATM != 0 || WISOATM != 0 || WISOATM_RESTART != 0 || FRAC_KINETIK != 0 || ISOLBM != 0 || \
      WISOLND != 0 || WISOLND_RESTART != 0 || ISOOCN != 0 || WISOOCN != 0 || ISOBERG != 0 || ISM != 0 || \
      SMB_TYP != 0 || CPLTYP != 0 || REFREEZING != 0 || DOWNSTS != 0 || DOWNSCALING != 0 || SHELFMELT != 0 || \
      CALVFLUX != 0 || HEATFWF != 0 || CONSEAU != 0 || HEATCALV != 0 || CYCC != 0 || OCYCC != 0 || \
      INTERACT_CYCC != 0 || OLDC14 != 0 || KC14 != 0 || KC14P != 0 || OXNITREUX != 0 || WINDINCC != 0 || \
      WINDS_ERA5 != 0 || O2ATM != 0 || N2OATM != 0 || FROG_CARBON != 0 || MEDUSA != 0 || sediment_loopback != 0 || \
      PATH != 0 || NEOD != 0 || BRINES != 0 || ARGON != 0 || OOISO != 0 || RAYLEIGH != 0 || IRON_LIMITATION != 0 || \
      CORAL != 0 || COASTAL != 0 || CEMIS != 0 || REMIN != 0 || REMIN_CACO3 != 0 || ARAG != 0 || SILICA != 0 || \
      BATHY != 0 || GEOTHERMAL != 0 || ABEL != 0 || PROGRESS != 0 || UNCORFLUX != 0 || APPLY_UNCORFWF != 0 || \
      F_PALAEO_FWF != 0 || LGMSWITCH != 0 || FRAZER_ARCTIC != 0 || WRAP_EVOL != 0 || CLIO_OUT_NEWGEN != 0 || \
      CTRL_FIRST_ITER != 0 || L_TEST != 3 || NIT_RAP != 0 || TIDEMIX != 0 || XSLOP != 1 || I_COUPL != 1 || \
      forced_winds != 0 || NC_BERG != 0 || NC_IMSK != 0 || CFC != 0 || NUDGING != 0 || PERTATMOS != 0 || \
      PERTOCEAN != 0 || LONG_SED_RUN != 0 || DOWN_T2M != 0 )
        std_config = .false.
#endif
        if ((is_record .or. is_replay) .and. .not. std_config) then
           write(*,*) "ocean_coupling: ocean_mode = ", trim(ocean_mode), " is only allowed for now with the standard"
           write(*,*) "                choixcomposantes.h / clio_switches.h / BC_switches.h / additional_flags.h"
           STOP 1
        endif

        if (is_record) then
           call execute_command_line("mkdir -p "//trim(ocean_bc_path))
        endif

        write(*,*) "ocean_coupling: ocean_mode = ", trim(ocean_mode), "   ocean_bc_path = ", trim(ocean_bc_path)

      END SUBROUTINE ocean_coupling_read_nml

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ocean_bc_update(ist)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Set the atmospheric ocean / sea-ice surface state for day ist of this run (ist = 0: initialisation call).
!  Replaces the direct calls to ec_oc2co in the coupler initialisation and in the daily loop.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        use comemic_mod,     only: irunlabel
        use ocean_bc_mod,    only: ocn_bc, ocean_bc_apply, ocean_bc_write, ocean_bc_read
        use OCEAN2COUPL_COM, only: ec_oc2co_gather

        INTEGER, INTENT(IN) :: ist
        INTEGER(ip)         :: kday

        ! state at the start of day max(ist,1): the init call and day 1 see the same (restart) ocean state
        kday = int(irunlabel,ip)*360_ip + int(max(ist-1,0),ip)

        if (is_replay) then
           call ocean_bc_read(kday, ocean_bc_path)
        else
           call ec_oc2co_gather()
           if (is_record) then
              if (ist == 0) then
                 ocn_bc_init = ocn_bc
              else if (ist == 1) then
                 if (any(ocn_bc /= ocn_bc_init)) then
                    write(*,*) "ocean_coupling: ocean state changed between initialisation and day 1;"
                    write(*,*) "                replay cannot reproduce this run. Stopping."
                    STOP 1
                 endif
                 deallocate(ocn_bc_init)
              endif
              call ocean_bc_write(kday, ocean_bc_path)
           endif
        endif

        call ocean_bc_apply()

      END SUBROUTINE ocean_bc_update

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      LOGICAL FUNCTION ocean_is_replay()
        ocean_is_replay = is_replay
      END FUNCTION ocean_is_replay

      END MODULE OCEAN_COUPLING_MOD
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! The End of All Things (op. cit.)
