!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2026 Didier M. Roche (a.k.a. dmr) | iLOVECLIM / FRATRES coding group

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!   Style sheet: v1.0.0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [ocean_coupling_mod]
!
!>     @brief Selects how the atmosphere gets its ocean / sea-ice boundary conditions: coupled, record, replay or
!!            climatology.
!
!      DESCRIPTION:
!>     ocean_mode is read from the optional group oceanctl of the global "namelist" file (absent group = coupled):
!!       "coupled" : ocn_bndcon gathered from CLIO every day (historical behaviour).
!!       "record"  : as coupled, and ocn_bndcon is written to <ocean_bndcon_path>ocean_bndcon_<year>.nc.
!!       "replay"  : ocn_bndcon read from <ocean_bndcon_path>ocean_bndcon_<year>.nc.
!!       "climatology" : ocn_bndcon read from the single 360-record file ocean_bndcon_file, cycled by day of year.
!!     In replay and climatology the ocean is prescribed: CLIO is initialised (its restart is read) but never stepped,
!!     and ec_co2oc, clio, the CLIO restart files and CLIO outputs are skipped.
!!     Replay of a record run (same restart, same executable) gives bit-identical atmosphere, land and coupler restarts
!!     (ROADMAP step 0a). Climatology is ROADMAP step 0b. For now all modes other than coupled are only allowed with
!!     the standard flag headers.
!!
!!     Contract (public entry points):
!!       ocean_coupling_read_nml(nml_unit) : read oceanctl from the open global namelist, check the configuration.
!!       ocean_bndcon_update(ist)          : set the ocean / sea-ice surface state for day ist (0 = initialisation).
!!       ocean_is_prescribed()             : .true. in replay and climatology (callers skip the ocean model and its files).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module ocean_coupling_mod

        use global_constants_mod, only: str_len, ip, dblp=>dp

        implicit none

        private

        public :: ocean_coupling_read_nml, ocean_bndcon_update, ocean_is_prescribed

! dmr&clo   Values of ocean_mode.
        character(len=*), parameter :: OCEAN_MODE_COUPLED = "coupled"   ! ocean model stepped, nothing saved
        character(len=*), parameter :: OCEAN_MODE_RECORD  = "record"    ! ocean model stepped, ocn_bndcon saved
        character(len=*), parameter :: OCEAN_MODE_REPLAY  = "replay"        ! ocean model frozen, recorded ocn_bndcon read
        character(len=*), parameter :: OCEAN_MODE_CLIM    = "climatology"   ! ocean model frozen, 360-day ocn_bndcon cycled

! dmr&clo   Namelist oceanctl and derived state.
        character(len=16)      :: ocean_mode        = OCEAN_MODE_COUPLED      ! coupled | record | replay | climatology
        character(len=str_len) :: ocean_bndcon_path = "outputdata/coupler/"   ! directory of the yearly ocean_bndcon files
        character(len=str_len) :: ocean_bndcon_file = ""                      ! climatology: the 360-record file
        logical                :: is_record         = .false.
        logical                :: is_replay         = .false.
        logical                :: is_clim           = .false.

        real(dblp), dimension(:,:,:), allocatable :: ocn_bndcon_init   ! record mode: init gather, checked against day 1

        namelist /oceanctl/ ocean_mode, ocean_bndcon_path, ocean_bndcon_file

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Namelist reading and configuration checks.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine ocean_coupling_read_nml(nml_unit)

          use, intrinsic :: iso_fortran_env, only: iostat_end

          integer(ip), intent(in) :: nml_unit   !< global namelist file, open and positioned after tstepctl

          integer(ip) :: ios
          logical     :: std_config

          read(nml_unit, nml=oceanctl, iostat=ios)
          if (ios == iostat_end) then
            ocean_mode = OCEAN_MODE_COUPLED
          else if (ios /= 0) then
            write(*,*) "ocean_coupling: error reading namelist group oceanctl, iostat = ", ios
            stop 1
          endif

          select case (trim(ocean_mode))
          case (OCEAN_MODE_COUPLED)
          case (OCEAN_MODE_RECORD)
            is_record = .true.
          case (OCEAN_MODE_REPLAY)
            is_replay = .true.
          case (OCEAN_MODE_CLIM)
            is_clim = .true.
            if (len_trim(ocean_bndcon_file) == 0) then
              write(*,*) "ocean_coupling: ocean_mode = climatology needs ocean_bndcon_file in namelist group oceanctl"
              stop 1
            endif
          case default
            write(*,*) "ocean_coupling: unknown ocean_mode = ", trim(ocean_mode), " (coupled | record | replay | climatology)"
            stop 1
          end select

          ! record / replay are restricted to the standard choixcomposantes.h, clio_switches.h, BC_switches.h and
          ! additional_flags.h: any flag differing from its standard value aborts (strict on purpose, for now)
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

          if ((is_record .or. is_replay .or. is_clim) .and. .not. std_config) then
            write(*,*) "ocean_coupling: ocean_mode = ", trim(ocean_mode), " is only allowed for now with the standard"
            write(*,*) "                choixcomposantes.h / clio_switches.h / BC_switches.h / additional_flags.h"
            stop 1
          endif

          if (is_record) then
            call execute_command_line("mkdir -p "//trim(ocean_bndcon_path))
          endif

          if (is_clim) then
            write(*,*) "ocean_coupling: ocean_mode = ", trim(ocean_mode), "   ocean_bndcon_file = ", trim(ocean_bndcon_file)
          else
            write(*,*) "ocean_coupling: ocean_mode = ", trim(ocean_mode), "   ocean_bndcon_path = ", trim(ocean_bndcon_path)
          endif

        end subroutine ocean_coupling_read_nml

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Daily ocean / sea-ice surface state for the atmosphere (replaces the direct calls to ec_oc2co).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine ocean_bndcon_update(ist)

          use comemic_mod,      only: irunlabel
          use ocean_bndcon_mod, only: ocn_bndcon, ocean_bndcon_apply, ocean_bndcon_write, ocean_bndcon_read,          &
                                      ocean_bndcon_read_clim
          use OCEAN2COUPL_COM,  only: ec_oc2co_gather

          integer(ip), intent(in) :: ist   !< day of this run [1..ntotday], 0 for the initialisation call

          integer(ip) :: kday

          ! state at the start of day max(ist,1): the initialisation call and day 1 see the same (restart) ocean state
          kday = int(irunlabel,ip)*360_ip + max(ist-1_ip, 0_ip)

          if (is_replay) then
            call ocean_bndcon_read(kday, ocean_bndcon_path)
          else if (is_clim) then
            call ocean_bndcon_read_clim(kday, ocean_bndcon_file)
          else
            call ec_oc2co_gather()
            if (is_record) then
              if (ist == 0) then
                ocn_bndcon_init = ocn_bndcon
              else if (ist == 1) then
                if (any(ocn_bndcon /= ocn_bndcon_init)) then
                  write(*,*) "ocean_coupling: ocean state changed between initialisation and day 1;"
                  write(*,*) "                replay cannot reproduce this run. Stopping."
                  stop 1
                endif
                deallocate(ocn_bndcon_init)
              endif
              call ocean_bndcon_write(kday, ocean_bndcon_path)
            endif
          endif

          call ocean_bndcon_apply()

        end subroutine ocean_bndcon_update

        function ocean_is_prescribed() result(prescribed)
          logical :: prescribed

          prescribed = is_replay .or. is_clim
        end function ocean_is_prescribed

      end module ocean_coupling_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
