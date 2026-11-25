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
!      MODULE: [iso_nml_mod]
!
!>     @brief Runtime settings of the atmospheric water isotopes: namelist group wisoatm_nml in the file wisoatm.nml
!!            (ROADMAP D6, step 2b).
!
!      DESCRIPTION:
!
!>     Reads, checks and reports the isotope settings, and turns the delta values given in the namelist into the
!!     isotopic ratios used by the model (rsmow = 1 for every species: a ratio is the isotopologue mass per water mass).
!!
!!     Delta conventions: inside the code deltas are plain fractions (order 1e-3, never multiplied by 1000); permil and
!!     per meg appear only in the namelist, converted once in iso_nml_read:
!!       d18o_*    delta 18O in permil (VSMOW)
!!       dexc_*    deuterium excess in permil:   dD = 8 d18O + dexc
!!       o17exc_*  17O excess in per meg:        ln(d17O + 1) = 0.528 ln(d18O + 1) + o17exc
!!     The 0.528 above is the reference slope of the 17O-excess definition, not a fractionation exponent.
!!
!!     Contract
!!       - iso_nml_read(adv_mode) once at initialisation, before the isotope fields are initialised; stops on a missing
!!         file, an unknown or not yet available setting, and unless adv_mode = 'fv' (O31).
!!       - iso_nml_write(unit) writes the settings to parameters.chk.
!!       - iso_ratios(d18o, dexc, o17exc) gives the ratios of all species for one set of deltas, as fractions.
!!       - After iso_nml_read: r_vap_init, r_land_evap, r_snow_evap and (ocean_ratio_mode = 'fixed') r_ocean_fixed hold
!!         the ratios of all species; iwater and iwat16 are 1.
!
!      REFERENCES: ROADMAP D6, D10, D30, D34, §3.3, §4; O31, O40.
!
!      AUTHORS: dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module iso_nml_mod

#if ( WISO == 1 )

        use global_constants_mod, only: str_len, ip, dblp=>dp, tK_zero_C
        use comatm,               only: nwisos, iwater, iwat16, iwat17, iwat18, iwat2h
        use iso_alphas,           only: ISO_FORM_LEN, iso_alphas_init, iso_alphas_report

        implicit none

        private

        public :: ISO_NML_FILE, iso_nml_read, iso_nml_write, iso_ratios
        public :: frac_mode, frac_evap, frac_lsc, frac_conv, iso_passive
        public :: ocean_ratio_mode, ocean_d18o_file, restart_read, restart_file, land_iso_mode
        public :: lbm_restart_read, lbm_restart_file, land_init_mode, land_init_file
        public :: d18o_soil_init, dexc_soil_init, o17exc_soil_init, d18o_snow_init, dexc_snow_init, o17exc_snow_init
        public :: r_vap_init, r_land_evap, r_snow_evap, r_ocean_fixed, dexc_ocean_frac, o17exc_ocean_frac
        public :: land_bucket, r_soil_init, r_snow_init
        public :: delta_check_min, delta_check_max
        public :: kin_wind, kin_k18_frac, kin_kd_ratio, cg_rv_gamma, bc_eq_lsc, bc_eq_conv
        public :: cond_temp, cond_phase, cond_tliq_k, cond_tice_k


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        character(len=*), parameter :: ISO_NML_FILE  = 'wisoatm.nml'   !< namelist file, in the run directory
        real(dblp),       parameter :: ISO_SLOPE_DEX = 8.0_dblp        !< dD = 8 d18O + d-excess
        real(dblp),       parameter :: ISO_SLOPE_17  = 0.528_dblp      !< reference slope of the 17O-excess definition
        real(dblp),       parameter :: ISO_PERMIL    = 1.0e-3_dblp     !< permil to fraction (namelist input only)
        real(dblp),       parameter :: ISO_PERMEG    = 1.0e-6_dblp     !< per meg to fraction (namelist input only)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Namelist variables (defaults)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr&clo   Fractionation (wiso:2d:0000): frac_mode = 0 none (passive tracers), 1 equilibrium, 2 equilibrium + kinetic,
!           for every process; frac_mode = -1: per process, frac_evap (ocean evaporation: 1 equilibrium with the ocean,
!           2 Craig-Gordon), frac_lsc (large-scale condensation) and frac_conv (convection) (1 equilibrium, 2 + kinetic
!           effect for snow). After iso_nml_read, frac_evap, frac_lsc and frac_conv hold the values in use.
        integer(ip)            :: frac_mode        = 0
        integer(ip)            :: frac_evap        = 0
        integer(ip)            :: frac_lsc         = 0
        integer(ip)            :: frac_conv        = 0

! dmr&clo   Initial vapour (replaces the datmini parameters; d-excess 0 as in the legacy code)
        real(dblp)             :: d18o_vap_init    = -12.0_dblp      !< permil
        real(dblp)             :: dexc_vap_init    =   0.0_dblp      !< permil
        real(dblp)             :: o17exc_vap_init  =   0.0_dblp      !< per meg

! dmr&clo   Fixed composition of the evaporation flux over land and snow (D10, until land coupling, step 3)
        real(dblp)             :: d18o_land_evap   = -12.0_dblp      !< permil
        real(dblp)             :: dexc_land_evap   =   0.0_dblp      !< permil
        real(dblp)             :: o17exc_land_evap =   0.0_dblp      !< per meg
        real(dblp)             :: d18o_snow_evap   = -12.0_dblp      !< permil
        real(dblp)             :: dexc_snow_evap   =   0.0_dblp      !< permil
        real(dblp)             :: o17exc_snow_evap =   0.0_dblp      !< per meg

! dmr&clo   Surface ocean composition: 'file' (d18O map, dD and 17O from dexc_ocean, o17exc_ocean) or 'fixed' (uniform)
        character(len=str_len) :: ocean_ratio_mode = 'file'
        character(len=str_len) :: ocean_d18o_file  = 'inputdata/isoatm/calculated_d18O_v1_1-T21.nc'
        real(dblp)             :: d18o_ocean       =   0.0_dblp      !< permil, 'fixed' only
        real(dblp)             :: dexc_ocean       =   0.0_dblp      !< permil, both modes
        real(dblp)             :: o17exc_ocean     =   0.0_dblp      !< per meg, both modes

! dmr&clo   wiso:3a:0001 Land isotopes: 'fixed' (evaporation over land and snow at the compositions above, the LBM
!           carries water only) or 'bucket' (isotopes in the soil water and the snow of the LBM, ROADMAP step 3a);
!           wiso:3a:0008 default 'bucket' (D70, run bf of the 3a evaluation, O63)
        character(len=str_len) :: land_iso_mode    = 'bucket'
!           wiso:3a:0002 'bucket': composition of the soil water and the snow at a cold start of the land isotopes (no
!           land isotope restart): the isotopes are set to these ratios times the water of the stores
        real(dblp)             :: d18o_soil_init   =  -8.0_dblp      !< permil
        real(dblp)             :: dexc_soil_init   =  10.0_dblp      !< permil
        real(dblp)             :: o17exc_soil_init =   0.0_dblp      !< per meg
        real(dblp)             :: d18o_snow_init   = -25.0_dblp      !< permil
        real(dblp)             :: dexc_snow_init   =  10.0_dblp      !< permil
        real(dblp)             :: o17exc_snow_init =   0.0_dblp      !< per meg
!           wiso:3a:0007 cold start from these values ('fixed') or from a map ('file', D69): netCDF on the model grid,
!           d18o_snow (required), dexc_snow, o17exc_snow, d18o_soil, dexc_soil, o17exc_soil (permil, per meg; optional),
!           the values above where a variable or a cell is missing. Made offline: tools/src/wisolbm_init
!           wiso:3a:0008 default 'file', with the map of run bf installed from inputs/isoatm (D70); ignored with 'fixed'
        character(len=str_len) :: land_init_mode   = 'file'
        character(len=str_len) :: land_init_file   = 'inputdata/isoatm/wisolbm_init.nc'

! dmr&clo   Restart of the isotopologues (netCDF, D34; available from ROADMAP 2c)
        logical                :: restart_read     = .false.
        character(len=str_len) :: restart_file     = 'startdata/wisoatm_restart.nc'
! dmr&clo   wiso:3a:0006 Restart of the land isotopes (land_iso_mode = 'bucket'; netCDF, D65), independent of the above
        logical                :: lbm_restart_read = .false.
        character(len=str_len) :: lbm_restart_file = 'startdata/wisolbm_restart.nc'

! dmr&clo   Fractionation factors (D35, iso_alphas): formulations by name, 17O exponents, snow supersaturation
        character(len=ISO_FORM_LEN) :: alpha_lv_form = 'majoube71'      !< liquid-vapour: majoube71 | horita94
        character(len=ISO_FORM_LEN) :: alpha_sv_form = 'majoube70_mn67' !< ice-vapour: majoube70_mn67 | ellehoj13
        character(len=ISO_FORM_LEN) :: diff_form     = 'merlivat78'     !< diffusivities: merlivat78 | cappa03
        real(dblp)             :: theta_lv17       = 0.529_dblp      !< alpha17 = alpha18**theta, liquid-vapour
        real(dblp)             :: theta_sv17       = 0.528_dblp      !< idem, ice-vapour
        real(dblp)             :: theta_diff17     = 0.518_dblp      !< D17/D = (D18/D)**theta (Barkan and Luz 2007)
        real(dblp)             :: ss_s0            = 1.0_dblp        !< supersaturation over ice S = ss_s0 - ss_s1 T(degC)
        real(dblp)             :: ss_s1            = 0.002_dblp      !< (Jouzel and Merlivat 1984; 0.002 to 0.004);
                                                                         !< 0.002 since wiso:2f:0003 (D59, 2f run S2)

! dmr&clo   wiso:2d:0002 Kinetic fractionation of ocean evaporation (frac_evap = 2, Craig-Gordon; O4)
!           kin_form = 'mj79': Merlivat and Jouzel (1979, eq. 9), R_E = (1 - k)(R_oc/alpha - h R_v)/(1 - h); k18 =
!           kin_k18 (smooth regime, read from MJ79 fig. 2; k18 is not given as a number in the paper); k_i = k18 eps_i/
!           eps_18 with eps_i = D/D_i - 1 from diff_form (MJ79: kD/k18 = eps_D/eps_18 = 0.88 with Merlivat 1978), for
!           17O always, for D unless kin_kd_ratio >= 0 is given; kin_wind = .true.: rough regime above 7 m/s,
!           k18 = 0.285 u10 + 0.82 permil (a fit of MJ79 fig. 2, source to verify; no wind dependence by default:
!           Bonne et al. 2019, Zannoni et al. 2022). 'diffusivity' (Mathieu and Bariac style) reserved.
        character(len=ISO_FORM_LEN) :: kin_form = 'mj79'
        logical                :: kin_wind         = .false.
        real(dblp)             :: kin_k18          = 6.0_dblp        !< permil, smooth regime of MJ79
! dmr&clo   wiso:2d:0003 Condensation (frac_lsc, frac_conv; wiso_cond_mod): temperature 't500' or 'temp4g'; phase 'ramp'
!           (liquid fraction linear between cond_tice and cond_tliq, degC) or 'water' (rain/snow of the water by tsurf,
!           Roche 2013)
        character(len=ISO_FORM_LEN) :: cond_temp  = 't500'
        character(len=ISO_FORM_LEN) :: cond_phase = 'ramp'
        real(dblp)             :: cond_tliq        = 0.0_dblp        !< degC, all liquid above
        real(dblp)             :: cond_tice        = -20.0_dblp      !< degC, all ice below
        real(dblp)             :: cond_tliq_k      = 273.15_dblp     !< the same in K (after reading)
        real(dblp)             :: cond_tice_k      = 253.15_dblp

! dmr&clo   wiso:2f:0002 Near-surface vapour ratio in the Craig-Gordon evaporation (O52): R_v = R_col alpha_lv(T_s)**gamma;
!           0 = the column ratio (D55); gamma = 1 for an exponential humidity profile distilled from the surface (idealised)
        real(dblp)             :: cg_rv_gamma      = 0.0_dblp

! dmr&clo   wiso:2f:0004 Below-cloud re-equilibration of rain with the vapour (O5, O54): fraction of the step's large-scale
!           (bc_eq_lsc) and convective (bc_eq_conv) rain brought to isotopic equilibrium with the vapour at tempsg, without
!           exchange of water (closed system per cell); 0: none. Snow does not exchange.
        real(dblp)             :: bc_eq_lsc        = 0.0_dblp
        real(dblp)             :: bc_eq_conv       = 0.0_dblp

        real(dblp)             :: kin_kd_ratio     = -1.0_dblp       !< kD / k18; < 0: eps_D/eps_18 from diff_form

! dmr&clo   Plausibility window of the delta-bound diagnostic (wiso_budget_mod, wiso:2c:0006); 17O from 18O with 0.528
        real(dblp)             :: d18o_check_min   = -100.0_dblp     !< permil
        real(dblp)             :: d18o_check_max   =   50.0_dblp     !< permil
        real(dblp)             :: dd_check_min     = -800.0_dblp     !< permil
        real(dblp)             :: dd_check_max     =  400.0_dblp     !< permil

        namelist /wisoatm_nml/ frac_mode, frac_evap, frac_lsc, frac_conv,                                             &
                               d18o_vap_init, dexc_vap_init, o17exc_vap_init,                                         &
                               d18o_land_evap, dexc_land_evap, o17exc_land_evap,                                      &
                               d18o_snow_evap, dexc_snow_evap, o17exc_snow_evap,                                      &
                               ocean_ratio_mode, ocean_d18o_file, d18o_ocean, dexc_ocean, o17exc_ocean,               &
                               land_iso_mode, d18o_soil_init, dexc_soil_init, o17exc_soil_init,                       &
                               d18o_snow_init, dexc_snow_init, o17exc_snow_init, land_init_mode, land_init_file,      &
                               restart_read, restart_file, lbm_restart_read, lbm_restart_file,                        &
                               alpha_lv_form, alpha_sv_form, diff_form, theta_lv17, theta_sv17, theta_diff17,         &
                               ss_s0, ss_s1, d18o_check_min, d18o_check_max, dd_check_min, dd_check_max,              &
                               kin_form, kin_wind, kin_k18, kin_kd_ratio, cond_temp, cond_phase, cond_tliq, cond_tice, &
                               cg_rv_gamma, bc_eq_lsc, bc_eq_conv

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Derived ratios (set by iso_nml_read)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp), dimension(nwisos) :: r_vap_init    = 1.0_dblp     !< initial vapour
        real(dblp), dimension(nwisos) :: r_land_evap   = 1.0_dblp     !< evaporation flux over land
        real(dblp), dimension(nwisos) :: r_snow_evap   = 1.0_dblp     !< evaporation (sublimation) flux over snow
        real(dblp), dimension(nwisos) :: r_ocean_fixed = 1.0_dblp     !< surface ocean, ocean_ratio_mode = 'fixed'
        logical                       :: land_bucket   = .false.      !< land_iso_mode = 'bucket' (wiso:3a:0002)
        real(dblp), dimension(nwisos) :: r_soil_init   = 1.0_dblp     !< soil water, cold start ('bucket')
        real(dblp), dimension(nwisos) :: r_snow_init   = 1.0_dblp     !< snow, cold start ('bucket')
        real(dblp)                    :: dexc_ocean_frac   = 0.0_dblp  !< dexc_ocean as a fraction ('file' mode)
        real(dblp)                    :: o17exc_ocean_frac = 0.0_dblp  !< o17exc_ocean as a fraction ('file' mode)
        real(dblp), dimension(nwisos) :: delta_check_min = -1.0_dblp  !< plausibility window, delta as a fraction
        real(dblp), dimension(nwisos) :: delta_check_max =  1.0_dblp  !< (iwat16: +-1e-6, the deviation from water)
        real(dblp)                    :: kin_k18_frac  = 6.0e-3_dblp  !< kin_k18 as a fraction (wiso:2d:0002)


      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Delta values to ratios of all species
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        pure function iso_ratios(d18o, dexc, o17exc) result(ratios)

          real(dblp), intent(in)        :: d18o                         !< delta 18O, fraction
          real(dblp), intent(in)        :: dexc                         !< deuterium excess, fraction
          real(dblp), intent(in)        :: o17exc                       !< 17O excess, fraction
          real(dblp), dimension(nwisos) :: ratios

          ratios(iwater) = 1.0_dblp
          ratios(iwat16) = 1.0_dblp
          ratios(iwat18) = 1.0_dblp + d18o
          ratios(iwat2h) = 1.0_dblp + ISO_SLOPE_DEX * d18o + dexc
          ratios(iwat17) = exp(ISO_SLOPE_17 * log(ratios(iwat18)) + o17exc)

        end function iso_ratios

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Reading and checks
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine iso_nml_read(adv_mode)

          character(len=*), intent(in) :: adv_mode                     !< water transport mode (moipar, D28)

          integer(ip)        :: unit_nml, ios
          logical            :: exists, ok
          character(len=256) :: msg

          inquire(file=ISO_NML_FILE, exist=exists)
          if (.not. exists) then
            write(*,*) "iso_nml_read: ", ISO_NML_FILE, " not found in the run directory (needed with WISO = 1)"
            stop 1
          endif

          open(newunit=unit_nml, file=ISO_NML_FILE, status='old', action='read', iostat=ios, iomsg=msg)
          if (ios /= 0) then
            write(*,*) "iso_nml_read: cannot open ", ISO_NML_FILE, ": ", trim(msg)
            stop 1
          endif
          read(unit_nml, nml=wisoatm_nml, iostat=ios, iomsg=msg)
          if (ios /= 0) then
            write(*,*) "iso_nml_read: error in group wisoatm_nml of ", ISO_NML_FILE, ": ", trim(msg)
            stop 1
          endif
          close(unit_nml)

          ! dmr&clo   isotope runs use the FV water transport (D30, O31)
          if (trim(adv_mode) /= 'fv') then
            write(*,*) "iso_nml_read: WISO = 1 requires adv_mode = 'fv' in namelist group moipar (found '",          &
                       trim(adv_mode), "')"
            stop 1
          endif

          ! dmr&clo   wiso:2d:0000 one switch for all processes, or per process with frac_mode = -1
          if (frac_mode < -1 .or. frac_mode > 2) then
            write(*,*) "iso_nml_read: frac_mode = ", frac_mode, " (expected -1, 0, 1 or 2)"
            stop 1
          endif
          if (frac_mode >= 0) then
            frac_evap = frac_mode
            frac_lsc  = frac_mode
            frac_conv = frac_mode
          endif
          call check_frac('frac_evap', frac_evap, 2)
          call check_frac('frac_lsc',  frac_lsc,  2)
          call check_frac('frac_conv', frac_conv, 2)
#if ( EVAPSI != 1 )
          ! dmr&clo   wiso:2d:0001 evaporation over sea ice (EVAPSI = 0): source and water bookkeeping to settle first (O45)
          if (frac_evap /= 0) then
            write(*,*) "iso_nml_read: frac_evap = ", frac_evap, " requires EVAPSI = 1 (no evaporation over sea ice, O45)"
            stop 1
          endif
#endif

          select case (trim(ocean_ratio_mode))
          case ('file', 'fixed')
          case default
            write(*,*) "iso_nml_read: unknown ocean_ratio_mode = ", trim(ocean_ratio_mode), " (file | fixed)"
            stop 1
          end select

          select case (trim(land_init_mode))
          case ('fixed', 'file')
          case default
            write(*,*) "iso_nml_read: unknown land_init_mode = ", trim(land_init_mode), " (fixed | file)"
            stop 1
          end select
          select case (trim(land_iso_mode))
          case ('fixed')
            if (lbm_restart_read) then
              write(*,*) "iso_nml_read: lbm_restart_read = .true. needs land_iso_mode = 'bucket'"
              stop 1
            endif
! dmr&clo   wiso:3a:0008 land_init_mode = 'file' is the default (D70): with 'fixed' it has nothing to set, ignored
            if (trim(land_init_mode) == 'file') then
              write(*,*) "iso_nml_read: land_iso_mode = 'fixed', land_init_mode = 'file' ignored"
            endif
          case ('bucket')
            land_bucket = .true.
          case default
            write(*,*) "iso_nml_read: unknown land_iso_mode = ", trim(land_iso_mode), " (fixed | bucket)"
            stop 1
          end select

! dmr&clo   wiso:2c:0002 restart_read: read by wiso_restart_mod at the end of ec_iatmphys

! dmr&clo   wiso:2d:0002 kinetic settings of the ocean evaporation
          select case (trim(kin_form))
          case ('mj79')
          case ('diffusivity')
            write(*,*) "iso_nml_read: kin_form = 'diffusivity' is not available yet (ROADMAP 2d, O4)"
            stop 1
          case default
            write(*,*) "iso_nml_read: unknown kin_form = ", trim(kin_form), " (mj79)"
            stop 1
          end select
          if (kin_k18 < 0.0_dblp .or. kin_k18 > 30.0_dblp .or. kin_kd_ratio > 2.0_dblp) then
            write(*,*) "iso_nml_read: kin_k18 = ", kin_k18, " permil (0 to 30), kin_kd_ratio = ", kin_kd_ratio,         &
                       " (< 0: from diff_form, or 0 to 2)"
            stop 1
          endif
          kin_k18_frac = kin_k18 * ISO_PERMIL
          if (bc_eq_lsc < 0.0_dblp .or. bc_eq_lsc > 1.0_dblp .or. bc_eq_conv < 0.0_dblp .or. bc_eq_conv > 1.0_dblp) then
            write(*,*) "iso_nml_read: bc_eq_lsc = ", bc_eq_lsc, ", bc_eq_conv = ", bc_eq_conv, " (0 to 1)"
            stop 1
          endif
          if (cg_rv_gamma < 0.0_dblp .or. cg_rv_gamma > 2.0_dblp) then
            write(*,*) "iso_nml_read: cg_rv_gamma = ", cg_rv_gamma, " (0 to 2)"
            stop 1
          endif

! dmr&clo   wiso:2d:0003 condensation settings
          if (cond_temp /= 't500' .and. cond_temp /= 'temp4g') then
            write(*,*) "iso_nml_read: unknown cond_temp = ", trim(cond_temp), " (t500, temp4g)"
            stop 1
          endif
          if (cond_phase /= 'ramp' .and. cond_phase /= 'water') then
            write(*,*) "iso_nml_read: unknown cond_phase = ", trim(cond_phase), " (ramp, water)"
            stop 1
          endif
          if (cond_tliq <= cond_tice) then
            write(*,*) "iso_nml_read: cond_tliq = ", cond_tliq, " must be above cond_tice = ", cond_tice
            stop 1
          endif
          cond_tliq_k = cond_tliq + tK_zero_C
          cond_tice_k = cond_tice + tK_zero_C

          call iso_alphas_init(alpha_lv_form, alpha_sv_form, diff_form, theta_lv17, theta_sv17, theta_diff17,         &
                               ss_s0, ss_s1, ok, msg)
          if (.not. ok) then
            write(*,*) "iso_nml_read: ", trim(msg)
            stop 1
          endif

          ! dmr&clo   namelist permil / per meg -> fractions, the only place where these units are converted
          r_vap_init    = iso_ratios(d18o_vap_init  * ISO_PERMIL, dexc_vap_init  * ISO_PERMIL, o17exc_vap_init  * ISO_PERMEG)
          r_land_evap   = iso_ratios(d18o_land_evap * ISO_PERMIL, dexc_land_evap * ISO_PERMIL, o17exc_land_evap * ISO_PERMEG)
          r_snow_evap   = iso_ratios(d18o_snow_evap * ISO_PERMIL, dexc_snow_evap * ISO_PERMIL, o17exc_snow_evap * ISO_PERMEG)
          r_soil_init   = iso_ratios(d18o_soil_init * ISO_PERMIL, dexc_soil_init * ISO_PERMIL, o17exc_soil_init * ISO_PERMEG)
          r_snow_init   = iso_ratios(d18o_snow_init * ISO_PERMIL, dexc_snow_init * ISO_PERMIL, o17exc_snow_init * ISO_PERMEG)
          r_ocean_fixed = iso_ratios(d18o_ocean     * ISO_PERMIL, dexc_ocean     * ISO_PERMIL, o17exc_ocean     * ISO_PERMEG)
          dexc_ocean_frac   = dexc_ocean   * ISO_PERMIL
          o17exc_ocean_frac = o17exc_ocean * ISO_PERMEG
          delta_check_min(iwat16) = -1.0e-6_dblp
          delta_check_max(iwat16) =  1.0e-6_dblp
          delta_check_min(iwat18) = d18o_check_min * ISO_PERMIL
          delta_check_max(iwat18) = d18o_check_max * ISO_PERMIL
          delta_check_min(iwat17) = exp(0.528_dblp*log(1.0_dblp + delta_check_min(iwat18))) - 1.0_dblp
          delta_check_max(iwat17) = exp(0.528_dblp*log(1.0_dblp + delta_check_max(iwat18))) - 1.0_dblp
          delta_check_min(iwat2h) = dd_check_min * ISO_PERMIL
          delta_check_max(iwat2h) = dd_check_max * ISO_PERMIL

          write(*,*) "iso_nml_read: water isotopes, frac_evap, frac_lsc, frac_conv = ", frac_evap, frac_lsc, frac_conv,  &
                     ", ocean_ratio_mode = ", trim(ocean_ratio_mode), ", land_iso_mode = ", trim(land_iso_mode)

        contains

          ! one per-process switch: 0, 1 or 2; values above vmax not available until their 2d step is in
          subroutine check_frac(name, val, vmax)
            character(len=*), intent(in) :: name
            integer(ip),      intent(in) :: val
            integer,          intent(in) :: vmax              !< highest value available (wiso:2d:0001: 1 for frac_evap)
            if (val < 0 .or. val > 2) then
              write(*,*) "iso_nml_read: ", name, " = ", val, " (expected 0, 1 or 2)"
              stop 1
            endif
            if (val > vmax) then
              write(*,*) "iso_nml_read: ", name, " = ", val, " is not available yet (at most ", vmax, ", ROADMAP 2d)"
              stop 1
            endif
          end subroutine check_frac

        end subroutine iso_nml_read

        !> No fractionation anywhere: every ratio is a mixture of the sources (strict check of wiso_budget_mod)
        pure function iso_passive() result(p)
          logical :: p
          p = frac_evap == 0 .and. frac_lsc == 0 .and. frac_conv == 0
        end function iso_passive

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Report to parameters.chk
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine iso_nml_write(unit_chk)

          integer(ip), intent(in) :: unit_chk                          !< open formatted unit (parameters.chk)

          write(unit_chk, '(a)')         '--- water isotopes (wisoatm.nml) ---'
          write(unit_chk, '(a,i0)')      'frac_mode        = ', frac_mode
          write(unit_chk, '(a,3(1x,i0))') 'in use: frac_evap, frac_lsc, frac_conv =', frac_evap, frac_lsc, frac_conv
          write(unit_chk, '(a,3f12.4)')  'vap_init  d18o, dexc (permil), o17exc (per meg) = ',                          &
                                         d18o_vap_init, dexc_vap_init, o17exc_vap_init
          write(unit_chk, '(a,3f12.4)')  'land_evap d18o, dexc (permil), o17exc (per meg) = ',                          &
                                         d18o_land_evap, dexc_land_evap, o17exc_land_evap
          write(unit_chk, '(a,3f12.4)')  'snow_evap d18o, dexc (permil), o17exc (per meg) = ',                          &
                                         d18o_snow_evap, dexc_snow_evap, o17exc_snow_evap
          write(unit_chk, '(a,a)')       'ocean_ratio_mode = ', trim(ocean_ratio_mode)
          write(unit_chk, '(a,a)')       'ocean_d18o_file  = ', trim(ocean_d18o_file)
          write(unit_chk, '(a,3f12.4)')  'ocean     d18o (fixed only), dexc (permil), o17exc (per meg) = ',             &
                                         d18o_ocean, dexc_ocean, o17exc_ocean
          write(unit_chk, '(a,a)')       'land_iso_mode    = ', trim(land_iso_mode)
          write(unit_chk, '(a,3f12.4)')  'soil_init d18o, dexc (permil), o17exc (per meg) = ',                          &
                                         d18o_soil_init, dexc_soil_init, o17exc_soil_init
          write(unit_chk, '(a,3f12.4)')  'snow_init d18o, dexc (permil), o17exc (per meg) = ',                          &
                                         d18o_snow_init, dexc_snow_init, o17exc_snow_init
          write(unit_chk, '(a,l1)')      'restart_read     = ', restart_read
          write(unit_chk, '(a,a)')       'restart_file     = ', trim(restart_file)
          write(unit_chk, '(a,a)')       'land_init_mode   = ', trim(land_init_mode)
          write(unit_chk, '(a,a)')       'land_init_file   = ', trim(land_init_file)
          write(unit_chk, '(a,l1)')      'lbm_restart_read = ', lbm_restart_read
          write(unit_chk, '(a,a)')       'lbm_restart_file = ', trim(lbm_restart_file)
          write(unit_chk, '(a,a)')       'alpha_lv_form    = ', trim(alpha_lv_form)
          write(unit_chk, '(a,a)')       'alpha_sv_form    = ', trim(alpha_sv_form)
          write(unit_chk, '(a,a)')       'diff_form        = ', trim(diff_form)
          write(unit_chk, '(a,3f8.4)')   'theta17 lv, sv, diff = ', theta_lv17, theta_sv17, theta_diff17
          write(unit_chk, '(a,2f8.4)')   'supersaturation s0, s1 (per degC) = ', ss_s0, ss_s1
          write(unit_chk, '(a,4f10.2)')  'check window d18o min max, dD min max (permil) = ',                           &
                                         d18o_check_min, d18o_check_max, dd_check_min, dd_check_max
          write(unit_chk, '(a,a,a,l1,a,f8.3,a,f8.4)') 'kinetic (frac_evap = 2): kin_form = ', trim(kin_form),         &
                                         ', kin_wind = ', kin_wind, ', kin_k18 (permil) = ', kin_k18,                    &
                                         ', kin_kd_ratio = ', kin_kd_ratio
          write(unit_chk, '(a,f7.3)') 'Craig-Gordon near-surface vapour: cg_rv_gamma = ', cg_rv_gamma
          write(unit_chk, '(a,f6.3,a,f6.3)') 'below-cloud re-equilibration of rain: bc_eq_lsc = ', bc_eq_lsc,             &
                                             ', bc_eq_conv = ', bc_eq_conv
          write(unit_chk, '(a,a,a,a,a,f7.2,a,f7.2)') 'condensation: cond_temp = ', trim(cond_temp),                   &
                                         ', cond_phase = ', trim(cond_phase), ', cond_tliq (degC) = ', cond_tliq,        &
                                         ', cond_tice (degC) = ', cond_tice
          call iso_alphas_report(unit_chk)

        end subroutine iso_nml_write

#endif

      end module iso_nml_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
