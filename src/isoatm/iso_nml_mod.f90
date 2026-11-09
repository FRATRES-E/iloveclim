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

#if ( WISOATM == 1 )

        use global_constants_mod, only: str_len, ip, dblp=>dp
        use comatm,               only: nwisos, iwater, iwat16, iwat17, iwat18, iwat2h
        use iso_alphas,           only: ISO_FORM_LEN, iso_alphas_init, iso_alphas_report

        implicit none

        private

        public :: ISO_NML_FILE, iso_nml_read, iso_nml_write, iso_ratios
        public :: frac_mode, ocean_ratio_mode, ocean_d18o_file, restart_read, restart_file
        public :: r_vap_init, r_land_evap, r_snow_evap, r_ocean_fixed, dexc_ocean_frac, o17exc_ocean_frac


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

! dmr&clo   Fractionation: 0 none (passive tracers), 1 equilibrium, 2 equilibrium + kinetic; only 0 until ROADMAP 2d
        integer(ip)            :: frac_mode        = 0

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

! dmr&clo   Restart of the isotopologues (netCDF, D34; available from ROADMAP 2c)
        logical                :: restart_read     = .false.
        character(len=str_len) :: restart_file     = 'startdata/wisoatm_restart.nc'

! dmr&clo   Fractionation factors (D35, iso_alphas): formulations by name, 17O exponents, snow supersaturation
        character(len=ISO_FORM_LEN) :: alpha_lv_form = 'majoube71'      !< liquid-vapour: majoube71 | horita94
        character(len=ISO_FORM_LEN) :: alpha_sv_form = 'majoube70_mn67' !< ice-vapour: majoube70_mn67 | ellehoj13
        character(len=ISO_FORM_LEN) :: diff_form     = 'merlivat78'     !< diffusivities: merlivat78 | cappa03
        real(dblp)             :: theta_lv17       = 0.529_dblp      !< alpha17 = alpha18**theta, liquid-vapour
        real(dblp)             :: theta_sv17       = 0.528_dblp      !< idem, ice-vapour
        real(dblp)             :: theta_diff17     = 0.518_dblp      !< D17/D = (D18/D)**theta (Barkan and Luz 2007)
        real(dblp)             :: ss_s0            = 1.0_dblp        !< supersaturation over ice S = ss_s0 - ss_s1 T(degC)
        real(dblp)             :: ss_s1            = 0.004_dblp      !< (Jouzel and Merlivat 1984; 0.002 to 0.004)

        namelist /wisoatm_nml/ frac_mode, d18o_vap_init, dexc_vap_init, o17exc_vap_init,                              &
                               d18o_land_evap, dexc_land_evap, o17exc_land_evap,                                      &
                               d18o_snow_evap, dexc_snow_evap, o17exc_snow_evap,                                      &
                               ocean_ratio_mode, ocean_d18o_file, d18o_ocean, dexc_ocean, o17exc_ocean,               &
                               restart_read, restart_file,                                                            &
                               alpha_lv_form, alpha_sv_form, diff_form, theta_lv17, theta_sv17, theta_diff17,         &
                               ss_s0, ss_s1

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Derived ratios (set by iso_nml_read)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp), dimension(nwisos) :: r_vap_init    = 1.0_dblp     !< initial vapour
        real(dblp), dimension(nwisos) :: r_land_evap   = 1.0_dblp     !< evaporation flux over land
        real(dblp), dimension(nwisos) :: r_snow_evap   = 1.0_dblp     !< evaporation (sublimation) flux over snow
        real(dblp), dimension(nwisos) :: r_ocean_fixed = 1.0_dblp     !< surface ocean, ocean_ratio_mode = 'fixed'
        real(dblp)                    :: dexc_ocean_frac   = 0.0_dblp  !< dexc_ocean as a fraction ('file' mode)
        real(dblp)                    :: o17exc_ocean_frac = 0.0_dblp  !< o17exc_ocean as a fraction ('file' mode)


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
            write(*,*) "iso_nml_read: ", ISO_NML_FILE, " not found in the run directory (needed with WISOATM = 1)"
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
            write(*,*) "iso_nml_read: WISOATM = 1 requires adv_mode = 'fv' in namelist group moipar (found '",          &
                       trim(adv_mode), "')"
            stop 1
          endif

          if (frac_mode /= 0) then
            write(*,*) "iso_nml_read: frac_mode = ", frac_mode, " is not available yet (only 0, ROADMAP 2d)"
            stop 1
          endif

          select case (trim(ocean_ratio_mode))
          case ('file', 'fixed')
          case default
            write(*,*) "iso_nml_read: unknown ocean_ratio_mode = ", trim(ocean_ratio_mode), " (file | fixed)"
            stop 1
          end select

! dmr&clo   wiso:2c:0002 restart_read: read by wiso_restart_mod at the end of ec_iatmphys

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
          r_ocean_fixed = iso_ratios(d18o_ocean     * ISO_PERMIL, dexc_ocean     * ISO_PERMIL, o17exc_ocean     * ISO_PERMEG)
          dexc_ocean_frac   = dexc_ocean   * ISO_PERMIL
          o17exc_ocean_frac = o17exc_ocean * ISO_PERMEG

          write(*,*) "iso_nml_read: water isotopes, frac_mode = ", frac_mode, ", ocean_ratio_mode = ",                    &
                     trim(ocean_ratio_mode)

        end subroutine iso_nml_read

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Report to parameters.chk
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine iso_nml_write(unit_chk)

          integer(ip), intent(in) :: unit_chk                          !< open formatted unit (parameters.chk)

          write(unit_chk, '(a)')         '--- water isotopes (wisoatm.nml) ---'
          write(unit_chk, '(a,i0)')      'frac_mode        = ', frac_mode
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
          write(unit_chk, '(a,l1)')      'restart_read     = ', restart_read
          write(unit_chk, '(a,a)')       'restart_file     = ', trim(restart_file)
          write(unit_chk, '(a,a)')       'alpha_lv_form    = ', trim(alpha_lv_form)
          write(unit_chk, '(a,a)')       'alpha_sv_form    = ', trim(alpha_sv_form)
          write(unit_chk, '(a,a)')       'diff_form        = ', trim(diff_form)
          write(unit_chk, '(a,3f8.4)')   'theta17 lv, sv, diff = ', theta_lv17, theta_sv17, theta_diff17
          write(unit_chk, '(a,2f8.4)')   'supersaturation s0, s1 (per degC) = ', ss_s0, ss_s1
          call iso_alphas_report(unit_chk)

        end subroutine iso_nml_write

#endif

      end module iso_nml_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
