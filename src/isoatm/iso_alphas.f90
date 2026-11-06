!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2021, 2026 Didier M. Roche (a.k.a. dmr) | iLOVECLIM / FRATRES coding group

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
!      MODULE: [iso_alphas]
!
!>     @brief Fractionation factors of the water isotopologues, with selectable formulations (ROADMAP D35, step 2b).
!
!      DESCRIPTION:
!
!>     Equilibrium factors are condensed phase over vapour (alpha > 1 for the heavy species), temperature in K.
!!     Species index k as in comatm (iwater, iwat16, iwat17, iwat18, iwat2h); iwater and iwat16 always give 1.
!!     17O is derived from 18O with the exponents theta: alpha17 = alpha18**theta (equilibrium) and
!!     D17/D = (D18/D)**theta_diff (diffusivities).
!!
!!     Formulations (names as given in wisoatm.nml):
!!       liquid-vapour  alpha_lv_form = 'majoube71'      Majoube (1971), 18O and D
!!                                      'horita94'       planned: Horita and Wesolowski (1994), coefficients to verify
!!       ice-vapour     alpha_sv_form = 'majoube70_mn67' Majoube (1970) for 18O, Merlivat and Nief (1967) for D
!!                                      'ellehoj13'      Ellehoj et al. (2013), coefficients from the abstract, to verify
!!       diffusivities  diff_form     = 'merlivat78'     Merlivat (1978)
!!                                      'cappa03'        Cappa et al. (2003)
!!       kinetic effect at snow formation (Jouzel and Merlivat 1984): supersaturation S = ss_s0 - ss_s1 * T(degC)
!!
!!     Corrections relative to the v21 code (wiso:2b:0002, sources checked with dmr, Nov 2026):
!!       - Majoube (1971) D liquid-vapour constant +52.612e-3 (the code had +56.612e-3: alpha_D 4 permil too large);
!!       - Majoube (1970) 18O ice-vapour constant 28.224e-3 (the code had 28.244e-3);
!!       - ice-vapour attributions: 18O Majoube (1970), D Merlivat and Nief (1967) (the code cited Majoube 1971a).
!!
!!     Contract
!!       - iso_alphas_init(...) once, before any use; returns ok = .false. and a message for an unknown or not yet
!!         available formulation or an invalid parameter (the caller stops).
!!       - alpha_lv, alpha_sv, alpha_k_snow, supersat are pure elemental; diff_ratio is pure elemental.
!!       - alpha_sv(T,k) * alpha_k_snow(T,k) is the effective ice-vapour factor of Jouzel and Merlivat (1984).
!
!      REFERENCES: ROADMAP D35, O37, O42; references in ROADMAP §8.
!
!      AUTHORS: dmr (v21 formulations), dmr&clo (dispatch layer, 2026)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module iso_alphas

#if ( WISOATM == 1 )

        use global_constants_mod, only: ip, dblp=>dp, tK_zero_C
        use comatm,               only: iwat17, iwat18, iwat2h

        implicit none

        private

        public :: ISO_FORM_LEN, iso_alphas_init, iso_alphas_report
        public :: alpha_lv, alpha_sv, alpha_k_snow, supersat, diff_ratio


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip), parameter :: ISO_FORM_LEN = 32                     !< length of a formulation name

        integer(ip), parameter :: LV_MAJOUBE71      = 1
        integer(ip), parameter :: SV_MAJOUBE70_MN67 = 1
        integer(ip), parameter :: SV_ELLEHOJ13      = 2
        integer(ip), parameter :: DIFF_MERLIVAT78   = 1
        integer(ip), parameter :: DIFF_CAPPA03      = 2

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Module state (set by iso_alphas_init)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip), save :: lv_code     = LV_MAJOUBE71
        integer(ip), save :: sv_code     = SV_MAJOUBE70_MN67
        integer(ip), save :: diff_code   = DIFF_MERLIVAT78
        real(dblp),  save :: theta_lv    = 0.529_dblp                   !< alpha17_lv = alpha18_lv**theta_lv
        real(dblp),  save :: theta_sv    = 0.528_dblp                   !< alpha17_sv = alpha18_sv**theta_sv
        real(dblp),  save :: theta_diff  = 0.518_dblp                   !< D17/D = (D18/D)**theta_diff
        real(dblp),  save :: ss_s0       = 1.0_dblp                     !< supersaturation S = s0 - s1 T(degC)
        real(dblp),  save :: ss_s1       = 0.004_dblp
        real(dblp),  save :: d18_ratio   = 0.9723_dblp                  !< D(H2 18O)/D(H2O)
        real(dblp),  save :: d2h_ratio   = 0.9755_dblp                  !< D(HDO)/D(H2O)


      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Set-up and report
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine iso_alphas_init(lv_form, sv_form, diff_form, th_lv, th_sv, th_diff, s0, s1, ok, msg)

          character(len=*), intent(in)  :: lv_form                     !< liquid-vapour formulation
          character(len=*), intent(in)  :: sv_form                     !< ice-vapour formulation
          character(len=*), intent(in)  :: diff_form                   !< diffusivity ratios
          real(dblp),       intent(in)  :: th_lv                       !< 17O exponent, liquid-vapour
          real(dblp),       intent(in)  :: th_sv                       !< 17O exponent, ice-vapour
          real(dblp),       intent(in)  :: th_diff                     !< 17O exponent, diffusivities
          real(dblp),       intent(in)  :: s0                          !< supersaturation at 0 degC
          real(dblp),       intent(in)  :: s1                          !< supersaturation slope, per degC
          logical,          intent(out) :: ok
          character(len=*), intent(out) :: msg

          ok  = .true.
          msg = ''

          select case (trim(lv_form))
          case ('majoube71')
            lv_code = LV_MAJOUBE71
          case ('horita94')
            ok = .false.
            msg = "alpha_lv_form = 'horita94' is planned, not available yet (coefficients to verify, O37)"
            return
          case default
            ok = .false.
            msg = "unknown alpha_lv_form = '" // trim(lv_form) // "' (majoube71)"
            return
          end select

          select case (trim(sv_form))
          case ('majoube70_mn67')
            sv_code = SV_MAJOUBE70_MN67
          case ('ellehoj13')
            sv_code = SV_ELLEHOJ13
          case default
            ok = .false.
            msg = "unknown alpha_sv_form = '" // trim(sv_form) // "' (majoube70_mn67 | ellehoj13)"
            return
          end select

          select case (trim(diff_form))
          case ('merlivat78')
            diff_code = DIFF_MERLIVAT78
            d18_ratio = 0.9723_dblp
            d2h_ratio = 0.9755_dblp
          case ('cappa03')
            diff_code = DIFF_CAPPA03
            d18_ratio = 0.9691_dblp
            d2h_ratio = 0.9839_dblp
          case default
            ok = .false.
            msg = "unknown diff_form = '" // trim(diff_form) // "' (merlivat78 | cappa03)"
            return
          end select

          if (th_lv <= 0.0_dblp .or. th_sv <= 0.0_dblp .or. th_diff <= 0.0_dblp .or. s0 <= 0.0_dblp) then
            ok = .false.
            msg = "theta_lv17, theta_sv17, theta_diff17 and ss_s0 must be positive"
            return
          endif

          theta_lv   = th_lv
          theta_sv   = th_sv
          theta_diff = th_diff
          ss_s0      = s0
          ss_s1      = s1

        end subroutine iso_alphas_init

        subroutine iso_alphas_report(unit_out)

          integer(ip), intent(in) :: unit_out                          !< open formatted unit

          write(unit_out, '(a,f8.4,a,f8.4)') 'alpha check, 20 degC: 18O lv = ', alpha_lv(tK_zero_C + 20.0_dblp, iwat18), &
                                             ', D lv = ', alpha_lv(tK_zero_C + 20.0_dblp, iwat2h)
          write(unit_out, '(a,f8.4,a,f8.4)') 'alpha check, -20 degC: 18O sv = ', alpha_sv(tK_zero_C - 20.0_dblp, iwat18), &
                                             ', D sv = ', alpha_sv(tK_zero_C - 20.0_dblp, iwat2h)
          write(unit_out, '(a,3f8.5)')       'D_i/D (17O, 18O, D) = ', diff_ratio(iwat17), diff_ratio(iwat18),            &
                                             diff_ratio(iwat2h)

        end subroutine iso_alphas_report

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Liquid-vapour equilibrium
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        pure elemental function alpha_lv(temp, k) result(alpha)

          real(dblp),  intent(in) :: temp                              !< temperature, K
          integer(ip), intent(in) :: k                                 !< species index (comatm)
          real(dblp)              :: alpha

          if (k == iwat18) then
            alpha = lv18(temp)
          else if (k == iwat17) then
            alpha = lv18(temp)**theta_lv
          else if (k == iwat2h) then
            alpha = lv2h(temp)
          else
            alpha = 1.0_dblp                                           ! iwater, iwat16
          endif

        end function alpha_lv

        pure elemental function lv18(temp) result(alpha)

          real(dblp), intent(in) :: temp
          real(dblp)             :: alpha

          ! dmr   Majoube (1971): ln alpha = 1137/T**2 - 0.4156/T - 2.0667e-3 (lv_code = LV_MAJOUBE71, the only one so far)
          alpha = exp(1137.0_dblp / temp**2 - 0.4156_dblp / temp - 2.0667e-3_dblp)

        end function lv18

        pure elemental function lv2h(temp) result(alpha)

          real(dblp), intent(in) :: temp
          real(dblp)             :: alpha

          ! dmr&clo   Majoube (1971): ln alpha = 24844/T**2 - 76.248/T + 52.612e-3 (v21 code had +56.612e-3)
          alpha = exp(24844.0_dblp / temp**2 - 76.248_dblp / temp + 52.612e-3_dblp)

        end function lv2h

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Ice-vapour equilibrium
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        pure elemental function alpha_sv(temp, k) result(alpha)

          real(dblp),  intent(in) :: temp                              !< temperature, K
          integer(ip), intent(in) :: k                                 !< species index (comatm)
          real(dblp)              :: alpha

          if (k == iwat18) then
            alpha = sv18(temp)
          else if (k == iwat17) then
            alpha = sv18(temp)**theta_sv
          else if (k == iwat2h) then
            alpha = sv2h(temp)
          else
            alpha = 1.0_dblp                                           ! iwater, iwat16
          endif

        end function alpha_sv

        pure elemental function sv18(temp) result(alpha)

          real(dblp), intent(in) :: temp
          real(dblp)             :: alpha

          select case (sv_code)
          case (SV_ELLEHOJ13)
            ! dmr&clo   Ellehoj et al. (2013): ln alpha = 0.0831 - 49.192/T + 8312.5/T**2 (to verify against the paper)
            alpha = exp(0.0831_dblp - 49.192_dblp / temp + 8312.5_dblp / temp**2)
          case default
            ! dmr&clo   Majoube (1970): ln alpha = 11.839/T - 28.224e-3 (v21 code had 28.244e-3)
            alpha = exp(11.839_dblp / temp - 28.224e-3_dblp)
          end select

        end function sv18

        pure elemental function sv2h(temp) result(alpha)

          real(dblp), intent(in) :: temp
          real(dblp)             :: alpha

          select case (sv_code)
          case (SV_ELLEHOJ13)
            ! dmr&clo   Ellehoj et al. (2013): ln alpha = 0.2133 - 203.10/T + 48888/T**2 (to verify against the paper)
            alpha = exp(0.2133_dblp - 203.10_dblp / temp + 48888.0_dblp / temp**2)
          case default
            ! dmr&clo   Merlivat and Nief (1967): ln alpha = 16288/T**2 - 0.0934 (constant to check, O42)
            alpha = exp(16288.0_dblp / temp**2 - 0.0934_dblp)
          end select

        end function sv2h

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Kinetic effect at snow formation (Jouzel and Merlivat 1984) and diffusivities
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        pure elemental function supersat(temp) result(s)

          real(dblp), intent(in) :: temp                               !< temperature, K
          real(dblp)             :: s

          s = ss_s0 - ss_s1 * (temp - tK_zero_C)

        end function supersat

        pure elemental function alpha_k_snow(temp, k) result(alpha)

          real(dblp),  intent(in) :: temp                              !< temperature, K
          integer(ip), intent(in) :: k                                 !< species index (comatm)
          real(dblp)              :: alpha

          real(dblp) :: s

          ! dmr   alpha_k = S / (1 + alpha_sv (S - 1) D/D_i); the effective factor is alpha_sv * alpha_k
          s = supersat(temp)
          alpha = s / (1.0_dblp + alpha_sv(temp, k) * (s - 1.0_dblp) / diff_ratio(k))

        end function alpha_k_snow

        pure elemental function diff_ratio(k) result(ratio)

          integer(ip), intent(in) :: k                                 !< species index (comatm)
          real(dblp)              :: ratio                             !< D_i / D(H2O)

          if (k == iwat18) then
            ratio = d18_ratio
          else if (k == iwat17) then
            ratio = d18_ratio**theta_diff
          else if (k == iwat2h) then
            ratio = d2h_ratio
          else
            ratio = 1.0_dblp                                           ! iwater, iwat16
          endif

        end function diff_ratio

#endif

      end module iso_alphas

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
