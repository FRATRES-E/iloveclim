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

!   Style sheet: v1.0.0 (standalone program, by analogy, D25)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      PROGRAM: [test_iso_alphas]
!
!>     @brief Pass/fail tests of the fractionation factors of src/isoatm/iso_alphas.f90 (ROADMAP 2b, D35).
!
!      DESCRIPTION:
!
!>     Reference values computed independently (python, from the published formulas) at -20, 0 and 20 degC, to 1e-12;
!!     textbook anchors at 20 degC (Majoube 1971: alpha18_lv ~ 1.0098, alphaD_lv ~ 1.085); structural checks:
!!     17O exponents, iwater and iwat16 give 1, S = 1 gives no kinetic effect, unknown names are refused.
!!     One line per check, exit code 1 on any failure.
!
!      AUTHORS: dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      program test_iso_alphas

        use global_constants_mod, only: ip, dblp=>dp, tK_zero_C
        use comatm,               only: iwater, iwat16, iwat17, iwat18, iwat2h
        use iso_alphas,           only: iso_alphas_init, alpha_lv, alpha_sv, alpha_k_snow, supersat, diff_ratio

        implicit none

        real(dblp), parameter :: TOL = 1.0e-12_dblp
        real(dblp), parameter :: TT(3) = [253.15_dblp, 273.15_dblp, 293.15_dblp]

        ! dmr&clo   python references (see DESCRIPTION), columns = TT
        real(dblp), parameter :: LV18(3) = [1.014132602567857_dblp, 1.011718982796772_dblp, 1.009793879292320_dblp]
        real(dblp), parameter :: LV2H(3) = [1.149225425386506_dblp, 1.112321652295485_dblp, 1.085031301017711_dblp]
        real(dblp), parameter :: SV18(3) = [1.018715723219610_dblp, 1.015233348238943_dblp, 1.012235719541133_dblp]
        real(dblp), parameter :: SV2H(3) = [1.174406305367668_dblp, 1.133041759168199_dblp, 1.100906976960226_dblp]
        real(dblp), parameter :: EL18(3) = [1.018663143654563_dblp, 1.014524243447619_dblp, 1.012095612027103_dblp]
        real(dblp), parameter :: EL2H(3) = [1.189902790288559_dblp, 1.133139028646432_dblp, 1.093477508536286_dblp]
        real(dblp), parameter :: AK18_M20 = 0.996476307155023_dblp    !< alpha_k 18O at -20 degC, S = 1.08, merlivat78

        integer(ip)        :: nfail, i
        logical            :: ok
        character(len=256) :: msg

        nfail = 0

        call init('majoube71', 'majoube70_mn67', 'merlivat78')
        do i = 1, 3
          call check_close('majoube71 lv 18O',      alpha_lv(TT(i), iwat18), LV18(i), TOL)
          call check_close('majoube71 lv D',        alpha_lv(TT(i), iwat2h), LV2H(i), TOL)
          call check_close('majoube70 sv 18O',      alpha_sv(TT(i), iwat18), SV18(i), TOL)
          call check_close('merlivat-nief67 sv D',  alpha_sv(TT(i), iwat2h), SV2H(i), TOL)
          call check_close('lv 17O = 18O**0.529',   log(alpha_lv(TT(i), iwat17)) / log(LV18(i)), 0.529_dblp, TOL)
          call check_close('sv 17O = 18O**0.528',   log(alpha_sv(TT(i), iwat17)) / log(SV18(i)), 0.528_dblp, TOL)
          call check_close('lv water = 1',          alpha_lv(TT(i), iwater), 1.0_dblp, 0.0_dblp)
          call check_close('lv iwat16 = 1',         alpha_lv(TT(i), iwat16), 1.0_dblp, 0.0_dblp)
          call check_close('sv iwat16 = 1',         alpha_sv(TT(i), iwat16), 1.0_dblp, 0.0_dblp)
        enddo

        ! dmr&clo   textbook anchors at 20 degC (Majoube 1971)
        call check_close('anchor 18O lv 20 degC', alpha_lv(tK_zero_C + 20.0_dblp, iwat18), 1.0098_dblp, 1.0e-4_dblp)
        call check_close('anchor D lv 20 degC',   alpha_lv(tK_zero_C + 20.0_dblp, iwat2h), 1.0850_dblp, 5.0e-4_dblp)

        ! dmr&clo   diffusivities and kinetic effect
        call check_close('merlivat78 D18/D',      diff_ratio(iwat18), 0.9723_dblp, 0.0_dblp)
        call check_close('merlivat78 DD/D',       diff_ratio(iwat2h), 0.9755_dblp, 0.0_dblp)
        call check_close('D17/D = (D18/D)**0.518 (v21: 0.9855)', diff_ratio(iwat17), 0.9723_dblp**0.518_dblp, TOL)
        call check_close('v21 value 0.9855',      diff_ratio(iwat17), 0.9855_dblp, 1.0e-4_dblp)
        call check_close('S(-20 degC) = 1.08',    supersat(tK_zero_C - 20.0_dblp), 1.08_dblp, TOL)
        call check_close('S(0 degC) = 1: alpha_k 18O = 1', alpha_k_snow(tK_zero_C, iwat18), 1.0_dblp, TOL)
        call check_close('alpha_k 18O at -20 degC', alpha_k_snow(tK_zero_C - 20.0_dblp, iwat18), AK18_M20, TOL)
        call check_true('alpha_k < 1 below 0 degC (18O, 17O, D)',                                                     &
                        all(alpha_k_snow(tK_zero_C - 20.0_dblp, [iwat17, iwat18, iwat2h]) < 1.0_dblp))
        call check_close('alpha_k water = 1',     alpha_k_snow(tK_zero_C - 20.0_dblp, iwater), 1.0_dblp, TOL)

        call init('majoube71', 'ellehoj13', 'cappa03')
        do i = 1, 3
          call check_close('ellehoj13 sv 18O',    alpha_sv(TT(i), iwat18), EL18(i), TOL)
          call check_close('ellehoj13 sv D',      alpha_sv(TT(i), iwat2h), EL2H(i), TOL)
        enddo
        call check_close('cappa03 D18/D',         diff_ratio(iwat18), 0.9691_dblp, 0.0_dblp)
        call check_close('cappa03 DD/D',          diff_ratio(iwat2h), 0.9839_dblp, 0.0_dblp)

        ! dmr&clo   refused settings
        call iso_alphas_init('horita94', 'majoube70_mn67', 'merlivat78', 0.529_dblp, 0.528_dblp, 0.518_dblp,           &
                             1.0_dblp, 0.004_dblp, ok, msg)
        call check_true('horita94 refused (planned)', .not. ok)
        call iso_alphas_init('majoube71', 'nonsense', 'merlivat78', 0.529_dblp, 0.528_dblp, 0.518_dblp,                &
                             1.0_dblp, 0.004_dblp, ok, msg)
        call check_true('unknown alpha_sv_form refused', .not. ok)
        call iso_alphas_init('majoube71', 'majoube70_mn67', 'merlivat78', 0.529_dblp, -0.528_dblp, 0.518_dblp,         &
                             1.0_dblp, 0.004_dblp, ok, msg)
        call check_true('negative theta refused', .not. ok)

        if (nfail > 0) then
          write(*,'(a,i0,a)') 'test_iso_alphas: ', nfail, ' FAILED'
          error stop 1
        endif
        write(*,'(a)') 'test_iso_alphas: all passed'

      contains

        subroutine init(lv, sv, df)

          character(len=*), intent(in) :: lv, sv, df

          call iso_alphas_init(lv, sv, df, 0.529_dblp, 0.528_dblp, 0.518_dblp, 1.0_dblp, 0.004_dblp, ok, msg)
          if (.not. ok) error stop 'test_iso_alphas: init failed'

        end subroutine init

        subroutine check_close(name, val, ref, tol)

          character(len=*), intent(in) :: name
          real(dblp),       intent(in) :: val, ref, tol

          if (abs(val - ref) <= tol) then
            write(*,'(a,a,es12.4)') 'PASS  ', name, abs(val - ref)
          else
            write(*,'(a,a,2f20.15)') 'FAIL  ', name, val, ref
            nfail = nfail + 1
          endif

        end subroutine check_close

        subroutine check_true(name, cond)

          character(len=*), intent(in) :: name
          logical,          intent(in) :: cond

          if (cond) then
            write(*,'(a,a)') 'PASS  ', name
          else
            write(*,'(a,a)') 'FAIL  ', name
            nfail = nfail + 1
          endif

        end subroutine check_true

      end program test_iso_alphas

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
