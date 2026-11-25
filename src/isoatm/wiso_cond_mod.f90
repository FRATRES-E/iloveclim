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
!      MODULE: [wiso_cond_mod]
!
!>     @brief Fractionation factor of the condensate, condensate / vapour (ROADMAP step 2d, wiso:2d:0003).
!
!      DESCRIPTION:
!
!>     The condensate of a cell forms at a condensation temperature T_c (cond_temp of wisoatm.nml):
!!       - 't500': 500 hPa, interpolated in log-pressure between temp4g (650 hPa) and temp2g (350 hPa) exactly as the
!!         water's saturation content (ec_detqmaxexact);
!!       - 'temp4g': 650 hPa.
!!     It is a mixture of liquid (fraction w) and ice, each in equilibrium with the vapour, so that
!!       alpha_eff = w alpha_lv(T_c) + (1 - w) alpha_ice(T_c),   alpha_ice = alpha_sv           (frac = 1)
!!                                                              alpha_ice = alpha_sv alpha_k(S) (frac = 2, Jouzel and
!!                                                              Merlivat 1984, S = ss_s0 - ss_s1 T_c(degC))
!!     with w (cond_phase):
!!       - 'ramp': w = (T_c - cond_tice)/(cond_tliq - cond_tice), clamped to [0, 1] (alpha_lv extrapolated below 0 degC
!!         for supercooled liquid);
!!       - 'water': the phase of the water's precipitation, w = 1 if tsurf >= tzero + 0.1 (rain), 0 otherwise (snow),
!!         as in Roche (2013); alpha at T_c.
!!     The condensate goes to rain and snow as the water's (melting does not change its composition; no below-cloud
!!     exchange, O5). iwat16: alpha = 1.
!!
!!     Contract
!!       - wiso_cond_temp(t2, t4): T_c (K) from temp2g, temp4g of a cell.
!!       - wiso_cond_liq(tc, ts): w from T_c and tsurf (K) of a cell.
!!       - wiso_cond_alpha(tc, w, k, frac): alpha_eff of species k; frac = 0 gives 1.
!
!      REFERENCES: ROADMAP D46, D47, D50, D56, O1. Jouzel, J. and Merlivat, L. (1984), J. Geophys. Res., 89, 11749-11757;
!                  Roche, D. M. (2013), Geosci. Model Dev., 6, 1481-1491.
!
!      AUTHORS: dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module wiso_cond_mod

#if ( WISO == 1 )

        use global_constants_mod, only: ip, dblp=>dp
        use comatm,               only: iwat16, rlogtl12

        implicit none

        private

        public :: wiso_cond_temp, wiso_cond_liq, wiso_cond_alpha

        real(dblp), parameter :: TS_RAIN   = 0.1_dblp                      !< rain if tsurf >= tzero + 0.1 (ec_moisture)

      contains

        !> Condensation temperature (K) of a cell
        function wiso_cond_temp(t2, t4) result(tc)

          use iso_nml_mod, only: cond_temp

          real(dblp), intent(in) :: t2                 !< temp2g, 350 hPa (K)
          real(dblp), intent(in) :: t4                 !< temp4g, 650 hPa (K)
          real(dblp)             :: tc

          real(dblp) :: alpha

          if (cond_temp == 't500') then
            ! as ec_detqmaxexact (atmphys_mod)
            alpha = (t2 - t4)*rlogtl12
            tc    = t4 + alpha*log(500e0_dblp/650e0_dblp)
          else
            tc    = t4
          endif

        end function wiso_cond_temp

        !> Liquid fraction of the condensate
        function wiso_cond_liq(tc, ts) result(w)

          use iso_nml_mod, only: cond_phase, cond_tliq_k, cond_tice_k
          use comphys,     only: tzero

          real(dblp), intent(in) :: tc                 !< condensation temperature (K)
          real(dblp), intent(in) :: ts                 !< tsurf of the cell (K), for cond_phase = 'water'
          real(dblp)             :: w

          if (cond_phase == 'water') then
            w = 0.0_dblp
            if (ts >= tzero + TS_RAIN) w = 1.0_dblp
          else
            w = max(0.0_dblp, min(1.0_dblp, (tc - cond_tice_k)/(cond_tliq_k - cond_tice_k)))
          endif

        end function wiso_cond_liq

        !> Fractionation factor condensate / vapour of species k
        function wiso_cond_alpha(tc, w, k, frac) result(alpha)

          use iso_alphas, only: alpha_lv, alpha_sv, alpha_k_snow

          real(dblp),  intent(in) :: tc                !< condensation temperature (K)
          real(dblp),  intent(in) :: w                 !< liquid fraction
          integer(ip), intent(in) :: k                 !< species (comatm)
          integer(ip), intent(in) :: frac              !< 0 none, 1 equilibrium, 2 + kinetic effect for the ice
          real(dblp)              :: alpha

          real(dblp) :: aice

          if (frac == 0 .or. k == iwat16) then
            alpha = 1.0_dblp
            return
          endif
          aice = 1.0_dblp
          if (w < 1.0_dblp) then
            aice = alpha_sv(tc, k)
            if (frac == 2) aice = aice*alpha_k_snow(tc, k)
          endif
          if (w > 0.0_dblp) then
            alpha = w*alpha_lv(tc, k) + (1.0_dblp - w)*aice
          else
            alpha = aice
          endif

        end function wiso_cond_alpha

#endif

      end module wiso_cond_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
