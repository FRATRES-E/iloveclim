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
!      MODULE: [wiso_evap_mod]
!
!>     @brief Isotopologue evaporation fluxes by surface type (ROADMAP step 2c, wiso:2c:0003).
!
!      DESCRIPTION:
!
!>     The isotopologue flux is the water flux at the composition of its source:
!!       - open ocean, and sea ice when it evaporates (EVAPSI = 0; D41): surface ocean ratios ratio_oceanatm
!!         (map or fixed, ocean_ratio_mode of wisoatm.nml); frac_evap = 1 (wiso:2d:0001): vapour in equilibrium
!!         with the surface ocean, ratio_oceanatm / alpha_lv(T_ocean), T_ocean the surface temperature of the
!!         ocean tile (that of qsurfn); frac_evap = 2 (wiso:2d:0002): Craig-Gordon with the closure of Merlivat and
!!         Jouzel (1979), R_E = (1 - k)(R_oc/alpha - h R_v)/(1 - h), h = q10/q_s the humidity of the bulk formula of
!!         the water flux, R_v the (column) vapour ratio, k18 from iso_nml_mod (smooth regime, or rough above 7 m/s
!!         with kin_wind), k_i = k18 eps_i/eps_18 with eps_i = D/D_i - 1 (MJ79 eq. 11 to first order, which gives
!!         their kD/k18 = 0.88; 17O consistently), kD optionally from kin_kd_ratio; computed as the flux
!!         E/(1 - h) (1 - k)(R_oc/alpha - h R_v), finite as h -> 1; a negative
!!         isotopologue flux (h R_v > R_oc/alpha) is set to 0 and counted (wisobounds_*, neclip); iwat16 follows
!!         the water. Fractionation is refused with EVAPSI = 0 (iso_nml_mod);
!!       - land (D10, D42): sublimation of snow at r_snow_evap, evaporation of the bucket at r_land_evap, fixed in
!!         wisoatm.nml until the land isotopes (step 3). The heat flux esice over the ice mask carries no water, hence
!!         no isotopologues.
!!
!!     Contract
!!       - wiso_evap_split(ev, i, j, e_ocean, e_snow, e_land[, t_ocean, q_s, q_a, u10]): ev(iwat16:iwat2h) of cell
!!         (i,j) from the parts of the water flux (m/s) by source; e_ocean + e_snow + e_land is the water flux
!!         evapn(i,j,nn,iwater) of the caller; where e_ocean /= 0, t_ocean (K) is required for frac_evap >= 1, and
!!         q_s, q_a (saturation and 10 m humidity of the bulk formula, e_ocean > 0 only if q_s > q_a) and u10 (m/s)
!!         for frac_evap = 2.
!
!      REFERENCES: ROADMAP D10, D36, D41, D42, D50, O4. Merlivat, L. and Jouzel, J. (1979), J. Geophys. Res., 84,
!                  5029-5033; Bonne et al. (2019), Nat. Commun., 10, 1632; Zannoni et al. (2022), J. Geophys. Res.
!                  Atmos. (no smooth/rough discontinuity observed).
!
!      AUTHORS: dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module wiso_evap_mod

#if ( WISOATM == 1 )

        use global_constants_mod, only: ip, dblp=>dp
        use comatm,               only: iwat16, iwat17, iwat18, iwat2h

        implicit none

        private

        public :: wiso_evap_split

        real(dblp), parameter :: MJ79_U_ROUGH = 7.0_dblp            !< m/s, smooth -> rough regime (MJ79)
        real(dblp), parameter :: MJ79_RA      = 0.285e-3_dblp       !< rough regime k18 = RA u10 + RB (fractions)
        real(dblp), parameter :: MJ79_RB      = 0.82e-3_dblp

      contains

        subroutine wiso_evap_split(ev, i, j, e_ocean, e_snow, e_land, t_ocean, q_s, q_a, u10)

          use iso_param_mod,   only: ratio_oceanatm
          use iso_nml_mod,     only: r_snow_evap, r_land_evap, frac_evap
          use iso_alphas,      only: alpha_lv
          use wiso_vap_mod,    only: wiso_vap_ratio_at
          use wiso_budget_mod, only: wiso_bud_count_eclip

          real(dblp),  intent(inout) :: ev(:)          !< (nwisos) evaporation of the cell, iwat16:iwat2h set
          integer(ip), intent(in)    :: i, j           !< cell
          real(dblp),  intent(in)    :: e_ocean        !< water flux from the surface ocean (m/s)
          real(dblp),  intent(in)    :: e_snow         !< water flux from snow (m/s)
          real(dblp),  intent(in)    :: e_land         !< water flux from the land bucket (m/s)
          real(dblp),  intent(in), optional :: t_ocean !< surface temperature of the ocean tile (K), frac_evap >= 1
          real(dblp),  intent(in), optional :: q_s     !< saturation humidity at t_ocean (bulk formula), frac_evap = 2
          real(dblp),  intent(in), optional :: q_a     !< 10 m humidity (bulk formula), frac_evap = 2
          real(dblp),  intent(in), optional :: u10     !< 10 m wind speed (m/s), frac_evap = 2

          integer(ip) :: k
          real(dblp)  :: h, ecg, eo

          if (frac_evap == 0 .or. e_ocean == 0.0_dblp) then
            ! passive: the flux at the composition of each source (wiso:2c:0003)
            do k = iwat16, iwat2h
              ev(k) = e_ocean*ratio_oceanatm(i, j, k) + e_snow*r_snow_evap(k) + e_land*r_land_evap(k)
            enddo
          else if (frac_evap == 2) then
            ! dmr&clo   wiso:2d:0002 Craig-Gordon / MJ79; e_ocean > 0 implies q_s > q_a, so h < 1
            if (.not. (present(t_ocean) .and. present(q_s) .and. present(q_a) .and. present(u10)))                    &
              error stop "wiso_evap_split: t_ocean, q_s, q_a, u10 needed for frac_evap = 2"
            h   = q_a/q_s
            ecg = e_ocean*q_s/(q_s - q_a)                       ! E/(1 - h)
            ev(iwat16) = e_ocean*ratio_oceanatm(i, j, iwat16) + e_snow*r_snow_evap(iwat16) + e_land*r_land_evap(iwat16)
            do k = iwat17, iwat2h
              eo = ecg*mj79_one_minus_k(u10, k)*(ratio_oceanatm(i, j, k)/alpha_lv(t_ocean, k)                         &
                                                 - h*wiso_vap_ratio_at(i, j, k))
              if (eo < 0.0_dblp) then
                eo = 0.0_dblp
                call wiso_bud_count_eclip(k)
              endif
              ev(k) = eo + e_snow*r_snow_evap(k) + e_land*r_land_evap(k)
            enddo
          else
            ! dmr&clo   wiso:2d:0001 frac_evap = 1: ocean flux in equilibrium with the surface ocean (alpha = 1 for iwat16)
            if (.not. present(t_ocean)) error stop "wiso_evap_split: t_ocean needed for frac_evap >= 1"
            do k = iwat16, iwat2h
              ev(k) = e_ocean*ratio_oceanatm(i, j, k)/alpha_lv(t_ocean, k) + e_snow*r_snow_evap(k)                 &
                    + e_land*r_land_evap(k)
            enddo
          endif

        end subroutine wiso_evap_split

        !> 1 - k of Merlivat and Jouzel (1979) for species k (iwat17:iwat2h): smooth regime k18 = kin_k18, or with
        !! kin_wind the rough regime above MJ79_U_ROUGH; k_i = k18 eps_i/eps_18, eps_i = D/D_i - 1 (diff_form);
        !! kD = kin_kd_ratio k18 if kin_kd_ratio >= 0
        pure function mj79_one_minus_k(u10, k) result(omk)

          use iso_nml_mod, only: kin_wind, kin_k18_frac, kin_kd_ratio
          use iso_alphas,  only: diff_ratio

          real(dblp),  intent(in) :: u10
          integer(ip), intent(in) :: k
          real(dblp)              :: omk

          real(dblp) :: k18, eps18

          k18 = kin_k18_frac
          if (kin_wind .and. u10 > MJ79_U_ROUGH) k18 = MJ79_RA*u10 + MJ79_RB
          eps18 = 1.0_dblp/diff_ratio(iwat18) - 1.0_dblp
          select case (k)
          case (iwat18)
            omk = 1.0_dblp - k18
          case (iwat2h)
            if (kin_kd_ratio >= 0.0_dblp) then
              omk = 1.0_dblp - kin_kd_ratio*k18
            else
              omk = 1.0_dblp - k18*(1.0_dblp/diff_ratio(iwat2h) - 1.0_dblp)/eps18
            endif
          case (iwat17)
            omk = 1.0_dblp - k18*(1.0_dblp/diff_ratio(iwat17) - 1.0_dblp)/eps18
          case default
            omk = 1.0_dblp
          end select

        end function mj79_one_minus_k

#endif

      end module wiso_evap_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
