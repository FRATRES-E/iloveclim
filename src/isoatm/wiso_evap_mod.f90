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
!>     frac_evap = 0 (passive, alpha = 1): the isotopologue flux is the water flux at the composition of its source:
!!       - open ocean, and sea ice when it evaporates (EVAPSI = 0; D41): surface ocean ratios ratio_oceanatm
!!         (map or fixed, ocean_ratio_mode of wisoatm.nml);
!!       - land (D10, D42): sublimation of snow at r_snow_evap, evaporation of the bucket at r_land_evap, fixed in
!!         wisoatm.nml until the land isotopes (step 3). The heat flux esice over the ice mask carries no water, hence
!!         no isotopologues.
!!     Equilibrium (frac_evap = 1) and Craig-Gordon (frac_evap = 2) come in step 2d.
!!
!!     Contract
!!       - wiso_evap_split(ev, i, j, e_ocean, e_snow, e_land): ev(iwat16:iwat2h) of cell (i,j) from the parts of the
!!         water flux (m/s) by source; e_ocean + e_snow + e_land is the water flux evapn(i,j,nn,iwater) of the caller.
!
!      REFERENCES: ROADMAP D10, D36, D41, D42.
!
!      AUTHORS: dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module wiso_evap_mod

#if ( WISOATM == 1 )

        use global_constants_mod, only: ip, dblp=>dp
        use comatm,               only: iwat16, iwat2h

        implicit none

        private

        public :: wiso_evap_split

      contains

        subroutine wiso_evap_split(ev, i, j, e_ocean, e_snow, e_land)

          use iso_param_mod, only: ratio_oceanatm
          use iso_nml_mod,   only: r_snow_evap, r_land_evap

          real(dblp),  intent(inout) :: ev(:)          !< (nwisos) evaporation of the cell, iwat16:iwat2h set
          integer(ip), intent(in)    :: i, j           !< cell
          real(dblp),  intent(in)    :: e_ocean        !< water flux from the surface ocean (m/s)
          real(dblp),  intent(in)    :: e_snow         !< water flux from snow (m/s)
          real(dblp),  intent(in)    :: e_land         !< water flux from the land bucket (m/s)

          integer(ip) :: k

          do k = iwat16, iwat2h
            ev(k) = e_ocean*ratio_oceanatm(i, j, k) + e_snow*r_snow_evap(k) + e_land*r_land_evap(k)
          enddo

        end subroutine wiso_evap_split

#endif

      end module wiso_evap_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
