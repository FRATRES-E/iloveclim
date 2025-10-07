!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Copyright 2024 Didier M. Roche (a.k.a. dmr)

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

      MODULE loveclim_transfer_mod

      USE global_constants_mod, ONLY: dblp=>dp, ip
      USE declars_mod,          ONLY: LT, JT, NOC_CBR, JX
      use global_constants_mod, only:days_year360d_i

      implicit none

      INTEGER, PARAMETER         :: KLSR = 0

      REAL(KIND=dblp),  dimension(JX)    :: ZX
      REAL(KIND=dblp),  dimension(JX)    :: mid_level

      REAL(KIND=dblp),  dimension(LT,JT,NOC_CBR), target   :: TM, SM    ! intialisee dans "fait_pointer_CC.f"
      REAL(KIND=dblp),  dimension(LT,JT,NOC_CBR)           :: DVOL      ! intialisee dans "fait_pointer_CC.f"
      REAL(KIND=dblp),  dimension(LT,JT,NOC_CBR)           :: DVOL_prev ! intialisee dans "fait_pointer_CC.f"

      INTEGER(KIND=ip), dimension(LT,JT,NOC_CBR)           :: MGT       ! intialisee dans "fait_pointer_CC.f"
      INTEGER(KIND=ip), dimension(LT,JT,NOC_CBR)           :: MGT_prev  ! intialisee dans "fait_pointer_CC.f"
      INTEGER(KIND=ip), dimension(LT,JT,NOC_CBR)           :: diff_MGT  ! intialisee dans "fait_pointer_CC.f"
      
      REAL(KIND=dblp),  dimension(LT,NOC_CBR)              :: SQRO2     ! intialisee dans "fait_pointer_CC.f"
      REAL(KIND=dblp),  dimension(LT,NOC_CBR), target      :: TM_surface! intialisee dans "fait_pointer_CC.f"
      REAL(KIND=dblp),  dimension(LT,NOC_CBR)              :: SABST_O   ! intialisee dans "fait_pointer_CC.f"
      REAL(KIND=dblp),  dimension(LT,NOC_CBR)              :: FRICE     ! intialisee dans "fait_pointer_CC.f"
      REAL(KIND=dblp),  dimension(LT,NOC_CBR)              :: WS_OC     ! intialisee dans "fait_pointer_CC.f"
      REAL(KIND=dblp),  dimension(LT,NOC_CBR)              :: ZMIX_OCN  ! intialisee dans "fait_pointer_CC.f"

      REAL(KIND=dblp),  dimension(JT)                      :: ZZ        ! intialisee dans "fait_pointer_CC.f"
      REAL(KIND=dblp)                                      :: OVOL      ! Calcule dans "fait_pointer_CC.f
      REAL(KIND=dblp)                                      :: OVOL_prev ! Calcule dans "fait_pointer_CC.f
      REAL(KIND=dblp)                                      :: total_area! Calcule dans "fait_pointer_CC.f
      

      REAL(KIND=dblp), dimension(LT,NOC_CBR,days_year360d_i) :: WIND_ERA5 ! intialisee dans "fait_pointer_CC.f"
      REAL(KIND=dblp), dimension(LT,JT,NOC_CBR)            :: IRON_LIM  ! intialisee dans "fait_pointer_CC.f"

      END MODULE loveclim_transfer_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
