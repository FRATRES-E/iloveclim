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

      MODULE MOD_OCEAN_BGC

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~         USE uuid_fort_wrap, ONLY: uuid_size
        USE global_constants_mod, ONLY: dblp=>dp, silp=>sp, ip, str_len, uuid_size

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        IMPLICIT NONE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   History
! dmr           Change from 0.0.0: created a proposal object derived type based on discussion with Guy
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        CHARACTER(LEN=5), PARAMETER :: version_mod ="0.0.1"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        CHARACTER(LEN=7), PARAMETER :: undef_str  = "NotSet"
        REAL(KIND=dblp),  PARAMETER :: undef_dblp = (-1.0_dblp)*HUGE(1.0_silp)
        REAL(KIND=silp),  PARAMETER :: undef_silp = (-1.0_silp)*HUGE(1.0_silp)
        INTEGER(KIND=ip), PARAMETER :: undef_ip   = (-1_ip)*HUGE(1_ip)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  dmr   Holder type for the ocean column biogeochemical properties
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        TYPE :: OCEAN_CHEM_COLUMN

          PRIVATE
          LOGICAL                   :: initialized =.false.
          INTEGER(kind=ip)          :: vert_size   = -1
          INTEGER(kind=ip)          :: eff_bottom  = -1

          ! --- WATER DISSOLVED CHEMISTRY CONTENT
          REAL(kind=dblp), dimension(:)  , allocatable :: OXYGEN
          REAL(kind=dblp), dimension(:)  , allocatable :: PHOSPHATE
          REAL(kind=dblp), dimension(:)  , allocatable :: NITRATE
          REAL(kind=dblp), dimension(:)  , allocatable :: ALKALINITY
          REAL(kind=dblp), dimension(:)  , allocatable :: DISS_INORGC
          REAL(kind=dblp), dimension(:)  , allocatable :: DISS_ORGC
          REAL(kind=dblp), dimension(:)  , allocatable :: DISS_ORGCS
          REAL(kind=dblp), dimension(:)  , allocatable :: PART_ORGC

          ! --- 


          CONTAINS
            PROCEDURE, PUBLIC  :: init => init_ocean_chem_column

        END TYPE OCEAN_CHEM_COLUMN


       CONTAINS


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE init_ocean_chem_column(this, column_size, OPT_effct_col)

       CLASS(OCEAN_CHEM_COLUMN)  , INTENT(OUT):: this
       INTEGER(kind=ip)          , INTENT(IN) :: column_size
       INTEGER(kind=ip), OPTIONAL, INTENT(IN) :: OPT_effct_col



         IF ( .NOT. this%initialized ) then

           IF (PRESENT(OPT_effct_col)) THEN
              this%eff_bottom = OPT_effct_col
           ENDIF

           this%vert_size = column_size

           ALLOCATE(this%OXYGEN(this%vert_size)) 
           ALLOCATE(this%PHOSPHATE(this%vert_size)) 
           ALLOCATE(this%NITRATE(this%vert_size)) 
           ALLOCATE(this%ALKALINITY(this%vert_size)) 
           ALLOCATE(this%DISS_INORGC(this%vert_size)) 
           ALLOCATE(this%DISS_ORGC(this%vert_size)) 
           ALLOCATE(this%DISS_ORGCS(this%vert_size)) 

           ! dmr --- Pre-initialization to HUGE value
           this%OXYGEN(:) = undef_dblp
           this%PHOSPHATE(:) = undef_dblp
           this%NITRATE(:) = undef_dblp
           this%ALKALINITY(:) = undef_dblp
           this%DISS_INORGC(:) = undef_dblp
           this%DISS_ORGC(:) = undef_dblp
           this%DISS_ORGCS(:) = undef_dblp

         ELSE ! I am called to initialize a variable that is already initialized, strange ...

            WRITE(*,*) "un-implemented [ABORT]"
            STOP 1
        ENDIF

        this%initialized = .TRUE.

       END SUBROUTINE init_ocean_chem_column

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      END MODULE MOD_OCEAN_BGC

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
