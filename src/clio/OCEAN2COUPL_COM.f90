!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Copyright for initial code Hugues Goosse, Thierry Fichefet, https://www.elic.ucl.ac.be/modx/index.php?id=289
!      Branched from version 1.2

!   Further development:
!   Copyright 2023 Didier M. Roche (a.k.a. dmr)

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

      MODULE OCEAN2COUPL_COM

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       USE global_constants_mod, ONLY: dblp=>dp, ip
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      IMPLICIT NONE
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   History
! dmr           Change from 0.0.0: Created a modular version from the initial no-module legacy code
! dmr           Change from 0.1.0: split into gather (here, CLIO -> ocn_bndcon on the T21 grid) and apply (ocean_bndcon_mod,
! dmr                              coupler side). detseaalb, initseaalb, oc2at, ec_shine moved to ocean_bndcon_mod.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      CHARACTER(LEN=5), PARAMETER :: version_mod ="0.2.0"
      

      PRIVATE
      PUBLIC :: ec_oc2co, ec_oc2co_gather
      
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Tentative list of variables communicated from ocean to atm:
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
      
      CONTAINS
      
      
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_oc2co(ist)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** communicate oceanic data to the coupler: gather (CLIO -> ocn_bndcon) then apply (ocn_bndcon -> surface state)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use ocean_bndcon_mod, only: ocean_bndcon_apply

      implicit none

      integer ist

      call ec_oc2co_gather()
      call ocean_bndcon_apply()

      return
      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_oc2co_gather()
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** interpolate the CLIO surface fields needed by the atmosphere onto the atmospheric grid (ocn_bndcon)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!! START_OF_USE_SECTION

      use ocean_bndcon_mod, only: ocn_bndcon, BNDCON_SST, BNDCON_SIC, BNDCON_TSI, BNDCON_HIC, BNDCON_HSN, oc2at

#if ( ISOOCN >= 2 )
      USE isoatm_mod, ONLY: ratio_oceanatm
      USE iso_param_mod, ONLY: ieau17, ieaud
      USE para0_mod, only: owisostrt, owisostop, oc2atwisoindx
#endif

      use para0_mod, only: imax, jmax
      use para_mod,  only:
      use bloc0_mod, only: ks2, scal
      use bloc_mod , only:
      use ice_mod  , only: albq, ts, hgbq, hnbq

!! END_OF_USE_SECTION

      implicit none

      integer ix,jy
#if ( ISOOCN >= 2 )
      INTEGER :: iz, iatmwiso
#endif
      real*8 zfld(imax,jmax)

! *** SST

      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = scal(ix,jy,ks2,1)
        enddo
      enddo
      call oc2at(zfld,ocn_bndcon(:,:,BNDCON_SST))

#if ( ISOOCN >= 2 )
! [NOTA] Loop here should be on the ocean part, hence owisostrt -> owisostop
      do iz = owisostrt, owisostop
! *** 17O, 18O, DH
      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = scal(ix,jy,ks2,iz)
        enddo
      enddo
      ! [NOTA] in the after-following call, the index should be oceanic
      !        need to convert the index from oceanic to atmospheric
      !        (only 17 -> 2H)
      call oc2atwisoindx(iz,iatmwiso)
      call oc2at(zfld,ratio_oceanatm(:,:,iatmwiso))
      enddo ! on iz, wisos ...
#endif

! *** ICE fraction

      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = 1.0-albq(ix,jy)
        enddo
      enddo
      call oc2at(zfld,ocn_bndcon(:,:,BNDCON_SIC))

! *** STI

      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = ts(ix,jy)
        enddo
      enddo
      call oc2at(zfld,ocn_bndcon(:,:,BNDCON_TSI))

! *** hic

      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = hgbq(ix,jy)
        enddo
      enddo
      call oc2at(zfld,ocn_bndcon(:,:,BNDCON_HIC))

! *** hsn

      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = hnbq(ix,jy)
        enddo
      enddo
      call oc2at(zfld,ocn_bndcon(:,:,BNDCON_HSN))

      return
      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      END MODULE OCEAN2COUPL_COM
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
