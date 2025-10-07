!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"
!dmr -- Adding the choice of components through the pre-processing options
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Copyright 2021 Didier M. Roche (a.k.a. dmr)

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
!
!      Auteur : Didier M. Roche 
!      Date   : 25 mai 2011
!      Derniere modification : 25 mai 2011, 23 juin 2021 (version 21)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE ec_isoatminit


       use global_constants_mod, only: ip, dblp=>dp
       
#if ( WISOATM == 1 )       

       use comsurf_mod, only: fractoc
       USE iso_param_mod, ONLY : dexcess, ratio_oceanatm, rsmow, delta_inv, fac_17Oexc
       USE lectNC, only: lect_2D_NC
       use comatm, only: nlat, nlon, nwisos, iwater, iwat17, iwat18, iwat2h

       IMPLICIT NONE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       REAL, DIMENSION(nlon, nlat) :: dataread
       REAL, DIMENSION(nlat, nlon) :: dataecb, Ones_Ecbilt = 1._dblp

       character(len=4), parameter :: varia="d18o"
       character(len=44), parameter:: pathtofile = "inputdata/isoatm/calculated_d18O_v1_1-T21.nc"
       INTEGER :: i,j

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Let's read the input file containing surface water d18O
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
       CALL lect_2D_NC(TRIM(pathtofile),dataread,var_nc = TRIM(varia)) 

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Reading the input file leads to an array indexed differently as ours. We thus reformat it. 
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       CALL format_nc_to_ec(dataread,dataecb,nlon,nlat)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       The dataset read has negative missing values that are incompatible with our formulation. Change this. 
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       WHERE(dataecb.LT.-1E10_dblp)
         dataecb = ABS(dataecb)
       ENDWHERE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Now dataecb contains a "delta 18O" read in pathtofile
!       Convert d18O to ratio 
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       ratio_oceanatm(:,:,iwat18) = delta_inv(Ones_Ecbilt, dataecb, iwat18)

! dmr
!       Then we assume a d-excess of 0 and a d17excess of zero
!     
!       d-excess is defined as : d-excess = dD - 8*d18O
!         choosing d-excess = 0.0 implies dD = 8 * d18O
!
!       --- 
!
!       17Oexcess is defined as : 
!           17Oexcess = ln(d17O/1000.+1.0) - 0.528*ln(d18O/1000.+1.0)
!         choosing 17Oexcess = 0.0 thus implies : 
!            R17O = ((d18O/1000.+1.0) **0.528) * R17smow
!         or : R17O = R18**0.528 * R17SMOW * R18SMOW**-0.528
! dmr

       ratio_oceanatm(:,:,iwat2h) = delta_inv(Ones_Ecbilt, 8._dblp*dataecb, iwat2h) ! (dataecb/1000.d0*8.d0 +1.d0)*rsmow(iwat2h)
       ratio_oceanatm(:,:,iwat17) = ratio_oceanatm(:,:,iwat18)**fac_17Oexc * (rsmow(iwat17) / rsmow(iwat18)**(fac_17Oexc)) ! ((dataecb/1000.d0+1.d0)**0.528d0)*rsmow(iwat17)

!~ #if (0) 
!~        DO i=1, iso_LAT
!~          DO j=1, iso_LON
!~            if (ratio_oceanatm(i,j,ieau18).LT.9999.d0)                   &
!~              WRITE(*,'(3G15.5)')                                        &
!~              (ratio_oceanatm(i,j,ieau17)/rsmow(ieau17)-1.d0)*1000.d0,   &
!~              (ratio_oceanatm(i,j,ieau18)/rsmow(ieau18)-1.d0)*1000.d0,   &
!~              (ratio_oceanatm(i,j,ieaud)/rsmow(ieaud)-1.d0)*1000.d0
!~ !           if (ratio_oceanatm(i,j,ieau18).LT.9999.d0)                   &
!~ !             WRITE(*,*) "d-excess: ",                                   &
!~ !           dexcess((ratio_oceanatm(i,j,ieaud)/rsmow(ieaud)-1.d0)*1000.d0&
!~ !               ,(ratio_oceanatm(i,j,ieau18)/rsmow(ieau18)-1.d0)*1000.d0)
!~          ENDDO
!~        ENDDO
 
!~ #endif


      call fix_input_field(fractoc,ratio_oceanatm,iwat17,iwat2h)


#endif

      contains
      
      
       SUBROUTINE fix_input_field(data_in,data_out,nzi,nzj)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!                            data_in, le champs de reference
!                            nx,ny ; la taille de tableau d'entree
!       Variables de sortie : 
!                            data_out, le champs a corriger
!-----|--1--------2---------3---------4---------5---------6---------7-|

!   INSERER ICI LES EVENTUELS "USE MODULE"

       IMPLICIT NONE

!   INSERER ICI LES DECLARATIONS DES VARIABLES D'ENTREE / SORTIE

       REAL, DIMENSION(:,:),   INTENT(IN)    :: data_in  ! nx, ny
       REAL, DIMENSION(:,:,:), INTENT(INOUT) :: data_out ! nx, ny, nz
 
!~        INTEGER, INTENT(IN) :: nx, ny, nz, nzi, nzj
       INTEGER, INTENT(IN) :: nzi, nzj

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

        REAL, PARAMETER :: undef = 10000.0D0
        REAL :: val_moy
        INTEGER :: nb_moy
        INTEGER :: i, j, k, l, ii, jj, kk, ll, m, n, mm, nn, z
        INTEGER :: nx, ny, nz
        INTEGER, PARAMETER :: s = 1
!-----|--1--------2---------3---------4---------5---------6---------7-|
!-----|--1--------2---------3---------4---------5---------6---------7-|

       nx = UBOUND(data_out,DIM=1)
       ny = UBOUND(data_out,DIM=2)
       nz = UBOUND(data_out,DIM=3)
       
       DO z = nzi, nzj

       DO i=1,nx
        DO j=1,ny

         IF ((data_in(i,j).NE.0.0d0).AND.(data_out(i,j,z).GE.undef))    &
         THEN

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       On effectue la moyenne des plus proches avec une valeur
!       pertinente (code repris de routageEAU.f90)
!-----|--1--------2---------3---------4---------5---------6---------7-|

          val_moy = 0.0d0
          nb_moy = 0

          ii = i
          jj = j
          m = -s
          n = s

          DO k=m, n, n
           DO l=m, n, n

              mm = ii+k
              nn = jj+l

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Attention, il ne faut pas etre circulaire en latitude !!! 
!-----|--1--------2---------3---------4---------5---------6---------7-|
              if (mm.LT.1) mm = 1  
              if (mm.GT.nx) mm = nx

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       En revanche il FAUT etre circulaire en longitude !!! 
!-----|--1--------2---------3---------4---------5---------6---------7-|

              if (nn.LT.1) nn = ny + nn
              if (nn.GT.ny) nn = nn - ny

              if (.NOT.(data_out(mm,nn,z).GT.undef)) then
                val_moy = val_moy + data_out(mm,nn,z)
                nb_moy = nb_moy + 1
              endif
           ENDDO
          ENDDO

          IF(nb_moy.NE.0) THEN
            data_out(i,j,z) = val_moy / nb_moy
          ELSE
            WRITE(*,*) "I can't fix location: ", i, j
            READ(*,*)
          ENDIF

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Fin de la partie "fix" du code
!-----|--1--------2---------3---------4---------5---------6---------7-|
         ENDIF

        ENDDO
       ENDDO

       ENDDO

       END SUBROUTINE fix_input_field
      

       END SUBROUTINE ec_isoatminit



! dmr --- The End of All Things (op. cit.)
