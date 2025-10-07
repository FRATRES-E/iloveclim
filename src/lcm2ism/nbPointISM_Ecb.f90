!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Cette routine sert a calculer le nombre de points ISM dans une
!      case du modele ECBilt
!      Cette routine est utilisee dans le transfert des variables ISM
!      sur la grille ECBilt
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche
!      Date   : 18 Juillet 2008
!      Derniere modification : 13 Janvier 2009
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#include "choixcomposantes.h"

#if ( DOWNSCALING == 2 || ISM > 1 )

       SUBROUTINE nbPointISM_Ecb(nx,ny,ecx,ecy,iEcb,jEcb,nbpoints,nmbp,coord_GrisEcb,coord_EcbGris, ec_lat_low, ec_lat_hig,        &
                                 ec_lon_low, ec_lon_hig)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  :
!        nx,ny      : dimension de la grille ISM
!        ecx, ecy   : dimension de la grille ECBilt
!        iEcb, jEcb : tableaux des indices ECBilt sur grille ISM
!       Variables de sortie :
!        nbpoints   : tableau du nombre de points ISM sur cases ECBilt
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       IMPLICIT NONE

       INTEGER, INTENT(IN) :: nx,ny,ecx,ecy,nmbp
       REAL, DIMENSION(nx,ny), INTENT(IN) :: iEcb, jEcb
       INTEGER, DIMENSION(ecx,ecy), INTENT(OUT) :: nbpoints

       INTEGER, DIMENSION(3,nx,ny), INTENT(OUT) :: coord_GrisEcb
       INTEGER, DIMENSION(2,nmbp,ecx,ecy), INTENT(OUT) :: coord_EcbGris

       INTEGER, INTENT(out) :: ec_lat_low, ec_lat_hig, ec_lon_low, ec_lon_hig

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       INTEGER i,j, ei,ej

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Initialisation
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       nbpoints(:,:) = 0
       coord_EcbGRis(:,:,:,:) = 0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       DO j=1,ny
         DO i=1,nx

!check           if (i.eq. 141.and.j.eq.131) WRITE(*,*) "PROBLEM :",          &
!check     &    iEcb(i,j),i,j
!check           if (FLOOR(iEcb(i,j)).GT.32) WRITE(*,*) "PROBLEM :",          &
!check     &    iEcb(i,j),i,j

           ei = MIN(ecx,FLOOR(iEcb(i,j)))
           ej = FLOOR(jEcb(i,j))

           nbpoints(ei,ej) = nbpoints(ei,ej) + 1

           coord_GrisEcb(1,i,j) = ei
           coord_GrisEcb(2,i,j) = ej
           coord_GrisEcb(3,i,j) = nbpoints(ei,ej)

           coord_EcbGRis(1,nbpoints(ei,ej),ei,ej) = i
           coord_EcbGRis(2,nbpoints(ei,ej),ei,ej) = j

         ENDDO
       ENDDO

       ec_lat_low = MINVAL(PACK(coord_GrisEcb(1,:,:),.true.),dim=1)
       ec_lat_hig = MAXVAL(PACK(coord_GrisEcb(1,:,:),.true.),dim=1)
       ec_lon_low = MINVAL(PACK(coord_GrisEcb(2,:,:),.true.),dim=1)
       ec_lon_hig = MAXVAL(PACK(coord_GrisEcb(2,:,:),.true.),dim=1)



       WRITE(*,*) "Extent of subGRID in main GRID :: ", ec_lat_low, ec_lat_hig, ec_lon_low, ec_lon_hig

       END SUBROUTINE nbPointISM_Ecb

#endif
