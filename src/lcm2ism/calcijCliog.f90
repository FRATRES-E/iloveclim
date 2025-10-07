!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine calcijCliog est piratee du couplage entre LOVECLIM
!       et les modele de calotte de glace (AG)ISM de P. Huybrechts
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : P. Huybrechts
!      Date   : 15 Juin 2008
!      Derniere modification : 15 Juin 2008
!-----|--1--------2---------3---------4---------5---------6---------7-|


!afq -- Adding the choice of components through the pre-processing opt.
#include "choixcomposantes.h"
!afq -- Adding the choice of components through the pre-processing opt.

       SUBROUTINE calcijCliog(nx,ny,latg,long,latC,lonC,iCg,jCg)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!         nx, ny     : dimension des tableaux sur la grille ISM
!         latg, long : tableaux des lat, lon sur la grille ISM
!         latC, lonC : tableaux des longitudes, latitudes des cases CLIO
!          sur la grille du modele CLIO
!       Variables de sortie : 
!         iCg, jCg   : tableaux de correspondance entre les deux grilles
!           comprenant les indices CLIO sur la grille ISM
!
!         --  initialement pour GISM : nx, ny = 165, 281
!         --  pour GRISLI Nord :       nx, ny = 241, 241
!
!-----|--1--------2---------3---------4---------5---------6---------7-|

#if ( SHELFMELT == 1 )
      use varsClio2ism_mod, only : tmsism
#endif         

      IMPLICIT NONE

      INTEGER, INTENT(in) :: nx,ny
      REAL, DIMENSION(nx,ny), INTENT(in) :: latg,long
      REAL, DIMENSION(122,65), INTENT(in) :: latC,lonC
      INTEGER, DIMENSION(nx,ny), INTENT(out) :: iCg,jCg

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

      INTEGER :: ki,kj
      REAL :: dx,dy,d1,d2,pi,fac,cf,lt,lg
      INTEGER :: i, j, k, l
      REAL :: lokal_lon
      real,dimension(122,65) :: loclonc ! afq
!      INTEGER :: mink, maxk
!-----|--1--------2---------3---------4---------5---------6---------7-|
!                         
!-----|--1--------2---------3---------4---------5---------6---------7-|

!      mink = 65 
!      maxk = 1

      loclonc(:,:)=lonc(:,:)
#if ( SHELFMELT == 1 )
      where (tmsism(:,:,20).eq.0.) 
         loclonC(:,:) = -1.
      end where
#endif


      pi=3.141592654
      fac=pi/180.
      do i=1,nx
       do j=1,ny
!dmr ki = point de depart du calcul de la grille = maximum
!dmr orig        ki=106
         ki=121
!dmr kj = point de depart du calcul de la grille = maximum
!dmr orig         kj=48
         kj=65
!dmr d1 sert plus loin a calculer un minimum, on y stocke au debut
!dmr la valeur maximum
         d1=2*pi
!dmr lt = latitude  de la case ISM i,j en radian
         lt=latg(i,j)*fac
!dmr lg = longitude de la case ISM i,j en radian
         lg=long(i,j)*fac
!dmr si la longitude est moins que pi => on rajoute 2 pi ?
!         if(lg.lt.pi) lg=lg+2*pi
         if(lg.lt.0.0) lg=lg+2*pi
         if(lg.gt.2*pi) lg=lg-2*pi
!dmr cf = cosinus de la latitude
         cf=cos(lt)
!dmr orig         do k=101,106
!dmr Rque : pour la grille GRISLI Nord, le minimum possible est 1  
         do k=1,122
!dmr orig           do l=48,57
!dmr Rque : pour la grille GRISLI Nord, le minimum possible est 39 
           do l=39,65
!dmr on ecrit tout sur [| 0, 360 |]
             if (loclonC(k,l).ge.0.) then
                if (loclonC(k,l).GT.2*pi) then
                   lokal_lon=loclonC(k,l)-2*pi
                else if (lonC(k,l).LT.0.0) then
                   lokal_lon=loclonC(k,l)+2*pi
                else
                   lokal_lon=loclonC(k,l)
                endif
!dmr Distance entre deux meridiens = cos(lat)*(lon1-lon2)
!             dx=cf*(lg-lonC(k,l))
                dx=cf*(lg-lokal_lon)
!dmr Distance entre deux paralleles = (lat1-lat2)
                dy=lt-latC(k,l) 
!dmr d2 = norme du vecteur (dx,dy) = distance vraie entre les deux
!dmr points de la grille
                d2=(dx*dx+dy*dy)**0.5
!dmr si d2 < d1 => d1=d2  <=> d1 est le minimum entre les deux points
!dmr on garde k et l dans ki et kj
                if(d2.le.d1) then
                   d1=d2
                   ki=k
                   kj=l
                endif
             endif ! afq, if loclonc is negative then we go to next point
           end do
         end do
!dmr la version originale conserve la valeur relative de la case Clio ?
!dmr orig         iCg(i,j)=ki-100
!dmr orig         jCg(i,j)=kj-47
         iCg(i,j)=ki
         jCg(i,j)=kj
!         if (ki.gt.maxk) maxk=kj
!         if (ki.lt.mink) mink=kj
       end do
      end do
!         print*, "min, max", maxk, mink
  
      END SUBROUTINE calcijCliog
