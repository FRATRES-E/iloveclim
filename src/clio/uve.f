!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:57 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:57 CET 2009

      SUBROUTINE uve
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!**AVANCEMENT EXPLICITE D UN PAS DE TEMPS BAROCLINE DE U et V**
!  avec taux advection verticale implicite. + termes visqueux en 1/Rt.
!    et taux diffusion verticale implicite.
!  advection suivant la verticale :  w(k) * [ V(k-1) + V(k) ] / 2
!  combinee avec version traitement visco.barotrope type "s" (cf notes)
!  Suprime termes visqueux-metriques - Calcule Tension de Fond.
! Cfcc [Cfc0] => option : Avec Coeff Coriolis Reciproque = -2.Omega.cos(phi)
!  modif : 19/08/96

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
!! END_OF_USE_SECTION


!dmr [NOEQUI]      dimension phizu(imax,jmax,kmax+1), phizv(imax,jmax,kmax+1)
!dmr [NOEQUI]      equivalence ( phizu(1,1,1) , phizzz(1,1,1,1) )
!dmr [NOEQUI]      equivalence ( phizv(1,1,1) , phizzz(1,1,1,2) )

      real(kind=dblp), dimension(imax,jmax):: phiaxu, phiaxv, phidxu
     &               , phidxv, phiayu, phiayv, phidyu, phidyv, detadx
     &               , detady, cuhe ,  cvhe

!--variables locales :
      real(kind=dblp), dimension(imax) ::  wud2dz
      real(kind=dblp) :: ccdif, ccexa, ccexd, ccgdx, ccgdy, ccxdif
     &                 , ccydif, corfac, cu, cud, cv, cvd, u2d2
     &                 , v1d2, ww3
      integer(kind=ip) ::  i, j, k
      
      
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Forcage lie a l'elevation ; Calcul de la Tension de Fond.       |
!-----------------------------------------------------------------------

      ccgdx = gpes * uns2dx
      ccgdy = gpes * uns2dy
      do j=ju1,ju2
       do i=iu1(j),iu2(j)
        detadx(i,j) = ccgdx * smx(i,j,3) *
     &              ((eta(i-1,j-1)-eta(i,j))+(eta(i-1,j)-eta(i,j-1)))
        detady(i,j) = ccgdy * smy(i,j,3) *
     &              ((eta(i-1,j-1)-eta(i,j))-(eta(i-1,j)-eta(i,j-1)))
       enddo
      enddo

      if (cdbot.ne.zero) then
!-----
!--Tension de fond : TauX,Y = Cdrag * |V(bot.)| * u,v(bot.) = Coeff * u,v(bot.)
!  calcul du Coeff. (dans avudz(-,-,1))(pour resolution implicite).
      do j=ju1,ju2
       do i=iu1(j),iu2(j)
         avudz(i,j,1) = sqrt( u(i,j,kfu(i,j))*u(i,j,kfu(i,j))
     &                      + v(i,j,kfu(i,j))*v(i,j,kfu(i,j)) ) * cdbot
         phifu(i,j) = -avudz(i,j,1)*u(i,j,kfu(i,j))
         phifv(i,j) = -avudz(i,j,1)*v(i,j,kfu(i,j))
       enddo
      enddo
!-----
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Calcul des Flux Verticaux Explicite ADvection et DiFfusion.    |
!----------------------------------------------------------------------

      if (.not.(txiadu.eq.one .and. txidfu.eq.one)) then

      ccexa = (1. - txiadu) * 0.125
      ccexd =  1. - txidfu
!--Debut de la 1ere boucle externe sur l'indice de niveau k :
      do k=ku1+1,ku2
!-----

!--calcul des flux verticaux d'advection et de diffusion de qqmvt :
!      ccexak = ccexa * unsdzw(k)
       do j=ju1,ju2
!- construction vitesses verticales (a un fact.mult. pres ) :
        do i=iu1(j),iu2(j)
!         wud2dz(i)= ccexak
          wud2dz(i)= ccexa
     &          *(w(i,j,k)+w(i-1,j,k)+w(i,j-1,k)+w(i-1,j-1,k))
        enddo
!- calcul des flux :
        do i=iu1(j),iu2(j)
          ccdif =  avudz(i,j,k) * ccexd
          phizzz(i,j,k,1) = tmu(i,j,k-1) * (
!    &            (dz(k)*u(i,j,k-1)+dz(k-1)*u(i,j,k))*wud2dz(i)
     &            ( u(i,j,k-1) + u(i,j,k) ) * wud2dz(i)
     &          + ( u(i,j,k-1) - u(i,j,k) ) * ccdif )
          phizzz(i,j,k,2) = tmu(i,j,k-1) * (
!    &       (dz(k)*v(i,j,k-1)+dz(k-1)*v(i,j,k))*wud2dz(i)
     &       ( v(i,j,k-1) + v(i,j,k) ) * wud2dz(i)
     &          + ( v(i,j,k-1) - v(i,j,k) ) * ccdif )
        enddo
       enddo
!--fin du calcul des flux verticaux d'advection et de diffusion.

!--Fin de la 1ere boucle externe sur l'indice de niveau k .
      enddo

!--calcul des flux V. Explicites termine.
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      ww3 = 0.
      ccxdif = ahu * unsdx
      ccydif = ahu * unsdy
      corfac = 2.*(1.-txifcu)
!--Debut de la 2nd  boucle externe sur l'indice de niveau k :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,j,u2d2,ccdif,phiaxu,phiaxv,phidxu,phidxv)
!$DIR SPP LOOP_PRIVATE(v1d2,phiayu,phiayv,phidyu,phidyv)
!$DIR SPP LOOP_PRIVATE(cud,cvd,cuhe,cvhe,cu,cv,ww3)
      do  k=ku2,ku1,-1

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Resolution Complete des Flux Horizontaux .        |
!---------------------------------------------------------

!**place reservee pour le calcul des taux de decentrement**

!--calcul des flux (advection "a", diffusion "d") dans les directions x et y :

      do j=ju1,ju2
       do i=iu1(j)-1,iu2(j)
        u2d2 = cmy(i,j,2) * 0.25 * (u(i,j,k) + u(i+1,j,k))
        ccdif = smxy(i,j,2) * ccxdif
        phiaxu(i,j) = u2d2  * (u(i,j,k)+u(i+1,j,k))
     &              + alphxu * abs(u2d2) * (u(i,j,k)-u(i+1,j,k))
        phiaxv(i,j) = u2d2  * (v(i,j,k)+v(i+1,j,k))
     &              + alphxv * abs(u2d2) * (v(i,j,k)-v(i+1,j,k))
        phidxu(i,j) = ccdif * (u(i,j,k)-u(i+1,j,k))
        phidxv(i,j) = ccdif * (v(i,j,k)-v(i+1,j,k))
        enddo
      enddo

      do  j=ju1-1,ju2
       do  i=iuf1(j),iuf2(j)
        v1d2 = cmx(i,j,1) * 0.25 * (v(i,j,k)+v(i,j+1,k))
        ccdif = cmxy(i,j,1) * ccydif
        phiayu(i,j) = v1d2  * (u(i,j,k)+u(i,j+1,k))
     &              + alphyu * abs(v1d2) * (u(i,j,k)-u(i,j+1,k))
        phiayv(i,j) = v1d2  * (v(i,j+1,k)+v(i,j,k))
     &              + alphyv * abs(v1d2) * (v(i,j,k)-v(i,j+1,k))
        phidyu(i,j) = ccdif * (u(i,j,k)-u(i,j+1,k))
        phidyv(i,j) = ccdif * (v(i,j,k)-v(i,j+1,k))
        enddo
       enddo

!--fin de calcul des flux horizontaux.

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Bilan de tous les Flux calcules Explicitement .   |
!---------------------------------------------------------

      do  j=ju1,ju2
       do  i=iu1(j),iu2(j)
        cud = smxy(i,j,3)*( unsdx*(phidxu(i-1,j)-phidxu(i,j))
     &                    + unsdy*(phidyu(i,j-1)-phidyu(i,j)) )
        cvd = smxy(i,j,3)*( unsdx*(phidxv(i-1,j)-phidxv(i,j))
     &                    + unsdy*(phidyv(i,j-1)-phidyv(i,j)) )

!- test 2nd terme corriolis : du/dt = - 2 . Omega . cos(Phi) . w
!fcc    ww3 = (w(i-1,j-1,k) + w(i,j,k)) + (w(i-1,j,k) + w(i,j-1,k))
!fcc &    + (w(i-1,j-1,k+1)+w(i,j,k+1)) + (w(i-1,j,k+1)+w(i,j-1,k+1))
        cuhe(i,j) = tmu(i,j,k) * ( cud
!fcc &      + fcucor(i,j) * ww3
     &      + smx(i,j,3)*uns2dx
     &          *((q(i-1,j-1,k)-q(i,j,k))+(q(i-1,j,k)-q(i,j-1,k)))
     &      + smxy(i,j,3)*( unsdx*(phiaxu(i-1,j)-phiaxu(i,j))
     &                    + unsdy*(phiayu(i,j-1)-phiayu(i,j))
     &                    - cmxdy(i,j)*u(i,j,k)*v(i,j,k)
     &                    + cmydx(i,j)*v(i,j,k)*v(i,j,k)     ))
        cvhe(i,j) = tmu(i,j,k) * ( cvd
!fcc &      + fcvcor(i,j) * ww3
     &      + smy(i,j,3)*uns2dy
     &          *((q(i-1,j-1,k)-q(i,j,k))-(q(i-1,j,k)-q(i,j-1,k)))
     &      + smxy(i,j,3)*( unsdx*(phiaxv(i-1,j)-phiaxv(i,j))
     &                    + unsdy*(phiayv(i,j-1)-phiayv(i,j))
     &                    - cmydx(i,j)*u(i,j,k)*v(i,j,k)
     &                    + cmxdy(i,j)*u(i,j,k)*u(i,j,k)     ))

        fub(i,j,k) = dz(k)*(cuhe(i,j)-cud)
        fvb(i,j,k) = dz(k)*(cvhe(i,j)-cvd)
        enddo
       enddo

      do  j=ju1,ju2
       do  i=iu1(j),iu2(j)
        cu = cuhe(i,j) + detadx(i,j) + corfac*fs2cor(i,j)*v(i,j,k)
     &          + unsdz(k)*(phizzz(i,j,k,1)-phizzz(i,j,k+1,1))
        cv = cvhe(i,j) + detady(i,j) - corfac*fs2cor(i,j)*u(i,j,k)
     &          + unsdz(k)*(phizzz(i,j,k,2)-phizzz(i,j,k+1,2))

        u(i,j,k) = u(i,j,k) + dtu * cu
        v(i,j,k) = v(i,j,k) + dtu * cv
        enddo
       enddo

!--Fin de la 2nd  boucle externe sur l'indice de niveau k .
      enddo

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine uve -
      end
