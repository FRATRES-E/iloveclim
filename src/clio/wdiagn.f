!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:57 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:57 CET 2009

      SUBROUTINE wdiagn
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- calcul de la vitesse verticale w depuis la surface jusqu'au fond.
!  avec compensation de la divergence barotrope -> dans w(-,-,1).
!  modif : 09/01/97

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
!! END_OF_USE_SECTION


!--variables locales equivalentes (pour raison de place memoire) :
      real(kind=dblp), dimension(imax):: phixj

      integer(kind=ip) :: i, j, k
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      write(66,*) 'wdiagn : W sans DIV(ub,vb)'

!--Debut de la boucle externe sur l'indice de latitude j :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,k,phixj,dzd2)
      do 400 j=js1,js2
!-----

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) calcul 2D ; Bilan de l'Eq. Continuite.                          |
!-----------------------------------------------------------------------

!--initialisation :
      do 100 i=1,imax
!       w(i,j,ks1) = 0.0
        w(i,j,ks2+1) = 0.0
 100  continue

!--Calcul de la divergence du transport barotrope (ub,vb) :
      do 130 i=is1(j),is2(j)+1
        phixj(i) = cmy(i,j,1) * (ub(i,j+1)+ub(i,j))
 130  continue
      do 150 i=is1(j),is2(j)
        w(i,j,1) = smxy(i,j,0) *
     &           ( uns2dx * (phixj(i)-phixj(i+1))
     &           + uns2dy * ( cmx(i,j  ,2)*(vb(i+1,j  )+vb(i,j  ))
     &                      - cmx(i,j+1,2)*(vb(i+1,j+1)+vb(i,j+1)) ))
 150  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Bilan de l'Eq. Continuite et Somme directement dans w .         |
!-----------------------------------------------------------------------

!--boucle sur l'indice de niveau, "k" :
      do 300 k=ks2,1+ks1,-1

      do 230 i=is1(j),is2(j)+1
        phixj(i) = cmy(i,j,1) * (u(i,j+1,k)+u(i,j,k))
 230  continue
!       phiy(i,j) = cmx(i,j,2) * (v(i+1,j,k)+v(i,j,k))

!     dzd2 = 0.5 * dz(k)
      do 250 i=is1(j),is2(j)
!       w(i,j,k) = w(i,j,k+1) + dz(k) * smxy(i,j,0) *
        w(i,j,k) = w(i,j,k+1) + dz(k) * ( smxy(i,j,0) *
!    &           ( uns2dx * (phix(i+1,j)-phix(i,j))
!    &           + uns2dy * (phiy(i,j+1)-phiy(i,j)) )
     &           ( uns2dx * (phixj(i+1)-phixj(i))
     &           + uns2dy * ( cmx(i,j+1,2)*(v(i+1,j+1,k)+v(i,j+1,k))
     &                      - cmx(i, j, 2)*(v(i+1, j, k)+v(i, j, k)) ))
     &        + w(i,j,1) / zw(kfs(i,j)) )
!       w(i,j,k+1) = w(i,j,k+1) * tms(i,j,k)
        w(i,j,k) = w(i,j,k) * tms(i,j,k)
 250  continue

!-----
 300  continue
!--fin boucle sur "k".

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Fin de la boucle externe sur l'indice de latitude j .
 400  continue

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine wdiagn -
      end
