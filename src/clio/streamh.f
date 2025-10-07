!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:55 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:55 CET 2009

      SUBROUTINE streamh(psi)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  calcule la fonction de courrant de transport horizontal (psi)
!  a partir du Flux Zonal integre sur la verticale Ub,
!   en commencant par le Nord vers le Sud.
!-- Methode : comme sur une grille C , psi au coins .
!  modif : 28/07/94

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use reper_mod
!! END_OF_USE_SECTION



!--dummy argument :
      real(kind=dblp), dimension(imax,jmax) :: psi

!--variables locale :
      integer(kind=ip), dimension(imax):: jnordb

!--- more locales
      real(kind=dblp) :: dxd2, dyd2, ps0, psd, psh
      integer(kind=ip):: i, ii, j, jj


!--initialisation :
      do j=1,jmax
       do i=1,imax
        psi(i,j) = 0.
       enddo
      enddo

!--Bordure Nord : pts ou la composante U = 0
      do i=1,imax
        jnordb(i) = jmax
        if (jsep(i).eq.jsepar+1) jnordb(i) = jsepar
      enddo

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  3 ) Calcul de la fonction de courant :
!-----------------------------------------------------------------------

!--valeur de la fonction de courrant (pts 3) sur la borbure Nord :
      jj = jnordb(1)
      psi(1,jj) = 0.
      psh = 0.
      dxd2 = 0.5 * dx
      do i=2,imax
!- calcul du flux a la frontiere Nord (integre vb.dx ) :
        jj = jnordb(i)
        ii = i - 1
        ps0 = dxd2 * vb(i,jj)
        psd = tms(ii,jj-1,ks2) * tms(ii,jj,ks2) * cmx(ii,jj,2) * ps0
        psh = psh + psd
        psi(i,jj) = psh
        psd = tms(i,jj-1,ks2)  * tms(i,jj,ks2)  * cmx(i,jj,2)  * ps0
        psh = psh + psd
      enddo

!--calcul de la fonction de courrant (pts 3)
!   en integrant (-ub1).(-dy) a partir de la vitesse ub au pts 1 :
      dyd2 = 0.5 * dy
      do i=ims1,ims2
       do j=jnordb(i),ju1,-1
          psi(i,j-1) = psi(i,j)
     &               + cmy(i,j-1,1) * dyd2 * (ub(i,j) + ub(i,j-1))
       enddo
      enddo

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  4 ) Raccord cyclique / entre grille .
!-----------------------------------------------------------------------

      if (ltest.ge.1) then
!--raccord cyclique pour psi :
      do j=jcl1,jcl2
        psi(ims1-1,j) = psi(ims2,j)
        psi(ims2+1,j) = psi(ims1,j)
      enddo
      endif

      return
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!- fin de la routine streamh -
      end
