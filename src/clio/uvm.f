!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:57 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:57 CET 2009

      SUBROUTINE uvm
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Resolution des termes verticaux de l'eq.qqmvt par appel a la routine "uvi".
!**AJUSTEMENT DE LA MOYENNE DES VITESSES HORIZONTALES**
!  modif : 08/01/96

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
!! END_OF_USE_SECTION


!--variables locales :
      real(kind=dblp), dimension(imax) :: ui, vi
      integer(kind=ip):: i, j, k
      real(kind=dblp) :: unsncl, vvber
      
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Calcul des moyennes de eta,ub,vb-spl , remplace eta, ub, vb :   |
!-----------------------------------------------------------------------

       j = js1
       do i=is1(j),is2(j)
        eta(i,j) = etaspl(i,j)
       enddo

!--Debut de la boucle externe sur l'indice de latitude j :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,k,ui,vi)
      jloop: do j=ju1,ju2
!-----

       do i=is1(j),is2(j)
        eta(i,j) = etaspl(i,j)
       enddo

       do i=iu1(j),iu2(j)
        ub(i,j) = ubspl(i,j)
        vb(i,j) = vbspl(i,j)
       enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Resolution Implicite des termes verticaux de l'eq. de qqmvt .   |
!-----------------------------------------------------------------------
!REFACTORING: the next line probably to be removed: already inside a "do j=ju1,ju2" loop.
!     do 280 j=ju1,ju2
        call uvi(j)
!REFACTORING: accordingly, the next line is redundant
 280  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Ajuste la vitesse Horizontale sur sa moyenne sur la Verticale . |
!-----------------------------------------------------------------------

      do i=iu1(j),iu2(j)
       ui(i)=0.0
       vi(i)=0.0
      enddo

      do k=ku1,ku2
       do i=iu1(j),iu2(j)
        ui(i) = ui(i) + dz(k)*u(i,j,k)
        vi(i) = vi(i) + dz(k)*v(i,j,k)
       enddo
      enddo
      do i=iu1(j),iu2(j)
        ui(i) = unshu(i,j)*(ub(i,j)-ui(i))
        vi(i) = unshu(i,j)*(vb(i,j)-vi(i))
      enddo
      do k=ku1,ku2
       do i=iu1(j),iu2(j)
        u(i,j,k) = tmu(i,j,k) * ( u(i,j,k) + ui(i) )
        v(i,j,k) = tmu(i,j,k) * ( v(i,j,k) + vi(i) )
       enddo
      enddo

!-----

      enddo jloop

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  6 ) raccord cyclique + raccord de grille pour eta, ub, vb, u, v .   |
!-----------------------------------------------------------------------

      if (ltest.eq.3) then
!--vitesse barocline nulle :
        vvber = unshu(iberp,jberp) * vb(iberp,jberp)
        do k=ks1,ks2
          u(iberp,jberp,k) = 0.
          v(iberp,jberp,k) =  tmu(iberp,jberp,k) * vvber
        enddo
      endif

      call raccord( eta, zero, 1, 0)
      call raccord( ub, zero, 1, 23)
      call raccord( vb, zero, 1, 35)
      call raccord( u(1,1,1), zero, kmax, 15)
      call raccord( v(1,1,1), zero, kmax, 27)

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  7 ) Computation of the average of the surface velocity              |
!-----------------------------------------------------------------------

      unsncl = one / DFLOAT(nclin)
      do j=ju1,ju2
        do i=iu1(j),iu2(j)
          umoy(i,j)=umoy(i,j)+u(i,j,ku2)*unsncl
          vmoy(i,j)=vmoy(i,j)+v(i,j,ku2)*unsncl
        enddo
      enddo

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine uvm -
      end
