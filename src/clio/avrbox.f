!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009

      SUBROUTINE avrbox(var, vinbox, volbox, coef,
     &       ibond1,ibond2, jlim1,jlim2, nboxmx,nbox, km, ns, nn99)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! calcule (dans vinbox) l'integrale sur chaque boite de la variable "var"
! km > 0 : Integrale et var 3 D. -- km = 0 : Integrale et var 2 D.
! N.B: Sans initialisation de vinbox (=> doit etre fait avant)
!      calcule du volume d'integration si ns = 0
!  modif : 17/08/97

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use newunit_clio_mod, only: mouchard_id

!! END_OF_USE_SECTION


!--dummy variables :

      integer :: nboxmx, km

      real(kind=dblp), dimension(imax,jmax,*) :: var
      real(kind=dblp), dimension(nboxmx,0:km) :: vinbox, volbox
      real(kind=dblp), dimension(jmax,*)      :: ibond1, ibond2
      real(kind=dblp), dimension(*)           :: coef, jlim1, jlim2     
                  
c~       dimension var(imax,jmax,*)
c~       dimension vinbox(nboxmx,0:km), volbox(nboxmx,0:km), coef(*)
c~       dimension ibond1(jmax,*), ibond2(jmax,*), jlim1(*), jlim2(*)


c~       real(kind=dblp) :: 
      integer(kind=ip):: nbox, ns, nn99, nb, k, j, i

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initialisation :
!-----------------------------------------------------------------------

      if (ns.eq.0) then
!- initialisation :
        do nb=1,nboxmx
         do k=0,km
          volbox(nb,k) = 0.0
         enddo
        enddo
      endif

      if (km.eq.0) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Integrale sur chacune des boites "nb" :
!-----------------------------------------------------------------------
      do nb=1,nbox
!-----

      if (ns.eq.0) then
!- Calcul des volumes :
        do j=jlim1(nb),jlim2(nb)
         do i=ibond1(j,nb),ibond2(j,nb)
           volbox(nb,0) = volbox(nb,0)
     &                  + coef(1) * ctmi(i,j,ks2,0)
         enddo
        enddo
      endif

      do j=jlim1(nb),jlim2(nb)
       do i=ibond1(j,nb),ibond2(j,nb)
        vinbox(nb,0) = vinbox(nb,0)
     &               + coef(1) * ctmi(i,j,ks2,0) * var(i,j,1)
       enddo
      enddo

!--fin de la boucle sur la boite "nb".
      enddo
!-----

      else
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Integrale sur chacun niveaux de chacune des boites "nb" :
!-----------------------------------------------------------------------
      do nb=1,nbox
       do k=1,km
!-----
      if (ns.eq.0) then
!- Calcul des volumes :
        do j=jlim1(nb),jlim2(nb)
         do i=ibond1(j,nb),ibond2(j,nb)
           volbox(nb,k) = volbox(nb,k)
     &                  + coef(k) * ctmi(i,j,k,0)
         enddo
        enddo
        volbox(nb,0) = volbox(nb,0) + volbox(nb,k)
      endif
!-----

      do j=jlim1(nb),jlim2(nb)
       do i=ibond1(j,nb),ibond2(j,nb)
        vinbox(nb,k) = vinbox(nb,k)
     &               + coef(k) * ctmi(i,j,k,0) * var(i,j,k)
        enddo
       enddo
      vinbox(nb,0) = vinbox(nb,0) + vinbox(nb,k)

!--fin des boucles sur la boite "nb" et le niveau "k".
       enddo
      enddo
!-----
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Impresion sur fichier "mouchard" :
      if (nn99.eq.2 .and. ns.eq.0) then
          write(mouchard_id,*)
        do nb=1,nbox
          write(mouchard_id,'(A,I3)') ' avrbox,'//
     &      ' Volume (m) [/dx*dy , surf -> fond] de la boite nb =', nb
          write(mouchard_id,'(1P5E12.5)') (volbox(nb,k),k=km,0,-1)
          write(mouchard_id,*)
        enddo
      endif

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine avrbox -
      end
