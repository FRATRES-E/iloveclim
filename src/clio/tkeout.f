!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:56 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:56 CET 2009

      SUBROUTINE tkeout(avuk,avsk)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
! transfert des Viscosites & Diffusivites Verticales depuis TKE...
!   vers les points wu & ws avec division par dz.
!  modif : 04/10/95

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
!! END_OF_USE_SECTION

!--dummmy variables :
      real(kind=dblp), dimension(imax,jmax,kmax) :: avuk, avsk
      
!--variables locales :
      real(kind=dblp) :: ccdz
      integer(kind=ip):: i, j, k


!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  1 ) Transfert des Diffusivites V.
!-----------------------------------------------------------------------

!--transfert (avec division par dzw) :
      do k=ks1+1,ks2
       do j=js1,js2
        do i=is1(j),is2(j)
         avsdz(i,j,k) = avsk(i,j,k) * unsdzw(k)
        enddo
       enddo
      enddo
      
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  2 )  Transfert des Viscosites V.
!-----------------------------------------------------------------------

!--raccord cyclique et autre (avuk) :
      call raccord(avuk(1,1,2), zero, kmax-1, 4)

      do k=ks1+1,ks2

!--moyenne sur 4 points (et division par dzw) :
       ccdz = 0.25 * unsdzw(k)
       do j=ju1,ju2
        do i=iu1(j),iu2(j)
         avudz(i,j,k) = ccdz *
     &   (avuk(i-1,j-1,k)+avuk(i,j,k)+avuk(i-1,j,k)+avuk(i,j-1,k))
        enddo
       enddo
       
      enddo
      return
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!- fin de la routine tkeout -
      end
