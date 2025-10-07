!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009

      SUBROUTINE iso_copy
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  modif : 07/03/99

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use datadc_mod

      use isoslope_mod
      use reper_mod
      use varno_mod

      use newunit_clio_mod, only: clio3_out_id

!! END_OF_USE_SECTION


!--- locales
      integer(kind=ip):: i, j, k
      real(kind=dblp) :: avi, avs


      write(clio3_out_id,'(A)')
     &      ' *** ISO_COPY: ppmodif: 20-03-97: gm90,class ***'

!- Copy Variables
      do k=1,kmax
       do j=1,jmax
        do i=1,imax
         avi=ai(k)*(c4x(i,j,k)**2+c4y(i,j,k)**2)+epsil
         avs=avsdz(i,j,k)+epsil
         avudz(i,j,k)=LOG10(avi)
         avsdz(i,j,k)=LOG10(avs)
         b(i,j,k)=c4x(i,j,k)
         bvf(i,j,k)=c4y(i,j,k)
        enddo
       enddo
      enddo
      call raccord(avsdz(1,1,1),zero,kmax,4)
      call raccord(avudz(1,1,1),zero,kmax,4)
      call raccord(b(1,1,1),zero,kmax,4)
      call raccord(bvf(1,1,1),zero,kmax,4)
      nvrl(nvravi) = 3
      nvrl(nvravs) = 3
      nvrl(nvrslx) = 3
      nvrl(nvrsly) = 3
      nvrl(nvras)= 0
      nvrl(nvrau)= 0
      nvrl(nvrb) = 0
      nvrl(nvrn2)= 0

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine iso_copy -
      end
