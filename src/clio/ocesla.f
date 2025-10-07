!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:50 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:50 CET 2009

      SUBROUTINE ocesla
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  Programme simulant un ocean statique agissant uniquement comme un reservoir
!  de chaleur. Ce programme n'est utilise que pour la mise au point!
!  modif : 01/04/94

!! START_OF_USE_SECTION

      use const_mod, only:

      use para0_mod, only: nsmax
      use para_mod , only:
      use bloc0_mod, only: scal, is1, is2, js1, js2, ks2, phiss
      use bloc_mod,  only:
      use global_constants_mod, only: dblp=>dp, ip
      use newunit_clio_mod, only: clio3_out_id
      
!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"

!! END_OF_INCLUDE_SECTION


! --- dmr Local variables arising from implicit none construct

      integer(ip) :: i,j,k

!        write(clio3_out_id,*)'avant',scal(imax-2,5,ks2,1),scal(imax-2,5,ks2,2)
!      &                    ,phiss(imax-2,5,1),phiss(imax-2,5,2)
      do 30 k=1,nsmax
         do 20 j=js1,js2
            do 10 i=is1(j),is2(j)
               scal(i,j,ks2,k)=scal(i,j,ks2,k) - phiss(i,j,k)
!              scal(i,j,ks2,k)=scal(i,j,ks2,k)
!    &                         -dts(ks2)*unsdz(ks2)*phiss(i,j,k)
 10         continue
 20      continue
 30   continue
!        write(clio3_out_id,*)'temp apres',scal(imax-2,5,ks2,1),
!    &              scal(imax-2,5,ks2,2)
      return
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!- fin de la routine slaboc -
      end
