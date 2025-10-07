!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:46 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:46 CET 2009

      SUBROUTINE flowsurf(titub, titvb)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  calcule le tranport horizontal dans la couche de surface :
!  remplace ub & vb (=transport total) par le transport partiel (surface).
!  modif : 29/12/94

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod

      use newunit_clio_mod, only: clio3_out_id
!! END_OF_USE_SECTION

!--dummy variables :
      character*(*) titub, titvb

!--local variables :
      character*19 cctit

      integer(kind=ip) :: kbase, lltit, i, j, k


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Definition des niveaux (surface) concernes .
!-----------------------------------------------------------------------

      write(clio3_out_id,*) 'Transport Surfacique (k2=kmax) : Input k1  ?'
      read(5,*)  kbase
      kbase = max(1,min(kbase,ks2))

      write(cctit,'(A3,I2,A14)') ' k=', kbase, ' -> Surf. Flow'
      lltit = min(6+len(cctit), len(titub))
      titub(7:lltit) = cctit
      lltit = min(6+len(cctit), len(titvb))
      titvb(7:lltit) = cctit

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Somme
!-----------------------------------------------------------------------

!--initialisation :
      do j=1,jmax
       do i=1,imax
        ub(i,j) = 0.
        vb(i,j) = 0.
       end do
      end do

      do k=kbase,ks2
       do j=1,jmax
        do i=1,imax
         ub(i,j) = ub(i,j) + dz(k) * u(i,j,k)
         vb(i,j) = vb(i,j) + dz(k) * v(i,j,k)
        end do
       end do
      end do

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      write(clio3_out_id,'(2A,I2,A,F8.2)') 'Ub,Vb : u,v integres ',
     &     'depuis z(kmax)=0 jusque z(', kbase, ') =', zw(kbase)

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine flowsurf -
      end
