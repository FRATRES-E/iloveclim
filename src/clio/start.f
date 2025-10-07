!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:54 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:54 CET 2009

      SUBROUTINE start
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  demarrage a partir de l'etat de repos :
!  modif : 24/09/99

!! START_OF_USE_SECTION

      use const_mod,only:

      use para0_mod,only: imax, jmax, nsmax
      use para_mod, only:
      use bloc0_mod,only: ks1, ks2, tpstot, scal, scal0, u, v, eta, ub, vb
      use bloc_mod, only: numit
      use ice_mod,  only: tbq, tfu, hgbq, albq, hnbq, ts, firg, fcsg
     &           , fleg, fsbbq, qstobq, xzo

      use global_constants_mod, only: dblp=>dp,ip

      use newunit_clio_mod, only: clio3_out_id
!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "ice.com"

!! END_OF_INCLUDE_SECTION

!--variables locales :

      integer(ip) :: i, j, k, n

      write(clio3_out_id,*) 'debut de start'

      tpstot = 0.0
      numit = 0

      do n=1,nsmax
       do k=ks1,ks2
        do j=1,jmax
         do i=1,imax
          scal(i,j,k,n) = scal0(k,n)
         enddo
        enddo
       enddo
      enddo

      do k=ks1,ks2
       do j=1,jmax
        do i=1,imax
         u(i,j,k) = 0.0
         v(i,j,k) = 0.0
        enddo
       enddo
      enddo

      do j=1,jmax
       do i=1,imax
        eta(i,j) = 0.0
        ub(i,j) = 0.0
        vb(i,j) = 0.0
       enddo
      enddo 

      do j=1,jmax
        do i=1,imax
!
!                        tfu: MELTING POINT OF SEA WATER.
!
            tfu(i,j)    = abs(273.15-0.0575*scal(i,j,ks2,2)+
     &                        1.710523e-03*sqrt(scal(i,j,ks2,2))**3-
     &                        2.154996e-04*scal(i,j,ks2,2)**2)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
            hgbq(i,j)   = 0.0
            albq(i,j)   = 1.0
            hnbq(i,j)   = 0.0
            ts(i,j)     = tfu(i,j)
            tbq(i,j,1)  = tfu(i,j)
            tbq(i,j,2)  = tfu(i,j)
            tbq(i,j,3)  = tfu(i,j)
            firg(i,j)   = 0.0
            fcsg(i,j)   = 0.0
            fleg(i,j)   = 0.0
            fsbbq(i,j)  = 0.0
            qstobq(i,j) = 0.0
            xzo(i,j)    = 0.001
        enddo
      enddo

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine start -
      end
