!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:47 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:47 CET 2009

      SUBROUTINE gm_flux(ilat)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- computation of the meridional heat and salt transport
!- !! we must have do this before modifying scal(i,j,k)=S or T

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip           

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use isoslope_mod
      use reper_mod
      use varno_mod
      
      use newunit_clio_mod, only: clio3_out_id      
      
!! END_OF_USE_SECTION


      real(dblp), dimension(imax,jmax,kmax) :: phix, phiy

      integer(kind=ip):: ns, ilat, k, j, i
            
      real(kind=dblp) :: xmult, u1, t1, v2, t2

      
      data ns /1/
      data xmult /4.1D-9/
      character*10 name

      write(clio3_out_id,'(2A,F5.0)') ' *** GM_FLUX: ppmodif: 14-04-97:',
     &   ' gm90,class *** : AITD=', aitd(kmax)

      call isoadvec

!-************************************************
!- GM90 ADVECTION FLUX (centered advection)
!-************************************************

      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            phix(i,j,k)=zero
            phiy(i,j,k)=zero
          enddo
        enddo
      enddo

      do k=1,kmax
        do j=js1,js2
          do i=ims1,ims2
            u1=uiso(i,j,k)
            t1=0.5*( scal(i,j,k,ns)+scal(i-1,j,k,ns) )
            phix(i,j,k)=ttm1(i,j,k)*cmy(i,j,1)*u1*t1
          enddo
        enddo
      enddo

      do k=1,kmax
        do j=js1,1+js2
          do i=ims1,ims2
            v2=viso(i,j,k)
            t2=0.5*( scal(i,j,k,ns)+scal(i,j-1,k,ns) )
            phiy(i,j,k)=ttm2(i,j,k)*cmx(i,j,2)*v2*t2
          enddo
        enddo
      enddo

      name='idl.gm'
      call iso_aver(ns,ilat,name,phix,phiy,xmult)

!-************************************************
!- OVERTURNING ADVECTION FLUX (center scheme)
!-************************************************
      do k=1,kmax
        do j=1,jmax
          do i=1,imax
            phix(i,j,k)=zero
            phiy(i,j,k)=zero
          enddo
        enddo
      enddo

      do k=1,kmax
        do j=js1,js2
          do i=ims1,ims2
            u1=uiso(i,j,k)
            phix(i,j,k)=ttm1(i,j,k)*cmy(i,j,1)*u1
          enddo
        enddo
      enddo

      do k=1,kmax
        do j=js1,1+js2
          do i=ims1,ims2
            v2=viso(i,j,k)
            phiy(i,j,k)=ttm2(i,j,k)*cmx(i,j,2)*v2
          enddo
        enddo
      enddo

      name='idl.overgm'
      call iso_aver(ns,ilat,name,phix,phiy,xmult)

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      end
