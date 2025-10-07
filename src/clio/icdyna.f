!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009

      SUBROUTINE icdyna
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  This routine calls ice dynamic and advcect sea ice properties.
!---
! Ccpl [Ccp0] => ligne specifique a la version avec [sans] couplage .
!---
!  modif : 28/09/99



      use para0_mod, only: imax, jmax
      use bloc0_mod, only: phiss, phivs, tmu, phisu, phisv, ust2s
     >             , is1, is2, scal, u, v, iu1, iu2, js1, js2, ks2, ks1
     >             , ku1, ku2, ju1, ju2, jeq
      use bloc_mod,  only: idyn
      use ice_mod,   only: tfu, alcd, alcr, hgbqp, hnbq, albq, hgbq
     >             , tairox, tairoy, sdvt
      use dynami_mod,only: uo, vo, hnm, hgm, ipo1i0, jpo1i0, ug, vg
     >             , npo1i0, sangvg, rhoco, cangvg
      use const_mod, only: zero, one, rho0



      implicit none

      real*8, dimension(5) ::  trest
      
!nb debut
      real*8 fluxbrines(imax,jmax)
      common / common_brines / fluxbrines
!nb fin

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! dmr ---     Local variables
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      integer :: i, j, k, n, ih

      real*8  :: sang, tairx, tairy, zmod, tglx, tgly, t11, t12, t21, t22
     &         , tmoy, tot, zustm


!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--1. Dynamics of sea ice-initialisation                               |
!-----------------------------------------------------------------------
!
      do j=js1,js2
        do i=is1(j),is2(j)
          tfu(i,j)   = abs(273.15-
     &                         0.0575*scal(i,j,ks2,2)+
     &                         1.710523e-3*sqrt(scal(i,j,ks2,2))**3-
     &                         2.154996e-4*scal(i,j,ks2,2)**2)

#if ( ISOOCN >= 1 )
          phiss(i,j,:) = 0.0 ! phiss with isoocn is on owatert size
#else
          phiss(i,j,0) = 0.0
          phiss(i,j,1) = 0.0
          phiss(i,j,2) = 0.0
#endif
!nb debut
!          fbrines(i,j) = 0.0
          fluxbrines(i,j)=0.0
!nb fin
          hgbqp(i,j)   = 0.0

        enddo
      enddo


       do k=ks1,ks2
        do j=js1,js2
          do i=is1(j),is2(j)
            phivs(i,j,k,1) = 0.0
          enddo
        enddo
       enddo


      if (idyn.eq.1) then
!
!     initialise array alcd for tracking rate of lead closure/opening
!     due to convervenge/divergence.
!     initialise array alcr for extra amount of lead opening/closing
!     due to shearing deformation.
!
      do j=1,jmax
        do i=1,imax
          alcd(i,j) = zero
          alcr(i,j) = zero
        enddo
      enddo
!
!--1.1. Mean ice and snow thicknesses.
!--------------------------------------
!
        do j=1,jmax
          do i=1,imax
            hnm(i,j) = (1.0-albq(i,j))*hnbq(i,j)
            hgm(i,j) = (1.0-albq(i,j))*hgbq(i,j)
!           uo(i,j)  = uost(i,j)*tmu(i,j,ku2)
!           vo(i,j)  = vost(i,j)*tmu(i,j,ku2)
            uo(i,j)  = u(i,j,ku2)*tmu(i,j,ku2)
            vo(i,j)  = v(i,j,ku2)*tmu(i,j,ku2)
!           uo(i,j)  = umoy(i,j)*tmu(i,j,ku2)
!           vo(i,j)  = vmoy(i,j)*tmu(i,j,ku2)
          enddo
        enddo 
!
!--1.2. Call to dynamics routine.
!---------------------------------
!
!-no-ice velocity point
        do n=1,npo1i0
         trest(n)=tmu(ipo1i0(n),jpo1i0(n),ks2)
         tmu(ipo1i0(n),jpo1i0(n),ks2)= 0.0
        enddo
!
!  Northern hemisphere.
!
! Cdzh dynamics according to Zhang and Hibler, 1997.
! Cdzr dynamics according to Zhang and Rothrock, 2000.
!
        call dynami_zh(+1)
!dzr    call dynami_zr(+1)
!
!  Southern hemisphere.
!
        call dynami_zh(-1)
!dzr    call dynami_zr(-1)
!
!-no-ice velocity point
        do n=1,npo1i0
         tmu(ipo1i0(n),jpo1i0(n),ks2)= trest(n)
        enddo
!
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--2. Computation of flux to the ocean.                                |
!-----------------------------------------------------------------------
!
        do j=ju1,ju2
           ih=sign(1,j-jeq)
           sang = real(ih)*sangvg
           do i=iu1(j),iu2(j)
!cp1        if (icoupl .ne. 0) then
!cp1         tairx = albq(i,j)*tairox(i,j)
!cp1 &              +albq(i-1,j)*tairox(i-1,j)
!cp1 &              +albq(i-1,j-1)*tairox(i-1,j-1)
!cp1 &              +albq(i,j-1)*tairox(i,j-1)
!
!cp1         tairy = albq(i,j)*tairoy(i,j)
!cp1 &              +albq(i-1,j)*tairoy(i-1,j)
!cp1 &              +albq(i-1,j-1)*tairoy(i-1,j-1)
!cp1 &              +albq(i,j-1)*tairoy(i,j-1)
!cp1        else
             tairx = (albq(i,j)+albq(i-1,j)
     &               +albq(i-1,j-1)+albq(i,j-1))*tairox(i,j)
!
             tairy = (albq(i,j)+albq(i-1,j)
     &               +albq(i-1,j-1)+albq(i,j-1))*tairoy(i,j)
!cp1        endif
             zmod  = sqrt((ug(i,j)-uo(i,j))**2+(vg(i,j)-vo(i,j))**2)
             tglx  = (4-albq(i,j)-albq(i-1,j)-albq(i-1,j-1)-albq(i,j-1))
!    &               *rhoco*zmod*(ug(i,j)-uo(i,j))
     &               *rhoco*zmod*(cangvg*(ug(i,j)-uo(i,j))-
     &                            sang*(vg(i,j)-vo(i,j)))
             tgly  = (4-albq(i,j)-albq(i-1,j)-albq(i-1,j-1)-albq(i,j-1))
!    &               *rhoco*zmod*(vg(i,j)-vo(i,j))
     &               *rhoco*zmod*(cangvg*(vg(i,j)-vo(i,j))+
     &                            sang*(ug(i,j)-uo(i,j)))
             phisu(i,j) = -(tairx+1.0*tglx)/(4*rho0)
             phisv(i,j) = -(tairy+1.0*tgly)/(4*rho0)
           enddo
         enddo  
!
        do j=js1,js2
           do i=is1(j),is2(j)
              t11=rhoco*((ug(i-1,j-1)-uo(i-1,j-1))**2
     &                  +(vg(i-1,j-1)-vo(i-1,j-1))**2)
              t12=rhoco*((ug(i-1,j)-uo(i-1,j))**2
     &                  +(vg(i-1,j)-vo(i-1,j))**2)
              t21=rhoco*((ug(i,j-1)-uo(i,j-1))**2
     &                  +(vg(i,j-1)-vo(i,j-1))**2)
              t22=rhoco*((ug(i,j)-uo(i,j))**2
     &                  +(vg(i,j)-vo(i,j))**2)
              tmoy=0.25*(t11+t12+t21+t22)
              tot=tmoy*(1-albq(i,j))
     &            +albq(i,j)*sqrt(tairox(i,j)**2+tairoy(i,j)**2)
              zustm=tot/rho0
              ust2s(i,j)=zustm*(one+sdvt(i,j))
            enddo
          enddo
      else
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--3. If no ice dynamics                                               |
!-----------------------------------------------------------------------
!
        do j=ju1,ju2
          do i=iu1(j),iu2(j)
!cp1        if (icoupl .ne. 0) then
!cp1           phisu(i,j) = -(tairox(i,j) + tairox(i-1,j)
!cp1 &                      + tairox(i-1,j-1) + tairox(i,j-1))/(rho0*4)
!cp1           phisv(i,j) = -(tairoy(i,j) + tairoy(i-1,j)
!cp1 &                      + tairoy(i-1,j-1) + tairoy(i,j-1))/(rho0*4)
!cp1        else
              phisu(i,j) = -(tairox(i,j))/(rho0)
              phisv(i,j) = -(tairoy(i,j))/(rho0)
!cp1        endif
          enddo
       enddo
!
        do j=js1,js2
           do i=is1(j),is2(j)
              tot=sqrt(tairox(i,j)**2+tairoy(i,j)**2)
              zustm=tot/rho0
              ust2s(i,j)=zustm*(one+sdvt(i,j))
           enddo
        enddo
!
      endif
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--4. Set umoy, vmoy at zero                                           |
!-----------------------------------------------------------------------

      do j=ju1,ju2
        do i=iu1(j),iu2(j)
!         umoy(i,j)=0.0
!         vmoy(i,j)=0.0
         enddo
      enddo

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine icdyna -
      end subroutine icdyna
