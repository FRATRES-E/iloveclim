!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:42 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:42 CET 2009

      SUBROUTINE adv(dt,ut,vt,s0,s1)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!
!    Avection methode Hibler.
!    Attention : pas conservatif au front de glace.
!    Pas toujours tres stable (instable si pad diffusion).
!

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use dynami_mod
!! END_OF_USE_SECTION


        real(kind=dblp), dimension(imax,jmax) :: ut, vt, s0, s1, div
     >                 , huc, hvc     
        real(kind=dblp) :: dt
        integer(kind=ip):: i, j, ip1, jp1, im1, jm1, jj


!
c~         dimension ut(imax,jmax),vt(imax,jmax),s0(imax,jmax)
c~         dimension s1(imax,jmax)
c~         dimension div(imax,jmax),huc(imax,jmax),hvc(imax,jmax)
!
        do j=js1-1,js2
           do i=is1(j)-1,is2(j)
              ip1=i+1
              huc(i,j)=(s1(i,j)/area(i,j)*dxc1(ip1,j)
     &                      +s1(ip1,j)/area(ip1,j)*dxc1(i,j))
     &                /(dxc1(ip1,j)+dxc1(i,j))
              huc(i,j)=.5*(s1(i,j)/area(i,j)+s1(ip1,j)/area(ip1,j))
           enddo
        enddo
!
        do j=js1-1,js2
           jp1=j+1
            do i=is1(j)-1,is2(j)
              hvc(i,j)=(s1(i,j)/area(i,j)*dxc2(i,jp1)
     &                       +s1(i,jp1)/area(i,jp1)*dxc2(i,j))
     &                 /(dxc2(i,j)+dxc2(i,jp1))
               hvc(i,j)=.5*(s1(i,j)/area(i,j)+s1(i,jp1)/area(i,jp1))
           enddo
        enddo
!
       do j=js1,js2
           jm1 = j-1
           do i=is1(j),is2(j)
              im1      = i-1
              div(i,j) = 2.0*
     &      (akappa(i,j,1,1)*(ut(i,j)*huc(i,j)-ut(im1,j)*huc(im1,j))-
     &        akappa(i,j,1,2)*(vt(i,jm1)*hvc(i,jm1)+vt(i,j)*hvc(i,j))-
     &        akappa(i,j,2,2)*(-vt(i,jm1)*hvc(i,jm1)+vt(i,j)*hvc(i,j))+
     &        akappa(i,j,2,1)*(ut(i,j)*huc(i,j)+ut(im1,j)*huc(im1,j)))
!             div(i,j)  = bkappa(i,j,1,1)*ut(i,j)*huc(i,j)
!     &                  -bkappa(i,j,1,2)*ut(im1,j)*huc(im1,j)
!     &                  +bkappa(i,j,2,2)*vt(i,j)*hvc(i,j)
!     &                  -bkappa(i,j,2,1)*vt(i,jm1)*hvc(i,jm1)
             s1(i,j)  = max(zero,s0(i,j)-area(i,j)*dt*div(i,j))
         enddo
       enddo
       do 200 jj=jcl1,jcl2
           s1(ims1-1,jj) = s1(ims2,jj)
           s1(ims2+1,jj) = s1(ims1,jj)
200    continue
!
       return
       end
