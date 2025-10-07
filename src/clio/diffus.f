!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009

      SUBROUTINE diffus(dt,difhx,difhy,fld0,fld1)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----
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


        real(kind=dblp), dimension(imax,jmax) :: difhx, difhy, fld0
     &                 , fld1, fld2, grxh, gryh, div, fact
     
        !dmr [FIX] -> variable div0 was used uninitialized all the time
        !          -> tentative fix with "save" and init to zero /
        !          Hugues G. says no "save"
        real(kind=dblp), dimension(imax,jmax) :: div0 = 0.0_dblp

        real(kind=dblp) :: dt, alfa, zindk, fldnw, amx
        integer(kind=ip):: j, i, k, jp1, ip1, jm1, im1, jj

c~         dimension difhx(imax,jmax),difhy(imax,jmax),
c~      &                  fld0(imax,jmax),fld1(imax,jmax)
c~ !
c~         dimension fld2(imax,jmax)
c~         dimension grxh(imax,jmax),gryh(imax,jmax),
c~      &            div0(imax,jmax),div(imax,jmax),fact(imax,jmax)
!
!  Fully implicit (alfa = 1) or CRANK-NICHOLSON (alfa = 0.5),
!  or fully  explicit (alfa = 0.0).
!
        alfa = 0.5
!
!  Relaxation constant.
!
        om = 0.5
!
        do j=1,jmax
          do i=1,imax
            fld1(i,j) = fld0(i,j)
            grxh(i,j) = 0.0
            gryh(i,j) = 0.0
          enddo
        enddo
!
      kloop:   do k=1,100
          zindk = real(1/k)
          do j=js1-1,js2
            jp1 = (j+1)
            do i=is1(j)-1,is2(j)
              ip1       = (i+1)
              grxh(i,j) = difhx(i,j)*(fld1(ip1,j)-fld1(i,j))
              gryh(i,j) = difhy(i,j)*(fld1(i,jp1)-fld1(i,j))
            enddo
          enddo

          do j=js1,js2
            jm1=j-1
            do i=is1(j),is2(j)
              im1=i-1
              fact(i,j) = bkappa(i,j,1,1)+bkappa(i,j,1,2)
     &                         +bkappa(i,j,2,2)+bkappa(i,j,2,1)        
!              fact(i,j) = 2.0*
!     &                    (difhx(i,j)*(akappa(i,j,1,1)+akappa(i,j,2,1))+
!     &                     difhx(im1,j)*(akappa(i,j,1,1)-akappa(i,j,2,1))+
!     &                     difhy(i,jm1)*(-akappa(i,j,1,2)+akappa(i,j,2,2))+
!     &                     difhy(i,j)*(akappa(i,j,2,2)+akappa(i,j,1,2)))
            enddo
          enddo
!
          do j=js1,js2
            jm1 = j-1
            do i=is1(j),is2(j)
              im1       = i-1
              div(i,j)  = bkappa(i,j,1,1)*grxh(i,j)
     &                          -bkappa(i,j,1,2)*grxh(im1,j)
     &                   +bkappa(i,j,2,2)*gryh(i,j)
     &                   -bkappa(i,j,2,1)*gryh(i,jm1)
!              div(i,j)  = 2.0*
!     &                    (akappa(i,j,1,1)*(grxh(i,j)-grxh(im1,j))+
!     &                     akappa(i,j,1,2)*(gryh(i,jm1)+gryh(i,j))+
!     &                     akappa(i,j,2,2)*(gryh(i,j)-gryh(i,jm1))+
!     &                     akappa(i,j,2,1)*(grxh(i,j)+grxh(im1,j)))
              div0(i,j) = zindk*div(i,j)+(1.0-zindk)*div0(i,j)
            enddo
          enddo

!
          do j=js1,js2
            do i=is1(j),is2(j)
              fldnw     = (fld0(i,j)+dt*
     &                           (alfa*(div(i,j)+fact(i,j)*fld1(i,j))+
     &                            (1.0-alfa)*div0(i,j)))/
     &                                (1.0+alfa*dt*fact(i,j))
              fld2(i,j) = fld1(i,j)+om*(fldnw-fld1(i,j))
            enddo
          enddo
!
          amx =  0.0
          do j=js1,js2
            do i=is1(j),is2(j)
              amx = max(amx,abs(fld2(i,j)-fld1(i,j)))
            enddo
          enddo

          do j=js1,js2
            do i=is1(j),is2(j)
              fld1(i,j) = fld2(i,j)
            enddo
          enddo

          do jj=jcl1,jcl2
             fld1(ims1-1,jj) = fld1(ims2,jj)
             fld1(ims2+1,jj) = fld1(ims1,jj)
          enddo
          
          if (amx.lt.2.0e-04) exit kloop
        enddo kloop
!
!       write(53,*) amx,k
        return
        end
