!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:45 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:45 CET 2009

      SUBROUTINE dynami_zh(ih)
!
! determine the velocity field of sea ice
!
! forcing: wind stress, water stress,
!           surface tilt
!
! ice-ice interaction: non-linear viscous-
!                       plastic law and bulk
!                       rheology
!
! method: Zhang and Hibler, 1997
!
! ih: hemispheric index (=1, NH; =2, SH)
! iu1: western most demarcation for velocity
!      points along a zonal belt (iu1.ge.2)
! iu2: eastern most demarcation for velocity
!      points along a zonal belt (iu2.le.imax-1)
!
! M.A. Morales Maqueda ?-?-?
! modified M.A. Morales Maqueda 31-12-1993
! modified H. Goosse 08-12-1994
! modified M.A. Morales Maqueda 09-02-2000
!
! Ccpl [Ccp0] => line pertaining to the
!                coupled (uncoupled) version
! Ccem        => concentric ellipse method
!                Hibler, 1979
! Ctem        => truncated ellipse method
!                Geiger et al., 1998
! Clds        => changes in ice concentration
!                associated with shearing
!


!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use ice_mod
      use dynami_mod
!! END_OF_USE_SECTION


      integer(kind=ip), parameter :: ntmax=1000

      real(kind=dblp), dimension(imax,jmax) :: presm, presh, gphs1
     >               ,gphs2, albqd, zmass, zcorl, usdtzm, uiner, viner
     >               ,uzcorl, vzcorl, tagnx, tagny, u0, v0, viszeta
     >               ,viseta, sigm11, sigm12, sigm22, c1m, c10, c1p, r1
     >               ,dvsg1, c2m, c20, c2p, r2, dvsg2, zaw, zbw, gamy
     >               ,up, vp, uaux=0.0_dblp, vaux=0.0_dblp, res, uvmsk
     &               ,zm_sc, za1ct
     >               ,za2ct, zb1, zb2, wstrn
      
      real(kind=dblp), dimension(ntmax) :: bety, gamx, betx, uu
      

!---  Other local variables

       integer(kind=ip):: ih, njmin, njmax, njmm1, imaxm1, j, i, jm1
     >                  , im1, iter, jp1, ip1, jter, is, ie



      real(kind=dblp) :: sang, ccgdx, ccgdy, usw, detadxi, detadyi
     >                  , zmabe, etapac, etaarc, tabe, pbe, pslbe, u1
     >                  , u2, ubering, zindu, usdtp, ze11, ze12, ze22
     >                  , ze21, trace, trace2, deter, zes2, delta, zes
     >                  , viseta1, aa, sigm1111, sigm2211, sigm1112
     >                  , sigm2212, sigm1121, sigm2221, sigm1122
     >                  , sigm2222, sigm1211, sigm1221, sigm1212
     >                  , sigm1222, ur, vr, zmod, zcd, bet, alp, gamma
     >                  , fact, resm, t1, t2, det, zden, zrug, zrvg

!
! define limits where to compute dynamics.
!
      if (ih.eq.1) then
        njmin = nlminn
        njmax = nlmaxn
      else
        njmin = nlmins
        njmax = nlmaxs
      endif
      njmm1  = njmin-1
      imaxm1 = imax-1
!
! initialise certain arrays. this could
! be done only once per run by introducing
! a flag and placing the arrays in a common
!
      do j=njmm1,njmax+1
        do i=1,imax
          uvmsk(i,j) = zero
          up(i,j)    = zero
          vp(i,j)    = zero
          zaw(i,j)   = one
          c10(i,j)   = zero
          r1(i,j)    = zero
          dvsg1(i,j) = zero
          c20(i,j)   = zero
          r2(i,j)    = zero
          dvsg2(i,j) = zero
        enddo
      enddo
!
! sign of turning angle for oceanic drag
!
      sang = real(ih)*sangvg
!
! ice mass, ice strength, and wind stress 
! at the center of the grid cells                                                
!
      do j=njmm1,njmax
        do i=1,imax
          zm_sc(i,j) = tms(i,j,ks2)
     &                 *(rhon*hnm(i,j)+rhog*hgm(i,j))
          wstrn(i,j) = exp(-c*albq(i,j))
          presm(i,j) = tms(i,j,ks2)
     &                 *pstarh*hgm(i,j)*wstrn(i,j)
!cp1     if (icoupl .eq. 0) then
          zb1(i,j)   = tms(i,j,ks2)*(one-albq(i,j))
          zb2(i,j)   = tms(i,j,ks2)*(one-albq(i,j))
!cp1     else
!cp1      zb1(i,j)   = tms(i,j,ks2)
!cp1                   *tenagx(i,j)*(one-albq(i,j))
!cp1      zb2(i,j)   = tms(i,j,ks2)
!cp1                   *tenagy(i,j)*(one-albq(i,j))
!cp1     endif
!
        enddo
      enddo

!
! lead area, mass, coriolis coefficients
! and wind stress at the corners of the 
! grid cells + momentum terms that are
! independent of the velocity field
!
      ccgdx = gpes*uns2dx
      ccgdy = gpes*uns2dy
      do 20 j=njmin,njmax
        jm1 = j-1
        do 10 i=iu1(j),iu2(j)
          im1        = i-1
!
          usw        = tmu(i,j,ku2)/
     &                 max( tms(i,j,ks2)*wght(i,j,2,2)
     &                     +tms(im1,j,ks2)*wght(i,j,1,2)
     &                     +tms(i,jm1,ks2)*wght(i,j,2,1)
     &                     +tms(im1,jm1,ks2)*wght(i,j,1,1),
     &                     zepsd2)
!
! leads area
!
         albqd(i,j)  = ( tms(i,j,ks2)*albq(i,j)*wght(i,j,2,2)
     &                  +tms(im1,j,ks2)*albq(im1,j)*wght(i,j,1,2)
     &                  +tms(i,jm1,ks2)*albq(i,jm1)*wght(i,j,2,1)
     &                  +tms(im1,jm1,ks2)*albq(im1,jm1)*wght(i,j,1,1))
     &                 *usw
!
! mass
!
          zmass(i,j) = ( zm_sc(i,j)*wght(i,j,2,2)
     &                  +zm_sc(im1,j)*wght(i,j,1,2)
     &                  +zm_sc(i,jm1)*wght(i,j,2,1)
     &                  +zm_sc(im1,jm1)*wght(i,j,1,1))
     &                 *usw
          uvmsk(i,j) = one-max(zero,sign(one,one-zmass(i,j)))
          zcorl(i,j) = zmass(i,j)*zfn(i,j)
!
! wind stress
!
!cp1     if (icoupl.eq.0) then
          tagnx(i,j) = (zb1(i,j)*wght(i,j,2,2)
     &                 +zb1(im1,j)*wght(i,j,1,2)
     &                 +zb1(i,jm1)*wght(i,j,2,1)
     &                 +zb1(im1,jm1)*wght(i,j,1,1))
     &                 *usw
     &                 *tenagx(i,j)
!
          tagny(i,j) = (zb2(i,j)*wght(i,j,2,2)
     &                 +zb2(im1,j)*wght(i,j,1,2)
     &                 +zb2(i,jm1)*wght(i,j,2,1)
     &                 +zb2(im1,jm1)*wght(i,j,1,1))
     &                 *usw
     &                 *tenagy(i,j)
!cp1     else
!cp1      tagnx(i,j) = (zb1(i,j)*wght(i,j,2,2)
!cp1 &                 +zb1(im1,j)*wght(i,j,1,2)
!cp1 &                 +zb1(i,jm1)*wght(i,j,2,1)
!cp1 &                 +zb1(im1,jm1)*wght(i,j,1,1))
!cp1 &                 *usw
!
!cp1      tagny(i,j) = (zb2(i,j)*wght(i,j,2,2)
!cp1 &                 +zb2(im1,j)*wght(i,j,1,2)
!cp1 &                 +zb2(i,jm1)*wght(i,j,2,1)
!cp1 &                 +zb2(im1,jm1)*wght(i,j,1,1))
!cp1 &                 *usw
!cp1     endif
!
! terms that are independent of the velocity field
! (surface tilt + wind stress)
!
         detadxi     = zmass(i,j)*ccgdx*smx(i,j,3)
     &                 *( (eta(im1,jm1)-eta(i,j))
     &                   +(eta(im1,j)-eta(i,jm1)))
         detadyi     = zmass(i,j)*ccgdy*smy(i,j,3)
     &                 *( (eta(im1,jm1)-eta(i,j))
     &                   -(eta(im1,j)-eta(i,jm1)))
         za1ct(i,j)  = tagnx(i,j)+detadxi
         za2ct(i,j)  = tagny(i,j)+detadyi
10      continue
20    continue
!
! cyclicity of velocity mask
!
      do j=njmin,njmax
        uvmsk(1,j)    = uvmsk(imaxm1,j)
        uvmsk(imax,j) = uvmsk(2,j)
      enddo
!
! initial conditions
!
      do j=njmin,njmax
        do i=1,imax
          ug(i,j) = uvmsk(i,j)*ug(i,j)
          vg(i,j) = uvmsk(i,j)*vg(i,j)
          u0(i,j) = ug(i,j)
          v0(i,j) = vg(i,j)
        enddo
      enddo
!
      if (ih.eq.1) then
!
! case of bering (compute velocity straight
!                 away!)
!
        zmabe   = ( rhon*hnm(iberp,jberp-1)
     &             +rhog*hgm(iberp,jberp-1)
     &             +rhon*hnm(iberp-1,jberp-1)
     &             +rhog*hgm(iberp-1,jberp-1)
     &             +rhon*hnm(ibera-1,jbera-1)
     &             +rhog*hgm(ibera-1,jbera-1)
     &             +rhon*hnm(ibera,jbera-1)
     &             +rhog*hgm(ibera,jbera-1))/4.0d0
        etapac  =  eta(iberp-1,jberp-3)+eta(iberp,jberp-3)
     &            +eta(iberp-1,jberp-2)+eta(iberp,jberp-2)
        etaarc  =  eta(ibera,jbera-3)+eta(ibera+1,jbera-3)
     &            +eta(ibera,jbera-2)+eta(ibera+1,jbera-2)
        detadyi =  zmabe*ccgdy*smy(iberp,jberp-1,3)
     &            *(etapac-etaarc)*0.25d0
        tabe    =  tenagy(iberp,jberp-1)+detadyi
        pbe     = ( presm(iberp,jberp-1)
     &             +presm(iberp-1,jberp-1)
     &             +presm(ibera,jbera-1)
     &             +presm(ibera-1,jbera-1))/4.0d0
        pslbe   =  pbe/(4.0d0*50.0d03)
        u1      =  max(zero,vo(iberp,jberp)
     &            +sign(one,tabe-pslbe)
     &             *sqrt(abs(tabe-pslbe)/rhoco))
        u2      =  min(zero,vo(iberp,jberp)
     &            +sign(one,tabe+pslbe)
     &             *sqrt(abs(tabe+pslbe)/rhoco))
        ubering =  (u1+u2)*50.0d0/125.0d0*sber
      endif
!
! solution of the momentum equation taking into
! account ice-ice interactions
!
      do 3000 iter=1,2*nbiter
!
        zindu = 0.5d0*mod(iter,2)
        usdtp = (one+2.0d0*zindu)*nbiter*usdt
!
! computation of free drift field for free slip
! boundary conditions
!
!       if (bound.gt.0.5d0) 
!    &    call frdrft(ih,usdtp,u0,v0,zmass,zcorl,tagnx,tagny)
!
! computation of viscosities
!
        do 40 j=njmm1,njmax
          jp1 = j+1
          do 30 i=1,imax
            ip1 = (i+1)-(imax-2)*(i/imax)
!
! rate of strain tensor
!
            ze11         = akappa(i,j,1,1)*
     &                      (ug(ip1,j)+ug(ip1,jp1)-(ug(i,j)+ug(i,jp1)))
     &                    +akappa(i,j,1,2)*
     &                      (vg(ip1,j)+vg(ip1,jp1)+vg(i,j)+vg(i,jp1))
            ze12         = akappa(i,j,2,2)*
     &                      (ug(i,jp1)+ug(ip1,jp1)-(ug(i,j)+ug(ip1,j)))
     &                    -akappa(i,j,2,1)*
     &                      (vg(i,j)+vg(ip1,j)+vg(i,jp1)+vg(ip1,jp1))
            ze22         = akappa(i,j,2,2)*
     &                      (vg(i,jp1)+vg(ip1,jp1)-(vg(i,j)+vg(ip1,j)))
     &                    +akappa(i,j,2,1)*
     &                      (ug(i,j)+ug(ip1,j)+ug(i,jp1)+ug(ip1,jp1))
            ze21         = akappa(i,j,1,1)*
     &                      (vg(ip1,j)+vg(ip1,jp1)-(vg(i,j)+vg(i,jp1)))
     &                    -akappa(i,j,1,2)*
     &                      (ug(ip1,j)+ug(ip1,jp1)+ug(i,j)+ug(i,jp1))
            trace        = ze11+ze22
            trace2       = trace*trace
            deter        = ze11*ze22-0.25d0*(ze12+ze21)**2
            zes2         = max(zero,trace2-4.0d0*deter)
            delta        = max(sqrt(trace2+zes2*usecc2),creepl)
            viszeta(i,j) = max(presm(i,j)/delta,zetamn)
!cem        viseta(i,j)  = viszeta(i,j)*usecc2
!cem        presh(i,j)   = presm(i,j)
            presh(i,j)   = delta*viszeta(i,j)
            zes          = max(zepsd2,sqrt(zes2))
            viseta1      = (presh(i,j)-viszeta(i,j)*trace)/zes
            viseta(i,j)  = min(viszeta(i,j)*usecc2,viseta1)
!
! Determination of stress tensor
!
            aa           = viszeta(i,j)*(ze11+ze22)
            sigm11(i,j)  = viseta(i,j)*(ze11-ze22)+aa
            sigm12(i,j)  = viseta(i,j)*(ze12+ze21)
            sigm22(i,j)  = viseta(i,j)*(ze22-ze11)+aa
30        continue
40      continue
!
! gradient of ice strength and some auxiliary arrays
!
        do j=njmin,njmax
          jm1 = j-1
          do i=iu1(j),iu2(j)
            im1        =  i-1
            gphs1(i,j) =  (alambd(i,j,2,2,2,1)-alambd(i,j,2,1,2,1))
     &                    *presh(i,jm1)
     &                   +(alambd(i,j,2,2,2,2)-alambd(i,j,2,1,2,2))
     &                    *presh(i,j)
     &                   -(alambd(i,j,2,2,1,1)+alambd(i,j,2,1,1,1))
     &                    *presh(im1,jm1)
     &                   -(alambd(i,j,2,2,1,2)+alambd(i,j,2,1,1,2))
     &                    *presh(im1,j)
            gphs2(i,j) = -(alambd(i,j,1,1,2,1)+alambd(i,j,1,2,2,1))
     &                    *presh(i,jm1)
     &                   -(alambd(i,j,1,1,1,1)+alambd(i,j,1,2,1,1))
     &                    *presh(im1,jm1)
     &                   +(alambd(i,j,1,1,2,2)-alambd(i,j,1,2,2,2))
     &                    *presh(i,j)
     &                   +(alambd(i,j,1,1,1,2)-alambd(i,j,1,2,1,2))
     &                    *presh(im1,j)
!
          usdtzm(i,j)  = usdtp*zmass(i,j)
          uiner(i,j)   = usdtzm(i,j)*u0(i,j)
          viner(i,j)   = usdtzm(i,j)*v0(i,j)
          uzcorl(i,j)  = ug(i,j)*zcorl(i,j)
          vzcorl(i,j)  = vg(i,j)*zcorl(i,j)
          enddo           
        enddo           

!
        do j=njmin,njmax
          jm1 = j-1
          do i=1,imax
            im1      = (i-1)+(imax-2)*(1/i)
!
! determination of c1m
!
            ze11     = -akappa(im1,jm1,1,1)
            sigm1111 = (viszeta(im1,jm1)+viseta(im1,jm1))*ze11
            sigm2211 = (viszeta(im1,jm1)-viseta(im1,jm1))*ze11
            ze11     = -akappa(im1,j,1,1)
            sigm1112 = (viszeta(im1,j)+viseta(im1,j))*ze11
            sigm2212 = (viszeta(im1,j)-viseta(im1,j))*ze11
            c1m(i,j) =  alambd(i,j,2,2,1,1)*sigm1111
     &                 +alambd(i,j,2,2,1,2)*sigm1112
     &                 +alambd(i,j,2,1,1,1)*sigm2211
     &                 +alambd(i,j,2,1,1,2)*sigm2212
!
! determination of c2m
!
            ze22     = -akappa(im1,jm1,2,2)
            sigm1111 = (viszeta(im1,jm1)-viseta(im1,jm1))*ze22
            sigm2211 = (viszeta(im1,jm1)+viseta(im1,jm1))*ze22
            ze22     = -akappa(i,jm1,2,2)
            sigm1121 = (viszeta(i,jm1)-viseta(i,jm1))*ze22
            sigm2221 = (viszeta(i,jm1)+viseta(i,jm1))*ze22
            c2m(i,j) =  alambd(i,j,1,1,2,1)*sigm2221
     &                 +alambd(i,j,1,1,1,1)*sigm2211
     &                 +alambd(i,j,1,2,1,1)*sigm1111
     &                 +alambd(i,j,1,2,2,1)*sigm1121
          enddo
        enddo
!
! determination of c1p and c2p
!
        do j=njmin,njmax
          jp1 = j+1
          do i=1,imax
            ip1      = (i+1)-(imax-2)*(i/imax)
            c1p(i,j) = c1m(ip1,j)
            c2p(i,j) = c2m(i,jp1)
          enddo
        enddo
!
! special boundary case, c2p(i,njmax)
!
        j = njmax
        do i=1,imax
          im1      = (i-1)+(imax-2)*(1/i)
          ze22     = akappa(im1,j,2,2)
          sigm1112 = (viszeta(im1,j)-viseta(im1,j))*ze22
          sigm2212 = (viszeta(im1,j)+viseta(im1,j))*ze22
          ze22     = akappa(i,j,2,2)
          sigm1122 = (viszeta(i,j)-viseta(i,j))*ze22
          sigm2222 = (viszeta(i,j)+viseta(i,j))*ze22
          c2p(i,j) = -alambd(i,j,1,1,2,2)*sigm2222
     &               -alambd(i,j,1,1,1,2)*sigm2212
     &               +alambd(i,j,1,2,1,2)*sigm1112
     &               +alambd(i,j,1,2,2,2)*sigm1122 
        enddo
!
        do 60 j=njmin,njmax
          jm1 = j-1
          do 50 i=iu1(j),iu2(j)
            im1 = i-1
!
! determination of c10
!
            ze12     =  akappa(im1,jm1,2,2)
            sigm1211 =  viseta(im1,jm1)*ze12
            ze12     =  akappa(i,jm1,2,2)
            sigm1221 =  viseta(i,jm1)*ze12
            ze12     = -akappa(im1,j,2,2)
            sigm1212 =  viseta(im1,j)*ze12
            ze12     = -akappa(i,j,2,2)
            sigm1222 =  viseta(i,j)*ze12
!
            c10(i,j) =  (alambd(i,j,1,1,1,1)-alambd(i,j,1,2,1,1))
     &                  *sigm1211
     &                 +(alambd(i,j,1,1,2,1)-alambd(i,j,1,2,2,1))
     &                  *sigm1221
     &                 -(alambd(i,j,1,1,2,2)+alambd(i,j,1,2,2,2))
     &                  *sigm1222
     &                 -(alambd(i,j,1,1,1,2)+alambd(i,j,1,2,1,2))
     &                  *sigm1212
     &                 -(c1m(i,j)+c1p(i,j))
!
! determination of c20
! 
            ze21     =  akappa(im1,jm1,1,1)
            sigm1211 =  viseta(im1,jm1)*ze21
            ze21     = -akappa(i,jm1,1,1)
            sigm1221 =  viseta(i,jm1)*ze21
            ze21     =  akappa(im1,j,1,1)
            sigm1212 =  viseta(im1,j)*ze21
            ze21     = -akappa(i,j,1,1)
            sigm1222 =  viseta(i,j)*ze21
!
            c20(i,j) = -(alambd(i,j,2,1,2,1)+alambd(i,j,2,2,2,1))
     &                  *sigm1221
     &                 -(alambd(i,j,2,1,2,2)+alambd(i,j,2,2,2,2))
     &                  *sigm1222
     &                 +(alambd(i,j,2,2,1,1)-alambd(i,j,2,1,1,1))
     &                  *sigm1211
     &                 +(alambd(i,j,2,2,1,2)-alambd(i,j,2,1,1,2))
     &                  *sigm1212
     &                 -(c2m(i,j)+c2p(i,j))
!
50        continue
60      continue
!
! mask c1m, c1p, c2m, and c2p
!
        do j=njmin,njmax
          do i=1,imax
            c1m(i,j) = c1m(i,j)*uvmsk(i,j)
            c1p(i,j) = c1p(i,j)*uvmsk(i,j)
            c2m(i,j) = c2m(i,j)*uvmsk(i,j)
            c2p(i,j) = c2p(i,j)*uvmsk(i,j)
          enddo
        enddo
!
! relaxation
!
        do j=njmm1,njmax
          do i=1,imax
            uaux(i,j) = ug(i,j)
            vaux(i,j) = vg(i,j)
          enddo
        enddo
!
        do 1000 jter=1,nbitdr 
!
          do j=njmin,njmax
            jm1 = j-1
            jp1 = j+1
            do i=iu1(j),iu2(j)
              im1        = i-1
              ip1        = i+1 
              up(i,j)    = uaux(i,j)
              vp(i,j)    = vaux(i,j)
!
! determination of right-hand side term of momentum equation
! (it is not possible to linearise the ice-ocean drag term,
! as doing so leads to important inaccuracies in the momentum
! balance when running the routine without pseudo-timesteps...)
!
              ur         = uaux(i,j)-uo(i,j)
              vr         = vaux(i,j)-vo(i,j)
              zmod       = sqrt(ur*ur+vr*vr)*(one-albqd(i,j))
              zcd        = rhoco*zmod
              zaw(i,j)   = one+uvmsk(i,j)*(usdtzm(i,j)+zcd*cangvg-one)
              r1(i,j)    = ( uiner(i,j)
     &                      +vzcorl(i,j)
     &                      +za1ct(i,j)
     &                      +zcd*( cangvg*uo(i,j)
     &                            -sang*(vo(i,j)-vaux(i,j)))
     &                      -gphs1(i,j))
     &                     *uvmsk(i,j)
              r2(i,j)    = ( viner(i,j)
     &                      -uzcorl(i,j)
     &                      +za2ct(i,j)
     &                      +zcd*( cangvg*vo(i,j)
     &                            +sang*(uo(i,j)-uaux(i,j)))
     &                      -gphs2(i,j))
     &                     *uvmsk(i,j)
!
! determination of divergence of stress tensor
!
              dvsg1(i,j) =  alambd(i,j,2,2,2,1)*sigm11(i,jm1)
     &                     +alambd(i,j,2,2,2,2)*sigm11(i,j)
     &                     -alambd(i,j,2,2,1,1)*sigm11(im1,jm1)
     &                     -alambd(i,j,2,2,1,2)*sigm11(im1,j)
     &                     -alambd(i,j,2,1,1,1)*sigm22(im1,jm1)
     &                     -alambd(i,j,2,1,2,1)*sigm22(i,jm1)
     &                     -alambd(i,j,2,1,1,2)*sigm22(im1,j)
     &                     -alambd(i,j,2,1,2,2)*sigm22(i,j)
     &                     +( alambd(i,j,1,2,2,1)
     &                       -alambd(i,j,1,1,2,1))*sigm12(i,jm1)
     &                     +( alambd(i,j,1,2,1,1)
     &                       -alambd(i,j,1,1,1,1))*sigm12(im1,jm1)
     &                     +( alambd(i,j,1,1,2,2)
     &                       +alambd(i,j,1,2,2,2))*sigm12(i,j)
     &                     +( alambd(i,j,1,1,1,2)
     &                       +alambd(i,j,1,2,1,2))*sigm12(im1,j)
              dvsg1(i,j) = ( dvsg1(i,j)
     &                      +c1m(i,j)*uaux(im1,j)
     &                      +c10(i,j)*uaux(i,j)
     &                      +c1p(i,j)*uaux(ip1,j))
     &                     *uvmsk(i,j)
!
              dvsg2(i,j) = -alambd(i,j,1,2,1,1)*sigm11(im1,jm1)
     &                     -alambd(i,j,1,2,2,1)*sigm11(i,jm1)
     &                     -alambd(i,j,1,2,1,2)*sigm11(im1,j)
     &                     -alambd(i,j,1,2,2,2)*sigm11(i,j)
     &                     -alambd(i,j,1,1,2,1)*sigm22(i,jm1)
     &                     -alambd(i,j,1,1,1,1)*sigm22(im1,jm1)
     &                     +alambd(i,j,1,1,2,2)*sigm22(i,j)
     &                     +alambd(i,j,1,1,1,2)*sigm22(im1,j)
     &                     +( alambd(i,j,2,1,1,1)
     &                       -alambd(i,j,2,2,1,1))*sigm12(im1,jm1)
     &                     +( alambd(i,j,2,1,2,1)
     &                       +alambd(i,j,2,2,2,1))*sigm12(i,jm1)
     &                     +( alambd(i,j,2,1,1,2)
     &                       -alambd(i,j,2,2,1,2))*sigm12(im1,j)
     &                     +( alambd(i,j,2,1,2,2)
     &                       +alambd(i,j,2,2,2,2))*sigm12(i,j)
              dvsg2(i,j) = ( dvsg2(i,j)
     &                      +c2m(i,j)*vaux(i,jm1)
     &                      +c20(i,j)*vaux(i,j)
     &                      +c2p(i,j)*vaux(i,jp1))
     &                     *uvmsk(i,j)
            enddo
          enddo
!
          if (ih.eq.1) then
!
! case of bering (ice drift at bering
! strait is used as boundary condition)
!
            i          = ibera
            j          = jbera-1
            dvsg2(i,j) = ( dvsg2(i,j)
     &                    -c2p(i,j)*vaux(i,j+1))
     &                   *uvmsk(i,j)
            i          = iberp
            j          = jberp-1
            dvsg2(i,j) = ( dvsg2(i,j)
     &                    -c2p(i,j)*vaux(i,j+1))
     &                   *uvmsk(i,j)
          endif
!
! solve tridiagonal system in the x-direction.
!
          do 70 j=njmin,njmax
            is         = iu1(j)
            ie         = iu2(j)
!
! cyclic/non-cyclic conditions will be applied 
! depending on whether, at a given latitude, the 
! domain spans a full zonal belt or not. 
!
            if (ie-is.lt.imax-3) then
!
! non-cyclic tridiagonal case
!
              bet        = c10(is,j)+zaw(is,j)
              uaux(is,j) = (r1(is,j)+dvsg1(is,j))/bet
              do i=is+1,ie
                gamx(i)   = c1p(i-1,j)/bet
                bet       = c10(i,j)+zaw(i,j)-c1m(i,j)*gamx(i)
                uaux(i,j) = ( r1(i,j)+dvsg1(i,j)
     &                       -c1m(i,j)*uaux(i-1,j))/bet
              enddo
              do i=ie-1,is,-1
                uaux(i,j) = uaux(i,j)-gamx(i+1)*uaux(i+1,j)
              enddo
            else
!
! cyclic tridiagonal case
!
! first tridiagonal step 
!
              alp        = c1p(ie,j)
              bet        = c1m(is,j)
              gamma      = -(c10(is,j)+zaw(is,j)) 
              betx(is)   = -(gamma+gamma)
              uaux(is,j) = (r1(is,j)+dvsg1(is,j))/betx(is)
              do i=is+1,ie-1
                gamx(i)   = c1p(i-1,j)/betx(i-1)
                betx(i)   = c10(i,j)+zaw(i,j)-c1m(i,j)*gamx(i)
                uaux(i,j) = ( r1(i,j)+dvsg1(i,j)
     &                       -c1m(i,j)*uaux(i-1,j))/betx(i)
              enddo
              gamx(ie)   = c1p(ie-1,j)/betx(ie-1)
              betx(ie)   = c10(ie,j)+zaw(ie,j)-alp*bet/gamma
     &                     -c1m(ie,j)*gamx(ie)
              uaux(ie,j) = ( r1(ie,j)+dvsg1(ie,j)
     &                      -c1m(ie,j)*uaux(ie-1,j))/betx(ie)
              do i=ie-1,is,-1
                uaux(i,j) = uaux(i,j)-gamx(i+1)*uaux(i+1,j)
              enddo
!
! second tridiagonal step
!
              uu(is) = -0.5d0
              do i=is+1,ie-1
                uu(i) = -c1m(i,j)*uu(i-1)/betx(i)
              enddo
              uu(ie) = (alp-c1m(ie,j)*uu(ie-1))/betx(ie)
              do i=ie-1,is,-1
                uu(i) = uu(i)-gamx(i+1)*uu(i+1)
              enddo
              fact = (uaux(is,j)+bet*uaux(ie,j)/gamma)
     &               /(one+uu(is)+bet*uu(ie)/gamma)
              do i=is,ie
                uaux(i,j) = uaux(i,j)-fact*uu(i)
              enddo
            endif
70        continue
!
! solve tridiagonal system in the y-direction
!
          do i=2,imaxm1
            bety(i)       = c20(i,njmin)+zaw(i,njmin)
            vaux(i,njmin) = (r2(i,njmin)+dvsg2(i,njmin))/bety(i)
          enddo
          do j=njmin+1,njmax
            do i=2,imaxm1
              gamy(i,j) = c2p(i,j-1)/bety(i)
              bety(i)   = c20(i,j)+zaw(i,j)-c2m(i,j)*gamy(i,j)
              vaux(i,j) = (r2(i,j)+dvsg2(i,j)-c2m(i,j)*vaux(i,j-1))
     &                    /bety(i)
            enddo
          enddo
          do j=njmax-1,njmin,-1
            do i=2,imaxm1
              vaux(i,j) = vaux(i,j)-gamy(i,j+1)*vaux(i,j+1)
            enddo
          enddo
!
! apply over/under-relaxation
!
          do j=njmin,njmax
            do i=2,imaxm1
              uaux(i,j) = up(i,j)+om*(uaux(i,j)-up(i,j))
              vaux(i,j) = vp(i,j)+om*(vaux(i,j)-vp(i,j))
            enddo
          enddo
!
          if (ih.eq.1) then
!
! case of bering
!
            uaux(ibera,jbera) = zero
            vaux(ibera,jbera) = -ubering
            uaux(iberp,jberp) = zero
            vaux(iberp,jberp) = ubering
          endif
!
! cyclicity
!
          do j=njmin,njmax
            uaux(1,j)    = uaux(imaxm1,j)
            uaux(imax,j) = uaux(2,j)
            vaux(1,j)    = vaux(imaxm1,j)
            vaux(imax,j) = vaux(2,j)
          enddo
!
! convergence test
!
          resm = zero
          do j=njmin,njmax
            do i=iu1(j),iu2(j)
              resm = max(resm,abs(uaux(i,j)-up(i,j)),
     &                        abs(vaux(i,j)-vp(i,j)))
            enddo
          enddo
          if(resm.lt.resl) go to 2000
!
! new stress tensor
!
          do j=njmm1,njmax
            jp1 = j+1
            do i=1,imax
              ip1  = (i+1)-(imax-2)*(i/imax)
!
! rate of strain tensor
!
              ze11 = akappa(i,j,1,1)*
     &               (uaux(ip1,j)+uaux(ip1,jp1)-(uaux(i,j)+uaux(i,jp1)))
     &              +akappa(i,j,1,2)*
     &               (vaux(ip1,j)+vaux(ip1,jp1)+vaux(i,j)+vaux(i,jp1))
              ze12 = akappa(i,j,2,2)*
     &               (uaux(i,jp1)+uaux(ip1,jp1)-(uaux(i,j)+uaux(ip1,j)))
     &              -akappa(i,j,2,1)*
     &               (vaux(i,j)+vaux(ip1,j)+vaux(i,jp1)+vaux(ip1,jp1))
              ze22 = akappa(i,j,2,2)*
     &               (vaux(i,jp1)+vaux(ip1,jp1)-(vaux(i,j)+vaux(ip1,j)))
     &              +akappa(i,j,2,1)*
     &               (uaux(i,j)+uaux(ip1,j)+uaux(i,jp1)+uaux(ip1,jp1))
              ze21 = akappa(i,j,1,1)*
     &               (vaux(ip1,j)+vaux(ip1,jp1)-(vaux(i,j)+vaux(i,jp1)))
     &              -akappa(i,j,1,2)*
     &               (uaux(ip1,j)+uaux(ip1,jp1)+uaux(i,j)+uaux(i,jp1))
!
              aa          = viszeta(i,j)*(ze11+ze22)
!
! is the relaxation of sigm really neccesary?
! the algorithm exhibits sometimes instabilities
! when the stress tensor is not relaxed. is
! this due to the fact that the scheme is a 9-point
! one? 
!
              sigm11(i,j) = sigm11(i,j)+om*( viseta(i,j)*(ze11-ze22)+aa
     &                                      -sigm11(i,j))
              sigm12(i,j) = sigm12(i,j)+om*( viseta(i,j)*(ze12+ze21)
     &                                      -sigm12(i,j))
              sigm22(i,j) = sigm22(i,j)+om*( viseta(i,j)*(ze22-ze11)+aa
     &                                      -sigm22(i,j))
            enddo
          enddo
!
1000    continue
2000    continue
!
! include coriolis term correction
!
        do j=njmin,njmax
          do i=iu1(j),iu2(j)
            ur        = uaux(i,j)-uo(i,j)
            vr        = vaux(i,j)-vo(i,j)
            zmod      = sqrt(ur*ur+vr*vr)*(one-albqd(i,j))
            zcd       = rhoco*zmod
            zbw(i,j)  = alpha*zcorl(i,j)+zcd*sang 
            t1        = zaw(i,j)*uaux(i,j)-zbw(i,j)*vg(i,j)
            t2        = zaw(i,j)*vaux(i,j)+zbw(i,j)*ug(i,j)
            det       = zaw(i,j)*zaw(i,j)+zbw(i,j)*zbw(i,j)
            zden      = sign(one,det)/max(zepsd2,abs(det))
     &                  *uvmsk(i,j)
            zrug      = (zaw(i,j)*t1+zbw(i,j)*t2)*zden
            zrvg      = (zaw(i,j)*t2-zbw(i,j)*t1)*zden
            ug(i,j)   = zindu*(ug(i,j)-zrug)+zrug
            vg(i,j)   = zindu*(vg(i,j)-zrvg)+zrvg
          enddo
        enddo
!
        if (ih.eq.1) then
!
! case of bering
!
          ug(ibera,jbera) = zero
          vg(ibera,jbera) = -ubering
          ug(iberp,jberp) = zero
          vg(iberp,jberp) = ubering
        endif
!
! cyclicity
!
        do j=njmin,njmax
          ug(1,j)    = ug(imaxm1,j)
          ug(imax,j) = ug(2,j)
          vg(1,j)    = vg(imaxm1,j)
          vg(imax,j) = vg(2,j)
        enddo
!
! initial conditions
!
        do j=njmin,njmax
          do i=iu1(j),iu2(j)
            u0(i,j) = 2.0d0*zindu*(u0(i,j)-ug(i,j))+ug(i,j)
            v0(i,j) = 2.0d0*zindu*(v0(i,j)-vg(i,j))+vg(i,j)
          enddo
        enddo
!
!       if (iter.eq.2.or.iter.eq.2*nbiter) write(79,*) ih,iter,resm,jter
!
3000  continue
!
!     normalised rate of energy disipation by 
!     shear and convergence ---> extra open-water production 
!
      do j=njmm1,njmax
        jp1 = j+1
        do i=1,imax
          ip1       = (i+1)-(imax-2)*(i/imax)
          ze11      = akappa(i,j,1,1)*
     &                 (ug(ip1,j)+ug(ip1,jp1)-(ug(i,j)+ug(i,jp1)))
     &               +akappa(i,j,1,2)*
     &                 (vg(ip1,j)+vg(ip1,jp1)+vg(i,j)+vg(i,jp1))
          ze12      = akappa(i,j,2,2)*
     &                 (ug(i,jp1)+ug(ip1,jp1)-(ug(i,j)+ug(ip1,j)))
     &               -akappa(i,j,2,1)*
     &                 (vg(i,j)+vg(ip1,j)+vg(i,jp1)+vg(ip1,jp1))
          ze22      = akappa(i,j,2,2)*
     &                 (vg(i,jp1)+vg(ip1,jp1)-(vg(i,j)+vg(ip1,j)))
     &               +akappa(i,j,2,1)*
     &                 (ug(i,j)+ug(ip1,j)+ug(i,jp1)+ug(ip1,jp1))
          ze21      = akappa(i,j,1,1)*
     &                 (vg(ip1,j)+vg(ip1,jp1)-(vg(i,j)+vg(i,jp1)))
     &               -akappa(i,j,1,2)*
     &                (ug(ip1,j)+ug(ip1,jp1)+ug(i,j)+ug(i,jp1))
          trace     = ze11+ze22
          trace2    = trace*trace
          deter     = ze11*ze22-0.25d0*(ze12+ze21)**2
          delta     = sqrt(trace2+(trace2-4.0d0*deter)*usecc2)
          alcr(i,j) = (1.0d0-max(0.0d0,sign(1.0d0,-(1.0d0-albq(i,j)))))*
     &                0.5*(delta-trace)*wstrn(i,j)
        enddo
      enddo
!
      return
      end
