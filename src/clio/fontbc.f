!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:46 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:46 CET 2009


      SUBROUTINE fontbc(kideb,kiut)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!
! This routine determines the time evolution of sea ice
! thickness, concentration and heat content due to the
! vertical and lateral thermodynamic accretion-ablation
! processes.
!
!---
! Ccpl [Ccp0] => ligne specifique a la version avec [sans] couplage .
!---

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip      

      use const_mod, only: zero, one, rho0

      use para0_mod, only: nbpt
      use para_mod, only: 
      use bloc0_mod, only: 
      use bloc_mod, only: icoupl
      use thermo_mod, only: tbqb, fstbqb, qstbqb, fratsb, fcsb, fleb
     &              , tsb, hnbqb, dvsbqb, dmnbqb, hgbqb, dvbbqb, dmgbqb
     &              , qlbqb, dvnbqb, dmgwib, albqb, fscbqb, qfvbqb
     &              , fltbqb, dvlbqb, cldqb, fsolgb, psbqb, tabqb 
     &              , vabqb, ratbqb, tfub, qabqb, hnpbqb, fbbqb
     &              , qlbbqb, thcmb

      use ice_mod, only: hndif, hgdif, tfsn, rcpn, tfsg, rcpg, emig
     >                 , hnzst, xkn, xkg, hth, xlg, hgmin, ddtb, swiqst
     >                 , stefan, beta, xln, rhon, parsub, xsn, rhog
     >                 , hmelt, alphs, cnscg, hakspl, nbits

      implicit none

      real(kind=dblp), dimension(nbpt)  ::  zep, zqsat, zemin, zrchu,
     &          ztfs, zksdh, ztbq, zfs, zab, zqmax, zindn, zindg, zfsup,
     &          zfocea, zffs, zdhg, zdhb, zalbqb, ykn, ykg, rcpdt, ztst
      real(kind=dblp), dimension(nbpt,2)::  qctbqb
      
!--- more locales
      integer(kind=ip), intent(in):: kideb, kiut
      integer(kind=ip):: ji, jt
      real(kind=dblp) :: zeps, zeps0, zeps1, zind1, zind23, zignn, zignm
     >                  , zig, zigm, xknw, xkgw, he, zhe, heshth, gdif
     >                  , zk, zdh, zhg, zst, zexpo, zfsab, zhq, rhoa, ce
     >                  , surthi, es, zssdqw, dzf, zf, dtsb, dtemax, umb
     >                  , zind3, zks, zkg, zd, zaux1, zaux2, zfsb, zkgb
     >                  , za1, za2, za3, zc1, zc2, zc3, zb1, zb2, zb3
     >                  , zd1, zd2, zd3, ztci, zti, zdhsm, zhf, zhn
     >                  , zqfont, zqf, ziqf, zdhcf, zhnf, zqnes, zqr
     >                  , ziqr, ztb, zf2, zihq, zidhb, zdhbf, zhgnew
     >                  , zhni, zihgnew, zihg, zdhgm, zdhnm, zhnfi
     >                  , zpc1, zpc2, zp1, zp2, dzrmh, zhnnew, quot
     >                  , tneq, zia, zian, zq, zqtot, ziq, dvn

      
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  1) Tolerance parameters.                                            |
!-----------------------------------------------------------------------
!
      zeps  = 1.0e-20
      zeps0 = 1.0e-13
      zeps1 = 1.0e-06
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  2) If tbqb(ji,1) > tfsn, tbqb(ji,1) = tfsn.                 |
!     If tbqb(ji,2/3) > tfsg, tbqb(ji,2/3) = tfsg.
!-----------------------------------------------------------------------
!
       do ji=kideb,kiut

        zind1        = max(zero,sign(one,hndif-hnbqb(ji)))
        zind23       = max(zero,sign(one,hgdif-hgbqb(ji)))
        qctbqb(ji,1) = max(zero,(tbqb(ji,1)-tfsn)*rcpn*hnbqb(ji))
     &                 *(1.0-zind1)
        qctbqb(ji,2) = (max(zero,(tbqb(ji,2)-tfsg)*rcpg
     &                *(hgbqb(ji)/2.0))+
     &                max(zero,(tbqb(ji,3)-tfsg)*rcpg
     &                *(hgbqb(ji)/2.0)))*
     &                (1.0-zind23)
        tbqb(ji,1)   = min(tfsn,tbqb(ji,1))
        tbqb(ji,2)   = min(tfsg,tbqb(ji,2))
        tbqb(ji,3)   = min(tfsg,tbqb(ji,3))
        
        enddo
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  3) Calculate some intermediate variables.                           |
!-----------------------------------------------------------------------
!
      do ji=kideb,kiut
        zemin(ji)  = emig
        zignn      = max(zero,sign(one,hnzst-hnbqb(ji)))
        zignm      = max(zero,sign(one,hndif-hnbqb(ji)))
        zignm      = max(zignm,zignn)
        zig        = max(zero,sign(one,-hnbqb(ji)))
        zigm       = max(zero,sign(one,hgdif-hgbqb(ji)))
        ztfs(ji)   = (1.0-zig)*tfsn+zig*tfsg
!
! Effective conductivity
!
        xknw       = xkn/(xkn+xkg)
        xkgw       = xkg/(xkn+xkg)
        he         = xknw*hgbqb(ji)+xkgw*hnbqb(ji)
        zhe        = max(zero,sign(one,2.0*he-hth))
        heshth     = he/hth
!
!       Linear surface temperature profile between he=0 and he=0.5*hth
!
!old    gdif       = (1.0-zhe)*heshth+zhe*0.5*(1.0+log(2.0*heshth))
!       gdif       = max(zeps,hakspl*log(hgbqb(ji)/(hakspl*hth)))
!
!       quadratic surface temperature profile between he=0 and he=0.5*hth
!
        gdif       = (1.0-zhe)*heshth*(2.0-heshth)
     &                +zhe*0.5*(1.5+log(2.0*heshth))
!       gdif       = min(gdif,1.5*one)
!
        gdif       = 1.0
        ykn(ji)    = gdif*xkn
        ykg(ji)    = gdif*xkg
!
        zk         = 2.0*((1.0-zignm)*ykn(ji)+2.0*zignm*ykg(ji))
        zdh        = (1.0-zignm)*hnbqb(ji)+
     &               zignm*((1.0+3.0*zigm)*hgbqb(ji)
     &               +4.0*ykg(ji)/ykn(ji)*hnbqb(ji))
        zksdh(ji)  = zk/zdh
        ztbq(ji)   = (1.0-zignm)*tbqb(ji,1)
     &               +zignm*(tbqb(ji,2)*(1.0-zigm)
     &               +tfub(ji)*zigm)
        enddo
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  4) Calculate zqmax, fstbqb, qstbqb AND zab.                         |
!-----------------------------------------------------------------------
!
     
      do ji=kideb,kiut
!
        zhg        = max(zero,sign(one,-hnbqb(ji)))
        zigm       = max(zero,sign(one,hgdif-hgbqb(ji)))
!
!                       GRENFELL AND MAYKUT, 1977.
!
        zst        = zhg*
     &               ((0.18+0.82*max(zero,1.0-(hgbqb(ji)/0.1)))
     &               *(1.0-cldqb(ji))+
     &               (0.35+0.65*max(zero,1.0-(hgbqb(ji)/0.1)))
     &               *cldqb(ji))
        zqmax(ji)  = max(zero,0.5*xlg*(hgbqb(ji)-hgmin))
        zexpo      = min(one,exp(-1.5*(hgbqb(ji)-0.1)))
        fstbqb(ji) = zst*fsolgb(ji)*zexpo
        zfsab      = zst*fsolgb(ji)*(1.0-zexpo)
        zhq        = zigm+(1.0-zigm)*max(zero,
     &               sign(one,qstbqb(ji)-zqmax(ji)))
        qstbqb(ji) = (qstbqb(ji)+(1.0-zhq)*zfsab*ddtb)*swiqst
        zab(ji)    = 1.0-zst*
     &               (zhq*zexpo+(1.0-zhq)*(swiqst+(1.0-swiqst)*zexpo))

        if (icoupl .eq. 0) then
          rhoa       = psbqb(ji)/(287.04*tabqb(ji))
          ce         = 1.75e-03
          zrchu(ji)  = rhoa*ce*vabqb(ji)
        else
          zrchu(ji)  = 0.0
        endif


        enddo
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  5) Calculate surface temperature (Newton-Raphson method).           |
!-----------------------------------------------------------------------
!
!NCAR: Newton-Raphson surface temperature.
!
      surthi = 0.10
      do ji=kideb,kiut
!       es         =  611.0*10.0**(9.5*(tsb(ji)-273.16)
!    &              /(tsb(ji)-7.66))
!       zqsat(ji)  =  (0.622*es)/(psbqb(ji)-(1.0-0.622)*es)
!       fcsb(ji)   =  zrchu(ji)*1004.0*(tsb(ji)-tabqb(ji))
!       fleb(ji)   =  zrchu(ji)*2.834e+06*(zqsat(ji)-qabqb(ji))
        rcpdt(ji)  =  (rcpn*min(hnbqb(ji),surthi)+
     &                 rcpg*max(surthi-hnbqb(ji),zero) )/ddtb
        ztst(ji)   =  tsb(ji)
     
      enddo
      
      do jt=1,nbits
        do ji=kideb,kiut

#if ( I_COUPL == 0 )

!!          if (icoupl .eq. 0) then
          fratsb(ji) =  zemin(ji)*(ratbqb(ji)-
     &                  stefan*tsb(ji)*tsb(ji)*tsb(ji)*tsb(ji))
          es         =  611.0*10.0**(9.5*(tsb(ji)-273.16)
     &                  /(tsb(ji)-7.66))
          zqsat(ji)  =  (0.622*es)/(psbqb(ji)-(1.0-0.622)*es)
          fcsb(ji)   =  zrchu(ji)*1004.0*(tsb(ji)-tabqb(ji))
          fleb(ji)   =  zrchu(ji)*2.834e+06*(zqsat(ji)-qabqb(ji))
          zssdqw     =  (zqsat(ji)*zqsat(ji)*psbqb(ji)/
     &                  (0.622*es))*alog(10.0)*
     &                  9.5*((273.16-7.66)/(tsb(ji)-7.66)**2)
          dzf        =  4.0*zemin(ji)*stefan*tsb(ji)
     &                  *tsb(ji)*tsb(ji)+
     &                  zrchu(ji)*(1004.0+2.834e+06*zssdqw)+
     &                  zksdh(ji)
     &                  +rcpdt(ji)
!!         else

#else
          fratsb(ji) =  ratbqb(ji)-(zemin(ji)*
     &                  stefan*tsb(ji)*tsb(ji)*tsb(ji)*tsb(ji))
          dzf        =  4.0*zemin(ji)*stefan*tsb(ji)
     &                  *tsb(ji)*tsb(ji)+
     &                  zksdh(ji)
     &                  +rcpdt(ji)+vabqb(ji)
!cp2 &                  -fderb(ji)
         fcsb(ji)    = 0.0
         fleb(ji)    = 0.0

!!          endif
#endif          
!! #endif
          zfs(ji)    =  zksdh(ji)*(ztbq(ji)-tsb(ji))
          zf         = -zab(ji)*fsolgb(ji)-fratsb(ji)
     &                  +fcsb(ji)+fleb(ji)-zfs(ji)
     &                  +rcpdt(ji)*(tsb(ji)-ztst(ji))
          dtsb       = -zf/dzf
!         dtsb       =  dtsb/3.0 
          dtemax     =  10.0
!         dtsb       =  min (max(dtsb/3.0,dtemax*-1.0),dtemax)
          dtsb       =  min (max(dtsb,dtemax*(-1.0)),dtemax)
          tsb(ji)    =  tsb(ji)+dtsb
!
!END NCAR: Newton-Raphson surface temperature.
!
       enddo
       enddo       

!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  6) Limitation of surface temperature and update of sensible         |
!     and latent heat fluxes.                                          |
!-----------------------------------------------------------------------
!
      do ji=kideb,kiut
        tsb(ji)    = min(ztfs(ji),tsb(ji))
!
!NCAR: Update surface fluxes.
!
#if ( I_COUPL == 0 )
!!       if (icoupl .eq. 0) then
        fratsb(ji) = zemin(ji)*(ratbqb(ji)-stefan*tsb(ji)*tsb(ji)
     &               *tsb(ji)*tsb(ji))
        es         = 611.0*10.0**(9.5*(tsb(ji)-273.16)/(tsb(ji)-7.66))
        zqsat(ji)  = (0.622*es)/(psbqb(ji)-(1.0-0.622)*es)
        fcsb(ji)   = zrchu(ji)*1004.0*(tsb(ji)-tabqb(ji))
        fleb(ji)   = zrchu(ji)*2.834e+06*(zqsat(ji)-qabqb(ji))
!!       else
#else
        fratsb(ji) = ratbqb(ji)-(zemin(ji)*stefan*tsb(ji)*tsb(ji)
     &               *tsb(ji)*tsb(ji))

!!       endif
#endif       
        zfs(ji)    = zksdh(ji)*(ztbq(ji)-tsb(ji))
!
!NCAR: Update surface fluxes.
!
      enddo
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  7)Calculate available heat for surface ablation.                    |
!-----------------------------------------------------------------------
!
       do ji=kideb,kiut
        zffs(ji)= fratsb(ji)-fleb(ji)-fcsb(ji)
     &            +zfs(ji)+zab(ji)*fsolgb(ji)
        zffs(ji)= max(zero,zffs(ji))
        zffs(ji)= zffs(ji)*max(zero,sign(one,tsb(ji)-ztfs(ji)))
       enddo
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  8) Calculate changes in internal temperature due to                 |
!     vertical diffusion processes.                                    |
!-----------------------------------------------------------------------
!
      do ji=kideb,kiut
!
!--8.1 Calculate intermediate variables.
!----------------------------------------
!
        umb        = 1.0-beta
        zind1      = max(zero,sign(one,hndif-hnbqb(ji)))
        zindn(ji)  = 1.0-zind1
        zind3      = max(zero,sign(one,hgdif-hgbqb(ji)))
        zindg(ji)  = zind1*zind3
        zks        = 2.0*ddtb*ykn(ji)/rcpn
        zkg        = 4.0*ddtb*ykg(ji)*(1.0-zindg(ji))/
     &               max(hgbqb(ji)*hgbqb(ji)*rcpg,zeps)
        zd         = ykn(ji)*hgbqb(ji)+2.0*ykg(ji)*hnbqb(ji)*zindn(ji)
        zaux1      = 2.0*zindn(ji)*zks*ykg(ji)/zd
        zaux2      = 2.0*zindn(ji)*zkg*ykn(ji)*hgbqb(ji)/zd
        zfsb       = ddtb*zfs(ji)
!
!--8.2. Fulfill the linear system matrix.
!-----------------------------------------
!
        zkgb       =  zkg*beta
        za1        = -zaux1*beta
        za2        = -zkgb
        za3        =  0.0
        zc1        =  0.0
        zc2        = -zaux2*beta
        zc3        = -zkgb
        zb1        =  zind1+hnbqb(ji)*zindn(ji)-za1
        zb2        =  1.0+zkgb-zc2
        zb3        =  1.0+3.0*zkgb
!
!--8.3. Fulfill the independent term vector.
!-------------------------------------------
!
        zd1        = hnbqb(ji)*zindn(ji)*tbqb(ji,1)
     &               -zindn(ji)*zfsb/rcpn+
     &               zind1*tsb(ji)+umb*zaux1*(tbqb(ji,2)-tbqb(ji,1))
        zd2        = tbqb(ji,2)-
     &               2.0*zind1*zfsb*(1.0-zindg(ji))/
     &               max(hgbqb(ji)*rcpg,zeps)+
     &               umb*(zaux2*(tbqb(ji,1)-tbqb(ji,2))+
     &                    zkg*(tbqb(ji,3)-tbqb(ji,2)))
        zd3        = tbqb(ji,3)+
     &               zkg*(2.0*tfub(ji)+umb*(tbqb(ji,2)-3.0*tbqb(ji,3)))
!
!--8.4. Solve the system (Gauss elimination method).
!----------------------------------------------------
!
        za1        = -za1/zb1
        zd1        =  zd1/zb1
        za2        = -za2/(zb2+zc2*za1)
        zd2        =  (zd2-zc2*zd1)/(zb2+zc2*za1)
        tbqb(ji,3) =  (zd3-zc3*zd2)/(zb3+zc3*za2)
        tbqb(ji,2) =  za2*tbqb(ji,3)+zd2
        tbqb(ji,1) =  za1*tbqb(ji,2)+zd1
        ztci       =  tbqb(ji,2)*(1.0-zindg(ji))+tfub(ji)*zindg(ji)
        zti        =
     &  (4.0*ykg(ji)*hnbqb(ji)*ztci+ykn(ji)*hgbqb(ji)
     &  *tsb(ji)*(1.0+3.0*zindg(ji)))/
     &  (4.0*ykg(ji)*hnbqb(ji)+ykn(ji)*hgbqb(ji)*(1.0+3.0*zindg(ji)))
        tbqb(ji,2) =  tbqb(ji,2)+zindg(ji)
     &                *((3.0*zti+tfub(ji))/4.0-tbqb(ji,2))
        tbqb(ji,3) =  tbqb(ji,3)+zindg(ji)
     &                *((zti+3.0*tfub(ji))/4.0-tbqb(ji,3))
        tbqb(ji,1) =  tbqb(ji,1)+zind1*(zti-tsb(ji))/2.0
!
       enddo
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  9.Take into account surface ablation and bottom accretion-abalation.|
!----------------------------------------------------------------------
!
      do ji=kideb,kiut
!
!--9.1. Surface ablation and update of snow thickness and qstbqb
!------------------------------------------------------------------
!
!  surface ablation
!
        zdhsm      = -(zffs(ji)*ddtb-qctbqb(ji,1))/xln
        zab(ji)    =  hnbqb(ji)
        zhf        =  hnbqb(ji)+hnpbqb(ji)+zdhsm
        zhn        =  1.0-max(zero,sign(one,-zhf))
        hnbqb(ji)  =  max(zero,zhf)
        dvsbqb(ji) =  (1.0-albqb(ji))*(hnbqb(ji)
     &                -zab(ji)-hnpbqb(ji))
        dvsbqb(ji) =  min(zero,dvsbqb(ji))
        dmnbqb(ji) =  rhon*dvsbqb(ji)
        zqfont     =  max(zero,-zhf)*xln
        zdhg(ji)   = -zqfont/xlg
        zqf        =  qstbqb(ji)+zdhg(ji)*xlg
        ziqf       =  max(zero,sign(one,zqf))
        zhq        =  max(zero,sign(one,qstbqb(ji)-zqmax(ji)))
        zdhg(ji)   =  zdhg(ji)+zhq*(ziqf*zdhg(ji)
     &              -(1.0-ziqf)*qstbqb(ji)/xlg)
        qstbqb(ji) =  qstbqb(ji)+zhq*(ziqf*zqf-qstbqb(ji))
        dvsbqb(ji) =  zdhg(ji)*(1.0-albqb(ji))
!
! If fleb is negative, snow condensates at the surface.
!
        zdhcf      =  hnbqb(ji)-parsub*fleb(ji)/(rhon*xsn)*ddtb
        hnbqb(ji)  =  max(zero,zdhcf)
        zhn        =  1.0-max(zero,sign(one,-hnbqb(ji)))
        zdhg(ji)   =  zdhg(ji)-max(zero,-zdhcf)*rhon/rhog
!
!  Internal temperature and qstbqb.
!
!
        tbqb(ji,1) = zhn*tbqb(ji,1)+(1.0-zhn)*tfub(ji)
        zhnf       = max(zero,sign(one,hnbqb(ji)-zab(ji)))
        tbqb(ji,1) = tbqb(ji,1)+
     &               (1.0-zab(ji)/max(hnbqb(ji),zeps))*
     &               (tsb(ji)-tbqb(ji,1))*zindn(ji)*zhnf
        zqnes      = (tfsg-tbqb(ji,2))*rcpg*(hgbqb(ji)/2.0)
        zqr        = qstbqb(ji)-zqnes
        ziqr       = max(zero,sign(one,zqr))
        tbqb(ji,2) = ziqr*tfsg+(1-ziqr)*(tbqb(ji,2)+
     &               qstbqb(ji)/(rcpg*(hgbqb(ji)/2)))
        qstbqb(ji) = ziqr*zqr*swiqst
!
!--9.2. Calculate bottom accretion-ablation and update qstbqb.
!--------------------------------------------------------------
!
        ztb        = tbqb(ji,3)*(1.0-zindg(ji))+tsb(ji)*zindg(ji)
        zf2        = 4.0*ykg(ji)*(tfub(ji)-ztb)/
     &               (hgbqb(ji)+zindg(ji)*(3.*hgbqb(ji)
     &                +4.0*ykg(ji)/ykn(ji)*hnbqb(ji)))
        zdhb(ji)   = ddtb*(zf2-fbbqb(ji)-qlbbqb(ji))/xlg
        zqf        = qstbqb(ji)+zdhb(ji)*xlg
        ziqf       = max(zero,sign(one,zqf))
        zihq       = max(zero,sign(one,qstbqb(ji)-zqmax(ji)))
        zidhb      = max(zero,sign(one,-zdhb(ji)))
        zdhb(ji)   = zdhb(ji)+zihq*zidhb*(ziqf*zdhb(ji)-
     &               (1.0-ziqf)*qstbqb(ji)/xlg)
     &               -qctbqb(ji,2)/xlg
        qstbqb(ji) = qstbqb(ji)+zihq*zidhb*(ziqf*zqf-qstbqb(ji))
        zdhbf      = max(hmelt,zdhb(ji))
        zfsup(ji)  = xlg*(1.0-albqb(ji))/albqb(ji)
     &               *(zdhbf-zdhb(ji))/ddtb
        zdhb(ji)   = zdhbf
        zhgnew     = hgbqb(ji)+zdhg(ji)+zdhb(ji)
        dvbbqb(ji) = (1.0-albqb(ji))*zdhb(ji)
!
!--9.3. Case of total ablation.
!--------------------------------
!
        zhni       =  hnbqb(ji)
        zihgnew    =  1.0-max(zero,sign(one,-zhgnew))
        zihg       =  max(zero,sign(one,-zhni))
        zdhgm      =  (1.0-zihgnew)*(zhgnew-qstbqb(ji)/xlg)
        zdhnm      =  (1.0-zihg)*zdhgm*rhog/rhon
        zhgnew     =  zihgnew*zhgnew
        dmgbqb(ji) =  dmgbqb(ji)+(1.0-albqb(ji))
     &                *(zhgnew-hgbqb(ji))*rhog
        zhnfi      =  zhni+zdhnm
        hnbqb(ji)  =  max(zero,zhnfi)
        dmnbqb(ji) =  dmnbqb(ji)+(1.0-albqb(ji))*(hnbqb(ji)-zhni)*rhon
        zfocea(ji) = -(zihg*zdhgm*xlg+(zhnfi-hnbqb(ji))*xln)/ddtb
        qstbqb(ji) =  zihgnew*qstbqb(ji)
!
!--9.4. Update internal temperature and ice thickness.
!-------------------------------------------------------
!
        tsb(ji)    =  zihgnew*tsb(ji)+(1.0-zihgnew)*tfub(ji)
        zidhb      =  max(zero,sign(one,-zdhb(ji)))
        zc1        = -zhgnew*0.5
        zpc1       =  min(0.5*one,-hgbqb(ji)*0.5-zdhg(ji))
        zc2        = -zhgnew
        zpc2       =  zidhb*zc2+(1.0-zidhb)*(-hgbqb(ji)-zdhg(ji))
        zp1        =  max(zpc1,zc1)
        zp2        =  max(zpc2,zc1)
        zep(ji)    =  tbqb(ji,2)
        tbqb(ji,2) =
     &            2.0*(-zp1*tbqb(ji,2)+(zp1-zp2)*tbqb(ji,3)
     &            +(zp2-zc1)*tfub(ji))*
     &            zihgnew/max(zhgnew,zeps)+(1.0-zihgnew)*tfub(ji)
        zp1        =  min(zpc1,zc1)
        zp2        =  min(zpc2,zc1)
        zp1        =  max(zc2,zp1)
        tbqb(ji,3) =  2.0*
     &                ((1.0-zidhb)*((zc1-zp2)*tbqb(ji,3)
     &                +(zp2-zc2)*tfub(ji))+
     &                 zidhb*((zc1-zp1)*zep(ji)
     &                 +(zp1-zc2)*tbqb(ji,3)))*
     &                zihgnew/max(zhgnew,zeps)+(1.0-zihgnew)*tfub(ji)
        hgbqb(ji)  =  zhgnew
!
       enddo
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  10) Surface accretion.                                              |
!-----------------------------------------------------------------------
!
      do ji=kideb,kiut
!
!  Archimedes principle:
!       zrmh       = (rhon*hnbqb(ji)+rhog*hgbqb(ji))/rho0
!
!  Lepparanta (1983):
!       dzrmh      = max(zero,(rhon*hnbqb(ji)+(rhog-rho0)
!    &               *hgbqb(ji))/(rhon+rho0-rhog))
        dzrmh      = max(zero,(rhon*hnbqb(ji)+(rhog-rho0)
     &               *hgbqb(ji))/(alphs*rhon+rho0-rhog))
!
!  New ice thickness.
!       zhgnew     = max(hgbqb(ji),zrmh)
!       zhnnew     = (rho0*zrmh-rhog*zhgnew)*rhoesn
!
! Lepparanta (1983):
!
        zhgnew     = max(hgbqb(ji),hgbqb(ji)+dzrmh)
!       zhnnew     = min(hnbqb(ji),hnbqb(ji)-dzrmh)
        zhnnew     = min(hnbqb(ji),hnbqb(ji)-alphs*dzrmh)
!
        zig        = 1.0-max(zero,sign(one,-zhgnew))
!
!  Compute new ice temperatures. snow temperature remains unchanged.
!
        quot       = (1.0-zig)+zig*min(one,hgbqb(ji)/max(zhgnew,zeps))
!       tneq       = cnscg*tbqb(ji,1)
!
!  Lepparanta (1983):
!
!       tneq       = cnscg*tbqb(ji,1)+(1.0-rhon/rhog)*tfub(ji)
        tneq       = alphs*cnscg*tbqb(ji,1)+
     &                   (1.0-alphs*(rhon/rhog))*tfub(ji)
!
        zep(ji)    = tbqb(ji,2)
!
!  Lepparanta (1983) (latent heat released during white ice formation
!  goes to the ocean -for lateral ablation-)
!       qlbqb(ji)  = qlbqb(ji)+dzrmh*(1.0-rhon/rhog)*xlg*(1.0-albqb(ji))
        qlbqb(ji)  = qlbqb(ji)+dzrmh*(1.0-alphs*(rhon/rhog))
     &                                      *xlg*(1.0-albqb(ji))
!
        tbqb(ji,2) = tneq-quot*quot*(tneq-tbqb(ji,2))
        tbqb(ji,3) = 2.0*tneq+quot*(tbqb(ji,3)
     &                   +zep(ji)-2.0*tneq)-tbqb(ji,2)
!
!  Changes in ice volume and ice mass.
!
        dvnbqb(ji) = (1.0-albqb(ji))*(zhgnew-hgbqb(ji))
        dmgwib(ji) = dmgwib(ji)+(1.0-albqb(ji))
     &                *(hnbqb(ji)-zhnnew)*rhon
!
!  Lepparanta (1983):
!
        dmgbqb(ji) = dmgbqb(ji)+(1.0-albqb(ji))
     &              *(zhgnew-hgbqb(ji))*rhog
        dmnbqb(ji) = dmnbqb(ji)
     &              -(1.0-albqb(ji))*(zhgnew-hgbqb(ji))*alphs*rhon
!    &              -(1.0-albqb(ji))*(zhgnew-hgbqb(ji))*rhon
!fd0    dmgbqb(ji) = dmgbqb(ji)+(1.0-albqb(ji))
!fd0 &              *(zhgnew-hgbqb(ji))*(rhog-alphs*rhon)
!    &              *(zhgnew-hgbqb(ji))*(rhog-rhon)
!
!  Actualize new snow and ice thickness.
!
        hnbqb(ji)  = zhnnew
        hgbqb(ji)  = zhgnew
      enddo
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  11. Lateral ablation.                                               |
!-----------------------------------------------------------------------
!
      do ji=kideb,kiut
        zab(ji)    = albqb(ji)
        zia        = 1.0-max(zero,sign(one,-hgbqb(ji)))
        zian       = 1.0-max(zero,sign(one,-hnbqb(ji)))
        albqb(ji)  = (1.0-zia)+zia*zab(ji)
        fscbqb(ji) = (1.0-zab(ji))*(1.0-thcmb(ji))*fstbqb(ji)
        qfvbqb(ji) = (zab(ji)*zfsup(ji)+(1.0-zia)*(1.0-zab(ji))*
     &               zfocea(ji))*ddtb
        qlbqb(ji)  = qlbqb(ji)+qfvbqb(ji)+(1.0-zia)*fscbqb(ji)*ddtb
        zq         = hnbqb(ji)*xln+hgbqb(ji)*xlg
        zqtot      = (1.0-albqb(ji))*zq
        ziq        = 1.0-max(zero,sign(one,qlbqb(ji)-zqtot))
        albqb(ji)  = ziq*(albqb(ji)+max(zero,qlbqb(ji)/max(zeps,zq)))+
     &               (1.0-ziq)
        fltbqb(ji) = ((1.0-zab(ji))*qstbqb(ji)-zqtot)/ddtb
!
!  Opening of leads: HAKKINEN & MELLOR.
!
        zalbqb(ji) = (albqb(ji)+max(zero,(-(zdhg(ji)+zdhb(ji))*hakspl*
     &           (1.0-zab(ji)))/max(zeps0,hgbqb(ji)
     &            +hnbqb(ji)*rhon/rhog)))
!
!  Opening of leads: HIbLER.
!
!       zalbqb(ji) = (albqb(ji)+max(zero,
!    &               ((-zdhg(ji)-zdhb(ji))*hibspl)/
!    &               max(zeps,hgbqb(ji)+hnbqb(ji)*rhon/rhog)))
!
!  Opening of leads: OLD.
!
!       zalbqb(ji) = albqb(ji)
!
        zalbqb(ji) = ziq*min(0.99*one,zalbqb(ji))+(1-ziq)
!
        tsb(ji)    = tsb(ji)+(1.0-ziq)*(tfub(ji)-tsb(ji))
        tbqb(ji,1) = tbqb(ji,1)+(1.0-ziq)*(tfub(ji)-tbqb(ji,1))
        tbqb(ji,2) = tbqb(ji,2)+(1.0-ziq)*(tfub(ji)-tbqb(ji,2))
        tbqb(ji,3) = tbqb(ji,3)+(1.0-ziq)*(tfub(ji)-tbqb(ji,3))
        dvlbqb(ji) = zia*(zab(ji)-albqb(ji))*hgbqb(ji)
        dmgbqb(ji) = dmgbqb(ji)+dvlbqb(ji)*rhog
        dvn        = zian*(zab(ji)-albqb(ji))*hnbqb(ji)
        dmnbqb(ji) = dmnbqb(ji)+dvn*rhon
        hnbqb(ji)  = ziq*hnbqb(ji)
        hnbqb(ji)  = ziq*hnbqb(ji)*(1.0-albqb(ji))
     &               /max(zeps,1.0-zalbqb(ji))
        hgbqb(ji)  = ziq*hgbqb(ji)*(1.0-albqb(ji))
     &               /max(zeps,1.0-zalbqb(ji))
        qstbqb(ji) = ziq*(1.0-zab(ji))*qstbqb(ji)
     &               /max(zeps,1.0-zalbqb(ji))
        albqb(ji)  = zalbqb(ji)
      enddo
!
      return
      end
