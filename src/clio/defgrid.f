!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009

      SUBROUTINE defgrid(nflag)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! initatilisation of the masque and of the array linked with the grid.
! If nflag=2 : write testing varaibles on the file "mouchard" ;
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  modif : 23/08/02

!! START_OF_USE_SECTION

#if ( SHELFMELT == 1 )
      ! afq -- a local land-ocean mask for the ISM:
      use varsClio2ism_mod, only : tmsism
#endif

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod,  only: cstmin, cstmax, radian, rterre, omega
     &                    , pi, zero, gpes, one, rho0, cpo, unsrt

      use para0_mod,  only: nsmax, ltest, jsepar, ncomax, nlpmax
     &                    , imax, jmax, kmax
      use para_mod,   only: nbsmax
      use bloc0_mod,  only: xslon, yslat, phifu, phifv, kniv, unshu
     &                    , hu, hs, hux, huy, xulon, yulat, xsedg, ysedg
     &                    , yulat, angle, yslatp, xslonp, xuedg, yuedg
     &                    , unsdzw, dzw, z, zw, scs, q2turb, scpme, phivs
     &                    , ims1, ims2, tms, tmu, ks1, ks2, js1, js2
     &                    , scal, scalr, ijudl, ijsdl, jdl1, jdl2, jeq
     &                    , ku1, ku2, dx, dy, unsdx, unsdy, jberpm
     &                    , unsdx, unsdy, jcl1, jcl2, iberp, ibera
     &                    , jbera, jberp, iberam, jberam, xulonp, yulatp
     &                    , is1, is2, unsdx, unsdy, iu1, iu2, iuf1, iuf2
     &                    , ju1, ju2, imu1, imu2, kfs, kfu, dtb, zurfow
     &                    , dts, dtu, isf1, isf2, afilt, slopmgm, aitd
     &                    , slopemax, avnu0, avnub, avv, ahh, ai, ahs
     &                    , alphgr, alphah, rifsmx, rifumx, msks, alphmi
     &                    , alphaz, msku, algrmn, avk0, avkb, ahe, ahu
     &                    , bering, unsdz, dz, iberpm, fss, phifs, phiss
     &                    , rappes, phimnx, bvf, avsdz, avudz, fqajc
     &                    , rappel, w, scal0, uns2dx, uns2dy
      use bloc_mod,   only: cmx, cmy, smx, smy, cmxy, smxy
     &                    , n1coin, n2coin, n3coin, n4coin
     &                    , i1coin, i2coin, i3coin, i4coin
     &                    , kslp, lslp, ijslp, nxslp, aire
     &                    , unsvol, ctmi, cmxdy, cmydx, covrai
     &                    , xang1, xang2, fub, fvb, phihhh, phizzz
     &                    , nrap, q, q2tmin, fs2cor, kfond, nxyslp
      use ice_mod,    only: tmics, tauc, reslum, areicebs, iicebern1
     &                    , iicebern2, jicebern1, jicebern2, areicebn
     &                    , iicebers1, iicebers2, jicebers1, jicebers2
      use dynami_mod, only: alambd, akappa, bkappa, area, dfhu, dfhv
     &                    , dxs1, gridsz, dxc1, dxc2, zindfa, ren
     &                    , dxs2, zfn, uvdif, wght, npo1i0, ipo1i0
     &                    , jpo1i0, zepsd1
      use reper_mod,  only: kbath, kbath2, iezon, iszon, jehsf, iehsf
     &                    , ishsf, jshsf, tithsf, nvhsf, ndhsf, xwpoln
     &                    , icheck, jcheck, kcheck, dxwi, dywj
     &                    , dlong, dlat, ylat1, xlon1, xaj1, yai1
     &                    , dxaj, dyai, jsep, kbath1, jnorth, xwi1
     &                    , ywj1

#if ( BATHY >= 2 || NC_BERG == 2 )
      use comemic_mod, only: iyear
      use palaeo_timer_mod, only: palaeo_year
      use update_clio_bathy_tools, only: la_date
#endif

#if ( BRINES >= 3 )
      use comemic_mod, only: iyear
      use palaeo_timer_mod, only: palaeo_year
      use brines_mod, only: la_date_brines
#endif

      use coord_mod, only: zlatt, zlont
      use newunit_clio_mod, only: clio3_out_id, mouchard_id

!! END_OF_USE_SECTION


      implicit none

!! START_OF_INCLUDE_SECTION

!! #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "ice.com"
! [SCRPTCOM] #include "reper.com"
! [SCRPTCOM] #include "dynami.com"

!! END_OF_INCLUDE_SECTION

!     include 'met.inc.com'

! ntypcn    Water type
! rrro      --
! dd1o       | Parameter for the abosorption of solar radiation in the ocean
! dd2o      __
! dx1       Relative length of a grid square (direction x). Gen. in radian
! dx2       Relative length of a grid square (direction y). Gen. in radian
! h1        Metric coefficient in the direction x (defined at the center
!           of the grid )
! h2        Metric coefficient in the direction y (defined at the center
!           of the grid )
! d2d1      Derivatite of h2 in the x direction (defined at the center)
! d1d2      Derivatite of h1 in the y direction (defined at the center)
! h1p       Idem h1 for the bottom left corner of the grid
! h2p       Idem h2 for the bottom left corner of the grid
! d2d1p     Idem d2d1 for the bottom left corner of the grid
! d1d2p     Idem d1d2 for the bottom left corner of the grid
! h1pp      Idem h1 for the middle up of the grid
! h2pp      Idem h2 for the middle up of the grid
!

!--- By reference variables
      integer(ip), intent(in) :: nflag

!--local variables :
      real(dblp), dimension(kmax) :: zbath, dzbath, unszw

      real(dblp), dimension(jmax) :: csxsa, csxua
!     dimension snxua(jmax)
      real(dblp), dimension(5)    :: zrr,zd1, zd2
      real(dblp), dimension(20)   :: tauco
!      real(dblp), dimension(imax,jmax) :: zlatt, zlont
      integer, dimension(10) :: ipt0v, jpt0v
      integer, dimension(imax, jmax) :: ntypcn
      real(dblp), dimension(imax,jmax) :: rrro
      real(dblp), dimension(imax,jmax) :: dd1o, dd2o

      real(dblp), dimension(imax,jmax) :: dx1, dx2,
     &          h1, h2, h1p, h2p, d1d2p, d2d1p, d1d2, d2d1
     &         ,h1pp, h2pp

!ph
!ph jsep2   aux. array to mask unneeded regions in rotated grid
!ph not valid CL30, only for CL15
      integer jsep2(83:imax),isep2(29:jmax)
!ph
      character(len=1)  ::  cc1(0:9)
      character(len=8)  ::  fmt1
      character(len=30) :: fmt, fmtrep

!nb
#if ( BATHY >= 2 )
      character*30 name_file
      character(len=5) :: charI
!!!      integer :: la_date
      real, dimension(imax,jmax,kmax) :: tms_old, tmu_old
#endif

!nb debut
      real(dblp), dimension(imax,jmax), save :: fluxbrines
c~       common / common_brines / fluxbrines
!nb fin
cnb grid
#if ( BATHY >=1 )
      real, dimension(4) :: test_wet
      real :: kbath_bas
#endif
!      common / coord /  zlont,zlatt

      data tauco  /6.6,6.6,7.0,7.2,7.1,6.8,6.5,6.6,7.1,7.6,
     &             6.6,6.1,5.6,5.5,5.8,5.8,5.6,5.6,5.6,5.6/
!    &             6.6,6.1,5.6,5.5,5.8,6.1,6.2,6.0,5.6,5.6/
!ph not valid CL30, only for CL15
      data jsep2 / 49, 50, 50, 50, 50, 49, 47, 46, 46, 47, 47, 47,
     &             56, 56, 56, 56, 63, 64, 65, 65, 65, 65, 65, 65,
     &             65, 64, 64, 61, 61, 60, 59, 46, 46, 46, 46, 46,
     &             46, 46, 46, 29 /

cnb isep2 not used: to be suppressed???
!      data isep2 /116,109,108,107,107,107,107,108,109,109,110,113,
!     &            118,120,121,121,110,109,112,112,110,110,110,110,
!     &            111,113,113,113,113,112,111,109,109,109,107,107,
!     &            105 /

!      data jsep2 / 49, 50, 50, 50, 50, 49, 47, 46, 46, 47, 47, 47,
!     &             47, 56, 56, 55, 55, 64, 65, 65, 65, 65, 64, 64,
!     &             63, 63, 60, 60, 59, 58, 45, 45, 45, 45, 45, 45,
!     &             45, 45, 45, 29 /
!ph

       integer(ip) :: kmaxp1, i, j, k, n, ns, nn, npt0v, ipm, jpm, im1
     >             , jm1, ip1, km, kkm, nn0vit, km1, nl, ij, nxslp0
     >             , nyslp0, llm, jj1, ii, jj2, ii2, ii1, nfrc, nv, jp1
     >             , indx, klim, indxp1, km3, kkm3
       real(dblp)  :: ysdeg, yudeg, yurad, yu0rad, ccsysw, sscysw, szphi
     >              , zlam, xsrad, ccmydx, sscyua, sscysa, yfac, usden
     >              , ysrad, xx4tms, volz, ccdxdy, ah, clat, alat
     >              , xurad, znivsp, znivin, dr, dl, effecta, ccsysa
     >              , yu2rad, yyalim, ccsyua, ssnyua, xu0m90, xudeg
     >              , xsdeg, ssnysw, ccsyuw, sscyuw, ccmxdy, ccfcor
     >              , ccfcu8, ssnyuw

       integer(ip) :: typeaux_dat_id, lat_dat_id, lon_dat_id, mask_defgrid_id
     >              , bath_om_id, bath_id
       
      kmaxp1 = kmax + 1

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initialisation of local and global variables (common) .         |
!-----------------------------------------------------------------------

      !write(*,*) 'TMS -1 ', tmu(56,3,14)

      do n=1,9
        cc1(n) = '-'
      enddo
      cc1(0) = '|'
      cc1(5) = '5'

      do j=1,jmax
       is1(j) = imax
       is2(j) = 1
       iu1(j) = imax
       iu2(j) = 1
       iuf1(j) = imax
       iuf2(j) = 1
       isf1(j) = imax
       isf2(j) = 1
      enddo 

      do k=1,kmax
        z(k)  = 0.0
        zw(k) = 0.0
        dz(k)  = 0.0
        dzw(k) = 0.0
        unsdz(k)  = 0.0
        unsdzw(k) = 0.0
      enddo
      z(kmax+1) = 0.0
      zw(kmax+1) = 0.0
      dzw(kmax+1) = 0.0
      unsdzw(kmax+1) = 0.0

      do j=1,jmax
       do i=1,imax
        hs(i,j) = 0.0
        huy(i,j) = 0.0
        hux(i,j) = 0.0
        hu(i,j) = 0.0
        unshu(i,j) = 0.0
!       kfs(i,j) = kmax
!       kfu(i,j) = kmax
        kniv(i,j,-1) = kmaxp1
        kniv(i,j, 0) = kmaxp1
        kniv(i,j, 1) = kmaxp1
cnb modified only outside of defgrid
      if (nflag.ne.3) then
        phifu(i,j)  = 0.0
        phifv(i,j)  = 0.0
      endif
cnb not modified outside of defgrid
        xang1(i,j)=0.0
        xang2(i,j)=1.0
!ph
!ph preset t- and u-grid longitudes and latitudes, and the
!ph respective gridcell corners
!ph
        xslon(i,j  )=90.0
        yslat(i,j  )=90.0
        xsedg(i,j,1)=90.0
        xsedg(i,j,2)=90.0
        xsedg(i,j,3)=90.0
        xsedg(i,j,4)=90.0
        ysedg(i,j,1)=90.0
        ysedg(i,j,2)=90.0
        ysedg(i,j,3)=90.0
        ysedg(i,j,4)=90.0
        xulon(i,j  )=90.0
        yulat(i,j  )=90.0
        xuedg(i,j,1)=90.0
        xuedg(i,j,2)=90.0
        xuedg(i,j,3)=90.0
        xuedg(i,j,4)=90.0
        yuedg(i,j,1)=90.0
        yuedg(i,j,2)=90.0
        yuedg(i,j,3)=90.0
        yuedg(i,j,4)=90.0
        angle(i,j)=0.0
!ph
        enddo
      enddo
      do j=1,jmax+1
        do i=1,imax+1
          xslonp(i,j)=90.0
          yslatp(i,j)=90.0
          xslonp(i,j)=90.0
          yslatp(i,j)=90.0
        enddo
      enddo

!dmr&nb --- update_bathy
      if (nflag.ne.3) then

      do ns=0,nsmax
       do j=1,jmax
        do i=1,imax
         fss(i,j,ns)  = 0.0
         phifs(i,j,ns)  = 0.0
         phiss(i,j,ns)  = 0.0
         rappes(i,j,ns) = 0.0
         phimnx(i,j,0,ns) = cstmin
         phimnx(i,j,1,ns) = cstmax
        enddo
       enddo
      enddo

!nb debut
       do j=1,jmax
        do i=1,imax
         fluxbrines(i,j) = 0.0
        enddo
       enddo
!nb fin


      do nn=1,6
       do j=1,jmax
        do i=1,imax
         phihhh(i,j,nn) = 0.0
        enddo
       enddo
      enddo

      do k=1,kmax
       do j=1,jmax
        do i=1,imax
         fub(i,j,k) = 0.0
         fvb(i,j,k) = 0.0
         q(i,j,k) = 0.0
         bvf(i,j,k) = 0.0
         avsdz(i,j,k) = 0.0
         avudz(i,j,k) = 0.0
         fqajc(i,j,k) = 0.0
         tms(i,j,k) = 0.0
         tmu(i,j,k) = 0.0
         rappel(i,j,k) = 0.0
        enddo
       enddo
      enddo 

      do k=1,kmax+1
       do j=1,jmax
        do i=1,imax
         w(i,j,k) = 0.0
         phizzz(i,j,k,1) = 0.0
         phizzz(i,j,k,2) = 0.0
        enddo
       enddo
      enddo

!--Initialisation of the scalairs,  forcings, and of the
!  turbuent kinetic energy
      do ns=1,nsmax
       do k=1,kmax
        nrap(k,ns) = 0
        do j=1,jmax
         do i=1,imax
          scalr(i,j,k,ns)=0.0
          scal(i,j,k,ns) = scal0(k,ns)
          phivs(i,j,k,ns) = 0.0
         enddo
        enddo
       enddo
      enddo

      do k=1,kmax
       do j=1,jmax
        do i=1,imax
         q2turb(i,j,k)=q2tmin
        enddo
       enddo
      enddo

!--Initialisation of the scalars associated with the freshwater flux
      do ns=1,nsmax
       do j=1,jmax
        do i=1,imax
         scs(i,j,ns) = scpme(ns)
        enddo
       enddo
      enddo

!dmr&nb --- update_bathy
      else ! if nflag == 3 -> bathy updating!

#if ( BATHY >= 2 )



       tms_old(:,:,:) = tms(:,:,:)
       tmu_old(:,:,:) = tmu(:,:,:)
       do k=1,kmax
        do j=1,jmax
         do i=1,imax

          tms(i,j,k) = 0.0
          tmu(i,j,k) = 0.0
         enddo
        enddo
       enddo

      !write(*,*) 'TMS 0 tms, tmu_old ', tmu(56,3,14), tmu_old(56,3,14)

#endif

      endif
!dmr&nb --- update_bathy

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Initialosarion of the constants ( resolution, grid ...)          |
!-----------------------------------------------------------------------

!--size (maximum)  of the basin (pts type scalar) :
      ims1 = 2
      ims2 = imax - 1
      js1 = 2
      js2 = jmax - 1
      ks1 = 1
      ks2 = kmax
!--others :
      ku2 = ks2
      ku1 = ks1

!--constantes linked with the resolution / axes /  Geography :
      call geogra(js1, js2, jeq, jdl1, jdl2, ijsdl, ijudl,
     &            iberp, ibera)
             !write(*,*) 'IBERP 1', iberp

      if (ltest.lt.1) then
!--basin type : rectangular
        jcl1 = js2
        jcl2 = js1
      else
!--basin with cyclic boundaries, definition of the linked areas :
        jcl1 = js1
        jcl2 = jeq
      endif

!--initialisation of the constants
      if (ltest.lt.0) then
        dx = 1000.
        dy = 1000.
      else
        dx = dlong * radian * rterre
        dy = dlat  * radian * rterre
      endif
      unsdx = 1.0 / dx
      unsdy = 1.0 / dy
      uns2dx = 0.5 / dx
      uns2dy = 0.5 / dy

!- Warning : Bering Stait closed => (iberp,a)=0,0. But not when living "defgrid"
      jberpm = jberp - 1
      iberpm = iberp - 1
      jberam = jbera - 1
      iberam = ibera - 1

!--initialisation of island with no surface :
!- Spitzberg :
      ipt0v(1) = 109
      jpt0v(1) =  55
!ic0  npt0v =  1
      npt0v =  0

!--point with ocean velocity but no sea-ice velocity
      npo1i0 = 1
      npo1i0 = 2
      ipo1i0(1) = 102
      jpo1i0(1) = 56
      ipo1i0(2) = 111
      jpo1i0(2) = 56


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) For each grid, initialisation of metric coefficients .   |
!-----------------------------------------------------------------------

!-------------------
! Conventions :    |
!-------------------
!  index 1 \ 2 <-> location in the direction x \ y
!    3eme indice <-> location on the grid :
!  0=Centre, 1=Side W x(i-1/2), 2=side S y(j-1/2), 3=corner SW x(i-1/2),y(j-1/2)
!-------------------
! Metric coeficients : cmx \ cmy = metric. coef. . direction x (h1) \ y (h2)
! Derivatives of the metriq. coef. : cmxdy \ cmydx = d(cmx)/dy  \ d(cmy)/dx located.3
!
! Other coeficients built from metric. coef. :
!  cmxy(i,j,0\3) = cmx(i,j,0\3) * cmy(i,j,0\3) = element de surface
!  cmxy(i,j,1\2) = cmx(i,j,1\2) / cmy(i,j,1\2)
!  smx = 1 / cmx ; smy = 1 / cmy ; smxy = 1 / cmxy
!-------------------
! loacl conventions :
!  cs <- cosine ; sc <- 1 / cosine ; sn <- sine [ou 1st letter doubled]
!  x\y_s\u_w\a : related to a
!   longitude (x) \ latitude (y) <-- differs from the convention metric coef.
!   for a scalar point (=Centre) (s) \ velocity (corner) (u) [default = s]
!   on the grid WW (w) \ AA (a)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--3.1 1st grid classical spheric , grid World (WW) :
!-------------------------------------------------------

!--initialisation of the geometric arrays cos,sin, etc ...
      do j=1,jsepar
        ysdeg = ylat1 + dlat * DFLOAT(j-1)
!DFG
        yudeg = ysdeg - 0.5 * dlat
!DFG
        yurad = (ysdeg - 0.5 * dlat) * radian
!-
        if (ltest.lt.0) then
!--Cartesian Grid :
          yu0rad = 45.0 * radian
          ccsysw = 1.0
          sscysw = 1.0
          ccsyuw = 1.0
          sscyuw = 1.0
          ccmxdy = 0.0
          ccfcor = omega * sin(yu0rad)
          ccfcu8 = -0.25 * omega * cos(yu0rad)
        else
!--Spherical grid:
          ccsysw = cos(ysdeg * radian)
          sscysw = 1.0 / ccsysw
          ccsyuw = cos(yurad)
          sscyuw = 1.0 / ccsyuw
          ssnyuw = sin(yurad)
          ssnysw = sin(ysdeg*radian)
          ccmxdy = -ssnyuw * unsrt
          ccfcor = omega * ssnyuw
          ccfcu8 = -0.25 * omega * ccsyuw
!--end of the difference between  cartesian/spheric .
        endif

!--metric : spherical grid long x lat :
        do i=1,imax
          zlatt(i,j)=yurad
!ph
!ph simple lon-lat-grid with longitude interval dlong and
!ph latitudinal width dlat. corners are indexed as follows:
!ph  1) SW  2) SE  3) NE  4) NW  of respective gridpoint
!ph momentum-gridpoint i,j is at SW corner of corresponding
!ph tracer point.
!ph
          yslat(i,j)=ysdeg
          ysedg(i,j,1)=yudeg
          ysedg(i,j,2)=yudeg
          ysedg(i,j,3)=yudeg+dlat
          ysedg(i,j,4)=yudeg+dlat
          yulat(i,j)=yudeg
          yuedg(i,j,1)=ysdeg-dlat
          yuedg(i,j,2)=ysdeg-dlat
          yuedg(i,j,3)=ysdeg
          yuedg(i,j,4)=ysdeg
          xsdeg = xlon1 + dlong * DFLOAT(i-1)
          xurad = (xsdeg - 0.5 * dlong) * radian
          xudeg = xsdeg - 0.5 * dlong
          xslon(i,j)=xsdeg
          xsedg(i,j,1)=xudeg
          xsedg(i,j,2)=xudeg+dlong
          xsedg(i,j,3)=xudeg+dlong
          xsedg(i,j,4)=xudeg
          xulon(i,j)=xudeg
          if ((xsdeg-dlong).le.0.0d0) then
            xuedg(i,j,1)=xsdeg-dlong+360.0d0
          else
            xuedg(i,j,1)=xsdeg-dlong
          endif
          xuedg(i,j,2)=xsdeg
          xuedg(i,j,3)=xsdeg
          if ((xsdeg-dlong).le.0.0d0) then
            xuedg(i,j,4)=xsdeg-dlong+360.0d0
          else
            xuedg(i,j,4)=xsdeg-dlong
          endif

          zlont(i,j)=xurad
!- coef :
          cmx(i,j,0) = ccsysw
          cmx(i,j,1) = ccsysw
          cmx(i,j,2) = ccsyuw
          cmx(i,j,3) = ccsyuw
          cmy(i,j,0) = 1.
          cmy(i,j,1) = 1.
          cmy(i,j,2) = 1.
          cmy(i,j,3) = 1.
!- inverse :
          smx(i,j,0) = sscysw
          smx(i,j,1) = sscysw
          smx(i,j,2) = sscyuw
          smx(i,j,3) = sscyuw
          smy(i,j,0) = 1.
          smy(i,j,1) = 1.
          smy(i,j,2) = 1.
          smy(i,j,3) = 1.
!- product, ratio of 2 coefs :
          cmxy(i,j,0) = ccsysw
          cmxy(i,j,1) = ccsysw
          cmxy(i,j,2) = ccsyuw
          cmxy(i,j,3) = ccsyuw
          smxy(i,j,0) = sscysw
          smxy(i,j,1) = sscysw
          smxy(i,j,2) = sscyuw
          smxy(i,j,3) = sscyuw
!- derivatives :
          cmxdy(i,j) = ccmxdy
          cmydx(i,j) = 0.
!- metric coefficients for sea ice dynamic:
          dx1(i,j)   = dlong*radian
          yu2rad     = (ysdeg + 0.5*dlat)*radian
          dx2(i,j)   = sin(yu2rad) -sin(yurad)
          h1(i,j)    = rterre*ccsysw
          h2(i,j)    = rterre/ccsysw
          d1d2(i,j)  = -rterre*sin(ysdeg*radian)*sscysw
          d2d1(i,j)  = 0.0
        enddo

!ph - rotation of grid (should yield zero for simple lon-lat grid)
        do i=1,imax-1
          angle(i,j)=atan2((yulat(i+1,j)-yulat(i,j)),
     &                     (xulon(i+1,j)-xulon(i,j)))
        enddo
        angle(imax,j)=angle(imax-1,j)

!--Coriolis factor  :
        do i=1,imax
          fs2cor(i,j) = ccfcor
!fcc      fcucor(i,j) = ccfcu8
!fcc      fcvcor(i,j) = 0.0
!- Sine (Geog. Lat.) :
          covrai(i,j) = ssnysw
        enddo
      enddo

      if (ltest.ge.3) then

!--3.3 North Atlantic + Arctic (AA), sherical grid N-S / E-W inversed
!-----------------------------------------------------------------------------

!       xudeg0 = xaj1  - 1.5 * dxaj
        xu0m90 = xaj1  - 1.5 * dxaj - 90.
        do j=jeq,jmax
!         xudeg = xudeg0 + dxaj * DFLOAT(j)
          xudeg = xu0m90 + dxaj * DFLOAT(j)
!- Caution : change of  sign ! csxua = -cos(Long_AA)
!         csxua(j) = -cos(xudeg * radian)
          csxua(j) = sin(xudeg * radian)
!fcc      snxua(j) = cos(xudeg * radian)
          csxsa(j) = sin((xudeg+0.5*dxaj)*radian)
        enddo
        yyalim = 90.0 - abs(dyai)
        do i=1,imax
          ysdeg = yai1 + dyai * DFLOAT(i-1)
          if (abs(ysdeg).ge.yyalim) then
            ccsysa = 1.
            yurad =  0.
            ysrad =  0.
            ccsyua = 1.
            ssnyua = 0.
          else
            ccsysa = cos(ysdeg * radian)
            yurad  = (ysdeg - 0.5 * dyai) * radian
            ysrad  = ysdeg * radian
            ccsyua = cos(yurad)
!fcc        ssnyua = sin(yurad)
          endif
          sscysa = 1.0 / ccsysa
          sscyua = 1.0 / ccsyua
          ccmydx =  sin(yurad) * unsrt
!--Metrique : 2nd spherical grid lat x long :
          do j=jsep(i),jmax
!- coef :
            cmx(i,j,0) = 1.
            cmx(i,j,1) = 1.
            cmx(i,j,2) = 1.
            cmx(i,j,3) = 1.
            cmy(i,j,0) = ccsysa
            cmy(i,j,1) = ccsyua
            cmy(i,j,2) = ccsysa
            cmy(i,j,3) = ccsyua
!- inverse :
            smx(i,j,0) = 1.
            smx(i,j,1) = 1.
            smx(i,j,2) = 1.
            smx(i,j,3) = 1.
            smy(i,j,0) = sscysa
            smy(i,j,1) = sscyua
            smy(i,j,2) = sscysa
            smy(i,j,3) = sscyua

!- product, ratio de 2 coefs :
            cmxy(i,j,0) = ccsysa
            cmxy(i,j,1) = sscyua
            cmxy(i,j,2) = sscysa
            cmxy(i,j,3) = ccsyua
            smxy(i,j,0) = sscysa
            smxy(i,j,1) = ccsyua
            smxy(i,j,2) = ccsysa
            smxy(i,j,3) = sscyua
!- derivatives :
            cmxdy(i,j) = 0.
            cmydx(i,j) = ccmydx
!- metric coefficients for sea ice dynamic:
            dx2(i,j)   = dlat*radian
            yu2rad     = (ysdeg + 0.5*dyai)*radian
            dx1(i,j)   = -sin(yu2rad) + sin(yurad)
            h1(i,j)    = rterre/ccsysa
            h2(i,j)    = rterre*ccsysa
            d1d2(i,j)  = 0.0
            d2d1(i,j)  = +rterre*sin(ysdeg*radian)*sscysa
            xsrad=(xaj1+dxaj*DFLOAT(j-1))*radian
            xurad=(xaj1-1.5*dxaj+dxaj*DFLOAT(j))*radian
            yslat(i,j)=-asin(cos(ysrad)*cos(xsrad))/radian
            yulat(i,j)=-asin(cos(yurad)*cos(xurad))/radian
            zlatt(i,j)=-asin(cos(yurad)*cos(xurad))
            szphi=(cos(yurad)*cos(xurad))
            zlam=atan2(cos(ysrad)*sin(xsrad),sin(ysrad))
            xslon(i,j)=(zlam+pi)/radian+69.0
            zlam=atan2(cos(yurad)*sin(xurad),sin(yurad))
            xulon(i,j)=(zlam+pi)/radian+69.0
            zlont(i,j)=zlam+pi+69.0*pi/180.0
!Cice        xang1(i,j)=(-sin(zphi)*cos(zlam))/max(zepsd1,cos(yurad))
            xang1(i,j)=(szphi*cos(zlam))/max(zepsd1,cos(yurad))
            xang2(i,j)=(sin(zlam))/max(zepsd1,cos(yurad))
          enddo
!--coef " between the two grids 2 grids :
          if (jsep(i).eq.jeq) then
            cmy(i,jeq,2) = 0.5 * ( 1. + cmy(i,jeq,0) )
            cmy(i,jeq,3) = 0.5 * ( 1. + cmy(i,jeq,1) )
            smy(i,jeq,2) = 1. / cmy(i,jeq,2)
            smy(i,jeq,3) = 1. / cmy(i,jeq,3)
            cmxy(i,jeq,2) = smy(i,jeq,2)
            cmxy(i,jeq,3) = cmy(i,jeq,3)
            smxy(i,jeq,2) = cmy(i,jeq,2)
            smxy(i,jeq,3) = smy(i,jeq,3)
            cmydx(i,jeq) = 0.5 * ccmydx
          endif
!--Coriolis factor :
          do j=jsep(i),jmax
            fs2cor(i,j) = omega * ccsyua * csxua(j)
!fcc        fcucor(i,j) = omega * snxua(j) * -0.25
!fcc        fcvcor(i,j) = omega * ssnyua * csxua(j) * 0.25
!- Sine (Geog. Lat) :
            covrai(i,j) = ccsysa * csxsa(j)
          enddo
         enddo ! this is the outer imax from l. 661
!ph
!ph set momentum gridpoints as corners of tracer cells
!ph and vice-versa in rotated grid
!ph
        do i=1,imax
          do j=jsep(i),jmax
            ipm=min(i+1,imax)
            jpm=min(j+1,jmax)
            xsedg(i,j,1)=xulon(i  ,j  )
            xsedg(i,j,2)=xulon(ipm,j  )
            xsedg(i,j,3)=xulon(ipm,jpm)
            xsedg(i,j,4)=xulon(i  ,jpm)
            ysedg(i,j,1)=yulat(i  ,j  )
            ysedg(i,j,2)=yulat(ipm,j  )
            ysedg(i,j,3)=yulat(ipm,jpm)
            ysedg(i,j,4)=yulat(i  ,jpm)
            ipm=max(i-1,1)
            jpm=max(j-1,1)
            xuedg(i,j,1)=xslon(ipm,jpm)
            xuedg(i,j,2)=xslon(i  ,jpm)
            xuedg(i,j,3)=xslon(i  ,j  )
            xuedg(i,j,4)=xslon(ipm,j  )
            yuedg(i,j,1)=yslat(ipm,jpm)
            yuedg(i,j,2)=yslat(i  ,jpm)
            yuedg(i,j,3)=yslat(i  ,j  )
            yuedg(i,j,4)=yslat(ipm,j  )
          enddo
        enddo
!ph
!ph avoid unsteady, huge coordinates in the outer space
!ph between the two grids
!ph
        do i=1,82
!L30
          yfac=0.25/float(jmax-jsep(i))
          jpm=jsep(i)-1
          do j=jsep(i),jmax
            xslon(i,j)=xslon(i,jpm)
            yslat(i,j)=yslat(i,jpm)+yfac*(j-jpm)
            xsedg(i,j,1)=xsedg(i,jpm,1)
            xsedg(i,j,2)=xsedg(i,jpm,2)
            xsedg(i,j,3)=xsedg(i,jpm,3)
            xsedg(i,j,4)=xsedg(i,jpm,4)
            ysedg(i,j,1)=yslat(i,j)-0.5*yfac
            ysedg(i,j,2)=yslat(i,j)-0.5*yfac
            ysedg(i,j,3)=yslat(i,j)+0.5*yfac
            ysedg(i,j,4)=yslat(i,j)+0.5*yfac
            xulon(i,j)=xulon(i,jpm)
            yulat(i,j)=yulat(i,jpm)+yfac*(float(j-jpm)-0.5)
            xuedg(i,j,1)=xuedg(i,jpm,1)
            xuedg(i,j,2)=xuedg(i,jpm,2)
            xuedg(i,j,3)=xuedg(i,jpm,3)
            xuedg(i,j,4)=xuedg(i,jpm,4)
            yuedg(i,j,1)=yulat(i,j)-0.5*yfac
            yuedg(i,j,2)=yulat(i,j)-0.5*yfac
            yuedg(i,j,3)=yulat(i,j)+0.5*yfac
            yuedg(i,j,4)=yulat(i,j)+0.5*yfac
          enddo
        enddo
!ph
!ph same procedure in the arctic ocean region
!ph
#if ( HRCLIO == 0 )
        do i=83,imax
          if (jsep2(i).ne.jmax) then
#elif ( HRCLIO == 1 )
!L15    do i=83,imax
!L30    THIS SHOULD BE MODIFIED FOR 30 LEVELS
        jsep2(imax)=jmax
        do i=imax,imax
          if (jsep2(i).ne.jmax) then
#endif
            yfac=0.25/float(jmax-jsep2(i))
            jpm=jsep2(i)-1
            do j=jsep2(i),jmax
              xslon(i,j)=xslon(i,jpm)
              yslat(i,j)=yslat(i,jpm)+yfac*(j-jpm)
              xsedg(i,j,1)=xsedg(i,jpm,1)
              xsedg(i,j,2)=xsedg(i,jpm,2)
              xsedg(i,j,3)=xsedg(i,jpm,3)
              xsedg(i,j,4)=xsedg(i,jpm,4)
              ysedg(i,j,1)=yslat(i,j)-0.5*yfac
              ysedg(i,j,2)=yslat(i,j)-0.5*yfac
              ysedg(i,j,3)=yslat(i,j)+0.5*yfac
              ysedg(i,j,4)=yslat(i,j)+0.5*yfac
              xulon(i,j)=xulon(i,jpm)
              yulat(i,j)=yulat(i,jpm)+yfac*(float(j-jpm)-0.5)
              xuedg(i,j,1)=xuedg(i,jpm,1)
              xuedg(i,j,2)=xuedg(i,jpm,2)
              xuedg(i,j,3)=xuedg(i,jpm,3)
              xuedg(i,j,4)=xuedg(i,jpm,4)
              yuedg(i,j,1)=yulat(i,j)-0.5*yfac
              yuedg(i,j,2)=yulat(i,j)-0.5*yfac
              yuedg(i,j,3)=yulat(i,j)+0.5*yfac
              yuedg(i,j,4)=yulat(i,j)+0.5*yfac
            enddo
          endif
        enddo

!ph
!ph angle of rotation of grid
!ph
        do i=1,imax-1
          do j=jsep(i),jmax
            angle(i,j)=atan2((yulat(i+1,j)-yulat(i,j)),
     &                       (xulon(i+1,j)-xulon(i,j)))
          enddo
        enddo
        i=imax
        do j=jsep(i),jmax
          angle(i,j)=angle(i-1,j)
        enddo

!-- end of the initialisation of the AA grid .
!-----
      endif
!ph
!ph second set of edges for the 2nd form of the 3-argument shade in
!ph ferret (the four point algorithm)
!ph
      do j=1,jmax
        do i=1,imax
          xslonp(i,j)=xsedg(i,j,1)
          yslatp(i,j)=ysedg(i,j,1)
          xulonp(i,j)=xuedg(i,j,1)
          yulatp(i,j)=yuedg(i,j,1)
        enddo
        xslonp(imax+1,j)=xsedg(imax,j,2)
        yslatp(imax+1,j)=ysedg(imax,j,2)
        xulonp(imax+1,j)=xuedg(imax,j,2)
        yulatp(imax+1,j)=yuedg(imax,j,2)
      enddo
      do i=1,imax
        xslonp(i,jmax+1)=xsedg(i,jmax,4)
        yslatp(i,jmax+1)=ysedg(i,jmax,4)
        xulonp(i,jmax+1)=xuedg(i,jmax,4)
        yulatp(i,jmax+1)=yuedg(i,jmax,4)
      enddo
      xslonp(imax+1,jmax+1)=xsedg(imax,jmax,3)
      yslatp(imax+1,jmax+1)=ysedg(imax,jmax,3)
      xulonp(imax+1,jmax+1)=xuedg(imax,jmax,3)
      yulatp(imax+1,jmax+1)=yuedg(imax,jmax,3)

      if (ltest.eq.3 .and. iberp.ne.0) then
!       smy(iberp-1,jberp,2) = 0.
!       smy(iberp , jberp,2) = 0.
!--Bering : local modification (AA) of the Metric Coefs :
        cmx(ibera, jbera,2) = cmx(iberpm,jberp,2)
        smx(ibera, jbera,2) = smx(iberpm,jberp,2)
        cmy(ibera, jbera,2) = cmy(iberpm,jberp,2)
        smy(ibera, jbera,2) = smy(iberpm,jberp,2)
        cmxy(ibera, jbera,2) = cmxy(iberpm,jberp,2)
        smxy(ibera, jbera,2) = smxy(iberpm,jberp,2)
!-
        cmx(iberam,jbera,2) = cmx(iberp, jberp,2)
        smx(iberam,jbera,2) = smx(iberp, jberp,2)
        smx(iberam,jbera,2) = smx(iberp, jberp,2)
        cmy(iberam,jbera,2) = cmy(iberp, jberp,2)
        smy(iberam,jbera,2) = smy(iberp, jberp,2)
        cmxy(iberam,jbera,2) = cmxy(iberp, jberp,2)
        smxy(iberam,jbera,2) = smxy(iberp, jberp,2)
!--Only for the diagnostics :
        cmxy(ibera,jbera,3) = cmxy(iberp,jberp,3)
        smxy(ibera,jbera,3) = smxy(iberp,jberp,3)
        cmxy(iberp, jberp,0) = cmxy(iberam,jberam,0)
        cmxy(iberpm,jberp,0) = cmxy(ibera, jberam,0)
        cmxy(ibera, jbera,0) = cmxy(iberpm,jberpm,0)
        cmxy(iberam,jbera,0) = cmxy(iberp, jberpm,0)
        smxy(iberp, jberp,0) = smxy(iberam,jberam,0)
        smxy(iberpm,jberp,0) = smxy(ibera, jberam,0)
        smxy(ibera, jbera,0) = smxy(iberpm,jberpm,0)
        smxy(iberam,jbera,0) = smxy(iberp, jberpm,0)
      endif

!--Modification of the metric  coeff. close to Gibraltar.
!     include 'defgrid.gibr.inc'

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Computation of metric coefficients needed in ice dynamics       |
!-----------------------------------------------------------------------
!
!
!--4.1. WEIGHTS FOR COMPUTING INTERPOLATED VALUES AT THE CORNERS
!       OF THE GRID SQUARES.
!----------------------------------------------------------------
!
      do j=2,jmax
        jm1 = j-1
        do i=1,imax
          im1           = (i-1)+(imax-2)*(1/i)
          usden         = 1.0/((dx1(i,j)+dx1(im1,j))
     &                    *(dx2(i,j)+dx2(i,jm1)))
          wght(i,j,1,1) = usden*dx1(i,j)*dx2(i,j)
          wght(i,j,1,2) = usden*dx1(i,j)*dx2(i,jm1)
          wght(i,j,2,2) = usden*dx1(im1,j)*dx2(i,jm1)
          wght(i,j,2,1) = usden*dx1(im1,j)*dx2(i,j)
        enddo
      enddo
!
!--4.2. metric coefficients at the grid squares corners,
!       and mid-sides points, and their derivatives (h1p,  h2p,
!       h1pp, h2pp, d1d2p,d2d1p), corilis coefficient (zfn)
!       and coefficients for the strain rate tensor (akappa).
!----------------------------------------------------------------
!
      do i=1,imax
        h1p(i,1)   = 0.0
        h2p(i,1)   = 1.0e+20
        d1d2p(i,1) = 1.0e+20
        d2d1p(i,1) = 0.0
!
        zfn(i,1)   = fs2cor(i,1)*2
      enddo
!
      do j=1,jmax-1
        do i=1,imax
          im1          = (i-1)+(imax-2)*(1/i)
          ip1          = (i+1)-(imax-2)*(i/imax)
          h1p(i,j+1)   = h1(i,j+1)*wght(i,j+1,2,2)+
     &                   h1(im1,j+1)*wght(i,j+1,1,2)+
     &                   h1(i,j)*wght(i,j+1,2,1)+
     &                   h1(im1,j)*wght(i,j+1,1,1)
          h2p(i,j+1)   = h2(i,j+1)*wght(i,j+1,2,2)+
     &                   h2(im1,j+1)*wght(i,j+1,1,2)+
     &                   h2(i,j)*wght(i,j+1,2,1)+
     &                   h2(im1,j)*wght(i,j+1,1,1)
          d1d2p(i,j+1) = 2.0*
     &                   (dx1(i,j)*(-h1(im1,j)+h1(im1,j+1))+
     &                    dx1(im1,j)*(-h1(i,j)+h1(i,j+1)))/
     &                   ((dx1(i,j)+dx1(im1,j))*(dx2(i,j)+dx2(i,j+1)))
          d2d1p(i,j+1) = 2.0*
     &                   (dx2(i,j+1)*(h2(i,j)-h2(im1,j))+
     &                    dx2(i,j)*(h2(i,j+1)-h2(im1,j+1)))/
     &                   ((dx1(i,j)+dx1(im1,j))*(dx2(i,j)+dx2(i,j+1)))
!
          h1pp(i,j)    = (dx2(i,j)*h1(i,j+1)+dx2(i,j+1)*h1(i,j))/
     &                   (dx2(i,j)+dx2(i,j+1))
          h2pp(i,j)    = (dx1(i,j)*h2(ip1,j)+dx1(ip1,j)*h2(i,j))/
     &                   (dx1(i,j)+dx1(ip1,j))
!
          zfn(i,j+1)   = 2*fs2cor(i,j+1)
        enddo
      enddo
      do i=1,imax
        ip1          = (i+1)-(imax-2)*(i/imax)
        h1pp(i,jmax) = 0.0
        h2pp(i,jmax) = (dx1(i,j)*h2(ip1,j)+dx1(ip1,j)*h2(i,j))/
     &                   (dx1(i,j)+dx1(ip1,j))
      enddo
!
      do j=1,jmax
        do i=1,imax
          ip1             = (i+1)-(imax-2)*(i/imax)
          akappa(i,j,1,1) = 1.0/(2.0*h1(i,j)*dx1(i,j))
          akappa(i,j,1,2) = d1d2(i,j)/(4.0*h1(i,j)*h2(i,j))
          akappa(i,j,2,1) = d2d1(i,j)/(4.0*h1(i,j)*h2(i,j))
          akappa(i,j,2,2) =  1.0/(2.0*h2(i,j)*dx2(i,j))
       enddo
      enddo
!
!--4.3. COEFFICIENTS FOR DIVERGENCE OF THE STRESS TENSOR.
!---------------------------------------------------------
!
      do j=2,jmax
        jm1 = j-1
        do i=1,imax
          im1=(i-1)+imax*(1/i)
          usden               =
     &           1.0/(h1p(i,j)*h2p(i,j)*(dx1(i,j)+dx1(im1,j))
     &           *(dx2(i,j)+dx2(i,jm1)))
          alambd(i,j,2,2,2,1) = usden*2.0*dx2(i,j)*h2(i,jm1)
          alambd(i,j,2,2,2,2) = usden*2.0*dx2(i,jm1)*h2(i,j)
          alambd(i,j,2,2,1,1) = usden*2.0*dx2(i,j)*h2(im1,jm1)
          alambd(i,j,2,2,1,2) = usden*2.0*dx2(i,jm1)*h2(im1,j)
          alambd(i,j,1,1,2,1) = usden*2.0*dx1(im1,j)*h1(i,jm1)
          alambd(i,j,1,1,1,1) = usden*2.0*dx1(i,j)*h1(im1,jm1)
          alambd(i,j,1,1,2,2) = usden*2.0*dx1(im1,j)*h1(i,j)
          alambd(i,j,1,1,1,2) = usden*2.0*dx1(i,j)*h1(im1,j)
          alambd(i,j,1,2,1,1) = usden*d1d2p(i,j)*dx2(i,j)*dx1(i,j)
          alambd(i,j,1,2,2,1) = usden*d1d2p(i,j)*dx2(i,j)*dx1(im1,j)
          alambd(i,j,1,2,1,2) = usden*d1d2p(i,j)*dx2(i,jm1)*dx1(i,j)
          alambd(i,j,1,2,2,2) = usden*d1d2p(i,j)*dx2(i,jm1)*dx1(im1,j)
          alambd(i,j,2,1,1,1) = usden*d2d1p(i,j)*dx2(i,j)*dx1(i,j)
          alambd(i,j,2,1,2,1) = usden*d2d1p(i,j)*dx2(i,j)*dx1(im1,j)
          alambd(i,j,2,1,1,2) = usden*d2d1p(i,j)*dx2(i,jm1)*dx1(i,j)
          alambd(i,j,2,1,2,2) = usden*d2d1p(i,j)*dx2(i,jm1)*dx1(im1,j)
        enddo
      enddo
!
!--4.4 METRICS FOR ADVECTION AND SCALAR DIFFUSION COEFFICIENTS.
!--------------------------------------------------------------
!
!  dxs1, dxs2: LENGHT OF THE GRID SQUARES SIDES.
!  dxc1, dxc2: LENGHT OF THE GRID SQUARES CENTRES.
!  area: SURFACE OF A GRID SQUARE.
!
      do j=1,jmax
        do i=1,imax
          dxs1(i,j)  = h1pp(i,j)*dx1(i,j)
          dxs2(i,j)  = h2pp(i,j)*dx2(i,j)
          dxc1(i,j)  = h1(i,j)*dx1(i,j)
          dxc2(i,j)  = h2(i,j)*dx2(i,j)
          area(i,j)  = dxc1(i,j)*dxc2(i,j)
        enddo
      enddo
!
      do j=2,jmax
       do i=1,imax
         im1=(i-1)+(imax-2)*(1/i)
         bkappa(i,j,1,1) = dxs2(i,j)/area(i,j)
         bkappa(i,j,1,2) = dxs2(im1,j)/area(i,j)
         bkappa(i,j,2,2) = dxs1(i,j)/area(i,j)
         bkappa(i,j,2,1) = dxs1(i,j-1)/area(i,j)
       enddo
      enddo
!
      do i=1,imax
         j=1
         bkappa(i,j,1,1) = dxs2(i,j)/area(i,j)
         bkappa(i,j,1,2) = dxs2(im1,j)/area(i,j)
         bkappa(i,j,2,2) = dxs1(i,j)/area(i,j)
         bkappa(i,j,2,1) = dxs1(i,1)/area(i,j)
      enddo
!--modification of the metric coefficients for Bering (Arctic)
      if (ltest.eq.3 .and. iberp.ne.0) then
        dxs1(ibera-1,jbera-1) = dxs1(iberp,jberp-1)
        dxs1(ibera,jbera-1)   = dxs1(iberp-1,jberp-1)
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Definition of the domain : reading of the bathymetry               |
!-----------------------------------------------------------------------

!--reading of the bathymetry file (thickness of each level (m) )
!        and of the map tha contains the number of levels)
      open(newunit=bath_om_id,
     &  file='inputdata/clio/bath.om',status='old'
     & ,form='UNFORMATTED')
      read (bath_om_id) zbath
      read (bath_om_id) dzbath
#if ( BATHY == 0 )
      read (bath_om_id) kbath
#endif
      close(bath_om_id)

c fl - reading the bathymetry file in vertical levels
c~ level-bath-0.5_t1_r0_0.4_T.txt
c~ level-bath-133.9_t1-21P_r0_0.4_T.txt
c~ level-bath-133.9_t1-21T_r0_0.4_T.txt
c~ level-bath-133.9_t21P_r0_0.4_T.txt
c~ level-bath-133.9_t21T_r0_0.4_T.txt

#if ( NC_BERG == 2 )
      la_date=palaeo_year-iyear
      write(*,*) 'The date for ice sheet update is ', la_date
#endif 

#if ( BATHY == 1 )
      open(newunit=bath_id,file='inputdata/clio/bath_new.txt',
     &status='old',form='FORMATTED')

#elif ( BATHY >= 2 )
      la_date=palaeo_year-iyear
      write(*,*) 'The date for bathy update is ', la_date

      write(charI,'(I5.5)'), la_date
      name_file ='inputdata/clio/bath_'//trim(charI)//'.txt'

      write(*,*) name_file

      open(newunit=bath_id,file=name_file,
     &status='old',form='FORMATTED')
!      open(unit=21,file='inputdata/clio/bath_new.txt',
!     &status='old',form='FORMATTED')
#endif

#if ( BRINES >= 3 )
      la_date_brines=palaeo_year-iyear
#endif

#if ( BATHY >=1 )
      do I=1,122,1
            read(bath_id,*) kbath(I,:)
      enddo
      close(bath_id)
#endif
             !write(*,*) 'KBATH 1 ', kbath(56-1,3-1),kbath(56,3)


!--initialisation of the depth scale :
      zw(ks2+1) = 0.0
      do k=ks2,ks1,-1
        z(k)  = -zbath(ks1+ks2-k)
        dz(k) = dzbath(ks1+ks2-k)
        zw(k)    = zw(k+1) - dz(k)
        unsdz(k) = 1.0 / dz(k)
        if (zw(k).ne.zero) unszw(k) = 1.0 / zw(k)
      enddo
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      do  k=ks1+1,ks2
        dzw(k) = z(k)  - z(k-1)
        unsdzw(k) = 1.0 / dzw(k)
      enddo

!dmr&nb --- update_bathy

      if(kfond.gt.0) then
!--flat bottom:
      do  j=1,jmax
       do i=1,imax
        kbath(i,j) = min(kfond,(kbath(i,j)*kmax))
       enddo
      enddo
      endif
             !write(*,*) 'KBATH 2 ', kbath(56-1,3-1),kbath(56,3)

!--closing the basin North and South (after save) :
      do i=1,imax
!dmr --- A priori kbath1 non utilise ...
        kbath1(i) = kbath(i,1)

        kbath(i,1) = 0 ! premiere ligne mise a zero au Sud
!dmr --- A priori kbath2 non utilise ...
        kbath2(i) = kbath(i,jmax)

        kbath(i,jmax) = 0 ! derniere ligne mise a zero au Nord
       enddo

!--closing the basin East and Ouest :
      do  j=1,jmax
        kbath(1,j) = 0
        kbath(imax,j) = 0
      enddo
             !write(*,*) 'KBATH 3 ', kbath(56-1,3-1),kbath(56,3)

!--closing between the 2 grids  (except at the Equator) :
      if (ltest.ge.2) then
        do i=1,imax
          if (jsep(i).ne.jeq) kbath(i,jsep(i)-1) = 0
        enddo
      endif
             !write(*,*) 'KBATH 4 ', kbath(56-1,3-1),kbath(56,3)

c nb modif kbath
#if ( BATHY >=1 )
c      do i=1,imax
c        do j=1,jmax
      do i=2,imax-1
        do j=2,jmax-1
c           test si au moins un carre autour
            test_wet(1) = min(kbath(i-1,j+1),kbath(i,j+1),kbath(i-1,j))!NW
            test_wet(2) = min(kbath(i,j+1),kbath(i+1,j+1),kbath(i+1,j))!NE
            test_wet(3) = min(kbath(i+1,j),kbath(i+1,j-1),kbath(i,j-1))!SE
            test_wet(4) = min(kbath(i,j-1),kbath(i-1,j-1),kbath(i-1,j))!SW
            if (kbath(i,j).gt.(test_wet(1))
     &        .and.(kbath(i,j).gt.test_wet(2))
     &        .and.(kbath(i,j).gt.test_wet(3))
     &        .and.(kbath(i,j).gt.test_wet(4))) then
              kbath_bas=maxval(test_wet(:)) !niveau kbath le plus bas permettant echange horizontal
c~               write(*,*), 'pas de carre, reduit kbath',
c~      &                       i,j,kbath(i,j) ,kbath_bas
              kbath(i,j)=kbath_bas

           endif
        enddo
      enddo
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  6 ) Initialisation of the mask and of the indexes of beginning/end of domain.
!-----------------------------------------------------------------------

!--6.1 Centre of the grid  (pt. type scalar) - land/sea mask  pts Scal. :
!---------------------------------------------------------------------------
             !write(*,*) 'KBATH 5 ', kbath(56-1,3-1),kbath(56,3)

!- initialisation of the depth and computation of the index of beginning
!  and end of basin  :
      do j=js1,js2
       is1(j) = imax
       is2(j) = 1
       do  i=ims1,ims2
         km = kbath(i,j)
         if(km.gt.0) then
           kkm = ks1 + ks2 - km
           kniv(i,j,0) = kkm
           hs(i,j) = -zw(kkm)
           if(is1(j).eq.imax) is1(j) = i
           is2(j) = i
           do  k=kkm,kmax
             tms(i,j,k) = 1.0
           enddo
         endif
        enddo
       enddo

      !write(*,*) 'TMS 1 ', tmu(56,3,14)


!--6.2 definition of the "openings" for the E & S boundaries :
!---------------------------------------------------------------

             !write(*,*) 'KBATH 6 ', kbath(56-1,3-1),kbath(56,3)
      if (ltest.ge.1) then
!--half-cycli correspondance for kbath :
      do  j=jcl1,jcl2
        kbath(ims1-1,j) = kbath(ims2,j) ! on copie la colonne ims2 = imax-1 dans la colonne ims1-1 = 1
      enddo
      endif

             !write(*,*) 'IBERP 2', iberp

             !write(*,*) 'KBATH 7 ', kbath(56-1,3-1),kbath(56,3)
      if (ltest.eq.3 .and. iberp.ne.0) then
!--demi correspondance  for kbath, Bering Strait :
        kbath(iberp, jberp) = kbath(iberam,jberam)
        kbath(iberpm,jberp) = kbath(ibera, jberam)
      endif

             !write(*,*) 'KBATH 8 ', kbath(56-1,3-1),kbath(56,3)

!--6.3 Corners of the grid (pt. type velocity) - land/sea mask for velocity :
!--------------------------------------------------------------------------

!- Computation of the depths and computation of the index of beginning end of basin :
      do  j=js1,js2
       iu1(j) = imax
       iu2(j) = 1
       do  i=ims1,ims2
!--for each point (i,j) :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
         km3 = min(kbath(i-1,j-1),kbath(i,j),kbath(i-1,j),kbath(i,j-1))
!         if ((i.eq.56) .and. (j.eq.3)) then
!             !write(*,*) 'TMS 1z ', km3, kbath(i-1,j-1),kbath(i,j),
!     >           kbath(i-1,j), kbath(i,j-1)
!         endif
!--
!--Suppresion  of the island with no surface  :
         nn0vit = 1
         do  n=1,npt0v
           nn0vit = min( nn0vit , abs(i-ipt0v(n)) + abs(j-jpt0v(n)) )
         enddo
!--
          !if ((i.eq.56) .and. (j.eq.3)) then
          !   write(*,*) 'TMS 1a ', nn0vit, km3
          !endif
         if(km3.gt.0 .and. nn0vit.eq.1) then
          !if ((i.eq.56) .and. (j.eq.3)) then
          !   write(*,*) 'TMS 1b ', tmu(56,3,14)
          !endif
           kkm3 = ks1 + ks2 - km3
           kniv(i,j,-1) = kkm3
           hu(i,j) = -zw(kkm3)
           unshu(i,j) = -unszw(kkm3)
           if(iu1(j).eq.imax) iu1(j) = i
           iu2(j) = i
           do  k=kkm3,kmax
             tmu(i,j,k) = 1.0
           !if ((i.eq.56) .and. (j.eq.3) .and. (k.eq.14)) then
           !   write(*,*) 'TMS 1c ', tmu(56,3,14)
           !endif
           enddo
         endif
!--end of the work for each point (i,j).
          enddo
         enddo

      !write(*,*) 'TMS 2 ', tmu(56,3,14)

!--6.4 definition of the  "openings" for  the W & N  boundaries (correspondance) :
!---------------------------------------------------------------------

      if (ltest.ge.1) then
      do j=jcl1,jcl2
!--half cyclic correspondance  kbath :
       kbath(ims2+1,j) = kbath(ims1,j) ! on copie la colonne ims1 = 2 dans la colonne ims2+1 = imax
!--cyclic correspondance for kniv, hu, unshu, tms et tmu :
       kniv(ims1-1,j,0) = kniv(ims2,j,0)
       kniv(ims2+1,j,0) = kniv(ims1,j,0)
       kniv(ims1-1,j,-1) = kniv(ims2,j,-1)
       kniv(ims2+1,j,-1) = kniv(ims1,j,-1)
       hu(ims1-1,j) = hu(ims2,j)
       hu(ims2+1,j) = hu(ims1,j)
       unshu(ims1-1,j) = unshu(ims2,j)
       unshu(ims2+1,j) = unshu(ims1,j)
       do k=ks1,ks2
        tms(ims1-1,j,k) = tms(ims2,j,k)
        tms(ims2+1,j,k) = tms(ims1,j,k)
        tmu(ims1-1,j,k) = tmu(ims2,j,k)
        tmu(ims2+1,j,k) = tmu(ims1,j,k)
       enddo
      enddo
      endif
!--

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Bering Strait, end of the closing :
      if (ltest.eq.3 .and. iberp.eq.0) then
!- if in the middle of 4 "land point"=> no connexion :
        do  j=2,jmax
          if (iberp.eq.0) then
            if (is1(j-1).gt.ims1 .and. is1(j).gt.ims1) then
              iberp = ims1
              jberp = j
            elseif (is2(j-1).lt.ims2 .and. is2(j).lt.ims2 ) then
              iberp = ims2
              jberp = j
            endif
          endif
        enddo
        if (iberp.eq.0) then
          write(clio3_out_id,*) 'STOP in the routine "defgrid" :'
          write(clio3_out_id,*) 'Problem in Bering closing !'
          stop
        endif
        ibera = iberp
        jbera = jberp
        jberpm = jberp - 1
        iberpm = iberp - 1
        jberam = jbera - 1
        iberam = ibera - 1
      endif
      if (ltest.eq.3 .and. iberp.ne.ibera) then
!--half correspondance  for kbath, Bering Strait :
        kbath(ibera, jbera) = kbath(iberpm,jberpm)
        kbath(iberam,jbera) = kbath(iberp, jberpm)
!--correspondance  Bering for kniv, hu, unshu, tms et tmu :
        kniv(iberp, jberp,0) = kniv(iberam,jberam,0)
        kniv(iberpm,jberp,0) = kniv(ibera, jberam,0)
        kniv(ibera, jbera,0) = kniv(iberpm,jberpm,0)
        kniv(iberam,jbera,0) = kniv(iberp, jberpm,0)
        kniv(ibera,jbera,-1) = kniv(iberp,jberp,-1)
        hu(ibera,jbera) = hu(iberp,jberp)
        unshu(ibera,jbera) = unshu(iberp,jberp)
        do  k=ks1,ks2
          tms(iberp, jberp,k) = tms(iberam,jberam,k)
          tms(iberpm,jberp,k) = tms(ibera ,jberam,k)
          tms(ibera, jbera,k) = tms(iberpm,jberpm,k)
          tms(iberam,jbera,k) = tms(iberp ,jberpm,k)
          tmu(ibera,jbera,k) = tmu(iberp,jberp,k)
        enddo
!--Correction of  "iu1(jberp),iu2(jberp)" (necessairy ?) :
!       j = jberp
!       do i=iberp+1,imax
!         if (iu1(j).eq.iberp.and.tmu(i,j,ks2).eq.one) iu1(j) = i
!       enddo
!       do i=iberp-1,1,-1
!         if (iu2(j).eq.iberp.and.tmu(i,j,ks2).eq.one) iu2(j) = i
!       enddo
!--Computation of  the coeff. "alpha" (-> bering = coeff for geostrophic control )
        bering = bering * gpes * hu(iberp,jberp)
     &                  * 0.5 / fs2cor(iberp,jberp)
        bering = max(zero, bering)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

      !write(*,*) "BERING IN DEFGRID == ", bering, iberp

      !write(*,*) 'TMS 3 ', tmu(56,3,14)

!--6.6 Corners of the grid (pt. type velocity) - mask sea+coast : kniv(1)
!--------------------------------------------------------------------------

!- computation of the maximum depth :
      do  j=js1,jmax
       do  i=ims1,imax
!--for each point (i,j) :
         km1 = max(kbath(i-1,j-1),kbath(i,j),kbath(i-1,j),kbath(i,j-1))
         kniv(i,j,1) = kmaxp1 - km1
       enddo
      enddo
!--
      if (ltest.ge.1) then
        do  j=jcl1,jcl2
!--cyclic boundary for kniv(1) :
          kniv(ims1-1,j,1) = kniv(ims2,j,1)
          kniv(ims2+1,j,1) = kniv(ims1,j,1)
        enddo
      endif
      if (ltest.eq.3 .and. iberp.ne.ibera) then
!--correspondance  Bering for kniv(1) :
        kniv(ibera,jbera,1) = kniv(iberp,jberp,1)
        kniv(ibera+1,jbera,1) = kniv(iberp-1,jberp,1)
        kniv(ibera-1,jbera,1) = kniv(iberp+1,jberp,1)
        do  ii=-1,1
          kniv(iberp+ii,jberp+1,1) = kmaxp1
        enddo
      endif




!--6.7 definition of the indexes for the computation of the fluxes :
!-------------------------------------------------------------------

!--initialisation of  isf1 and isf2 :
      do  j=js1,1+js2
       isf1(j) = min(is1(j-1),is1(j))
       isf2(j) = max(is2(j-1),is2(j))
      enddo

!--initialisation of iuf1 and iuf2 , huy and hux :
      do j=1,jmax - 1
        iuf1(j) = min(iu1(j),iu1(j+1))
        iuf2(j) = max(iu2(j),iu2(j+1))
        do  i=1,imax - 1
          huy(i,j) = min(hu(i,j),hu(i,j+1))
          hux(i,j) = min(hu(i,j),hu(i+1,j))
        enddo
       enddo

!--initialisation of imu1 et imu2 , ju1 et ju2 ,
!- correction (empty lign) of is2, iu2, isf2 & iuf2 :
      imu1 = imax
      imu2 = 1
      ju1 = jmax
      ju2 = 1
      do  j=1,jmax
        imu1 = min(imu1, iu1(j))
        imu2 = max(imu2, iu2(j))
        if (iu1(j) .gt.iu2(j) ) then
          iu2(j) = ims2
        else
          if (ju1.eq.jmax) ju1 = j
          ju2 = j
        endif
        if (is1(j) .gt.is2(j) ) is2(j) = ims2
        if (isf1(j).gt.isf2(j)) isf2(j) = ims2
        if (iuf1(j).gt.iuf2(j)) iuf2(j) = ims2
       enddo

!--6.9 initialisation of the indexes indicating the bottom : kfs & kfu
!-----------------------------------------------------------------------

!--the level corresponding to the bottom : kfs (pts Sclal), kfu (pts Vites.)
      do  j=1,jmax
       do i=1,imax
        kfs(i,j) = min( kmax, kniv(i,j,0) )
        kfu(i,j) = min( kmax, kniv(i,j,-1))
        enddo 
       enddo
#if ( SHELFMELT == 1 )
      tmsism(:,:,:) = tms(:,:,:)
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  7 ) Other Indices - Surface for Integration over the whole domain:        |
!-----------------------------------------------------------------------

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--7.1 Initialisation of the arrays indicating the corners :

      do  k=ks1,ks2
       n1coin(k) = 0
       n2coin(k) = 0
       n3coin(k) = 0
       n4coin(k) = 0
       do  j=js1,js2
        ij = (j-1) * imax
        do  i=is1(j),is2(j)
          xx4tms =  tms(i-1,j,k) + tms(i+1,j,k)
     &            + tms(i,j-1,k) + tms(i,j+1,k)
          if (xx4tms.eq.2.0d0 .and. tms(i,j,k).eq.one ) then
            if (tmu(i,j,k).eq.one) then
              n = 1 + n1coin(k)
              if (n.gt.ncomax) goto 1710    ! goto error handling for n > ncomax
              i1coin(n,k) = i + ij
              n1coin(k) = n
            elseif (tmu(i+1,j,k).eq.one) then
              n = 1 + n2coin(k)
              if (n.gt.ncomax) goto 1710    ! goto error handling for n > ncomax
              i2coin(n,k) = i + ij
              n2coin(k) = n
            elseif (tmu(i,j+1,k).eq.one) then
              n = 1 + n3coin(k)
              if (n.gt.ncomax) goto 1710    ! goto error handling for n > ncomax
              i3coin(n,k) = i + ij
              n3coin(k) = n
            elseif (tmu(i+1,j+1,k).eq.one) then
              n = 1 + n4coin(k)
              if (n.gt.ncomax) goto 1710    ! goto error handling for n > ncomax
              i4coin(n,k) = i + ij
              n4coin(k) = n
            endif
          endif
        enddo
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--7.3 Initialisation of the arrays indcating the slopes:

!- direction X :
      nl = 0
      do k=ks1+1,ks2
       do  j=js1,js2
        do  i=is1(j),is2(j)+1
         if (tms(i-1,j,k).ne.one .or. tms(i,j,k).ne.one) cycle
         if (tms(i-1,j,k-1).eq.one .and. tms(i,j,k-1).eq.zero) then
           nl = nl + 1
           if (nl.gt.nlpmax) goto 1730      ! goto error handling for nl > nlpmax
           ijslp(nl) = i - 1 + (j-1) * imax
           kslp(nl) = k
           lslp(nl) = 1
         endif
         if (tms(i-1,j,k-1).eq.zero .and. tms(i,j,k-1).eq.one) then
           nl = nl + 1
           if (nl.gt.nlpmax) goto 1730      ! goto error handling for nl > nlpmax
           ijslp(nl) = i + (j-1) * imax
           kslp(nl) = k
           lslp(nl) = -1
         endif
        enddo
       enddo
      enddo

      nxslp = nl

!- direction Y :
      do  k=ks1+1,ks2
       do  j=js1,js2+1
        do  i=isf1(j),isf2(j)
         if (tms(i,j-1,k).ne.one .or. tms(i,j,k).ne.one) cycle
         if (tms(i,j-1,k-1).eq.one .and. tms(i,j,k-1).eq.zero) then
           nl = nl + 1
           if (nl.gt.nlpmax) goto 1730      ! goto error handling for nl > nlpmax
           ijslp(nl) = i + (j-2) * imax
           kslp(nl) = k
           lslp(nl) = imax
         endif
         if (tms(i,j-1,k-1).eq.zero .and. tms(i,j,k-1).eq.one) then
           nl = nl + 1
           if (nl.gt.nlpmax) goto 1730      ! goto error handling for nl > nlpmax
           ijslp(nl) = i + (j-1) * imax
           kslp(nl) = k
           lslp(nl) = -imax
         endif
        enddo
       enddo
      enddo
      nxyslp = nl

      nxslp0 = nxslp
      nyslp0 = nl - nxslp
      if (kfond.gt.-2) then
        nxslp  = 0
        nxyslp = 0
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--7.5 Initialisation of the Surfaces for the  Integration over the whole domain :

!--initialisation :
      do  llm=0,1
       do  k=1,kmax
        do  j=1,jmax
         do  i=1,imax
           ctmi(i,j,k,llm) = 0.
         enddo
        enddo
       enddo
      enddo

!--Point type Scalar :
      do  k=1,kmax
       do  j=js1,js2
        do  i=is1(j),is2(j)
          ctmi(i,j,k,0) = tms(i,j,k) * cmxy(i,j,0)
        enddo
       enddo
      enddo

!--Point type velocity :
      do  k=1,kmax
       do  j=ju1,ju2
        do  i=iu1(j),iu2(j)
          ctmi(i,j,k,1) = tmu(i,j,k) * cmxy(i,j,3)
        enddo
       enddo
      enddo
!     if (ltest.eq.3 .and. iberp.ne.ibera) then
!--Bering : ajoute 1 pt vitesse : <- si necessaire ?
!       i = iberp
!       j = jberp
!       do 775 k=1,kmax
!         ctmi(i,j,k,1) = tmu(i,j,k) * cmxy(i,j,3)
!775    continue
!     endif

!--computation of the ocean "volume" (in fact : volume /(dx*dy)) in metres :
      volz = 0.0
      do  j=js1,js2
        do  i=is1(j),is2(j)
          volz = volz + cmxy(i,j,0)*hs(i,j)
        enddo
       enddo
      unsvol = 1.0 / volz
!--end of the computation of the  volume.
      ccdxdy = dx * dy / 1.0d+12
      do  j=1,jmax
        do  i=1,imax
!           aire(i,j) = ctmi(i,j,ks2,0)*dx*dy/1.0d+12
            aire(i,j) = ctmi(i,j,ks2,0) * ccdxdy
        enddo
      enddo

! surface for mean w
        zurfow=0.0
        do j=js1,js2
         do i=is1(j),is2(j)
           zurfow=zurfow+aire(i,j)*tms(i,j,ks2)
         enddo
        enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  8 ) Checking - writing of parametres and indices of the domain .       |
!-----------------------------------------------------------------------

      if(nflag.ne.2) goto 899
!--checking
      write(mouchard_id,*) 'dtb, dtu, dts(k) :'
      write(mouchard_id,*)  dtb, dtu, dts
      write(mouchard_id,*) 'dx, dy, pi :'
      write(mouchard_id,*)  dx, dy, pi
      write(mouchard_id,*) 'icheck, jcheck, kcheck :'
      write(mouchard_id,*)  icheck, jcheck, kcheck

      if (kcheck.eq.0 .and. icheck.ge.1 .and. icheck.le.imax) then
!-----
!- Metric Coef. : writing for  i=icheck :
      do  jj1=1,jmax,10
      jj2 = min(jj1+9,jmax)
      write(mouchard_id,'(3(A,I3))') 'cmx(0), i=', icheck,
     &                      ' de j=', jj1, ' a j=', jj2
      write(mouchard_id,*) (cmx(icheck,j,0),j=jj1,jj2)
      enddo
      do  jj1=1,jmax,10
      jj2 = min(jj1+9,jmax)
      write(mouchard_id,'(3(A,I3))') 'cmy(0), i=', icheck,
     &                      ' de j=', jj1, ' a j=', jj2
      write(mouchard_id,*) (cmy(icheck,j,0),j=jj1,jj2)
      enddo
      do  jj1=1,jmax,10
      jj2 = min(jj1+9,jmax)
      write(mouchard_id,'(3(A,I3))') 'cmx(3), i=', icheck,
     &                      ' de j=', jj1, ' a j=', jj2
      write(mouchard_id,*) (cmx(icheck,j,3),j=jj1,jj2)
      enddo
      do  jj1=1,jmax,10
      jj2 = min(jj1+9,jmax)
      write(mouchard_id,'(3(A,I3))') 'cmy(3), i=', icheck,
     &                      ' de j=', jj1, ' a j=', jj2
      write(mouchard_id,*) (cmy(icheck,j,3),j=jj1,jj2)
      enddo
!-----
      endif !! on kcheck, jcheck

      if (kcheck.eq.0 .and. jcheck.ge.1 .and. jcheck.le.jmax) then
!-----
!- Metric coef. : wrting for j=jcheck :
      do  ii1=1,imax,10
      ii2 = min(ii1+9,imax)
      write(mouchard_id,'(3(A,I3))') 'cmx(0), j=', jcheck,
     &                      ' de i=', ii1, ' a i=', ii2
      write(mouchard_id,*) (cmx(i,jcheck,0),i=ii1,ii2)
      enddo
      do  ii1=1,imax,10
      ii2 = min(ii1+9,imax)
      write(mouchard_id,'(3(A,I3))') 'cmy(0), j=', jcheck,
     &                      ' de i=', ii1, ' a i=', ii2
      write(mouchard_id,*) (cmy(i,jcheck,0),i=ii1,ii2)
      enddo

      do  ii1=1,imax,10
      ii2 = min(ii1+9,imax)
      write(mouchard_id,'(3(A,I3))') 'cmx(3), j=', jcheck,
     &                      ' de i=', ii1, ' a i=', ii2
      write(mouchard_id,*) (cmx(i,jcheck,3),i=ii1,ii2)
      enddo
      do  ii1=1,imax,10
      ii2 = min(ii1+9,imax)
      write(mouchard_id,'(3(A,I3))') 'cmy(3), j=', jcheck,
     &                      ' de i=', ii1, ' a i=', ii2
      write(mouchard_id,*) (cmy(i,jcheck,3),i=ii1,ii2)
      enddo
!-----
      endif !! on kcheck, jcheck

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      write(mouchard_id,*) 'z :'
      write(mouchard_id,*)  z
      write(mouchard_id,*) 'dz :'
      write(mouchard_id,*)  dz
      write(mouchard_id,*) 'zw :'
      write(mouchard_id,*)  zw
      write(mouchard_id,*) 'dzw :'
      write(mouchard_id,*)  dzw

      write(mouchard_id,*) 'ahu, ahe :'
      write(mouchard_id,*)  ahu, ahe
      write(mouchard_id,*) 'ahs :'
      write(mouchard_id,*)  ahs
      write(mouchard_id,*) 'ai  :'
      write(mouchard_id,*)  ai
      write(mouchard_id,*) 'slopemax  :'
      write(mouchard_id,*)  slopemax
      write(mouchard_id,*) 'aitd  :'
      write(mouchard_id,*)  aitd
      write(mouchard_id,*) 'slopmgm  :'
      write(mouchard_id,*)  slopmgm
      write(mouchard_id,*) 'afilt,ahh,avv  :',afilt,ahh,avv
      write(mouchard_id,*) 'avnub :'
      write(mouchard_id,*)  avnub
      write(mouchard_id,*) 'avnu0 :'
      write(mouchard_id,*)  avnu0
      write(mouchard_id,*) 'avkb :'
      write(mouchard_id,*)  avkb
      write(mouchard_id,*) 'avk0 :'
      write(mouchard_id,*)  avk0
      write(mouchard_id,*) 'rifumx, rifsmx :'
      write(mouchard_id,*)  rifumx, rifsmx
      write(mouchard_id,*) 'alphah :'
      write(mouchard_id,*)  alphah
      write(mouchard_id,*) 'alphgr :'
      write(mouchard_id,*)  alphgr
      write(mouchard_id,*) 'algrmn :'
      write(mouchard_id,*)  algrmn
      write(mouchard_id,*) 'alphmi :'
      write(mouchard_id,*)  alphmi
      write(mouchard_id,*) 'alphaz :'
      write(mouchard_id,*)  alphaz

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      write(mouchard_id,'(A)') 'Indices linked with the domain :'
      write(mouchard_id,*) 'jsep(i) :'
      write(mouchard_id,'(20I4)') (jsep(i),i=1,imax)
      write(mouchard_id,*) 'jnorth(i) :'
      write(mouchard_id,'(20I4)') (jnorth(i),i=1,imax)

      write(mouchard_id,*) 'dxwi, dywj, xwi1, ywj1 :'
      write(mouchard_id,*)  dxwi, dywj, xwi1, ywj1
      write(mouchard_id,*) 'dxaj, dyai, xaj1, yai1, xwpoln :'
      write(mouchard_id,*)  dxaj, dyai, xaj1, yai1, xwpoln
      write(mouchard_id,*) 'jsepar, jeq, ijsdl, ijudl :'
      write(mouchard_id,*)  jsepar,jeq,ijsdl,ijudl
      write(mouchard_id,*) 'jcl1, jcl2, jdl1, jdl2 :'
      write(mouchard_id,*)  jcl1, jcl2, jdl1, jdl2
!-----
      write(mouchard_id,'(A,2I4)') 
     & 'Detroits : pour fichier "evolu" et "*.x" :'
     &                  //' nvhsf,ndhsf =', nvhsf, ndhsf
      do  nv=1,max(nvhsf,ndhsf)
        write(mouchard_id,'(A,4I4)') tithsf(nv), ishsf(nv), iehsf(nv),
     &                      jshsf(nv), jehsf(nv)
      enddo
      write(mouchard_id,*)
!-----
      write(mouchard_id,*) 'bering :', bering
      write(mouchard_id,*) 'iberp, jberp, ibera, jbera :'
      write(mouchard_id,*)  iberp, jberp, ibera, jbera
      write(mouchard_id,*) 'cmx,cmy(iberp-1,jberp,2) :'
      write(mouchard_id,*)  cmx(iberp-1,jberp,2), cmy(iberp-1,jberp,2)
      write(mouchard_id,*) 'cmx,cmy(iberp,jberp,2) :'
      write(mouchard_id,*)  cmx(iberp,jberp,2), cmy(iberp,jberp,2)
      write(mouchard_id,*) 'cmx,cmy(ibera-1,jbera,2) :'
      write(mouchard_id,*)  cmx(ibera-1,jbera,2), cmy(ibera-1,jbera,2)
      write(mouchard_id,*) 'cmx,cmy(ibera,jbera,2) :'
      write(mouchard_id,*)  cmx(ibera,jbera,2), cmy(ibera,jbera,2)

      do  n=1,npt0v
        write(mouchard_id,*) 'Ajout ile surf=0 en (i,j)=', ipt0v(n), jpt0v(n)
      enddo 
      write(mouchard_id,*) 'n1coin(k) :'
      write(mouchard_id,'(20I4)') (n1coin(k),k=1,kmax)
      write(mouchard_id,*) 'n2coin(k) :'
      write(mouchard_id,'(20I4)') (n2coin(k),k=1,kmax)
      write(mouchard_id,*) 'n3coin(k) :'
      write(mouchard_id,'(20I4)') (n3coin(k),k=1,kmax)
      write(mouchard_id,*) 'n4coin(k) :'
      write(mouchard_id,'(20I4)') (n4coin(k),k=1,kmax)
!     do 892 k=1,kmax
!       write(mouchard_id,*) 'Coins 1, k=', k
!       write(mouchard_id,'(20I5)') (i1coin(n,k),n=1,n1coin(k))
!       write(mouchard_id,*) 'Coins 2, k=', k
!       write(mouchard_id,'(20I5)') (i2coin(n,k),n=1,n2coin(k))
!       write(mouchard_id,*) 'Coins 3, k=', k
!       write(mouchard_id,'(20I5)') (i3coin(n,k),n=1,n3coin(k))
!       write(mouchard_id,*) 'Coins 4, k=', k
!       write(mouchard_id,'(20I5)') (i4coin(n,k),n=1,n4coin(k))
!892  continue
      write(mouchard_id,*) 'nxyslp, nXslope, nYslope, kfond :'
      write(mouchard_id,*) nxyslp, nxslp0, nyslp0, kfond

      write(mouchard_id,*) 'ju1, ju2, imu1, imu2 :'
      write(mouchard_id,*)  ju1, ju2, imu1, imu2
      write(mouchard_id,'(A)') 'Indices Debut(=1)/Fin(=2) :'
      write(mouchard_id,'(A)') '  j | is1,is2   iu1,iu2  iuf1,iuf2'
     &          //' isf1,isf2     (iszon,iezon)(0/1/2/3) :'
      do  j=jmax,1,-1
        write(mouchard_id,'(I3,1X,A1,4(2I4,2X),3(2I4,1X),2I4)') j, '|'
     &   , is1(j),  is2(j),  iu1(j),  iu2(j)
     &   , iuf1(j), iuf2(j), isf1(j), isf2(j)
     &   , (iszon(j,n),iezon(j,n),n=0,nbsmax)
      enddo 

!-- writing of the file  kbath (=sum[tm(k)]) :
!      <-- transfered in the routine  "local"

      write(mouchard_id,*) 'array kfs(i,j) :'
      if (kmax.le.15) then
        nfrc = 125
        fmt1 = 'Z1'
      else
        nfrc = 41
        fmt1 = 'I3'
      endif
      write(fmtrep,'(A,I3,A)') '(',iabs(nfrc),'A1)'
      do  ii1=1,imax,nfrc
       ii2=min(ii1+nfrc-1,imax)
       write(fmt,'(A,I3,2A)') '(',(ii2-ii1+1),fmt1,',A1,I3)'
       write(mouchard_id,*) 'portion de i=',ii1,' a i=',ii2
       write(mouchard_id,fmtrep) (cc1(mod(i,10)),i=ii1,ii2)
       do  j=jmax,1,-1
       write(mouchard_id,fmt) (kfs(i,j),i=ii1,ii2),'|',j
       enddo
      enddo
      write(mouchard_id,*) 'volz, volume-real :'
      write(mouchard_id,*)  volz, volz*dx*dy

 899  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  9 ) Needed in ice dynamic (comprise tms).                           |
!-----------------------------------------------------------------------
!
!--9.1. Ice thickness and advection masks.
!-------------------------------------------
!
      do  j=1,jmax
        do  i=1,imax
          zindfa(i,j) = tms(i,j,ks2)
        enddo
       enddo

!  Diffusion coefficients.
!
       ah = (uvdif/ren)*gridsz
       do  j=1,jmax-1
        jp1 = (j+1)
        do  i=1,imax
          ip1       = (i+1)-(imax-2)*(i/imax)
          dfhu(i,j) = tms(i,j,ks2)*tms(ip1,j,ks2)*ah
          dfhv(i,j) = tms(i,j,ks2)*tms(i,jp1,ks2)*ah
        enddo
       enddo
!     write (89,*) ah
!     write (89,*) dfhv
!

cnb tauc not modified elsewhere
      if (nflag.ne.3) then
!--9.2. Determine cloud optical depths as a FUNCTION of
!       latitude (CHOU ET AL., 1981).
!------------------------------------------------------
!
      do  j=js1,js2
         do  i=is1(j),is2(j)
         alat    = asin(covrai(i,j))/radian
         clat    = (95.0-alat)/10.0
         indx    = 1+int(clat)
         indxp1  = indx+1
         dl      = clat-int(clat)
         dr      = 1.0-dl
         tauc(i,j) = dr*tauco(indx)+dl*tauco(indxp1)
         enddo
       enddo 
!ic0  return
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! 10 ) Solar absorption at the various levels .                        |
!-----------------------------------------------------------------------

      open(newunit=typeaux_dat_id,
     &     file='inputdata/clio/typeaux.dat',status='old')
      do j=1,jmax
!        read(typeaux_dat_id,'(122(I3))') (ntypcn(i,j),i=1,imax)
         read(typeaux_dat_id,*) (ntypcn(i,j),i=1,imax)
      enddo
      close(typeaux_dat_id)

      data zrr /0.58,0.62,0.67,0.77,0.78/
      data zd1 /0.35,0.60,1.00,1.50,1.40/
      data zd2 /23.0,20.0,17.0,14.0,7.90/
      do  j=1,jmax
        do  i=1,imax
          rrro(i,j) = zrr(ntypcn(i,j))
          dd1o(i,j) = zd1(ntypcn(i,j))
          dd2o(i,j) = zd2(ntypcn(i,j))
         enddo
       enddo

      do  j=1,jmax
        do  i=1,imax
          do k=1,kmax
            reslum(i,j,k) = 0.0
          enddo
          znivin = 0.0
          klim=max(11,kfs(i,j))
          do  k=klim,kmax
c~             znivsp =fotr(rrro(i,j),dd1o(i,j),dd2o(i,j),-zw(k+1))
            znivsp = rrro(i,j)*exp(max(-500.d0,zw(k+1)/dd1o(i,j)))
     &              +(1.0-rrro(i,j))*exp(max(-500.d0,zw(k+1)/dd2o(i,j)))

            reslum(i,j,k) = (znivsp - znivin)*tms(i,j,k)
     &                      *dts(k)/dts(ks2)*unsdz(k)/(rho0*cpo)
            znivin = znivsp
          enddo
         enddo
       enddo



!DFG start
!
! set horizontal tracer and momentum grid masks
!
      do j=1,jmax
        hs(1,j)=hs(imax-1,j)
        hs(imax,j)=hs(2,j)
        hu(1,j)=hu(imax-1,j)
        hu(imax,j)=hu(2,j)
        do i=1,imax
          msks(i,j)=1
          if (hs(i,j).eq.0) msks(i,j)=0
          msku(i,j)=1
          if (hu(i,j).eq.0) msku(i,j)=0
          angle(i,j)=angle(i,j)*float(msku(i,j))
        enddo
      enddo
!DFG end

!output of latitude and longitude
      open (newunit=lat_dat_id,file='outputdata/ocean/lat.dat')
      do j=1,jmax
#if ( HRCLIO == 0 )
         write(lat_dat_id,'(122( F10.5))' ) (zlatt(i,j),i=1,imax)
#elif ( HRCLIO == 1 )
         write(lat_dat_id,'(242( F10.5))' ) (zlatt(i,j),i=1,imax)
#endif
      enddo
      write(lat_dat_id,*)
      close (lat_dat_id)
      open (newunit=lon_dat_id,file='outputdata/ocean/long.dat')
      do j=1,jmax
#if ( HRCLIO == 0 )
         write(lon_dat_id,'(122( F10.5))' ) (zlont(i,j),i=1,imax)
#elif ( HRCLIO == 1 )
         write(lon_dat_id,'(242( F10.5))' ) (zlont(i,j),i=1,imax)
#endif
      enddo
      write(lon_dat_id,*)
      close (lon_dat_id)
      open (newunit=mask_defgrid_id,file='outputdata/ocean/mask.dat')
      do j=1,jmax
#if ( HRCLIO == 0 )
         write(mask_defgrid_id,'(122( I1))' ) (msks(i,j),i=1,imax)
#elif ( HRCLIO == 1 )
         write(mask_defgrid_id,'(242( I1))' ) (msks(i,j),i=1,imax)
#endif
      enddo
      write(mask_defgrid_id,*)
      close (mask_defgrid_id)

!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! 11 ) Ice shelves and icebergs .                                         |
!-----------------------------------------------------------------------
! melting of iceberg in the North
!
      iicebern1=98
      iicebern2=115
      jicebern1=46
      jicebern2=55
      areicebn=0.0
      do i=iicebern1,iicebern2
       do j=jicebern1,jicebern2
        areicebn=areicebn+area(i,j)*tms(i,j,ks2)
       enddo
      enddo

      iicebers1=2
      iicebers2=121
      jicebers1=1
      jicebers2=9
      areicebs=0.0
      do i=iicebers1,iicebers2
       do j=jicebers1,jicebers2
        areicebs=areicebs+area(i,j)*tms(i,j,ks2)
       enddo
      enddo

      write (mouchard_id,*) 'surface iceberg melting N,S',areicebn,areicebs

! Definition of the ice-shelves
!     effecta=10E3
      effecta=5E3
      do i=1,imax
       do j=1,jmax
         tmics(i,j)=0.0
       enddo
      enddo
! 1. Rhonne Ice-shelf
      do i=93,103
        tmics(i,2)=tms(i,2,ks2)*dxs1(i,2)*effecta
      enddo
      do i=93,96
        tmics(i,3)=tms(i,3,ks2)*dxs1(i,2)*effecta
      enddo
! 2. East-Weddel 1
      do i=104,106
        tmics(i,3)=tms(i,3,ks2)*dxs1(i,3)*effecta
      enddo
      do i=107,109
        tmics(i,4)=tms(i,4,ks2)*dxs1(i,4)*effecta
      enddo
! 3. East-Weddel 2
      do i=110,114
        tmics(i,4)=tms(i,4,ks2)*dxs1(i,4)*effecta
      enddo
      do i=115,117
        tmics(i,5)=tms(i,5,ks2)*dxs1(i,5)*effecta
      enddo
! 4. Larsen
      do j=4,5
        tmics(93,j)=tms(93,j,ks2)*dxs2(93,j)*effecta
      enddo
! 5. Aimery
      do i=14,17
        tmics(i,5)=tms(i,5,ks2)*dxs1(i,5)*effecta
      enddo
! 6. Ross
      do i=48,61
        tmics(i,2)=tms(i,2,ks2)*dxs1(i,2)*effecta
      enddo
! 7. Getz
      do i=73,77
        tmics(i,3)=tms(i,3,ks2)*dxs1(i,3)*effecta
      enddo
! 8. George VI
      do i=84,87
        tmics(i,4)=tms(i,4,ks2)*dxs1(i,4)*effecta
      enddo

      !write(*,*) 'TMS 4 ', tmu(56,3,14)
!#if ( BATHY == 2 )

!      if (nflag.eq.3) then
!      write(*,*) 'TMS 5 ', tmu(56,3,14),tms_old(56,3,14),
!     >            abs(tms(56,3,14)-tms_old(56,3,14))
!      do i=1,imax
!      do j=1,jmax
!      do k=1,kmax
!         if (abs(tms(i,j,k)-tms_old(i,j,k)).gt.0.0) then
!             write(*,*) "AAARGH ... tms, i,j,k", tms(i,j,k),
!     >                  tms_old(i,j,k), i,j,k
!         endif
!         if (abs(tmu(i,j,k)-tmu_old(i,j,k)).gt.0.0) then
!             write(*,*) "AAARGH ... tmu, i,j,k", tmu(i,j,k),
!     >                  tmu_old(i,j,k), i,j,k
!         endif
!      enddo
!      enddo
!      enddo
!      endif

!#endif

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! 13 ) Pathological cases: error handling                              |
!-----------------------------------------------------------------------

 1710 continue
      write(clio3_out_id,'(A,2I5)') 'Stop in routine "defgrid" :'
     &    //' ncomax too small ! (k,ncomax)= ', k, ncomax
      stop

 1730 continue
      write(clio3_out_id,'(A,4I5)') 'Stop in routine "defgrid" :'
     &    //' nlpmax too small ! (i,j,k,nlpmax)= ', i,j,k, ncomax
      stop



!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- end of the routine defgrid -
      end subroutine defgrid

