!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:40 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:40 CET 2009

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_initecbilt
!-----------------------------------------------------------------------
! *** initialisation of ECBilt
!-----------------------------------------------------------------------


      USE comatm,  only: nlat, nlon
      use comphys, only: iscencel
      use comemic_mod, only: fini
      
c~ #if ( ISOATM >= 1 )
c~       use comsurf_mod, only: fractoc
c~ #endif
      use comunit


c~ #if ( ISOATM >= 1 )
c~       USE isoatm_mod, ONLY : ratio_oceanatm
c~ #endif

      use ECBiltTimer_mod, only: ec_inimdldate

!~ --- [NOTA] temporary addition. I do not see the point of calling celest through there and not through a proper init
      use atmrad_mod, only: celest
      use newunit_mod, only: coef_dat_id, berg_dat_id, gauss_asc_id, win_dat_id,
     &               sum_dat_id, lwrref_dat_id , lwrcoef_dat_id, swrref_dat_id,
     &               swrcoef_dat_id, namelistecbilt_id, GHG_dat_id,
     &               VOLC_dat_id, SUL_dat_id, OZONE_dat_id, mbcs2_cor_id,
     &               scenario2Xco2_dat_id, TSI_RM_dat_id,
     &               book_id, ocbasin_id, ocheattr_id, runoff_id,info_id,
     &               error_id, ipcc_id
#if ( CEMIS == 1 )
     &               , carbon_emission_dat_id
#endif     
      implicit none

! *** open files


      logical :: success

! --- start of former openatinpfiles.h
c-----------------------------------------------------------------------
c *** open statements of atmosphere input files:
c *** units 90-99 are reserved for initial states
c *** units 20 and higher are reserved for other parts of ecbilt
c-----------------------------------------------------------------------

c *** initialisation data

      open(newunit=coef_dat_id,file='inputdata/coef.dat',
     &     status='old',form='unformatted')
      open(newunit=berg_dat_id,file='inputdata/berg.dat',
     &     status='old',form='unformatted')
      open(newunit=gauss_asc_id,file='inputdata/gauss.asc',
     &     status='old',form='formatted')
      open(newunit=win_dat_id,file='inputdata/win.dat',
     &                             status='old',form='formatted')
      open(newunit=sum_dat_id,file='inputdata/sum.dat',
     &                   status='old',form='formatted')

!      open(iuo+14,file='inputdata/ocbas.dat',
!     &                   status='old',form='formatted')

      open(newunit=lwrref_dat_id,file='inputdata/lwrref.dat'
     &      ,form='unformatted')

      open(newunit=lwrcoef_dat_id,file='inputdata/lwrcoef.dat'
     &      ,form='unformatted')

      open(newunit=swrref_dat_id,file='inputdata/swrref.dat'
     &      ,form='unformatted')

      open(newunit=swrcoef_dat_id,file='inputdata/swrcoef.dat'
     &      ,form='unformatted')

      open(newunit=namelistecbilt_id
     &      ,file = 'namelistecbilt',status='old',form='formatted')

      open(newunit=GHG_dat_id,file='inputdata/GHG.dat')
      open(newunit=TSI_RM_dat_id,file='inputdata/TSI-RM.dat')
      open(newunit=VOLC_dat_id,file='inputdata/VOLC.dat')
      open(newunit=SUL_dat_id,file='inputdata/SUL.dat',form='unformatted')
      open(newunit=OZONE_dat_id,file='inputdata/OZONE.dat')

      open(newunit=mbcs2_cor_id,file='inputdata/mbcs2_cor')
      open(newunit=scenario2Xco2_dat_id,file='inputdata/scenario2Xco2.dat')

#if ( CEMIS == 1 )
      open(newunit=carbon_emission_dat_id,
     &    file='inputdata/carbon_emission.dat')
#endif

! --- end of former openatinpfiles.h

      call ec_iniparameterat


! --- start of former openatoutfiles.h

c-----------------------------------------------------------------------
c *** open statements of atmosphere output files:
c *** units 90-99 are reserved for initial states
c *** units 1-20 and above 40  are reserved for other parts of ecbilt
c-----------------------------------------------------------------------


c *** open data output files

      open(newunit=book_id,file='outputdata/globals/book'//fini
     &          ,form='formatted')
      open(newunit=ocbasin_id, file='outputdata/atmos/ocbasin'//fini,
     &           form='formatted')
      open(newunit=ocheattr_id, file='outputdata/atmos/ocheattr'//fini,
     &           form='unformatted')
#if( ROUTEAU == 0 )
      open(newunit=runoff_id, file='outputdata/land/runoff'//fini,
     &           form='formatted')
#endif
c *** open info file

      open(newunit=info_id,file='outputdata/globals/info'//fini
     &          ,form='formatted')

c *** open error file

      open(newunit=error_id,file='outputdata/globals/error'//fini
     &          ,form='formatted')

      open(newunit=ipcc_id,file='outputdata/globals/ipcc'//fini
     &          ,form='formatted')

! --- end of former openatoutpfiles.h


#if ( ISM == 1 )
      if (flgism) call ec_initism
#endif

!******change H. Renssen, 25-02-2003
      if (iscencel.eq.1) success = celest(.true.) !~ dmr [NOTA] .true. == initialization done in that call !!
!      if (iscencel.eq.2) call bretagnon
!******end of change

      call ec_inierror
      call ec_inimdldate
      call ec_iatmpar
      call ec_iatmdyn
      call ec_iatmphys
      call ec_inioutparat
      call ec_inioutparat2

      call ec_atmstate
#if ( CLAQUIN == 1 )
!dmr --- Ajout pour le forcage de Claquin et al., 2003
      print*, "Initialisation Claquin"
      call rad_forc_add(1)
!dmr --- Ajout pour le forcage de Claquin et al., 2003
#endif
#if ( WISOATM == 1 )
      call ec_isoatminit
      ! moved to ec_isoatminit
c~       call fix_input_field(fractoc,ratio_oceanatm,nlat,nlon,neauiso,
c~      &                     ieau17,ieaud)
#endif
      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_iatmpar


      USE comatm, only: pi, radius, rgas, plevel, tlevel, om, p0, rdtime
     &  , fzero, grav, phi, dt, dtime, alogpl2tl2, alogtl12, alogtl1pl2
     &  , rlogtl12, dtt, nlat, nlon, nvl, nsh2, dp, cosfi, sinfi, tanfi
      USE comdyn, only:
      use comphys, only: alphad, alphas, alphav, cc1, cc2, cc3, cpair, cvair
     &  , cwater, dp0, dp1, dp2, dtqmi, dtqmj, dtqmk, gamd, gamma, potfac1
     &  , potfac2, rdtqmi, rdtqmj, rdtqmk, rkappa, rlatfus, rlatsub, rlatvap
     &  , roair, rowat, sboltz, tqmimin, tqmjmin, tqmkmin, tzero
      use comemic_mod, only: iatm, fracto
      use comrunlabel_mod, only: irunlabelf
      use comsurf_mod, only: fractoc, fractn, epss, nld, noc, nse
      use comunit, only: iuo
      use newunit_mod, only: gauss_asc_id
      
      implicit none

      real*8  ari(nlat/2)
      real*8  rj,rw,dumwei

      integer ilat,k,ird,i,j

! *** constants
! *** rlatvap: latent heat of condensation in J/kg
! *** rlatsub: latent heat of sublimation in J/kg
! *** rlatfus: latent heat of fusion in J/kg

      pi=4d0*datan(1d0)
      fzero=sin(0.25*pi)
!     let op om is 2 keer de hoeksnelheid van de aarde !!!!
      om=4d0*pi/(24d0*3600d0)
      grav=9.8
      rgas=287.
      radius=6.37e+6
      rowat=1000.
      roair=1.25
      rlatvap=2.5e+06
      rlatsub=2.8e+06
      rlatfus=0.3e+06
      sboltz=5.67e-08
      cwater=4180.
      cpair=1004.
      cvair=717.
      gamma=cpair/cvair
      rkappa=(gamma-1.)/gamma

      plevel(1)=2.0d+4
      plevel(2)=5.0d+4
      plevel(3)=8.0d+4

      tlevel(1)=3.5d+4
      tlevel(2)=6.5d+4

      dp=plevel(2)-plevel(1)

      dp0=2.0d+4
      dp1=3.0d+4
      dp2=5.0d+4

      p0=1d5

      rlogtl12=1d0/log(tlevel(1)/tlevel(2))
      alogtl12=log(tlevel(1)/tlevel(2))
      alogtl1pl2=log(tlevel(1)/plevel(2))
      alogpl2tl2=log(plevel(2)/tlevel(2))

      potfac1=(tlevel(1)/p0)**rkappa
      potfac2=(tlevel(2)/p0)**rkappa

      gamd=grav/cpair
      tzero=273.15d0
      alphad=roair*cpair
      alphas=roair*rlatsub
      alphav=roair*rlatvap

! *** constants in clausius clapeyron relation

      cc1=0.662*611.2
      cc2=17.67
      cc3=29.66

! *** definitions for tabel of qmax values used in iatmphys
! *** i corresponds to temperature at 650 hPa
! *** j corresponds to temperature difference between ground and 650 hPa
! *** k corresponds to temperature difference between 650 and 350 hPa

      tqmimin=200d0
      dtqmi=2d0
      rdtqmi=1d0/dtqmi
      tqmjmin=-10d0
      dtqmj=2d0
      rdtqmj=1d0/dtqmj
      tqmkmin=5d0
      dtqmk=2d0
      rdtqmk=1d0/dtqmk

! *** time step of the atmospheric model:
! *** dt    : fraction of one day
! *** dtime : in seconds
! *** dtt   : dimensionless

      dt     = 1d0/dble(iatm)
      dtime  = dt*(24d0*3600d0)
      rdtime = 1d0/dtime
      dtt    = dt*pi*4d0

! *** gauss points and weights

      ilat=nlat/2
  10  continue
        read(gauss_asc_id,220,end=15) rj,rw
        ird=int(rj)
        if (ird.eq.ilat) then
          do i=1,ird
            read(gauss_asc_id,220) ari(i),dumwei
          enddo
          goto 20
        else
          goto 10
        endif
  15    continue
        call ec_error(4)
  20  continue

      close(gauss_asc_id)

      do i=1,ilat
        phi(i)=-ari(ilat+1-i)
        phi(ilat+i)=ari(i)
      enddo

      do i=1,nlat
        phi(i)=asin(phi(i))
      enddo
c~ #if ( PMIP == 1 )
c~ !mab:
c~       open(1010,file='phi_lat.dat',form='unformatted',status='replace')
c~       write(1010) (phi(i),i=1,nlat)
c~       close(1010)
c~ #endif

      do i=1,nlat
        cosfi(i)=cos(phi(i))
        sinfi(i)=sin(phi(i))
        tanfi(i)=tan(phi(i))
      enddo

! *** land/sea/sea-ice fraction


      do j=1,nlon
        do i=1,nlat
          fractoc(i,j)=fracto(i,j)
          fractn(i,j,noc)=fractoc(i,j)
          fractn(i,j,nld)=1-fractoc(i,j)
          fractn(i,j,nse)=0.0
        enddo
      enddo

!      do j=nlat,1,-1
!         write (6, '(i3,x,64(i1))') j, (int(fracto(j,i)), i = 1, nlon)
!      enddo

![moved to comsurf_mod as parameter]      epss=1d-10


220   format(f18.10,f17.10)

      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_iniparameterat
!-----------------------------------------------------------------------
! *** initialisation of parameters from the namelist
!-----------------------------------------------------------------------


      USE comatm
      USE comdyn
      USE comphys
      use comemic_mod, only: irunlabel, isatfor
      use comrunlabel_mod, only: irunlabelf
      use comsurf_mod
      use comunit, only: iuo
      use comatfor, only: nafyear, nbsatfor


      use global_constants_mod, only: solar_constant
      use newunit_mod, only: namelistecbilt_id, parameterschk_id
      
      implicit none


      NAMELIST /runatctl/ iadyn,iaphys,ipert,initfield,initdate
      NAMELIST /dispar/   tdis,addisl,addish,trel,tdif,idif,initfield,initdate
      NAMELIST /dfmpar/  rrdef1,rrdef2,h0
      NAMELIST /moipar/  ihavm,ivavm,imsink,tdifq,gpm500,relhmax,
     *                   hmoisr,umoisr,rainmax
      NAMELIST /cloudpar/ relhcrit,relhfac
      NAMELIST /forpar/  iartif,ipvf1,ipvf2,ipvf3,ipvf4,ipvf5
      NAMELIST /radpar/    solarc,iradcloud,iscenghg,isghgstrt,
     *                     iscentsi,istsistrt,iscenvol,isvolstrt,
     *                     iscensul,issulstrt,isceno3,iso3strt ,
     *                     iscencel,iens,numens,emisoc,emisse,emisld,
     *                     albin,albis,albice,alphd,alphdi,alphs,cgren,
     *                     albocef,facttsi,bup,AMPWIR,EXPIR,HPROFTROP,
     *                     HPROFEQ,HPROFW,AMPEQIR,HPROFAN,AMPANIR,
     *                     eccf,oblf,omwebf,AMPANIR2,HPROFAN2,
     *                     irunlabeloffset,iscenyear, cdrag, evfac
     *                   , mag_alpha
      NAMELIST /satfor/   isatfor,nbsatfor,nafyear,iclimflux
      NAMELIST /fluxpar/ cdrag,cwdrag,dragan,dragla,uv10rfx,uv10m,
     *                   uv10rws,ndayws
!dmr @-@ iceb0
     *    , uv10rwv
!dmr @-@ iceb0
      NAMELIST /fluxcorw/ corAN,corPN,corAC,corID,corAS,corPS,
     *                   corAA


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! intgrtn parameter:                                                   C
! iadyn:      with (1) or without (0) atmospheric dynamics.            C
! iaphys:     with (1) or without (0) atmospheric physics.             C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! disspt parameter:                                                    C
! tdis:       Ekman dissipation time scale atlower level               C
! addisl:     parameter for computing dissipation time scale for land. C
! addish:     parameter for computing dissipation time scale over      C
!             mountains                                                C
! trel:       temperature relaxation time scale.                       C
! tdif:       damping time scale of hyperviscosity of the smallest     C
!             waves at all levels                                      C
! idif:       determines scale-selectivity of hyperviscosity: power of C
!             Laplace operator on PV                                   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! dfmtion parameter:                                                   C
! rrdef1:     rossby radius of deformation layer 1 (200 - 500 hPa).    C
! rrdef2:     rossby radius of deformation layer 2 (500 - 800 hPa).    C
! h0:         mountain scale height.                                   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! moisture Parameter:                                                  C
! ihavm:      with (1) or without (0) horizontal divergence of moistureC
! ivavm:      with (1) or without (0) vertical divergence of moisture. C
! imsink:     with (1) or without (0) the source and sink of moisture. C
! tdifq:      diffusion time scale for moisture in days                C
! gpm500  :   mean 500 hPa geopotential height [m] used in the         C
!             calculation of the groundpressure and temperature        C
! relhmax:    maximum relative humidity before saturation.             C
! hmoisr:     reduction factor of mountain heights in order to tune    C
!             the amount of water that is allowed to pass a            C
!             topographic barier                                       C
! umoisr:     reduction factor of 800 hPa winds used in advecting the  C
!             moisture                                                 C
! rainmax:    maximum rate of dynamic and convective rain in m/s.      C
!             if the maximum is reached by both type of rain           C
!             oversaturation occurs.                                   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Surface flux Parameters:                                             C
! cdrag:      coefficient in sensible and latent air-surface heat flux C
! cwdrag:     coefficient in winds stress                              C
! uv10rfx:    reduction factor of 800 hPa winds to obtain              C
!             10 mtr wind used in the definition of the surface        C
!             heat fluxes                                              C
!             used in SUBROUTINE wind10 of atmdyn.f                    C
! uv10m:      minimum value of 10 mtr wind used in all surface fluxes  C
!             used in SUBROUTINE wind10 of atmdyn.f                    C
! uv10rws:    reduction factor of 800 hPa winds to obtain              C
!             10 mtr wind used in the definition of the wind stress    C
!             to drive the ocean currents                              C
!             used in SUBROUTINE wind10 of atmdyn.f                    C
! uv10rwv:    JONO: reduction factor 800hPa to 10m to obtain windvectorC
!             to drive the icebergs                                    C
! dsnm:       maximum land snow depth                                  C
! bmoism:     maximum bottom moisture.                                 C
! ndayws:     averaging period of wind stresses in days                C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Flux correction parameters                                           C
! corAN       flux correction in the North Atlantic                    C
! corPN       flux correction in the North Pacific                     C
! corAC       flux correction in the Arctic                            C
! corID       flux correction in the Indian                            C
! corAS       flux correction in the South Atlantic                    C
! corPS       flux correction in the South Pacific                     C
! corAA       flux correction in the Southern Ocean                    C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! forcing parameter:                                                   C
! iartif:     with (1) or without (0) artificial forcing               C
! ipvf1 :     with (1) or without (0) diabatic heating                 C
! ipvf2 :     with (1) or without (0) advection of f by divergent wind C
! ipvf3 :     with (1) or without (0) stretching term                  C
! ipvf4 :     with (1) or without (0) advection of zeta by divergent   C
!             wind                                                     C
! ipvf5 :     with (1) or without (0) advection of temperature by      C
!                                  divergent wind                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! tcont parameter:                                                     C
! solarc:     solar constant.                                          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! satfor parameter:                                                    C
! isatfor:    if (1) saves nafyear of atmospheric data to disk to be   C
!             used to drive the ocean in uncoupled mode                C
! nbsatfor:   first year in the integrations to start saving nafyears  C
!             of data                                                  C
! nafyear:    number of years of data to save                          C
! iclimflux:  if (1) daily climatological sst and heatflux are output  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      iadyn    = 1
      iaphys   = 1
      ipert    = 0
      initfield= 1
      initdate = 0

      tdis   = 3.0
      addisl = 0.5
      addish = 0.5
      trel   = 50.0
      tdif   = 1.0
      idif   = 4

      rrdef1 = .110
      rrdef2 = .070
      h0     = 3.

      ihavm  = 1
      ivavm  = 1
      imsink = 1
      tdifq  = 5d0
      gpm500 = 5500
      relhmax= 1.0
      hmoisr = 1d0
      umoisr = 0.8
      rainmax= 1e-5

      cdrag  = 1.4e-03
      cwdrag = 2.0e-03
      dragan = 10.0
      dragla = 10.0
      uv10rfx= 0.8
      uv10m  = 4.
      uv10rws= 0.8
      ndayws = 30
!dmr @-@ iceb0
      uv10rwv= 0.8
!dmr @-@ iceb0
      corAN  = -0.075d0
      corPN  =  1.0d0
      corAC  = -0.20d0
      corID  =  0.0d0
      corAS  = -0.075d0
      corPS  =  0.d0
      corAA  =  0.d0

      relhcrit = 0.5d0
      relhfac  = 1.0d0


      iartif = 0
      ipvf1  = 1
      ipvf2  = 1
      ipvf3  = 1
      ipvf4  = 1
      ipvf5  = 1

!     solarc = 1353

! dmr      solarc = 1365
      solarc = solar_constant
      solarm=solarc

! dmr [NOTA]: somehow, solarc variable seems useless, it is recomputed in atmphys0.f from solarm
! dmr [DELETE] => solarc as a global variable, move to local variable in ec_solar

! dmr moved the solar constant to the global_constants module


      iradcloud= 0
      iscenghg = 0
      irunlabeloffset = 3000
      iscenyear = 0
      isghgstrt= 1990
      iscentsi = 0
      istsistrt=1500
      facttsi =1.0
      iscenvol = 0
      isvolstrt=1500
      iscensul = 0
      issulstrt=1850
      isceno3 = 0
      iso3strt=1860
      iscencel = 0
      eccf   = 0.016724
      oblf   = 23.446
      omwebf = 102.04
      iens=0
      bup=0.13
      numens=1
      emisoc=1.0
      emisse=1.0
      emisld=1.0
      albin = 0.43
      albis = 0.43
      albice = 0.43
      alphd  = 0.70
      alphdi = 0.62
      alphs  = 0.53
      cgren = 0.04
      albocef =1.0
      AMPWIR    = 1.0
      AMPEQIR   = 1.0
      EXPIR     = 0.3333
      HPROFTROP = 1.0
      HPROFEQ   = 1.0
      HPROFW    = 1.0
      HPROFAN   = 1.0
      AMPANIR   = 1.0
      HPROFAN2  = 1.0
      AMPANIR2  = 1.0

      mag_alpha = 1.0
      evfac=1d0

      isatfor = 0
      nbsatfor= 0
      nafyear = 0
      iclimflux = 0


      read(namelistecbilt_id, NML = runatctl)
      read(namelistecbilt_id, NML = dispar)
      read(namelistecbilt_id, NML = dfmpar)
      read(namelistecbilt_id, NML = moipar)
      read(namelistecbilt_id, NML = cloudpar)
      read(namelistecbilt_id, NML = forpar)
      read(namelistecbilt_id, NML = radpar)
      read(namelistecbilt_id, NML = satfor)
      read(namelistecbilt_id, NML = fluxpar)
      read(namelistecbilt_id, NML = fluxcorw)

#if ( F_PALAEO >= 1 )
! dmr&aurel : those two are read in radpar above
! --- in case of F_PALAEO we ignore possible use conflicts
      iscenghg = 1
      iscencel = 1
#endif


      open(newunit=parameterschk_id,file="parameters.chk"
     >    ,form="formatted")

      write(parameterschk_id, 900) 'iaphys   =', iaphys
      write(parameterschk_id, 900) 'iadyn    =', iadyn
      write(parameterschk_id, 900) 'ipert    =', ipert
      write(parameterschk_id, 900) 'initfield=', initfield
      write(parameterschk_id, 900) 'initdate =', initdate

      irunlabelf=irunlabel-initdate
      write(parameterschk_id, 900) 'irunlabelf=', irunlabelf

      write(parameterschk_id, 910) 'tdis     =', tdis
      write(parameterschk_id, 910) 'addisl   =', addisl
      write(parameterschk_id, 910) 'addish   =', addish
      write(parameterschk_id, 910) 'trel     =', trel
      write(parameterschk_id, 910) 'tdif     =', tdif
      write(parameterschk_id, 900) 'idif     =', idif

      write(parameterschk_id, 910) 'h0       =', h0
      write(parameterschk_id, 910) 'rrdef1   =', rrdef1
      write(parameterschk_id, 910) 'rrdef2   =', rrdef2

      write(parameterschk_id, 910) 'relhcrit =', relhcrit
      write(parameterschk_id, 910) 'relhfac  =', relhfac

      write(parameterschk_id, 910) 'cdrag    =', cdrag
      write(parameterschk_id, 910) 'cwdrag   =', cwdrag
      write(parameterschk_id, 910) 'dragan   =', dragan
      write(parameterschk_id, 910) 'dragla   =', dragla
      write(parameterschk_id, 910) 'uv10rfx  =', uv10rfx
      write(parameterschk_id, 910) 'uv10m    =', uv10m
      write(parameterschk_id, 910) 'uv10rws  =', uv10rws
      write(parameterschk_id, 900) 'ndayws   =', ndayws
!dmr @-@ iceb0
      write(parameterschk_id, 910) 'uv10rwv  =', uv10rwv
!dmr @-@ iceb0
      write(parameterschk_id, 910) 'corAN    =', corAN
      write(parameterschk_id, 910) 'corPN    =', corPN
      write(parameterschk_id, 910) 'corAC    =', corAC
      write(parameterschk_id, 910) 'corID    =', corID
      write(parameterschk_id, 910) 'corAS    =', corAS
      write(parameterschk_id, 910) 'corPS    =', corPS
      write(parameterschk_id, 910) 'corAA    =', corAA


      write(parameterschk_id, 900) 'ihavm    =', ihavm
      write(parameterschk_id, 900) 'ivavm    =', ivavm
      write(parameterschk_id, 900) 'imsink   =', imsink
      write(parameterschk_id, 910) 'tdifq    =', tdifq
      write(parameterschk_id, 910) 'gpm500   =', gpm500
      write(parameterschk_id, 910) 'relhmax  =', relhmax
      write(parameterschk_id, 910) 'hmoisr   =', hmoisr
      write(parameterschk_id, 910) 'umoisr   =', umoisr
      write(parameterschk_id, 910) 'rainmax  =', rainmax


      write(parameterschk_id, 900) 'iartif   =', iartif
      write(parameterschk_id, 900) 'ipvf1    =', ipvf1
      write(parameterschk_id, 900) 'ipvf2    =', ipvf2
      write(parameterschk_id, 900) 'ipvf3    =', ipvf3
      write(parameterschk_id, 900) 'ipvf4    =', ipvf4
      write(parameterschk_id, 900) 'ipvf5    =', ipvf5

      write(parameterschk_id, 910) 'solarc   =', solarc
      write(parameterschk_id, 900) 'iradcloud=', iradcloud
      write(parameterschk_id, 900) 'iscenghg =', iscenghg
      write(parameterschk_id, 900) 'isghgstrt=', isghgstrt
      write(parameterschk_id, 900) 'iscentsi =', iscentsi
      write(parameterschk_id, 900) 'istsistrt=', istsistrt
      write(parameterschk_id, 910) 'facttsi  =', facttsi
      write(parameterschk_id, 900) 'iscenvol =', iscenvol
      write(parameterschk_id, 900) 'isvolstrt=', isvolstrt
      write(parameterschk_id, 900) 'iscensul =', iscensul
      write(parameterschk_id, 900) 'issulstrt=', issulstrt
      write(parameterschk_id, 910) 'bup      =', bup

      write(parameterschk_id, 900) 'isatfor  =', isatfor
      write(parameterschk_id, 900) 'nbsatfor =', nbsatfor
      write(parameterschk_id, 900) 'nafyear  =', nafyear
      write(parameterschk_id, 900) 'iclimflux=', iclimflux

      call flush(parameterschk_id)
      emisn(noc)=emisoc
      emisn(nse)=emisse
      emisn(nld)=emisld


900   format(a12,1x,i6)
910   format(a12,1x,e12.5)

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_inioutparat
!-----------------------------------------------------------------------


      USE comatm,  only: nlat, nlon
      use comdiag
      use comemic_mod, only:
      use comunit, only: iuo


      use newunit_mod, only: namelistecbilt_id
      
      implicit none


      integer  ts(20),t(20),tstrat(20),cp(20),lsp(20),pp(20),evap(20),
     *    eminp(20),q(20),r(20),tcc(20),u(20),v(20),uv10(20),ageu(20),
     *         agev(20),omega(20),psi(20),qgpv(20),gh(20),chi(20),
     *         sp(20),shf(20),lhf(20), hforc(20),vforc(20),ssr(20),
     *         tsr(20),ttr(20),str(20),albs(20),albp(20),ustress(20),
     *         vstress(20),cdragw(20), cdragv(20),bm(20),sdl(20),
     *         runoffo(20), runoffl(20),dumt1(20),dumt2(20),dumu1(20),
     *         dumu2(20),snow(20),t2m(20),hic(20),swrs(20), lwrs(20) !mohr
#if ( ISOATM >= 1 )
     *        ,pp18(20)
     *        ,sn18(20)
     *        ,pp17(20)
     *        ,sn17(20)
     *        ,ppd(20)
     *        ,snd(20)
     *        ,q18(20)
     *        ,q17(20)
     *        ,qd(20)
#endif
#if ( ISOATM >= 2 )
     *        ,bm18(20)
     *        ,sdl18(20)
     *        ,bm17(20)
     *        ,sdl17(20)
     *        ,bmd(20)
     *        ,sdld(20)
#endif
      integer i

!!mohr
!! swrs --> shortwave radiation at surface
!! lwrs --> longwave radiation at surface

      NAMELIST /outatctl/ixout,ioutdaily,ioutyearly,meantype,
     &                    meanyl,meantot,ifrendat
      NAMELIST /wratpar/ ts,t,tstrat,t2m,cp,lsp,pp,snow,evap,eminp,q,
     * r,tcc,u,v,uv10,ageu,agev,omega,psi,qgpv,gh,chi,sp,shf,lhf,
     * hforc,vforc,ssr,tsr,ttr,str,swrs,lwrs,albs,albp,ustress,vstress,cdragw,
     * cdragv,bm,sdl,runoffo, runoffl,hic,dumt1,dumt2,dumu1,dumu2
#if ( ISOATM >= 1 )
     *,pp18 ,sn18 ,pp17 ,sn17 ,ppd ,snd, q18, q17, qd
#endif
#if ( ISOATM >= 2 )
     *,bm18 ,sdl18 ,bm17 ,sdl17 ,bmd ,sdld
#endif
      NAMELIST /flgout/irad
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ixout:        output frequency for instantanous field in days        C
! meantype      output monthly mean fields.                            C
! meantype      output seasonal mean fields.                           C
! meantot       computes whole period monthly or seasonal mean fields. C
! meanyl        computes yearly monthly or seansonal mean fields.      C
! ioutdaily     output instantanous fields.                            C
! ioutyearly    output yearly mean fields.                             C
! ts:           output surface temperature.                            C
! t:            output temperature.                                    C
! u             output wind component U.                               C
! v:            output wind component V.                               C
! om:           output wind component omega.                           C
! psi:          output stream FUNCTION.                                C
! shf:          output surface sensible heat flux.                     C
! lhf:          output surface latent heat flux.                       C
! lsp:          output large scale precipitation.                      C
! cp:           output convective precipitation.                       C
! q:            output specific humidity.                              C
! r:            output relative humidity.                              C
! ageu:         output ageostrophic wind component U.                  C
! agev:         output ageostrophic wind component V.                  C
! ssr:          output surface solar radiation (downward).             C
! tsr:          output top solar radiation (downward).                 C
! ttr:          output top thermal radiation (upward).                 C
! str:          output surface thermal radiation (upward).             C
! pp:           output total percipitation.                            C
! evap:         output evaporation.                                    C
! eminp:        output evaporation minus precipitation.                C
! albs:         output surface albedo.                                 C
! albp:         output planetary albedo.                               C
! ustress:      output u wind stress.                                  C
! vstress:      output v wind stress.                                  C
! runoffo:      output runoff over ocean.                              C
! runoffl:      output runoff over land.                               C
! sdl:          output snow depth over land.                           C
! chi:          output velocity potential.                             C
! sp:           output surface pressure                                C
! uv10:         output wind magnitude at 10 meter height.              C
! cdragw:       output drag coefficient of wind.                       C
! cdragv:       output drag coefficient of sensible heat flux.         C
! richar:       output richardson number.                              C
! qgpv:         output quasi geostrophic potential vorticity.          C
! gh:           output geopotential height.                            C
! tcc:          output total cloud cover.                              c
! dumt1:        free output variable at 350 hPa, 650 hPa and 1000 hPa  c
! dumt2:        free output variable at 350 hPa, 650 hPa and 1000 hPa  c
! dumu1:        free output variable at 200 hPa, 500 hPa and 800 hPa   c
! dumu2:        free output variable at 200 hPa, 500 hPa and 800 hPa   c
! mohr
! swrs:         output shortwave radiation at surface                  c
! lwrs:         output longwave radiation at surface                   c
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!** Here following are default values for above mentioned parameters.
!** These parameters can be updated in the namelist.

      itel=0
      instcount=0

      ixout  = 30
      ioutdaily = 0
      ioutyearly = 0
      meantype = 1
      meantot  = 0
      meanyl   = 1
      irad     = 1

      do i = 1, 20
        ts(i)      = 0
        t(i)       = 0
        tstrat(i)  = 0
        t2m(i)     = 0
        cp(i)      = 0
        lsp(i)     = 0
        pp(i)      = 0
        snow(i)    = 0
        evap(i)    = 0
        eminp(i)   = 0
        q(i)       = 0
        r(i)       = 0
        tcc(i)     = 0
        u(i)       = 0
        v(i)       = 0
        uv10(i)    = 0
        ageu(i)    = 0
        agev(i)    = 0
        omega(i)   = 0
        psi(i)     = 0
        qgpv(i)    = 0
        gh(i)      = 0
        chi(i)     = 0
        sp(i)      = 0
        shf(i)     = 0
        lhf(i)     = 0
        hforc(i)   = 0
        vforc(i)   = 0
        ssr(i)     = 0
        tsr(i)     = 0
        ttr(i)     = 0
        str(i)     = 0
        albs(i)    = 0
        albp(i)    = 0
        ustress(i) = 0
        vstress(i) = 0
        cdragw(i)  = 0
        cdragv(i)  = 0
        bm(i)      = 0
        sdl(i)     = 0
        runoffo(i) = 0
        runoffl(i) = 0
        hic(i)     = 0
        dumt1(i)   = 0
        dumt2(i)   = 0
        dumu1(i)   = 0
        dumu2(i)   = 0
        swrs(i)    = 0 !mohr
        lwrs(i)    = 0 !mohr
#if ( ISOATM >= 1 )
        pp18(i)    = 0
        sn18(i)    = 0
        pp17(i)    = 0
        sn17(i)    = 0
        ppd(i)     = 0
        snd(i)     = 0
        q18(i)     = 0
        q17(i)     = 0
        qd(i)      = 0
#endif
#if ( ISOATM >= 2 )
        bm18(i)= 0
        sdl18(i)   = 0
        bm17(i)= 0
        sdl17(i)   = 0
        bmd(i) = 0
        sdld(i)    = 0
#endif

      enddo

      read(namelistecbilt_id, NML = outatctl)
      read(namelistecbilt_id, NML = wratpar)
      read(namelistecbilt_id, NML = flgout)

      do i = 1, 20
        newts(i)     = ts(i)
        newt(i)      = t(i)
        newtstrat(i) = tstrat(i)
        newt2m(i)    = t2m(i)
        newcorain(i) = cp(i)
        newdyrain(i) = lsp(i)
        newtorain(i) = pp(i)
        newsnow(i)   = snow(i)
        newevap(i)   = evap(i)
        neweminp(i)  = eminp(i)
        newrmoisg(i) = q(i)
        newrelhum(i) = r(i)
        newtcc(i)    = tcc(i)
        newu(i)      = u(i)
        newv(i)      = v(i)
        newuv10(i)   = uv10(i)
        newdivu(i)   = ageu(i)
        newdivv(i)   = agev(i)
        newomega(i)  = omega(i)
        newpsi(i)    = psi(i)
        newqgpv(i)   = qgpv(i)
        newgh(i)     = gh(i)
        newchi(i)    = chi(i)
        newsp(i)     = sp(i)
        newhflux(i)  = shf(i)
        neweflux(i)  = lhf(i)
        newhforg(i)  = hforc(i)
        newvforg(i)  = vforc(i)
        newssr(i)    = ssr(i)
        newtsr(i)    = tsr(i)
        newttr(i)    = ttr(i)
        newstr(i)    = str(i)
        newalbs(i)   = albs(i)
        newalbp(i)   = albp(i)
        newustress(i)= ustress(i)
        newvstress(i)= vstress(i)
        newcdragw(i) = cdragw(i)
        newcdragv(i) = cdragv(i)
        newbmoisg(i) = bm(i)
        newsdl(i)    = sdl(i)
        newrunoffo(i)= runoffo(i)
        newrunoffl(i)= runoffl(i)
        newhic(i)    = hic(i)
        newdumt1(i)  = dumt1(i)
        newdumt2(i)  = dumt2(i)
        newdumu1(i)  = dumu1(i)
        newdumu2(i)  = dumu2(i)
        newswrs(i)   = swrs(i)     !mohr
        newlwrs(i)   = lwrs(i)     !mohr

#if ( ISOATM >= 1 )
        newpp18(i)   = pp18(i)
        newsn18(i)   = sn18(i)
        newpp17(i)   = pp17(i)
        newsn17(i)   = sn17(i)
        newppd(i)    = ppd(i)
        newsnd(i)    = snd(i)
        newrmoisg18(i) = q18(i)
        newrmoisg17(i) = q17(i)
        newrmoisgd(i)  = qd(i)
#endif
#if ( ISOATM >= 2 )
        newbmoisg18(i) = bm18(i)
        newsdl18(i)  = sdl18(i)
        newbmoisg17(i) = bm17(i)
        newsdl17(i)  = sdl17(i)
        newbmoisgd(i)= bmd(i)
        newsdld(i)   = sdld(i)
#endif
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_inioutparat2
!-----------------------------------------------------------------------


       use comatm
       use comdiag
       use comemic_mod, only:
       use comunit


      implicit none

      character*60 part1
      integer      totvar(80,20)
      integer      i,l
      integer :: outp_atmosparam_id


      NAMELIST /outatctl/ixout,ifrendat
      NAMELIST /flgout/irad
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ixout:        output frequency for instantanous field in days        C
! meantype      output monthly mean fields.                            C
! meantype      output seasonal mean fields.                           C
! meantot       computes whole period monthly or seasonal mean fields. C
! meanyl        computes yearly monthly or seansonal mean fields.      C
! ioutdaily     output instantanous fields.                            C
! ioutyearly    output yearly mean fields.                             C
! ts:           output surface temperature.                            C
! t:            output temperature.                                    C
! u             output wind component U.                               C
! v:            output wind component V.                               C
! om:           output wind component omega.                           C
! psi:          output stream function.                                C
! shf:          output surface sensible heat flux.                     C
! lhf:          output surface latent heat flux.                       C
! lsp:          output large scale precipitation.                      C
! cp:           output convective precipitation.                       C
! q:            output specific humidity.                              C
! r:            output relative humidity.                              C
! ageu:         output ageostrophic wind component U.                  C
! agev:         output ageostrophic wind component V.                  C
! ssr:          output surface solar radiation (downward).             C
! tsr:          output top solar radiation (downward).                 C
! ttr:          output top thermal radiation (upward).                 C
! str:          output surface thermal radiation (upward).             C
! pp:           output total percipitation.                            C
! evap:         output evaporation.                                    C
! eminp:        output evaporation minus precipitation.                C
! albs:         output surface albedo.                                 C
! albp:         output planetary albedo.                               C
! ustress:      output u wind stress.                                  C
! vstress:      output v wind stress.                                  C
! runoffo:      output runoff over ocean.                              C
! runoffl:      output runoff over land.                               C
! sdl:          output snow depth over land.                           C
! chi:          output velocity potential.                             C
! sp:           output surface pressure                                C
! uv10:         output wind magnitude at 10 meter height.              C
! cdragw:       output drag coefficient of wind.                       C
! cdragv:       output drag coefficient of sensible heat flux.         C
! richar:       output richardson number.                              C
! qgpv:         output quasi geostrophic potential vorticity.          C
! gh:           output geopotential height.                            C
! tcc:          output total cloud cover.                              c
! dumt1:        free output variable at 350 hPa, 650 hPa and 1000 hPa  c
! dumt2:        free output variable at 350 hPa, 650 hPa and 1000 hPa  c
! dumu1:        free output variable at 200 hPa, 500 hPa and 800 hPa   c
! dumu2:        free output variable at 200 hPa, 500 hPa and 800 hPa   c
! mohr
! swrs:         shortwave radiation at surface
! lwrs:         longwave  radiation at surface
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!** Here following are default values for above mentioned parameters.
!** These parameters can be updated in the namelist.

      itel=0
      instcount=0

      ixout  = 30
      ioutdaily = 0
      ioutyearly = 0
      meantype = 0
      meantot  = 0
      meanyl   = 0
      irad     = 1
      thirdd(:)= 0

      open(newunit=outp_atmosparam_id,file='outp_atmos.param')
      read(outp_atmosparam_id,'(/,/,/,/,/,/,/,/)')
      read(outp_atmosparam_id,*) numtotvar,fill_value,missing_value
      read(outp_atmosparam_id,*)

      do i=1,numtotvar

         read(outp_atmosparam_id,"(A)") part1
         nametotvar(i,1)=trim(part1)
         read(outp_atmosparam_id,*) (nametotvar(i,l),l=2,5)
         read(outp_atmosparam_id,*) (newtotvar(i,l),l=1,7)
         do l=1,6
            newtotvar(i,l)=newtotvar(i,l+1)
         enddo

         IF ( newtotvar(i,4)==1 ) ioutyearly = 1
         IF ((newtotvar(i,3)==1).OR.(newtotvar(i,2)==1) ) meantype = 1
         IF ( newtotvar(i,2)==1 ) meanyl     = 1
         IF ( newtotvar(i,3)==1 ) meantot    = 1
         IF ( newtotvar(i,1)==1 ) ioutdaily  = 1

         IF ((newtotvar(i,2)==1).OR.(newtotvar(i,3)==1)
     &                          .OR.(newtotvar(i,4)==1) ) THEN
         SELECT CASE ( nametotvar(i,5) )
         CASE ("T2")
            thirdd(1)=1
         CASE ("T3")
            thirdd(2)=1
         CASE ("T4")
            thirdd(3)=1
         CASE ("U3")
            thirdd(4)=1
         CASE ("N")
            l=0
         CASE DEFAULT
            call ec_error(123)
         END SELECT

         END IF

      enddo

      close(outp_atmosparam_id)

      return
      end

