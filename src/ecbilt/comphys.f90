!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module est la conversion en module FORTRAN90/95 du fichier
!       comatm.h
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche (pour la conversion)
!      Date   : 15 Mai 2018
!      Derniere modification : 16 Mai 2018
!-----|--1--------2---------3---------4---------5---------6---------7-|
!dmr -- Added optional components choice - Tue Jun 14 12:07:40 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 12:07:40 CEST 2011
!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** File:     comphys.h
! *** Contents: Common declarations for physical part of atmospheric
! ***           model of ECbilt
!      common /ec_cgamma/ gamgr,gamean,gamvar,gamsq
!      solarc:      solar constant.
!      common /ec_ctemag/ tempsg,temp2g,temp4g,tempsgn
!      tempsg:      mean temperature at 10 meters height.
!      tempsgn:     temperature at 10 meters height for each surface
!                   type.
!      temp2g:      temperature at 350 mb.
!      temp4g:      temperature at 650 mb.
!
!      common /ec_ctempm/ tempm,thform,temp2gm,temp4gm,
!                     tempsgm,tsurfm
!      tempm:       global mean temperature.
!      thform:      global mean
!      temp2gm:     global mean temperature at 350 mb.
!      temp4gm:     global mean temperature at 650 mb.
!      tempsgm:     global mean temperature at 10 meters height.
!      tsurfm:      global mean temperature at surface.
!
!      common /ec_cpot/   potm,pot2g,pot4g
!      potm:        global mean potential temperature.
!      pot2g:       potential temperature at 350 mb.
!      pot4g:       potential temperature at 650 mb.
!
!      common /ec_csurf/  tsurf,lsmask,sst,sstday,tsurfn,fractn
!      tsurf:       mean temperature at surface.
!      lsmask:      land sea mask.
!      sst:         sea surface tmeperature.???????????????
!      sstday:      sea surface temperature each day.?????
!      tsurfn       surface temperature for each surface type
!      fractn       fraction of each surface type
!      fracto       fraction of land/ocean
!
!      common /ec_sunfr/  solarc,Q0,solardref
!      solarc:      sun constant.
!
!      common /ec_cswrad/ hesw1,hesw2,hesws,albes,albea,albeaw,albeas,
!      *                abso1,abso2,albesn,heswsn
!      hesw1:       solar radiation heating rate.
!      hesw2:       solar radiation heating rate.
!      hesws:       mean solar radiation heating rate at the surface.
!      albes:       mean albedo of the surface
!      albea:       solar radiation reflective coefficient.
!      albeaw:      albedo of land in winter
!      albeas:      albedo of land in summer
!      abso1:       solar radiation absorbtion coefficient.
!      abso2:       solar radiation absorbtion coefficient.
!      albesn:      albedo of each surface type
!      heswsn:      solar radiation heating rate for each surface type
!
!      common /ec_cpar1/  rowat,roair,cdrag,cpair,cvair,rlatvap,rlatcon,
!      *                sboltz,rlatsub,rlatfus,cwater,gamma,rkappa,ps,
!      *                aa,bb,gamd,tzero,alphad,alphav,alphas,
!      *                dragan,dragla,cdragvn,epss,cwdrag
!      rowat:       density of water.
!      roair:       mean air desity at sea level.
!      cdrag:       coefficient in sensible and latent air-sea heat flux.
!      cwdrag:      coefficient in wind stress.
!      cdragv:      coef. in sen. en lat. heat flux depending on roughness
!                   length and stability (Richardson number)
!      richar:      richardson number
!      dragane:     rotation in the comp. of the wind stress as a function
!                   of latitude
!      cpair:       specific heat of dry air at constant pressure.
!      cvair:       specific heat of dry air at constant volum.
!      rlatvap:     latent heat uptake due to evaporation.
!      rlatcon:     latent heat release due to condensition.
!      sboltz:      stefan-boltzmann constant.
!      rlatsub:     latent heat of sublimation.
!      rlatfus:     latent heat of fusion.
!      cwater:      4180
!      gamma:       cpair/cvair
!      rkappa:      =(gamma-1)/gamma
!      ps:          surface pressure.
!      aa:          =(350/1000)**rkappa.
!      bb:          =(650/1000)**rkappa
!      gamad:
!      tzero:       =273.15
!      epss :       = treshold for surface computation
!      alphad:      =roair*cdrag*cpair
!      alphav:      =roair*cdrag*rlatvap
!      alphas:      =roair*cdrag*rlatsub
!      albes:       surface albedo.
!      dragan:      maximum turning angle in the computation of wind stress
!      dragla:      latitude below which the turning angle is reduced
!      cdragvn:     drag coefficeint for each surface type (see cdragv)

!      common/ fluxcore/ corAN,corPN,corAC,corID,corAS,corPS,corAA
!      corAN       flux correction in the North Atlantic
!      corPN       flux correction in the North Pacific
!      corAC       flux correction in the Arctic
!      corID       flux correction in the Indian
!      corAS       flux correction in the South Atlantic
!      corPS       flux correction in the South Pacific
!      corAA       flux correction in the Southern Ocean


!      common /ec_clwrad/ ulrad1,ulrad2,ulrads,dlrads,dlradsn,ulradsn
!      ulrad2:      upwards long wave radiation at ?.
!      ulrad1:      upwards long wave radiation at ?.
!      ulrads:      mean upwards long wave radiation at the surface.
!      ulradsn:     upwards long wave radiation at the surface.
!                   for each surface type
!      dlrads:      mean downwards long wave radiation at the surface.
!      dlradsn:     downwards long wave radiation at the surface
!                   for each surface type
!
!      common /ec_csflux/  Hflux(nlat,nlon),Eflux(nlat,nlon),dEflux(nlat,nlon),
!                       hfluxn(nlat,nlon,ntyps),eflux(nlat,nlon,ntyps)
!      Hflux:       mean sensible heat flux at the surface.
!      Eflux:       mean latent heat flux at the surfac.
!      dEflux:      latent heat flux at the surfac?????????????????.
!      hfluxn:      sensible heat flux for each surface type
!      efluxn:      latent heat flux for each surface type
!
!      common /ec_crain/  corain,dyrain,torain
!      corain;      convective rain.
!      dyrain:      dynamic rain.
!      torain:      total rain.
!
!      common /ec_cmoist/ qsats,qsat4,relhum,qg
!      qsats:
!      qsat4
!      relhum:      relative humidity.
!      qg
!
!      common /ec_moisturew/ relhmax, ihavm, ivavm, imsink
!      relhmax:     maximum relative humidity before saturation.
!      ihavm:       with (1) or without (0) horizontal divergence of moisture.
!      ivavm:       with (1) or without (0) vertical divergence of moisture.
!      imsink:      with (1) or without (0) the source or sink of moisture.
!
!      common /ec_cmois/  rmoiss,rmoisg,precip,evap,dcmoisg,evapn
!      rmoisg:      specific humidity.
!      precip:      (not found)
!      evap:        evaporation.
!      evapn:       evaporation for each surface type
!      dcmoisg:
!      common /ec_conv/   imoisr,imoisc,tetacr,gams,teta,convn
!      imoisr
!      imoisc
!      tetacr
!      gams
!      teta:        0.5*(pot2g - pot4g)
!      convn
!
!      common /ec_cthfor/ thforg1,thforg2,thfor1,thfor2,vhforg1,vhforg2
!      thforg1
!      thforg2
!      thfor1
!      thfor2
!      vhforg1:     heating force for 350 mb.
!      vhforg2:     heating force for 650 mb.
!
!      common /ec_cvorfor/ vfor1,vfor2,vfor3,vforg1,vforg2,vforg3
!      vfor1
!      vfor2
!      vfor3
!      vforg1:      vorticity forcing at 200 mb.
!      vforg2:      vorticity forcing at 500 mb.
!      vforg3:      vorticity forcing at 800 mb.
!
!      common /ec_cfluxpar/ evfac,evfaca
!      evfac:       maximum evaporation factor over land.
!      ecfaca:      actual evaporation factor over land.
!
!      common /ec_calbedo/ albice,albsea,albsnow,albland,albseaw,albseas,
!     *                abstow,abstos
!      albice:      albedo of ice.
!      albsea:      albedo of sea.
!      albsnow:     albedo of snow.
!      albland:     albedo of land.
!
!      common /ec_cstype/ noc,nse,nld
!      noc:       type number of the ocean
!      nse:       type number of the sea ice
!      nld:       type number of the land
!-----------------------------------------------------------------------

      MODULE comphys

      use global_constants_mod, only: dblp=>dp, silp=>sp, ip
#if ( HOURLY_RAD > 0 )
      use global_constants_mod, only: days_year360d_i
#endif
      use comatm, only: nlat, nlon, nsh2, nvl, nwisos

!~ #if ( ISOATM >= 1 )
!~       use iso_param_mod, only: neauiso
!~ #endif

      IMPLICIT NONE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      integer(ip), parameter :: iqmtab=50,jqmtab=20,kqmtab=20

      real(dblp), dimension(nlat,nlon) :: gamgr, temp0g, temp2g, temp4g
      real(dblp), dimension(0:2)       :: tempm
      real(dblp), dimension(nlat)      :: Q0
      !dmr [NOTA] dummy init for solardref in the following
      real(dblp)                       :: solarc, solardref = 1.0_dblp,ecc,obl,omweb,ecc2,so,perh, rowat,roair,cdrag,cwdrag,cpair  &
                ,cvair      &
                ,rlatvap,rlatcon,sboltz,rlatsub,rlatfus,cwater,gamma,rkappa,potfac1,potfac2,gamd,tzero,alphad,alphav,alphas        &
                ,corAN,corPN,corAC,corID,corAS,corPS,corAA
      real(dblp), dimension(nlat,nlon) :: cdragw,richar
      real(dblp), dimension(nlat)      :: dragane
      real(dblp), dimension(49)        :: chsea, cdsea, chlan, cdlan
      real(dblp), dimension(nlat,nlon) :: tetacr, gams, teta,thforg1,thforg2, thforg0, vhforg0, vhforg1=0.0_dblp, vhforg2=0.0_dblp &
                                        , relhum         &
                                        , vforg1, vforg2, vforg3
      real(dblp), dimension(nsh2)      :: vfor1, vfor2, vfor3

      real(dblp) :: evfac

!~ #if ( ISOATM >= 1 )
!~       real(dblp), dimension(nsh2,neauiso)     :: rmoiss
!~       real(dblp), dimension(nlat,nlon,neauiso):: rmoisg

!~ #else
      real(dblp), dimension(nsh2,nwisos)      :: rmoiss
      real(dblp), dimension(nlat,nlon,nwisos) :: rmoisg
!~ #endif

      real(dblp), dimension(nlat,nlon)        :: qmount

      real(dblp)                              :: dp2,dp1,dp0,tdifq,gpm500,relhmax,hmoisr,umoisr,rainmax
      integer(ip)                             :: ihavm, ivavm, imsink

!~ #if ( ISOATM >= 1 )
!~       real(dblp), dimension(nlat,nlon,neauiso):: corain, dyrain, cormois, torain, cosnow, dysnow, tosnow
!~ #else
      real(dblp), dimension(nlat,nlon,nwisos) :: corain, dyrain, torain, cosnow, dysnow, tosnow, cormois
!~ #endif
      real(dblp), dimension(nlat,nlon)                 :: co2
      real(dblp), dimension(0:iqmtab)                  :: tqmi
      real(dblp), dimension(0:jqmtab)                  :: tqmj
      real(dblp), dimension(0:kqmtab)                  :: tqmk
      real(dblp), dimension(0:iqmtab,0:jqmtab,0:kqmtab):: qmtabel
      real(dblp), dimension(nlat,nlon)                 :: tcc, tccd
      real(dblp), dimension(nlat,nlon,nvl)             :: dumt1,dumt2,dumu1,dumu2
      real(dblp), dimension(nlat,nlon)                 :: uv10,uvw10,utot10,vtot10

!dmr @-@ iceb0
! JONO march 2004
! analogue to uv10rws and uv10() introducing
! uv10rwv =  interpolation factor 800hPa to 10 meter
! default = 0.8 , to obtain:
! {utot(nlat,nlon),vtot()} the windvector used to force the icebergs
!dmr @-@ iceb0

      integer(ip):: ndayws,iradcloud,iscenghg,isghgstrt,iscentsi,iscenvol,istsistrt,isvolstrt,iscensul,issulstrt,isceno3,iso3strt  &
                   ,iscencel,iens,numens,irunlabeloffset,iscenyear

      real(dblp) :: cc1,cc2,cc3,tqmimin,dtqmi,rdtqmi,tqmjmin,dtqmj,rdtqmj,tqmkmin,dtqmk,rdtqmk
      real(dblp) :: facttsi,bup,eccf,oblf,omwebf,relhcrit,relhfac,emisoc,emisse,emisld,albin,albis,albice,alphd,alphdi,alphs,cgren
      real(dblp) :: dragan,dragla,uv10rfx,uv10m,uv10rws,uv10rwv

!***  common /ec_swrscheme/costref,salbref,swrref,swrcost,swrsalb,dayfr(nlat),
!***          kosz(nlat),solarf(nlat)
!***  costref contains regional averaged cosine of the zenith angle
!***  salbref contains regional averaged surface albedo from ????????
!***  swrref  contains reference short wave radiation fluxes from
!***          radiation calculated using KRCM with NCEP humidity, ipcc95
!***          greenhousegascontrations, isccpd2 cloud climatology and
!***          ECHAM4 ozone climatology
!***  swrcost linear sensitivity coefficient for SWR changes due to anomalies
!***          from the cosine of the zenith angle wrt the reference value
!***  swrsalb linear sensitivity coefficient for SWR changes due to anomalies
!***          from the surface albedo for clear sky and 3rd order polynomial fit
!***          for unity overcast fluxes.
!*** index 1: flux levels
!*** index 2: regions
!*** index 3: month's
!*** index 4: 0 clear sky, 1 cloudy sky, except for surface albedo: 1, 2,3
!***          correspond to coefficients of 3rd order polynomial fit
!*** region classification and flux level definition same as LWR


      real(dblp), dimension(27,12)           :: costref, salbref
      real(dblp), dimension(8,27,12,0:1)     :: swrref, swrcost
      real(dblp), dimension(8,27,12,0:3)     :: swrsalb

      real(dblp), dimension(nlat)            :: dayfr = 0.5_dblp, kosz=1.0_dblp, solarf= 1350.0_dblp ! [NOTA] dummy fix with the inits
      real(dblp), dimension(nlat,nlon)       :: dso4 = 0.0_dblp, tas1 ! [NOTA] dummy fix with the inits
      real(dblp), dimension(3000,nlat,nlon)  :: sulopt

!***  common lwrscheme/tncep,qancep,ghgipcc,lwrref,lwrt,lwrts,lwrqa,lwrghg,irn

!***  tncep contains reference vertical temperature profiles from the
!***  ncep reanalysis (1982-1994) for 12 month's, 27 regions and 19 levels:
!***  10,20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000,p2m,psurf
!***  except when surface pressure lower than 1000, then p2m is inserted
!***  at the appropriate position as given by ipl (see below)

!***  qancep:total precipitable water content below 500 mb from ncep, reduced
!***  in order to tune the LW fluxes (minus 15 %)

!***  ghgref: green house gas concentrations 1990 from ipcc '92:
!           1: pCO2-CO2 conc. (ppmv)
!           2: pCH4-methane conc. (ppbv)
!           3: pN2O-nitrous-oxide (ppbv)
!
!           pCFC-concentrations CFCs (pptv) per CFC type (pptv):
!           4: pCFC(1) - CFC-11
!           5: pCFC(2) - CFC-12
!           6: pCFC(3) - CFC-113
!           7: pCFC(4) - CFC-114
!           8: pCFC(5) - CFC-115
!
!           pHCFC-concentrations HCFCs (pptv) per type (pptv):
!           9: pHCFC(1) - HCFC-22
!          10: pHCFC(2) - HCFC-123
!          11: pHCFC(3) - HCFC-124
!          12: pHCFC(4) - HCFC-125
!          13: pHCFC(5) - HCFC-134A
!          14: pHCFC(6) - HCFC-141B
!          15: pHCFC(7) - HCFC-142B
!          16: pHCFC(8) - HCFC-143A
!          17: pHCFC(9) - HCFC-152A
!
!          18: pCTC-concentration Carbon TetraChloride (CCl4) (pptv)
!          19: pMCF-concentration Methyl Chloroform (CH3CCl3) (pptv)

!***  O3echam4(kg/kg): based on ECHAM4 climatology; pressure levels in file
!***  for each season: djf-mam-jja-son
!***  ccisccp: total cloud cover monthly climatology: derived from isccp, used
!***  in the linearisation procedure of KRCM, tuned
!***  (middle and high clouds -15%) for KRCM to give 'observed'  LWR fluxes.

!***  lwrref: reference long wave radiation fluxes, for four month's, jan,
!***  apr, july, oct, calculated with KRCM with
!***  data from several sources: temperature and humidity from ncep, ground
!***  pressure and cloud cover climatology from ISCCP, ozone from ECHAM4:
!***  values of the first index represent:  1: OLR,  2: upward at 200 mb,
!***  3: upward at 500 hPa, 4: upward at surface, 5: downward at 200 mb,
!***  6: downward at 500 mb, 7: downward at surface
!***  and fourth index corresponds to clearsky (0) and unity overcast (1).
!***  lwrrefc: same as lwrref but corrected for systematic difference of
!***  linearised scheme wrt KRCM for NCEP dataset 1982-1994
!***  lwrt: first partial derivative of lwr fluxes wrt temperature at different
!***  levels (see tncep and lwrref)
!***  lwrts: 4th order polynomial fit of lwr flux dependence on surface temperature
!***  lwrqa: first partial derivative of lwr fluxes wrt total precipitable water
!***  content, distributed according to  ncep mean vertical profiles
!***  lwrghg: fits of lwr flux dependence on several ghg concentrations
!***  irn: index of regions used for definition of the reference profiles and fits
!***  27: Green Land area
!***  26: Rocky Mountain area
!***  25: Himalaya area
!***  24: Andes area
!***  23: Antarctica area
!***  22: zonal mean 75N-90N land
!***  21: zonal mean 75N-90N sea
!***  20: zonal mean 60N-75N land
!***  19: zonal mean 60N-75N sea
!***  18: zonal mean 45N-60N land
!***  17: zonal mean 45N-60N sea
!***  16: zonal mean 30N-45N land
!***  15: zonal mean 30N-45N sea
!***  14: zonal mean 15N-30N land
!***  13: zonal mean 15N-30N sea
!***  12: zonal mean 15S-15N land
!***  11: zonal mean 15S-15N sea
!***  10: zonal mean 30S-15S land
!***   9: zonal mean 30S-15S sea
!***   8: zonal mean 45S-30S land
!***   7: zonal mean 45S-30S sea
!***   6: zonal mean 60S-45S land
!***   5: zonal mean 60S-45S sea
!***   4: zonal mean 75S-60S land
!***   3: zonal mean 75S-60S sea
!***   2: zonal mean 90S-75S land
!***   1: zonal mean 90S-75S sea
!***  ipl: index of pressure level in tncep containing 2mtr
!***       temperature
!***  pisccp: surface pressure anual mean, region averaged

      real(dblp), dimension(19,27,12) :: tncep
#if ( F_PALAEO >= 1 )
      real(dblp), dimension(20,798001):: ghgscen
#else
      real(dblp), dimension(20,10100) :: ghgscen
#endif
      real(dblp), dimension(27,12)    :: qancep
      real(dblp), dimension(19)       :: ghgipcc

      integer(ip)                     :: y1scenghg,nyscenmax


      real(dblp), dimension(18,nlat,nlon,2):: dtemp
      real(dblp), dimension(27,12)         :: z500ncep
      real(dblp), dimension(27)            :: rlogts
      real(dblp), dimension(27)            :: pisccp
      real(dblp), dimension(20)            :: ghg
      real(dblp), dimension(17)            :: rlogtl
      real(dblp), dimension(17)            :: pncep


      real(dblp), dimension(32,64,12)     :: ccisccp ![NOTA] is 32, 64 == nlat, nlon ??
      real(dblp), dimension(2,27,12)      :: tncep12
      real(dblp), dimension(7,27,4,0:1)   :: lwrref
      real(dblp), dimension(7,27,4,0:1,2) :: lwrflux
      real(dblp), dimension(7,18,27,4,0:1):: lwrt
      real(dblp), dimension(7,4,27,4,0:1) :: lwrts
      real(dblp), dimension(7,27,4,0:1)   :: lwrqa
      real(dblp), dimension(7,4,27,4,0:1) :: lwrqts
      real(dblp), dimension(7,19,27,4,0:1):: lwrghg

      integer(ip), dimension(nlat,nlon,2) :: irn
      integer(ip), dimension(27)          :: ipl

      real(dblp)                          :: AMPWIR,AMPEQIR,EXPIR,HPROFTROP,HPROFEQ,HPROFW,HPROFAN,AMPANIR,HPROFAN2,AMPANIR2 &
                                           , mag_alpha
      real(dblp), dimension(10100)        :: solartsi
      real(dblp), dimension(10100,12,4)   :: solarvol
      real(dblp)                          :: solarm
      real(dblp), dimension(nlat)         :: solarcl

      real(dblp), dimension(nlat,nlon)    :: fswdsfc,fswusfc

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#if ( HOURLY_RAD > 0 )
! dmr in the following, the first dimension is the number of steps per day
      integer(ip), parameter                       :: n_steps_day = 6
      real(dblp), dimension(n_steps_day,nlat)      :: solflux_hourly
#endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



      END MODULE comphys
