!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:29 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:29 CET 2009

c- issu de la routine "TVM" ; modif (JMC) : 08/10/00; inclu
c- dans emic (driess) :30/08/2002
*********************************************************************
*      Terrestrial Vegetation annual Model
*
*  Purpose: Computation of forest/grass/desert ratio;
*           npp and living and soil biomass
*
*  By: V.Brovkin
*  Modifyed by: A.Ganopolski
*  Last modification: 17.11.97
*  Started cleanup, dmr, 2019-10-17
**********************************************************************
c- nyvegt = frequence (in Years) of "call veget"
c      if = 0 --> return flgveg=False and do not "call veget" anymore.
c- nwrveg = frequence (in Years) of ASCII output (0 => no output)
c---------
c-First call of "veget" (iyear=0) : read "veget.par"
c kveget < 0  => compute initial vegetation using T,Pre,GDD0 from "veget.init"
c kveget > 0  => read vegetation on "veget.rest"
c   after 1rst call : kveget =abs(kveget)  and other param. are kept constant
c-In Running Loop (iyear > 0) :
c kveget = 0  => no effect on Atmospheric Model (Albedo not changed)
c         but write T_moy,Pre_moy & GDD0 on file "veget.init" & "veget.outp"
c kveget = 1  => synchronous run of vegetation model
c kveget > 1  => do "kveget" iterations of vegetation to get faster equilibrium
c----------
c   prcmin= minimun daily precip (in m) available for vegetation in warm area.
c bmtdry = Soil-Water threshold(m) to compute "tpsdry"=Nb_days/yr of dry soil C.
c tmxdry : above this limit(Nb_days/yr), soil_dryness start to affect vegetation
c dtrdry = time interval(days) for linear transition [tmxdry,tmxdry+dtrdry] ;
c  If tpsdry(=Nb_day_dry) > tmxdry+dtrdry => Precip(Veget) is reduced by rpfdry
c rpfdry = reduction factor applied to Precip(Veget) if dry condition satisfied

*********************************************************************
c input data: annual mean temperature in degr. Celc. - ave_t
c             annual mean precipitation, mm/year - ave_pr
c             growing degree days above 0, degr. Celc - gdd0
c             CO2 concentration in the atmosphere, ppm - co2
c
c output data: st(i,j) - forest's share in cell i,j
c              sg(i,j) - grass share in cell i,j
c              sd(i,j) - desert's share in cell i,j
c              snlt(i,j) -needle leaved trees ratio
c              st(i,j)+sg(i,j)+sd(i,j)=1, 0<= snlt(i,j)<= 1
c              plai(i,j) - annual average of leaf area index, m2/m2
c              anup(i,j) - annual uptake of carbon in cell, Gt -> kg/m2 ???
c
c COMMON block described in buffer.inc
*********************************************************************

      SUBROUTINE veget(ist,jst,dtime,epss,patmCO2,fracgr,darea,tempgr)
**************************************************************************
c~ [WEIRD] apparently unused ... 2019-10-16
c~ #if ( IMSK == 1 )
c~ cdmr --- Added for the LGM ... : prevent vegetation over the ice-sheets !
c~       USE input_icemask, only:
c~ cdmr --- Added for the LGM ... : prevent vegetation over the ice-sheets !
c~ #endif

      use global_constants_mod, only: dblp=>dp, ip

#if ( CYCC == 2 )
        USE veget_iso,   only: b1g13, b1g14, b1t13, b1t14, b2g13, b2g14
     &               , b2t13, b2t14, b3g13, b3g14, b3t13, b3t14, b4g13
     &               , b4g14, b4t13, b4t14

        USE carbone_co2, ONLY: new_run_c
#if ( KC14 == 1 )
        USE C_res_mod, ONLY: cav_la14_b, cav_la_b
#endif
#endif
#if ( ISOATM >= 1 )
        USE iso_param_mod
#endif

#if( PF_CC > 0 )
        USE LECT_2D_NC__genmod
#endif

      use comland_mod, only: bmoisg,bmoism,rainf,snowf, albsnow,albland,
     &  forestfr,alblbm

      use comatm,      only: nlat, nlon, iwater
      use comphys,     only: tosnow, torain
      use comdiag,     only: irad
      use comemic_mod, only: iyear, nyears, flgveg, iatm, nstpyear, nwrskip
      use veget_mod,   only: anup, ave_pr, ave_pr05, ave_t, b1, b1g, b1t
     &             , b2, b2g, b2t, b3, b3g, b3t, b4, b4g, b4t, betag
     &             , betat, co2ghg, fco2veg, gamm2, gdd0, i0dfor, ieqveg
     &             , ieqvegwr, iscendef, iveg, ivegstrt, lat, lon, mdfor
     &             , ndfor, numvegvar, sd, sdr, sg, sgr, snlt, st
     &             , st_const, str, veg_fill_value, veg_missing_value
     &             , phi, namevegvar, newvegvar, farea, snltr, sd_const
     &             , stock
      use comrunlabel_mod, only: irunlabelf
      use veget_sub_mod,   only: ccstat, ccdyn, ccstatR, ccdynR, initcpar

      use newunit_mod, only: info_id
      
      implicit none

      include 'netcdf.inc'
      real*8 veg_albd
*********************************************************************
!      real*8    albsnow(nlat),albland(nlat,nlon,4)
!#if ( ISOATM >= 2 )
!      real*8    bmoisg(nlat,nlon,5),forestfr(nlat,nlon)
!#else
!      real*8    bmoisg(nlat,nlon),forestfr(nlat,nlon)
!#endif
!      real*8    bmoism(nlat,nlon)
C     real*8    rs(nlat,nlon)
!#if ( ISOATM >= 2 )
!      real*8    alblbm(nlat,nlon),rainf(nlat,nlon,5),snowf(nlat,nlon,5)
!#else
!      real*8    alblbm(nlat,nlon),rainf(nlat,nlon),snowf(nlat,nlon)
!#endif
!      real*8    alblbmR(nlat,nlon),alblandR(nlat,nlon,4)
!      real*8    alblandismR(nlat,nlon,4,3),forestfrR(nlat,nlon)
!      common /ec_lbmbmois/ bmoisg,bmoism,rainf,snowf
!      common /ec_lbmcalbedo/albsnow,albland,forestfr,alblbm
!      common /ec_lbmcalbedo2/alblandR,forestfrR,alblbmR,alblandismR

c--dummy variables :
c- input (except 1rst call) :
      integer ist,jst, istep
c~ #if ( WRAP_EVOL == 2 )
c~       integer iyearvegout
c~ #endif
      real*8 dtime, epss
      real*8 fracgr(nlat,nlon), darea(nlat), tempgr(nlat,nlon)
c- output :
c~ [WEIRD] not used ...
c~ #if ( WRAP_EVOL == 2 )
c~       character*6 endyear
c~       character*3 endday
c~       logical existe
c~ #endif
c--local variables saved (in common) from one call to the other :
      logical, save         ::  flgwrt
      
      integer(kind=ip),save ::  kveget, nyvegt, nwrveg, ncumvg, iyr0vg
      real(kind=dblp), save ::  tmxdry, dtrdry, rpfdry, unsdry, prcmin
     &               , bmtdry
      
      integer(kind=ip),dimension(4,0:1), save :: ns0inv


      integer ieq

c~       real*8 temveg(nlat,nlon), gd0veg(nlat,nlon),
c~      &     prcveg(nlat,nlon,2), tpsdry(nlat,nlon)
c~ cdmr --- Ajout du GDD5
c~       real*8 gd5veg(nlat,nlon)
c~ cdmr --- Ajout du GDD5

      real(kind=dblp), dimension(nlat,nlon), save :: temveg,gd0veg
     &               , tpsdry, gd5veg, freezeI, thawI, prcday = 0.0d0

      real(kind=dblp), dimension(nlat,nlon,2), save :: prcveg



      real*4  outdata(64,32)

c~       real*8 prcmin, bmtdry, prcday(nlat,nlon)
      real*8 fareal(nlat,nlon)
      
      character(len=15), save ::  titveg
      character(len=15) :: titv00

!      common / com0veg / kveget, nyvegt, nwrveg, ncumvg,
!     &                   iyr0vg, ns0inv(4,0:1)
!      common / com1veg / tmxdry,dtrdry,rpfdry,unsdry,
!     &                   temveg,gd0veg,prcveg,tpsdry
!cdmr --- Ajout du GDD5
!     &                   , gd5veg,
!cdmr --- Ajout du GDD5
!     &                   freezeI, thawI

c~       common / com2veg / prcmin, bmtdry, prcday
c~       common / comcveg / titveg


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c--local variables :
      integer i,j,k,ns,ii,kk2,nlatd2, iyrloc, ireg
      real*8 zero, one, xx, xxalb, zlai, ddc, bmtd00
      real*8 dzero, ttemp, prcrit, zmoy1, zmoy2, zmoy3, zmoy4
      real*8, SAVE :: albet(4), albeg(4), albed(4), albegc(4), albedc(4)
      integer istart,igroen
      integer ios,ios2
      double precision patmCO2

!dmr --- added by Huan Li to force the vegetation cover to the CMIP LUH dataset
!dmr --- regardless of what is (still) computed by VECODE
#if ( VEG_LUH == 1 )
      REAL*8, DIMENSION(nlat,nlon) :: st_luh
      REAL*8, DIMENSION(nlat,nlon) :: sg_luh
      REAL*8, DIMENSION(nlat,nlon) :: sd_luh
#endif

cdmr --- Ajout du GDD5
      real*8 dcinq
cdmr --- Ajout du GDD5

      character*60 part1
      character(len=256) :: output_filename
      real*8 ari(nlat/2)
      real*8 rj,rw,dumwei
      integer ilat,ird, status


#if ( PF_CC > 0 )
       REAL, DIMENSION(nlat,nlon) :: dataread
       REAL, DIMENSION(nlat,nlon) :: dataecb


       character(len=4), parameter :: varia="PF"
       character(len=200), parameter ::
     & pathtofile = "inputdata/permafrost/Only_Permafrost-aoT.nc"

#endif

!--- IDs for files
       integer :: vegetpar_id, gaussasc_id, outp_vegetparam_id
     >          , vegetdat_id, vegetrest_id, luh2_vegtxt_id
     >          , vegetinit_id, veget_accurest_id, freezeIctl_id
     >          , freezeIdat_id, permafctl_id, permafdat_id


c- Albedo of tree, grass, dessert for Winter,Spring,Summer and Fall:
c     data albet / 0.13 , 0.13 , 0.13 , 0.13 /
c     data albeg / 0.20 , 0.20 , 0.20 , 0.20 /
c     data albed / 0.33 , 0.33 , 0.33 , 0.33 /
c- Albedo of grass over cold area (Steppe) = albeg + albegc
c     data albegc / -.06 , -.04 , -.02 , -.04 /
c- Albedo of bright sand desert (Sahara) = albed + albedc
c     data albedc / 0.07 , 0.07 , 0.07 , 0.07 /
c--For output :
      data flgwrt / .TRUE. /

c~ #if ( WRAP_EVOL == 2 )
c~       logical veg_act
c~ #endif


       logical :: returnValue


      zero = 0.
      one  = 1.
      dzero = 0.d0
cdmr --- Ajout du GDD5
      dcinq=5.d0
cdmr --- Ajout du GDD5

c~ #if ( WRAP_EVOL == 2 )
c~       if (iyear.eq.0 .and. ist.eq.1 .and. jst.eq.1 .and. iday.eq.30 ) then
c~              iyearvegout = 0
c~       elseif (iyear.eq. 0 ) then
c~              iyearvegout = 1
c~       else
c~              iyearvegout = iyear
c~       endif
c~ #endif
c-----

c~ #if ( WRAP_EVOL == 2 )
c~       if(initialization.eqv..true.) then
c~ #else
      if (iyear.eq.0) then
c~ #endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  0 ) Read "veget.par" and Define parameter for this routine.         |
c-----------------------------------------------------------------------

c--Initialization of Cumul.array

c--ouverture de "veget.uptake" pour le calcul du flux de carbon
c--entre l atmosphere et la biomasse

c      open(iveg+3,file='veget.uptake',status='unknown')

c- Read Vegetation parameter :
      open(newunit=vegetpar_id,file='veget.par',status='old')
      read(vegetpar_id,*)
      read(vegetpar_id,*)
      read(vegetpar_id,'(A)') titveg
      read(vegetpar_id,*)
      read(vegetpar_id,*) kveget, nyvegt, nwrveg
      read(vegetpar_id,*)
      read(vegetpar_id,*) prcmin, bmtdry, tmxdry, dtrdry, rpfdry
      read(vegetpar_id,*)
      read(vegetpar_id,*) ieqveg,ieqvegwr
      read(vegetpar_id,*)
      read(vegetpar_id,*) iscendef ,ivegstrt
      read(vegetpar_id,*)
      read(vegetpar_id,*) (albet(i),i=1,4)
      read(vegetpar_id,*)
      read(vegetpar_id,*) (albeg(i),i=1,4)
      read(vegetpar_id,*)
      read(vegetpar_id,*) (albed(i),i=1,4)
      read(vegetpar_id,*)
      read(vegetpar_id,*) (albegc(i),i=1,4)
      read(vegetpar_id,*)
      read(vegetpar_id,*) (albedc(i),i=1,4)
      read(vegetpar_id,*)
      read(vegetpar_id,*) gamm2
      read(vegetpar_id,*)
      read(vegetpar_id,*) betag,betat
      read(vegetpar_id,*)
      read(vegetpar_id,*) fco2veg
      close(vegetpar_id)

!dmr --- Following section has been added by Huan Li (ported into the trunk 2016-12-14)
!dmr     It fixes a small bug in the netcdf output of the vegetation model

chli--- calculating latitude based on Gauss points,replacing iuo+7 by 111, commented the WRAP_EVOL syntax

      open(newunit=gaussasc_id, file='inputdata/gauss.asc',status='old',
     &     form='formatted')
      rewind(gaussasc_id)
      ilat=nlat/2
   10 continue
        read(gaussasc_id,220,end=15) rj,rw
        ird=int(rj)
        if (ird.eq.ilat) then
          do i=1,ird
            read(gaussasc_id,220) ari(i),dumwei
          enddo
          goto 20
        else
          goto 10
        endif
   15 continue
C     call ec_error(4)
   20 continue

      do i=1,ilat
        phi(i)=-ari(ilat+1-i)
        phi(ilat+i)=ari(i)
      enddo

      do i=1,nlat
        phi(i)=asin(phi(i))
      enddo

 220  format(f18.10,f17.10)
      close(gaussasc_id)

!dmr --- End of added section from Huan Li

c~ #if ( WRAP_EVOL == 2 )
c~ C *** gauss points and weights
c~       rewind(iuo+7)
c~       ilat=nlat/2
c~    10 continue
c~         read(iuo+7,220,end=15) rj,rw
c~         ird=int(rj)
c~         if (ird.eq.ilat) then
c~           do i=1,ird
c~             read(iuo+7,220) ari(i),dumwei
c~           enddo
c~           goto 20
c~         else
c~           goto 10
c~         endif
c~    15 continue
c~ C     call ec_error(4)
c~    20 continue

c~       do i=1,ilat
c~         phi(i)=-ari(ilat+1-i)
c~         phi(ilat+i)=ari(i)
c~       enddo

c~       do i=1,nlat
c~         phi(i)=asin(phi(i))
c~       enddo

c~  220  format(f18.10,f17.10)

c~ #endif
c--- reading outp_veget.param, !mohr
c--- here some variables for netcdf file are defined
      open(newunit=outp_vegetparam_id,file='outp_veget.param')
      read(outp_vegetparam_id,'(/,/,/,/,/,/,/,/)')
      read(outp_vegetparam_id,*) numvegvar,veg_fill_value
     >                          ,veg_missing_value
      read(outp_vegetparam_id,*)

      do i=1,numvegvar

         read(outp_vegetparam_id,"(A)") part1
         namevegvar(i,1)=trim(part1)
         read(outp_vegetparam_id,*) (namevegvar(i,k),k=2,5)
         read(outp_vegetparam_id,*) (newvegvar(i,k),k=1,7)
         do k=1,6
            newvegvar(i,k)=newvegvar(i,k+1)
         enddo

C         If ( newvegvar(i,4)==1 ) ioutyearly = 1
C         IF ((newvegvar(i,3)==1).OR.(newvegvar(i,2)==1) ) meantype = 1
C         IF ( newvegvar(i,3)==1 ) meanyl     = 1
C         IF ( newvegvar(i,2)==1 ) meantot    = 1
C         IF ( newvegvar(i,1)==1 ) ioutdaily  = 1

      enddo

      close(outp_vegetparam_id)

c-AM
C beta_i divided once  for all by log(2)
      betat=betat/log(2.0)
      betag=betag/log(2.0)

c-AM (2008)
c***  read deforestation
      i0dfor=0
      ndfor=0
#if ( WRAP_EVOL == 2 )
c-YSD (2011)
c***  Netcdfisation
      if (iscendef.eq.1) then
        status=nf_open("inputdata/VEGET.nc", nf_nowrite, ireg) !ouvre le fichier
        status=nf_inq_dimid(ireg, 'time', i) !recupère l'id du temps
        status=nf_inq_dimlen(ireg, i, ndfor) !recupère à partir de l'id, le nombre de pas de temps
        if(ndfor.gt.mdfor) STOP 'Please adjust  mdfor  in veget.h'
        status=nf_inq_varid(ireg, 'time', i) !recupère l'id du temps
        status=nf_get_vara_real(ireg, i, (/1/), (/ndfor/), VegetTime)
        ivegstrt=int(VegetTime(1)) !la 1ere valeur du temps donne la date de debut du forcage
        status=nf_inq_varid(ireg, "vegfrac", j) !récupère l'id de la variable Sul
        status=nf_get_vara_real(ireg, j, (/1,1,1/), (/64,32,ndfor/), farea) !charge les valeurs dans la variable sulopt
        ivegstrt=max(i0dfor,ivegstrt)
        write(info_id,*) 'scen Veget start=',ivegstrt, "AD nbline=",ndfor
      endif
#else
      if (iscendef.eq.1) then
        open(newunit=vegetdat_id,file='inputdata/VEGET.dat'
     >      ,form='unformatted')
        read(vegetdat_id) ndfor,i0dfor
        if(ndfor.gt.mdfor) STOP 'Please adjust  mdfor  in veget.h'
        if(ndfor.le.0) STOP 'Wrong HEADER in VEGET.dat'
c
        do k=1,ndfor
            read(vegetdat_id) fareal
            farea(:,:,k)=fareal(:,:)
        enddo
        close(vegetdat_id)
c*** if deforestation=VEGET.dat ref year should not be less than i0dfor
        ivegstrt=max(i0dfor,ivegstrt)
      endif
#endif
c
c-AM (2008)
c  in case of constant vegetation distribution
      if (iscendef.eq.-1) i0dfor=ivegstrt

c- check param :
      unsdry = 0.0
      if (dtrdry.gt.1.0e-6) unsdry = 1.0 / dtrdry

      if (nyvegt.eq.0) then
        nyvegt = nyears + 1
        kveget = abs(kveget)
        flgveg = .FALSE.
        return
      endif
      nwrveg = (nwrveg/nyvegt) * nyvegt
      if (nwrveg.eq.0) then
        nwrveg = nyears + 1
        flgwrt = .FALSE.
      endif
      iyr0vg = 0

#if ( PF_CC > 0 )
!dmr --- so far permafrost is not computed interactively. Thus, I read
!it ...

       CALL lect_2D_NC(TRIM(pathtofile),dataread,var_nc = TRIM(varia))

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Reading the input file leads to an array indexed differently as
!      ours. We thus reformat it.
!      NOTA: Not necessary here since it is read directly as lat,lon
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       CALL format_nc_to_ec(dataread,dataecb,nlon,nlat)

       dataecb = dataread

       WHERE(dataecb.LT.0.0)
         dataecb = 0.0
       ENDWHERE

       WRITE(*,*) "Fin lecture permafrost", MAXVAL(dataecb),
     & MINVAL(dataecb)
#endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Before the 1rst time step of the run : Setup Initial Vegetation.|
c-----------------------------------------------------------------------

c--reverse season index in South.H. = ns0inv(ns,0)
      do ns=1,4
        ns0inv(ns,0) = 1 + mod(ns+1,4)
        ns0inv(ns,1) = ns
      enddo

c- Initialisation of Vegetation Parameters :
      if (kveget.ne.0) returnValue=initcpar()

      if (kveget.gt.0) then
c- Read Vegetation from file "veget.rest".
      open(newunit=vegetrest_id,file='startdata/veget.rest',status='old'
     >     , form='unformatted')
      read(vegetrest_id) iyr0vg,bmtd00,titv00
      read(vegetrest_id) st
      read(vegetrest_id) sg
      read(vegetrest_id) sd
      read(vegetrest_id) snlt
      read(vegetrest_id) temveg
      read(vegetrest_id) b1t
      read(vegetrest_id) b1g
      read(vegetrest_id) b2t
      read(vegetrest_id) b2g
      read(vegetrest_id) b3t
      read(vegetrest_id) b3g
      read(vegetrest_id) b4t
      read(vegetrest_id) b4g
      read(vegetrest_id) b1
      read(vegetrest_id) b2
      read(vegetrest_id) b3
      read(vegetrest_id) b4
      read(vegetrest_id) anup
#if ( CYCC == 2 )
      IF (new_run_c.EQ.0) THEN
        read(vegetrest_id) B1T13
        read(vegetrest_id) B1G13
        read(vegetrest_id) B2T13
        read(vegetrest_id) B2G13
        read(vegetrest_id) B3T13
        read(vegetrest_id) B3G13
        read(vegetrest_id) B4T13
        read(vegetrest_id) B4G13
        read(vegetrest_id) B1T14
        read(vegetrest_id) B1G14
        read(vegetrest_id) B2T14
        read(vegetrest_id) B2G14
        read(vegetrest_id) B3T14
        read(vegetrest_id) B3G14
        read(vegetrest_id) B4T14
        read(vegetrest_id) B4G14
#if ( KC14 == 1 )
        read(vegetrest_id) cav_la14_b
        read(vegetrest_id) cav_la_b
#endif
      ELSE

        CALL ccstat_isotope()

      ENDIF
#endif
cdmr --- Ajout dans le cas de valeurs aberrantes dans le fichier de
cdmr --   redemarrage

       WHERE(b4t.LT.0.0)
        b4t = 0.0
       ENDWHERE

       WHERE(b4g.LT.0.0)
        b4g = 0.0
       ENDWHERE

       WHERE(b4.LT.0.0)
        b4 = 0.0
       ENDWHERE

!dmr --- following section added by Huan Li to force the vegetation to LUH2 dataset
#if ( VEG_LUH == 1 )

c read the fraction of tree, grass and desert from LUH2 dataset
      open(newunit=luh2_vegtxt_id,file = 'startdata/luh2_veg.txt'
     >    ,status='old', form="formatted")
      do i=1,nlat
        do j=1,nlon
           read(luh2_vegtxt_id,*) st_luh(i,j), sg_luh(i,j),sd_luh(i,j)
        enddo
      enddo
      close(luh2_vegtxt_id)

c- pass fraction of tree, grass and desert from LUH2 dataset to vecode
      st(:,:)= st_luh(:,:)
      sg(:,:)= sg_luh(:,:)
      sd(:,:)= sd_luh(:,:)
#endif

c
c iscendef = -1 : veget does not evolve (maintained at its initial distribution)
c          = +1 : land-use scenario
c          =  0 : no constraint
c
!dmr >>> Ajout pour compilation -C -fpe0
      ios2 = 0
!dmr <<<
      if ((iscendef.eq.1).or.(iscendef.eq.-1)) then
c Reference vegetation (note no need to save sd_const <- AM)
c ivegstrt define the reference year for anomalies of cropland
       read(vegetrest_id,end=99,iostat=ios) st_const
   99  if ((ios.ne.0).or.(irunlabelf.eq.ivegstrt)) st_const(:,:)=st(:,:)
       sd_const(:,:)=sd(:,:)
       read(vegetrest_id,end=98,iostat=ios2)str
       read(vegetrest_id,end=98,iostat=ios2)sgr
       read(vegetrest_id,end=98,iostat=ios2)sdr
      endif
   98 if ((ios2.ne.0).or.(iscendef.ne.1)) then
       sdr(:,:)=sd(:,:)
       sgr(:,:)=sg(:,:)
       str(:,:)=st(:,:)
      endif
      close(vegetrest_id)
!dmr ### [DELETE] -- to avoid the creation of a fort.67 file ...      
!dmr ???      write(67,'(A)') 'Initial Vegetation <- Read file "veget.rest".'
!dmr ???      write(67,'(A,I6,A,F6.3,2A)') ' year:', iyr0vg,
!dmr ???     &        ' Dry.Soil Lim=', bmtd00, ' ; title: ', titv00
!dmr ### [DELETE] 
      snltr(:,:)=snlt(:,:)

      do i=1,nlat
        do j=1,nlon
          stock(i,j)=b1(i,j)+b2(i,j)+b3(i,j)+b4(i,j)
c~ [DEPRECATED] these variables where intended for the LOCH carbon cycle
c~           stockloch(i,j)=fracgr(i,j)*darea(i)*stock(i,j)*1E-12
c~           anuploch(i,j)=fracgr(i,j)*darea(i)*anup(i,j)*1E-12
c~           anuploch(i,j)=fco2veg*anuploch(i,j)
        enddo
      enddo
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-----
      elseif (kveget.lt.0) then
c- Read from file=veget.init : An_Mean_Temp, gdd0 & An_Mean_Precip.
      open(newunit=vegetinit_id,file='startdata/veget.init',status='old'
     >    ,form='unformatted')
      read(vegetinit_id) temveg
      read(vegetinit_id) gd0veg
chli---read gd5veg from file veget.init
      read(vegetinit_id) gd5veg
      read(vegetinit_id) ((prcveg(i,j,1),i=1,nlat),j=1,nlon)
      read(vegetinit_id) ((prcveg(i,j,2),i=1,nlat),j=1,nlon)
      read(vegetinit_id) tpsdry
      read(vegetinit_id) iyrloc, bmtd00, titv00
      close(vegetinit_id)
!dmr ### [DELETE] -- to avoid the creation of a fort.67 file ... 
!dmr ???      write(67,'(A)') 'Initial Veg. computed from file "veget.init" :'
!dmr ???      write(67,'(A,I6,A,F6.3,2A)') ' year:', iyrloc,
!dmr ???     &        ' Dry.Soil Lim=', bmtd00, ' ; title: ', titv00
!dmr ### [DELETE]       
c-----
      endif

      else
c------------------------------------------------
c- during time integration : <--> kveget > or = 0


*************************************************************************
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) current time step : cum. atmospheric fields used for vegetation :
c------------------------------------------------

      istep=(ist-1)*iatm+jst

      ncumvg = ncumvg + 1
      do j=1,nlon
       do i=1,nlat
         ttemp = tempgr(i,j)-273.15d0
         temveg(i,j)=temveg(i,j)+ttemp
c- torain is in m/s ; dtime = length of 1 iter, in second
c~ #if ( ISOATM >= 1 )
c~        prcday(i,j)=prcday(i,j)+dtime*(torain(i,j,ieau)+tosnow(i,j,ieau))
c~ #else
         prcday(i,j)=prcday(i,j)+dtime*(torain(i,j,iwater)
     &              +tosnow(i,j,iwater))
c~ #endif
         gd0veg(i,j)=gd0veg(i,j)+max(ttemp,dzero)
cdmr --- Ajout du GDD5
         gd5veg(i,j)=gd5veg(i,j)+max(ttemp-dcinq,dzero)
cdmr --- Ajout du GDD5
         freezeI(i,j)=freezeI(i,j)+abs(min(ttemp,dzero))
         thawI(i,j)=thawI(i,j)+max(ttemp,dzero)
       enddo
      enddo

c- last iter. of the day :
      if (mod(istep,iatm).eq.0) then
        prcrit = abs(prcmin)
        do j=1,nlon
         do i=1,nlat
c- total precipitation :
           prcveg(i,j,1) = prcveg(i,j,1)+prcday(i,j)
c- precipitation above the daily threshold "prcmin" :
           if (prcday(i,j).ge.prcrit)
     &       prcveg(i,j,2)=prcveg(i,j,2)+prcday(i,j)
           prcday(i,j) = 0.
c- cumul Nb_days of dry soil conditions :
c~ #if ( ISOATM >= 2 )
c~            if (bmoisg(i,j,ieau).le.bmtdry)
c~      &        tpsdry(i,j) = tpsdry(i,j) + 1.
c~ #else
           if (bmoisg(i,j,iwater).le.bmtdry)
     &        tpsdry(i,j) = tpsdry(i,j) + 1.
c~ #endif
         enddo
        enddo
      endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c~ #if ( WRAP_EVOL == 2 )
c~ ! SDubinkina, needed for restart every N month, N<12
c~       if (mod(irunlabeld*iatm+istep,nstpyear*nyvegt).eq.0) then
c~ #else
      if (mod(istep,nstpyear*nyvegt).eq.0) then
c~ #endif
c-----------------------------------------------------------
c--Last time step of the year : compute Anual_Mean Climato :
        zmoy1 = 0.
        if (ncumvg.ge.1) zmoy1 = 1.d0 / DBLE(ncumvg)
        zmoy2 = zmoy1 * 360.
c- conversion : hauteur cumulee (m) -> precip. (mm/yr)
        zmoy3 = zmoy2 * DBLE(iatm*1000)
c- conversion Nb_day/yr
        zmoy4 = zmoy2 * DBLE(iatm)
c- precip for veget: Sum precip > prcmin (if prcmin > 0) + Reduction fct(tpsdry)
        kk2 = 2
        if (prcmin.lt.dzero) kk2 = 1
!dmr ### [DELETE] -- to avoid the creation of a fort.67 file ...         
!dmr ???        write(67,'(A,I2,1P4E14.6)') ' veget:k2,zmoy1,2,3,4=', kk2,
!dmr ???     &                               zmoy1, zmoy2, zmoy3, zmoy4
!dmr ### [DELETE]
cmdr --- For CTRL conditions MAXVAL freezeI is on the order of 86 000
        do j=1,nlon
         do i=1,nlat
           temveg(i,j) = temveg(i,j) * zmoy1
           gd0veg(i,j) = gd0veg(i,j) * zmoy2
cdmr --- Ajout du GDD5
           gd5veg(i,j) = gd5veg(i,j) * zmoy2
cdmr --- Ajout du GDD5
           freezeI(i,j) = freezeI(i,j) * zmoy2
           thawI(i,j) = thawI(i,j) * zmoy1
           prcveg(i,j,1) = prcveg(i,j,1) * zmoy3
           prcveg(i,j,2) = prcveg(i,j,2) * zmoy3
           tpsdry(i,j) = tpsdry(i,j) * zmoy4
         enddo
        enddo


#if( 1 )
      open(newunit=freezeIctl_id,file='outputdata/atmos/freezeI.ctl')
      write(freezeIctl_id,fmt="('dset   ^freezeI.dat')")
      write(freezeIctl_id,fmt="('options big_endian')")
      write(freezeIctl_id,fmt="('undef ',1p,e12.4)") -1.0e20
      write(freezeIctl_id,fmt="('title ECBILT orography')")
      write(freezeIctl_id,fmt="('xdef ',i3,' linear ',2f7.3)") 64,0.00,5.625
      write(freezeIctl_id,fmt="('ydef ',i3,' levels')") 32
      write(freezeIctl_id,
     >             fmt="(' -85.7606 -80.2688 -74.7445 -69.2130 -63.6786')")
      write(freezeIctl_id,
     >             fmt="(' -58.1430 -52.6065 -47.0696 -41.5325 -35.9951')")
      write(freezeIctl_id,
     >             fmt="(' -30.4576 -24.9199 -19.3822 -13.8445 -8.30670')")
      write(freezeIctl_id,fmt="(' -2.76890 2.76890 8.30670 13.8445 19.3822')")
      write(freezeIctl_id,fmt="('  24.9199 30.4576 35.9951 41.5325 47.0696')")
      write(freezeIctl_id,fmt="('  52.6065 58.1430 63.6786 69.2130 74.7445')")
      write(freezeIctl_id,fmt="('  80.2688 85.7606')")
      write(freezeIctl_id,fmt="('zdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
      write(freezeIctl_id,fmt="('tdef ',i4,' linear 1jan0001  1YR')") 1
      write(freezeIctl_id,fmt="('vars 1')")
      write(freezeIctl_id,fmt="('fi       1  99 Freeze Index ECBILT')")
      write(freezeIctl_id,fmt="('endvars')")
      close(freezeIctl_id)
      do i=1,nlon
        do j=1,nlat
          outdata(i,j)=freezeI(j,i)
        enddo
      enddo
      open(newunit=freezeIdat_id,CONVERT='BIG_ENDIAN'
     >    ,file='outputdata/atmos/freezeI.dat'
     >    ,form='unformatted'
     >    ,access='direct',recl=Size(outdata)*Kind(outdata(1,1)))
      write(freezeIdat_id,REC=1) outdata
      close(freezeIdat_id)

#endif

c- Reduce Precip by 0 -> rpfdry if Nb_Day_Dry(=tpsdry) > tmxdry -> tmxdry+dtrdry
        do j=1,nlon
         do i=1,nlat
           prcveg(i,j,2) = prcveg(i,j,kk2) *
     &      (one-rpfdry*min(one,max(zero,(tpsdry(i,j)-tmxdry)*unsdry)))
         enddo
        enddo

      else
c-----------------------------------------------------------
c-Not [last time_step of the year] => return directly (no call veget.routines).
c~ #if ( WRAP_EVOL == 2 )
c~ !         return
c~ ! SDubinkina, needed for restart every N month, N<12
c~ c-Not [last time_step of the year] => go to the saving files part
c~ c-Since it could be the last time_step of the calculation
c~         goto 25
c~ #else
        return
c~ #endif
      endif
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

c~ #if ( WRAP_EVOL == 2 )
c~ ! SDubinkina, needed for restart every N month, N<12
c~       veg_act=.false.
c~       if ((iyearvegout.eq.0).and.(mod(irunlabeld*iatm+istep,nstpyear*nyvegt).eq.0)) veg_act=.true.
c~       do 281 ieq=1,ieqveg

c~ !       if (kveget.lt.0 .or. iyear*kveget.ne.0)  then
c~ ! SDubinkina, needed for restart every N month, N<12
c~       if (kveget.lt.0 .or. iyearvegout*kveget.ne.0 .or. veg_act) then
c~ #else
      do ieq=1,ieqveg
      if (kveget.lt.0 .or. iyear*kveget.ne.0) then
c~ #endif
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3 ) Call the Vegetation Routines :
c-----------------------------------------------------------------------

c- co2 enrichment
c-AM      co2ghg = ghg(1)
          co2ghg = patmCO2

c...    SPATIAL LOOP

!      WRITE(*,*) "iter - veg=", ieq
      do j=1,nlon
       do i=1,nlat
#if( PF_CC > 0 )
!          fr_ndx   = freezeI(i,j)
!dmr --- Version de l'index invente par ckc, entre 0 et 1
          fr_ndx   = SQRT(freezeI(i,j))/(SQRT(freezeI(i,j))
     &    +SQRT(thawI(i,j)))
          if ( dataecb(i,j).LT.1000.) then
            pf_percent(i,j) = dataecb(i,j)/1000.
          else
            pf_percent(i,j) = 1.
          endif
#endif
        if (fracgr(i,j).gt.epss) then
c- Transfert var. for Veget. model (& modif Precip_2) :
          lat = i
          lon = j
          ave_t    = temveg(i,j)
          gdd0     = gd0veg(i,j)

          ave_pr   = prcveg(i,j,1)
          ave_pr05 = prcveg(i,j,2)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c... calculation of dynamics of carbon pools and forest ratio

          if (kveget.lt.0) then
            if (irad.eq.1) returnValue = ccstatR(fracgr,darea)
            returnValue = ccstat(fracgr,darea)
          else
            do k=1,kveget
              if (irad.eq.1) returnValue = ccdynR(fracgr,darea)
              returnValue = ccdyn(fracgr,darea)
            enddo
          endif
        endif
       enddo
      enddo
#if( PF_CC > 0 )
#if( 1 )
      open(newunit=permafctl_id,file='outputdata/atmos/permaf.ctl')
      write(permafctl_id,fmt="('dset   ^permaf.dat')")
      write(permafctl_id,fmt="('options big_endian')")
      write(permafctl_id,fmt="('undef ',1p,e12.4)") -1.0e20
      write(permafctl_id,fmt="('title ECBILT orography')")
      write(permafctl_id,fmt="('xdef ',i3,' linear ',2f7.3)") 64,0.00,5.625
      write(permafctl_id,fmt="('ydef ',i3,' levels')") 32
      write(permafctl_id,fmt="(' -85.7606 -80.2688 -74.7445 -69.2130 -63.6786')")
      write(permafctl_id,fmt="(' -58.1430 -52.6065 -47.0696 -41.5325 -35.9951')")
      write(permafctl_id,fmt="(' -30.4576 -24.9199 -19.3822 -13.8445 -8.30670')")
      write(permafctl_id,fmt="(' -2.76890 2.76890 8.30670 13.8445 19.3822')")
      write(permafctl_id,fmt="('  24.9199 30.4576 35.9951 41.5325 47.0696')")
      write(permafctl_id,fmt="('  52.6065 58.1430 63.6786 69.2130 74.7445')")
      write(permafctl_id,fmt="('  80.2688 85.7606')")
      write(permafctl_id,fmt="('zdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
      write(permafctl_id,fmt="('tdef ',i4,' linear 1jan0001  1YR')") 1
      write(permafctl_id,fmt="('vars 1')")
      write(permafctl_id,fmt="('pf       1  99 Permafrost Depth ECBILT')")
      write(permafctl_id,fmt="('endvars')")
      close(permafctl_id)
      do i=1,nlon
        do j=1,nlat
          outdata(i,j)=pf_percent(j,i)
        enddo
      enddo
      open(newunit=permafdat_id,CONVERT='BIG_ENDIAN'
     >    ,file='outputdata/atmos/permaf.dat'
     >    ,form='unformatted'
     >    ,access='direct',recl=Size(outdata)*Kind(outdata(1,1)))
      write(permafdat_id,REC=1) outdata
      close(permafdat_id)
#endif
#endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

      if (kveget.ne.0) then

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Compute albedo of land (except over Greenland & Antartica )
c      Compute bmoism en fonction de veget
c-----------------------------------------------------------------------

c- Antarctica Ice i: 1--> 5 [90.S -> 61.875 S]
c- Greenland Ice : continent & i >= 26 + |j-57| (inclus)
      nlatd2 = nlat / 2
      istart=6

C      bmoism(i,j)=bmoisg(i,j)
C      rs(i,j)=0.

      do i=1,nlat
         do j=1,nlon
C           bmoism(i,j)=bmoisg(i,j)
            bmoism(i,j)=0.15
C           rs(i,j)=0.
         enddo
      enddo

#if ISM == 1
chg allow tundra on Antarctica
c      if (flgisma) istart=6
      if (flgisma) istart=1
chg --
#endif

!dmr --- following section added by Huan Li to force the vegetation to LUH2 dataset
#if ( VEG_LUH == 1 )
c- pass fraction of tree, grass and desert from LUH2 dataset to vecode
      st(:,:)= st_luh(:,:)
      sg(:,:)= sg_luh(:,:)
      sd(:,:)= sd_luh(:,:)
#endif

!dmr --- Added to fix a underflow BUG
      WHERE(sd.LT.(10.d0*TINY(sd(1,1))))
        sd = 0.d0
      ENDWHERE
      WHERE(sg.LT.(10.d0*TINY(sg(1,1))))
        sg = 0.d0
      ENDWHERE
      WHERE(st.LT.(10.d0*TINY(st(1,1))))
        st = 0.d0
      ENDWHERE
      WHERE(sdr.LT.(10.d0*TINY(sdr(1,1))))
        sdr = 0.d0
      ENDWHERE
      WHERE(str.LT.(10.d0*TINY(str(1,1))))
        str = 0.d0
      ENDWHERE
      WHERE(sgr.LT.(10.d0*TINY(sgr(1,1))))
        sgr = 0.d0
      ENDWHERE
!dmr --- Added to fix a underflow BUG


      do j=1,nlon
       do i=istart,nlat
      igroen=26+abs(j-57)

c~ [DEPRECATED]
c~ #if ( ISM == 1 )
c~       if (flgismg) igroen=nlat+1
c~ #endif
c~ #if ( WRAP_EVOL == 2 )
c~ #if LGM == 1
c~         if ( fracgr(i,j).gt.epss) then
c~ #else
c~         if ( fracgr(i,j).gt.epss .and. i.lt.igroen ) then
c~ #endif
c~ #else
cdmr ###        if ( fracgr(i,j).gt.epss .and. i.lt.igroen ) then
cdmr --- Suppress ! We need to call it over the ice-sheet now !
        if ( fracgr(i,j).gt.epss ) then ! .and. i.lt.igroen ) then
c~ #endif
          ii =i / nlatd2
          if (ii.eq.2) ii=1
c- transition -> Tundra & Steppe (xx=1) : T_ann entre 4 et 0 deg.C
          xx = (4. - temveg(i,j)) * 0.25
          xx = min(1.,max(0.,xx))
c--Bright sand desert (Sahara) : + albedc :
          ddc = 0.
Csah1     if (  (i.ge.19.and.i.le.21) .and. (j.le.10.or.j.ge.62)
Csah1&     .or. (i.eq.22.and.j.le.7) ) ddc = 1.
          if ( (i.ge.19.and.i.le.23).and.(j.le.12.or.j.ge.62) ) ddc=1.
          if ( i.eq.23 .and. j.ge.5 .and. j.le.12 ) ddc = 0.5
c--test Albedo funct(LAI) :
C_2       zlai = blai(i,j,1)*st(i,j)+blai(i,j,2)*sg(i,j)
C_2       xxalb = veg_albd(i,zlai,snlt(i,j),temveg(i,j),ddc)
          do ns=1,4
            xxalb = st(i,j)*albet(ns)
     &            + sd(i,j)*( albed(ns) + ddc*albedc(ns) )
     &            + sg(i,j)*( albeg(ns) + xx*albegc(ns) )
c!Attention pour calcul FR, creer bmoismR, rsR et adapter ds lbm
           bmoism(i,j)=((st(i,j)*(1-snlt(i,j)))*0.25)
     &       + (sd(i,j)*0.1)+(sg(i,j)*0.15)+(st(i,j)*snlt(i,j)*0.25)
clbm2      bmoism(i,j)=((st(i,j)*(1-snlt(i,j)))*acwt*zrt/1020.)
clbm2&            + (sd(i,j)*acwd*zrd/1020.)
clbm2&            + (sg(i,j)*acwg*zrg/1020.)
clbm2&            + (st(i,j)*snlt(i,j)*acwn*zrn/1020.)
clbm1      rs(i,j)=((st(i,j)*(1-snlt(i,j)))*rst)+(sd(i,j)*rsd)
clbm1&            +(sg(i,j)*rsg)+(st(i,j)*snlt(i,j)*rsn)

c--Test Eq.Rain.Forest Low Albedo :
C_1       if (i.ge.15.and.i.le.18) xxalb = .09
          if (i.eq.16.or.i.eq.17) xxalb = .09
          if (i.eq.15.or.i.eq.18) xxalb = 0.5*xxalb + .045
C_2       do 430 ns=1,4
            albland(i,j,ns0inv(ns,ii)) = xxalb

c~ [DEPRECATED]
c~ #if ( ISM == 1 )
c~             if (flgism) then
c~              alblandism(i,j,ns0inv(ns,ii),tu) = xxalb
c~              alblandism(i,j,ns0inv(ns,ii),ca) = 0.55
c~ chg+ modified ice albedo
c~              if (i.lt.7) alblandism(i,j,ns0inv(ns,ii),ca) = albsnow(i)
c~              if (i.lt.7) alblandismR(i,j,ns0inv(ns,ii),ca) = albsnow(i)
c~ chg-
c~             endif
c~ #endif

!if ( IMSK == 1 )
!cdmr --- Added to get the icesheet properly set with reference to the icemask, even if not flgism
!           if (icemask(i,j).gt.0.9) then
!             st(i,j)=0.0 ! no trees on ice-sheets
!             sg(i,j)=0.0 ! no grass on ice-sheets
!             sd(i,j)=1.0 ! real desert on ice-sheets
!             snlt(i,j)=0.0 ! for consistency (not useful, but ...)
!           endif
!dmr --- Added to get the icesheet properly set with reference to the icemask, even if not flgism
!endif

          enddo
          forestfr(i,j) = st(i,j)
        endif
       enddo
      enddo

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Ouput on ASCII file.
c-----------------------------------------------------------------------

      if (ieqveg.ne.1) then
c~ #if ( WRAP_EVOL == 2 )
c~       if((mod(ieq,ieqvegwr).eq.0).or.((initialization.eqv..true.).and.(ieq.eq.1))) then
c~ #else
c      if((mod(ieq,ieqvegwr).eq.0).or.((iyear.eq.0).and.(ieq.eq.1))) then
      if((mod(ieq,ieqvegwr).eq.0)) then
c~ #endif
         flgwrt=.TRUE.
      else
         flgwrt=.FALSE.
      endif
      endif

c- Compute & write Global & Zonal mean Var. ; Write 2.D vegetation Var.
      if (flgwrt) call veget_wr(
     &             kveget, nyvegt, nwrveg,iyr0vg, iyear, nyears,
     &             temveg,gd0veg,prcveg,tpsdry,
     &             prcmin,bmtdry, epss,fracgr,darea,
cdmr --- Ajout du GDD5
     &             gd5veg,
cdmr --- Ajout du GDD5
     &             titveg)
      enddo
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  6 ) Prepare the next run (write restart) or the next year (reset).  |
c-----------------------------------------------------------------------

c~ #if ( WRAP_EVOL == 2 )
c~ ! SDubinkina, needed for restart every N month, N<12
c~       if (istep.eq.(ntotday*iatm)) then
c~ #else
!      if (iyear+nyvegt.gt.nyears) then
cdmr -- changed conditions to get regular restart as for atmosphere
      if ((iyear.ne.0).and.(mod(istep,nstpyear).eq.0)) then
        if (mod(iyear,nwrskip).eq.0.or.iyear.eq.nyears) then

c~ #endif
c--Last call : Write restart files veget.init & veget.rest
        open(newunit=vegetinit_id,file='veget.init',status='unknown'
     &             ,form='unformatted')
        write(vegetinit_id) temveg
        write(vegetinit_id) gd0veg
chli-- write gdd5 into file veget.init
        write(vegetinit_id) gd5veg
chli-- write gdd5 into file veget.init
        write(vegetinit_id) ((prcveg(i,j,1),i=1,nlat),j=1,nlon)
        write(vegetinit_id) ((prcveg(i,j,2),i=1,nlat),j=1,nlon)
        write(vegetinit_id) tpsdry
        write(vegetinit_id) iyear+iyr0vg, bmtdry, titveg
        close(vegetinit_id)

        if (kveget.ne.0) then
        open(newunit=vegetrest_id,file='veget.rest',status='unknown'
     &             ,form='unformatted')

        write(vegetrest_id) iyear+iyr0vg, bmtdry, titveg
        write(vegetrest_id) st
        write(vegetrest_id) sg
        write(vegetrest_id) sd
        write(vegetrest_id) snlt
        write(vegetrest_id) temveg
        write(vegetrest_id) b1t
        write(vegetrest_id) b1g
        write(vegetrest_id) b2t
        write(vegetrest_id) b2g
        write(vegetrest_id) b3t
        write(vegetrest_id) b3g
        write(vegetrest_id) b4t
        write(vegetrest_id) b4g
        write(vegetrest_id) b1
        write(vegetrest_id) b2
        write(vegetrest_id) b3
        write(vegetrest_id) b4
        write(vegetrest_id) anup
#if ( CYCC == 2 )
        write(vegetrest_id) B1T13
        write(vegetrest_id) B1G13
        write(vegetrest_id) B2T13
        write(vegetrest_id) B2G13
        write(vegetrest_id) B3T13
        write(vegetrest_id) B3G13
        write(vegetrest_id) B4T13
        write(vegetrest_id) B4G13
        write(vegetrest_id) B1T14
        write(vegetrest_id) B1G14
        write(vegetrest_id) B2T14
        write(vegetrest_id) B2G14
        write(vegetrest_id) B3T14
        write(vegetrest_id) B3G14
        write(vegetrest_id) B4T14
        write(vegetrest_id) B4G14
#if ( KC14 == 1 )
        write(vegetrest_id) cav_la14_b
        write(vegetrest_id) cav_la_b
#endif
#endif
        if (iscendef.eq.1) then
          write(vegetrest_id) st_const
          write(vegetrest_id) str
          write(vegetrest_id) sgr
          write(vegetrest_id) sdr
        else
          write(vegetrest_id) st
        endif

        close(vegetrest_id)
!dmr ### [DELETE]
!        close(iveg+7)   !dmr -> this seems to relate to outputdata/vegetation/veget.zav
!                       !dmr    i see no good reason to have that here ...
!dmr ### [DELETE]

!dmr ### [DELETE] -- to avoid the creation of a fort.67 file ...        
!dmr ???        write(67,'(2A,I8,A,F8.3)') ' Last call: write "veget.init"',
!dmr ???     &      ' & "veget.rest" ; yr=', iyear+iyr0vg, ' ; CO_2=', co2ghg
!dmr ### [DELETE]    
      else
!dmr ### [DELETE] -- to avoid the creation of a fort.67 file ...               
!dmr ???        write(67,'(2A,I8,A,F8.3)') ' Last call: write "veget.init"',
!dmr ???     &                     ' ; yr=', iyear+iyr0vg, ' ; CO_2=', co2ghg
!dmr ### [DELETE]    
      endif
c-----
      endif
c~ #if ( WRAP_EVOL == 0 || WRAP_EVOL == 1 )
      endif
c~ #endif
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Before starting a new year, reset Cum.Var (Climato) at Zero :
      ncumvg = 0
      do j=1,nlon
       do i=1,nlat
         temveg(i,j) = 0.0
         gd0veg(i,j) = 0.0
cdmr --- Ajout du GDD5
         gd5veg(i,j) = 0.0
cdmr --- Ajout du GDD5
         prcveg(i,j,1) = 0.0
         prcveg(i,j,2) = 0.0
         tpsdry(i,j) = 0.0
         freezeI(i,j)= 0.0
         thawI(i,j)= 0.0
       enddo
      enddo
      kveget = abs(kveget)
#if ( WRAP_EVOL == 2 )
! SDubinkina, needed for restart every N month, N<12
c--read Cum.Var (Climato) from a file :

      if(initialization.eqv..true.) then
        open(newunit=veget_accurest_id,file='veget_accu.rest'
     >      ,status='old',form='unformatted',iostat=ios)
      if (ios.eq.0) then
          read(veget_accurest_id) ncumvg
          read(veget_accurest_id) temveg
          read(veget_accurest_id) prcday
          read(veget_accurest_id) gd0veg
chli---read the variable gd5veg
          read(veget_accurest_id) gd5veg
          read(veget_accurest_id) ((prcveg(i,j,1),i=1,nlat),j=1,nlon)
          read(veget_accurest_id) ((prcveg(i,j,2),i=1,nlat),j=1,nlon)
          read(veget_accurest_id) tpsdry
      endif
      close(veget_accurest_id)

      endif
! SDubinkina, needed for restart every N month, N<12
c-Not [last time_step of the year] and [last time_step of the year]
   25 continue
      if ((istep.eq.(ntotday*iatm)).and.(initialization.eqv..false.)) then
c--Last call : Write restart file veget_accu.rest

        open(newunit=veget_accurest_id,file='veget_accu.rest'
     >      ,status='unknown',form='unformatted')
        write(veget_accurest_id) ncumvg
        write(veget_accurest_id) temveg
        write(veget_accurest_id) prcday
        write(veget_accurest_id) gd0veg
chli---write gd5veg into file veget_accu.test
        write(veget_accurest_id) gd5veg
        write(veget_accurest_id) ((prcveg(i,j,1),i=1,nlat),j=1,nlon)
        write(veget_accurest_id) ((prcveg(i,j,2),i=1,nlat),j=1,nlon)
        write(veget_accurest_id) tpsdry
        write(veget_accurest_id) iyearvegout+iyr0vg, bmtdry, titveg
        close(veget_accurest_id)

      else
c-Not [last time_step of the run] => return without saving veget_accu.rest
	return
      endif
#endif

      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine veget -
      end subroutine veget
