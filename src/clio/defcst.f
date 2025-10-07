!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009

      SUBROUTINE defcst(nn99)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Initialisation of the various physical variables and parameters
!  nn99, : 1 => open "mouchard" and nn99 <- 2
!   0 => if kfond=-1,-3 , open "mouchard" and nn99 <- 2
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  modif : 23/08/02

!! START_OF_USE_SECTION

      use const_mod, only: cpo, cstmax, cstmin, epsil, gpes, omega, one,
     &    pi, radian, rho0, rterre, svrdrp, unsrt, untour, yeaday, zero



      use para0_mod, only: ltest, imax, jmax, kmax, nsmax, jsepar
      use para_mod,  only: nvmax

      use bloc0_mod, only: ahrap, ajcmix, alphxu, alphxv, alphyu, alphyv
     &  , avsn2, bering, bheat, ccfmn, ccfmx, cdbot, ibera, iberp,
     &    ihyster, jbera, jberp, nitrun, nsav, qavs, txifcb,txifcu,
     &    vcor, xslop, zfluxm, zfluxms
     &  , avk0, avkb, avnub, slopmgm, rifsmx, avnu0, aitd, slopemax, ai
     &  , ahs, rifumx, ahu, ahh, ahe, afilt, avv, dts, xfreshsv, dtu
     &  , dtb, scssv, scpme, zfluxms, zfluxm, xslop, vcor, txifcu
     &  , txifcb, scal0, alphah, alphmi, alphaz, txeflx, alphgr, algrmn
     &  , txiads, txiadu, txidfs, txidfu, bering_prev


      use bloc_mod, only: ghmax, ghmin, icoupl, icoutp, idyn, itau_slow,
     &    kajul, kavsmx, kfond, kforc, kinput, koutpu, kstart, lstab,
     &    mdforc, nclin, nclmoy, ninfo, nitrap, nsewfr, nsplaj, nsplit,
     &    ntmoy, numit, nwa, nwjl, nwm, nwtal, nwtest, nwtoa
     &  , nwtom, q2tmin, sqrghm, varfor, vkappa, vlmin, zlotur

      use datadc_mod, only: nvrajc, nvral, nvras, nvrau, nvret, nvrfq,
     &    nvrfw, nvrhg, nvrhn, nvrmom, nvrqs, nvrtbq, nvrtgx, nvrtgy,
     &    nvrtke, nvrts, nvru, nvrub, nvrug, nvrum, nvrv, nvrvb, nvrvg,
     &    nvrvm, nvrxzo

      use ice_mod, only: alphs, amax, beta, cevap, cnscg, ddtb, emig,
     &    exld, hakdif, hakspl, hgdif, hglim, hgmin, hibspl, hmelt,
     &    hndif, hnzst, hth, nbits, parlat, parsub, rcpg, rcpn, rhoesn,
     &    rhog, rhon, sglace, stefan, swiqst, tfsg, tfsn, too, uscomi,
     &    vkarmn, xkg, xkn, xlg, xln, xsn, zemise, acrit, hgcrit

      use dynami_mod, only: alpha, angvg, bound, c, cangvg, creepl, cw,
     &    dm, gridsz, iameth, nbitdf, nbitdr, nbiter, nlmaxn, nlmaxs,
     &    nlminn, nlmins, om, pstarh, ren, resl, rhoco, rhoco2, sangvg,
     &    sber, usdt, usecc2, uvdif, zepsd1, zepsd2, zetamn

      use reper_mod,  only: icheck, jcheck, kcheck, unstyr, unitfx, rapp0, rapp1
     &                   , yforc, filcor, scalwr
      use varno_mod,  only: nvrl, ltyp

#if ( ISOOCN >= 1 )
      use para0_mod, only: owatert, isoocn_restart
#endif

      use global_constants_mod, only: dblp=>dp, ip


#if ( BATHY >= 2 )
      use update_clio_bathy_tools, only: bering_nb, bering_date,
     >       bering_value, la_date,
     >       kamax_nb, kamax_date, kamax_value
      use comcoup_mod, only: kamax
#endif


      use ipcc_output_mod, only: tmc0

      use newunit_clio_mod, only: clio3_out_id, mouchard_id

!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

! [IMPLCTNONE] #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "ice.com"
! [SCRPTCOM] #include "dynami.com"
! [SCRPTCOM] #include "reper.com"
! [SCRPTCOM] #include "varno.com"

!! END_OF_INCLUDE_SECTION

!--local variables :
      real(dblp), dimension(kmax) :: coef2, coef3, coef4, coef5, coef6,
     &                               coef7, coef8, coef9, coef10,coef11
!     dimension zrr(5),zd1(5),zd2(5)
!ic0  dimension acrit(2), hgcrit(2)


      integer(ip) :: i, jj, k, n, nn, nninfo, ns, nv, nn99
      real(dblp)  :: ecc, pstar

      character*70 line
#if ( ISOOCN >= 1 )
      real  tmpzfluxm,tmpzfluxms(owatert)
#else
      real  tmpzfluxm,tmpzfluxms
#endif
      integer ios
!--instructions "data" :
! [SCRPTCOM] #include "datadc.com"

!~       real*8 moc,tmc,tmc0,tsurfmean,cland,thex
!~       common/IPCC_out2/moc,tmc,tmc0,tsurfmean,cland,thex

#if ( BATHY >= 2 )
       real distance_before, distance_after
#endif

      integer(ip) :: run_param_id, fresh_for_dat_id, thermo_param_id,
     &               dynami_param_id, correcw_dat_id
      integer(kind=ip)  :: lok_nsmax = 13
        
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Reading of the parameter of the run ; arrays of associated parameters |
!-----------------------------------------------------------------------

       open(newunit=run_param_id,file='run.param',status='old')
!- header 1 :
        read(run_param_id,*)
        read(run_param_id,*)
!- Checking at one particular point :
        icheck = 0
        jcheck = 0
        kcheck = 0
        read(run_param_id,'(A)') line
        if (line(14:14).ne.' ') read(line(15:),*) icheck, jcheck, kcheck
!ray    read(run_param_id,'(15x,3i5)') icheck, jcheck, kcheck
!- header 2 :
      do 100 n=1,4
        read(run_param_id,*)
 100  continue
!- option of the code :
!---
! kfond > 0 => Nb_niv.=kfond (flat); < -2 => DownSloping , -1,-3 => mouchard
! lstab > 1 => AjConv by pairs, lstab times ; -1 total mixing; 0 : nothing
!        -3 =>  permutation & mixing of a fraction ajcmix
      read(run_param_id,*)
      read(run_param_id,*) kfond, lstab, nsewfr, bering, ajcmix, xslop
!- initial value (by level) for each scalar :
      print*, "in defcst, read nsmax = ", nsmax, lok_nsmax
      if (lok_nsmax.GT.nsmax) then
         lok_nsmax = nsmax
      endif
      do 105 n=1,lok_nsmax
        read(run_param_id,*)
        read(run_param_id,*) (scal0(k,n),k=kmax,1,-1)
 105  continue
!- conversion Celcius -> Kelvin :
      do 110 k=1,kmax
        scal0(k,1) = scal0(k,1) + 273.15d0
 110  continue
!- restoring :
      nitrap = 0
      ahrap = 0.
      read(run_param_id,*)
      read(run_param_id,*) unstyr, (rapp0(k),k=0,(lok_nsmax+2)), ahrap, nitrap
      read(run_param_id,*)
      read(run_param_id,*) (rapp1(k),k=kmax,1,-1)
!- time steps :
      read(run_param_id,*)
      read(run_param_id,*) dts(kmax), dtu, dtb
      read(run_param_id,*)
      read(run_param_id,*) (coef2(k),k=kmax,1,-1)
!- Parameters associated with the coupling with an atmospheric model
! icoupl=0 : not coupled
! icoupl=1 : following of a preceding run
! icoupl=2 : set time to zero
! icoupl=3 : Beginning of a coupled run with starting state
! Coupling type
! icoutp=1 : LMD 5bis
! icoutp=2 : LMD Z
! icoutp=3 : ECBILT
! icoutp=4 : (to be continued)

      read(run_param_id,*)
      read(run_param_id,*) icoupl, icoutp, itau_slow
!-----
!- nb-iterations :
      read(run_param_id,*)
      read(run_param_id,*) nitrun, nsplit, nsplaj
      read(run_param_id,*)
      read(run_param_id,*) nclin, nclmoy
!- option start/read/write :
      read(run_param_id,*)
      read(run_param_id,*) kstart, kinput, koutpu
!- frequ. writing :
      read(run_param_id,*)
      read(run_param_id,*) nsav, ninfo, ntmoy !dmr ### ,ntout
!- option for the sea-ice ...
!       nwtest=0, WRITE OUTPUTS OF YEARS:
!                              nw(..), 2*nw(..), ...
!       nwtest=1, WRITE OUTPUTS OF YEARS:
!                              1, nw(..)+1, 2*nw(..)+1, ...
      read(run_param_id,*)
      read(run_param_id,*) nwjl, nwm, nwa, nwtoa, nwtom
      read(run_param_id,*)
      read(run_param_id,*) nwtest, nwtal

!- Bottom drag coefficient
      read(run_param_id,*)
!     read(run_param_id,*) cdbot
      read(run_param_id,*) cdbot, txifcb, txifcu

!- diffusivity & visc. H. :
      read(run_param_id,*)
      read(run_param_id,*) ahs(kmax), ahu, ahe
      read(run_param_id,*)
      read(run_param_id,*) (coef3(k),k=kmax,1,-1)
!- diffusivity & visc. V. :
      read(run_param_id,*)
      read(run_param_id,*) avnub(kmax), avnu0(kmax), rifumx
      read(run_param_id,*)
      read(run_param_id,*) avkb(kmax),  avk0(kmax),  rifsmx
      read(run_param_id,*)
      read(run_param_id,*) (coef4(k),k=kmax,1,-1)
      read(run_param_id,*)
      read(run_param_id,*) (coef5(k),k=kmax,1,-1)
      read(run_param_id,*)
      read(run_param_id,*) (coef6(k),k=kmax,1,-1)
      read(run_param_id,*)
      read(run_param_id,*) (coef7(k),k=kmax,1,-1)
!- upwind rate :
      read(run_param_id,*)
      read(run_param_id,*) alphxu, alphxv, alphyu, alphyv
      read(run_param_id,*)
      read(run_param_id,*) alphah(1), alphah(2),
     &           (alphgr(ns),ns=1,lok_nsmax)
     &         , (algrmn(ns),ns=1,lok_nsmax)
      read(run_param_id,*)
      read(run_param_id,*) (alphmi(k),k=kmax,1,-1)
      read(run_param_id,*)
      read(run_param_id,*) (alphaz(k),k=kmax,1,-1)
!- reading of the parameters for the formulation of Avs(N2)
      read(run_param_id,*)
      read(run_param_id,*) kavsmx, qavs, avsn2, ccfmn, ccfmx
!- rate Implic. / Explic.
      read(run_param_id,*)
      read(run_param_id,*) txiadu, txiads, txidfu, txidfs, txeflx(1:lok_nsmax)
!- modif forcing :
      read(run_param_id,*)
      read(run_param_id,'(A)') line
!ray  read(run_param_id,'(I4,6X,A40)') mdforc,filcor
!- which forcing files read :
! kforc : (unite-> constant, dizaine-> saisonnier, centaine ...)
      read(run_param_id,*)
      read(run_param_id,*) kforc, (yforc(ns),ns=0,lok_nsmax)

!- value of precipitations and value at surface for the scalars
      read(run_param_id,*)
      read(run_param_id,*) (scpme(ns),ns=1,lok_nsmax)
      read(run_param_id,*)
      read(run_param_id,*) (scssv(ns),ns=1,lok_nsmax)

!- constants and parameters for the turbulence model.
      read(run_param_id,*)
      read(run_param_id,*) q2tmin,zlotur,vlmin,varfor,kajul
      vkappa  = 0.4d0
      ghmax  = 0.0233d0
      ghmin  = -0.28d0
      sqrghm = sqrt(-ghmin)

!- volume correction + bottom heat flux + hysteresis or not
      read(run_param_id,*)
      read(run_param_id,*) vcor,bheat,ihyster

!- Parameters for isoslope.f, isodiffu.f, isoadvec.f
      read(run_param_id,*)
      read(run_param_id,*) ai(kmax),slopemax(kmax)
      read(run_param_id,*)
      read(run_param_id,*) (coef8(k),k=kmax,1,-1)
      read(run_param_id,*)
      read(run_param_id,*) (coef9(k),k=kmax,1,-1)
      read(run_param_id,*)
      read(run_param_id,*) aitd(kmax),slopmgm(kmax),afilt,ahh,avv
      read(run_param_id,*)
      read(run_param_id,*) (coef10(k),k=kmax,1,-1)
      read(run_param_id,*)
      read(run_param_id,*) (coef11(k),k=kmax,1,-1)

      close(run_param_id)

#if ( BATHY >= 2 )
!nb read Bering scenario
      call read_bering
      bering=bering_value(1)
      do i=2,bering_nb
           if ((la_date.le.bering_date(i-1))
     >        .and.(la_date.ge.bering_date(i))) then
              distance_before=bering_date(i-1)-la_date
              distance_after=la_date-bering_date(i)
              if (distance_before.le.distance_after) then
                  bering=bering_value(i-1)
              else
                  bering=bering_value(i)
              endif
           elseif (la_date.le.bering_date(bering_nb)) then
              bering=bering_value(bering_nb)
           endif
       enddo
       write(*,*) 'Bering value at year ', la_date, '  updated to ', bering


!nb read kamax scenario
      call read_kamax
      kamax=kamax_value(1)
      !print*, 'la_date', la_date
      do i=2,kamax_nb
           !print*, 'kamax_date',kamax_date(i-1),kamax_date(i)
           if ((la_date.le.kamax_date(i-1))
     >        .and.(la_date.ge.kamax_date(i))) then
              distance_before=kamax_date(i-1)-la_date
              distance_after=la_date-kamax_date(i)
              if (distance_before.le.distance_after) then
                  kamax=kamax_value(i-1)
              else
                  kamax=kamax_value(i)
              endif
            elseif(la_date.le.kamax_date(kamax_nb)) then
              kamax=kamax_value(kamax_nb)
           endif
       enddo
       write(*,*) 'kamax value at year ', la_date, '  updated to ',kamax
#endif

cnb keep bering value in memory
       bering_prev=bering


      xfreshsv(:)=0.0
      if (ihyster.ne.0) then
        open(newunit=fresh_for_dat_id,file='inputdata/fresh_for.dat')
         do i=1,5000
           read(fresh_for_dat_id,*) jj,xfreshsv(i)
         enddo
        close(fresh_for_dat_id)
      endif



      icheck = min(imax,icheck)
      jcheck = min(jmax,jcheck)
      kcheck = min(kmax,kcheck)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--1.2 computation of the other parameters :
!--------------

      do 120 k=1,kmax
        dts(k) = coef2(k) * dts(kmax)
        ahs(k) = coef3(k) * ahs(kmax)
        avnub(k) = coef4(k) * avnub(kmax)
        avnu0(k) = coef5(k) * avnu0(kmax)
        avkb(k)  = coef6(k) * avkb(kmax)
        avk0(k)  = coef7(k) * avk0(kmax)
        ai(k)=coef8(k)*ai(kmax)
        slopemax(k)=coef9(k)*slopemax(kmax)
        aitd(k)=coef10(k)*aitd(kmax)
        slopmgm(k)=coef11(k)*slopmgm(kmax)
 120  continue
!- unchanged values unchanged at the bottom :
        avnub(1) = coef4(1)
        avnu0(1) = coef5(1)
        avkb(1)  = coef6(1)
        avk0(1)  = coef7(1)
!-
      if(nsplit.le.nsplaj.or.nsplaj.lt.0) then
        write(clio3_out_id,*) 'STOP in "defcst", nsplit & nsplaj =',
     &             nsplit, nsplaj
        stop
      endif
      if(nsewfr.eq.0) nsewfr = 1
      nn = 500 * nsav
      if (nn.ge.1.and.nn.le.numit) then
        write(clio3_out_id,*) 'STOP because Pb for parameter "nsav" :'
     &           //' too frequent output !'
        stop
      endif
      if ( koutpu.ge.2 ) then
!- test for the computation of the averages of Ajust.Conv. :
        if (ninfo.eq.0) then
          write(clio3_out_id,*) 'Defcst : modif of the parameter "ninfo" :'
          write(clio3_out_id,*) 'Old , New :', ninfo, nsav
          ninfo = nsav
        elseif ( mod(nsav,ninfo).ne.0 ) then
          nninfo = abs(ninfo)
          do 130 nn=nninfo,nsav
            if (mod(nsav,nn).eq.0) goto 135
 130      continue
          nn = nsav
 135      continue
          nn = sign(nn,ninfo)
          write(clio3_out_id,*) 'Defcst : modif of the parameter "ninfo" :'
          write(clio3_out_id,*) 'Old , New :', ninfo, nn
          ninfo = nn
        endif
      endif

!--Bering Strait opening:
      ibera = 0
      iberp = 0
      jbera = jmax
      jberp = jsepar
      if ( ltest.eq.3 .and. bering.ne.0. ) iberp = 1
      !write(*,*) 'IBERP 10 ', iberp, ltest, bering

!--Modificarion of the forcaging if needed :
      read(line,*) mdforc
      if (mdforc.ne.0) read(line,'(10x,A40)') filcor

!--Opening of the file. "mouchard" (=unit 99) and modification of "nn99" for the outputs
      if (nn99.eq.0 .and. mod(kfond,2).eq.-1 ) nn99 = 1
      if(nn99.eq.1) then
        nn99 = 2
        open(newunit=mouchard_id,
     &       file='outputdata/globals/mouchard',status='unknown')
      endif

!--Initialisation of the arrsy "nvrl" (depending on  "koutpu"),
!   nvrl(nv)=1 => write the variable "nv" in the ouput/restart file.
      nvrl(nvret ) = 1
      nvrl(nvrub ) = 1
      nvrl(nvrvb ) = 1
      nvrl(nvru  ) = 1
      nvrl(nvrv  ) = 1
      nvrl(nvrtke) = 1
      do 170 ns=1,nsmax
        nvrl(ns  ) = 1
        if (koutpu.ge.2) nvrl(ns+nvrfw) = 1
 170  continue
      if (koutpu.ge.2) then
        nvrl(nvrfw ) = 1
        nvrl(nvrajc) = 1
!       nvrl(nvrb ) = 1
!       nvrl(nvrn2 ) = 1
        nvrl(nvras ) = 1
        nvrl(nvrau ) = 1
      endif
      nvrl(nvrum)  = 1
      nvrl(nvrvm)  = 1
      nvrl(nvrhg)  = 1
      nvrl(nvrfq)  = 1
      nvrl(nvrqs)  = 1
      nvrl(nvrhn)  = 1
      nvrl(nvral)  = 1
      nvrl(nvrts)  = 1
      nvrl(nvrug)  = 1
      nvrl(nvrvg)  = 1
      nvrl(nvrtbq) = 1
      nvrl(nvrxzo) = 1
      nvrl(nvrtgx) = 1
      nvrl(nvrtgy) = 1
      nvrl(nvrmom) = 1
      do 180 nv=1,nvmax
        if (ltyp(nv).eq.99) nvrl(nv) = 0
 180  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Initiailisation of   the numerical/physical constants  |
!-----------------------------------------------------------------------

      zero = 0.d0
      one = 1.d0
      epsil = 0.1d-9
      cstmin = -1.d20
      cstmax =  1.d20
      untour = 360.d0

      pi = 4.0 * atan(one)
      radian = pi / 180.0
      omega = 2.0 * pi / 86164.0
      rterre = 6371000.0
      unsrt = 1.0 / rterre
      gpes = 9.80d0
      rho0 = 1030.0
      svrdrp = 1.0d-6
      cpo = 4002.
!--number of day in a year :
!ic0  yeaday = 365.25d0
      yeaday = 365.0
      yeaday = 360.0

!--used in redforc, informe : change unit(model)  -> standard  unit
      unitfx(0) = 1.
      do 290 ns=1,nsmax
        scalwr(ns) = 0.
        unitfx(ns) = 1.
 290  continue
!- unitfx(0) = Year(s) ;
!- unitfx(1) = rho.Cp : Fx. en W/m2 ; unitfx(2) = Year / 34.7 g/l : Fx. en m/y
      unitfx(0) = yeaday * 86400.0
      unitfx(1) = -4.1d+06
      unitfx(2) =  9.1d+05

!- initialisation of the bottom stress  <- moved to "defgrid"

      if (nn99.eq.2) write(mouchard_id,*) 'year : yeaday =', yeaday

!ic0  return

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Values linked with the sea-ice thermodynamics   |
!-----------------------------------------------------------------------

!--3.1 Physical constants.
!--------------------------

      tfsn=273.15d0
      tfsg=273.05d0
      xkn=0.22d0
      xkg=2.034396d0
      rcpn=6.9069e+05
      rcpg=1.8837e+06
      xlg=300.33e+06
      xln=110.121e+06
      xsn=2.8e+06
      rhog=900.0
      rhon=330.0
!
!     emig=0.97d0
!     emig=1.0d0
!     sglace=4.
      sglace=6.
!
      rhoesn=1/rhon
!
!--3.2. Physical constants for the heat fluxes.
!-----------------------------------------------
!
      stefan=5.67e-08
!
      too=273.16d0
      vkarmn=0.4d0
      cevap=2.5e+06
      zemise=0.97d0
      zemise=1.00d0
!
!--3.3.Physical parameters.
!--------------------------
!
      open(newunit=thermo_param_id,
     &     file='inputdata/clio/thermo.param', status='old')
      read(thermo_param_id,*)
      read(thermo_param_id,*)
      read(thermo_param_id,*)
      read(thermo_param_id,*)
      read(thermo_param_id,*)
      read(thermo_param_id,*) hmelt
      read(thermo_param_id,*)
      read(thermo_param_id,*) acrit(1)
      read(thermo_param_id,*)
      read(thermo_param_id,*) acrit(2)
      read(thermo_param_id,*)
      read(thermo_param_id,*) hgcrit(1)
      read(thermo_param_id,*)
      read(thermo_param_id,*) hgcrit(2)
      read(thermo_param_id,*)
      read(thermo_param_id,*) emig
!
!--3.4. Numerical parameters.
!----------------------------
!
      read(thermo_param_id,*)
      read(thermo_param_id,*) hgmin
      read(thermo_param_id,*)
      read(thermo_param_id,*) hndif
      read(thermo_param_id,*)
      read(thermo_param_id,*) hgdif
      read(thermo_param_id,*)
      read(thermo_param_id,*) hglim
      read(thermo_param_id,*)
      read(thermo_param_id,*) amax
      read(thermo_param_id,*)
      read(thermo_param_id,*) swiqst
!
      uscomi=1.0/(1.0-amax)
!
!--3.5. Numerical caracteristics.
!---------------------------------
!
      read(thermo_param_id,*)
      read(thermo_param_id,*)
      read(thermo_param_id,*) beta
      read(thermo_param_id,*)
      read(thermo_param_id,*) ddtb
      read(thermo_param_id,*)
      read(thermo_param_id,*) nbits
      read(thermo_param_id,*)
      read(thermo_param_id,*) parlat
      read(thermo_param_id,*)
      read(thermo_param_id,*) hakspl
      read(thermo_param_id,*)
      read(thermo_param_id,*) hibspl
      read(thermo_param_id,*)
      read(thermo_param_id,*) exld
      read(thermo_param_id,*)
      read(thermo_param_id,*) hakdif
      read(thermo_param_id,*)
      read(thermo_param_id,*) hth
      read(thermo_param_id,*)
      read(thermo_param_id,*) hnzst
      read(thermo_param_id,*)
      read(thermo_param_id,*) parsub
      read(thermo_param_id,*)
      read(thermo_param_id,*) alphs
!
      xkn=hakdif*xkn
      xkg=hakdif*xkg
      if ((hndif.gt.100.0).or.(hgdif.gt.100.0)) then
        cnscg = 0.0
      else
        cnscg = rcpn/rcpg
      endif
!
      close(thermo_param_id)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Values for ice-dynamic model.                                   |
!-----------------------------------------------------------------------
!
!--3.1. begin reading file of parameters.
!-----------------------------------------
!
      open(newunit=dynami_param_id,
     &      file='inputdata/clio/dynami.param',status='old')
!
      read(dynami_param_id,*)
      read(dynami_param_id,*)
      read(dynami_param_id,*)
      read(dynami_param_id,*)
      read(dynami_param_id,*)
!
!--3.2. numerical characteristics.
!---------------------------------
!
      read(dynami_param_id,*)
      read(dynami_param_id,*) idyn
      read(dynami_param_id,*)
      read(dynami_param_id,*) zepsd1
      read(dynami_param_id,*)
      read(dynami_param_id,*) zepsd2
      read(dynami_param_id,*)
      read(dynami_param_id,*) nlminn
      read(dynami_param_id,*)
      read(dynami_param_id,*) nlmaxs
      nlmaxn=jmax-1
      nlmins=3
      read(dynami_param_id,*)
      read(dynami_param_id,*) usdt
      read(dynami_param_id,*)
      read(dynami_param_id,*) alpha
      read(dynami_param_id,*)
      read(dynami_param_id,*) bound
      read(dynami_param_id,*)
      read(dynami_param_id,*) dm
      read(dynami_param_id,*)
      read(dynami_param_id,*) nbitdf
      read(dynami_param_id,*)
      read(dynami_param_id,*) nbiter
      read(dynami_param_id,*)
      read(dynami_param_id,*) nbitdr
      read(dynami_param_id,*)
      read(dynami_param_id,*) om
      read(dynami_param_id,*)
      read(dynami_param_id,*) resl
!
!--3.3. physical parameters.
!---------------------------
!
      read(dynami_param_id,*)
      read(dynami_param_id,*) cw
      read(dynami_param_id,*)
      read(dynami_param_id,*) angvg
      read(dynami_param_id,*)
      read(dynami_param_id,*) pstar
      read(dynami_param_id,*)
      read(dynami_param_id,*) c
      read(dynami_param_id,*)
      read(dynami_param_id,*) zetamn
      read(dynami_param_id,*)
      read(dynami_param_id,*) creepl
      read(dynami_param_id,*)
      read(dynami_param_id,*) ecc
      read(dynami_param_id,*)
      read(dynami_param_id,*) uvdif
      read(dynami_param_id,*)
      read(dynami_param_id,*) ren
      read(dynami_param_id,*)
      read(dynami_param_id,*) gridsz
      read(dynami_param_id,*)
      read(dynami_param_id,*) iameth
!
      close(dynami_param_id)
!
      usecc2 = 1.0/(ecc*ecc)
      rhoco  = rho0*cw
      rhoco2 = rhoco*rhoco
      angvg  = angvg*radian
      sangvg = sin(angvg)
      cangvg = cos(angvg)
      pstarh = pstar/2.0
      sber   = 1.-max(zero,sign(one,0.1-bering))
!
!- 3.4 Mean value of w

      zfluxm=0.0
      zfluxms=0.0
      if (vcor.gt.zero) then
        open(newunit=correcw_dat_id,
     &       file='startdata/correcw.dat',form='formatted')
        ios=0
        do while (ios==0)
#if ( ISOOCN >= 1 )
        if (isoocn_restart.NE.0) then
!dmr --- Pour les isotopes, j'ai besoin d'un format specifiant
         read(correcw_dat_id,'(6E25.18)',end=974,iostat=ios) tmpzfluxm,tmpzfluxms
         zfluxm=tmpzfluxm
         zfluxms=tmpzfluxms ! 2 is salt there ...
        else
         read(correcw_dat_id,*,end=974,iostat=ios) tmpzfluxm,tmpzfluxms(2)
         zfluxm=tmpzfluxm
         zfluxms(2)=tmpzfluxms(2) ! 2 is salt there ...
        endif
#else
! S'il n'y a qu'une valeur sur la ligne => tmc & ios>0 => va ligne 974
        read(correcw_dat_id,*,end=974,iostat=ios) tmpzfluxm,tmpzfluxms
        zfluxm=tmpzfluxm
        zfluxms=tmpzfluxms
#endif
       enddo
974    if (ios.NE.0)tmc0=tmpzfluxm
       close(correcw_dat_id)
       write(clio3_out_id,*) 'correction vol applied',vcor,zfluxm,zfluxms
      endif

!--3.4. Geostrophic velocities
!------------------------------
!
!     open (89,file='/u14/grpastr/hgs/data/geost/ugeost.mig',
!    &      form='unformatted',status='old')
!     read(89) uost
!     close(89)
!     open (89,file='/u14/grpastr/hgs/data/geost/vgeost.mig',
!    &      form='unformatted',status='old')
!     read(89) vost
!     close(89)

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- end of the routine defcst -
      end
