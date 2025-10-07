!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009

      subroutine init_CLIO(irunlabelclio)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  modif : 23/08/02


!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

! [DEPRECATED]      use comemic_mod, only: flgicb
      use const_mod

#if ( PATH >= 1 )
      use path_mod, only: path_init, scalstart_path, scalend_path
#endif

#if ( NEOD >= 1 )
      use neodymium_mod, only: neodymium_init, 
     & scalstart_neodymium, scalend_neodymium, init_source_netcdf,
     & coeff_init
#endif


      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use ice_mod
      use dynami_mod

! --- BdB 05-2019: add first year in clio for writing time in output
      use comemic_mod, only: init_year_clio

!dmr for updating of the tracer & speed masks
      use update_clio_bathy_tools, only: get_clio_masks_fnm
#if ( BATHY >= 1 )
!      use update_clio_bathy_tools, only: get_clio_masks_fnm, read_clio_masks
     >         , read_clio_masks
     >         , update_diff_masks
      use global_constants_mod, only: str_len
      use unix_like_libs, only: is_file
#endif


       use mchd99_mod, only: nn99
       use clio_control_mod, only: ktvar,ntrmax,mixage, yrsec, daysec
     >                           , dtsd2, unsplt

       use scale_mod, only: scale_init

#if ( BRINES == 2 )
      use global_constants_mod, only: str_len
      use brines_mod, only: read_frac2D
#elif ( BRINES == 3 )
      use global_constants_mod, only: str_len
      use brines_mod, only: read_frac2D, la_date_brines
#elif ( BRINES == 4 )
      use brines_mod, only: read_frac1D, la_date_brines, get_brines_frac
      use global_constants_mod, only: str_len
#endif

      use newunit_clio_mod, only: clio3_out_id

!dmr --- Added the testing of to_from_CLIO
      use to_and_from_clio, only: get_indexes_C, get_lonlat_C
!dmr --- Added the testing of to_from_CLIO

      USE OCEAN2COUPL_COM, only: initseaalb

!! END_OF_USE_SECTION



      integer(kind=ip):: jflag, nn2t, nn3t, numcpl
      real(kind=dblp) :: xday

!dmr @-@ iceb0
!mb added this to read correct restart file res******.om
      integer irunlabelclio
!nb debut
!      integer use_brines
#if ( BRINES >= 1 )
      real*8 fluxbrines(imax,jmax)
      common / common_brines / fluxbrines
#endif
#if ( BRINES >= 2 )
      character(len=str_len) :: frac_filename 
#endif
#if ( BRINES == 3 )
      character(len=5) :: charI
#endif
!nb fin

!--local variables :
      integer(kind=ip), dimension(imax,8) :: irn
      integer(kind=ip), dimension(jmax,8) :: jrn      
      integer(kind=ip) :: iyear
c~       dimension irn(imax,8), jrn(jmax,8)
      character*6 chf

#if ( BATHY >= 1 )
!dmr for updating of the tracer & speed masks
      character(len=str_len) :: file_path
      logical                :: logical_elmt
!dmr
#endif

!dmr --- Added the testing of to_from_CLIO
      integer :: input_i = 103, input_j = 50, indexes_C(2)
      real, dimension(2) :: lon_lat_C
!dmr --- Added the testing of to_from_CLIO


      open (newunit=clio3_out_id,file='outputdata/globals/clio3.out')
      write(clio3_out_id,*) 'clio3 : Adv. Alterne X,Y (nsewfr)'
!ic0  write(clio3_out_id,*) 'om : Adv. Alterne X,Y (nsewfr)'

      call foroutp(irn,jrn)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1) Preparation of the run.                                          |
!-----------------------------------------------------------------------

      jflag = 0
      nn99 = 0
      call defcst(nn99)
#if ( BATHY >= 1 )
!dmr --- reading the previous tms, tmu if it exists in the startdata
      file_path='startdata/'//trim(get_clio_masks_fnm(irunlabelclio))
!      file_path=trim(get_clio_masks_fnm(irunlabelclio))
      write(*,*) 'read_clio_mask', file_path
      if (is_file(trim(file_path))) then 
         call read_clio_masks(file_path)
      else
         write(*,*) '!! no restart for tms and tmu'
      endif
!dmr
#endif

      call defgrid(nn99)

#if ( BATHY >= 1 )
!dmr --- updating the differences in tms/tmu masks
      logical_elmt = update_diff_masks(tms,tmu)
!dmr
#endif
      call redforc(nn99)
      call initseaalb
!--Zonaly uniform, time dependant forcing :
      ktvar = abs(kforc) / 100

!--kstart = 0 start from routine start   (ocean at rest, no ice)
!--       = 1                    redrun* (follow up of a run, same conditions)
!--       > 1 transition :
!--       = 2                    staoc*  (ocean not at rest, sea ice prescribed)
!--kinput = 0/2 (nn2t=0) restart from the binairy file (*b)
!--       = 1/3 (nn2t=1)                  NetCDF       (*c)
!--koutpu = 0/2 (nn3t=0) output on a binary file (*b)
!--       = 1/3 (nn3t=1)             NetCDF  (*c)

      nn2t = mod(kinput,2)
      nn3t = mod(koutpu,2)
      ntrmax = int(dts(ks2)/ddtb)
! mab: dts(ks2)=86400, ddtb=86400
      if ((ntrmax .ne. 1) .and. (icoupl.ge.1)) then
          write(clio3_out_id,*) 'problem with ntrmax in  coupled mode ',ntrmax
          stop
      endif
      if (kstart.eq.0) then
!- case kstart = 0 :
        call start
      elseif (kstart.eq.2) then
        if (nn2t.eq.0) then
          call staocb(kinput,'resto.om')
        else
!         call staocc(kinput,'resto.ncdf')
        endif
      else
        if (nn2t.eq.0) then
!MB added irunlabel to restart filename
           write(*,*) "CLIO restart type:", kinput, nn99
           write(chf,'(I6.6)') irunlabelclio
           call redrunb(kinput,nn99,'startdata/res'//chf//'.om')
        else
!ncd      call redrunc(kinput,nn99,'rest.ncdf')
          write(clio3_out_id,*) 'om : Version without NetCDF '//
     &               '=> reading rest.ncdf impossible !'
          stop
        endif
      endif

! mab: icoupl .eq. 2 means "set time to zero"
      if (icoupl.ge.2)  then
        numit = 0
        tpstot = 0.0
      endif
      nstart = numit + 1
      nlast  = numit + nitrun
! mab: yeaday=360,defined in clio/sources/defcst
      yrsec  = yeaday*86400.
      daysec  = 86400.
      dtsd2 = 0.5*dts(ks2)
      iyear = 1+int(tpstot/yrsec)
! --- BdB 05-2019: initialise the first year, needed for writing time
      init_year_clio = 1+int(tpstot/yrsec)
      xday = mod(tpstot,yrsec)/daysec
!      nsav  = nlast
      write(clio3_out_id,'(A,I8,A,F7.2,A,I10)') ' start of the run ; year=',
     &     iyear,' +', xday, ' days ; iteration=', nstart
      if (ninfo.eq.0) ninfo = nlast + 1
      if (nsav.eq.0)  nsav  = nlast + 1
      mixage = max(min(lstab+lstab,2),-lstab)
      unsplt = 1.0 / DFLOAT(nsplit - nsplaj)

!sai  call splsclr1
      call etat(0, nn99)
      call informe(nn99)
      if (nstart.eq.1) call inforun(nn99)
!     call conti3d

!- Initialisation of the coupling
      if (icoupl.ge.1) then
!cp1     call mytid_ocn
!cp1     call ocn_rfm_cpl
         numcpl=0
!cp2     call inicmo2(nitrun,freqcpl,int(dts))
!         call inicmo3(nitrun,freqcpl,int(dts))
!         write(clio3_out_id,*) ' ocn recv initi forc from cpl, itr :', nstart
      endif

#if ( ISOOCN >= 1 )
      call init_isoocn
#endif



#if ( PATH >= 1 )
        call path_init(scal(:,:,:,scalstart_path:scalend_path))
#endif

#if ( NEOD >= 1 )
        call neodymium_init(scal(:,:,:,scalstart_neodymium:scalend_neodymium))
        call init_source_netcdf
        call coeff_init
#endif



        call scale_init()


#if ( BRINES == 2 )
      frac_filename='inputdata/clio/frac.nc'
      call read_frac2D(frac_filename)
#elif( BRINES == 3 )
      write(*,*) 'The date for brines update is ', la_date_brines

      write(charI,'(I5.5)'), la_date_brines
      frac_filename ='inputdata/clio/frac_clio_'//trim(charI)//'.nc'

      write(*,*) frac_filename
      call read_frac2D(frac_filename)
#elif( BRINES == 4 )
      call read_frac1D
      call get_brines_frac
#endif


!dmr --- Added the testing of to_from_CLIO
!dmr --- 2022-05-05
      write(*,*) "given, i, j          ", input_i, input_j
      
      lon_lat_C = get_lonlat_C(input_i, input_j)
      
      write(*,*) "obtained, lon, lat ==", lon_lat_C(1), lon_lat_C(2)

      indexes_C = get_indexes_C(lon_lat_C(1), lon_lat_C(2))
             
             
      write(*,*) "obtained, i, j ==    ", indexes_C(1), indexes_C(2)
!dmr --- Added the testing of to_from_CLIO



      end subroutine init_CLIO

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


#if ( F_PALAEO >= 1 && UNCORFLUX >= 1 )
      subroutine init_CLIO_ormen()


!! START_OF_USE_SECTION

      use zonesOrmenfwf_mod, only: init_zonesOrmenfwf_mod
      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use dynami_mod
!! END_OF_USE_SECTION

!! START_OF_INCLUDE_SECTION

#include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "dynami.com"

!! END_OF_INCLUDE_SECTION

      call init_zonesOrmenfwf_mod(area)

      end subroutine init_CLIO_ormen
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      SUBROUTINE CLIO(istep,irunlabelclio,ntotday)

! mab:clio is called in emic.f once a day at the last atmospheric step; istep=i
! i=1,ntotday

      use global_constants_mod, only: dblp=>dp, ip

!dmr [DEPRECATED]      use comemic_mod, only: flgicb
      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use ice_mod
      use dynami_mod
      use scale_mod, only: scale

      use update_clio_bathy_tools, only: get_clio_masks_fnm, write_clio_masks
!#if ( BATHY >= 1 )
!dmr --- added for bathymetric ouput writing
!      use update_clio_bathy_tools, only: get_clio_masks_fnm, write_clio_masks
!#endif

      use engtur_mod, only: engtur

      use mchd99_mod, only: nn99
      use clio_control_mod, only: ntrmax,mixage,dtsd2,yrsec,daysec
     >                           ,unsplt

      use newunit_clio_mod, only: clio3_out_id

#if ( ICEBERG > 0 ) 
      use iceberg_mod, only: write_iceberg_restartdata 
#endif

#if ( BRINES >= 1 )
      use brines_mod, only : brines
#endif


! APW ? -> dmr : what is the use of liftoff ?
!      integer liftoff
! APW ?
      character*6 fchnum
      integer irunlabelclio
c~       common /clio_control_int/ ktvar,ntrmax,mixage
c~       common /clio_control_fp/ dtsd2,yrsec,daysec,unsplt

! APW ? -> dmr : what is the use of liftoff ?
!      common /dbug_numit_com/ liftoff
! APW


      integer(kind=ip):: nsew, ither, nuclin, nn3t, ncfch
      real(kind=dblp) :: xjour, xjour1, ccsplt
      integer(kind=ip) :: iyear
      integer :: istep, ntotday
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2) Daily do loop.                                                   |
!-----------------------------------------------------------------------
! APW ? -> dmr : what is the use of liftoff ?
!        if (liftoff.eq.0) liftoff=numit
! APW ?

        numit = numit + 1
! modif for  the counting of the interations in couple mode
        nsew = mod(numit,nsewfr)
!--The date corresponds to the middle of the time step.
!--iyear is the year, the first year is the year 1.
!--The integer part of xjour is the day of the year(1-yeaday) and the
!--non-integer one is the fraction of this day.
! mab: is dtsd2 half a day? (86400 / 2)
        tpstot = tpstot + dtsd2
        iyear  = 1+int(tpstot/yrsec)
        xjour  = mod(tpstot,yrsec)/daysec +1.0
        xjour1 = mod(tpstot+dts(ks2),yrsec)/daysec +1.0
        tpstot = tpstot + dtsd2

!-- 2.1. start of an iteration :
!-------------------------------

!ic0    if (ktvar.ne.0) call tvforc(nn99)
!sai    call splsclr2(jflag,iyear,xjour)
#if ( I_COUPL == 0 )
!        if (icoupl.eq.0) then
           call forcat(iyear,xjour)
#elif ( I_COUPL == 3 )           
!        elseif ((icoupl.eq.3) .and. (numit.eq.nstart)) then
        if (numit.eq.nstart) then
!cp2       call rfmcpl2(numcpl)
           call forcat(iyear,xjour)
        endif
!       elseif (icoupl.ge.1) then
!        else
#endif
!cp2       call rfmcpl2(numcpl)
!           call rfmcpl3(numcpl)
!        endif
        call icdyna
        call icdadv(xjour)

!cfc    call cfc_flx(nn99)

        do ither=1,ntrmax
           call thersf(ither,istep)
!dmr&afq           call thersf(ither,ntrmax,istep,irunlabelclio)

!cp1       if (icoupl.ge.1.and.ither.lt.ntrmax) then
!cp1          call ocn_s2_cpl
!cp1          call ocn_rfm_cpl
!cp1       endif
!nb debut
#if ( BRINES >= 1 )
      call brines
#endif
!nb fin
        enddo
!       call ocesla
!tk0    call vdiffu(nn99)
        call diftur
        call engtur
        call conti3d
        call flucor
        call isoslope
        if (aitd(kmax).gt.zero) call isoadvec
        call scale(nsew,nn99)
        call scali
        call etat(mixage, nn99)
        do 250 nuclin=1,nclin
          call uve
          call barot
          ccsplt = 0.0
          do 200 numspl=1,nsplit
            if (numspl.gt.nsplaj) ccsplt = unsplt
            if (ahe.eq.zero) then
              call uvb0et(ccsplt)
            else
              call uvbfet(ccsplt)
            endif
 200      continue
!         call uvi   ! --> ds uvm
          call uvm
 250    continue
!       write(mouchard_id,*) 'clio',toticesm

!--partie ICEBERG de l'ITERATION
!--------------------

#if ( ICEBERG > 0 )
        call iceberg(istep)
	if ( mod(istep,ntotday).eq.0) then 
          call write_iceberg_restartdata 
        endif 
#endif

!--end of an iteration

!-- 2.2. Outputs.
!-----------------
!dmr ###        if (ntout.eq.1) then
        call outnet(iyear,xjour,xjour1)
!dmr ###        else
!dmr ###        call outave(iyear,xjour,xjour1)
!dmr ###        endif
        if (mod(numit,ninfo).eq.0 .or. ntmoy.ge.1) call inforun(nn99)

!#if ( ICEBERG == 1 )
!         call iceb_out
!#endif

        if (mod(numit-nstart+1,nsav).eq.0 .or.(numit-nstart+1).eq.nlast)
     >   then
! APW ? -> dmr : what is the use of liftoff ?
!        if (mod(numit-liftoff,nsav).eq.0) THEN
!          write(mouchard_id,*)'APW clio.f:(mod(numit-liftoff,nsav).eq.0)',
!     &    numit-liftoff,nsav
! APW ? -> dmr : what is the use of liftoff ?
!--definition of the root f or the name of the restart file:
!           indicf = (nlast - numit) / nsav
!           ncfch = max(1,indicf)
!           ncfch = 4 - int( log10( DFLOAT(ncfch) ) )
           write(fchnum,'(I6.6)') irunlabelclio
!--   vertical velocity "w" n agreement with (u,v) :
           if (koutpu.ge.2) call conti3d

!--   save of teh results
           nn3t = mod(koutpu,2)
           if (nn3t.eq.1) then
!     ncd        call savrunc(koutpu,nn99,'res'//fchnum(ncfch:)//'.ncdf')
              write(clio3_out_id,*) 'om : Version without NetCDF '//
     &   '=> writing of  res'//fchnum(ncfch:)//'.ncdf impossible !'
!              call savrunb(koutpu,nn99,'res'//fchnumch:)//'.om')
              call savrunb(koutpu,nn99,'res'//fchnum//'.om')
#if ( ISOOCN >= 1 )
              call write_isoocn
#endif
!     dmr --- Fix for the output
!              num_clio_res = fchnum(ncfch:)
!              WRITE(*,*) "num_clio_res : ", num_clio_res, fchnum(ncfch:)
!     dmr --- End fix
           else
!              call savrunb(koutpu,nn99,'res'//fchnum(ncfch:)//'.om')
              call savrunb(koutpu,nn99,'res'//fchnum//'.om')
!     dmr --- Fix for the output
!              num_clio_res = fchnum(ncfch:)
!              WRITE(*,*) "num_clio_res : ", num_clio_res, fchnum(ncfch:)
!     dmr --- End fix
#if ( ISOOCN >= 1 )
              call write_isoocn
#endif

!nb always write mask for restart
!#if ( BATHY >= 1 )
!dmr --- Tentative positionning of the tracer mask file writing
          call write_clio_masks(trim(get_clio_masks_fnm(irunlabelclio)))
!dmr ---
!#endif

           endif
           write(clio3_out_id,*) 'deriv:',deriv
           write(clio3_out_id,*) 'day:',xjour,' ------ ','year:',iyear
!--   buffeur associed with evolu emptied (No=90) :
           if (numit+ninfo.le.nlast) call flush(90)
        endif

!-receive  and/or  send forcing in coupled mode
!        if (icoupl.ge.1) then
!cp1       call ocn_s2_cpl
!cp2       call s2cpl2(numcpl)
!           call s2cpl3(numcpl)
!           write(clio3_out_id,*) 'sst and ice has been sent to cpl', numit
!cp1       call ocn_rfm_cpl
!cp1       write(clio3_out_id,*) 'fluxes have been received by sioclm', numit
!        endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3) End of the run.                                                  |
!-----------------------------------------------------------------------

!      xday = xjour - 1. + dtsd2/daysec
!      write(clio3_out_id,'(A,I8,A,F7.2,A,I10)') '  end of the run ; year=', iyear,
!     &    ' +', xday, ' days ; iteration=', numit
!      if(nn99.eq.2) close(99)

!-deconnection of sioclm from PVM in coupled mode
!      if (icoupl.ge.1) then
!cp1     call ocn_stop_ocn
!cp2     call quitcmo2
!         call quitcmo3
!      endif

!      stop
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- end of the routine clio -
      end
