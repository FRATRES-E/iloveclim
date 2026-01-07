!     dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:27 CET 2009
#include "choixcomposantes.h"
!     dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:27 CET 2009


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
! dmr  MAIN PROGRAM FOR THE iLOVECLIM COUPLED EARTH SYSTEM MODEL
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|

      program emic

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
! dmr  By reference variables and functions
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|

      use comatm, only: nlat, nlon, darea, dtime
      use global_constants_mod, only: days_year360d_i

! flgveg    (flag vegetation, i.e. vecode = on or off)
! iatm      (numbers of atmospheric time steps in one day)
! iday      (counts the day during the run)
! imonth    (counts the months during the run)
! ilan      (number of land timesteps in one day)
! irunlabel (CLIO identifier)
! iyear     (counts the years during the run
! ntotday   (total number of days in proposed integration)

      use comemic_mod, only: flgveg, iatm, iday, imonth, ilan, irunlabel, iyear
     &                     , ntotday, nstpyear

#if ( ISM >= 2 )
      use comemic_mod, only: nwrskip, day
#endif

! --- BdB 05-2019: added new_year_atm and new_year_veg
      use comemic_mod, only: new_year_atm, new_year_veg, time_in_years
     &                     , current_int_atm, current_int_veg

      use comemic_mod, only: pretty_print_exec_time

! dmr comsurf is needed for the variables:
! epss (epsilon, small?)
! nld  (n land = number for land in mixed arrays)
! noc  (n ocean, idem)
! nse  (n sea ice, idem)
      use comsurf_mod, only: noc, nse, nld, epss, fractn, tempsgn


#if ( ISM == 2 || ISM == 3 )
!     dmr FLAG AJOUT GRISLI
!     mab:input_timerCplGRIS defines the frequency of coupling between
!     loveclim couples once a day, and if freqdecoup=1 it is coupled
!     if it s bigger than 1 it s decoupled
      USE input_timerCplGRIS
!     mab:input_flagsGRIS defines kind of calotte (0=no,1+variable but
!     forced,2=interactive) nord (=2)/sud(=0)
      USE input_flagsGRIS
!     dmr FLAG AJOUT GRISLI
#endif

#if ( DOWNSTS == 1 || DOWNSCALING == 2 || ISM >= 2 )
      use sgout_mass_balance_mod, only: initakkuVars, akkuVars,
     &       sgout_smb
      use transfer_subgrid_ecb_mod, only: subgrid_ecb_wrapper
      use transfer_ecb_subgrid_mod, only: transfer_wrapper_subgrid
#endif

#if ( DOWNSTS == 1)
      use vertDownsc_mod  , only: create_surftemp_var_d
#endif

! BdB 02-2019: in coherence with calls below (for F_PALAEO == 1),
!              remove ISM > 0 option when F_PALAEO >= 1
#if ( ISM == 2 || ISM == 3 || F_PALAEO == 1 || DOWNSCALING == 2 )
      use ecbilt_topography, only: ec_topo
#endif

#if ( CYCC >= 2 )
      use carbone_co2, only: PA0_C, PA_C
#endif

#if ( KC14 == 1 )
      USE mod_sync_time, ONLY: KENDY ! FLAG for end_of_year is true
#endif
#if ( F_PALAEO >= 1 )
      use palaeo_timer_mod, only: reload_topo, palaeo_timer
#endif
#if ( BATHY >= 2 )
      use update_clio_bathy_tools, only: reload_bathy, la_date
#endif
#if ( PATH >= 1 )
      use path_mod, only: path_step
#endif
#if ( NEOD > 0 )
      use neodymium_mod, only: neodymium_step
#endif


#if ( MEDUSA == 1 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
!     @mohr  Adding some libraries to use MEDUSA -
!            New version with sediments coupled
!            Last modification: 2016-10-25
!--- dmr 2020-06-02 proper declaration of variables for MEDUSA

      USE flux_to_sediments_mod,   only: fait_pointer_o2S, flux_o2s_ini
     >                                  ,sedfluxave
      USE flux_from_sediments_mod, only: flux_s2o_ini, update_flx_from_sediments
      USE medusa_wrapper,          only: i_medusastep, ayears_medusacurrent
     >                  , ayears_medusaini, ndays_medusastep, ayears_medusastep
     >                  , medusa_timer, nmedusasteps_ncoutfrequency
     >                  , medusa_wrap_init, medusa_wrap_step, medusa_wrap_close
     >                  , medusa_wrap_prep
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
      USE flux_from_sediments_mod, only : restart_flx_from_sediments
#endif /* On MEDUSA */

#if ( CARAIB > 0 || CARAIB_FORC_W > 0 )
      use ec_co2ca, only: climvars_acc, reset_climvars
#endif

#if ( CARAIB > 0 )
      use ec_co2ca, only: co2_for_CARAIB, sync_coupler_caraib
     >                  , orbitals_for_CARAIB
      use caraib_cpl, only: global_loop, end_run
#endif

      use ipcc_output_mod, only: ipcc_output

#if ( PROGRESS == 1 )
      use basic_libs,      only: progress
#endif

      use initcoupledmodel_mod,  only: ec_initemic, init_coupled_components
      use emic_write_state_mod,  only: ec_writestate
      use atmos_composition_mod, only: set_PGA_CO2, get_PGA_CO2, get_PCO2_REF
     &                          , atmos_carbon_update, get_PCO2_VEG
#if ( OCYCC == 1 )
      use ocycc_main,            only: ocycc_step
#endif

#if ( CLM_INDICES >= 1 )
      USE CLIMATE_INDICES_MOD, ONLY: GLOBAL_RE_INIT, GLOBAL_FINALIZE
#endif

#if ( CLIO_OUT_NEWGEN == 1 )
      USE GRID_IO_NC, ONLY: grid_io_reinit=>GLOBAL_RE_INIT
     &                    , daily_io_nc=>DAILYSTEP_IO_NC
     &                    , nbyearsinfile
#endif
!!    USE OMP_LIB

#if ( forced_winds == 1 )
      use offline_wind_forcing, only: update_winds_off
#endif

#if ( WINDS_ERA5 == 1 )
      use WINDFORC_CLIOERA5, only : get_daily_ERA5_CLIO_WINDSSTRENGTH
      use bloc_mod, only: normUV_ERA5
#endif

#if ( IRON_LIMITATION == 1 )
      use IRONLIM_PISCES, only : get_PISCES_CLIO_IRONLIM
      use bloc_mod, only: LFe_PISCES
#endif

      USE COUPL2OCEAN_COM, only: ec_co2oc
      USE OCEAN2COUPL_COM, only: ec_oc2co

#if ( FROG_EXP > 0)
      use main_lib_FROG, only: INITIALIZE_VAMP, GET_COUPLING_STEP
     >                         , STEPFWD_VAMP, INITIALIZE_CARBON_STOCK

      use CPL2FROG_mod,    only: INIT_CPL2VAMP, GET_VAMPVARS
     &                          , DAILY_UPDATE_VAMPVARS 
     &                          , RESET_VAMPVARS_TIMER
#endif


      use landmodel_mod, only: ec_la2co, ec_co2la, ec_lbm, ec_lae2co
c~      >                       , ec_sumfluxland


      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
! dmr  Local variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|

      integer i,j,k
      logical :: result_function_call

! dmr progname is used for pretty printing of the progress bar ...
#if ( PROGRESS == 1 )
      character*10 progname
#endif

#if ( MEDUSA == 1 )
      LOGICAL              :: l_ncout = .true. ! dmr --- added a default value 2020-06-02
#endif

#if ( FROG_EXP > 0)
      logical :: well_done
#endif

!cnb caraib CYCC
!#if ( CARAIB > 0 )
!      integer pix_car
!      parameter (pix_car=622)
!      common /veccarb/ stock_ysoilr(pix_car)
!      real*4  stock_ysoilr
!      common /inidata/ ylit_ini(2,pix_car)
!     >                ,yhum_ini(pix_car)
!      real*4 ylit_ini,yhum_ini
!      integer ngt
!#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|
!       Main code of the program starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8--|

      call pretty_print_exec_time("Start Execution")

!$OMP PARALLEL
      WRITE(*,*) "ECHO OMP ..."
!$OMP END PARALLEL

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|
! dmr   Call of the different components of init => place in a general INIT routine?
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|

cnb try to call first to have the date t update bathy
#if ( F_PALAEO >= 1 && ( BATHY >= 2 || NC_BERG == 2 ) )
         write(*,*) "call palaeo_timer with init=2", irunlabel
         CALL palaeo_timer(2,irunlabel,n_days=ntotday)
#endif

      result_function_call = ec_initemic()

      result_function_call = init_coupled_components(iday,imonth)

#if ( forced_winds == 1 )
      result_function_call = update_winds_off()
#endif

#if ( WINDS_ERA5 == 1 )
      normUV_ERA5(:,:,:)=get_daily_ERA5_CLIO_WINDSSTRENGTH()
#endif

#if ( IRON_LIMITATION == 1 )
      LFe_PISCES(:,:,:)=get_PISCES_CLIO_IRONLIM()
#endif



#if ( MEDUSA == 1 )
!-----|--1--------2---------3---------4---------5---------6---------7-|
!     dmr   Initialization of the grid sizes for CLIO and MEDUSA
!     mohr  Coupling the full medusa model 2016-10-25
!-----|--1--------2---------3---------4---------5---------6---------7-|

      write(*,*) 'MEDUSA OPTION ACTIVATED' !mohr

      !tsup = 0.d0
      !temps = 0.d0   ! in days, double precision
      !timestep = 3600 ! in days, integer 10 years
      !timestep = 360
      !tinf = 0.d0

      CALL fait_pointer_O2S()       ! From FLUX_TO_SEDIMENTS_MOD:
                                    !  - connects *_ma pointers from MBIOTA_MOD
                                    !    onto their targets in MARINE_BIO_MOD
                                    !    and LOVECLIM_TRANSFER_MOD

!     ! initialize the accumulator for the fluxes ocean to sediment to zeros

      CALL flux_o2s_ini()           ! From FLUX_TO_SEDIMENTS_MOD:
                                    !  - sets to zero all the *_mave
                                    !    variables from MBIOTA_MOD
      CALL flux_s2o_ini()           ! From FLUX_FROM_SEDIMENTS_MOD:
                                    !  - sets to zero all the *_sed2oc
                                    !    and *_loopback variables
                                    !    from MBIOTA_MOD

      CALL medusa_wrap_init()       ! From MEDUSA_WRAPPER:
                                    !  - opens Medusa's ERR, LOG and DBG files
                                    !  - calls INIT_FILELIST_MEDUSA
                                    !    (init list of  NC files)
                                    !  - calls INIT_TIMECONTROL_MEDUSA
                                    !    (init. time control vairbales and
                                    !    characteristic interval lengths)
                                    !  - calls SETUP_ILOVECLIM_MEDUSA_XCHANGE
                                    !    (prepares data for SEAFLOOR_SETUP)
                                    !  - calls SEAFLOOR_SETUP
                                    !  - calls InitEquilibParameters (from mod_equilibcontrol.F)
                                    !  - calls InitProcessParameters (from mod_processcontrol.F)
                                    !  - calls SETUP_SEDCORE_SYSTEM (from mod_sedcore.F)
                                    !  - calls setup_o2s (from mod_iloveclim_o2s.F)
                                    !    (allocates arrays in MOD_ILOVECLIM_O2S)
                                    !  - calls setup_s2o (from mod_iloveclim_s2o.F)
                                    !    (allocates arrays in MOD_ILOVECLIM_S2O)
                                    !  - calls ini_fluxes_s2o (from mod_iloveclim_s2o.F)
                                    !    (sets arrays in MOD_ILOVECLIM_S2O to zero)


      CALL medusa_wrap_prep()       ! From MEDUSA_WRAPPER:
                                    !  - calls InitSeafloorFrom*File*
                                    !    and calls SEDIMENT_TO_OCEAN if required
                                    !  - calls ini_fluxes_o2s (from mod_iloveclim_o2s.F)
                                    !    (sets arrays in MOD_ILOVECLIM_O2S to zero)
                                    !  - calls OCEAN_TO_SEDIMENT
                                    !  - calls OPEN_NCFILES_MEDUSA(ayears_medusaini)
                                    !    (creates NC files and writes initial state)

                                    ! Initialize time control counters and
                                    ! clocks for Medusa

!     !the following variables are declared in medusa_wrapper.F
      i_medusastep         = 0
      ayears_medusacurrent = ayears_medusaini
      ndays_medusastep     = MIN(ndays_medusastep, ntotday)
      ayears_medusastep    = DBLE(ndays_medusastep)/360.0D+00


      write(*,*) "vars MEDUSA == =", i_medusastep, ayears_medusacurrent
     >                             , ndays_medusastep, ayears_medusastep

!nb read restart for sediment fluxes
      call restart_flx_from_sediments(1)

#endif

#if ( CARAIB > 0 )
      call global_loop(0,sync_coupler_caraib)
!      call fait_pointer_CARAIB
#endif

#if (DOWNSTS == 1 || DOWNSCALING == 2)
      call Init_subgrid
#endif

!END OF INITIALIZATION

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|
! dmr  Init of the palaeo _timer forcing used to apply scenarios
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|

#if ( F_PALAEO >= 1 )
         write(*,*) "call palaeo_timer, init", irunlabel
         CALL palaeo_timer(1,irunlabel,n_days=ntotday)
#endif

#if ( F_PALAEO >= 1 && UNCORFLUX >= 1 )
         CALL init_CLIO_ormen()
#endif

#if ( F_PALAEO == 1 )
         CALL load_topo_masq(.TRUE.) ! parameter given is init
         CALL ec_masq(.FALSE.)       ! not an init
         CALL ec_topo                ! update topography
#if (DOWNSCALING == 2 )
         call subgrid_ecb_wrapper
         call ec_masq(.FALSE.)
         call ec_topo
#endif
#elif ( F_PALAEO == 3)
         call load_topo_masq_hr(.TRUE.) ! afq -- we read a prescribed topography
         CALL ec_topo
         call subgrid_ecb_wrapper
         call ec_masq(.FALSE.)
         call ec_topo
#elif ( ISM >= 2  || DOWNSCALING == 2 )
         call ec_topo                ! update topography, afq (rmount etc.)
         call subgrid_ecb_wrapper
         call ec_masq(.FALSE.)
         call ec_topo
#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|
! dmr&afq  Tentative placement of the CLIO masq call
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|

         CALL clio_masq

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|
! afq  Initialization for the subgrid downscaling
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|

#if ( DOWNSTS == 1 || DOWNSCALING == 2 )
         CALL initakkuVars
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|
! dmr  Initialization for GRISLI and accordingly the coupling procedures
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|

#if ( ISM == 2 || ISM == 3 )

! mab: 0=without model of the ice sheet
!      1=use agISM  (without starting it)
!      2=utilize and start GRISLI
!      3=ice-sheet coupled but fixed

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|
!     dmr   Ajout de l initialisation du couplage ECV - GRISLI (LUDUS)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|

      WRITE(*,*) "INIT COUPLAGE ECV - GRISLI"

      freqsortie = nwrskip*freqdecoup

      IF ((nord_GRIS.GE.1).OR.(sud_GRIS.GE.1)) THEN

!-----|--1--------2---------3---------4---------5---------6---------7-|
!     dmr   Ajout de l initialisation du modele GRISLI
!-----|--1--------2---------3---------4---------5---------6---------7-|
         WRITE(*,*) "INIT GRISLI"
         IF (nord_GRIS.GE.1) CALL GRISLI_INIT

!     dmr La ligne suivante lit la topo GRISLI et le masque calotte GRISLI
!     dmr dans un fichier externe. A reutiliser pour un forcage offline ?
!     dmr         CALL input_VAR_GRISLI

         WRITE(*,*) "INTERPOL GRISLI => L"

! dmr --- added subgrid_ecb_wrapper call back. Do not replace Interp yet.
         call subgrid_ecb_wrapper()

!     mab: in this subroutine the new icemask is defined (1,if there is ice,
!     otherwise 0) the icemask is also updated in regard of the ocean (there
!     no ice possible and a differentiation between greenland and antarctica
!     is made. if the flag is set to .TRUE. then a parametrised mask is used)
         CALL ec_masq(.FALSE.)
         CALL ec_topo

       WRITE(*,*) "===================================="
       WRITE(*,'(a,i8,a)') "Coupling timestep iLOVECLIM to GRISLI is ",
     &                    timCplGRISday, "  days"
       WRITE(*,*) "===================================="
       ENDIF
#endif

#if ( CARAIB > 0 || CARAIB_FORC_W > 0 )
       call reset_climvars()
!cnb and init of terrestrial biosphere carbon in CARAIB
!       open(25,file='/home/climwork/textier/caraib-git/wkdir/
!     >testCARAIBnewISO_1000y_v5/biomass1000.res',
!     >       form='unformatted')
!       do ngt = 1, pix_car
!           call read_init(2,ngt) !iread=0 -> seulement npft0
!           stock_ysoilr(ngt) = ylit_ini(1,ngt)+ylit_ini(2,ngt)
!     >                         +yhum_ini(ngt)
!           write(*,*) 'carbon veget dans emic', ylit_ini(1,ngt),
!     >        ylit_ini(2,ngt), yhum_ini(ngt)
!       enddo
!       close(25)
#endif


!cnb init of carbon cycle moved here after change of mask and topo
#if ( CYCC >= 2 )
!-----|--1--------2---------3---------4---------5---------6---------7-|
!     dmr   Initialization of atm. CO2 in in Coupled carbon cycle mode
!-----|--1--------2---------3---------4---------5---------6---------7-|
      CALL ECO2(0,fractn(1,1,nld),darea)
      WRITE(*,*) "emic.f (II) : PA0_C, PA_C", PA0_C, PA_C
#endif

! --- BdB 05-2019: initialise new_year_atm and new_year_veg
      new_year_atm = 0
      new_year_veg = 0
      current_int_atm = 0
      current_int_veg = 0


#if ( CLM_INDICES >= 1 )
      CALL GLOBAL_RE_INIT()
#endif


#if ( FROG_EXP > 0 )
      well_done = INITIALIZE_VAMP()
      if (well_done) then
        WRITE(*,*) "FROG INITIALIZATION COMPLETE"
      endif
      well_done = INIT_CPL2VAMP(GET_COUPLING_STEP())

      call DAILY_UPDATE_VAMPVARS()
      well_done = INITIALIZE_CARBON_STOCK(GET_VAMPVARS())
      call RESET_VAMPVARS_TIMER
#endif

#if ( CLIO_OUT_NEWGEN == 1 )
!dmr&ec First initialization of first file, initialize with index zero
      CALL grid_io_reinit(0)
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|
! dmr  Here all the inits are finished. We start the main integration
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|

!     *** OCEAN-SEAICE >>

! dmr  i = the counter of days in the integration
      do i=1,ntotday

#if ( PROGRESS == 1 )
         if(i.eq.1) write(*,*)
         write(progname,'(A,I6)') 'YEAR',int(irunlabel+(i/360))
         call progress(progname,i-1,ntotday-1)
#endif

!     mab: oceanic data to coupler
        call ec_oc2co(i)

!*** PERTURBATION OF THE ATMOSPHERE VIA tsurfn
#if ( PERTATMOS == 1 )
        if (i.eq.1) call ec_perttsurfn
#endif

!     ***   ATMOSPHERE >>

! dmr j = the counter for atmospheric timesteps
         do j=1,iatm

            call ec_update(i,j)
            call ec_co2at
!#if ( CARAIB > 0 || CARAIB_FORC_W > 0 )
!            call climvars_acc(modulo(i-1,360)+1,iatm,j)
!#endif
            call ec_at2co
            call ec_fluxes(noc)
            call ec_fluxes(nse)


!     *** integrate atmosphere

            call ec_ecbilt(i,j)


!     ***     LAND >>
!     mab: ilan=1
            do k=1,ilan
               call ec_la2co
               call ec_fluxes(nld)
               call ec_co2la

!     *** integrate land

               call ec_lbm(i,j,k)
               call ec_lae2co
               call ec_sumfluxland(k)
            enddo               ! boucle sur k, ilan

!     ***     << LAND

            call ec_sumfluxocean(i,j)

c~ [UNUSED] // code from LLN's iLOVECLIM 1.3 ... maintain ?
c~ !*** NUDGING
c~ #if ( NUDGING == 1 )
c~           call nudging(i,j)
c~ #endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!     *** OCEAN-SEAICE >>
!-----|--1--------2---------3---------4---------5---------6---------7-|

          if (j.eq.iatm) then ! dmr end of the day !

#if ( FROG_EXP > 0)
              call DAILY_UPDATE_VAMPVARS()
#endif


#if ( CARAIB > 0 || CARAIB_FORC_W > 0 )
              call climvars_acc(modulo(i-1,360)+1,iatm,j)
#endif


#if ( OCYCC == 1 )
               result_function_call = ocycc_step(iday,imonth)              ! step in ocycc daily: get the bio and non-bio!
#endif
               call ec_co2oc(i)                                            ! transmit fields coupler to ocean

c~ [UNUSED] // code from LLN's iLOVECLIM 1.3 ... maintain ?
c~ #if ( PERTOCEAN == 1 )
c~ !             !*** PERTURBATION OF THE OCEAN VIA scal
c~ ! !         if (i.eq.1) call ec_pertscal(i)
c~ !           call ec_pertscal(i)
c~             if (i.ne.1) call ec_pertscal(i)
c~             !this is just to write the scal file
c~ #endif

               call clio(i,irunlabel+iyear,ntotday)                        ! integrate ocean-seaice for a daily step (physics)

!#if ( OCYCC == 1 )
!               call sync_lcm_ocycc(1)                                      ! synchronisation des champs lcm (CLIO) -> OCYCC
!#endif

               result_function_call = atmos_carbon_update()                ! update the atmospheric pCO2

          endif ! j == iatm

!-----|--1--------2---------3---------4---------5---------6---------7-|
!     *** << OCEAN-SEAICE
!     *** >> VEGETATION
!-----|--1--------2---------3---------4---------5---------6---------7-|

          if (flgveg) then ! VECODE

            call veget(i,j,dtime,epss,get_PCO2_VEG(),fractn(1,1,nld),
     &         darea,tempsgn(1,1,nld))

          endif               ! flgveg

!     ***     << VEGETATION

      result_function_call = IPCC_output(i,j)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|
! dmr  Starts business for the downscaling
!      This section of the code is called once per atmospheric timestep
!       e.g. every four hours.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|

! afq : we may want to use the downscaling without the ice sheet model...
#if ( DOWNSTS == 1 || ISM >= 2 )
          ! for now (4/7/16) we can't run the ice sheet without _d downscaling
          call create_surftemp_var_d()
#endif

#if ( DOWNSTS == 1 || DOWNSCALING == 2 || ISM >= 2 )
! prepares annual outputs on the subgrid + if ISM>=2 computes the SMB
          call sgout_smb()

!afq subgrid           IF ((nord_GRIS.GE.2).OR.(sud_GRIS.GE.2)) THEN
          call akkuVars(0,iatm)
!afq subgrid           ENDIF
#endif

        enddo                  ! boucle j, 1=>iatm

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|
! tex Call CARAIB and atmospheric oxygen isotopic composition main routines
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|

!cnb call caraib here (before computes PCO2)
#if (CARAIB > 1)
      if (mod(i,days_year360d_i).eq.0) then
         call co2_for_CARAIB
         call orbitals_for_CARAIB
         call global_loop(i/days_year360d_i,sync_coupler_caraib)
#if ( CYCC > 1 )
         call send_caraib2lbm
#endif
         call frac_acc()

#if (OXYISO > 0)
         call d18O_iso()
#endif

#if (D17ISO > 0)
         call D17_iso()
#endif

#if (WAXISO > 0)
         call dDwax_iso()
#endif

         call reset_climvars()
      endif
#elif ( CARAIB == 1 )
      if (mod(i,days_year360d_i).eq.0) then
          call reset_climvars()
      endif
#elif ( CARAIB == 0 && CARAIB_FORC_W > 0 )
      if (mod(i,days_year360d_i).eq.0) then
          call reset_climvars()
      endif
#endif

#if ( FROG_EXP > 0)
      if (mod(i,days_year360d_i).eq.0) then
        !!!! FROG
        well_done = STEPFWD_VAMP(GET_VAMPVARS())
        WRITE(*,*) "CALLED FROG !!!!"
        call RESET_VAMPVARS_TIMER
      endif

#endif



!     ***     << CARAIB & OXYISO & D17ISO & WAXISO

#if ( KC14 == 1 )
      IF ( KENDY.EQ.1) THEN
        CALL C14ATM_DP
      ENDIF
#endif

#if ( CYCC >= 2 )
         CALL ECO2(1,fractn(1,1,nld),darea)
#endif

!     ***   << ATMOSPHERE

#if ( ISM >= 2 || DOWNSTS == 1 || DOWNSCALING == 2 )
!     dmr >> ICESHEETS, GRISLI + afq >> subgrid downscaling (outputs)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!     Interpolation des variables climat ECBilt & CLIO sur la grille
!     ISM lue dans l'Init_L2GRISLI
!-----|--1--------2---------3---------4---------5---------6---------7-|

#if ( ISM < 2 )
         cpl:if (mod(i,360).EQ.0.0) then ! il est temps de sortir les champs subgrid
#else
         cpl:if (mod(i,timCplGRISday).EQ.0.0) then ! il est temps de coupler

!     dmr GRISLI est force ou couple
            IF ((nord_GRIS.GE.1).OR.(sud_GRIS.GE.1)) THEN

               WRITE(*,*) "Je couple a GRISLI", FLOOR(day)+1, timCplGRISday
#endif
               CALL akkuVars(1,iatm)
! -- afq: this is commented as only ann. fields are interpolated now --->
!!     mab: nwrskip=number of years between writing the model state to disk
!               if (mod(iyear,nwrskip).eq.0.or.(iyear.eq.nyears)) then
!!     mab(lcm2ism/sources/) 1 means interpolating AND writing results out
!                  CALL Interp_ECBilt_GRISLI(1)
!               else
!                  CALL Interp_ECBilt_GRISLI(0)
!               endif
!     -- afq            <---
               call transfer_wrapper_subgrid

!     CALL Interp_CLIO_GRISLI

!-----|--1--------2---------3---------4---------5---------6---------7-|
!     Appel du modele GRISLI NORD
!-----|--1--------2---------3---------4---------5---------6---------7-|
#if ( ISM == 2 )
               IF (nord_GRIS.GE.1) CALL ISM_NORD(timCplGRISyr)
#endif

#if ( ISM >= 2 )
            endif               ! sur nord_GRIS ou sud_GRIS
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!     Reinitialisation des variables d'accumulation
!-----|--1--------2---------3---------4---------5---------6---------7-|
            CALL initakkuVars

#if ( ISM >= 2 )
!     dmr GRISLI est couple
            IF ((nord_GRIS.GE.2).OR.(sud_GRIS.GE.2)) THEN

! dmr  Added call to the subgrid_ecb_wrapper. Will replace Interp in the end.
               call subgrid_ecb_wrapper()
!     dmr   Aggregation des champs GRISLI sur ECBilt
!     mab(lcm2ism/sources/): interpolation only for GRISLI NORD???
!     contains commented part about calving!!!

!     dmr   Mise a jour du masque glaciaire ECBilt
               CALL ec_masq(.FALSE.)
!     dmr   Mise a jour de la topographie ECBilt
               CALL ec_topo

            ENDIF
#endif
         ENDIF cpl              ! sur la periode couplage

!     dmr << ICESHEETS, GRISLI
!     mab: end of if(ism==2)
#endif


!     ***     >> MEDUSA

#if ( MEDUSA == 1 )
!-----|--1--------2---------3---------4---------5---------6---------7-|
!     @mohr: calling medusa inside the main loop !2015-06-15 !!
!     trying 2015-06-30 17:05
!     now trying new version 2016-01-25
!     updated now to version 2016-10-25
!dmr --- updated and simplified 2021-11-29
!-----|--1--------2---------3---------4---------5---------6---------7-|

         IF (i == ntotday) THEN
                                    ! On the last day of the simulation experiment:
           medusa_timer = 0         !  - always call Medusa

         ELSE

           medusa_timer =  MOD(i, ndays_medusastep)

         ENDIF


c~          write(*,*) "Called sedfluxave ==", medusa_timer
c~      >             , ndays_medusastep, i

                                    ! Collect the values from OCYCC for sediment
         CALL sedfluxave(medusa_timer, ndays_medusastep)
                                    ! From FLUX_TO_SEDIMENTS_MOD
                                    !  - accumulates the *_ma arrays (from MBIOTA_MOD)
                                    !    into the *_mave arrays (from MBIOTA_MOD)
                                    !  - if medusa_timer==0:
                                    !    + divides *_mave by ndays_medusastep
                                    !    + does the speciation of DIC
                                    !    + copies contents into *_mafond arrays
                                    !      (from MBIOTA_MOD)
                                    !    + calls xchange_fluxes_o2s
                                    !      (from MOD_ILOVECLIM_O2S): transfers the
                                    !      *_mafond arrays from MBIOTA_MOD into
                                    !      MOD_ILOVECLIM_O2S
                                    !    + calls flux_o2s_ini


         IF (medusa_timer == 0) THEN

           i_medusastep = i_medusastep + 1

                                    ! Override the default l_ncout = .FALSE.
                                    ! whenever the output frequency is matched.
                                    ! If i == ntotday, l_ncout is already .TRUE.
           IF (MOD(i_medusastep, nmedusasteps_ncoutfrequency) == 0) THEN
              l_ncout = .TRUE.
           ENDIF

           WRITE(*,*) 'MEDUSA OPTION ACTIVATED FOR DAY ', i
           WRITE(*,*) 'ayears_medusacurrent',  ayears_medusacurrent
           WRITE(*,*) '& ', ayears_medusacurrent + ayears_medusastep


           CALL  medusa_wrap_step(i_medusastep, ayears_medusacurrent,
     &          ayears_medusastep, l_ncout)
                                    ! From MEDUSA_WRAPPER
                                    !  - argument values:
                                    !     + istep = i_medusastep
                                    !     + atime0 = ayears_medusacurrent
                                    !     + datime = ayears_medusastep
                                    !     + l_write_nc: NetCDF output request
                                    !  - calls OCEAN_TO_SEDIMENT (from MOD_ILOVECLIM_O2S)
                                    !    (transfers data from MOD_ILOVECLIM_O2S
                                    !     to MOD_SEAFLOOR_CENTRAL)
                                    !  - checks outcome from OCEAN_TO_SEDIMENT
                                    !  - calls SOLVSED_ONESTEP to solve diagenesis eqns
                                    !  - calls REACLAY_X_CORELAY to regularize core layers
                                    !  - calls WRITERES_NCFILES_MEDUSA
                                    !    if l_write_nc==.TRUE.
                                    !  - calls SEDIMENT_TO_OCEAN
                                    !  - calls xchange_fluxes_s2o
                                    ! Update Medusa's time counter [yr]
          ayears_medusacurrent = ayears_medusaini + DBLE(i)/360.0D+00

          IF (i < ntotday) THEN     ! If we have not yet reached the last day,
                                    ! limit the length of the next Medusa step
                                    ! to the maximum number of days left.
            ndays_medusastep = MIN((ntotday - i), ndays_medusastep)
            ayears_medusastep = DBLE(ndays_medusastep)/360.0D+00
          ENDIF

          call update_flx_from_sediments() ! update the values of the fluxes towards CLIO

          WRITE(*,*) 'MEDUSA CALC FINISHED DAY ',i,ayears_medusacurrent

         end if

c~ #endif /* LONG_SED_RUN*/
#endif /* MEDUSA */


!     ***      << MEDUSA

!END OF  MAIN LOOP new version 2016-10-25 !mohr

!     ***      << MEDUSA


#if ( F_PALAEO >= 1 )
         CALL palaeo_timer(0,(int(irunlabel+(i/360)))) ! current time for timer is based on irunlabel
#endif

!nb modify bathymetry
#if ( BATHY == 3 )
         if (mod(i,reload_bathy).EQ.0) then
           write(*,*) 'call change bathy'
           CALL change_grid
         endif
#endif

#if ( F_PALAEO == 1 )
         if (mod(i,reload_topo).EQ.0) then
           CALL load_topo_masq(.FALSE.)
           CALL ec_masq(.FALSE.)
           CALL ec_topo
         endif
#endif
#if ( F_PALAEO == 3 )
         if (mod(i,reload_topo).EQ.0) then
            call subgrid_ecb_wrapper
            call ec_masq(.FALSE.)
            call ec_topo
         endif
#endif

!vm#if ( KC14 == 1 )
!vm      IF ( KENDY.EQ.1) THEN
!vm        CALL C14ATM_DP
!vm      ENDIF
!vm#endif

#if ( CLIO_OUT_NEWGEN == 1 )
       call daily_io_nc()
       if (mod(i,360*nbyearsinfile).eq.0) then
          CALL grid_io_reinit(i/(360*nbyearsinfile))
       endif
#endif


!     dmr --- Ecriture de l etat de redemarrage si necessaire.
         result_function_call = ec_writestate(i,ntotday)

      enddo                     ! boucle sur i, 1=>ntotday

!-----|--1--------2---------3---------4---------5---------6---------7-|
!     A partir d'ici la simulation est finie !!!
!     Fermeture et récupération des ordures (garbage_collection)
!-----|--1--------2---------3---------4---------5---------6---------7-|


!-----|--1--------2---------3---------4---------5---------6---------7-|
!     Fermeture du fichier 'C_reservoirs.txt'
!-----|--1--------2---------3---------4---------5---------6---------7-|

#if ( CYCC >= 2 )
      CALL out_cycc(-2,fractn(1,1,nld),darea)
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!     Fermeture de MEDUSA
!-----|--1--------2---------3---------4---------5---------6---------7-|

#if (MEDUSA == 1)
      CALL medusa_wrap_close()
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!     Fermeture de CARAIB
!-----|--1--------2---------3---------4---------5---------6---------7-|

#if ( CARAIB > 1 )
      call end_run()
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       FINALIZING THE CLIMATE INDICES
!-----|--1--------2---------3---------4---------5---------6---------7-|

#if ( CLM_INDICES >= 1 )
      CALL GLOBAL_FINALIZE
#endif


!-----|--1--------2---------3---------4---------5---------6---------7-|
!     Finalize the screen printing
!-----|--1--------2---------3---------4---------5---------6---------7-|

      call pretty_print_exec_time("End Execution")

! --- BdB 05-2019: deallocate array of output years
      deallocate(time_in_years)

      call ec_error(999)
      end program emic

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      The End of All Things (op. cit.)
!-----|--1--------2---------3---------4---------5---------6---------7-|
