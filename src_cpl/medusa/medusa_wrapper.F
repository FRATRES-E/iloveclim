#ifdef FN_THISFILE
#undef FN_THISFILE
#endif
#define FN_THISFILE "medusa_wrapper.F"
#ifndef __LINE__
#define __LINE__ 0
#endif

#define TEST 0
!-----|--1---------2---------3---------4---------5---------6---------7-|
!
!-----|--1---------2---------3---------4---------5---------6---------7-|

      MODULE MEDUSA_WRAPPER

!-----|--1---------2---------3---------4---------5---------6---------7-|
!      Ce module contient les définition pour initialiser et appeler
!       MEDUSA en mode couplé avec iLOVECLIM.
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche & Guy Munhoven
!      Date   : 05 janvier 2014
!      Modified by: Santiago Moreira, Didier M. Roche & Guy Munhoven
!      Derniere modification : 31 janvier 2022 (Guy Munhoven)
!-----|--1---------2---------3---------4---------5---------6---------7-|

!-----|--1---------2---------3---------4---------5---------6---------7-|

      USE mod_medinterfaces

!-----|--1--------2---------3---------4---------5---------6---------7-|


      CHARACTER(LEN=*), PARAMETER ::
     &  fn_thisfile     = FN_THISFILE
      CHARACTER(LEN=*), PARAMETER ::
     &  fmt_info_a      = '("['//fn_thisfile//']: ", A)' ,
     &  fmt_info_ai     = '("['//fn_thisfile//']: ", A, I0)',
     &  fmt_infolin_ia  = '("['//fn_thisfile//':", I0, "]: ", A)',
     &  fmt_infolin_iai = '("['//fn_thisfile//':", I0, "]: ", A, I0)'


                                    ! Common variables
      INTEGER, SAVE    :: n_columns = -1

                                    ! Time control
      DOUBLE PRECISION :: ayears_medusaini
      INTEGER          :: ndays_medusastep
      INTEGER          :: nmedusasteps_ncoutfrequency

      DOUBLE PRECISION :: ayears_medusacurrent, ayears_medusastep
      INTEGER          :: i_medusastep = 0
      INTEGER          :: medusa_timer = 0


      PRIVATE :: INIT_TIMECONTROL_MEDUSA

      ! [GM 31JAN2022]
      ! Auxiliary arrays for mass balance calculations
      !  - allocated in medusa_wrap_init()
      !  - used in medusa_wrap_prep() and medusa_wrap_step(...)
      !  - deallocated in medusa_wrap_close()
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:), SAVE
     &                 :: total_xm_0, total_xm_1,
     &                    topflx_xm_01, botflx_xm_01, total_rea_xm_01


      CONTAINS


!-----|--1--------2---------3---------4---------5---------6---------7-|
      SUBROUTINE MEDUSA_WRAP_INIT
!-----|--1--------2---------3---------4---------5---------6---------7-|
!     Cette subroutine permet d'initialiser MEDUSA. Appelee une fois
!     en debut de simulation. 
!     (dans l'environnement logiciel LUDUS)
!
!     C'est une reprise du 00Coupling_to_LVCLLSCE.howto de Guy Munhoven
!     reformatté en SUBROUTINE et dans un MODULE. 
!
!     Auteur : Didier M. Roche & Guy Munhoven
!     Modifie par: Santiago Moreira
!     Date   : 05 janvier 2014
!     Derniere modification : 31 janvier 2022 (Guy Munhoven)
!    
!-----|--1---------2---------3---------4---------5---------6---------7-|

      USE mod_defines_medusa
      USE mod_basicdata_medusa
      USE mod_gridparam,            ONLY: idnt, idnw
      USE mod_indexparam
      USE mod_seafloor_central,     ONLY: seafloor_setup
      USE mod_iloveclim_o2s,        ONLY: setup_o2s
      USE mod_iloveclim_s2o,        ONLY: setup_s2o, ini_fluxes_s2o
!      USE mod_iloveclim_setup,      ONLY:
!     &                                setup_iloveclim_medusa_xchange,
!     &                                iloveclim_nix, iloveclim_njy
      USE mod_iloveclim_setup

      USE mod_files_medusa

      USE mod_sedcore,              ONLY: SETUP_SEDCORE_SYSTEM
      USE mod_equilibcontrol,       ONLY: InitEquilibParameters
      USE mod_processcontrol,       ONLY: InitProcessParameters


      IMPLICIT NONE


      INTEGER            :: nix_cnpoh, njy_cnpoh

      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE
     &  :: dcnpoh_c_host
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE
     &  :: dcnpoh_n_host
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE
     &  :: dcnpoh_p_host
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE
     &  :: dcnpoh_o_host
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE
     &  :: dcnpoh_h_host
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE
     &  :: dcnpoh_remin_o2_host


      INTEGER :: iflag

! End generic definitions

      WRITE(*,*) ' MEDUSA_WRAP_INIT: EXECUTING! '
! 1. Open the log and diagnostic files

      OPEN(jp_stderr, FILE="medusa.err")
      OPEN(jp_stdlog, FILE="medusa.log")

c~ #ifdef DEBUG
c~       IF (jp_stddbg /= jp_stderr) OPEN(jp_stddbg, FILE="medusa.dbg")
c~       WRITE(jp_stddbg, fmt_infolin_ia) (__LINE__), 'Starting'
c~       CALL DEBUG_MEDUSA_MAIN_INDICES
c~ #endif

                                    ! Initialize list of requested files
                                    ! in MOD_FILES_MEDUSA
      CALL INIT_FILELIST_MEDUSA()
                                    ! Initialize the time control for
                                    ! Medusa (this MODULE):
                                    !  - ayears_medusaini: initial time [yr]
                                    !    (may be overriden later by
                                    !    InitSeafloorFromNetCDFFiles called
                                    !    in MEDUSA_WRAP_PREP below)
                                    !  - ndays_medusastep: Medusa time-step
                                    !    length [day]
                                    !  - nmedusasteps_ncoutfrequency: output
                                    !    frequency (time interval between
                                    !    records, in numbers of Medusa
                                    !    time steps)
      CALL INIT_TIMECONTROL_MEDUSA()



! 2. Complete configuration information in MOD_ILOVECLIM_SETUP
!    See details of the setup process in fait_pointer_MED_OCN.F
      CALL setup_iloveclim_medusa_xchange ! from mod_iloveclim_setup

      IF (nomcompo > 0) THEN
        nix_cnpoh = iloveclim_nix   ! Effective number of seafloor grid
        njy_cnpoh = iloveclim_njy   ! elements to consider for the
                                    ! stoechiometric data of OrgMatter_CNP
                                    ! components:
                                    !  - one MAY use
                                    !    nix_cnpoh = 1 AND njy_cnpoh = 1
                                    !    if every OrgMatter_CNP component
                                    !    has a uniform composition for all
                                    !    seafloor grid elements;
                                    !  - one MUST use
                                    !    nix_cnpoh = nix_host AND
                                    !    njy_cnpoh = nix_host
                                    !    if one or more OrgMatter_CNP
                                    !    components have a non-uniform
                                    !    composition!

#ifdef DEBUG
        print *, "x and y", nix_cnpoh, njy_cnpoh
#endif

        ALLOCATE(dcnpoh_c_host(nomcompo, nix_cnpoh, njy_cnpoh))
        ALLOCATE(dcnpoh_n_host(nomcompo, nix_cnpoh, njy_cnpoh))
        ALLOCATE(dcnpoh_p_host(nomcompo, nix_cnpoh, njy_cnpoh))
        ALLOCATE(dcnpoh_o_host(nomcompo, nix_cnpoh, njy_cnpoh))
        ALLOCATE(dcnpoh_h_host(nomcompo, nix_cnpoh, njy_cnpoh))
        ALLOCATE(dcnpoh_remin_o2_host(nomcompo, nix_cnpoh, njy_cnpoh))


        dcnpoh_c_host(ioo_om, 1:nix_cnpoh, 1:njy_cnpoh) 
     &    =  iloveclim_om_c(:,:)
        dcnpoh_n_host(ioo_om, 1:nix_cnpoh, 1:njy_cnpoh)
     &    =  iloveclim_om_n(:,:)
        dcnpoh_p_host(ioo_om, 1:nix_cnpoh, 1:njy_cnpoh)
     &    =  iloveclim_om_p(:,:)
        dcnpoh_o_host(ioo_om, 1:nix_cnpoh, 1:njy_cnpoh)
     &    =  0.0D+00
        dcnpoh_h_host(ioo_om, 1:nix_cnpoh, 1:njy_cnpoh)
     &    =  0.0D+00

        dcnpoh_remin_o2_host(ioo_om, 1:nix_cnpoh, 1:njy_cnpoh)
     &    =  iloveclim_om_ro2(:,:)


                                    ! From MOD_SEAFLOOR_CENTRAL.F
        CALL SEAFLOOR_SETUP(iloveclim_ocngrid_mask,
     &                              iloveclim_ocngrid_lon,
     &                              iloveclim_ocngrid_lat,
     &                              iloveclim_ocngrid_surf,
     &                              dcnpoh_c_host, dcnpoh_n_host,
     &                              dcnpoh_p_host, dcnpoh_o_host,
     &                              dcnpoh_h_host, dcnpoh_remin_o2_host,
     &                              n_columns)

        DEALLOCATE(dcnpoh_c_host)
        DEALLOCATE(dcnpoh_n_host)
        DEALLOCATE(dcnpoh_p_host)
        DEALLOCATE(dcnpoh_o_host)
        DEALLOCATE(dcnpoh_h_host)
        DEALLOCATE(dcnpoh_remin_o2_host)

      ELSE
                                    ! If there are no OrgMatter_CNP
                                    ! components in the model, the
                                    ! dcnpoh* arguments may not be used.

                                    ! From MOD_SEAFLOOR_CENTRAL.F
        CALL SEAFLOOR_SETUP(iloveclim_ocngrid_mask,
     &                              iloveclim_ocngrid_lon,
     &                              iloveclim_ocngrid_lat,
     &                              iloveclim_ocngrid_surf,
     &                              n_columns = n_columns)

      ENDIF


#ifdef DEBUG
      PRINT *, 'The number of n_columns is : ', n_columns
#endif


! 3. Initialise equilibrium and process parameter values
      CALL InitEquilibParameters    ! From mod_equilibcontrol.F
      CALL InitProcessParameters    ! From mod_processcontrol.F

! 4. Set up MOD_SEDCORE

      CALL SETUP_SEDCORE_SYSTEM(cfn_ncin_sedcore, cfn_ncout_sedcore)


! 5. Set up MOD_ILOVECLIM_O2S
! 5.1. Allocate arrays
      CALL setup_o2s                ! From MOD_ILOVECLIM_O2S

! 6. Set up MOD_ILOVECLIM_S2O
! 6.1. Allocate arrays

      CALL setup_s2o                ! From MOD_ILOVECLIM_S2O

! 6.2. Provisionally clear arrays (set to zero)
!      (may be reset by call of InitSeafloorFromNetCDFFiles
!       in medusa_wrap_prep below)
      CALL ini_fluxes_s2o           ! From MOD_ILOVECLIM_S2O


! 7. Prepare work arrays for the mass balance calculations
      ! [GM 31JAN2022]
      ALLOCATE(total_xm_0(ncompo,n_columns))
      ALLOCATE(total_xm_1(ncompo,n_columns))
      ALLOCATE(topflx_xm_01(ncompo,n_columns))
      ALLOCATE(botflx_xm_01(ncompo,n_columns))
      ALLOCATE(total_rea_xm_01(ncompo,n_columns))


#ifdef DEBUG
      WRITE (*,*) ' MEDUSA_WRAP_INIT: EXECUTED! '
#endif
      WRITE (*,*) ' MEDUSA_WRAP_INIT: EXECUTED! '


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE MEDUSA_WRAP_INIT
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
      SUBROUTINE INIT_TIMECONTROL_MEDUSA()
!-----------------------------------------------------------------------

      USE mod_defines_medusa
      USE mod_execontrol_medusa, ONLY: ABORT_MEDUSA
      USE mod_logunits


      IMPLICIT NONE

      NAMELIST /nml_timecontrol/ ayears_medusaini,
     &                           ndays_medusastep,
     &                           nmedusasteps_ncoutfrequency

      !!INTEGER,          INTENT(OUT) :: kflag


      INTEGER :: iflag
      INTEGER :: medusa_tsi_unit


      ! Standard I/O related data
      ! -------------------------

      CHARACTER(LEN=*), PARAMETER :: c_fmtinf_a = 
     &  '("[[INIT_TIMECONTROL_MEDUSA]]: ", A)'

      CHARACTER(LEN=*), PARAMETER :: c_fmterr_a = 
     &  '("[[INIT_TIMECONTROL_MEDUSA]] error: ", A)'

#ifdef DEBUG
      CHARACTER(LEN=*), PARAMETER :: c_fmtdbg_a = 
     &  '("[[INIT_TIMECONTROL_MEDUSA]] debug: ", A)'
#endif      


      iflag = RESERVE_LOGUNIT(medusa_tsi_unit)

      IF (iflag /= 0) THEN
        WRITE(jp_stderr,c_fmterr_a)
     &     'No more free units available -- aborting!'
        CALL ABORT_MEDUSA()
#ifdef DEBUG

      ELSE

        WRITE(jp_stddbg,'()')
        WRITE(jp_stddbg,c_fmtdbg_a, ADVANCE='NO')
     &                              'Assigning MEDUSA_TSI unit '
        WRITE(jp_stddbg,'(I0)') medusa_tsi_unit
#endif

      ENDIF

                                    ! Set default values
      ayears_medusaini            = 0.0D+00
      ndays_medusastep            = -1
      nmedusasteps_ncoutfrequency = -1

      OPEN(medusa_tsi_unit, FILE='medusa.tsi')

!      !READ(medusa_tsi_unit, NAMELIST=nml_timecontrol)
      READ(medusa_tsi_unit, NML=nml_timecontrol)

      CLOSE(medusa_tsi_unit)

      iflag = FREE_LOGUNIT(medusa_tsi_unit)
      IF (iflag /= 0) THEN
        WRITE(jp_stderr,c_fmterr_a)
     &    'Unexpected error while freeing MEDUSA_TSI unit ' //
     &    ' -- continuing!'
      ELSE
          medusa_tsi_unit   = -1
      ENDIF


      IF (ndays_medusastep <= 0) THEN
        WRITE(jp_stderr, c_fmterr_a)
     &     '"ndays_medusastep" must be positive -- resetting to 360'
        ndays_medusastep = 360
      ENDIF

      IF (nmedusasteps_ncoutfrequency <= 0) THEN
        WRITE(jp_stderr, c_fmterr_a)
     &     '"nmedusasteps_ncoutfrequency" must be positive -- ' //
     &     'resetting to 1'
        nmedusasteps_ncoutfrequency = 1
      ENDIF


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE INIT_TIMECONTROL_MEDUSA
!-----------------------------------------------------------------------




!-----|--1---------2---------3---------4---------5---------6---------7-|
      SUBROUTINE MEDUSA_WRAP_PREP()
!-----|--1---------2---------3---------4---------5---------6---------7-|


      USE mod_defines_medusa,       ONLY: jp_stderr
      USE mod_execontrol_medusa,    ONLY: ABORT_MEDUSA
      USE mod_seafloor_init
      USE mod_files_medusa,         ONLY: OPEN_NCFILES_MEDUSA
      USE mod_files_medusa,         ONLY: cfn_nmlin_init,cfn_ncin_init

      USE mod_files_medusa,         ONLY: cfn_ncin_flx, ctitle_ncfiles

      USE mod_iloveclim_s2o,        ONLY: SEDIMENT_TO_OCEAN
      USE mod_iloveclim_o2s,        ONLY: ini_fluxes_o2s,
     &                                    OCEAN_TO_SEDIMENT

      ! [GM 31JAN2022]
      ! Added USE clause to perform vertical integration of
      ! sediment profiles in REACLAY
      USE mod_zone_mass_totals


      IMPLICIT NONE


                                    ! time in years, not necessarily zero!
                                    ! If Medusa is initialized from 
                                    !   a NetCDF file, the time of the last
                                    !   record of that time is retrieved
                                    !   and used to override the input value
                                    !   of atime0
                                    ! else
                                    !   the input atime0 is used as is.

      INTEGER :: iflag
      INTEGER :: i_request


! END OF DECLARATIONS

      WRITE (*,*) ' MEDUSA_WRAP_PREP: EXECUTING! '


                                    ! Initialise the sediment columns,
                                    ! either from the namelist, or from
                                    ! the last time record of a NetCDF
                                    ! REACLAY file previously generated.
                                    ! (includes the initialisation of
                                    ! sediment-to-ocean fluxes,
                                    ! either from scratch (0) or
                                    ! from a previously generated NetCDF
                                    ! FLX file, via <cfn_ncin_flx> arg.)

      IF (cfn_ncin_init /= "/dev/null") THEN
                                    ! From MOD_SEAFLOOR_INIT

c~         write(*,*) 'value of atime before...', ayears_medusaini
c~         write(*,*) 'using cfn_ncinit = ', trim(cfn_ncin_init)
c~         write(*,*) 'using cfn_flx = ', trim(cfn_ncin_flx)
        !!ctitle_ncfiles

        CALL InitSeafloorFromNetCDFFiles(cfn_reaclay = cfn_ncin_init,
     &                                   cfn_flx     = cfn_ncin_flx,
     &                                   atime = ayears_medusaini)

        write(*,*) 'value of atime after...', ayears_medusaini


                                    ! Update the fluxes in S2O if they
                                    ! have been read in from <cfn_ncin_flx>
        IF (cfn_ncin_flx /= "/dev/null") THEN
          CALL SEDIMENT_TO_OCEAN(iflag)
        ENDIF

      ELSEIF (cfn_nmlin_init /= "/dev/null") THEN
                                    ! From MOD_SEAFLOOR_INIT
        CALL InitSeafloorFromNamelistFile(cfn_nmlin_init)

      ELSE

        WRITE(jp_stderr,'()')
        WRITE(jp_stderr,'("Error: no valid initialisation file given:")')
        WRITE(jp_stderr,'(" - cfn_nmlin_init = """, A, """")')
     &                              TRIM(cfn_nmlin_init)
        WRITE(jp_stderr,'(" - cfn_ncin_init = """, A, """")')
     &                              TRIM(cfn_ncin_init)
        WRITE(jp_stderr,'("Aborting!")')
        CALL ABORT_MEDUSA()

      ENDIF


      CALL ini_fluxes_o2s           ! Not necessary, but no harm either
      CALL OCEAN_TO_SEDIMENT(i_request = 0, rbflag = iflag)



                                    ! Create the NetCDF files that have
                                    ! been requested and by giving atime0
                                    ! request to store the initial
                                    ! state of the model
      CALL OPEN_NCFILES_MEDUSA(ayears_medusaini)


c~       print *, "NETCDF FILES CREATED and initial state written"


                                    ! Initialisation complete, we can
                                    ! loop back to iLOVECLIM so that it
                                    ! perfoms his steps ...
                                    ! Run iLoveclim over one sediment time-step

c~       WRITE(*,*) ' MEDUSA_WRAP_PREP: EXECUTED!'
c~       WRITE(*,*) 'INFO: ', ctitle_ncfiles


      ! [GM 31JAN2022]
      ! Initial component contents in the REACLAY parts of the surface
      ! sediment (vertically integrated, for each component, for each column)
      CALL ZONE_MASS_TOTALS(total_xm = total_xm_0)



      RETURN


!-----|--1---------2---------3---------4---------5---------6---------7-|
      END SUBROUTINE MEDUSA_WRAP_PREP
!-----|--1---------2---------3---------4---------5---------6---------7-|




!-----|--1---------2---------3---------4---------5---------6---------7-|
      SUBROUTINE MEDUSA_WRAP_STEP(istep, atime0, datime, l_write_nc)
!-----|--1---------2---------3---------4---------5---------6---------7-|
!      Cette subroutine permet d'effectuer un pas de temps MEDUSA. 
!       (dans l'environnement logiciel LUDUS)
!
!      C'est une reprise du 00Coupling_to_LVCLLSCE.howto de Guy Munhoven
!       reformatté en SUBROUTINE et dans un MODULE. 
!
!      Auteur : Didier M. Roche & Guy Munhoven
!      Date   : 05 janvier 2014
!      Derniere modification : 31 janvier 2022 (Guy Munhoven)
!-----|--1--------2---------3---------4---------5---------6---------7-|

      USE mod_defines_medusa,       ONLY: jp_stderr
      USE mod_execontrol_medusa,    ONLY: ABORT_MEDUSA
      USE mod_iloveclim_o2s,        ONLY: OCEAN_TO_SEDIMENT
      USE mod_iloveclim_s2o,        ONLY: SEDIMENT_TO_OCEAN
      USE mod_files_medusa,         ONLY: WRITERES_NCFILES_MEDUSA
      USE mod_sedcore,              ONLY: REACLAY_X_CORELAY

      ! [GM 31JAN2022]
      ! Added USE clauses to perform vertical integration of
      ! sediment profiles in REACLAY and to do the mass balance
      ! calculations. MOD_MASS_BALANCES also provides the
      ! TYPE(MASSBAL_TERMS) definition.
      USE mod_indexparam,           ONLY: jc_to_io, jf_to_io
     &                                  , nsolut, nsolid
      USE mod_seafloor_central,     ONLY: GET_BOUNDARY_CONDS,
     &                                    GET_BOUNDARY_FLUXES
      USE mod_zone_mass_totals
      USE mod_mass_balances


      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: atime0
      INTEGER,          INTENT(IN)  :: istep
      DOUBLE PRECISION, INTENT(IN)  :: datime
      LOGICAL, INTENT(IN), OPTIONAL :: l_write_nc


      INTEGER :: i_request, n_trouble, i_flag
      LOGICAL :: lloc_write_nc
      DOUBLE PRECISION              :: atime


      ! [GM 31JAN2022]
      ! Auxiliary arrays and variables required for mass balance calculations
      DOUBLE PRECISION, DIMENSION(nsolut) :: wcflx
      DOUBLE PRECISION, DIMENSION(nsolid) :: wfflx, bfflx
      TYPE(MASSBAL_TERMS) :: g_bulk_0, g_bulk_1
      TYPE(MASSBAL_TERMS) :: g_bulk_topflx, g_bulk_botflx, g_bulk_rea
      TYPE(MASSBAL_TERMS) :: g_orgm_0, g_orgm_1
      TYPE(MASSBAL_TERMS) :: g_orgm_topflx, g_orgm_botflx, g_orgm_rea
      INTEGER :: i_column


! END GENERIC DEFINITIONS


      WRITE(*,*) 'MEDUSA_WRAP_STEP: EXECUTING in step: ', istep 

#if (TEST > 0)
      print *, "i_medusastep.............", i_medusastep
      print *, "ayears_medusacurrent.....", ayears_medusacurrent
      print *, "ndays_medusastep.........", ndays_medusastep
      print *, "ayears_medusastep........", ayears_medusastep
#endif

      IF (PRESENT(l_write_nc)) THEN
        lloc_write_nc =  l_write_nc
      ELSE
        lloc_write_nc = .FALSE.
      ENDIF

      i_flag = 0

                                    ! Before doing calculations, let's
                                    ! update the grids in MEDUSA using
                                    ! the corresponding ones in the
                                    ! ocean model... 
! !mohr
      PRINT*, "value of istep", istep
      PRINT*, "value of irequest", i_request
      PRINT*, "value of i_flag", i_flag

                                    ! the following first part has to be
                                    ! done in the init part, ensuring
                                    ! that the clio fluxes are provided
                                    ! to o2s first!!
![??? GM] The model must work with i_request=0 (i.e., must not give i_flag=1)
![??? GM] else, something is not correctly done.
      i_request = 0 ! Mandatorily read new dataset of boundary conditions
      i_request = 1

![??? GM] reading the next two blocks of comments, I think
![??? GM] that something must have gone wrong somewhere
![??? GM] i_flag = 1 after OCEAN_TO_SEDIMENT completes is
![??? GM] actually a fatal error!
c~       print *, "i_request, i_flag", i_request, i_flag
c~       print *, "n_datasets_read", n_datasets_read
c~       print *, " mod_iloveclim_o2s_setup",  mod_iloveclim_o2s_setup
c~       print *, "iloveclim_o2s_id",iloveclim_o2s_id

      CALL OCEAN_TO_SEDIMENT(i_request = i_request, rbflag = i_flag)

c~       print *, "Nope, I was wrong, after", i_request, i_flag
c~       print *, "CAREFUL!WHY i_flag becomes 1 after OCEAN_TO_SEDIMENT?"
c~ !      !i_flag = 0 !I had to force this value... !mohr [XXX]
      
      SELECT CASE(i_flag)
        CASE(0)
            CONTINUE
        CASE(-1)
          IF(i_request /= 0) THEN
            WRITE(jp_stderr,fmt_infolin_ia) (__LINE__),
     &        'Warning: iterating with boundary conditions from '//
     &        'previous time step.'
          ELSE
            WRITE(jp_stderr,fmt_infolin_ia) (__LINE__),
     &        'Missing mandatory boundary conditions in '//
     &        'ocean_to_sediment. Aborting!'
            CALL ABORT_MEDUSA()
          ENDIF
        CASE(1)
          WRITE(jp_stderr,fmt_infolin_ia) (__LINE__),
     &        'Fatal error (code 1) in ocean_to_sediment. Aborting!'
          CALL ABORT_MEDUSA()
        CASE DEFAULT
          WRITE(jp_stderr,fmt_infolin_ia) (__LINE__),
     &        'Unknown error in ocean_to_sediment. Aborting!'
          CALL ABORT_MEDUSA()
      END SELECT

c~       print *,"flag 001 m_step, before solvsed_onestep"

c~       print *,    "atime0",     atime0
c~       print *,    "datime",     datime
c~       print *,    "n_columns",  n_columns
c~       print *,    "i_flag",     i_flag
c~       print *,    "n_trouble",  n_trouble

      i_flag = 0

      CALL SOLVSED_ONESTEP(atime0, datime, n_columns,
     &                       i_flag, n_trouble)

      atime = atime0 + datime

c~       print *,"flag 002 m_step", atime
      
c~       print *,    "atime0",     atime0
c~       print *,    "datime",     datime
c~       print *,    "n_columns",  n_columns
c~       print *,    "i_flag",     i_flag
c~       print *,    "n_trouble",  n_trouble


                                    ! [GM 31JAN2022]
                                    ! Retrieve global mass balance terms:
                                    !
                                    ! - vertical integration for the end
                                    !   of the time step, both contents
                                    !   reaction rates over the time step
                                    !   (NB: reaction rates and fluxes are
                                    !   for the entire tiem step, but due to
                                    !   the implicit time approach, registered
                                    !   for the end of the time step).
      CALL ZONE_MASS_TOTALS(total_xm = total_xm_1, 
     &                              total_rea_xm = total_rea_xm_01)
      total_rea_xm_01(:,:) = total_rea_xm_01(:,:) !*datime

                                    ! - retrieve the fluxes accross the
                                    !   top and bottom boundaries of REACLAY
      DO i_column = 1, n_columns
        CALL GET_BOUNDARY_CONDS(i_column = i_column, gbcflag = i_flag,
     &                              wfflx = wfflx)
        topflx_xm_01(jf_to_io(:),i_column) = wfflx(:) !*datime
          
        CALL GET_BOUNDARY_FLUXES(i_column = i_column, gbfflag = i_flag,
     &                              wcflx = wcflx, bfflx = bfflx)
        topflx_xm_01(jc_to_io(:),i_column) = wcflx(:) !*datime
        botflx_xm_01(jf_to_io(:),i_column) = bfflx(:) !*datime
        botflx_xm_01(jc_to_io(:),i_column) = 0.0D+00
      ENDDO

      IF (datime > 0.0D+00) THEN
           topflx_xm_01(:,:)    = topflx_xm_01(:,:)*datime
           botflx_xm_01(:,:)    = botflx_xm_01(:,:)*datime
           total_rea_xm_01(:,:) = total_rea_xm_01(:,:)*datime
      ELSE
           total_xm_0=total_xm_1
      ENDIF


                                    ! - produce mass balances for the
                                    !   required conservation properties
                                    !   (C, Alk,, O2, PO4), integrate them
                                    !   globally and analyse them
      CALL MASS_BALANCES_GLOBAL(    total_xm_0, total_xm_1,
     &                              topflx_xm_01, botflx_xm_01,
     &                              total_rea_xm_01, n_columns,
     &                              g_bulk_0, g_bulk_1, g_bulk_topflx,
     &                              g_bulk_botflx, g_bulk_rea,
     &                              g_orgm_0, g_orgm_1, g_orgm_topflx,
     &                              g_orgm_botflx, g_orgm_rea,
     &                              l_write2log=.TRUE.)


       WRITE(*,'(A,15(1X,D12.5))') "CTRL_MASS_BALANCES_BULK", 
     &   g_bulk_topflx%c/1.E12, g_bulk_topflx%alk/1.E12,
     &   g_bulk_topflx%po4/1.E12,
     &   g_bulk_botflx%c/1.E12, g_bulk_botflx%alk/1.E12,
     &   g_bulk_botflx%po4/1.E12,
     &   g_bulk_rea%c/1.E12, g_bulk_rea%alk/1.E12,
     &   g_bulk_rea%po4/1.E12,
     &   g_bulk_0%c/1.E12, g_bulk_0%alk/1.E12, g_bulk_0%po4/1.E12,
     &   g_bulk_1%c/1.E12, g_bulk_1%alk/1.E12, g_bulk_1%po4/1.E12
       WRITE(*,'(A,15(1X,D12.5))') "CTRL_MASS_BALANCES_ORGANIC", 
     &   g_orgm_topflx%c/1.E12, g_orgm_topflx%alk/1.E12,
     &   g_orgm_topflx%po4/1.E12,
     &   g_orgm_botflx%c/1.E12, g_orgm_botflx%alk/1.E12,
     &   g_orgm_botflx%po4/1.E12,
     &   g_orgm_rea%c/1.E12, g_orgm_rea%alk/1.E12,
     &   g_orgm_rea%po4/1.E12,
     &   g_orgm_0%c/1.E12, g_orgm_0%alk/1.E12, g_orgm_0%po4/1.E12,
     &   g_orgm_1%c/1.E12, g_orgm_1%alk/1.E12, g_orgm_1%po4/1.E12

      ! Upon completion of the previous call, the following mass balance
      ! terms can be retrieved (units are moles or equivalents for stocks,
      ! and moles or equivalents per time-step for the fluxes):
      !  - g_bulk_0%alk : global alkalinity stock in the surface sediment
      !                   at the beginning of the Medusa time step
      !                   just completed, porewater + solids,
      !                   organically + inorganically carried
      !  - g_bulk_0%c   : idem for carbon (organic and inorganic)
      !  - g_bulk_0%o2  : idem for O2, incl. negative O2 content of OrgM
      !  - g_bulk_0%po4 : idem for phosphate (solids only)
      !
      !  - g_bulk_1%xyz : same as above, at the end of the time step
      !
      !  - g_bulk_topflx%xyz : top fluxes for the above, over the
      !                   whole time step
      !
      !  - g_bulk_botflx%xyz : bottom fluxes for the above, over the
      !                   whole time step
      !
      !  - g_bulk_rea%xyz : total reaction rate for the above, over the
      !                   whole time step
      !
      !  - g_orgm_0%xyz
      !    g_orgm_1%xyz
      !    g_orgm_botflx%xyz : same as above, but carried by for each
      !                   of the xyz carried by organic matter only.
      !
      ! Information can now be extracted from the ten structures and
      ! transferred elsewhere if required



                                    !  - and finally prepare the mass balance
                                    !    calculations of the next step.
      total_xm_0(:,:) = total_xm_1(:,:)

      ! [GM 31JAN2022] End of additions for mass balance calculations



      CALL REACLAY_X_CORELAY(atime)

c~       print *,"flag 003 writting NetCDF files?", lloc_write_nc
      IF (lloc_write_nc) CALL WRITERES_NCFILES_MEDUSA(atime)


      CALL SEDIMENT_TO_OCEAN(i_flag)


      WRITE(*,*) 'MEDUSA_WRAP_STEP:  timestep: ', atime0, atime


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE MEDUSA_WRAP_STEP
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
      SUBROUTINE MEDUSA_WRAP_CLOSE
!-----------------------------------------------------------------------


      USE mod_files_medusa,         ONLY: CLOSE_RESFILES_MEDUSA
      USE mod_defines_medusa,       ONLY: jp_stderr, jp_stdlog
#ifdef DEBUG
      USE mod_defines_medusa,       ONLY: jp_stddbg
#endif


      IMPLICIT NONE


#ifdef DEBUG
      WRITE(jp_stddbg,fmt_infolin_ia) 
     &   (__LINE__), 'Starting MEDUSA_WRAP_CLOSE'
#endif
      WRITE(*,*) ' MEDUSA_WRAP_CLOSE: EXECUTING. END OF MEDUSA' 


      CALL CLOSE_RESFILES_MEDUSA

! Deallocate work arrays for the mass balance calculations
      ! [GM 31JAN2022]
      DEALLOCATE(total_xm_0)
      DEALLOCATE(total_xm_1)
      DEALLOCATE(topflx_xm_01)
      DEALLOCATE(botflx_xm_01)
      DEALLOCATE(total_rea_xm_01)


#ifdef DEBUG
!-----Close debug file if still open------------------------------------
      IF (jp_stddbg /= jp_stderr) CLOSE(jp_stddbg)
      WRITE(jp_stddbg, fmt_infolin_ia) 
     &   (__LINE__), 'Ending MEDUSA_WRAP_CLOSE'
#endif
      CLOSE(jp_stderr)
      CLOSE(jp_stdlog)

      WRITE(*,*) ' MEDUSA_WRAP_CLOSE: EXECUTED. END OF MEDUSA' 


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE MEDUSA_WRAP_CLOSE
!-----------------------------------------------------------------------

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      The End of All Things (op. cit.)
!-----|--1--------2---------3---------4---------5---------6---------7-|
      END MODULE MEDUSA_WRAPPER
!-----|--1--------2---------3---------4---------5---------6---------7-|
