!=======================================================================
      MODULE MOD_ILOVECLIM_O2S
     
!=======================================================================
!!!!!!

      IMPLICIT NONE


      INTEGER, SAVE, PRIVATE :: n_grid_seafloor   = -1
      INTEGER, SAVE, PRIVATE :: nix               = -1
      INTEGER, SAVE, PRIVATE :: njy               = -1

                                    ! Number of boundary condition sets
                                    ! previously read in:
                                    ! - not even initialised: -1
      INTEGER, SAVE :: iloveclim_o2s_id = -1


                                    ! Has this module been setup
                                    ! (arrays allocated) ?
      LOGICAL, SAVE, PRIVATE :: mod_iloveclim_o2s_setup = .FALSE.


      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_temp
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_sali
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_dbsl


      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_co2
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_co3
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_hco3
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_oxyg


      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calc
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_clay


      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_tracer01 !mohr


#ifdef WITH_C13
                                    ! C-13 isotope related arrays
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_oc13

      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgmC13
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calcC13
#endif

#ifdef WITH_C14
                                    ! C-14 isotope related arrays
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_oc14

      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgmC14
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calcC14
#endif

#ifdef WITH_O18
                                    ! O-18 isotope related arrays
                                    ! (for usage of O-18 as a passive
                                    ! tracer in calcite shells)
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calcO18
#endif


      CONTAINS


!-----------------------------------------------------------------------
      SUBROUTINE setup_o2s
!-----------------------------------------------------------------------
!-----|--1---------2---------3---------4---------5---------6---------7-|
!     Cette subroutine contient l'allocation des tableaux d'échange
!     iLOVECLIM vers MEDUSA au moment de l'initialisation. Elle est
!     appelée une fois au début de simulation
!     (dans l'environnement logiciel LUDUS)
!
!     Auteur : Didier M. Roche & Guy Munhoven
!     Date   : 05 fevrier 2014
!     Dernière modification : 21 octobre 2016 (Guy Munhoven)
!-----|--1---------2---------3---------4---------5---------6---------7-|

      USE mod_defines_medusa,       ONLY: jp_stderr
      USE mod_execontrol_medusa,    ONLY: ABORT_MEDUSA
      USE mod_seafloor_central,     ONLY: IJ_COLUMNS_USED,
     &                                    N_COLUMNS_USED
!      USE mod_seafloor_central

      IMPLICIT NONE

      CHARACTER(LEN=*), PARAMETER :: c_fmterr_a =
     &  '("[SETUP_ILOVECLIM_O2S/setup_o2s] error: ", A)'


      IF (.NOT. mod_iloveclim_o2s_setup) THEN

        CALL IJ_COLUMNS_USED(nix, njy)
        CALL N_COLUMNS_USED(n_grid_seafloor)

        IF (n_grid_seafloor /= -1) THEN

          ALLOCATE(iloveclim_temp(nix,njy))
          ALLOCATE(iloveclim_sali(nix,njy))
          ALLOCATE(iloveclim_dbsl(nix,njy))

          ALLOCATE(iloveclim_co2(nix,njy))
          ALLOCATE(iloveclim_co3(nix,njy))
          ALLOCATE(iloveclim_hco3(nix,njy))
          ALLOCATE(iloveclim_oxyg(nix,njy))

          ALLOCATE(iloveclim_orgm(nix,njy))
          ALLOCATE(iloveclim_calc(nix,njy))
          ALLOCATE(iloveclim_clay(nix,njy))

          ALLOCATE(iloveclim_tracer01(nix,njy))

#ifdef WITH_C13
          ALLOCATE(iloveclim_oc13(nix,njy))

          ALLOCATE(iloveclim_orgmC13(nix,njy))
          ALLOCATE(iloveclim_calcC13(nix,njy))
#endif

#ifdef WITH_C14
          ALLOCATE(iloveclim_oc14(nix,njy))
          ALLOCATE(iloveclim_orgmC14(nix,njy)) 
          ALLOCATE(iloveclim_calcC14(nix,njy))
#endif

#ifdef WITH_O18
          ALLOCATE(iloveclim_calcO18(nix,njy))
#endif

          mod_iloveclim_o2s_setup = .TRUE.

        ELSE

          WRITE(jp_stderr,c_fmterr_a)
     &      'MOD_SEAFLOOR_CENTRAL not yet set up -- aborting!'
          CALL ABORT_MEDUSA()

        ENDIF

      ELSE

        WRITE(jp_stderr,c_fmterr_a)
     &         'MOD_ILOVECLIM_O2S already set up -- aborting!'
        CALL ABORT_MEDUSA()

      ENDIF


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE setup_o2s
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
      SUBROUTINE ini_fluxes_o2s
!-----------------------------------------------------------------------


      USE mod_defines_medusa,       ONLY: jp_stderr
      USE mod_execontrol_medusa,    ONLY: ABORT_MEDUSA

      USE mod_iloveclim_setup,      ONLY: iloveclim_ocngrid_dbsl

      IMPLICIT NONE


      CHARACTER(LEN=*), PARAMETER :: c_fmterr_a =
     &  '("[MOD_ILOVECLIM_O2S/ini_fluxes_s2o] error: ", A)'

      IF (.NOT. mod_iloveclim_o2s_setup) THEN
        WRITE(jp_stderr,c_fmterr_a)
     &    'MOD_ILOVECLIM_O2S not yet set up -- aborting!'
        CALL ABORT_MEDUSA()
      ENDIF


      iloveclim_temp(:,:) = 0.0D+00
      iloveclim_sali(:,:) = 0.0D+00
      iloveclim_dbsl(:,:) = iloveclim_ocngrid_dbsl(:,:)

      iloveclim_oxyg(:,:) = 0.0D+00
      iloveclim_co2(:,:)  = 0.0D+00
      iloveclim_co3(:,:)  = 0.0D+00
      iloveclim_hco3(:,:) = 0.0D+00

      iloveclim_clay(:,:) = 0.0D+00
      iloveclim_orgm(:,:) = 0.0D+00
      iloveclim_calc(:,:) = 0.0D+00

      iloveclim_tracer01(:,:) = 0.0D+00

#ifdef WITH_C13
      iloveclim_oc13(:,:)    = 0.0D+00

      iloveclim_orgmC13(:,:) = 0.0D+00
      iloveclim_calcC13(:,:) = 0.0D+00
#endif

#ifdef WITH_C14
      iloveclim_oc14(:,:)    = 0.0D+00

      iloveclim_orgmC14(:,:) = 0.0D+00
      iloveclim_calcC14(:,:) = 0.0D+00
#endif

#ifdef WITH_O18
      iloveclim_calcO18(:,:) = 0.0D+00
#endif


      iloveclim_o2s_id = iloveclim_o2s_id + 1

      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE ini_fluxes_o2s
!-----------------------------------------------------------------------

c~ #define VALUE_DEBUG 1

c~ #if ( VALUE_DEBUG == 1)
c~       SUBROUTINE safe_gard(tabel,min_v, max_v, label)
        
c~         DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: tabel
c~         DOUBLE PRECISION,                 INTENT(in) :: min_v, max_v
c~         CHARACTER(len=*),                 INTENT(in) :: label
        
c~         DOUBLE PRECISION :: tmin_v, tmax_v
      
c~         tmin_v = MINVAL(tabel)
c~         tmax_v = MAXVAL(tabel)
        
c~         if ((tmin_v.LT.min_v).OR.(tmax_v.GT.max_v)) then
c~           WRITE(*,*)
c~           WRITE(*,*) "==== TROUBLESOME VALUES IN OCEAN TO SEDIMENT ===="
c~           WRITE(*,*) "        ARRAY is "//TRIM(label)// "   !!!!"
c~           WRITE(*,*) "min", tmin_v, min_v
c~           WRITE(*,*) "max", tmax_v, max_v
c~         endif
        
c~       END SUBROUTINE safe_gard
c~ #endif



!-----------------------------------------------------------------------
       SUBROUTINE OCEAN_TO_SEDIMENT(i_request, rbflag)
!-----------------------------------------------------------------------

      ! Loads a new complete set of forcing arrays into
      ! MOD_ILOVECLIM_COUPLED, using data in MOD_ILOVECLIM_O2S.
      ! MOD_ILOVECLIM_O2S data are expected to be expressed in the same
      ! units as those adopted for iLoveclim. They are being converted
      ! here as required.

      USE mod_seafloor_wdata, ONLY: WDATA_CONTAINER
      USE mod_defines_medusa, ONLY: jp_stderr
      USE mod_execontrol_medusa, ONLY: ABORT_MEDUSA

      USE mod_indexparam
      USE mod_materialcharas

!      USE mod_seafloor_central

      USE mod_seafloor_central, ONLY:
     &                              COLUMN_N2IJ, SAVE_BOUNDARY_CONDS,
     &                              GET_MATERIALCHARAS, COLUMN_AREA4N

! =====================
! Variable declarations
! =====================

      IMPLICIT NONE

! Argument list variables
! -----------------------

      INTEGER, INTENT(IN)
     &   :: i_request               ! Which type of request ?
                                    ! possible values for <i_request>
                                    !  0: mandatorily read new data
                                    !  1: accept missing new data, proceed with
                                    !     currently active set
                                    ! -1: close open data files, and reset
                                    !     all counters to non-initialized

      INTEGER, INTENT(OUT)
     &   :: rbflag                  ! Return flag: possible values for <rbflag> reserved so far
                                    ! -1: no new data set read in
                                    !     (only if i_request /= 0)
                                    !  0: everything went OK
                                    !  1: mandatory read failed

! Local Variables
! ---------------


      INTEGER       :: i, j, n
      INTEGER, SAVE :: n_datasets_read     =  0
      INTEGER, SAVE :: iloveclim_o2s_previd = -1

      TYPE(WDATA_CONTAINER) :: wdata
!!      DOUBLE PRECISION :: wsalin, wdbsl, wtmpdc
      DOUBLE PRECISION, DIMENSION(nsolut) :: wconc
      DOUBLE PRECISION, DIMENSION(nsolid) :: wfflx

      INTEGER :: iflag

      CHARACTER(LEN=*), PARAMETER :: c_fmterr_a = 
     &  '("[MOD_ILOVECLIM_O2S/OCEAN_TO_SEDIMENT] error: ", A)'



! End of Declarations
! ===================
c~       print*, "SUBROUTINE OCEAN_TO_SEDIMENT, flag 001"
c~       print*, "IDs === v000 === ", iloveclim_o2s_previd, n_datasets_read
c~       print*, mod_iloveclim_o2s_setup, iloveclim_o2s_id
c~       print*, "i_request ==", i_request
      
      IF (n_datasets_read == 0) THEN

         IF (.NOT. mod_iloveclim_o2s_setup) THEN
            WRITE(jp_stderr,c_fmterr_a)
     &         'MOD_ILOVECLIM_O2S not yet set up - aborting'
            CALL ABORT_MEDUSA()
         ENDIF


         IF (iloveclim_o2s_id == -1) THEN
            WRITE(jp_stderr,c_fmterr_a)
     &         'Data not yet loaded into MOD_ILOVECLIM_O2S - aborting'
            CALL ABORT_MEDUSA()
         ENDIF

      ENDIF

c~ #ifdef DEBUG
c~       print*, "SUBROUTINE OCEAN_TO_SEDIMENT, flag 002"
c~       write(*,*) "IDs ===", iloveclim_o2s_id, iloveclim_o2s_previd
c~ #endif

      IF (iloveclim_o2s_id <= iloveclim_o2s_previd) THEN
                                    ! No new data available!
         IF(i_request == 0) THEN

            rbflag = 1
            RETURN

         ELSE
                                    ! No more *new* data! However, a complete
            rbflag = -1             ! set may previously have been stored.
                                    ! Use it, but issue warning!
         ENDIF

      ENDIF

c~ #ifdef DEBUG
c~       print*, "SUBROUTINE OCEAN_TO_SEDIMENT, flag 003"
c~ #endif


      IF (i_request == -1) THEN

        n_datasets_read = 0
        iloveclim_o2s_previd = -1
        iloveclim_o2s_id = -1

        rbflag = 0

        RETURN

      ENDIF

                                    ! OK, everything cleared now:
                                    ! proceed to loading the data
                                    ! into MOD_SEAFLOOR_CENTRAL.
c~ #ifdef DEBUG
!      print*, "ngrid_seafloor",n_grid_seafloor
!      print*, "SUBROUTINE OCEAN_TO_SEDIMENT, flag 004"
!      print*, "orgm_correc = 0.001, calc_correc = 0.002" 
c~ #endif

c~ #if ( VALUE_DEBUG == 1 )
c~       call safe_gard(iloveclim_temp(:,:),-5.0,40.0
c~      >            , "iloveclim_temp") ! [degC]
c~       call safe_gard(iloveclim_sali(:,:),0.0,45.0
c~      >            , "iloveclim_sali")  ! [-]
c~       call safe_gard(iloveclim_co3(:,:)*1.0D+03,0.0,0.3
c~      >            , "iloveclim_co3")
c~       call safe_gard(iloveclim_hco3(:,:)*1.0D+03,1.9,2.5
c~      >            , "iloveclim_hco3")
c~       call safe_gard(iloveclim_co2(:,:)*1.0D+03,0.0,0.12
c~      >            , "iloveclim_co2")        
c~       call safe_gard(iloveclim_oxyg(:,:)*1.0D-03,0.0,0.4
c~      >            , "iloveclim_oxyg")       
c~       call safe_gard(iloveclim_clay(:,:),0.0019,0.0021
c~      >            , "iloveclim_clay")       
c~       call safe_gard(iloveclim_calc(:,:)*mol_calc*1.E12
c~      > ,0.0,0.08  , "iloveclim_calc")       
c~       call safe_gard((iloveclim_orgm(:,:)/om_c)*mol_om*1.E12
c~      > ,0.0,0.035 , "iloveclim_orgm")
c~       WRITE(*,*) "mol_s values ... = ", mol_calc, mol_om, om_c
c~ #endif   


      DO n = 1, n_grid_seafloor

        CALL COLUMN_N2IJ(i_column = n, iflag = iflag, ix = i, jy = j)

        CALL GET_MATERIALCHARAS(i_column = n, iflag = iflag)

c~         CALL COLUMN_AREA4N(i_column = n, iflag = iflag
c~      &                   , sfc_area=surface_area)
     
        wdata%wtmpc   = iloveclim_temp(i,j) ! [degC]
        wdata%wsalin = iloveclim_sali(i,j)  ! [-]
        wdata%wdbsl  = iloveclim_dbsl(i,j)  ! [m]

       ! wtmpdc = iloveclim_temp(i,j)  ! [degC]
       ! wsalin = iloveclim_sali(i,j)  ! [-]
       ! wdbsl  = iloveclim_dbsl(i,j)  ! [m]


                                    ! Now store the information in the
                                    ! respective MOD_SEAFLOOR_CENTRAL arrays

                                    ! - solutes

                                    ! According to mbiota.F in Climber:
                                    ! units for PO4, NO3, O2, DOC : micromol/kg
                                    ! units for DIC and ALK: mol/kg

!       ! print*, "iloveclim_hco3", iloveclim_hco3
!       ! print*, "iloveclim_oxyg", iloveclim_oxyg

        wconc(ic_co3)  = iloveclim_co3(i,j)  * 1.0D+03 ! from  mol/kg to mol/m3
        wconc(ic_hco3) = iloveclim_hco3(i,j) * 1.0D+03 ! from  mol/kg to mol/m3
        wconc(ic_co2)  = iloveclim_co2(i,j)  * 1.0D+03 ! from  mol/kg to mol/m3

 


        wconc(ic_o2)   = iloveclim_oxyg(i,j) * 1.0D-03 ! from µmol/kg to mol/m3

c~ #ifdef DEBUG
c~         print *, "ISOTOPES in mod_iloveclim_o2s.F ", iloveclim_oc13(1,1)
c~ #endif
#ifdef WITH_C13
        wconc(ic_dic_c13)  = iloveclim_oc13(i,j) * 1.0D+03 ! [??? GM] units to be checked
#endif
#ifdef WITH_C14
        wconc(ic_dic_c14)  = iloveclim_oc14(i,j) * 1.0D+03 ! [??? GM] units to be checked
#endif



c~        print*, "mod_iloveclim_o2s.F , checking particle fluxes..."
c~        print*, "clay_ini =", iloveclim_clay(i,j)
![XXX] fluxes in iloveclim_XXX are /day, except for clay (/yr)
                                    ! - solids
                                    ! Clay: units are [kg/m2/yr] in iLoveclim,
                                    ! and [kg/m2/yr] in Medusa.
        wfflx(if_clay) = iloveclim_clay(i,j)

                                    ! Calcite : units are [Tmol C/m2/yr] in iLoveclim,
                                    ! to be converted to [kgC/m2/yr] here
                                    ! molC -> kg calc => * mol_calc
                                    ! Tkg -> kg => *1E12
        wfflx(if_calc) = iloveclim_calc(i,j) *mol_calc * 1.E12
     
                                    ! POC : units are [Tmol C/m2/yr] in iLoveclim,
                                    ! to be converted to [kg OrgMatter/m2/yr]
                                    ! molC -> mol.m.o => / om_c
                                    ! mol.m.o -> kg o.m. => * mol_om
                                    ! Tkg -> kg => *1E12
                                    ! here (not [kg C/m2/yr]!)
        wfflx(if_om)   = (iloveclim_orgm(i,j)/om_c) * mol_om * 1.E12


#ifdef WITH_C13
        wfflx(if_om_c13)   = (iloveclim_orgmC13(i,j)/om_c) * mol_om ! [??? GM] units to be checked
        wfflx(if_calc_c13) = (iloveclim_calcC13(i,j) * mol_calc     ! [??? GM] units to be checked
#endif

#ifdef WITH_C14
        wfflx(if_om_c14)   = (iloveclim_orgmC14(i,j)/om_c) * mol_om ! [??? GM] units to be checked
        wfflx(if_calc_c14) = (iloveclim_calcC14(i,j) * mol_calc     ! [??? GM] units to be checked
#endif

#ifdef WITH_O18
        wfflx(if_calc_o18) = (iloveclim_calcO18(i,j) * mol_calc     ! [??? GM] units to be checked
#endif

!        CALL SAVE_BOUNDARY_CONDS(i_column = n, gbcflag = iflag,
!     &                              wdbsl = wdbsl,
!     &                              wtmpdc = wtmpdc,
!     &                              wsalin = wsalin,
!     &                              wconc = wconc(:),
!     &                              wfflx = wfflx(:))

        CALL SAVE_BOUNDARY_CONDS(i_column = n, gbcflag = iflag,
     &                              wdata = wdata,
     &                              wconc = wconc(:),
     &                              wfflx = wfflx(:))
  
      ENDDO


                                    ! Final book-keeping:

                                    ! Increment the counter and the
                                    ! dataset-ID tracker now that the
                                    ! boundary condition set is complete.
      n_datasets_read = n_datasets_read + 1
      iloveclim_o2s_previd = iloveclim_o2s_id
    
      IF (rbflag /= -1) rbflag = 0  ! Clear <rbflag> only if it has not
                                    ! been set to -1 (do not clear the warning)

      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE OCEAN_TO_SEDIMENT
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      SUBROUTINE xchange_fluxes_o2s
     &  (temp_mafond, salt_mafond,
     &   oco2_mafond, oco3_mafond, ohco3_mafond,
     &   ooxy_mafond, ono3_mafond, opo4_mafond,
     &   clay_mafond, TPP_mafond, caco3_mafond,
     &   oc13_mafond, TPPC13_mafond, caco3C13_mafond,
     &   oc14_mafond, TPPC14_mafond, caco3C14_mafond,
     &   caco3O18_mafond, tracer01_mafond)

!-----------------------------------------------------------------------
!-----|--1--------2---------3---------4---------5---------6---------7-|
!     Cette subroutine contient l'initialisation des tableaux d'échange
!     de MEDUSA vers iLOVECLIM au moment de l'initialisation. Elle est
!     appelée une fois au début de simulation
!
!     ATTENTION : devrait etre possible de relire le restart de la
!     simulation precedente ???
!
!     (dans l'environnement logiciel LUDUS)
!
!     Auteur : Didier M. Roche & Guy Munhoven & Santiago Moreira
!     Date   : 12 fevrier 2014
!     Derniere modification : 12 avril 2022 (Guy&dmr)
!-----|--1--------2---------3---------4---------5---------6---------7-|

      DOUBLE PRECISION, DIMENSION(:,:),
     & INTENT(IN) :: TPP_mafond, caco3_mafond, temp_mafond,
     &               salt_mafond, ooxy_mafond,
     &               clay_mafond, tracer01_mafond,
     &               oco2_mafond, oco3_mafond, ohco3_mafond,
     &               ono3_mafond, opo4_mafond,
     &               oc13_mafond, oc14_mafond,
     &               TPPC13_mafond, caco3C13_mafond,
     &               TPPC14_mafond, caco3C14_mafond,
     &               caco3O18_mafond
     
        integer :: lbC, ubC, HPlb, HPub
        
! 8. Prepare one sediment step

! 8.1 Copy the integrated boundary fluxes and averaged boundary conditions
!     into MOD_ILOVECLIM_O2S

!-----|--1--------2---------3---------4---------5---------6---------7-|
!     CLIO_variable dimensions are:
!       CLIO_variable(Lon, Lat, Levels, nb_tracers) == (imax, jmax, kmax, nsmax)
!       ... and OCYCC_variable are:
!       OCYCC_variable(Lat, Levels, Lon) == (LT, JT, NOC_CBR)

!     When transfering variables from one grid to the other at same depth and same
!     tracer number, you need to transpose in lat, lon

!     Note a Benets : moreover the depth are inverted as follow:
!     first layer in OCYCC at surface and in the deep in CLIO

!     To transfer a variable in therefore follows:

!     FORALL (k=1:kmax)
!       odoc(:,(kmax+1-k),:) = TRANSPOSE(scal(:,:,k,3))
!     END FORALL

!-----|--1---------2---------3---------4---------5---------6---------7-|
! dmr Hence, the transpose included below for the _mafond variables that
!      are on the OCYCC grid and the iloveclim_ that are on the CLIO grid

!dmr --- 2021-09-22
!    --- Carbon variables size has changed due to removal of Halo Points in the Carbon grid.
!    --- Hence the restoring of Halo Points here to be coherent with CLIO

      lbC  = LBOUND(iloveclim_temp,dim=1)+1
      ubC  = UBOUND(iloveclim_temp,dim=1)-1
      HPlb = LBOUND(iloveclim_temp,dim=1)
      HPub = UBOUND(iloveclim_temp,dim=1)

      iloveclim_temp(lbC:ubC,:)   = TRANSPOSE(temp_mafond(:,:))
      iloveclim_temp(HPlb,:) = iloveclim_temp(ubC,:)
      iloveclim_temp(HPub,:) = iloveclim_temp(lbC,:)
            
      iloveclim_sali(lbC:ubC,:)   = TRANSPOSE(salt_mafond(:,:))
      iloveclim_sali(HPlb,:) = iloveclim_sali(ubC,:)
      iloveclim_sali(HPub,:) = iloveclim_sali(lbC,:)
            
                                    ! DBSL is constant in time
      !iloveclim_dbsl(:,:) ... left unchanged

      iloveclim_co2(lbC:ubC,:)    = TRANSPOSE(oco2_mafond(:,:))
      iloveclim_co2(HPlb,:) = iloveclim_co2(ubC,:)
      iloveclim_co2(HPub,:) = iloveclim_co2(lbC,:)
            
      iloveclim_co3(lbC:ubC,:)    = TRANSPOSE(oco3_mafond(:,:))
      iloveclim_co3(HPlb,:) = iloveclim_co3(ubC,:)
      iloveclim_co3(HPub,:) = iloveclim_co3(lbC,:)
            
      iloveclim_hco3(lbC:ubC,:)   = TRANSPOSE(ohco3_mafond(:,:))
      iloveclim_hco3(HPlb,:) = iloveclim_hco3(ubC,:)
      iloveclim_hco3(HPub,:) = iloveclim_hco3(lbC,:)      

      iloveclim_oxyg(lbC:ubC,:)   = TRANSPOSE(ooxy_mafond(:,:))
      iloveclim_oxyg(HPlb,:) = iloveclim_oxyg(ubC,:)
      iloveclim_oxyg(HPub,:) = iloveclim_oxyg(lbC,:)

      iloveclim_clay(lbC:ubC,:)   = TRANSPOSE(clay_mafond(:,:))
      iloveclim_clay(HPlb,:) = iloveclim_clay(ubC,:)
      iloveclim_clay(HPub,:) = iloveclim_clay(lbC,:)
      
      
      iloveclim_orgm(lbC:ubC,:)   = TRANSPOSE(TPP_mafond(:,:))
      iloveclim_orgm(HPlb,:) = iloveclim_orgm(ubC,:)
      iloveclim_orgm(HPub,:) = iloveclim_orgm(lbC,:)
            
      iloveclim_calc(lbC:ubC,:)   = TRANSPOSE(caco3_mafond(:,:))
      iloveclim_calc(HPlb,:) = iloveclim_calc(ubC,:)
      iloveclim_calc(HPub,:) = iloveclim_calc(lbC,:)

#ifdef WITH_C13
      iloveclim_oc13(lbC:ubC,:)    = TRANSPOSE(oc13_mafond(:,:))
      iloveclim_oc13(HPlb,:) = iloveclim_oc13(ubC,:)
      iloveclim_oc13(HPub,:) = iloveclim_oc13(lbC,:)      

      iloveclim_orgmC13(lbC:ubC,:) = TRANSPOSE(TPPC13_mafond(:,:))
      iloveclim_orgmC13(HPlb,:) = iloveclim_orgmC13(ubC,:)
      iloveclim_orgmC13(HPub,:) = iloveclim_orgmC13(lbC,:)  
            
      iloveclim_calcC13(lbC:ubC,:) = TRANSPOSE(caco3C13_mafond(:,:))
      iloveclim_calcC13(HPlb,:) = iloveclim_calcC13(ubC,:)
      iloveclim_calcC13(HPub,:) = iloveclim_calcC13(lbC,:)      
#endif

#ifdef WITH_C14
      iloveclim_oc14(lbC:ubC,:)    = TRANSPOSE(oc14_mafond(:,:))
      iloveclim_oc14(HPlb,:) = iloveclim_oc14(ubC,:)
      iloveclim_oc14(HPub,:) = iloveclim_oc14(lbC,:)     
      
      iloveclim_orgmC14(lbC:ubC,:) = TRANSPOSE(TPPC14_mafond(:,:))
      iloveclim_orgmC14(HPlb,:) = iloveclim_orgmC14(ubC,:)
      iloveclim_orgmC14(HPub,:) = iloveclim_orgmC14(lbC,:)
       
      iloveclim_calcC14(lbC:ubC,:) = TRANSPOSE(caco3C14_mafond(:,:))
      iloveclim_calcC14(HPlb,:) = iloveclim_calcC14(ubC,:)
      iloveclim_calcC14(HPub,:) = iloveclim_calcC14(lbC,:)
            
#endif

#ifdef WITH_O18
      iloveclim_calcO18(lbC:ubC,:) = TRANSPOSE(caco3O18_mafond(:,:))
      iloveclim_calcO18(HPlb,:) = iloveclim_calcO18(ubC,:)
      iloveclim_calcO18(HPub,:) = iloveclim_calcO18(lbC,:)      
#endif


                                    ! Non reactive tracer for testing purposes
      iloveclim_tracer01(lbC:ubC,:) = TRANSPOSE(tracer01_mafond(:,:))
      iloveclim_tracer01(HPlb,:) = iloveclim_tracer01(ubC,:)
      iloveclim_tracer01(HPub,:) = iloveclim_tracer01(lbC,:)  
      
      iloveclim_o2s_id = iloveclim_o2s_id + 1 
            
      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE xchange_fluxes_o2s
!-----------------------------------------------------------------------


!=======================================================================
      END MODULE MOD_ILOVECLIM_O2S
!=======================================================================
