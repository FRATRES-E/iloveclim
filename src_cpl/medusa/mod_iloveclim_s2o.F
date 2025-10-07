#define LOOP_ORG_MATTER
#define LOOP_CALCITE
#define LOOPBACK_LBTOT
! If "LOOPBACK_LBTOT" is #defined, the loopback fluxes are fed into the
! surface ocean layer, homogeneously distributed over the whole ocean area.
! if "LOOPBACK_LBTOT" is not #defined, a local distribution is used. 
!=======================================================================
      MODULE MOD_ILOVECLIM_S2O
!=======================================================================

      IMPLICIT NONE

      INTEGER, SAVE, PRIVATE :: n_grid_seafloor   = -1
      INTEGER, SAVE, PRIVATE :: nix               = -1
      INTEGER, SAVE, PRIVATE :: njy               = -1

                                    ! Has this module been setup
                                    ! (arrays allocated) ?

      LOGICAL, SAVE, PRIVATE :: mod_iloveclim_s2o_setup = .FALSE.

                                    ! Number of results sets
                                    ! previously written out:
                                    ! - not even initialised: -1

      INTEGER, SAVE          :: iloveclim_s2o_id = -1


!-----------------------------------------------------------------------
! Return fluxes to the ocean (top diffusive solute fluxes)
!-----------------------------------------------------------------------

      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_dic_tf
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_alk_tf
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_oxyg_tf
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_no3_tf
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_po4_tf

#ifdef WITH_C13
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_dic13_tf
#endif

#ifdef WITH_C14
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_dic14_tf
#endif


!-----------------------------------------------------------------------
! Burial fluxes out of the bottom of the sediment mixed layer
! (bottom advective solid fluxes)
!-----------------------------------------------------------------------
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm_bf
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calc_bf
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_clay_bf

#ifdef WITH_C13
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm13_bf
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calc13_bf
#endif

#ifdef WITH_C14
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm14_bf
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calc14_bf
#endif

#ifdef WITH_O18
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calc18_bf
#endif


!-----------------------------------------------------------------------
! Loop-back fluxes (if any)
!-----------------------------------------------------------------------

#ifdef LOOP_ORG_MATTER
! Save here the organic matter loss through the sediment bottom
! for reinjection by rivers flowing from the continent ("Loop-Back")
! We do this directly in terms of the solute constituents.

      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm_oxyg_lb  ! O2 to be consumed
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm_dic_lb   ! DIC produced
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm_alk_lb   ! ALK consumed
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm_no3_lb   ! NO3 produced
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm_po4_lb   ! PO4 produced

      DOUBLE PRECISION, SAVE :: iloveclim_orgm_oxyg_lbtot
      DOUBLE PRECISION, SAVE :: iloveclim_orgm_dic_lbtot
      DOUBLE PRECISION, SAVE :: iloveclim_orgm_alk_lbtot
      DOUBLE PRECISION, SAVE :: iloveclim_orgm_no3_lbtot
      DOUBLE PRECISION, SAVE :: iloveclim_orgm_po4_lbtot

#  ifdef WITH_C13
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm_dic13_lb ! DIC-13 produced

      DOUBLE PRECISION, SAVE :: iloveclim_orgm_dic13_lbtot
#  endif

#  ifdef WITH_C14
                                    ! The C14 loop-back flux should be zero
                                    ! (dead carbon)
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm_dic14_lb ! DIC-14 produced

      DOUBLE PRECISION, SAVE :: iloveclim_orgm_dic14_lbtot
#  endif

#endif

#ifdef LOOP_CALCITE
! Save here the calcite loss through the sediment bottom
! for reinjection by rivers flowing from the continent ("Loop-Back")
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calc_lb

      DOUBLE PRECISION, SAVE :: iloveclim_calc_lbtot

#  ifdef WITH_C13
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calc13_lb      

      DOUBLE PRECISION, SAVE :: iloveclim_calc13_lbtot
#  endif

#  ifdef WITH_C14
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calc14_lb

      DOUBLE PRECISION, SAVE :: iloveclim_calc14_lbtot
#  endif

#endif


! total area
      DOUBLE PRECISION, SAVE :: sfc_area_tot


!-----------------------------------------------------------------------
! Total amounts of materials in the sediment columns
!-----------------------------------------------------------------------

      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm_ts
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calc_ts
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_clay_ts

#ifdef WITH_C13
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm13_ts
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calc13_ts
#endif

#ifdef WITH_C14
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_orgm14_ts
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calc14_ts
#endif

#ifdef WITH_O18
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), ALLOCATABLE
     &   :: iloveclim_calc18_ts
#endif


      CONTAINS


!-----------------------------------------------------------------------
      SUBROUTINE setup_s2o
!-----------------------------------------------------------------------
!-----|--1---------2---------3---------4---------5---------6---------7-|
!     Cette subroutine contient l'allocation des tableaux d'échange
!     MEDUSA vers iLOVECLIM au moment de l'initialisation. Elle est
!     appelée une fois au début de simulation
!     (dans l'environnement logiciel LUDUS)
!
!     Auteur : Didier M. Roche & Guy Munhoven
!     Date   : 05 fevrier 2014
!     Dernière modification : 20 octobre 2016 (Guy Munhoven)
!-----|--1--------2---------3---------4---------5---------6---------7-|
 
      USE mod_defines_medusa,       ONLY: jp_stderr
      USE mod_execontrol_medusa,    ONLY: ABORT_MEDUSA
      USE mod_seafloor_central,     ONLY: IJ_COLUMNS_USED,
     &                                    N_COLUMNS_USED


      IMPLICIT NONE


      CHARACTER(LEN=*), PARAMETER :: c_fmterr_a =
     &  '("[MOD_ILOVECLIM_SETUP/setup_s2o] error: ", A)'


      IF (.NOT. mod_iloveclim_s2o_setup) THEN

        CALL IJ_COLUMNS_USED(nix, njy)
        CALL N_COLUMNS_USED(n_grid_seafloor)

        IF (n_grid_seafloor /= -1) THEN

          ALLOCATE(iloveclim_dic_tf(nix,njy))
          ALLOCATE(iloveclim_alk_tf(nix,njy))
          ALLOCATE(iloveclim_oxyg_tf(nix,njy))
          ALLOCATE(iloveclim_no3_tf(nix,njy))
          ALLOCATE(iloveclim_po4_tf(nix,njy))

          ALLOCATE(iloveclim_clay_bf(nix,njy))
          ALLOCATE(iloveclim_orgm_bf(nix,njy))
          ALLOCATE(iloveclim_calc_bf(nix,njy))

#ifdef WITH_C13
          ALLOCATE(iloveclim_dic13_tf(nix,njy))
#endif

#ifdef WITH_C13
          ALLOCATE(iloveclim_dic14_tf(nix,njy))
#endif


#ifdef LOOP_ORG_MATTER
          print*, "[MOD_ILOVECLIM_O2S] OrgMatter Loopback ENabled"
#  ifdef LOOPBACK_LBTOT
          print*, " - loopback fluxes distributed homogeneously " //
     &               "over the whole ocean surface"
#  else
          print*, " - loopback fluxes distributed somehow(?) " //
     &               "over the ocean surface"
#  endif
          ALLOCATE(iloveclim_orgm_dic_lb(nix,njy))
          ALLOCATE(iloveclim_orgm_alk_lb(nix,njy))
          ALLOCATE(iloveclim_orgm_oxyg_lb(nix,njy))
          ALLOCATE(iloveclim_orgm_no3_lb(nix,njy))
          ALLOCATE(iloveclim_orgm_po4_lb(nix,njy))

#  ifdef WITH_C13
          ALLOCATE(iloveclim_orgm_dic13_lb(nix,njy))
#  endif
#  ifdef WITH_C14
          ALLOCATE(iloveclim_orgm_dic14_lb(nix,njy))
#  endif
#else
          print*, "[MOD_ILOVECLIM_O2S] OrgMatter Loopback DISabled"
#endif

#ifdef LOOP_CALCITE
          print*, "[MOD_ILOVECLIM_O2S] Calcite Loopback ENabled"
#  ifdef LOOPBACK_LBTOT
          print*, " - loopback fluxes distributed homogeneously " //
     &               "over the whole ocean surface"
#  else
          print*, " - loopback fluxes distributed somehow(?) " //
     &               "over the ocean surface"
#  endif
          ALLOCATE(iloveclim_calc_lb(nix,njy))
#  ifdef WITH_C13
          ALLOCATE(iloveclim_calc13_lb(nix,njy))
#  endif
#  ifdef WITH_C14
          ALLOCATE(iloveclim_calc14_lb(nix,njy))
#  endif
#else
          print*, "[MOD_ILOVECLIM_O2S] Calcite Loopback DISabled"
#endif


          ALLOCATE(iloveclim_clay_ts(nix,njy))
          ALLOCATE(iloveclim_orgm_ts(nix,njy))
          ALLOCATE(iloveclim_calc_ts(nix,njy))

#ifdef WITH_C13
          ALLOCATE(iloveclim_orgm13_ts(nix,njy))
          ALLOCATE(iloveclim_calc13_ts(nix,njy))
#endif

#ifdef WITH_C14
          ALLOCATE(iloveclim_orgm14_ts(nix,njy))
          ALLOCATE(iloveclim_calc14_ts(nix,njy))
#endif

#ifdef WITH_O18
          ALLOCATE(iloveclim_calc18_ts(nix,njy))
#endif

          mod_iloveclim_s2o_setup = .TRUE.

        ELSE

          WRITE(jp_stderr,c_fmterr_a)
     &      'MOD_SEAFLOOR_CENTRAL not yet set up -- aborting!'
          CALL ABORT_MEDUSA()

        ENDIF

      ELSE

        WRITE(jp_stderr,c_fmterr_a)
     &         'MOD_ILOVECLIM_S2O already set up -- aborting!'
        CALL ABORT_MEDUSA()

      ENDIF


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE setup_s2o
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
      SUBROUTINE ini_fluxes_s2o
!-----------------------------------------------------------------------
!-----|--1---------2---------3---------4---------5---------6---------7-|
!     Cette subroutine contient l'initialisation des tableaux d'échange
!     de MEDUSA vers iLOVECLIM au moment de l'initialisation. Elle est
!     appelée une fois au début de simulation.
!     (dans l'environnement logiciel LUDUS)
!     Ne peut être appelé qu'après set-up du sédiment.
!
!     Auteur : Guy Munhoven
!     Date   : 11 octobre 2016
!     Dernière modification : 20 octobre 2016 (Guy Munhoven)
!-----|--1---------2---------3---------4---------5---------6---------7-|

      USE mod_defines_medusa
      USE mod_execontrol_medusa,    ONLY: ABORT_MEDUSA
      USE mod_indexparam,           ONLY: nsolut, nsolid


      IMPLICIT NONE


      INTEGER :: n
      INTEGER :: iflag

      CHARACTER(LEN=*), PARAMETER :: c_fmterr_a =
     &  '("[MOD_ILOVECLIM_S2O/ini_fluxes_s2o] error: ", A)'



      IF (.NOT. mod_iloveclim_s2o_setup) THEN
        WRITE(jp_stderr,c_fmterr_a)
     &    'MOD_ILOVECLIM_S2O not yet set up -- aborting!'
        CALL ABORT_MEDUSA()
      ENDIF

                                    ! Set the sediment-to-ocean
                                    ! fluxes in MOD_ILOVECLIM_S2O for
                                    ! iLoveclim to run for one sediment
                                    ! step to zero
      iloveclim_dic_tf(:,:)        = 0.0D+00
      iloveclim_alk_tf(:,:)        = 0.0D+00
      iloveclim_oxyg_tf(:,:)       = 0.0D+00
      iloveclim_no3_tf(:,:)        = 0.0D+00
      iloveclim_po4_tf(:,:)        = 0.0D+00

      iloveclim_orgm_bf(:,:)       = 0.0D+00
      iloveclim_calc_bf(:,:)       = 0.0D+00
      iloveclim_clay_bf(:,:)       = 0.0D+00


#ifdef WITH_C13
      iloveclim_dic13_tf(:,:)      = 0.0D+00

      iloveclim_orgm13_bf(:,:)     = 0.0D+00
      iloveclim_calc13_bf(:,:)     = 0.0D+00
#endif


#ifdef WITH_C14
      iloveclim_dic14_tf(:,:)      = 0.0D+00

      iloveclim_orgm14_bf(:,:)     = 0.0D+00
      iloveclim_calc14_bf(:,:)     = 0.0D+00
#endif


#ifdef WITH_O18
      iloveclim_calc18_bf(:,:)     = 0.0D+00
#endif


#ifdef LOOP_ORG_MATTER
      iloveclim_orgm_dic_lb(:,:)   = 0.0D+00
      iloveclim_orgm_alk_lb(:,:)   = 0.0D+00
      iloveclim_orgm_oxyg_lb(:,:)  = 0.0D+00
      iloveclim_orgm_no3_lb(:,:)   = 0.0D+00
      iloveclim_orgm_po4_lb(:,:)   = 0.0D+00

      iloveclim_orgm_oxyg_lbtot    = 0.0D+00
      iloveclim_orgm_dic_lbtot     = 0.0D+00
      iloveclim_orgm_alk_lbtot     = 0.0D+00
      iloveclim_orgm_no3_lbtot     = 0.0D+00
      iloveclim_orgm_po4_lbtot     = 0.0D+00

#  ifdef WITH_C13
      iloveclim_orgm_dic13_lb(:,:) = 0.0D+00
      iloveclim_orgm_dic13_lbtot   = 0.0D+00
#  endif

#  ifdef WITH_C14
      iloveclim_orgm_dic14_lb(:,:) = 0.0D+00
      iloveclim_orgm_dic14_lbtot   = 0.0D+00
#  endif
#endif


#ifdef LOOP_CALCITE
      iloveclim_calc_lb(:,:)       = 0.0D+00

      iloveclim_calc_lbtot         = 0.0D+00

#  ifdef WITH_C13
      iloveclim_calc13_lb(:,:)     = 0.0D+00
      iloveclim_calc13_lbtot       = 0.0D+00
#  endif

#  ifdef WITH_C14
      iloveclim_calc14_lb(:,:)     = 0.0D+00
      iloveclim_calc14_lbtot       = 0.0D+00
#  endif
#endif


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE ini_fluxes_s2o
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
      SUBROUTINE SEDIMENT_TO_OCEAN(kflag)
!-----------------------------------------------------------------------

      ! Returns model outcome to the arrays in this module, carrying
      ! out unit conversions appropriate for iLoveclim
      ! so that the results can be directly used in the
      ! appropriate subroutines of the biogeochemistry module
      ! Can only be called after the sediment has been initialised
      ! (i.e., the sediment state arrays in MOD_SEAFLOOR_CENTRAL
      ! have been initialised, else tsolids cannot be calculated).

      USE mod_defines_medusa, ONLY: jp_stderr
      USE mod_execontrol_medusa, ONLY: ABORT_MEDUSA
      USE mod_indexparam
      USE mod_materialcharas
      USE mod_seafloor_central, ONLY:
     &                              COLUMN_IJ2N, GET_BOUNDARY_FLUXES,
     &                              GET_MATERIALCHARAS, COLUMN_AREA4N
      USE mod_zone_mass_totals, only:
     &                              REACLAY_TOTALS

      IMPLICIT NONE

      INTEGER, INTENT(OUT)                :: kflag


      DOUBLE PRECISION, DIMENSION(nsolut) :: wcflx
      DOUBLE PRECISION, DIMENSION(nsolid) :: bfflx
      DOUBLE PRECISION, DIMENSION(nsolid) :: tsolids
      DOUBLE PRECISION                    :: sfc_area

      INTEGER :: i,j,n, iflag, sumfluxfile

      CHARACTER(LEN=*), PARAMETER :: c_fmterr_a =
     &  '("[MOD_ILOVECLIM_SETUP/SEDIMENT_TO_OCEAN] error: ", A)'


      IF (.NOT. mod_iloveclim_s2o_setup) THEN
        WRITE(jp_stderr,c_fmterr_a)
     &    'MOD_ILOVECLIM_S2O not yet set up -- aborting!'
        CALL ABORT_MEDUSA()
      ENDIF

#ifdef DEBUG
      print*, "SUBROUTINE SEDIMENT_TO_OCEAN"
#endif

                                    ! Prepare loop-back totals:
                                    ! always reste these to zero

       sfc_area_tot = 0.0D+00       ! Reset total area to zero 

#ifdef LOOP_ORG_MATTER
                                    ! Reset the total loop-back fluxes
                                    ! to zero if loop-back is active
      iloveclim_orgm_oxyg_lbtot   = 0.0D+00
      iloveclim_orgm_dic_lbtot    = 0.0D+00
      iloveclim_orgm_alk_lbtot    = 0.0D+00
      iloveclim_orgm_no3_lbtot    = 0.0D+00
      iloveclim_orgm_po4_lbtot    = 0.0D+00

#  ifdef WITH_C13
      iloveclim_orgm_dic13_lbtot  = 0.0D+00
#  endif

#  ifdef WITH_C14
      iloveclim_orgm_dic14_lbtot  = 0.0D+00
#  endif

#endif

#ifdef LOOP_CALCITE
      ! Reset the total Loop-back fluxes to zero if loop-back is active
      iloveclim_calc_lbtot        = 0.0D+00

#  ifdef WITH_C13
      iloveclim_calc13_lbtot      = 0.0D+00
#  endif

#  ifdef WITH_C14
      iloveclim_calc14_lbtot      = 0.0D+00
#  endif

#endif

    

      DO i = 1, nix
!      DO i = 2, nix-1 !test remove halo

         DO j = 1, njy

          CALL COLUMN_IJ2N(i, j, iflag, n)

          IF (iflag == 0) THEN      ! if there exists n for (i,j)


            CALL GET_BOUNDARY_FLUXES(i_column = n, gbfflag = iflag,
     &                              wcflx = wcflx(:),
     &                              bfflx = bfflx(:))

            CALL GET_MATERIALCHARAS(i_column = n, iflag = iflag)

            CALL COLUMN_AREA4N(i_column = n, iflag = iflag,
     &                              sfc_area = sfc_area)

            sfc_area_tot = sfc_area_tot + sfc_area


                                    ! - Solutes

                                    ! Units:
                                    !  - wcflx [mol/m2/yr]
                                    !  - iloveclim_dic_tf [mol/grid_elt/yr] --- dmr&gm 2021-11-17
            iloveclim_dic_tf(i,j) =
     &        (wcflx(ic_co3) + wcflx(ic_hco3) + wcflx(ic_co2))*sfc_area

                                    ! Units:
                                    !  - wcflx [mol/m2/yr]
                                    !  - iloveclim_alk_tf [eq/grid_elt/yr] ---[gm dixit: unit is mol equivalent H+]
            iloveclim_alk_tf(i,j) =
     &        (wcflx(ic_hco3) + wcflx(ic_co3) + wcflx(ic_co3))*sfc_area

                                    ! Units:
                                    !  - wcflx [mol/m2/yr]
                                    !  - iloveclim_oxyg_tf [mol/grid_elt/yr]
            iloveclim_oxyg_tf(i,j) =
     &         wcflx(ic_o2)*sfc_area


            ! GM: Currently, NO3 and PO4 are not explicitly considered in
            ! the configuration of Medusa coupled to iLoveclim
            ! We therefore derive the NO3 and PO4 fluxes from the
            ! O2 flux. Notice that when O2 is taken up, NO3 and PO4 are
            ! released. Therefore the minus-sign. Units are here the
            ! same as for O2 (currently [mol/column/yr])

                                    ! NO3 derived from O2 consumption rate
                                    !  - iloveclim_no3_tf [mol/grid_elt/yr]
            iloveclim_no3_tf(i,j) =
     &         -iloveclim_oxyg_tf(i,j)/om_remin_o2*om_n


                                    ! PO4 derived from O2 consumption
                                    !  - iloveclim_po4_tf [mol/grid_elt/yr]
            iloveclim_po4_tf(i,j) =
     &         -iloveclim_oxyg_tf(i,j)/om_remin_o2*om_p


                                    ! - Solids

                                    ! Units:
                                    !  - bfflx [kg/m2/yr]
                                    !  - iloveclim_clay_bf [kg/grid_elt/yr]
            iloveclim_clay_bf(i,j) =
     &         bfflx(if_clay)*sfc_area

                                    ! Units:
                                    !  - bfflx [kg/m2/yr]
                                    !  - iloveclim_calc_bf [molC/grid_elt/yr]
            iloveclim_calc_bf(i,j) =
     &        (bfflx(if_calc)/mol_calc)*sfc_area

                                    ! Units:
                                    !  - bfflx [kg/m2/yr]
                                    !  - iloveclim_orgm_bf [molC/grid_elt/yr]
            iloveclim_orgm_bf(i,j) =
     &        (bfflx(if_om)/mol_om*om_c)*sfc_area



            CALL REACLAY_TOTALS(i_column = n, iflag = iflag,
     &                              tsolids = tsolids(:))

                                    ! Units:
                                    !  - tsolids [kg/m2]
                                    !  - iloveclim_clay_ts [kg/grid_elt]
            iloveclim_clay_ts(i,j) =
     &          tsolids(if_clay)*sfc_area

                                    ! Units:
                                    !  - tsolids [kg/m2]
                                    !  - iloveclim_calc_ts [molC/grid_elt]
            iloveclim_calc_ts(i,j) =
     &         (tsolids(if_calc)/mol_calc)*sfc_area

                                    ! Units:
                                    !  - tsolids [kg/m2]
                                    !  - iloveclim_orgm_ts [molC/grid_elt]
            iloveclim_orgm_ts(i,j) =
     &         (tsolids(if_om)/mol_om*om_c)*sfc_area

#ifdef WITH_C13
                                    ! Units:
                                    !  - tsolids [kg/m2]
                                    !  - iloveclim_calc13_ts [molC/grid_elt]
            iloveclim_calc13_ts(i,j) =
     &         (tsolids(if_calc_c13)/mol_calc)*sfc_area

                                    ! Units:
                                    !  - tsolids [kg/m2]
                                    !  - iloveclim_orgm_ts [molC/grid_elt]
            iloveclim_orgm13_ts(i,j) =
     &         (tsolids(if_om_c13)/mol_om*om_c)*sfc_area
#endif

#ifdef WITH_C14
                                    ! Units:
                                    !  - tsolids [kg/m2]
                                    !  - iloveclim_calc13_ts [molC/grid_elt]
            iloveclim_calc14_ts(i,j) =
     &         (tsolids(if_calc_c14)/mol_calc)*sfc_area

                                    ! Units:
                                    !  - tsolids [kg/m2]
                                    !  - iloveclim_orgm_ts [molC/grid_elt]
            iloveclim_orgm14_ts(i,j) =
     &         (tsolids(if_om_c14)/mol_om*om_c)*sfc_area
#endif

#ifdef WITH_O18
                                    ! Units:
                                    !  - tsolids [kg/m2]
                                    !  - iloveclim_calc13_ts [molC/grid_elt]
            iloveclim_calc18_ts(i,j) =
     &         (tsolids(if_calc_o18)/mol_calc)*sfc_area
#endif



#ifdef LOOP_ORG_MATTER
            ! GM: iloveclim_orgm_dic_lb/tot should be *added* to the atmosphere
            ! (amount of CO2 produced during fossil organic matter oxidation
            ! on land)
                                    ! Units:
                                    !  - bfflx(if_om) [kg OrgMatter/m2/yr], source
                                    !  - iloveclim_orgm_dic_lb [molC/m2/yr], source
            iloveclim_orgm_dic_lb(i,j)
     &        = bfflx(if_om)/mol_om*om_c

                                    !  - iloveclim_orgm_dic_lbtot [molC/yr], source
            iloveclim_orgm_dic_lbtot
     &        = iloveclim_orgm_dic_lbtot
     &          + iloveclim_orgm_dic_lb(i,j)*sfc_area


            ! GM: Please check in the ocean carbon cycle model if the alkalinity
            ! release into the water during organic matter production is 
            ! taken into account. If not, set iloveclim_orgm_alk_lb(i,j) to zero!
            ! iloveclim_orgm_alk_lb/tot must be *subtracted* from the surface ocean,
            ! following some distribution (reflecting riverine input)
                                    ! Units:
                                    !  - bfflx(if_om) [kg OrgMatter/m2/yr], source
                                    !  - iloveclim_orgm_alk_lb [molO2/m2/yr], SINK
 
            ! OCYCC does not consider alkalinity release during photosynthesis
!            iloveclim_orgm_alk_lb(i,j) = 0.0D+00
            ! If it gets this one day, please uncomment the next two lines
            iloveclim_orgm_alk_lb(i,j)
     &        = -bfflx(if_om)/mol_om*(om_p + om_n)

                                    !  - iloveclim_orgm_alk_lbtot [molC/yr], sink
            iloveclim_orgm_alk_lbtot
     &        = iloveclim_orgm_alk_lbtot
     &          + iloveclim_orgm_alk_lb(i,j)*sfc_area


            ! GM: O2 to be considered only if iLoveclim has a complete
            ! oxygen balance. Must be set to zero if, e.g., the atmosphere
            ! is considered an infinite reservoir.
                                    ! Units:
                                    !  - bfflx(if_om) [kg OrgMatter/m2/yr], source
                                    !  - iloveclim_orgm_oxyg_lb [molO2/m2/yr], SINK
            iloveclim_orgm_oxyg_lb(i,j)
     &        = -bfflx(if_om)/mol_om*om_remin_o2

                                    !  - iloveclim_orgm_oxyg_lbtot [molO2/yr], sink
            iloveclim_orgm_oxyg_lbtot
     &        = iloveclim_orgm_oxyg_lbtot
     &          + iloveclim_orgm_oxyg_lb(i,j)*sfc_area


            ! GM: iloveclim_orgm_no3_lb/tot should be *added* to the surface ocean
            ! (amount of net nitrification)
                                    ! Units:
                                    !  - bfflx(if_om) [kg OrgMatter/m2/yr], source
                                    !  - iloveclim_orgm_no3_lb [molNO3/m2/yr], source
            iloveclim_orgm_no3_lb(i,j)
     &        = bfflx(if_om)/mol_om*om_n

                                    !  - iloveclim_orgm_no3_lbtot [molNO3/yr], source
            iloveclim_orgm_no3_lbtot
     &        = iloveclim_orgm_no3_lbtot
     &          + iloveclim_orgm_no3_lb(i,j)*sfc_area


            ! GM: iloveclim_orgm_po4_lb/tot should be *added* to the surface ocean
            ! (amount of net phosphate erosion on land),
            ! following some distribution (reflecting riverine input from weathering)
                                    ! Units:
                                    !  - bfflx(if_om) [kg OrgMatter/m2/yr], source
                                    !  - iloveclim_orgm_po4_lb [molPO4/m2/yr], source
            iloveclim_orgm_po4_lb(i,j)
     &        = bfflx(if_om)/mol_om*om_p

                                    !  - iloveclim_orgm_po4_lbtot [molPO4/yr], source
            iloveclim_orgm_po4_lbtot
     &        = iloveclim_orgm_po4_lbtot
     &          + iloveclim_orgm_po4_lb(i,j)*sfc_area
#endif



#ifdef LOOP_CALCITE
            ! GM: iloveclim_calc_lb/tot should be used as follows:
            !  - add 2 * iloveclim_calc_lb/tot as DIC (sic!) to the surface ocean
            !  - add 2 * iloveclim_calc_lb/tot as ALK to the surface ocean
            !  - take out 1 * iloveclim_calc_lb/tot as CO2 from the atmosphere
            ! These repsectively represent the amounts of DIC and ALK
            ! produced, and of atmospheric CO2 consumed during carbonate
            ! weathering on land, following  CaCO3 + CO2 + H2O -> 2 HCO3

                                    ! Units:
                                    !  - bfflx(if_calc) [kg Calcite/m2/yr], source
                                    !  - iloveclim_calc_lb [molC/m2/yr], source
            iloveclim_calc_lb(i,j)
     &        = bfflx(if_calc)/mol_calc*calc_c

                                    !  - iloveclim_calc_lbtot [molC/yr], source
            iloveclim_calc_lbtot
     &        = iloveclim_calc_lbtot
     &          + iloveclim_calc_lb(i,j)*sfc_area
#endif



          ENDIF

        ENDDO

      ENDDO


#if ( 1 ) 

       open(newunit=sumfluxfile
     &        ,file='outputdata/sediments/sum_fluxes.txt'
     &        ,status='old',action='write', form='formatted'
     &        ,position="append")
     
       write(sumfluxfile,'(10F20.12)') 
     &       SUM(iloveclim_dic_tf) /1.E12,
     &       SUM(iloveclim_alk_tf) /1.E12,
     &       SUM(iloveclim_oxyg_tf)/1.E12,
     &       SUM(iloveclim_po4_tf) /1.E12,
     &       SUM(iloveclim_calc_bf)/1.E12,
     &       SUM(iloveclim_orgm_bf)/1.E12
       close(sumfluxfile)

#endif


      iloveclim_s2o_id = iloveclim_s2o_id + 1

      kflag = 0


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE SEDIMENT_TO_OCEAN
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
      SUBROUTINE xchange_fluxes_s2o
     &  (ODOCS_sed2oc, ODIC_sed2oc, OALK_sed2oc,
     &   OO2_sed2oc, ONO3_sed2oc, OPO4_sed2oc,
     &   ODIC13_sed2oc, ODIC14_sed2oc,
     &   orgm_dic_loopback, orgm_alk_loopback,
     &   orgm_oxyg_loopback, orgm_no3_loopback, orgm_po4_loopback,
     &   calc_loopback,
     &   orgm_dic13_loopback, calc13_loopback,
     &   orgm_dic14_loopback, calc14_loopback)
!-----------------------------------------------------------------------

      IMPLICIT NONE
!-----|--1---------2---------3---------4---------5---------6---------7-|
!     This subroutine prepares the distribution of dissolved fluxes
!     arising from ocean-sediment exchange as derived from Medusa
!     calculations, and returns them in arrays suitable for usage
!     in the carbon model, on 2D horizontal iLoveclim grids.
!
!     Two classes of such fluxes are considered
!
!      1. diffusive transport of dissolution and remineralization
!         products across the ocean-sediment interface, exchanged
!         with the bottom ocean:
!          - ODOCS_sed2oc (provided, but set to 0 for the time being)
!          - ODIC_sed2oc
!          - OALK_sed2oc
!          - OO2_sed2oc
!          - ONO3_sed2oc
!          - OPO4_sed2oc
!         These "O<solute>_sed2oc" are derived from the corresponding
!         "iloveclim_<solute>_tf" arrays (from this module), which have
!         been previously initialized by the SEDIMENT_TO_OCEAN subroutine
!         (also from this module).
!
!      2. solid burial fluxes, across the bottom-of-the-sediment,
!         expressed in terms of equivalent solute fluxes that would
!         arise from complete dissolution and remineralization
!         of these solids, and that are provided for (optional) usage
!         as riverine inputs (and atmospheric consumption and release
!         fluxes) to the surface ocean:
!          - orgm_dic_loopback
!          - orgm_alk_loopback
!          - orgm_no3_loopback
!          - orgm_po4_loopback
!          - orgm_oxyg_loopback
!          - calc_loopback [molC/gridelt/yr]
!         For each solid (orgm, calc), the corresponding group of arrays
!         is set to 0 if the respective pre-processor token
!         (LOOP_ORG_MATTER, LOOP_CALCITE) is #undefined.
!         If the respective pre-processor token (LOOP_ORG_MATTER, LOOP_CALCITE)
!         is #defined, these "<solid>_<solute>_loopback" arrays are
!         currently derived as follows:
!          - if LOOPBACK_LBTOT is #defined they are derived from the
!            corresponding "iloveclim_<solid>_<solut>_lbtot" arrays
!            (from this module), which have been previously initialized
!            by the SEDIMENT_TO_OCEAN subroutine (also from this module)
!          - if LOOPBACK_LBTOT is #undefined, they are obtained by
!            remapping "iloveclim_<solid>_<solut>_lb"
!
!         **************************************************************
!         * NONE OF THE LOOPBACK OPTIONS IS REALISTIC AND THEY MUST    *
!         * ONLY BE USED FOR DEVELOPMENT (MASS CONSERVATION) PURPOSES. *
!         * THEY ABSOLUTELY HAVE TO BE REPLACED BY A DECENT RIVERINE   *
!         * INPUT DISTRIBUTION SCHEME.                                 *
!         **************************************************************
!
!     Auteur : Santiago Moreira (mohr) & Guy Munhoven
!     Date   : 2015-08-26
!     Derniere modification :  17 decembre 2024 (Guy Munhoven)
!     References: ----
!-----|--1---------2---------3---------4---------5---------6---------7-|
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT)
     &  :: ODOCS_sed2oc, ODIC_sed2oc, OALK_sed2oc,  OO2_sed2oc,
     &     ONO3_sed2oc, OPO4_sed2oc

      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT)
     &  :: orgm_dic_loopback, orgm_alk_loopback,
     &     orgm_no3_loopback, orgm_po4_loopback, orgm_oxyg_loopback

      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT)
     &  :: calc_loopback

      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT)
     &  :: ODIC13_sed2oc, orgm_dic13_loopback, calc13_loopback

      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT)
     &  :: ODIC14_sed2oc, orgm_dic14_loopback, calc14_loopback

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: XChange_Array

      INTEGER :: lbnd, ubnd
      ! [XXX] Check if All set up (not eq to -1)

#ifdef DEBUG
      print*,"SUBROUTINE xchange_fluxes_s2o"
#endif

      lbnd = LBOUND(iloveclim_dic_tf,DIM=1)+1
      ubnd = UBOUND(iloveclim_dic_tf,DIM=1)-1

      ALLOCATE(XChange_Array(                                          
     & LBOUND(iloveclim_dic_tf,DIM=2):UBOUND(iloveclim_dic_tf,DIM=2)   
     &,LBOUND(iloveclim_dic_tf,DIM=1):UBOUND(iloveclim_dic_tf,DIM=1))  
     &        )
      
      XChange_Array = TRANSPOSE(iloveclim_dic_tf(:,:))
      ODIC_sed2oc   = XChange_Array(:,lbnd:ubnd)
      print*, "XChange fluxes ODIC_sed2oc ... min/max" , 
     &  MINVAL(ODIC_sed2oc(:,:)), MAXVAL(ODIC_sed2oc(:,:))
      XChange_Array = TRANSPOSE(iloveclim_alk_tf(:,:))
      OALK_sed2oc   = XChange_Array(:,lbnd:ubnd)
      
      XChange_Array = TRANSPOSE(iloveclim_oxyg_tf(:,:))
      OO2_sed2oc    = XChange_Array(:,lbnd:ubnd)
      
      ODOCS_sed2oc  = 0.0D+00       !there is no value of a DOC tf flux now
      
      XChange_Array = TRANSPOSE(iloveclim_no3_tf(:,:))
      ONO3_sed2oc   = XChange_Array(:,lbnd:ubnd)
      
      XChange_Array = TRANSPOSE(iloveclim_po4_tf(:,:))
      OPO4_sed2oc   = XChange_Array(:,lbnd:ubnd)

#ifdef WITH_C13
      ODIC13_sed2oc = TRANSPOSE(iloveclim_dic13_tf(:,:))(:,lbnd:ubnd)
#else
      ODIC13_sed2oc = 0.0D+00
#endif
#ifdef WITH_C14
      ODIC14_sed2oc = TRANSPOSE(iloveclim_dic14_tf(:,:))(:,lbnd:ubnd)
#else
      ODIC14_sed2oc = 0.0D+00
#endif


#ifdef LOOP_ORG_MATTER

#  ifdef LOOPBACK_LBTOT
!test
!      iloveclim_orgm_dic_lbtot=6540838384951.56
!      iloveclim_orgm_alk_lbtot=-1058156037428.06
!      iloveclim_orgm_po4_lbtot=62244472789.8858

      orgm_dic_loopback(:,:)   = iloveclim_orgm_dic_lbtot / sfc_area_tot
      orgm_alk_loopback(:,:)   = iloveclim_orgm_alk_lbtot / sfc_area_tot
      orgm_oxyg_loopback(:,:)  = iloveclim_orgm_oxyg_lbtot /sfc_area_tot
      orgm_no3_loopback(:,:)   = iloveclim_orgm_no3_lbtot / sfc_area_tot
      orgm_po4_loopback(:,:)   = iloveclim_orgm_po4_lbtot / sfc_area_tot

!      write(*,*) 'TEST SED loopback', iloveclim_orgm_dic_lbtot, 
!     &  iloveclim_orgm_alk_lbtot, iloveclim_orgm_po4_lbtot

#    ifdef WITH_C13
      orgm_dic13_loopback(:,:) = iloveclim_orgm_dic13_lbtot/sfc_area_tot
#    else
      orgm_dic13_loopback(:,:) = 0.0D+00
#    endif

#    ifdef WITH_C14
      orgm_dic14_loopback(:,:) = iloveclim_orgm_dic14_lbtot/sfc_area_tot
#    else
      orgm_dic14_loopback(:,:) = 0.0D+00
#    endif

#  else

      orgm_dic_loopback(:,:)   = TRANSPOSE(iloveclim_orgm_dic_lb(:,:))(:,lbnd:ubnd) ! modifier les parentheses (attention halo points removed)
      orgm_alk_loopback(:,:)   = TRANSPOSE(iloveclim_orgm_alk_lb(:,:))(:,lbnd:ubnd)
      orgm_oxyg_loopback(:,:)  = TRANSPOSE(iloveclim_orgm_oxyg_lb(:,:))(:,lbnd:ubnd)
      orgm_no3_loopback(:,:)   = TRANSPOSE(iloveclim_orgm_no3_lb(:,:))(:,lbnd:ubnd)
      orgm_po4_loopback(:,:)   = TRANSPOSE(iloveclim_orgm_po4_lb(:,:))(:,lbnd:ubnd)

#    ifdef WITH_C13
      orgm_dic13_loopback(:,:) = TRANSPOSE(iloveclim_orgm_dic13_lb(:,:))(:,lbnd:ubnd)
#    else
      orgm_dic13_loopback(:,:) = 0.0D+00
#    endif

#    ifdef WITH_C14
      orgm_dic14_loopback(:,:) = TRANSPOSE(iloveclim_orgm_dic14_lb(:,:))(:,lbnd:ubnd)
#    else
      orgm_dic14_loopback(:,:) = 0.0D+00
#    endif

#  endif

#else

      orgm_oxyg_loopback(:,:)  =  0.0D+00
      orgm_dic_loopback(:,:)   =  0.0D+00
      orgm_alk_loopback(:,:)   =  0.0D+00
      orgm_no3_loopback(:,:)   =  0.0D+00
      orgm_po4_loopback(:,:)   =  0.0D+00

      orgm_dic13_loopback(:,:) =  0.0D+00

      orgm_dic14_loopback(:,:) =  0.0D+00

#endif


#ifdef LOOP_CALCITE

#  ifdef LOOPBACK_LBTOT
!test
!      iloveclim_calc_lbtot=3336628474912.36

      calc_loopback(:,:)   = iloveclim_calc_lbtot / sfc_area_tot
!      write(*,*) 'TEST SED loopback calc', iloveclim_calc_lbtot

#    ifdef WITH_C13
      calc13_loopback(:,:) = iloveclim_calc13_lbtot / sfc_area_tot
#    else
      calc13_loopback(:,:) = 0.0D+00
#    endif

#    ifdef WITH_C14
      calc14_loopback(:,:) = iloveclim_calc14_lbtot / sfc_area_tot
#    else
      calc14_loopback(:,:) = 0.0D+00
#    endif

#  else

      calc_loopback(:,:)   = TRANSPOSE(iloveclim_calc_lb(:,:))(:,lbnd:ubnd)

#    ifdef WITH_C13
      calc13_loopback(:,:) = TRANSPOSE(iloveclim_calc13_lb(:,:))(:,lbnd:ubnd)
#    else
      calc14_loopback(:,:) = 0.0D+00
#    endif

#    ifdef WITH_C14
      calc14_loopback(:,:) = TRANSPOSE(iloveclim_calc14_lb(:,:))(:,lbnd:ubnd)
#    else
      calc14_loopback(:,:) = 0.0D+00
#    endif

#  endif

#else

      calc_loopback(:,:)   = 0.0D+00

      calc13_loopback(:,:) = 0.0D+00

      calc14_loopback(:,:) = 0.0D+00

#endif


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE xchange_fluxes_s2o
!-----------------------------------------------------------------------


!=======================================================================
      END MODULE MOD_ILOVECLIM_S2O
!=======================================================================
