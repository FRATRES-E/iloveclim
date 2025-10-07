!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"
!dmr -- Adding the choice of components through the pre-processing options

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!=======================================================================
      module flux_from_sediments_mod
!=======================================================================
#if ( MEDUSA > 0 )

      USE declars_mod, ONLY: LT, JT, NOC_CBR

      IMPLICIT NONE


      CONTAINS


!-----------------------------------------------------------------------
      SUBROUTINE flux_s2o_ini()
!-----------------------------------------------------------------------

      USE mbiota_mod, ONLY: ODOCS_sed2oc, ODIC_sed2oc, OALK_sed2oc,    &
                            OO2_sed2oc, ONO3_sed2oc, OPO4_sed2oc,      &
                            ODIC13_sed2oc, ODIC14_sed2oc
      USE mbiota_mod, ONLY: orgm_dic_loopback, orgm_alk_loopback,      &
                            orgm_oxyg_loopback,                        &
                            orgm_no3_loopback, orgm_po4_loopback,      &
                            calc_loopback,                             &
                            orgm_dic13_loopback, calc13_loopback,      &
                            orgm_dic14_loopback, calc14_loopback

      IMPLICIT NONE


#ifdef DEBUG
      print*, "flux_s2o_ini"
#endif

      ODOCS_sed2oc  = 0.0D+00
      ODIC_sed2oc   = 0.0D+00
      OALK_sed2oc   = 0.0D+00
      OO2_sed2oc    = 0.0D+00
      ONO3_sed2oc   = 0.0D+00
      OPO4_sed2oc   = 0.0D+00
      ODIC13_sed2oc = 0.0D+00
      ODIC14_sed2oc = 0.0D+00

      orgm_dic_loopback   = 0.0D+00
      orgm_alk_loopback   = 0.0D+00
      orgm_no3_loopback   = 0.0D+00
      orgm_po4_loopback   = 0.0D+00
      orgm_oxyg_loopback  = 0.0D+00

      orgm_dic13_loopback = 0.0D+00
      orgm_dic14_loopback = 0.0D+00

      calc_loopback       = 0.0D+00

      calc13_loopback     = 0.0D+00
      calc14_loopback     = 0.0D+00


      RETURN


!-----------------------------------------------------------------------
      END SUBROUTINE flux_s2o_ini
!-----------------------------------------------------------------------

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----------------------------------------------------------------------
      SUBROUTINE update_flx_from_sediments()
!-----------------------------------------------------------------------

      USE mbiota_mod, ONLY: ODOCS_sed2oc, ODIC_sed2oc, OALK_sed2oc,    &
                            OO2_sed2oc, ONO3_sed2oc, OPO4_sed2oc,      &
                            ODIC13_sed2oc, ODIC14_sed2oc

      USE mbiota_mod, ONLY: orgm_dic_loopback, orgm_alk_loopback,      &
                            orgm_no3_loopback, orgm_po4_loopback,      &
                            orgm_oxyg_loopback,                        &
                            calc_loopback,                             &
                            orgm_dic13_loopback, calc13_loopback,      &
                            orgm_dic14_loopback, calc14_loopback
      
      USE mod_iloveclim_s2o, ONLY: xchange_fluxes_s2o


      IMPLICIT NONE


      CALL flux_s2o_ini()
      
      CALL xchange_fluxes_s2o(ODOCS_sed2oc, ODIC_sed2oc, OALK_sed2oc, &
           OO2_sed2oc, ONO3_sed2oc, OPO4_sed2oc, &
           ODIC13_sed2oc, ODIC14_sed2oc, &
           orgm_dic_loopback, orgm_alk_loopback, &
           orgm_oxyg_loopback, orgm_no3_loopback, orgm_po4_loopback, &
           calc_loopback, &
           orgm_dic13_loopback, calc13_loopback, &
           orgm_dic14_loopback, calc14_loopback)            

      RETURN

!-----------------------------------------------------------------------
      END SUBROUTINE update_flx_from_sediments
!-----------------------------------------------------------------------

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----------------------------------------------------------------------
      SUBROUTINE flux_from_sediments()
!-----------------------------------------------------------------------

      USE declars_mod, ONLY: LT, JT, NOC_CBR
      USE loveclim_transfer_mod, ONLY: DVOL, MGT, OVOL
      USE loveclim_transfer_mod, ONLY: SQRO2, ZZ, total_area    !area, depth of cells
      USE mod_iloveclim_s2o, ONLY: xchange_fluxes_s2o

      USE marine_bio_mod, ONLY: oalk, odic, oo2, odocs, ono3, opo4
#ifdef WITH_C13
      USE marine_bio_mod, ONLY: odic13
#endif
#ifdef WITH_C14
      USE marine_bio_mod, ONLY: odic14
#endif

      USE mbiota_mod, ONLY: ODOCS_sed2oc, ODIC_sed2oc, OALK_sed2oc,    &
                            OO2_sed2oc, ONO3_sed2oc, OPO4_sed2oc,      &
                            ODIC13_sed2oc, ODIC14_sed2oc,              &
                            summary_flux_S2O_dic                       

      USE mbiota_mod, ONLY: orgm_dic_loopback, orgm_alk_loopback,      &
                            orgm_no3_loopback, orgm_po4_loopback,      &
                            orgm_oxyg_loopback,                        &
                            calc_loopback,                             &
                            orgm_dic13_loopback, calc13_loopback,      &
                            orgm_dic14_loopback, calc14_loopback

      USE mbiota_mod, ONLY: riverine_dic_input_tot, &
                            riverine_alk_input_tot, &
                            riverine_O2_input_tot, &
                            riverine_NO3_input_tot, &
                            riverine_PO4_input_tot 
#ifdef WITH_C13
      USE mbiota_mod, ONLY: riverine_dic13_input_tot
#endif
#ifdef WITH_C14
      USE mbiota_mod, ONLY: rriverine_dic14_input_tot
#endif


      USE mod_iloveclim_setup, ONLY: kfs_fond


      IMPLICIT NONE


      INTEGER :: i, n
      INTEGER :: klb                ! indox of layer for receiving loopback fluxes
      INTEGER :: kbot               ! vertical index of bottom layer

      DOUBLE PRECISION :: DMASS     ! water mass of the current cell [kg/grid_elt]
                                    ! DMASS = DVOL [m^3/grid_elt] * density [kg/m^3])

      DOUBLE PRECISION :: area      !surface

                                    ! time step length applied for the
                                    ! the integration of the sediment and
                                    ! weathering fluxes [yr] 
      DOUBLE PRECISION, PARAMETER :: dt_sedweaflx = 1.0d+00/360.0d+00

                                    ! seawater density [kg/m^3] provisional fix
                                    ! TO BE RETRIEVED FROMA GLOBAL FUNCTION
      !DOUBLE PRECISION, PARAMETER :: rho_sw = 1.028d+03
      DOUBLE PRECISION, PARAMETER :: rho_sw = 1.000d+03 ! compatibility value


! check precision with Guy
      REAL :: riverine_dic_input
      REAL :: riverine_alk_input
      REAL :: riverine_O2_input
      REAL :: riverine_NO3_input
      REAL :: riverine_PO4_input

#ifdef WITH_C13
      REAL :: riverine_dic13_input
#endif
#ifdef WITH_C14
      REAL :: riverine_dic14_input
#endif


      DO i = 1, LT
        DO n = 1, NOC_CBR

          IF (MGT(i,1,n) == 1) THEN

            kbot = INT(kfs_fond(i,n))
            DMASS = DVOL(i, kbot, n) * rho_sw

            ! All O<XYZ>_sed2oc fluxes (sediment-to-ocean fluxes across the SWI)
            !  - are in mol/grid_elt/yr
            !  - follow Medusa sign conventions, i.e., are positive *into*
            !    the sediment, and must therefore revrt sign when
            !    considered going into the ocean
            !
            ! DIC and ALK concentrations OCYCC are in in mol/kg:
            !          / DMASS [kg/grid_elt] :  grid_elt^-1 -> kg^-1
            !          * dt_sedweaflx: yr^-1 -> time_step^-1
            ! DOCS, PO4, NO3 and O2 are in µmol/kg
            !          * 1.0d+06 : mol -> µmol

            ODOCS(i, kbot, n) = ODOCS(i, kbot, n) &
              - ODOCS_sed2oc(i,n)/DMASS * 1.0d+06 * dt_sedweaflx

            ODIC(i, kbot, n) = ODIC(i, kbot, n) &
              - ODIC_sed2oc(i,n)/DMASS * dt_sedweaflx

            OALK(i, kbot, n) = OALK(i, kbot, n) &
              - OALK_sed2oc(i,n)/DMASS * dt_sedweaflx

            OO2(i, kbot, n) = OO2(i, kbot, n) &
              - OO2_sed2oc(i,n)/DMASS * 1.0d+06 * dt_sedweaflx

            ONO3(i, kbot, n) =   ONO3(i, kbot, n) &
              - ONO3_sed2oc(i,n)/DMASS * 1.0d+06 * dt_sedweaflx

            OPO4(i, kbot, n) =   OPO4(i, kbot, n) &
              - OPO4_sed2oc(i,n)/DMASS * 1.0d+06 * dt_sedweaflx

#ifdef WITH_C13
            ODIC13(i, kbot, n) =   ODIC13(i, kbot, n) &
              - ODIC13_sed2oc(i,n)/DMASS * dt_sedweaflx
#endif

#ifdef WITH_C14
            ODIC14(i, kbot, n) =   ODIC14(i, kbot ,n) &
              - ODIC14_sed2oc(i,n)/DMASS * dt_sedweaflx
#endif

          END IF

        END DO
      END DO


! Dealing with the riverine influxes at the surface

!#if ( sediment_loopback == 1 )
!      riverine_dic_input_tot = 0.0
!#endif

      DO i = 1, LT
        DO n = 1, NOC_CBR

          IF (MGT(i,1,n) == 1) THEN

            ! Units:
            !   orgm_dic_loopback [molC/m^2/yr]
            !   orgm_alk_loopback [eq/m^2/yr]
            !   orgm_dic_loopback [molC/m^2/yr]
            !   calc_loopback [molCaCO3/m^2/yr] = [molC/m^2/yr]

            ! Unlike O<XYZ>_sed2oc fluxes, <solid>_<solute>_loopback
            ! fluxes are in units of mol/m^2/yr and thus have to be
            ! multiplied by the grid-element surface area to use
            ! them per grid element.

            area = sqro2(i,n)


            ! klb: index of the injection layer for loopback fluxes,
            ! safe-guarded to lie within the water column.
            ! Please select by uncommenting.

            klb = 1  ! Topmost layer (first below surface)
            !klb = MIN(2, INT(kfs_fond(i,n))) ! Second layer below surface
            !klb = MAX(INT(kfs_fond(i,n))-1, 1) ! Second layer above the seafloor
            !klb = MAX(INT(kfs_fond(i,n)), 1) ! Bottom layer (first above sea-floor sediment)


            DMASS = DVOL(i, klb, n) * rho_sw

!if we use the loopback : fluxes are computed from what goes into sediments
#if ( sediment_loopback == 1 )
            riverine_dic_input=(orgm_dic_loopback(i,n) + calc_loopback(i,n))*area/DMASS
            !riverine_dic_input_tot=riverine_dic_input_tot+riverine_dic_input
            riverine_alk_input=(orgm_alk_loopback(i,n) + calc_loopback(i,n)*2.0D+00)*area/DMASS
            riverine_O2_input=orgm_oxyg_loopback(i,n)*area/DMASS*1.0d+06
            riverine_NO3_input=orgm_no3_loopback(i,n)*area/DMASS*1.0d+06
            riverine_PO4_input=orgm_po4_loopback(i,n)*area/DMASS*1.0d+06
#ifdef WITH_C13
            riverine_dic13_input=(orgm_dic13_loopback(i,n) + calc13_loopback(i,n))*area/DMASS
#endif
#ifdef WITH_C14
            riverine_dic14_input=(orgm_dic14_loopback(i,n) + calc14_loopback(i,n))*area/DMASS
#endif

!else the fluxes from rivers are fixed, with uniform repartition
#else
            !riverine_dic_input_total total
            riverine_dic_input=riverine_dic_input_tot*area/total_area
            riverine_alk_input=riverine_alk_input_tot*area/total_area
            riverine_O2_input=riverine_O2_input_tot*area/total_area
            riverine_NO3_input=riverine_NO3_input_tot*area/total_area
            riverine_PO4_input=riverine_PO4_input_tot*area/total_area
#ifdef WITH_C13
            riverine_dic13_input=riverine_dic13_input_tot*area/total_area
#endif
#ifdef WITH_C14
            riverine_dic14_input=riverine_dic14_input_tot*area/total_area
#endif

#endif

!tbd            !ODIC(i, klb, n) = ODIC(i, klb, n) &
!tbd            !  + (orgm_dic_loopback(i,n) + calc_loopback(i,n))*area/DMASS * dt_sedweaflx
            ODIC(i, klb, n) = ODIC(i, klb, n) &
              +  riverine_dic_input* dt_sedweaflx
!tbd
            !OALK(i, klb, n) =   OALK(i, klb, n) &
!tbd            !  + (orgm_alk_loopback(i,n) + calc_loopback(i,n)*2.0D+00)*area/DMASS * dt_sedweaflx
            OALK(i, klb, n) =   OALK(i, klb, n) &
              + riverine_alk_input * dt_sedweaflx

!tbd            OO2(i, klb, n) = OO2(i, klb, n) &
!tbd              + orgm_oxyg_loopback(i,n)*area/DMASS*1.0d+06 * dt_sedweaflx
            OO2(i, klb, n) = OO2(i, klb, n) &
              + riverine_O2_input * dt_sedweaflx

!tbd            ONO3(i, klb, n) = ONO3(i, klb, n) &
!tbd              + orgm_no3_loopback(i,n)*area/DMASS*1.0d+06 * dt_sedweaflx
            ONO3(i, klb, n) = ONO3(i, klb, n) &
              + riverine_NO3_input * dt_sedweaflx

!tbd            OPO4(i, klb, n) =  OPO4(i, klb, n) &
!tbd              + orgm_po4_loopback(i,n)*area/DMASS*1.0d+06 * dt_sedweaflx
            OPO4(i, klb, n) =  OPO4(i, klb, n) &
              + riverine_PO4_input * dt_sedweaflx

#ifdef WITH_C13
!tbd            ODIC13(i, klb, n) = ODIC13(i, klb, n) &
!tbd              + (orgm_dic13_loopback(i,n) + calc13_loopback(i,n))*area/DMASS * dt_sedweaflx
            ODIC13(i, klb, n) = ODIC13(i, klb, n) &
              + riverine_dic13_input * dt_sedweaflx
#endif

#ifdef WITH_C14
!tbd            ODIC14(i, klb, n) = ODIC14(i, klb, n) &
!tbd              + (orgm_dic14_loopback(i,n) + calc14_loopback(i,n))*area/DMASS * dt_sedweaflx
            ODIC14(i, klb, n) = ODIC14(i, klb, n) &
              + riverine_dic14_input * dt_sedweaflx
#endif

          END IF

        END DO
      END DO

      !write(*,*) 'riverine_dic_input_tot', riverine_dic_input_tot


      ! ??? CALL flux_s2o_ini()


      summary_flux_S2O_dic = SUM(ODIC_sed2oc(:,:))/1.D12



      RETURN

!-----------------------------------------------------------------------
      END SUBROUTINE flux_from_sediments
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE restart_flx_from_sediments(choix)
!-----------------------------------------------------------------------
!nb to write and read restart fluxes from sediments

      USE mbiota_mod, ONLY: ODOCS_sed2oc, ODIC_sed2oc, OALK_sed2oc,    &
                            OO2_sed2oc, ONO3_sed2oc, OPO4_sed2oc,      &
                            ODIC13_sed2oc, ODIC14_sed2oc

      USE mbiota_mod, ONLY: orgm_dic_loopback, orgm_alk_loopback,      &
                            orgm_no3_loopback, orgm_po4_loopback,      &
                            orgm_oxyg_loopback,                        &
                            calc_loopback,                             &
                            orgm_dic13_loopback, calc13_loopback,      &
                            orgm_dic14_loopback, calc14_loopback

      USE mbiota_mod, ONLY: riverine_dic_input_tot, &
                            riverine_alk_input_tot, &
                            riverine_O2_input_tot, &
                            riverine_NO3_input_tot, &
                            riverine_PO4_input_tot 
#ifdef WITH_C13
      USE mbiota_mod, ONLY: riverine_dic13_input_tot
#endif
#ifdef WITH_C14
      USE mbiota_mod, ONLY: rriverine_dic14_input_tot
#endif

      USE loveclim_transfer_mod, ONLY: DVOL, MGT, SQRO2
      USE declars_mod, ONLY: LT, NOC_CBR


      IMPLICIT NONE

      INTEGER :: choix

      INTEGER :: fich_num, nrecl
      CHARACTER*12, PARAMETER :: fich_res_name="rest_sed.dat"
      CHARACTER*22, PARAMETER ::                                       &
                              fich_res_name_old="startdata/rest_sed.dat"

      DOUBLE PRECISION, PARAMETER :: rho_sw = 1.000d+03 ! compatibility value
      INTEGER :: i,n, klb
      DOUBLE PRECISION :: area      !surface
      DOUBLE PRECISION :: DMASS     ! water mass of the current cell [kg/grid_elt]
                                    ! DMASS = DVOL [m^3/grid_elt] * density [kg/m^3])

! check precision with Guy
      REAL :: riverine_dic_input
      REAL :: riverine_alk_input
      REAL :: riverine_O2_input
      REAL :: riverine_NO3_input
      REAL :: riverine_PO4_input

#ifdef WITH_C13
      REAL :: riverine_dic13_input
#endif
#ifdef WITH_C14
      REAL :: riverine_dic14_input
#endif

      LOGICAL :: existe
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Determination d'un numero de fichier libre ...
!-----|--1--------2---------3---------4---------5---------6---------7-|

       fich_num=298
       print*, "init_mb_fich", fich_num
       existe=.TRUE.

       DO WHILE (existe)
        INQUIRE(fich_num,OPENED=existe)
        print*, "init_mb_fich", fich_num, existe
        fich_num=fich_num+1 
       END DO 


!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Cas ou l'on ecrit le restart
!-----|--1--------2---------3---------4---------5---------6---------7-|
       IF (choix.EQ.0) THEN


!       Determination de la taille de l'ecriture
!       ----------------------------------------
       nrecl = 8*SIZE(ODOCS_sed2oc)+10*SIZE(orgm_dic_loopback)
       nrecl = nrecl * KIND(ODOCS_sed2oc)

!      write restart
!      -------------
        OPEN(UNIT=fich_num, FILE=fich_res_name, STATUS='unknown',       &
              ACCESS='direct', RECL=nrecl, ACTION='write')

        WRITE(UNIT=fich_num, REC=1)                                    &
                            ODOCS_sed2oc, ODIC_sed2oc, OALK_sed2oc,    &
                            OO2_sed2oc, ONO3_sed2oc, OPO4_sed2oc,      &
                            ODIC13_sed2oc, ODIC14_sed2oc,              &
                            orgm_dic_loopback, orgm_alk_loopback,      &
                            orgm_no3_loopback, orgm_po4_loopback,      &
                            orgm_oxyg_loopback,                        &
                            calc_loopback,                             &
                            orgm_dic13_loopback, calc13_loopback,      &
                            orgm_dic14_loopback, calc14_loopback

        CLOSE(UNIT=fich_num)



!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Cas ou l'on lit le restart
!-----|--1--------2---------3---------4---------5---------6---------7-|
       ELSE IF (choix.EQ.1) THEN

!       Determination de la taille de l'ecriture
!       ----------------------------------------
       nrecl = 8*SIZE(ODOCS_sed2oc)+10*SIZE(orgm_dic_loopback)
       nrecl = nrecl * KIND(ODOCS_sed2oc)

!     read restart
!     ------------
      WRITE(*,*) "Lecture du redemarrage carbone a partir du fichier : "
      WRITE(*,*) fich_res_name_old

        OPEN(UNIT=fich_num, FILE=fich_res_name_old, STATUS='old',       &
              ACCESS='direct', RECL=nrecl, ACTION='read')

        READ(UNIT=fich_num,REC=1)                                      &
                            ODOCS_sed2oc, ODIC_sed2oc, OALK_sed2oc,    &
                            OO2_sed2oc, ONO3_sed2oc, OPO4_sed2oc,      &
                            ODIC13_sed2oc, ODIC14_sed2oc,              &
                            orgm_dic_loopback, orgm_alk_loopback,      &
                            orgm_no3_loopback, orgm_po4_loopback,      &
                            orgm_oxyg_loopback,                        &
                            calc_loopback,                             &
                            orgm_dic13_loopback, calc13_loopback,      &
                            orgm_dic14_loopback, calc14_loopback
                     
        CLOSE(UNIT=fich_num)

      !check this is the same as in flux_from_sediments
      klb = 1  ! Topmost layer (first below surface)

      !computes total fluxes from rivers
      !---------------------------------
      !write(*,*) 'restart before riverine_dic_input_tot', riverine_dic_input_tot
      riverine_dic_input_tot = 0.0
      riverine_alk_input_tot = 0.0
      riverine_O2_input_tot = 0.0
      riverine_NO3_input_tot = 0.0
      riverine_PO4_input_tot = 0.0
#ifdef WITH_C13
            riverine_dic13_input_tot=0.0
#endif
#ifdef WITH_C14
            riverine_dic14_input_tot=0.0
#endif

      DO i = 1, LT
        DO n = 1, NOC_CBR

          IF (MGT(i,1,n) == 1) THEN
            area = sqro2(i,n)
            DMASS = DVOL(i, klb, n) * rho_sw
            riverine_dic_input=(orgm_dic_loopback(i,n) + calc_loopback(i,n))*area/DMASS
            riverine_dic_input_tot=riverine_dic_input_tot+riverine_dic_input

            riverine_alk_input=(orgm_alk_loopback(i,n) + calc_loopback(i,n)*2.0D+00)*area/DMASS
            riverine_alk_input_tot=riverine_alk_input_tot+riverine_alk_input

            riverine_O2_input=orgm_oxyg_loopback(i,n)*area/DMASS*1.0d+06
            riverine_O2_input_tot=riverine_O2_input_tot+riverine_O2_input

            riverine_NO3_input=orgm_no3_loopback(i,n)*area/DMASS*1.0d+06
            riverine_NO3_input_tot=riverine_NO3_input_tot+riverine_NO3_input

            riverine_PO4_input=orgm_po4_loopback(i,n)*area/DMASS*1.0d+06
            riverine_PO4_input_tot=riverine_PO4_input_tot+riverine_PO4_input

#ifdef WITH_C13
            riverine_dic13_input=(orgm_dic13_loopback(i,n) + calc13_loopback(i,n))*area/DMASS
            riverine_dic13_input_tot=riverine_dic13_input_tot+riverine_dic13_input
#endif
#ifdef WITH_C14
            riverine_dic14_input=(orgm_dic14_loopback(i,n) + calc14_loopback(i,n))*area/DMASS
            riverine_dic14_input_tot=riverine_dic14_input_tot+riverine_dic14_input
#endif

           ENDIF
         ENDDO
      ENDDO

     !write(*,*) 'restart after riverine_dic_input_tot', riverine_dic_input_tot
     !write(*,*) 'riverine_alk_input_tot', riverine_alk_input_tot
     !write(*,*) 'riverine_O2_input_tot', riverine_O2_input_tot
     !write(*,*) 'riverine_NO3_input_tot', riverine_NO3_input_tot
     !write(*,*) 'riverine_PO4_input_tot', riverine_PO4_input_tot


      ENDIF

      RETURN

!-----------------------------------------------------------------------
      END SUBROUTINE restart_flx_from_sediments
!-----------------------------------------------------------------------



#endif

!=======================================================================
      END MODULE flux_from_sediments_mod
!=======================================================================

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
