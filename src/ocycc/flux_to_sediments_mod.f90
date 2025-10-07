!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!=======================================================================
       MODULE flux_to_sediments_mod
!=======================================================================

#if ( MEDUSA > 0 )

       USE mod_iloveclim_setup, ONLY: kfs_fond

       IMPLICIT NONE

       CONTAINS

!-----------------------------------------------------------------------
       SUBROUTINE fait_pointer_O2S
!-----------------------------------------------------------------------

! Guy Munhoven (20 October 2016)
! Please notice that
! (1) nothing is required here for NO3 and PO4, since the sediment-to-ocean
!     return fluxes of these two are directly derived from the oxygen
!     consumption flux, and do not depend on the concentration at the
!     boundary;
! (2) nothing is required here for oxygen isotopes, as they are not
!     influenced by biogeochemistry;
! (3) although the cabron isotope boundary conditions are included, they
!     could possibly be left out to simplify , and all the isotopic components
!     released nto porewaters returned to the water column above.


       USE loveclim_transfer_mod,   ONLY: tm,      sm
       USE mbiota_mod,              ONLY: temp_ma, salt_ma

       USE marine_bio_mod,          ONLY: odic,    oalk,    oo2
       USE mbiota_mod,              ONLY: odic_ma, oalk_ma, ooxy_ma

#ifdef WITH_C13
       USE marine_bio_mod,          ONLY: oc13    ! [??? GM] -- > the names are correct
       USE mbiota_mod,              ONLY: oc13_ma ! or the names which correspond
#endif

#ifdef WITH_C14
       USE marine_bio_mod,          ONLY: oc14    ! [??? GM]
       USE mbiota_mod,              ONLY: oc14_ma ! or the names which correspond
#endif


       IMPLICIT NONE


       temp_ma => tm
       salt_ma => sm
       odic_ma => odic
       oalk_ma => oalk
!nb 
       ooxy_ma => oo2(:,:,:)

#ifdef WITH_C13
       oc13_ma => oc13    ! [??? GM] -- > correct!
#endif

#ifdef WITH_C14
       oc14_ma => oc14    ! [??? GM] --> correct!
#endif

       !INCLUDING NO3 and PO4 in sediments
       !ono3_ma => ono3
       !opo4_ma => opo4

       RETURN


!-----------------------------------------------------------------------
       END SUBROUTINE fait_pointer_O2S
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
       SUBROUTINE sedfluxave(med_timer, timestep) !nb sediments
!-----------------------------------------------------------------------
!-----|--1---------2---------3---------4---------5---------6---------7-|
!
!
!-----|--1---------2---------3---------4---------5---------6---------7-|

      USE declars_mod, ONLY: LT, JT, NOC_CBR
      USE loveclim_transfer_mod, ONLY: MGT, mid_level, SQRO2

      USE mod_iloveclim_setup, ONLY: kfs_fond
      USE mod_iloveclim_o2s,   ONLY: xchange_fluxes_o2s


      USE mbiota_mod, ONLY: temp_ma, temp_mave, temp_mafond
      USE mbiota_mod, ONLY: salt_ma, salt_mave, salt_mafond

      USE mbiota_mod, ONLY: ooxy_ma, ooxy_mave, ooxy_mafond
      USE mbiota_mod, ONLY: oalk_ma, oalk_mave, ooxy_mave
      USE mbiota_mod, ONLY: odic_ma, odic_mave
      USE mbiota_mod, ONLY:                     oco2_mafond, &
                                                oco3_mafond, &
                                                ohco3_mafond

      USE mbiota_mod, ONLY:  clay_ma,  clay_mave,  clay_mafond
      USE mbiota_mod, ONLY:   TPP_ma,   TPP_mave,   TPP_mafond
      USE mbiota_mod, ONLY: caco3_ma, caco3_mave, caco3_mafond

      USE mbiota_mod, ONLY: ono3_mafond, opo4_mafond
      USE mbiota_mod, ONLY: tracer01_ma, tracer01_mave, tracer01_mafond
      USE mbiota_mod, ONLY: caco3_mabot


#ifdef WITH_C13
                                    ! *c13_ma and *c13_mave are only required
                                    ! if C13 is traced
      USE mbiota_mod, ONLY:     oc13_ma,     oc13_mave
      USE mbiota_mod, ONLY:   TPPC13_ma,   TPPC13_mave
      USE mbiota_mod, ONLY: caco3C13_ma, caco3C13_mave
#endif
                                    ! *c13_mafond are always required
      USE mbiota_mod, ONLY:                                  oc13_mafond
      USE mbiota_mod, ONLY:                                TPPC13_mafond
      USE mbiota_mod, ONLY:                              caco3C13_mafond


#ifdef WITH_C14
                                    ! *c14_ma and *c14_mave are only required
                                    ! if C14 is traced
      USE mbiota_mod, ONLY:     oc14_ma,     oc14_mave
      USE mbiota_mod, ONLY:   TPPC14_ma,   TPPC14_mave
      USE mbiota_mod, ONLY: caco3C14_ma, caco3C14_mave
#endif
                                    ! *c14_mafond is always required
      USE mbiota_mod, ONLY:                                  oc14_mafond
      USE mbiota_mod, ONLY:                                TPPC14_mafond
      USE mbiota_mod, ONLY:                              caco3C14_mafond


      USE mbiota_mod, ONLY: summary_flux_O2S_calc, summary_flux_O2S_orgm

!nb 
      use const_mod, only: gpes, rho0


#ifdef WITH_O18
                                    ! *O18_ma and *O18_mave are only required
                                    ! if O18 is traced
      USE mbiota_mod, ONLY: caco3O18_ma, caco3O18_mave
#endif
                                    ! *O18_mafond is always required
      USE mbiota_mod, ONLY:                              caco3O18_mafond
      use carbonate_speciation_mod, only: incche


      IMPLICIT NONE


      INTEGER, INTENT(IN) :: med_timer
      INTEGER, INTENT(IN) :: timestep ! length of medusa time step in days

      INTEGER  :: i, j, n, kyr, k
      REAL     :: sCO2, xpCO2, xCO2, xHCO3, xCO3, z_h
      REAL     :: p_bar

      ! The clay flux would normally be set together with the other _ma
      ! fluxes. Temporarily it is set here to constant value
!      clay_ma(:,:,:)  = 2.0D-03 ! kg m**-2 yr**-1
!nb and gm modified clay flux, with z depth = -mid_level(j)
      !DO n = 1, NOC_CBR
      !  DO i = 1, LT
          DO j = 1, JT
      clay_ma(:,j,:)  = 2.0D-03 + 10.0**                                &
         (-2.4+mid_level(j)/1250.)*0.1*2650 ! kg m**-2 yr**-1
       !write(*,*) 'clay',j, mid_level(j), clay_ma(1,j,8)
          ENDDO
       !write(*,*) 'clay', clay_ma(1,:,8)

!      print*, "checking timestep in flux_to_sediments...", timestep, med_timer !!-- correct! 

      DO n = 1, NOC_CBR
        DO i = 1, LT
          DO j = 1, JT

            temp_mave(i,j,n) =  temp_mave(i,j,n) +  temp_ma(i,j,n) ! replacing the previous value when
            salt_mave(i,j,n) =  salt_mave(i,j,n) +  salt_ma(i,j,n) ! mod(NYR,NYRave) is not zero

            odic_mave(i,j,n) =  odic_mave(i,j,n) +  odic_ma(i,j,n)
            oalk_mave(i,j,n) =  oalk_mave(i,j,n) +  oalk_ma(i,j,n)
            ooxy_mave(i,j,n)  =  ooxy_mave(i,j,n) +  ooxy_ma(i,j,n)

#ifdef WITH_C13
            oc13_mave(i,j,n)     =     oc13_mave(i,j,n) +     oc13_ma(i,j,n)
            TPPC13_mave(i,j,n)   =   TPPC13_mave(i,j,n) +   TPPC13_ma(i,j,n) !this variable may exist already
            caco3C13_mave(i,j,n) = caco3C13_mave(i,j,n) + caco3C13_ma(i,j,n) !with different name
#endif

#ifdef WITH_C14
            oc14_mave(i,j,n)     =     oc14_mave(i,j,n) +     oc14_ma(i,j,n)
            TPPC14_mave(i,j,n)   =   TPPC14_mave(i,j,n) +   TPPC14_ma(i,j,n) !this variable may exist already
            caco3C14_mave(i,j,n) = caco3C14_mave(i,j,n) + caco3C14_ma(i,j,n) !with different name
#endif

#ifdef WITH_O18
            caco3O18_mave(i,j,n) = caco3O18_mave(i,j,n) + caco3O18_ma(i,j,n) !this variable does not yet exist
#endif

          ENDDO
        ENDDO
      ENDDO

!mohr

!             averaging


      IF (med_timer == 0) THEN

        write(*,*) "timestep in MEDUSA ave ===", timestep
        
            !! variables are already accumulated per timestep in MAPHOT.f
            !! i.e. in maphot, caco3_ma = caco3_ma + new_caco3, with new_caco3 in Tmol.m-2.day-1
            
            clay_mave(:,:,:) =   clay_ma(:,:,:)
            TPP_mave(:,:,:)  =   TPP_ma(:,:,:)   ! SUM(TPP,t=1,nb_timestep) -> [TPP_ma]
            caco3_mave(:,:,:)=   caco3_ma(:,:,:) ! this saves the value every timestep days
        
        
        
! --- dmr                   [TODO] check thoroughly the accumulation of fluxes below for values ...
        DO i = 1, LT
          DO j = 1, JT
            DO n = 1, NOC_CBR

              temp_mave(i,j,n)     =     temp_mave(i,j,n)/timestep
              salt_mave(i,j,n)     =     salt_mave(i,j,n)/timestep

              odic_mave(i,j,n)     =     odic_mave(i,j,n)/timestep
              oalk_mave(i,j,n)     =     oalk_mave(i,j,n)/timestep
              ooxy_mave(i,j,n)     =     ooxy_mave(i,j,n)/timestep

              clay_mave(i,j,n)     =     clay_mave(i,j,n) 
              TPP_mave(i,j,n)      =      TPP_mave(i,j,n) ! [TPP_MAVE]   -> Tmols.m-2.timestep-1
              caco3_mave(i,j,n)    =    caco3_mave(i,j,n) ! [CACO3_MAVE] -> Tmols.m-2.timestep-1

#ifdef WITH_C13
              oc13_mave(i,j,n)     =     oc13_mave(i,j,n)/timestep
              TPPC13_mave(i,j,n)   =   TPPC13_mave(i,j,n)/timestep
              caco3C13_mave(i,j,n) = caco3C13_mave(i,j,n)/timestep
#endif

#ifdef WITH_C14
              oc14_mave(i,j,n)     =     oc14_mave(i,j,n)/timestep
              TPPC14_mave(i,j,n)   =   TPPC14_mave(i,j,n)/timestep
              caco3C14_mave(i,j,n) = caco3C14_mave(i,j,n)/timestep
#endif

#ifdef WITH_O18
              caco3O18_mave(i,j,n) = caco3O18_mave(i,j,n)/timestep
#endif

            ENDDO
          ENDDO
        ENDDO


        DO i = 1, LT
          DO n = 1, NOC_CBR

            IF (MGT(i,1,n) == 1) THEN   ! si OCEAN en surface !!

              temp_mafond(i,n) = temp_mave(i,INT(kfs_fond(i,n)),n)
              salt_mafond(i,n) = salt_mave(i,INT(kfs_fond(i,n)),n)


              IF ( salt_mafond(i,n) < 1.E-3 ) THEN
                WRITE(*,*) "sal_fond === ", salt_mafond(i,n), i,n, kfs_fond(i,n)
                STOP
              ENDIF

        
!!              call set_values_constant()

!nb added new routine from Guy
              ! pressure in bar, rho0 should be replaced by rho
!dmr --- [NOTA] j is not defined here, so j+1 is 22 => overflow
!    ---        I think it should be replaced by INT(kfs_fond(i,n))

! ---              p_bar= rho0 * gpes * (-1) * mid_level(j+1)*1e-5
              p_bar= rho0 * gpes * (-1) * mid_level(INT(kfs_fond(i,n)))*1e-5

              !CALL incche(temp_mafond(i,n),salt_mafond(i,n),           &
              CALL incche(temp_mafond(i,n)+273.15,salt_mafond(i,n),    &
                      p_bar,                                           &
                      odic_mave(i,INT(kfs_fond(i,n)),n),               &
                      oalk_mave(i,INT(kfs_fond(i,n)),n),               &
                      sCO2, xpCO2, xCO2, xHCO3, xCO3, z_h)

              ooxy_mafond(i,n)  =  ooxy_mave(i,INT(kfs_fond(i,n)),n)
              oco2_mafond(i,n)  =  xCO2
              oco3_mafond(i,n)  =  xCO3
              ohco3_mafond(i,n) = xHCO3

!~               IF ( INT(kfs_fond(i,n)).GT.6 ) then
!~                  WRITE(*,*) "Happens in i,n", caco3_mave(i,int(kfs_fond(i,n)),n)
!~                  READ(*,*)
!~               ENDIF
                 

              !! FLUXES of TPP_mafond and caco3_mafond are per year 
              TPP_mafond(i,n)   =   TPP_mave(i,INT(kfs_fond(i,n)),n)
              caco3_mafond(i,n) = caco3_mave(i,int(kfs_fond(i,n)),n)

              clay_mafond(i,n)  =  clay_mave(i,int(kfs_fond(i,n)),n)

              ono3_mafond(i,n)  = 0.0D0 !NO3 and PO4 are not considered specifically
              opo4_mafond(i,n)  = 0.0D0 !in the OCYCC model. They are set to ZERO

#ifdef WITH_C13
              TPPC13_mafond(i,n)   =   TPPC13_mave(i,INT(kfs_fond(i,n)),n)
              caco3C13_mafond(i,n) = caco3C13_mave(i,kfs_fond(i,n),n)
#else
              TPPC13_mafond(i,n)   =  0.0D+00
              caco3C13_mafond(i,n) =  0.0D+00
#endif

#ifdef WITH_C14
              TPPC14_mafond(i,n)   =   TPPC14_mave(i,INT(kfs_fond(i,n)),n)
              caco3C14_mafond(i,n) = caco3C14_mave(i,kfs_fond(i,n),n)
#else
              TPPC14_mafond(i,n)   =  0.0D+00
              caco3C14_mafond(i,n) =  0.0D+00
#endif

#ifdef WITH_O18
              caco3O18_mafond(i,n) = caco3O18_mave(i,kfs_fond(i,n),n)
#else
              caco3O18_mafond(i,n) =  0.0D+00
#endif


            ELSE

              temp_mafond(i,n) =  0.0D0
              salt_mafond(i,n) =  0.0D0
              ooxy_mafond(i,n)  = 0.0D0
              oco2_mafond(i,n)  = 0.0D0
              oco3_mafond(i,n)  = 0.0D0
              ohco3_mafond(i,n) = 0.0D0
              TPP_mafond(i,n)   = 0.0D0
              caco3_mafond(i,n) = 0.0D0
              clay_mafond(i,n)  = 0.0D0
              ono3_mafond(i,n)  =   0.0D0  
              opo4_mafond(i,n)  =   0.0D0    
       
#ifdef WITH_C13
              TPPC13_mafond(i,n)   =  0.0D+00
              caco3C13_mafond(i,n) =  0.0D+00
#endif

#ifdef WITH_C14
              TPPC14_mafond(i,n)   =  0.0D+00
              caco3C14_mafond(i,n) =  0.0D+00
#endif

#ifdef WITH_O18
              caco3O18_mafond(i,n) =  0.0D+00
#endif


            END IF

          END DO
        END DO


!!      !! turning all fluxes to zero for debugging               
!!        temp_mafond  =  0.0
!!        salt_mafond  =  0.0
!!        oco2_mafond  =  0.0
!!        oco3_mafond  =  0.0
!!        ohco3_mafond =  0.0
!!        ooxy_mafond  =  0.0
!!        ono3_mafond  =  0.0
!!        opo4_mafond  =  0.0
!!        clay_mafond  =  0.0
!!        TPP_mafond   =  0.0
!!        caco3_mafond =  0.0
!!        oc13_mafond  =  0.0 
!!        TPPC13_mafond   = 0.0
!!        caco3C13_mafond = 0.0
!!        oc14_mafond     = 0.0
!!        TPPC14_mafond   = 0.0
!!        caco3C14_mafond = 0.0
!!        caco3O18_mafond = 0.0
!!        tracer01_mafond = 0.0

!!        print*, "set fluxes to sediments constant"
!!        call set_values_constant()
!!        print*, "set fluxes to sediments to zero"
!!        call  flux_o2s_ini()


!!        print*, "checking fluxes to sediments..."
!!        print*, "temperature...  ", temp_mafond(33:36,103:106)
!!        print*, "salinity...     ", salt_mafond(33:36,103:106)
!!        print*, "nit / phos ...  ", ono3_mafond(33:36,103:106), opo4_mafond(33:36,103:106)
!!        print*, "oxygen...       ", ooxy_mafond(33:36,103:106)
!!        print*, "oco2...         ", oco2_mafond(33:36,103:106)
!!        print*, "ohco3...        ", ohco3_mafond(33:36,103:106)
!!        print*, "oco3...         ", oco3_mafond(33:36,103:106)
!!        print*, "Clay...         ", clay_mafond(33:36,103:106)
!!        print*, "OrgMatt...      ", TPP_mafond(33:36,103:106)

        summary_flux_O2S_calc = SUM(caco3_mafond(:,:)*SQRO2(:,:))
        summary_flux_O2S_orgm =   SUM(TPP_mafond(:,:)*SQRO2(:,:))      

        CALL xchange_fluxes_o2s (temp_mafond, salt_mafond, &
             oco2_mafond, oco3_mafond, ohco3_mafond, &
             ooxy_mafond, ono3_mafond, opo4_mafond, &
             clay_mafond, TPP_mafond, caco3_mafond, &
             oc13_mafond, TPPC13_mafond, caco3C13_mafond, &
             oc14_mafond, TPPC14_mafond, caco3C14_mafond, &
             caco3O18_mafond, tracer01_mafond)
             
        CALL flux_o2s_ini() !after averaging we set values to zero again
        CALL reset_fluxes_o2s_tozero

      ENDIF ! on med_timer == 0

      RETURN

!-----------------------------------------------------------------------
      END SUBROUTINE sedfluxave
!-----------------------------------------------------------------------




!-----------------------------------------------------------------------
       SUBROUTINE flux_o2s_ini() !nb sediments
!-----------------------------------------------------------------------

       USE declars_mod, ONLY: LT, JT, NOC_CBR

       USE mbiota_mod, ONLY: temp_mave, salt_mave
       USE mbiota_mod, ONLY: odic_mave, oalk_mave, ooxy_mave
       USE mbiota_mod, ONLY: ono3_mave, opo4_mave
       USE mbiota_mod, ONLY: tracer01_mave, clay_mave, TPP_mave, caco3_mave
       USE mbiota_mod, ONLY: caco3_mabot

#ifdef WITH_C13
       USE mbiota_mod, ONLY: oc13_mave, TPPC13_mave, caco3C13_mave
#endif

#ifdef WITH_C14
       USE mbiota_mod, ONLY: oc14_mave, TPPC14_mave, caco3C14_mave
#endif

#ifdef WITH_O18
       USE mbiota_mod, ONLY: caco3O18_mave
#endif

       USE mbiota_mod
       USE mod_iloveclim_setup, ONLY: kfs_fond


       IMPLICIT NONE

!-----|--1--------2---------3---------4---------5---------6---------7-|
       INTEGER n, i, j
!-----|--1--------2---------3---------4---------5---------6---------7-|

       FORALL (n=1:NOC_CBR, i=1:LT, j=1:JT)

         temp_mave(i,j,n)     = 0.0D+00
         salt_mave(i,j,n)     = 0.0D+00

         odic_mave(i,j,n)     = 0.0D+00
         oalk_mave(i,j,n)     = 0.0D+00
         ooxy_mave(i,j,n)     = 0.0D+00

         ono3_mave(i,j,n)     = 0.0D+00
         opo4_mave(i,j,n)     = 0.0D+00

         clay_mave(i,j,n)     = 0.0D+00
         TPP_mave(i,j,n)      = 0.0D+00
         caco3_mave(i,j,n)    = 0.0D+00

#ifdef WITH_C13
         oc13_mave(i,j,n)     = 0.0D+00
         TPPC13_mave(i,j,n)   = 0.0D+00
         caco3C13_mave(i,j,n) = 0.0D+00
#endif

#ifdef WITH_C14
         oc14_mave(i,j,n)     = 0.0D+00
         TPPC14_mave(i,j,n)   = 0.0D+00
         caco3C14_mave(i,j,n) = 0.0D+00
#endif

#ifdef WITH_O18
         caco3O18_mave(i,j,n) = 0.0D+00
#endif

       END FORALL
       
       FORALL (n=1:NOC_CBR, i=1:LT)

!dmr --- Not the case for this version where we do not accumulate fluxes
!         at top
!         TPP_mas(i,n)=0
!         caco3_mas(i,n)=0

         !  temp_ma(i,j,int(kfs_fond(i,n))) =  0.0D0
         !  salt_ma (i,j,int(kfs_fond(i,n))) =  0.0D0
         !  ooxy_ma(i,j,int(kfs_fond(i,n)))  = 0.0D0

         !TPP_ma(i,int(kfs_fond(i,n)),n)   = 0.0D0
         !caco3_ma(i,int(kfs_fond(i,n)),n) = 0.0D0

        ! TPP_ma = 0.0D+00
        ! caco3_ma = 0.0D+00
         
         temp_mafond(i,n) =  0.0D0
         salt_mafond(i,n) =  0.0D0
         ooxy_mafond(i,n)  = 0.0D0
         oco2_mafond(i,n)  = 0.0D0
         oco3_mafond(i,n)  = 0.0D0
         ohco3_mafond(i,n) = 0.0D0
         !INCLUDING NO3 and PO4 in the sediment     
         ono3_mafond(i,n)  = 0.0D0
         opo4_mafond(i,n) = 0.0D0


         TPP_mafond(i,n)   = 0.0D0
         caco3_mafond(i,n) = 0.0D0
         caco3_mabot(i,n) = 0.0d0


       END FORALL

       RETURN

!-----------------------------------------------------------------------
       END SUBROUTINE flux_o2s_ini
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------                                                                              
       SUBROUTINE reset_fluxes_o2s_tozero !nb sediments                                                                                        
!-----------------------------------------------------------------                                                                               

       USE declars_mod, ONLY: LT, JT, NOC_CBR

!       USE mbiota_mod, ONLY: temp_mave, salt_mave                                                                                                   
!       USE mbiota_mod, ONLY: odic_mave, oalk_mave, ooxy_mave                                                                                        
!       USE mbiota_mod, ONLY: tracer01_mave, clay_mave, TPP_mave, caco3_mave                                                                         

       USE mbiota_mod
       USE mod_iloveclim_setup, ONLY: kfs_fond

       IMPLICIT NONE

!-----|--1--------2---------3---------4---------5---------6---------7-|                                                                              
       INTEGER n, i, j
!-----|--1--------2---------3---------4---------5---------6---------7-|                                                                              
       print*, "initalizing sediment values to ZERO, 2015-10-02"!mohr                                                                                


       FORALL (n=1:NOC_CBR, i=1:LT, j=1:JT)

         temp_mave(i,j,n) = 0.0D+00
         salt_mave(i,j,n) = 0.0D+00

         odic_mave(i,j,n) = 0.0D+00
         oalk_mave(i,j,n) = 0.0D+00
         ooxy_mave(i,j,n) = 0.0D+00

         ono3_mave(i,j,n) = 0.0D+00
         opo4_mave(i,j,n) = 0.0D+00

         clay_mave(i,j,n) = 0.0D+00
         clay_mave(i,j,n) = 0.0D+00
         TPP_mave(i,j,n)  = 0.0D+00
         caco3_mave(i,j,n)= 0.0D+00
         !  temp_ma(i,j,int(kfs_fond(i,n))) =  0.0D0                                                                                                 
         !  salt_ma (i,j,int(kfs_fond(i,n))) =  0.0D0                                                                                                 
         !  ooxy_ma(i,j,int(kfs_fond(i,n)))  = 0.0D0                                                                                                 

       END FORALL

         TPP_ma(:,:,:)   = 0.0D0
         caco3_ma(:,:,:) = 0.0D0

       FORALL (n=1:NOC_CBR, i=1:LT)

         temp_mafond(i,n) =  0.0D0
         salt_mafond(i,n) =  0.0D0
         ooxy_mafond(i,n)  = 0.0D0
         oco2_mafond(i,n)  = 0.0D0
         oco3_mafond(i,n)  = 0.0D0
         ohco3_mafond(i,n) = 0.0D0
         !INCLUDING NO3 and PO4 in the sediment                                                                                                    
         ono3_mafond(i,n)  = 0.0D0
         opo4_mafond(i,n) = 0.0D0
         TPP_mafond(i,n)   = 0.0D0
         caco3_mafond(i,n) = 0.0D0

         
       END FORALL

       RETURN
!-----------------------------------------------------------------------                                                                           
       END SUBROUTINE reset_fluxes_o2s_tozero
!-----------------------------------------------------------------------    

!-----------------------------------------------------------------------
      SUBROUTINE set_values_constant()
!-----------------------------------------------------------------------

      USE declars_mod, ONLY: LT, NOC_CBR
      USE loveclim_transfer_mod, ONLY: MGT

      USE mbiota_mod, ONLY: TPP_mafond, caco3_mafond,                  &
                temp_mafond, salt_mafond, ooxy_mafond, clay_mafond,    &
                clay_mafond, tracer01_mafond,                          &
                oco2_mafond, oco3_mafond, ohco3_mafond

      USE mbiota_mod, ONLY: ono3_mafond, opo4_mafond
      USE mod_iloveclim_setup, ONLY: kfs_fond

      IMPLICIT NONE


      INTEGER :: i, j, n, kyr, k

      DO i = 1, LT
        DO n = 1, NOC_CBR

          IF ( MGT(i,kfs_fond(i,n),n) == 1 ) THEN

            TPP_mafond(i,n)    =  12.0D+00 ! [mumolC cm**-2 yr**-1] ???
            caco3_mafond(i,n)  =  10.0D+00 ! [mumolC cm**-2 yr**-1] ???

            temp_mafond(i,n)   =   2.0D+00 ! [°C]
            salt_mafond(i,n)   =  35.0D+00 ! [-]

            ooxy_mafond(i,n)   = 100.0D+03 ! [µmol kg**-1]

            oco2_mafond(i,n)   = 0.030D-03 ! [mol kg**-1]
            ohco3_mafond(i,n)  = 2.230D-03 ! [mol kg**-1]
            oco3_mafond(i,n)   = 0.090D-03 ! [mol kg**-1]

            ono3_mafond(i,n)   = 0.0D-03   ! [mumol kg**-1]
            opo4_mafond(i,n)   = 0.0D-03   ! [mumol kg**-1]
            clay_mafond(i,n)   = 2.0D-03   ! [kg m**-2 yr**-1]

            tracer01_mafond(i,n) = 1.0D+00

#ifdef WITH_C13
            oc13_mafond(i,n)     = 0.0D+00
            TPPC13_mafond(i,n)   = 0.0D+00
            caco3C13_mafond(i,n) = 0.0D+00
#endif

#ifdef WITH_C14
            TPPC14_mafond(i,n)   = 0.0D+00
            caco3C14_mafond(i,n) = 0.0D+00
            oc14_mafond(i,n)     = 0.0D+00
#endif

#ifdef WITH_O18
            caco3O18_mafond(i,n) = 0.0D+00
#endif

          END IF

        END DO
      END DO

      PRINT*, "ocean values for MEDUSA set constant"

!----------------------------------------------------------------------- 
      END SUBROUTINE set_values_constant
!-----------------------------------------------------------------------
#endif

!=======================================================================
       END MODULE flux_to_sediments_mod
!=======================================================================
