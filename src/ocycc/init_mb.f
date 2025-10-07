!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
!********************************************************************
                    SUBROUTINE INIT_MB
!********************************************************************
!              OCEAN MARINE BIO INITIALISATION
!
!  By J. Bendtsen
!  Last modification: 01.02.99
!********************************************************************
      use global_constants_mod, only: dblp=>dp, ip           

      USE declars_mod, ONLY: JT, LT, NOC_CBR
      USE loveclim_transfer_mod, ONLY:

      USE marine_bio_mod, ONLY: OPO4, ONO3, ODOC, ODOCS, ODIC, OSI, OO2, OALK,
     &       OPOC, OC13, OC14
      USE marine_bio_mod, ONLY:OALK_ini, ODOCS_ini, ONO3_ini, OPO4_ini, OSI_ini
      
      USE mbiota_mod, ONLY: PHYTO_M, ZOO_M
      USE carbone_co2, ONLY: c14rstd
      USE iso_dioxygen_mod, ONLY: iair

#if ( OOISO == 1 )
      USE iso_dioxygen_mod, ONLY: iair16, iair17, iair18
#endif

      IMPLICIT NONE

      integer(kind=ip) :: i,j,n,km
      REAL(kind=dblp), PARAMETER ::
     &     PHYTO_ini    = 0.01_dblp,
     &     ZOO_ini      = 0.01_dblp,
     &     OO2_ini      = 250._dblp,    ! mumol/kg
!     &     OO2_iso01    = 249.4_dblp,   ! mumol/kg
!     &     OO2_iso02    = 0.09612_dblp, ! mumol/kg
!     &     OO2_iso03    = 0.51327_dblp, ! mumol/kg
     &     OO2_iso02    = 0.0959_dblp, ! mumol/kg
     &     OO2_iso03    = 0.5121_dblp, ! mumol/kg
     &     ODIC_ini     = 2257.77165E-06_dblp,
     &     ODOC_ini     = 0.0_dblp,
     &     OPOC_ini     = 0.0_dblp,
     &     OC13_cst_ini = -0.1_dblp,
     &     OC14_cst_ini = 0.85_dblp !vm coeff 85 donne par Toggweiler, 1989
!********************************************************************


!dmr --- REFACTORING: streamlining the ini, no use to have several different
!                     ways of getting the values

!...1) definition of initial nutrient, oxygen and carbon distributions

      PHYTO_M(:,:,:)= PHYTO_ini
      ZOO_M(:,:,:)  = ZOO_ini

      OPO4(:,:,:)   = OPO4_ini             ! mumol/kg
      ONO3(:,:,:)   = ONO3_ini             ! mumol/kg
      OSI(:,:,:)    = OSI_ini

      OO2(:,:,:,iair)  = OO2_ini            ! mumol/kg
#if ( OOISO == 1 )
      OO2(:,:,:,iair17)  = OO2_iso02
      OO2(:,:,:,iair18)  = OO2_iso03
      OO2(:,:,:,iair16)  = OO2(:,:,:,iair)-OO2(:,:,:,iair18)-OO2(:,:,:,iair17)
#endif 

      OALK(:,:,:)  = OALK_ini

#if ( OXNITREUX == 1 )
      ON2O(:,:,:)  = ON2O_ini
#endif
      ODIC(:,:,:)  = ODIC_ini
      ODOCS(:,:,:) = ODOCS_ini
      ODOC(:,:,:)  = ODOC_ini
      OPOC(:,:,:)  = OPOC_ini

      OC13(:,:,:) = OC13_cst_ini*ODIC(:,:,:)
      OC14(:,:,:) = OC14_cst_ini*c14rstd*ODIC(:,:,:)

      return
      END SUBROUTINE INIT_MB
