!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:50 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:50 CET 2009
! module de declaration de variables
! marine_bio
! remplace marine_bio.inc
! C. Dumas 29-01-2007
!  Verifie pour C2.5 07-12-2007 C. Dumas & D. Roche

module marine_bio_mod

USE declars_mod
use para0_mod, ONLY: NISOO2

implicit none
! 15 chiffres significatifs et exposant a 3 chiffres :
INTEGER, parameter :: pr=selected_real_kind(15,3)


REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: OALK, ODIC
REAL(kind=8), dimension(LT,JT,NOC_CBR,NISOO2), target :: OO2
REAL, dimension(LT,JT,NOC_CBR) :: CM, C14M, OPO4, ONO3, OSI,  ODIC_diff = 0.0d0
REAL, dimension(LT,JT,NOC_CBR) :: ODOC, OC13, OC14, ODOCS, ODOC13,ODOCS13
REAL, dimension(LT,JT,NOC_CBR) :: OPOC ! OPOC is always zero, it should be removed but kept for restart
!### REAL, dimension(LT,JT,NOC_CBR) :: OPOC_flux,caco3_flux !nb flux entrant
#if ( OXNITREUX == 1 )
REAL, dimension(LT,JT,NOC_CBR) :: ON2O !vm n2o
#endif
!vm&dmr --- Ajout de la production nette en variable globale
DOUBLE PRECISION, dimension(LT,JT,NOC_CBR) :: phyto_prod

!cnb - sediments
REAL*8, dimension(LT,JT,NOC_CBR) :: OCO3 !nb useless, to be deleted: OCO3SAT(LT,JT,NOC_CBR)
!nb useless, to deleted? REAL*8, dimension(LT,JT,NOC_CBR) :: OCAL(LT,JT,NOC_CBR)

REAL :: pco2_diag

!* SURFACE
!**********
#if ( KC14 == 1 )
REAL :: FDIC !vm  FDIC = flux carbon ocn->atm 
REAL :: FOAC14 !vm  FOAC14 = flux C14 ocn->atm
#endif
REAL, dimension(LT,NOC_CBR) :: FOPO4, FONO3, FOSI, FOALK, FODIC, FODOC
REAL, dimension(LT,NOC_CBR, NISOO2) :: FOO2
REAL, dimension(LT,NOC_CBR) :: FOC13, FOCO2, FOC14, FODOCS, FODOC13 
REAL, dimension(LT,NOC_CBR) :: FODOCS13, oxCO2 = 0.0d0, oxpCO2
REAL, dimension(LT,NOC_CBR) :: osCO2, oxHCO3, oxCO3
#if ( OXNITREUX == 1 )
REAL, dimension(LT,NOC_CBR) :: FON2O !vm: flux n2o
#endif

!* PRODUCTION
REAL(pr) :: eher, zinges, ecan, sigma_md, sigma_m

!REFACTORING DONE: replaced Oeta(:,1:5) by OetaC*, OetaN* and OetaO2*
!  POM related arrays are 3D, DOM related ones remain 1D
!  as they do not seem to have to fullfil any global constraints.
!  Also removed commented parts for clarity.
REAL(pr), dimension(LT, JT, NOC_CBR) :: OetaC_POMoxid, OetaN_POMoxid, OetaO2_POMoxid
REAL(pr), dimension(LT, JX, NOC_CBR) :: OetaC_POMrain, OetaN_POMrain
REAL(pr), dimension(JT) :: OetaC_DOMoxid_1D, OetaN_DOMoxid_1D, OetaO2_DOMoxid_1D
REAL(pr), dimension(NOC_CBR, LT) :: OetaC_POMsedin, OetaN_POMsedin
REAL(pr), dimension(NOC_CBR, LT) :: OetaO2_POMsedin


!tbd REAl(pr), dimension(JT,5) :: Oeta_floor
REAL :: oc13frac, oc13fair, oc13bio
!cnb - sediment model
REAL(pr) :: ter_car_flux, ter_car_total


!* Initial conditions for conservative tracers
REAL :: OPO4_ini, ONO3_ini, OSI_ini, OALK_ini, ODOCS_ini
#if ( OXNITREUX == 1 )
REAL :: ON2O_ini !vm n2o
#endif

!* dmr --- Ajout pour coherence entre mbiota et maphot : profondeur de production & du niveau 500
!* dmr !!! Attention ces declarions sont modeles dependantes
!* dmr ---
!* dmr Dans CLIMBER : J500 = 7 et JPROD = 2
!* dmr --- Pour CLIO : 
!cnb INTEGER, PARAMETER :: J500 = 10, JPROD = 6
INTEGER, PARAMETER :: J500 = 12, JPROD = 6

REAL(pr) :: global_dic_1
REAL(pr) :: global_po4_1
REAL(pr) :: global_no3_1
REAL(pr) :: global_alk_1
REAL(pr) :: global_doc_1
REAL(pr) :: global_docs_1
REAL(pr) :: global_o2_1
REAL(pr) :: global_oc13_1
REAL(pr) :: global_oc14_1

#if ( BATHY >= 1 )
!nb global value for correction
REAL(pr) :: global_dic
REAL(pr) :: global_po4
REAL(pr) :: global_no3
REAL(pr) :: global_alk
REAL(pr) :: global_oc13
REAL(pr) :: global_oc14
REAL(pr) :: tot_alk_prev
REAL(pr) :: tot_phos_prev
REAL(pr) :: tot_nit_prev
#endif

end module marine_bio_mod
