!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:58 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:58 CET 2009
! module de declaration de variables
! C_RESERVOIRS
! remplace C_RESERVOIRS.inc
! N. Bouttes 04-12-2007
!  Verifie pour C2.5 04-03-2008 N. Bouttes & D. Roche

!* Contient les reservoirs de carbone initiaux et dynamiques
!* et le c13 initial
!* pour les oceans et land

module C_res_mod

! USE declar_mod


use file_libs, only: fileDescriptor

implicit none
!dmr --- Pour essai de compilation ...
INTEGER, PARAMETER :: NCO2max=64924 !230000
!dmr --- Pour essai de compilation ...

!dmr&nb Numero pour le fichier C_reservoir.txt
INTEGER :: c_res_fich
#if ( KC14 == 1 )
INTEGER :: c14_res_fich, c_flux_fich, dc14_fich, c14_bilan_fich
#endif

#if ( CORAL ==1 )
INTEGER :: coral_res_fich, c_nino_fich
#endif

#if ( COASTAL ==1 )
INTEGER :: coastal_res_fich
#endif

type(fileDescriptor), public,save :: c13_res_fich, c_ocean_fich, fPOC_fich, fCAL_fich

REAL*8 :: ca_oc_ini, ca_la_ini
REAL*8 :: ca_oc_rest, ca_la_rest
REAL*8 :: alk_oc_rest
REAL*8 :: cav_oc=0.0d0, cav_la=0.0d0, cav_oc2
REAL*8 :: ca13_oc_ini, ca13_la_ini
REAL*8 :: ca13_oc_rest, ca13_la_rest
REAL*8 :: ca13_at_ini, dc13at_ini
REAL*8 :: ca_oc_vol,coc_odoc,coc_odocs, coc_odoc13, coc_odocs13 !*8 ?
REAL*8 :: cav_oc13,cav_la13 !*8 ?
REAL*8 :: EMIS_CUM
REAL*8 :: EMIS_C13_CUM
REAL*8 :: CAV_OC_P = 0.0d0
REAL*8 :: CAV_LA_P = 0.0d0
REAL*8 :: c13atm_rest
!cnb REAL*8 :: c13atm, c13init, c14init
REAL*8 :: c13atm, c13init

!dmr&vm --- Pour le14C ...
#if ( KC14 == 1 )

REAL(KIND=8) :: ca14_oc_ini, ca14_la_ini, cav_oc14_b, cav_oc14=0.0d0, cav_oc_b &
, cav_la14_b, cav_la14=0.0d0, cav_la_b, ca14_oc_rest, ca14_la_rest
!dmr&vm --- Fluxes
REAL(KIND=8) :: FC12OA,FC12LA,FC14OA, FC14LA

#endif
!* pour le co2 cycle

INTEGER :: NCO2
REAL*8, DIMENSION(NCO2max) ::  TPSCO2, PCO2M
end module C_res_mod
