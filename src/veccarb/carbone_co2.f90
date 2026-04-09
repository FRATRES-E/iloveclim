!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module sert à declarer les variables co2 du cycle du
!       carbone
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche, Nathaelle Bouttes
!      Date   : 14 decembre 2009
!      Derniere modification : 09 octobre 2012, Didier M. Roche,
!      Veronique Mariotti
!-----|--1--------2---------3---------4---------5---------6---------7-|

       MODULE carbone_co2

!   INSERER ICI LES EVENTUELS "USE MODULE"

!cnb
#if ( COMATM == 1 )
      use comatm, only: nlat, nlon
#endif

       IMPLICIT NONE

#if ( COMATM == 0 )
#include "veget.h"
#include "comsurf.h"
#endif


       INTEGER, PARAMETER :: new_run_c = 1

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Declaration des CO2 initiaux et courant
!-----|--1--------2---------3---------4---------5---------6---------7-|
       REAL :: PA_C, PA0_C, PA_C_D
       REAL, PARAMETER :: C14DEC = 1.d0/8240.d0, c14rstd = 1.176d-12
       REAL :: C14ATM, C14ATM0, C14ATM_rest
#if ( KC14 == 1 )
       REAL :: PRODC14, MPRODC14, LIMC14MASSE, N14LIM
#endif

! modified the 02/12/2026 EA for the production file in ec14.f

#if ( KC14P == 1)
       INTEGER, PARAMETER :: NC14max=2
       REAL, dimension(NC14max) :: PC14M
       INTEGER :: KTIME
       INTEGER :: n
       INTEGER :: NC14
       INTEGER :: NYR
       INTEGER :: NYR0
       INTEGER  :: NYR01
       INTEGER  :: NYRSR
       REAL :: PC14VAR
       REAL, dimension(NC14max) :: TPSC14

#endif
#if ( O2ATM == 1 )
       REAL :: PA_O, PA0_O
#endif


#if ( CEMIS == 1 )
       integer, parameter :: nb_emis= 450 !nb of lines with data in emission file
       !integer, parameter :: nb_emis= 110 !nb of lines with data in emission file
       real, dimension(nb_emis) :: cemis
#endif

       END MODULE carbone_co2
