!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce sous-programme permet l'initialisation du routage de l'eau à
!       partir de la topographie et du masque terre/ocean dans ECBilt
!      ... dans l'environnement logiciel LUDUS
!
!     NOTA BENE
!     =========
!     Il est rendu obligatoire à cause de l'incompatibilité foncière
!      entre les inclusifs comland.h comdyn.h qui ne me permettent pas
!      un accès à la fois à fractl et à rmount  !!!
!
!
!      Auteur : Didier M. Roche
!      Date   : 04 décembre 2010
!      Derniere modification : idem
!-----|--1--------2---------3---------4---------5---------6---------7-|

      SUBROUTINE ec_inirunoff_2

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree :
!       Variables de sorties :
!-----|--1--------2---------3---------4---------5---------6---------7-|

!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|

      USE routageEAU_mod, ONLY: nbexuts, mask_lnd, exuti, exutj, eni,
     & enj, river_basin
#if ( COMATM == 1 )
      USE comatm
      USE comdyn, only: rmount
#endif


        IMPLICIT NONE

!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|
#if ( COMATM == 0 )
#include "comatm.h"
! Pour avoir la topographie rmount
#include "comdyn.h"
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!        Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|

      CALL routageEAU(nlat,nlon,rmount,mask_lnd,nbexuts,exuti,exutj,
     &                 eni,enj)

      CALL out_routageEAU(nlat,nlon,eni,enj,river_basin)

      END SUBROUTINE ec_inirunoff_2
