!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:55 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:55 CET 2009
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module contient les variables du routage interactif de l'eau
!       (dans l'environnement logiciel LUDUS)
!   
!      Auteur : Didier M. Roche 
!      Date   : 04 décembre 2010
!      Derniere modification : idem
!-----|--1--------2---------3---------4---------5---------6---------7-|

       MODULE routageEAU_mod

       use taillesGrilles,  only: ilat=>iEcb, jlon=>jEcb

#if ( F_PALAEO == 3 )
       use taillesGrilles,  only: grislat=>sgnx, grislon=>sgny 
#endif

       IMPLICIT NONE

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Pas beau mais je n'ai pas le choix tant qu'ECBilt n'est pas en
!        modules (vrai pour lbm aussi)
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       INTEGER, PARAMETER, PRIVATE :: ilat = 32, jlon=64
       INTEGER, PARAMETER, PRIVATE :: nb_max_rbas = 13
!        INTEGER, PARAMETER, PRIVATE :: grislat = 241, grislon=241

       
       
       INTEGER, PARAMETER :: nbexuts=100000

       INTEGER, DIMENSION(nbexuts) :: exuti, exutj
       INTEGER, DIMENSION(ilat,jlon) :: eni, enj
       INTEGER, DIMENSION(ilat,jlon) :: mask_lnd
       INTEGER, DIMENSION(ilat,jlon) :: river_basin
#if ( F_PALAEO == 3 )
       INTEGER, DIMENSION(nbexuts) :: exutiSG, exutjSG
       INTEGER, DIMENSION(grislat,grislon) :: eniSG, enjSG
#endif

       REAL, PARAMETER :: epsi_lon = 1.0D-10

       REAL, DIMENSION(nb_max_rbas) :: flux_per_basin = 0.0d0
!!!### #if ( WATGRISCONS == 1 )
!!!###        REAL, DIMENSION(ilat,jlon) :: negat_precip
!!!### #endif

       END MODULE routageEAU_mod
