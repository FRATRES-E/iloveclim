!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine <squelette> sert a ???
!       et caetera, et caetera
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche, Thibaut Caley 
!      Date   : 02 mai 2012    
!      Derniere modification : ~
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE write_isolbm

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : ~
!       Variables de sortie : ~
!-----|--1--------2---------3---------4---------5---------6---------7-|
#if ( WISOATM == 1 && ISOLBM >= 1 )

       use comatm, only: iwater, nwisos

       use comland_mod, only: bmoisg, runofo, dsnow

       IMPLICIT NONE

       integer :: wisolbm_restartdat_id

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

!-----|--1--------2---------3---------4---------5---------6---------7-|
!-----|--1--------2---------3---------4---------5---------6---------7-|


      OPEN(newunit=wisolbm_restartdat_id,FILE='wisolbm_restart.dat'
     >    ,STATUS='unknown',FORM='unformatted')

        WRITE(wisolbm_restartdat_id) bmoisg(:,:,iwater+1:nwisos),
     >              runofo(:,:,iwater+1:nwisos),
     >              dsnow(:,:,iwater+1:nwisos)

         CLOSE(wisolbm_restartdat_id)
#endif
       END SUBROUTINE write_isolbm
