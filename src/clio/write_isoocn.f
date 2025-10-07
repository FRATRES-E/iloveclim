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
!      Derniere modification : 25 septembre 2018 (dmr)
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE write_isoocn

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : scal ocean tracers, of ns_max size
!       Variables de sortie : ~
!-----|--1--------2---------3---------4---------5---------6---------7-|
#if ( ISOOCN >= 1 )

!! START_OF_USE_SECTION

      use para0_mod, only: owisostrt, owisostop
      use bloc0_mod, only: scal

      use global_constants_mod, only: dblp=>dp, ip

!! END_OF_USE_SECTION

      implicit none
      
      integer :: wisoocn_restartdat_id

!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"

!! END_OF_INCLUDE_SECTION

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

         OPEN(newunit=wisoocn_restartdat_id,FILE='wisoocn_restart.dat'
     >       ,STATUS='unknown',FORM='unformatted')

         WRITE(wisoocn_restartdat_id) scal(:,:,:,owisostrt:owisostop)
         CLOSE(wisoocn_restartdat_id)

#endif
       END SUBROUTINE write_isoocn
