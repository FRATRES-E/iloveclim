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

       SUBROUTINE write_isoatm

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : ~
!       Variables de sortie : ~
!-----|--1--------2---------3---------4---------5---------6---------7-|
#if ( ISOATM >= 1 )

       USE iso_param_mod, ONLY : ieau, neauiso


#if ( COMATM == 1 )
      USE comatm
      USE comphys
#endif

      implicit none

#if ( COMATM == 0 )
#include "comatm.h"
#include "comphys.h"
#endif

      integer :: wisoatm_restartdat_id

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

!-----|--1--------2---------3---------4---------5---------6---------7-|
!-----|--1--------2---------3---------4---------5---------6---------7-|


      OPEN(newunit=wisoatm_restartdat_id,FILE='wisoatm_restart.dat'
     >    ,STATUS='unknown',FORM='unformatted')

       WRITE(wisoatm_restartdat_id) rmoisg(:,:,ieau+1:neauiso)
     >              ,torain(:,:,ieau+1:neauiso)
     >              ,tosnow(:,:,ieau+1:neauiso)

         CLOSE(wisoatm_restartdat_id)
#endif
       END SUBROUTINE write_isoatm
