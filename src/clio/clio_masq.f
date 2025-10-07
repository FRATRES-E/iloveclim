!dmr -- Ajout du choix optionnel des composantes - Tue Dec 15 16:33:22 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Tue Dec 15 16:33:22 CET 2009

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       This subroutine computes an integer mask on the ocean grid for
!        the iceberg coupling
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche , Marianne Buegelmayer
!      Date   : 25 Aout 2011
!      Derniere modification :
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE clio_masq

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  :
!       Variables de sortie :
!-----|--1--------2---------3---------4---------5---------6---------7-|

!! START_OF_USE_SECTION

       USE varsCliotemp_mod, ONLY: ocn_mask

!       IMPLICIT NONE

       use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
!! END_OF_USE_SECTION

!! START_OF_INCLUDE_SECTION


#include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"

!! END_OF_INCLUDE_SECTION


!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       ocn_mask = CEILING(tms(:,:,kmax))

       END SUBROUTINE clio_masq
