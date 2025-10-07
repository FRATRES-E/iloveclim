!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module est un portage en FORTRAN 90 du module CLIO initial
!       cree en FORTRAN 77.
!      (dans l'environnement logiciel LUDUS)
!
!      Auteur : ??, Didier M. Roche
!      Date   : ??
!      Derniere modification : 05 Avril 2017
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module const_mod

      use global_constants_mod, only: dp

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  fichier "const.com", incorpore par instruction 'include' dans les programmes
!   (et les routines des programmes) :
!      OM, CLASS, FEEDOM, GRADP, STATIS, DIFBIN, TROMNC, UNIBIN, TRSBATH.
!  modif : 06/02/98
!  modif : 01/03/04 <- declarations explicites pour inclusion dans LOCH

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "cstfix"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      real(kind=dp) ::  cpo, cstmin, cstmax, epsil, gpes, omega, one, pi, radian, rho0, rterre, svrdrp, unsrt, untour, yeaday, zero


      end module const_mod
!--fin du fichier "const.com"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
