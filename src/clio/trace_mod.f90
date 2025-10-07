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
!      Derniere modification : 19 Aout 2014
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module trace_mod

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  bloc "trace.com" : incorpore par instruction 'include' dans les routines
!      liees aux traceurs : cfc_flx.f forcng.f
!  modif : 12/04/99

!--blocs common :
!   Constants needed for solubility computations of CFCs in seawater
!       from Warner and Weiss (1985, Deep Sea Res., equation 6).
!   Note that F is in mol.kg-1.atm-1 following table 5 of WW(1985).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      real(kind=8), parameter ::                    &
                    a1cfc11=-232.0411d0             &
                  , a2cfc11=322.5546d0              &
                  , a3cfc11=120.4956d0              &
                  , a4cfc11=-1.39165d0              &
                  , b1cfc11=-.146531d0              &
                  , b2cfc11=0.093621d0              &
                  , b3cfc11=-.0160693d0

      real(kind=8), parameter ::                    &
                    a1cfc12=-220.2120d0             &
                  , a2cfc12=301.8695d0              &
                  , a3cfc12=114.8533d0              &
                  , a4cfc12=-1.39165d0              &
                  , b1cfc12=-.147718d0              &
                  , b2cfc12=0.093175d0              &
                  , b3cfc12=-.0157340d0

!       cfc11atm = variable holding time-interpolated SH CFC-11
!       cfc12atm = variable holding time-interpolated SH CFC-12


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "cfcmhe1"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(kind=8), dimension(2) :: atmcfc11, atmcfc12
        real(kind=8), dimension(66,2) :: cfc11, cfc12

      end module trace_mod

!--fin du fichier "trace.com"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
