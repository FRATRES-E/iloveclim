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

      module varno_mod

      use para_mod, only: nvmax

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! originally written with IMPLICIT INTEGER (I-N), REAL (A-H, O-Z)
! transferred as such in fortran90
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! fichier "varno.com" : incorpore par instruction 'include' dans les routines :
!  CLASS, GRADP, STATIS, TROMNC, FLTROM,
!    defcst, savrunb, redrunb, savrunc, redrunc, foroutp,
!    ncdfout, rdclpar, process, scadhv, meridflu, moyen, local, binout
!  (commun a toutes les routines de sortie (output) de resultats)
!   inclus apres "type.com", "para.com", "bloc.com".
!------------------
! Chaque variable est reperee par un numero "nv" (nv > 0),
!  nv = -1,0 pour fichiers (et autre) regroupant les moyennes longitudinales
!  nv = -2 pour le fichier regroupant les bilans globaux et par niveaux.
!  nv de -10-nsmax a -10 pour transports meridiens de masse (-10) et de scalaire
!-----
! La localisation de la variable sur la maille est consignee
!  dans le tableau "ltyp", croissant dans le sens trigo., =99 => non defini
!     0 a 3 : verticalement centre ; 4 a 7 interface entre 2 niveaux ;
!     8 a 11 : sans localisation verticale ( eta , ub , vb ...)
!     + 12 si 1ere composantE vectorrielle, + 24 si 2nd composante.
!  titres et formats : lus dans le fichier "class.par"
!------------------
!  modif : 30/01/98
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "indeks"
!--common definissant les correspondances (variables, niveaux associes) :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer, dimension(0:nvmax) :: krl1, krlm, kvar2d, ltyp, nvrl

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "cvname"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      character(len=3), dimension(-2:nvmax) :: titcv

      end module varno_mod
!--fin du fichier "varno.com"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
