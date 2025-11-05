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

      module para_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  fichier "para.com" : incorpore par instruction 'include' dans les programmes
!   (et les routines des programmes) :
!      OM, CLASS, FEEDOM, GRADP, STATIS, DIFBIN, TROMNC, UNIBIN, TRSBATH.
! contient la plupart des parametres definissant les tableaux courrants :
!--ocean mondial, 3x3 deg, version 2 Sous-grilles reunies en une seule.
!-----
! Ctke [Ctk0] => ligne specifique a la version avec [sans] TKE .
! Cice [Cic0] => ligne specifique a la version avec [sans] glace marine .
!-----
!  modif : 01/05/98
!
! Reorganisation : 1/03/2004; A. Mouchet pour couplage avec Loch
! -> cree para0.com inclus dans para.com et dans declarations Loch
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-parametres lies a la taille du domaine; freq. des donnees; etc.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      USE para0_mod, ONLY: kmax, nsmax

      implicit none

      private :: kmax, nsmax

!--parametres donnant le nombre de niveaux dans la glace
      integer, parameter :: nkb0 = 3

!--parametres fixant le nombre de zones (= nb bassins), nb de detroits :
! [indispensable pour inclusion du fichier  reper.com]
      integer, parameter :: nbsmax = 3 , nhsfmx = 10, nltmax = 4

!--parametres pour les sorties sur fichier "evolu" par la routine "informe" :
! [indispensable pour inclusion du fichier  reper.com]
!dmr -- faux a priori ...      parameter ( ninfmx = 30 + nhsfmx + 2*kmax + 20)
      integer, parameter :: ninfmx = 36 + nhsfmx + nsmax*kmax + 20, nchinf = 5 , nchsep = nchinf+2

!--parametres pour couplage : Nombre de tableaux envoyes et recus par l'ocean :
      integer, parameter :: ntocn = 4 , ntatm = 10

!--parametre indiquant le rang "k" (ds le tableau general),
!   du 1er niveau occupe par la variable:
! [indispensable pour inclusion du fichier  vareq.com]

      integer, parameter :: krlu  =  0
      integer, parameter :: krlfw = krlu - nsmax - 5
      integer, parameter :: krlfc = krlfw + 1
      integer, parameter :: krlfs = krlfw + 2
      integer, parameter :: krlfs3= krlfw + 3
      integer, parameter :: krlfs4= krlfw + 4
      integer, parameter :: krlps = krlu - 4
      integer, parameter :: krlet = krlu - 3
      integer, parameter :: krlub = krlu - 2
      integer, parameter :: krlvb = krlu - 1
      integer, parameter :: krlv  = krlu + kmax
      integer, parameter :: krlt  = krlu + kmax*2
      integer, parameter :: krls  = krlu + kmax*3
      integer, parameter :: krls3 = krlu + kmax*4
      integer, parameter :: krls4 = krlu + kmax*5
      integer, parameter :: krlb  = krlu + kmax*(2+nsmax)
      integer, parameter :: krln2 = krlu + kmax*(3+nsmax)
      integer, parameter :: krlas = krlu + kmax*(4+nsmax)
      integer, parameter :: krlau = krlu + kmax*(5+nsmax)
      integer, parameter :: krlw  = krlu + kmax*(6+nsmax)
      integer, parameter :: krltke= krlu + kmax*(7+nsmax) + 1
      integer, parameter :: krlajc= krlu + kmax*(8+nsmax) + 2

!  Avu=>Avi S=>Psx T=>Psy U=>Ugm V=>Vgm W=>Wgm B=>Slx N2=>Sly
      integer, parameter :: krlavi=krlau
      integer, parameter :: krlavs=krlas
      integer, parameter :: krlslx=krlb
      integer, parameter :: krlsly=krln2
      integer, parameter :: krlpsx=krlt
      integer, parameter :: krlpsy=krls
      integer, parameter :: krlugm=krlu
      integer, parameter :: krlvgm=krlv
      integer, parameter :: krlwgm=krlw

      integer, parameter :: krlusl= krlajc + kmax
      integer, parameter :: krlvsl= krlusl + nsmax   + 2
      integer, parameter :: krlhac= krlusl + nsmax*2 + 4
      integer, parameter :: krleac= krlhac + 1

!ic0  parameter( krlhg = krlajc )
      integer, parameter :: krlhg = krleac + 1
      integer, parameter :: krlfq = krlhg  + 1
      integer, parameter :: krlqs = krlhg  + 2
      integer, parameter :: krlal = krlhg  + 3
      integer, parameter :: krlhn = krlhg  + 4
      integer, parameter :: krlts = krlhg  + 5
      integer, parameter :: krlug = krlhg  + 6
      integer, parameter :: krlvg = krlhg  + 7
!ic0  parameter( kvsice = 0 )
      integer, parameter :: kvsice = 8

      integer, parameter :: krlvaf= krlt
      integer, parameter :: krlvdf= krlvaf + kmax
      integer, parameter :: krlaxt= krlvaf + kmax*2
      integer, parameter :: krlayt= krlvaf + kmax*3
      integer, parameter :: krlhat= krlvaf + kmax*4
      integer, parameter :: krlhdt= krlvaf + kmax*5
      integer, parameter :: krlvat= krlvaf + kmax*6
      integer, parameter :: krlvdt= krlvaf + kmax*7 + 1

!--parametres lies a la definition du tableau general utilise pour les sorties:
! [indispensable pour inclusion du fichier  var??.com]
      integer, parameter :: ltymax = 11
      integer, parameter :: krlmin = -5 - nsmax
      integer, parameter :: nvmax = 99 , nv3dmx = 9+nsmax
!- nv2dmx = Nb. Var. 2D (fixe+Ns+Ice) ; kv2dmx = Nb. Niv. reserves pour var. 2D
      integer, parameter :: nv2dmx = 9+nsmax+15 , kv2dmx = 11 + 3*nsmax &
             + kvsice
      integer, parameter :: krlmax = krlmin + 1 + nv3dmx*kmax + kv2dmx

      end module para_mod
!--fin du fichier "para.com"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
