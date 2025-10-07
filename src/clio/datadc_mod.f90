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

      module datadc_mod

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! fichier "datadc.com" (anciennement dans fichier "varno.com")
!  regroupant les instructions "DATA" ; inclus apres TOUTES les declarations !
! incorpore par instruction 'include' dans les routines :
!  CLASS, STATIS, FLTROM, defcst, savrunb, redrunb, savrunc, redrunc, foroutp,
!     ncdfout, rdclpar, process, meridflu, local, binout
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      integer ::  nvrt, nvrs, nvrs3, nvrs4, nvrfw, nvrfc, nvrfs, nvrfs3, nvrfs4, nvreac, nvret, nvrub, nvrvb, nvrps, nvrhac, nvrajc&
             , nvrb, nvru, nvrv, nvrtke, nvrw, nvrn2, nvras, nvrau, nvrusl, nvrvsl

!- variables pour Sea Ice :
      integer :: nvrhg, nvrfq, nvrqs, nvral, nvrhn, nvrts

!- variables pour Bilan detaille :
      integer :: nvraxt, nvrayt, nvrhat, nvrhdt, nvrvat, nvrvdt, nvrvaf, nvrvdf

!- variables pour Sea Ice :
      integer :: nvrum, nvrvm, nvrug, nvrvg, nvrtbq, nvrxzo, nvrtgx, nvrtgy, nvrmom

!- variables pour diffus. Isopyc. ou G.&M.W. :
      integer :: nvravi, nvravs, nvrslx, nvrsly, nvrpsx, nvrpsy, nvrugm, nvrvgm, nvrwgm

!--definition du numero de chaque variable (pour + de transparence) :

      data  nvrt   /  1 /
      data  nvrs   /  2 /
      data  nvrs3  /  3 /
      data  nvrs4  /  4 /
      data  nvrfw  / 10 /
      data  nvrfc  / 11 /
      data  nvrfs  / 12 /
      data  nvrfs3 / 13 /
      data  nvrfs4 / 14 /

      data  nvreac / 20 /
      data  nvret  / 21 /
      data  nvrub  / 22 /
      data  nvrvb  / 23 /
      data  nvrps  / 24 /
      data  nvrhac / 25 /

      data  nvrajc / 30 /
      data  nvrb   / 31 /
      data  nvru   / 32 /
      data  nvrv   / 33 /
      data  nvrtke / 34 /
      data  nvrw   / 35 /
      data  nvrn2  / 36 /
      data  nvras  / 37 /
      data  nvrau  / 38 /
      data  nvrusl / 40 /
      data  nvrvsl / 41 /

!- variables pour Sea Ice :
      data  nvrhg  / 50 /
      data  nvrfq  / 51 /
      data  nvrqs  / 52 /
      data  nvral  / 53 /
      data  nvrhn  / 54 /
      data  nvrts  / 55 /

!- variables pour Bilan detaille :
      data  nvraxt / 60 /
      data  nvrayt / 61 /
      data  nvrhat / 62 /
      data  nvrhdt / 63 /
      data  nvrvat / 64 /
      data  nvrvdt / 65 /
      data  nvrvaf / 66 /
      data  nvrvdf / 67 /

!- variables pour Sea Ice :
      data  nvrum  / 70 /
      data  nvrvm  / 71 /
      data  nvrug  / 72 /
      data  nvrvg  / 73 /
      data  nvrtbq / 74 /
      data  nvrxzo / 75 /
      data  nvrtgx / 76 /
      data  nvrtgy / 77 /
      data  nvrmom / 80 /

!- variables pour diffus. Isopyc. ou G.&M.W. :
      data  nvravi / 90 /
      data  nvravs / 91 /
      data  nvrslx / 92 /
      data  nvrsly / 93 /
      data  nvrpsx / 94 /
      data  nvrpsy / 95 /
      data  nvrugm / 96 /
      data  nvrvgm / 97 /
      data  nvrwgm / 98 /

      end module datadc_mod
!--fin du fichier "datadc.com"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
