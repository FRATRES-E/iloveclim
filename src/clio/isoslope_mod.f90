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

      module isoslope_mod

      use para0_mod, only: imax, jmax, kmax

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! originally written with IMPLICIT INTEGER (I-N), REAL (A-H, O-Z)
! transferred as such in fortran90
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  fichier "isoslope.com" : incorpore par instruction 'include' dans :
!    isoadvec isodiffu isoslope scali [iso_copy iso_flux isofilter]
! (commun a toutes les routines de diffus. Isopycnale et Genk & McWillliams)
! regroupe isoslope.com + isoadvec.com de P.P.M. (vers. du 23/03/97 & 16/10/96)
!-----
!  modif : 25/05/99

!  ppmodif : 23-03-97 : iso,gm90  (<- fichier isoslope.com d'origine)
!------------------------------------------------------------------
! GLOBAL VARIABLES for ISOPYCNAL DIFFUSION
!------------------------------------------------------------------
!              B-GRID faces (V=velocity) (S=scalar)
!------------------------------------------------------------------
!        WEST      SOUTH   MIDDLE   WEST     SOUTH     MIDDLE
!        + + +     + + +   + + +    + + +    + + +     V + V
!        + 1 3     3 2 +   1 0 +    V S V    V S V     + S +
!        + 5 +     + 6 +   + 2 +    + + +    + + +     V + V
!------------------------------------------------------------------
!  ai(k) = isopycnal mixing coefficient (m^2/s) (added to ahh)
!  slopemax(k) = maximum slope of isopycnals for ISO  (iroutine=1)
!  slopmgm(k)  = maximum slope of isopycnals for GM90 (iroutine=2)
!  ttm1=mask pt 1, ttm2=mask pt 2
!  c1x=drox/droz,c2y=droy/droz  where (droz=bvf)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "isopycnale"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax,kmax)  :: ttm1, ttm2, c1x, c2y, c4x, c4y ! dmr [DELETED], dscalz
!      real(kind=8), dimension(imax,jmax,kmax+1):: fisoz ! dmr [DELETED], fisoz

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  ppmodif : 16-10-96 : gm90  (<- fichier isoadvec.com ajoute ici)
!  aitd(k) = isopycnal thickness coefficient (m^2/s)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "isobolus"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax,kmax)  :: c4xgm, c4ygm

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "isobolus2"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax,kmax)  :: uiso, viso
      real(kind=8), dimension(imax,jmax,kmax+1):: wiso

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "indexiso"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer, dimension(imax)            :: jsdom1, jsdom2

      end module isoslope_mod
!--fin du fichier "isoslope.com"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
