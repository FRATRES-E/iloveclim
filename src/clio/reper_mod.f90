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

      module reper_mod

      use para0_mod, only: imax, jmax, kmax, nsmax
      use para_mod,  only: nchsep, nbsmax, nhsfmx, ninfmx


      implicit none

      private :: nchsep


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  fichier "reper.com" : incorpore par instruction 'include' dans les routines :
!     CLASS, GRADP, FEEDOM, UNIBIN, barot,
!     informe, defcst, defgrid, geogra, redforc,
!     ncdfout, moyen, streamh, streamv, meridflu,
!     checkwst, local, binout, sepgl, defgl, unigl.
! (commun a toutes les routines de traitement des resultats et de
!  definition du domaine) ; inclus apres "type.com", "para.com".
!-----
!  modif : 29/03/99

!--blocs common :

!-----
!- variables that must be transfered from defcst to defgrid and redforc
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "kdefini"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer, dimension (imax)      :: kbath1, kbath2
      integer, dimension (imax,jmax) :: kbath

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "kdefini"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8)                      :: unstyr
      real(kind=8), dimension(kmax)     :: rapp1
      real(kind=8), dimension(0:nsmax)  :: yforc, unitfx
      real(kind=8), dimension(0:nsmax+2) :: rapp0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "cfilcor"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      character(len=40) :: filcor

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "jinform"
!-----
!- variables utilisees dans routine "informe" pour sortie sur fich. "evolu" :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer                    :: nvinfo,nvhsf,nferme,icheck,jcheck, kcheck, nocean, ksud, jmsud, knor, jmnor
      integer, dimension (ninfmx):: ktsum


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "sinform"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8)                    :: zmdeta
      real(kind=8), dimension(nsmax)  :: scalwr
      real(kind=8), dimension(ninfmx) :: vinfor, vinfom


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "cinform"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      character(len=30)                        :: fmtw
      character(len=nchsep), dimension(ninfmx) :: titvar

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "distances"
!-----
!- variables utilisees dans la definition du domaine :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8) :: dlong, dlat, dniv, xlon1, ylat1, zniv1, xalon1, yalat1, dxaj, dyai, xaj1, yai1, xwpoln, dxwi, xwi1, dywj, ywj1

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "morceaux"
!- indices delimitant les grilles, zones, bassins et detroits :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer                                  :: ndhsf, jmvlat
      integer, dimension (imax)                :: jsep, jnorth
      integer, dimension (nhsfmx)              :: ishsf,iehsf,jshsf,jehsf
      integer, dimension (0:nbsmax)            :: jsbas, jebas, jezon
      integer, dimension (imax,jmax)           :: jmaxl
      integer, dimension (imax,jmax,0:1)       :: jgeogr
      integer, dimension (jmax,0:nbsmax)       :: iszon,iezon,isbas
      integer, dimension (-imax:imax+imax)     :: icl
      integer, dimension (jmax,0:kmax,0:nbsmax):: iebas


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "vraielat"
!- variables et coeff. utilisee pour Moyenne_Zonale(Vraie_Latitude)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax)    :: yvrlat, ageogr, bgeogr
      real(kind=8), dimension(imax,jmax,0:1):: rgeogr

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "cclieu"
!- titres associes :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      character(len=8 )                      :: ttvlat
      character(len=20), dimension(nhsfmx)   :: tithsf
      character(len=10), dimension(0:nbsmax) :: titzon

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----
!--variables utilisees pour regrouper les resultats sur a grille Globale :
!     common / globgrid /
!    & tmgl(imax,jmax,kmax,2), csgl(imax,jmax,4), surfgl(imax,jmax,2)

!     common / limglob /
!    & kfgl(imax,jmax), imdl1(jmax), imdl2
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      end module reper_mod
!--fin du fichier "reper.com"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
