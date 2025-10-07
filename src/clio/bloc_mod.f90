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

      module bloc_mod

      use para0_mod, only: imax, jmax, kmax, nsmax, nlpmax, ncomax, nrpmax
      use para_mod,  only: nbsmax, nltmax, ntatm, ntocn
      use global_constants_mod, only: days_year360d_i      

      implicit none


      public
! originally written with IMPLICIT INTEGER (I-N), REAL (A-H, O-Z)
! transferred as such in fortran90

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  fichier "bloc.com" : incorpore par instruction 'include' dans les routines :
!      OM, CLASS, FEEDOM, GRADP, STATIS, DIFBIN, TROMNC, UNIBIN,
!      start, flucor, scale, slopez, slopes, scadew, scadns, scali,
!      uve, uvi, barot, uvb0et, uvbfet, uvm,
!      informe, defcst, defgrid, redforc, correct,
!      conti3d, etat, vdiffu, alph2dc, alphdkc, alphdec, raccord,
!      savrunb, redrunb, savrunc, redrunc,
!      ncdfout, moyen, streamv, meridflu, scadhv, checkwst, streamh, stream1h,
!      vague, local, binout, defgl, unigl, foroutp, lisstab, flowsurf.
!  inclus apres "type.com", "para.com".
!  modif : 06/10/99
!  modif : 02/02/04 <- LOCH : 'bloc0.com'
!  modif : 01/03/04 <- LOCH
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "iteration"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!--Tableaux de variables en evolution (+ transfert entre routines) :
      integer :: numit, numspl, ninstb, nclin, nclmoy

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "dynam2"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!--Tableaux de flux et facteurs d'evolution et variables en evolution :
      real(kind=8), dimension(nsmax)                      :: deriv
      real(kind=8), dimension(nlpmax)                     :: alpslp, uvcslp
      real(kind=8), dimension(imax,jmax)                  :: etaspl, ubspl, vbspl, umoy, vmoy
      real(kind=8), dimension(imax,jmax,6)       , target :: phihhh
      real(kind=8), dimension(imax,jmax,kmax)             :: q = 0.0d0, fub,fvb, vlturb, avqdz, tm2tur, avuloc = 0.0d0
      real(kind=8), dimension(imax,jmax,kmax+1,2), target :: phizzz

!sai  common / saison /
!sai &  tmens(imax,jmax,nmois),smens(imax,jmax,nseas),
!sai &  txmens(imax,jmax,nmois),tymens(imax,jmax,nmois),
!sai &  d2tmns(imax,jmax,nmois),d2smns(imax,jmax,nseas),
!sai &  d2txms(imax,jmax,nmois),d2tyms(imax,jmax,nmois),
!sai &  flxss(imax,jmax,nsmax,nmois),d2fxss(imax,jmax,nsmax,nmois),
!sai &  flxsur(imax,jmax,nsmax)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "icouplage"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8) :: master
      integer      :: icoupl, icoutp, itau_slow


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "tbforcing"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax,ntocn):: tfmocn
      real(kind=8), dimension(imax,jmax,ntatm):: ttoocn

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "cfmetric"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax)      :: cmxdy, cmydx, fs2cor, aire, covrai, xang1, xang2
      real(kind=8), dimension(imax,jmax,0:3)  :: smx, smy, cmxy, smxy
      real(kind=8), dimension(imax,jmax,0:3), target :: cmx, cmy

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "surfvol"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8)                               :: zvols, zvolo, zvolv, zvolw, zsurf, unsvol
      real(kind=8), dimension(kmax)              :: zsurfs, zsurfo, zsurfv
      real(kind=8), dimension(0:nltmax)          :: zvolsla, zvolola
      real(kind=8), dimension(0:nbsmax)          :: zvolsba, zvoloba
      real(kind=8), dimension(kmax,0:nltmax)     :: zsurfsla, zsurfola
      real(kind=8), dimension(kmax,0:nbsmax)     :: zsurfsba, zsurfoba
      real(kind=8), dimension(imax,jmax,kmax,0:1):: ctmi

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "localise"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer                              :: nslp(nlpmax), nslpdw, nxslp, nxyslp
      integer, dimension(kmax)             :: n1coin, n2coin, n3coin, n4coin
      integer, dimension(nlpmax)           :: ijslp, kslp, lslp
      integer, dimension(kmax,nsmax)       :: nrap
      integer, dimension(ncomax,kmax)      :: i1coin, i2coin, i3coin, i4coin
      integer, dimension(nrpmax,kmax,nsmax):: ijrap

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "lerun"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer :: nstart, nlast, nsplaj, nsplit, ninfo, nwjl, nwtal, nwm, nwa , nwtest, idyn, nwtoa, nwtom, lstab, nsewfr, kstart   &
             , kinput, koutpu, nitrap, ntmoy, kfond, kforc, mdforc, kavsmx !dmr ### , ntout

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "coetur"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8) :: vkappa, q2tmin, ghmax, ghmin, zlotur, vlmin, varfor, sqrghm
      integer      :: kajul

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "cerun"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      character(len=6) :: refexp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Additional stuff
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      real(kind=8), dimension(imax,jmax):: zmix_CC
      real(kind=8), dimension(imax,jmax,days_year360d_i) :: normUV_ERA5 = 0.0d0 ! Winds ERA5
      real(kind=8), dimension(imax,jmax,kmax) :: LFe_PISCES = 0.0d0 ! Iron Limitation PISCES

      end module bloc_mod

!--fin du fichier "bloc.com"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
