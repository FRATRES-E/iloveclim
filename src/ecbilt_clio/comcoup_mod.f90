!dmr -- Added optional components choice - Tue Jun 14 12:07:40 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 12:07:40 CEST 2011

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module est la conversion en module FORTRAN90/95 du fichier
!       comcoup.h
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche (pour la conversion)
!      Date   : 06 juin 2018
!      Derniere modification : 06 juin 2018, 17 juin 2022
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module comcoup_mod

      use global_constants_mod, only: dblp=>dp, silp=>sp, ip
      use comatm, only: nlat, nlon, nwisos

!~ #if ( ISOATM >=1 )
!~       use iso_param_mod, only: neauiso
!~ #endif
!~ #if ( ISOOCN >= 1 )
!~       use para0_mod, only: owatert
!~ #endif

      implicit none

      private :: nlat, nlon

!23456789012345678901234567890123456789012345678901234567890123456789012
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** File:     comcoup.h
! *** Contents: Common declarations for coupling module of ECbilt
!dmr @-@ iceb0
! JONO march 2004 introducing wind vector to drive icebergs;
! 'coupling sum' of total wind at 10 meter;
! sumu10(),sumv10() = sum of utot10(),vtot10() = utot3,vtot3 * uv10rwv
!dmr @-@ iceb0
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


#if ( HRCLIO == 0 )
      integer(ip), parameter :: ijatm=nlon*nlat, ijocn=122*65, komax=17
#if ( BATHY >=2 )
      integer(ip)            :: kamax
#else
      integer(ip)            :: kamax=14
!      integer(ip)            :: kamax=13 !from 7.5ka
#endif

!---dmr [SUPPRESS]      integer(ip), parameter :: kamax2=14
#elif ( HRCLIO == 1 )
      integer(ip), parameter :: ijatm=nlon*nlat, ijocn=242*128,komax=17
      integer(ip)            :: kamax=38
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      integer(ip)            :: iobclint,iobtropt
!---dmr [KAMAX_VAR]      integer(ip), dimension(ijatm,kamax) :: indo2a
      integer(ip), dimension(:,:), ALLOCATABLE :: indo2a
!---dmr [SUPPRESS]      integer(ip), dimension(ijatm,kamax2) :: indo2a2
      integer(ip), dimension(ijocn,komax) :: inda2o

      real(dblp),  dimension(nlat,nlon)   :: cohesws, clhesws, cohesw0, clhesw0, cohesw1, clhesw1, cohesw2, clhesw2
      real(dblp),  dimension(nlat,nlon)   :: coulrads,clulrads,coulrad0,clulrad0,coulrad1,clulrad1,coulrad2,clulrad2
      real(dblp),  dimension(nlat,nlon)   :: codlrads,cldlrads,coeflux,cleflux,cohflux,clhflux

!~ #if ( ISOATM >= 1 )
!~       real(dblp),  dimension(nlat,nlon,neauiso):: coevap, clevap
!~ #else
      real(dblp),  dimension(nlat,nlon,nwisos)   :: coevap = 0.0_dblp , clevap = 0.0_dblp
!~ #endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(dblp),  dimension(nlat,nlon)   :: winstua, winstva, sumohfx, sumoswr, sumihfx, sumiswr

!~ #if ( ISOATM >= 2 && ISOOCN >= 1 )
!~       real(dblp),  dimension(nlat,nlon,owatert):: sumofwf, sumisno
!~ #else
      real(dblp),  dimension(nlat,nlon,nwisos)   :: sumofwf, sumisno
!~ #endif
      real(dblp),  dimension(nlat,nlon)   :: sumty, sumtx, sumuv10, couptcc=0.5_dblp

!~ #if ( ISOATM >= 2 )
!~       real(dblp),  dimension(nlat,nlon,neauiso):: couprf, coupsf
!~ #else
      real(dblp),  dimension(nlat,nlon,nwisos)   ::  couprf = 0.0_dblp, coupsf = 0.0_dblp
!~ #endif
      real(dblp),  dimension(nlat,nlon)   :: samix

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~ #if ( ISOATM >= 2 && ISOOCN >= 1 )
!~       real(dblp),  dimension(nlat,nlon,owatert):: sumrl, sumro
!~ #else
      real(dblp),  dimension(nlat,nlon,nwisos)   :: sumrl = 0.0_dblp, sumro = 0.0_dblp
!~ #endif
      real(dblp)                          :: sumhsn,sumhss,sumohsn,sumohss

!~ #if ( ISOATM >= 2 && ISOOCN >= 1 )
!~       real(dblp),  dimension(nlat,nlon,owatert):: couprunl, coupruno
!~ #else
      real(dblp),  dimension(nlat,nlon,nwisos)   :: couprunl = 0.0_dblp, coupruno = 0.0_dblp
!~ #endif

#if ( ICEBERG == 2 && ISM != 2 )
      real(dblp),  dimension(nlat,nlon)   :: coupiceb,sumiceb
#endif

      real(dblp) :: couphsnn,couphsns
      real(dblp),  dimension(nlat,nlon)   :: sumicof, sumpress

!---dmr [KAMAX_VAR]      real(dblp),  dimension(ijatm,kamax) :: wo2a
      real(dblp), dimension(:,:), ALLOCATABLE :: wo2a
!---dmr [SUPPRESS]      real(dblp),  dimension(ijatm,kamax2) :: wo2a2
      real(dblp),  dimension(ijocn,komax) :: wa2o

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



!!! [DEPRECATED] 2019-09-27
!!!### #if ( WATGRISCONS == 1 )
!!!### !mdr --- added water conservation between ECBilt and GRISLI
!!!###       real(dblp),  dimension(nlat,nlon)   :: coupruncorr, sumuncorr, sumouncorr
!!!### 
!!!### #endif
!!! [DEPRECATED] 2019-09-27
!dmr @-@ iceb0
      real(dblp),  dimension(nlat,nlon)   :: sumu10, sumv10
!dmr @-@ iceb0


!~       common /ec_rcoup/winstua,winstva,sumtx,sumty,                              &
!~      &      couprf,coupsf,couptcc,samix,                                         &
!~      &      cohesws,clhesws,cohesw0,clhesw0,cohesw1,clhesw1,                     &
!~      &      cohesw2,clhesw2,coulrads,clulrads,coulrad0,clulrad0,                 &
!~      &      coulrad1,clulrad1,coulrad2,clulrad2,codlrads,cldlrads,               &
!~      &      coeflux,cleflux,cohflux,clhflux,coevap,clevap,                       &
!~      &      sumrl,sumro,couprunl,coupruno,couphsnn,couphsns,                     &
!~      &      sumohfx,sumoswr,sumihfx,sumiswr,sumofwf,sumisno,                     &
!~      &      sumhsn,sumhss,sumohsn,sumohss,sumicof,sumuv10,sumpress               &
!~ #if( UNCORRUNOFF == 1 )
!~      &      ,coupgismelt,sumgismelt,sumogismelt                                  &
!~ #endif
!~ #if ( WATGRISCONS == 1 )
!~      &      ,coupruncorr,sumuncorr,sumouncorr                                    &
!~ #endif

!~ !dmr @-@ iceb0
!~      &      ,sumu10, sumv10
!~ !dmr @-@ iceb0

!~       common /ec_icoup/ iobtropt,iobclint,indo2a,inda2o

!~       common /ec_gridinfo/ wo2a,wa2o

      end module comcoup_mod
