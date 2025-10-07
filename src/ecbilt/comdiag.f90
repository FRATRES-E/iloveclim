!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module est la conversion en module FORTRAN90/95 du fichier
!       comatm.h
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche (pour la conversion)
!      Date   : 16 Mai 2018
!      Derniere modification : 16 Mai 2018
!-----|--1--------2---------3---------4---------5---------6---------7-|
!dmr -- Added optional components choice - Tue Jun 14 12:07:40 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 12:07:40 CEST 2011
!23456789012345678901234567890123456789012345678901234567890123456789012
!-----------------------------------------------------------------------
! *** File:     comdiag.h
! *** Contents: Common declarations for diagnostics of ECbilt
!-----------------------------------------------------------------------
! *** COMMON  /writectl/ixout,ioutdaily,ioutyearly,itel,minterv,
!                      meantype,meantot,meanyl,ifrendat,instcount
!     ixout:   output frequency for instantanous field in days
!     ioutdaily: output (1) or not (0) instantanous fields.
!     ioutyearly: output (1) or not (0) yearly mean fields.
!     itel:    counts the number of days in the output interval
!     minterv: counter used to compute monthly or seasonal mean.
!     meantype = 1 monthly mean
!     meantype = 2 seasonal mean.
!     meantot = 1 output total integration period monthly or seasonal
!               mean
!     meanyl  = 1 output yearly monthly or seasonal mean.
!     ifrendat: frequency of writing restart data in days.
!     instcount: counter used for output instantaneous fields.
!     irad =1 if calculation of radiative forcings
!-----------------------------------------------------------------------

      MODULE comdiag

      use global_constants_mod, only: dblp=>dp, silp=>sp, ip

      use comatm, only: nlat, nlon, nvl

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer(ip) ::     ixout,ioutdaily,ioutyearly,itel,minterv,ntotdays, meantype,meantot,meanyl,ifrendat,instcount,irad

!~       common /ec_writectl/ixout,ioutdaily,ioutyearly,itel,minterv,
!~      &                 meantype,meantot,meanyl,ifrendat,instcount,irad
!~ #if ( PMIP == 1 )

      integer(ip)                        :: numtotvar
      integer(ip), dimension(4)          :: thirdd
      real(dblp)                         :: fill_value, missing_value
      real(dblp), dimension(80,20)       :: newtotvar
      character(len=60), dimension(80,6) :: nametotvar

!~       common /ec_netcdf/newtotvar,fill_value,missing_value,numtotvar,
!~      &                  nametotvar,thirdd
!~ #endif

      integer(ip), dimension(20)         :: newts, newt2m, newt, newtstrat, newu, newv, newpsi, newhforg, newvforg, newhflux       &
                                          , newdyrain, newcorain, newrmoisg, newrelhum, newdivv, newssr, newtsr, newstr, newttr    &
                                          , newbmoisg, newtorain, newevap, neweminp, newalbs, newustress, newvstress, newsdl       &
                                          , newrunoffo, newrunoffl, neweflux, newdivu, newuv10, newchi, newsp, newalbp, newgh      &
                                          , newcdragw, newcdragv, newqgpv, newtcc, newdumt1, newdumt2, newdumu1, newdumu2          &
                                          , newomega, newsnow, newhic, newswrs, newlwrs
#if ( ISOATM >= 1 )
      integer(ip), dimension(20)         :: newpp18, newsn18, newpp17, newsn17, newppd, newsnd, newrmoisg18, newrmoisg17, newrmoisgd
#endif
#if ( ISOATM >= 2 )
      integer(ip), dimension(20)         :: newbmoisg18, newsdl18, newbmoisg17, newsdl17, newbmoisgd, newsdld
#endif

      real(dblp), dimension(nlat,nlon)   :: s1u200, s2u200, s1u500, s2u500, s1u800, s2u800, s1temp4g, s2temp4g, s1temp2g, s2temp2g &
                                          , s1tstrat, s2tstrat, s1tsurf, s2tsurf, s1tempsg, s2tempsg, s1t2m, s2t2m, s1v200, s2v200 &
                                          , s1v500, s2v500, s1v800, s2v800, s1vhforg1, s2vhforg1, s1vhforg2, s2vhforg2, s1vforg1   &
                                          , s2vforg1, s1pground, s2pground, s1vforg2,s2vforg2, s1vforg3,s2vforg3, s1dyrain         &
                                          , s2dyrain, s1corain, s2corain, s1snow, s2snow, s1Hflux, s2Hflux, s1Eflux, s2Eflux       &
                                          , s1hesw, s2hesw, s1hesws, s2hesws, s1ulrad1, s2ulrad1, s1nlrads, s2nlrads, s1swrs       &
                                          , s2swrs, s1lwrs, s2lwrs, s1bmoisg, s2bmoisg, s1rmoisgw3, s2rmoisgw3, s1relhum, s2relhum &
                                          , s1torain, s2torain, s1evap, s2evap, s1eminp, s2eminp, s1albes, s2albes, s1albep        &
                                          , s2albep,  s1winstu1, s2winstu1, s1winstv1, s2winstv1, s1dsnow, s2dsnow, s1hic, s2hic   &
                                          , s1runofo, s2runofo, s1runofl, s2runofl, s1uv10, s2uv10, s1cdragw, s2cdragw, s1cdragv   &
                                          , s2cdragv, s1tcc, s2tcc

      real(dblp), dimension(nlat,nlon,nvl)::s1omeg, s2omeg, s1psi, s2psi, s1udivg, s2udivg, s1vdivg, s2vdivg, s1chi, s2chi, s1gh   &
                                          , s2gh, s1qgpv, s2qgpv, s1dt1, s2dt1, s1dt2, s2dt2, s1du1, s2du1, s1du2, s2du2

#if ( ISOATM >= 1 )
      real(dblp), dimension(nlat,nlon)   :: s1torain18, s2torain18, s1tosnow18, s2tosnow18, s1torain17, s2torain17, s1tosnow17     &
                                          , s2tosnow17, s1toraind, s2toraind, s1tosnowd, s2tosnowd, s1rmoisg18, s2rmoisg18         &
                                          , s1rmoisg17, s2rmoisg17, s1rmoisgd, s2rmoisgd, s1dsnow18, s2dsnow18, s1dsnow17          &
                                          , s2dsnow17,  s1dsnowd, s2dsnowd, s1bmoisg18, s2bmoisg18, s1bmoisg17, s2bmoisg17         &
                                          , s1bmoisgd, s2bmoisgd
#endif

      real(dblp), dimension(nlat,nlon,12) :: sxu200, syu200, sxu500, syu500, sxu800, syu800, sxomeg3, syomeg3, sxomeg2, syomeg2    &
                                           , sxomeg1, syomeg1, sxtemp4g, sytemp4g, sxtemp2g, sytemp2g, sxtstrat, sytstrat, sxtsurf &
                                           , sytsurf, sxtempsg, sytempsg, sxt2m, syt2m, sxv200, syv200, sxv500, syv500, sxv800     &
                                           , syv800, sxgrpsi3, sygrpsi3, sxgrpsi2, sygrpsi2, sxgrpsi1, sygrpsi1, sxpground         &
                                           , sypground, sxvhforg1, syvhforg1, sxvhforg2, syvhforg2, sxvforg1, syvforg1, sxvforg2   &
                                           , syvforg2, sxvforg3, syvforg3, sxudivg3, syudivg3, sxudivg2, syudivg2, sxudivg1        &
                                           , syudivg1, sxvdivg3, syvdivg3, sxvdivg2, syvdivg2, sxvdivg1, syvdivg1, sxdyrain        &
                                           , sydyrain, sxcorain, sycorain, sxhflux, syhflux, sxeflux, syeflux, sxhesw, syhesw      &
                                           , sxhesws, syhesws, sxulrad1, syulrad1, sxnlrads, synlrads, sxswrs, syswrs, sxlwrs      &
                                           , sylwrs, sxbmoisg, sybmoisg, sxrmoisgw3, syrmoisgw3, sxrelhum, syrelhum, sxtorain      &
                                           , sytorain, sxsnow, sysnow, sxevap, syevap, sxeminp, syeminp, sxalbes, syalbes, sxalbep &
                                           , syalbep, sxwinstu1, sywinstu1, sxwinstv1, sywinstv1, sxdsnow, sydsnow, sxhic, syhic   &
                                           , sxrunofo, syrunofo, sxrunofl, syrunofl, sxuv10, syuv10, sxchi3, sychi3, sxchi2        &
                                           , sychi2, sxchi1, sychi1, sxgh3, sygh3, sxgh2, sygh2, sxgh1, sygh1, sxqgpv3, syqgpv3    &
                                           , sxqgpv2, syqgpv2, sxqgpv1, syqgpv1, sxcdragw, sycdragw, sxcdragv, sycdragv, sxtcc     &
                                           , sytcc, sxdt13, sydt13, sxdt12, sydt12, sxdt11, sydt11, sxdt23, sydt23, sxdt22, sydt22 &
                                           , sxdt21, sydt21, sxdu13, sydu13, sxdu12, sydu12, sxdu11, sydu11, sxdu23, sydu23, sxdu22&
                                           , sydu22, sxdu21, sydu21
#if ( ISOATM >= 1 )
      real(dblp), dimension(nlat,nlon,12) :: sxtorain18, sytorain18, sxtosnow18, sytosnow18, sxtorain17, sytorain17, sxtosnow17    &
                                           , sytosnow17, sxtoraind , sytoraind , sxtosnowd , sytosnowd , sxrmoisg18, syrmoisg18    &
                                           , sxrmoisg17, syrmoisg17, sxrmoisgd , syrmoisgd
#endif
#if ( ISOATM >= 2 )
      real(dblp), dimension(nlat,nlon,12) :: sxdsnow18, sydsnow18, sxdsnow17, sydsnow17, sxdsnowd, sydsnowd, sxbmoisg18,sybmoisg18 &
                                           , sxbmoisg17,sybmoisg17,sxbmoisgd, sybmoisgd
#endif

      integer(ip), parameter             :: nbasa = 8
      real(dblp),  dimension(nbasa)      :: precan, evapan, runoan, arocbasa
      real(dblp),  dimension(nlat,nbasa) :: hfmeanan
      integer(ip), dimension(nlat,nlon)  :: iocbasa


      END MODULE comdiag
