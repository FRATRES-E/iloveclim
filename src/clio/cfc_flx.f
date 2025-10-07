!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009

#define CFC_FLX 0

      SUBROUTINE cfc_flx(nn99)
#if ( CFC_FLX == 0 )
        integer, intent(in) :: nn99
#endif
#if ( CFC_FLX >= 1 )
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  This routine compute the air-sea exchange of CFC 11 & 12.
!  (Rq : initialement incorpore a "icdyna.Fom" et separe le 08/11/98)
!  modif : 22/04/99
!---

!! START_OF_USE_SECTION

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use ice_mod
      use dynami_mod
      use trace_mod
      use newunit_clio_mod, only: clio3_out_id
      
!! END_OF_USE_SECTION

!! START_OF_INCLUDE_SECTION

#include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "ice.com"
! [SCRPTCOM] #include "dynami.com"
! [SCRPTCOM] #include "trace.com"

!! END_OF_INCLUDE_SECTION

! Also define Coefficients for Schmidt numbers (Wanninkhof, JGR, 1992)
      parameter (a11=4039.8, b11=264.7, c11=8.2552, d11=0.10359)
      parameter (a12=3713.2, b12=243.3, c12=7.5879, d12=0.095215)

!----------------------------------------------------------------------
!   Set up CFC-fluxes from the atmosphere into the ocean, taking into
!  account the solubility of CFCs in seawater (Warner and Weiss, 1985).
!----------------------------------------------------------------------

      do 20 j=js1,js2
         do 10 i=is1(j),is2(j)
! Set integer variable to indicate whether SH or NH:
        if(j.lt.jeq) then
            ihemcfc = 1
        else
            ihemcfc = 2
        endif

! Set up CFC-11 and CFC-12 solubilities (tsol,ssol used for T,S
!       with T in degrees K, and S in parts per thousand):

!         mon = nint(tmonth)
!         tsol = sstobs(i,j,mon) + 273.0
!         ssol = (salobs(i,j,mon) * 1000.) + 35.0

          tsol = scal(i,j,ks2,1)
          ssol = scal(i,j,ks2,2)

          cfcsol11 = exp(a1cfc11 + a2cfc11*(100./tsol) +
     $          a3cfc11*log(tsol/100.) + a4cfc11*((tsol/100.)**2) +
     $          ssol*(b1cfc11 + b2cfc11*(tsol/100.) +
     $          b3cfc11*((tsol/100.)**2)) )
          cfcsol12 = exp(a1cfc12 + a2cfc12*(100./tsol) +
     $          a3cfc12*log(tsol/100.) + a4cfc12*((tsol/100.)**2) +
     $          ssol*(b1cfc12 + b2cfc12*(tsol/100.) +
     $          b3cfc12*((tsol/100.)**2)) )

!  Test output some CFC solubilities (near the dateline):
!         if(first .and. .not. mxpas2) then
!            if(i.eq.48 .and. j.eq.2) write(stdout,9088)
!9088         format(/,'Some examples of CFC solubilities (as a test):'
!    $,//,'   Temperature   Salinity   CFC-11*100  CFC-12*1000')
!            if((kmt(i,j).ge.1) .and. (i.eq.48))
!    $   write(stdout,9089) t(i,1,jc,nm,1), ssol,
!    $                          cfcsol11*100., cfcsol12*1000.
!         endif
!9089    format(2f11.3, 2f13.3)

!      if (i.eq.100.and.j.eq.50) then
!        write(120,*) tsol,ssol,cfcsol11,cfcsol12
!      endif

! Use a variety of techniques to get the ultimate "piston velocities":


!       cfcrestore = restore with e-folding time of 20 days (for example)
!               (i.e., surface level is restored towards atmospheric
!               concentration * solubility, factoring in a time scale
!               of 1/(20 days)

!  Value for the restoring
!      scalr(i,j,ks2,3) =cfcsol11*atmcfc11(ihemcfc)
!      scalr(i,j,ks2,4) =cfcsol12*atmcfc12(ihemcfc)
!      phiss(i,j,3) =0.0
!      phiss(i,j,4) =0.0

!       cfcwinds = Use complete Wanninkhof formulation with schmidt number
!                  and wind speed dependencies.  The wind speed dependent
!                  gas exchange is derived from E&K wind speed data (note
!                  interpolation over land and ice may give weird
!                  GFDL_prep_data values).

!       cfcliss = Use complete Liss and Merlivat (1986) formulation
!                  with schmidt number and wind speed dependencies.

!  Set up gammacfc to follow Wanninkhof (1992, JGR) formula:
        xsst = scal(i,j,ks2,1)-273.15


!  The Schmidt numbers (sc11, sc12) are slightly different (7% or so)
!    depending on which CFC is being considered:
        sc11 = a11 - b11*xsst + c11*(xsst**2) - d11*(xsst**3)
        sc12 = a12 - b12*xsst + c12*(xsst**2) - d12*(xsst**3)

! ifndef cfcliss
!  Wanninkhof's formula:
        xk11 = 0.31 * (vabq(i,j)**2) * sqrt(clio3_out_id0/sc11)
        xk12 = 0.31 * (vabq(i,j)**2) * sqrt(clio3_out_id0/sc12)
! endif

! ifdef cfcliss
!  Use Liss and Merlivat (1986) formulation instead ifdef cfcliss
!       if(vabq(i,j).le.3.6) xk11 = 0.17 * vabq(i,j)
!       if(vabq(i,j).gt.3.6 .and. vabq(i,j).le.13.)
!    $                  xk11 = 2.85 * vabq(i,j) - 9.65
!       if(vabq(i,j).gt.13.) xk11 = 5.9 * vabq(i,j) - 49.3
!       xk12 = xk11
!       if(vabq(i,j).le.3.6) then
!          xk11 = xk11 * ((sc11/600)**(-2./3.))
!          xk12 = xk12 * ((sc12/600)**(-2./3.))
!       else
!          xk11 = xk11 * ((sc11/600)**(-1./2.))
!          xk12 = xk12 * ((sc12/600)**(-1./2.))
!       endif
! endif

! Then the k (cm/hour) into time-scale (days) for level 1:
! Computation of the flux
!       t11(i,j) = (zw(1)/xk11) / 24.0
!       t12(i,j) = (zw(1)/xk12) / 24.0
        flcfc11  = xk11/(3600.*100.)*
     $    (cfcsol11*atmcfc11(ihemcfc) - scal(i,j,ks2,3))
        flcfc12  = xk12/(3600.*100.)*
     $    (cfcsol12*atmcfc12(ihemcfc) - scal(i,j,ks2,4))
        phiss(i,j,3) = - flcfc11 * dts(ks2) * unsdz(ks2)
        phiss(i,j,4) = - flcfc12 * dts(ks2) * unsdz(ks2)
!------------------------------------------------------------------------
! NB: Now limit the air-sea CFC flux in direct proportion to the
!     percentage coverage of sea-ice:
!------------------------------------------------------------------------
        phiss(i,j,3) = phiss(i,j,3) * albq(i,j)
        phiss(i,j,4) = phiss(i,j,4) * albq(i,j)

!  Then just add the CFC fluxes to the surface level CFC-11 and CFC-12....
 10      continue
 20   continue

      return
#endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine cfc_flx -
      end
