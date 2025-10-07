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

      module ice_mod

      use para0_mod, only: imax, jmax, kmax
      use para_mod,  only: nbsmax, nkb0

#if ( ISOOCN >= 1 )
      use para0_mod, only: owatert
#endif

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! originally written with IMPLICIT INTEGER (I-N), REAL (A-H, O-Z)
! transferred as such in fortran90
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  bloc "ice.com" : 'include' in the routines linked with the ice
!  modif : 25/09/98
!  modif : 18/02/05 AnneM
!  modif : 04/10/11 DidierMR

! tfsn      Melting point temperature of the snow
! tfsg      Melting point temperature of the ice
! xkn       Conductivity of the snow
! xkg       Conductivity of the ice
! rcpn      Density times specific heat for the snow
! rcpg      Density times specific heat for the ice
! xlg       Latent heat of fusion of the ice
! xln       Latent heat of fusion of the snow
! rhog      Density of the ice
! rhon      Density of the snow
! emig      Emissivity of the ice
! sglace    Salinity of the ice
! hmelt     Maximum melting at the bottom
! acrit(2)  Minimum fraction for leads
! hgcrit(2) Ice thickness for lateral accretion
! hgmin     Ice thickness corr. to max. energy stored in brine pocket
! hndif     Computation of temp. in snow or not
! hgdif     Computation of temp. in ice or not
! hglim     Minimum ice thickness
! amax      Maximum lead fraction
! uscomi    =1.0/(1.0-amax)
! beta      Numerical caracteritic of the scheme for diffusion in ice
! ddtb      Time step for ice thermodynamics (s)
! swiqst    Energy stored in brine pocket or not
! parlat    Percentage of energy used for lateral ablation
! hakspl    Slope of distr. for Hakkinen-Mellor's lateral melting
! hibspl    Slope of distribution for Hibler's lateral melting
! exld      Exponent for leads-closure rate
! hakdif    Coefficient for diffusions of ice and snow
! hth       Threshold thickness for comp. of eq. thermal conductivity
! hnzst     Thickness of the surf. layer in temp. computation
! parsub    Switch for snow sublimation or not
! alphs     Used to take into account settling during snow-ice formation
! cnscg     ratio  rcpn/rcpg
! nbits     Number of time steps in Newton -Raphson procedure
! stefan    Stefan-Boltzman constant
! xsn       Sublimation heat for the snow
! too       temperature of the triple point for water
! vkarmn    von Karman constant
! cevap     Latent heat of evaporation of water
! zemise    Emissivity of water
! rhoesn    1/rhon
! firg      IR flux over the ice (only used for outputs)
! fcsg      Sensible heat flux over the ice (only used for outputs)
! fleg      Latent heat flux over the ice (only used for outputs)
! dvosbq    Variation of volume at surface (only used for outputs)
! dvobbq    Variation of ice volume at the bottom ice (only used for outputs)
! dvolbq    Total variation of ice volume (only used for outputs)
! dvonbq    Lateral Variation of ice volume (only used for outputs)
! ts        Surface temperature of the ice
! tfu       Melting point temperature of sea water
! hnbq      Snow thickness
! hgbq      Ice thickness
! hgbqp     Ice production/melting
! albq      Leads fraction
! qstobq    Energy stored in the brine pockets
! fbbq      Heat flux at the ice base
! tbq       Temperature inside the ice/snow layer
! dmnbq     Variation of snow mass
! dmgbq     Variation of ice mass
! qlbq      heat balance of the lead (or of the open ocean)
! qcmbq     Energy needed to bring the ocean surface layer until its freezing
!           point (at a factor 2)
! thcm      part of the solar energy used in the lead heat budget
! fstrbq    Solar flux transmitted trough the ice
! ffltbq    Array linked with the max heat contained in brine pockets (?)
! fscmbq    Linked with the solar flux below the ice (?)
! fsbbq     Also linked with the solar flux below the ice (?)
! qfvbq     Array used to store energy in case of toral lateral ablation (?)
! xzo       rugosity of the ice (no more used)
! dmgwi     Variation of the mass of snow ice
! psbq      Surface air pressure
! tabq      Surface air temperature
! qabq      Surface air humidity
! vabq      Surface wind velocity
! hnplbq    Snow precipitation
! fevabq    Evaporation flux
! fsolcn    Solar flux at the ocean surface
! fsolg     Solar flux at the ice surface
! flecn     Latent heat flux at the ocean surface
! fcscn     Sensible heat flux at the ocean surface
! tenagx    Wind stress at the ice surface (x)
! tenagy    Wind stress at the ice surface (y)
! albege    Albedo of the snow or ice (only for outputs)
! tairox    Wind stress at the ocean surface (x)
! tairoy    Wind stress at the ocean surface (y)
! ratbqg    Longwave downward radiation flux over the ice
! ratbqo    Longwave downward radiation flux over the ocean
! cloud     Cloud fraction
! tdew      Air relative humidity
! albecn    Albedo of the ocean (only for outputs)
! tauc      Cloud optical depth
! runoff    river runoff
! sdvt      u*^2/(Stress/density)
! fcm1      Solar flux at the ocean surface or the ocean-ice interface
! fcm2      Non-solar flux at the ocean surface
! fwat      Freshwater flux (change of definition between the routines)
! reslum    Relative absorption of solar radiation in each ocean level
! alct      lead opening/closure rate due to thermodynamics 
! alcd      lead opening/closure rate due to convergence/divergence 
! alcr      lead opening/closure rate due to shear motion if Clds is activated
! hgcol     Collection thickness of sea ice in leads when it is variable (Cvhg enabled)
! ficebergn Heat due to iceberg melting in the Northern Hemisphere
! ficebergs Heat due to iceberg melting in the Southern Hemisphere
! iiceber, jiceber index for location of iceberg melting
! areiceb   surface for iceberg melting
! tmics     Mask for ice shelve melting
! toticesm  total melting due to ice shelve melting
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "cstbg"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer                                    :: nbits
      real(kind=8)                               :: tfsn, tfsg, xkn, xkg, rcpn, rcpg, xlg, xln, rhog, rhon, emig, sglace, hmelt    &
                  , hgmin, hndif, hgdif, hglim, amax, uscomi, beta, ddtb, swiqst, parlat, hakspl, hibspl, exld, hakdif, hth, hnzst &
                  ,parsub, alphs, cnscg
      real(kind=8), dimension(2)                 :: acrit, hgcrit

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "fluxsf"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8)                               :: stefan, xsn, too, vkarmn, cevap, zemise, rhoesn

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "comdia"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax)         :: firg = 0.0d0, fcsg = 0.0d0, fleg = 0.0d0, dvosbq, dvobbq, dvolbq, dvonbq

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "comban"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax)         :: ts, tfu, hnbq, hgbq, hgbqp = 0.0d0, qstobq, fbbq = 0.0d0, dmnbq, dmgbq, qlbq, &
                    qcmbq, thcm, fstrbq, ffltbq, fscmbq, fsbbq, qfvbq, xzo, dmgwi, alct, alcd, alcr, hgcol
      real(kind=8), dimension(imax,jmax,nkb0)    :: tbq

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "comban2"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax)         :: albq

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "comfor"
!-uwind et vwind en common pour les icebergs
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      real(kind=8), dimension(imax,jmax)         :: psbq, tabq = 0.0d0, qabq = 0.0d0, vabq, fevabq, fsolcn = 0.0d0, fsolg, flecn& 
                                                  , fcscn & 
                                                  , tenagx &
                                                  , tenagy, albege = 0.0d0, ratbqg, ratbqo, cloud, tdew, uwind, vwind           &
                                                  , albecn = 0.0d0, tauc     &
                                                  , runoff, sdvt=0.0d0
      real(kind=8), dimension(imax,jmax)         :: tairox = 0.0d0, tairoy = 0.0d0 ! [NOTA] dummy init

#if ( ISOOCN >= 1 )
      real(kind=8), dimension(imax,jmax,owatert) :: hnplbq
#else
      real(kind=8), dimension(imax,jmax)         :: hnplbq
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "comca"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax)         :: fcm2
#if ( ISOOCN >= 1 ) 
      real(kind=8), dimension(imax,jmax,owatert) :: fwat
#else
      real(kind=8), dimension(imax,jmax)         :: fwat, fwruno
#endif
      real(kind=8), dimension(imax,jmax,0:kmax+1):: reslum

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "comca2"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax)         :: fcm1 = 0.0d0


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "comiceb"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8)                               :: ficebergn,ficebergs, areicebn, areicebs, toticesm
      integer                                    :: iicebern1, iicebern2, jicebern1, jicebern2, iicebers1, iicebers2,jicebers1     &
             , jicebers2
      real(kind=8), dimension(imax,jmax)         :: tmics

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "comstr"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=4), dimension(57,kmax+1,0:nbsmax):: vwx

!age  common / sorage/ agen(imax,jmax),ageg(imax,jmax)
!cp2  common / fderice / fder(imax,jmax)
      end module ice_mod
!
!--fin du fichier "ice.com"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
