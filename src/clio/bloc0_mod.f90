!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module est un portage en FORTRAN 90 du module CLIO initial
!       cree en FORTRAN 77.
!      (dans l'environnement logiciel LUDUS)
!
!      Auteur : ??, Didier M. Roche
!      Date   : ??, 19 Aout 2014
!      Derniere modification : 19 Juin 2018
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module bloc0_mod

      use para0_mod, only: imax, jmax, kmax, nsmax, nrpmax

#if ( ISOOCN >= 1 )
      use para0_mod, only: owatert
#endif

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  fichier "bloc0.com" : incorpore dans bloc.com,loch.com
!  creation : 10/02/2004 <- LOCH

!-- Bloc atm2loch

      real(kind=8), dimension(imax, jmax)        :: temp_vent,temp_press

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     real(kind=8)                                :: dx, dy, unsdx, unsdy, uns2dx, uns2dy
!!!     real(kind=8), dimension(kmax)               :: dz, unsdz
     real(kind=8), dimension(kmax)               :: unsdz
     real(kind=8), dimension(kmax),target               :: dz
     real(kind=8), dimension(kmax+1)             :: z, zw, dzw, unsdzw
     real(kind=8), dimension(imax,jmax)          :: huy, hux, hu, unshu
     real(kind=8), dimension(imax,jmax)          :: xslon, xulon, yslat, yulat
     real(kind=8), dimension(imax,jmax)          :: hs, angle
     real(kind=8), dimension(imax,jmax,4)        :: xsedg, xuedg, ysedg, yuedg
     real(kind=8), dimension(imax+1,jmax+1)      :: xslonp, xulonp, yslatp, yulatp
!!!     real(kind=8), dimension(imax,jmax,kmax)     :: tms, tmu
     real(kind=8), dimension(imax,jmax,kmax)     :: tmu
     real(kind=8), dimension(imax,jmax,kmax),target     :: tms
!DFG

#if ( OCYCC == 1 )
!dmr --- Ajout pour recuperation des champs OCYCC autres que scalaires (pour sorties)
     real(kind=8), dimension(imax,jmax)           :: oxpco2_clio = 0.0d0, TPPma_clio = 0.0d0, CACO3ma_clio = 0.0d0
     real(kind=8), dimension(imax,jmax,kmax)      :: OPOCFlux_clio
#endif

#if ( OOISO == 1 )
     real(kind=8), dimension(imax,jmax,kmax)      :: phyto_clio, prodO2_clio, respO2_clio, reminO2_clio, NCP_clio
#endif

      integer, dimension(imax,jmax)               :: msks, msku

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     real(kind=8)                                          :: tpstot
     real(kind=8), dimension(imax,jmax)                    :: daeta
     real(kind=8), dimension(imax,jmax)                    :: eta,ub,vb
     real(kind=8), dimension(imax,jmax,kmax)               :: u, v,b=0.0d0, bvf, avsdz, avudz, fqajc
     real(kind=8), dimension(imax,jmax,kmax+1)             :: w,q2turb
     real(kind=8), dimension(imax,jmax,0:nsmax)            :: fss
     real(kind=8), dimension(imax,jmax,kmax,nsmax)         :: scal

#if ( CORAL == 1 )
     real(kind=8), dimension(imax,jmax,kmax)         :: coral_area_clio
     real(kind=8), dimension(imax,jmax,kmax)         :: coral_prod_clio
     real(kind=8), dimension(imax,jmax,kmax)         :: coral_mass_clio
     real(kind=8), dimension(imax,jmax,kmax)         :: omega_arag_clio
     real(kind=8), dimension(imax,jmax,kmax)         :: oco3_clio
     real(kind=8), dimension(imax,jmax,kmax)         :: tau_bleach_clio
     real(kind=8), dimension(imax,jmax,kmax)         :: DHW_clio
     real(kind=8), dimension(imax,jmax,kmax)         :: PH_clio
#endif
#if ( OCYCC == 1 )
      real(kind=8), dimension(imax, jmax, kmax) :: fPOC_flx_clio, fCAL_flx_clio
#endif
#if ( REMIN == 1 )
     real(kind=8), dimension(imax,jmax,kmax)         :: kremin_clio
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "forcing"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     real(kind=8)                                          :: spvr, ahrap
     real(kind=8), dimension(imax,jmax)                    :: phifu, phifv, phisu, phisv, ust2s = 0.0d0, ust2b
     real(kind=8), dimension(kmax,nsmax)                   :: scal0
     real(kind=8), dimension(nsmax)                        :: scpme, scssv
     real(kind=8), dimension(imax,jmax,nsmax)              :: scs
     real(kind=8), dimension(imax,jmax,kmax)               :: rappel
     real(kind=8), dimension(nrpmax,kmax,nsmax)            :: rapint
     real(kind=8), dimension(imax,jmax,0:nsmax)            :: rappes,phifs,phiss
     real(kind=8), dimension(imax,jmax,kmax,nsmax)         :: scalr, phivs
     real(kind=8), dimension(imax,jmax,0:1,0:nsmax)        :: phimnx

#if( APPLY_UNCORFWF == 1 )
     real(kind=8), dimension(imax,jmax)                    :: w_uncor
#endif
#if( UNCORRUNOFF == 1 )
     real(kind=8), dimension(imax,jmax)                    :: fwatgis
#endif
!---dmr&afq Moved to varsCONSEAU_mod.f90
!--- #if ( CONSEAU == 1 ) /* Water conservation flag between ECBilt and GRISLI */
!---      real(kind=8), dimension(imax,jmax)                    :: fwatconseau
!---      real(kind=8)                                          :: trendconseau
!--- #endif

#if ( F_PALAEO_FWF == 2 )
     real(kind=8), dimension(imax,jmax)                    :: fwatconseau
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "lerun2"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer :: nitrun,nsav


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "limites"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer                           :: ks1, ks2, ku1, ku2, ims1, ims2, js1, js2, imu1, imu2, ju1, ju2, jcl1, jcl2, jeq, jdl1   &
             , jdl2, ijsdl, ijudl, iberp, jberp, ibera, jbera, iberpm, jberpm, iberam, jberam
      integer, dimension(jmax)          :: is1, is2, iu1, iu2, isf1, isf2, iuf1, iuf2
      integer, dimension(imax,jmax,-1:1):: kniv

      integer, dimension(imax,jmax), target :: kfs
      integer, dimension(imax,jmax) :: kfu

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "runpara"(meters)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8)                  :: dtu, dtb, cdbot, ahu, ahe, alphxu, alphxv, alphyu, alphyv, afilt, ahh, avv, rifsmx, rifumx  &
                  , qavs, avsn2, ccfmn, ccfmx, txiadu, txiads, txidfu, txidfs, txifcb, txifcu, bering, ajcmix, xslop, bering_prev

      real(kind=8), dimension(2)    :: alphah
      real(kind=8), dimension(kmax) :: slopemax, ahs, ai, aitd, slopmgm, alphmi, alphaz, avnu0, avnub, avk0, avkb
      real(kind=8), dimension(nsmax):: alphgr, algrmn, txeflx


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "runpara"(meters) [bis]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(kmax)     :: dts

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "fluxsel"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      real(kind=8), dimension(imax,jmax):: frapsav


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "volcor"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer                      :: ihyster
      real(kind=8)                 :: zurfow,zflux0,zfluxm,vcor,bheat
      real(kind=8), dimension(2500):: xfreshsv

#if ( ISOOCN >= 1 && WISOATM == 1 )
      real(kind=8), dimension(owatert)   :: zflux0s,zfluxms
#else
      real(kind=8)                 :: zflux0s,zfluxms
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!--fin du fichier "bloc0.com"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      end module bloc0_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
