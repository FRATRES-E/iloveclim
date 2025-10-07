!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2024 Didier M. Roche (a.k.a. dmr)

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:50 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:50 CET 2009


      MODULE mod_bgc_fluxes_OCNATM

       USE global_constants_mod, ONLY: dblp=>dp, ip, days_year360d_i
       USE declars_mod,          ONLY: NOC_CBR, JT, LT
       USE iso_dioxygen_mod,     ONLY: iair
       USE para0_mod,            ONLY: NISOO2

       IMPLICIT NONE

       PRIVATE

       PUBLIC :: OCN_BIO_FLUXES

#if ( WINDINCC == 0 )
      REAL(kind=dblp)  :: fco2ex
      parameter (fco2ex=6.12e-5) ! m/s
#endif
       CONTAINS


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE OCN_BIO_FLUXES(TM, SM, FRICE, OO2, O2_sat_thistime,                                      &
                         oxpCO2, osCO2, oxCO2, oxHCO3, oxCO3, ODIC, FOC14, FOALK, FODOC, FODOC13, FODOCS, & 
                         FODOCS13, FONO3, FOO2, FOPO4, FODIC, FOC13, OC14, OALK, ODOC, ODOC13, ODOCS,     &
                         ODOCS13, ONO3, OPO4, OC13, ON2O)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!  By: J. Bendtsen
!  Last modification: 2024.10.29 (dmr,ecl)
!
!  Purpose: main routine for ocean biology and biogeochemistry fluxes to/from the atmosphere
!  Modified version by jwyang @Oct 2020
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       use marine_bio_mod, only: JPROD
       use marine_bio_mod, only: FDIC, FOAC14, ODIC_diff

      IMPLICIT NONE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr [nb,tb] By reference variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!  variables as INTENT(IN)
      
      ! dmr Variables from the physics: temperature, salinity, sea-ice
      real(kind=dblp), dimension(LT,JT,NOC_CBR), intent(in) :: TM            ! from loveclim_transfer_mod
      real(kind=dblp), dimension(LT,JT,NOC_CBR), intent(in) :: SM            ! from loveclim_transfer_mod
      real(kind=dblp), dimension(LT, NOC_CBR), intent(in)   :: FRICE         ! from loveclim_transfer_mod
      
      ! Used by Compute_oxy_fluxes and Compute_oxy_atm_fluxes SUBROUTINE 
      real(kind=dblp), dimension(LT, JPROD, NOC_CBR), intent(in) :: O2_sat_thistime ! from O2SAT_mod
       
      ! Used by Compute_pCO2_fluxes SUBROUTINE 
      real, dimension(LT, NOC_CBR), intent(in):: oxpCO2                  ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(in):: osCO2                   ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(in):: oxCO2                   ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(in):: oxHCO3                  ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(in):: oxCO3                   ! from marine_bio_mod

!  variables as INTENT(INOUT)
      real(kind=dblp), dimension(LT, JT, NOC_CBR, NISOO2), INTENT(inout) :: OO2
      real(kind=8), dimension(LT, JT, NOC_CBR), intent(inout)            :: ODIC ! from marine_bio_mod

      !dmr --- Surface fluxes (FO**)
      real, dimension(LT, NOC_CBR), intent(inout):: FOALK                ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(inout):: FODOC                ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(inout):: FODOC13              ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(inout):: FODOCS               ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(inout):: FODOCS13             ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(inout):: FONO3                ! from marine_bio_mod
      real, dimension(LT, NOC_CBR, NISOO2), intent(inout):: FOO2         ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(inout):: FOPO4                ! from marine_bio_mod
      
      real, dimension(LT, JT, NOC_CBR), intent(inout):: ODOC             ! from marine_bio_mod
      real, dimension(LT, JT, NOC_CBR), intent(inout):: ODOC13           ! from marine_bio_mod
      real, dimension(LT, JT, NOC_CBR), intent(inout):: ODOCS            ! from marine_bio_mod
      real, dimension(LT, JT, NOC_CBR), intent(inout):: ODOCS13          ! from marine_bio_mod
      real, dimension(LT, JT, NOC_CBR), intent(inout):: ONO3             ! from marine_bio_mod
      real, dimension(LT, JT, NOC_CBR), intent(inout):: OPO4             ! from marine_bio_mod      
      
      ! Use by Compute_pCO2_fluxes SUBROUTINE 
      real, dimension(LT, NOC_CBR), intent(inout):: FOC14                ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(inout):: FOC13                ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(inout):: FODIC                ! from marine_bio_mod
      real, dimension(LT, JT, NOC_CBR), intent(inout):: OC14             ! from marine_bio_mod
      real, dimension(LT, JT, NOC_CBR), intent(inout):: OALK             ! from marine_bio_mod
      real, dimension(LT, JT, NOC_CBR), intent(inout):: OC13             ! from marine_bio_mod      
      
      ! Used by Compute_oxnitrous_flux SUBROUTINE
      real, dimension(LT, JT, NOC_CBR), intent(inout), OPTIONAL:: ON2O                ! from marine_bio_mod     
      
! dmr --- These are the local variables
      INTEGER(kind=ip) :: iiso
      REAL(kind=dblp), dimension(LT,JT,NOC_CBR) :: ODIC_avt1, ODIC_avt2


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ecl   Appel des subroutines
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
       CALL Compute_oxy_flux( TM, FRICE, OO2, O2_sat_thistime, FOO2 ) 

if ( PRESENT(ON2O) ) then
       CALL Compute_oxnitrous_flux( TM, SM, FRICE, ON2O )
endif

       CALL Compute_pcO2_flux(TM, FRICE, oxpCO2, osCO2, oxCO2, oxHCO3, oxCO3,  &
                                           ODIC, FOC14, FODIC, FOC13, OC14, OC13 )
                                                 
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Application des flux de surface aux tableaux ad hoc : Subroutine CC_surface_flux
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       ODIC_avt1 = ODIC

       CALL CC_surface_flux(FOC14,    OC14,    LT,JT,NOC_CBR)
       CALL CC_surface_flux(FOPO4,    OPO4,    LT,JT,NOC_CBR)
       CALL CC_surface_flux(FONO3,    ONO3,    LT,JT,NOC_CBR)
       CALL CC_surface_flux(FOALK,    OALK,    LT,JT,NOC_CBR)
       CALL CC_surface_flux(FODIC,    ODIC,    LT,JT,NOC_CBR)
       CALL CC_surface_flux(FODOC,    ODOC,    LT,JT,NOC_CBR)
       CALL CC_surface_flux(FODOC13,  ODOC13,  LT,JT,NOC_CBR)
       CALL CC_surface_flux(FOC13,    OC13,    LT,JT,NOC_CBR)
       CALL CC_surface_flux(FODOCS,   ODOCS,   LT,JT,NOC_CBR)
       CALL CC_surface_flux(FODOCS13, ODOCS13, LT,JT,NOC_CBR)     
         
       do iiso=1,NISOO2
         CALL CC_surface_flux(FOO2(:,:,iiso), OO2(:,:,:,iiso), LT,JT,NOC_CBR)
         !WRITE(*,*) "d18O flux",(OO2(84,1,29,4)/OO2(84,1,29,2)/(2005.2E-6)-1)*1000 
         !WRITE(*,*) "d17O flux",(OO2(84,1,29,3)/OO2(84,1,29,2)/(379.9E-6_dblp)-1)*1000

         !WRITE(*,*) "D17 flux",1000000.*(LOG(1.+((OO2(84,1,29,3)/OO2(84,1,29,2)/(379.9E-6_dblp)-1)*1000)/1000.) &
         !                    -0.518*LOG(1.+((OO2(84,1,29,4)/OO2(84,1,29,2)/(2005.2E-6)-1)*1000)/1000.)) 

         !READ(*,*)
       enddo

       ODIC_avt2 = ODIC

       ! dmr --  Bricolage pour un vrai calcul du pCO2 a partir des flux ocean
       ODIC_diff = ODIC_diff + (ODIC_avt2-ODIC_avt1)

       
      END SUBROUTINE OCN_BIO_FLUXES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      FUNCTION O2_transfer_velocity(tm_cell,surf_wind) result(kg_O2)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!     DESCRIPTION : Function for gaz transfer velocity. 
!     Schmitto2 is the Schmidt number for oxygen.
!     O2_factor is the air sea transfer velocity.
!     X_conv is a constant factor to convert the piston velocity from [cm/hr] to [m/s].

!     Schmidt number is calculated according Keeling et al in 1998, Eq. 12,
!     p. 141-163 : Schmitto2 = 1638 - 81.83T + 1.483T**2 - 0.008004T**3 
!     The temperature (T) is in degrees Celsius. 

!     O2_factor (Kgo2) is calculated according Wanninkhof, 1992. 
!     Kgo2 (cm/h) = 0.39 * uav**2 (schmitto2/660)**(-0.5) 
!     uav (m/s) is the long term average wind speed. 
!     0.39 is a constant adjusted to match the large-scale mass balance
!     constraints for radiocarbon. 

        REAL(kind=dblp), INTENT(in) :: tm_cell, surf_wind
        REAL(kind=dblp), PARAMETER  :: Xconv=1._dblp/3.6e+05 
        REAL(kind=dblp)             :: schmitto2, kg_O2

 
        schmitto2=1638.0-81.83*tm_cell+1.483*tm_cell**2-0.008004*tm_cell**3 
        kg_O2=Xconv*0.337*surf_wind**2*(schmitto2/660._dblp)**(-0.5)

! TEST SIMU - Utiliser dans Nicholson, 2012 : 
!        schmitto2=1920.4-135.6*tm_cell+5.2122*tm_cell**2-0.10939*tm_cell**3+0.00093777*tm_cell**4
!        kg_O2=Xconv*0.251*surf_wind**2*(schmitto2/660._dblp)**(-0.5)

! [NOTA] drm&ec ---  Other equation version found in OOISO == 0 or OOISO == 1:
! The schimdt_O2 calculations below do not significantly change the oxygen (and
! isotope) results. This leads to a very slight decrease in oxygens. 
!~ # if ( OOISO == 0 )
!~         schmitto2=1953.4-128.0*tm_cell+3.9918*tm_cell**2 -0.050091*tm_cell**3
!~         schmitto2=1638.0-81.83*tm_cell+1.483*tm_cell**2-0.008004*tm_cell**3 
!~         kgo2(i,n)=(0.3*ws*ws+2.5*(0.5246+ttc*(0.016256+ttc*0.00049946)))*sqrt(660./schmitto2)
!~ # else 
!~         schmitto2=1920.4-135.6*tm_cell+5.2122*tm_cell**2-0.10939*tm_cell**3+0.00093777*tm_cell**4
!~         kg_O2=Xconv*0.251*surf_wind**2*(schmitto2/660._dblp)**(-0.5)
!~ # endif


      END FUNCTION O2_transfer_velocity


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       SUBROUTINE Compute_oxy_flux( TM, FRICE, OO2, O2_sat_thistime, FOO2)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Cette routine sert a appliquer les flux de surface Océan-Atmosphère appliqué pour l'oxygène lorsque 
!      le flage OOISO ==1 ou 0, WINDINCC == 0 ou 1 
!      vm --- New parametrization for O2, based on PISCES model
!
!      Auteur : J. Bendtsen
!      Date   : 15 juillet 2010
!      Derniere modification : Formation subroutine le 30/10/2024 (dmr, ecl)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       use loveclim_transfer_mod, only: MGT, SQRO2, WIND_ERA5
       use marine_bio_mod, only: JPROD 
       USE loveclim_transfer_mod, only: DVOL

#if ( OOISO == 1 )
       use iso_dioxygen_mod, only: O2_transfer_velocity_iso
       use iso_dioxygen_mod, only: iair16, nairiso
#endif       

      IMPLICIT NONE 
      
      real(kind=dblp), dimension(LT,JT,NOC_CBR), intent(in) :: TM        ! from loveclim_transfer_mod
      real(kind=dblp), dimension(LT, NOC_CBR), intent(in)   :: FRICE     ! from loveclim_transfer_mod
      real(kind=dblp), dimension(LT, JPROD, NOC_CBR), intent(in) :: O2_sat_thistime ! from O2SAT_mod
      !real(kind=dblp), dimension(LT, JT, NOC_CBR), intent(in) :: O2_sat_thistime ! from O2SAT_mod
      
      real(kind=dblp), dimension(LT, JT, NOC_CBR, NISOO2), INTENT(inout) :: OO2
      real, dimension(LT, NOC_CBR, NISOO2), intent(inout) :: FOO2              ! from marine_bio_mod


      ! Local variable
      INTEGER(kind=ip)  :: I,N
      REAL(kind=dblp)   :: O2_sat, O2_dif, norm_wind, temp_loc, O2_factor
      REAL(kind=dblp), dimension(NISOO2) :: kg_times_O2dif
      REAL(kind=dblp), parameter :: wind=4.7 ! m/s cnb

      INTEGER(kind=ip)  :: iiso   ! Loop index
      REAL(kind=dblp) :: d18O, d17O       

#define LIMIT_OCEAN_TEMP 0

    long_lp:  do n=1,NOC_CBR   ! dmr&ecl --- Begin SPATIAL loop
    lati_lp:   do i=1,LT

        if (MGT(i,1,n).eq.1) then !MGT

! dmr&ec --- Different options for wind forcing of molecular oxygen exchange       
#if ( WINDINCC == 2 )
          norm_wind = fco2ex !dmr&ec - fco2ex has a strange value of 6.E-5, not coherent with wind at 4.7 m.s-1
#elif ( WINDINCC == 1 )
         norm_wind = WS_OC(i,n)
#else
         norm_wind = wind
#endif

! ec --- If we choose the average wind at the ocean surface: we can apply a
! global wind of 4.7 m/s, or we can apply a wind varying with latitude = ERA5.
#if ( WINDS_ERA5 == 1 )
         norm_wind = WIND_ERA5(i,n,days_year360d_i)
#endif

! dmr&ec --- Different options for temp forcing of molecular oxygen exchange       
#if ( LIMIT_OCEAN_TEMP == 1 )
          temp_loc=min(35.,TM(i,1,n))
#else
          temp_loc=TM(i,1,n)
#endif

! Compute gas exchange coefficient for O2 from Waninkhof and Keeling equations 
          O2_factor = O2_transfer_velocity(temp_loc,norm_wind)
          O2_sat = O2_sat_thistime(i,1,n)
          O2_dif = OO2(i,1,n,iair) - O2_sat 
          kg_times_O2dif(iair) = (-1) * O2_factor * O2_dif

#if ( OOISO == 1 )
          if ( kg_times_O2dif(iair) == 0.0 ) then
            kg_times_O2dif(iair+1:NISOO2) = 0.0 
          else           
            kg_times_O2dif(iair+1:NISOO2) = O2_transfer_velocity_iso(temp_loc,O2_factor,O2_sat, OO2(i,1,n,:))
            !kg_times_O2dif(iair+1:NISOO2) = BRICOLE_O2dif(O2_dif, kg_times_O2dif(iair), OO2(i,1,n,:))
          endif
#endif

          DO iiso=iair,NISOO2
            FOO2(i,n,iiso)=(kg_times_O2dif(iiso)*(1-FRICE(i,n))*SQRO2(i,n))
          ENDDO  
         

         ! FOO2(i,n,iair:nairiso) = 0.0
           
          d18O = (FOO2(i,n,4)/FOO2(i,n,2)/(2005.2E-6)-1)*1000 
          !WRITE(*,*) "FLux O2 d18O",d18O
          d17O = (FOO2(i,n,3)/FOO2(i,n,2)/(379.9E-6_dblp)-1)*1000
          !WRITE(*,*) "FLux O2 d17O",d17O
          !WRITE(*,*) "FLux O2 D17",1000000.*(LOG(1.+d17O/1000.)-0.518*LOG(1.+d18O/1000.))
          !WRITE(*,*) "FLux ",FOO2(i,n,1),FOO2(i,n,2), FOO2(i,n,3), FOO2(i,n,4) 
          !READ(*,*)

        endif !MGT
         
      enddo lati_lp
      enddo long_lp
      
       END SUBROUTINE Compute_oxy_flux
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       SUBROUTINE Compute_oxnitrous_flux( TM, SM, FRICE, ON2O )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Cette routine sert a appliquer les flux de surface Océan-Atmosphère appliqué pour les oxydes nitreux 
!      Ne peut s'appliquer que lorsque le flag OXNITREUX == 1
!      Fonctionne avec WINDINCC == 0 ou 1
!      Toute la partie NITROUS OXIDE est basee sur la parametrisation du modele PISCES
!
!      Auteur : J. Bendtsen
!      Date   : 15 juillet 2010
!      Derniere modification : Formation subroutine le 30/10/2024 (dmr, ecl)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       use loveclim_transfer_mod, only: ZZ, MGT

       IMPLICIT NONE 

       ! dmr Variables from the physics: temperature, salinity, sea-ice
       REAL(kind=dblp), dimension(LT,JT,NOC_CBR), intent(in) :: TM        ! from loveclim_transfer_mod
       REAL(kind=dblp), dimension(LT,JT,NOC_CBR), intent(in) :: SM        ! from loveclim_transfer_mod
       REAL(kind=dblp), dimension(LT, NOC_CBR), intent(in)   :: FRICE     ! from loveclim_transfer_mod
       REAL, dimension(LT, JT, NOC_CBR), intent(inout):: ON2O                ! from marine_bio_mod

      ! Local variable
      INTEGER(kind=ip)                                             :: I,N
      REAL(kind = dblp)                                            :: schmittn2o, n2oatm, ws_n2o, ttc_n2o
      REAL(kind = dblp), dimension(LT,NOC_CBR)                     :: kgn2o
      REAL(kind = dblp)                                            :: zenx,znx,nx0,nx1,nx2,nx22,nx3,nx4,nx5
      REAL(kind = dblp)                                            :: fldn2o, flun2o
      REAL(kind = dblp)                                            :: tkel, sal, qtt, qtt2, zqtt

!  1) Set volumetric solubility constants for n2o in ml/l (Weiss & Price, 1980) ------------------------------------------------------
      nx0=-165.8806
      nx1=222.8743
      nx2=92.0792
      nx22=-1.48425
      nx3=-0.056235
      nx4=0.031619
      nx5=-0.0048472

!  2) compute gaz exchange coefficient (Wanninkopf(1992) equation 8 (with chemical enhancement, in cm/h) ---------

      do i=1,LT  ! ec -- Time Loop Begin 1
       do n=1,NOC_CBR
        ttc_n2o=min(39.,TM(i,1,n))
        schmittn2o=2301.1-151.1*ttc_n2o+4.7364*ttc_n2o**2-0.059431*ttc_n2o**3
        
# if ( WINDINCC == 0 )
        ws_n2o=fco2ex           !vm: variable "wind stress" prescrite
# else
        !dmr&vm real wind ...
        ws_n2o= WS_OC(i,n)
# endif

        kgn2o(i,n)=(0.3*ws_n2o*ws_n2o+2.5*(0.5246+ttc_n2o*(0.016256+ttc_n2o*0.00049946)))*sqrt(660./schmittn2o)
        kgn2o(i,n)=kgn2o(i,n)/100./3600.*(1-FRICE(i,n))      
       enddo
      enddo ! End of time loop 1
      
!  3) calculating the N2O solubility -------------------------------------------------------------------------------------------------------------

      do i=1,LT  ! Time loop begin 2
       do n=1,NOC_CBR
        tkel=TM(i,1,n)+273.15   !on met la temperature en Kelvin
        qtt=tkel*0.01
        qtt2=qtt*qtt
        sal=SM(i,1,n)+(1.-MGT(i,1,n))*35.
        zqtt=log(qtt)
        znx=nx0+nx1/qtt+nx2*zqtt+nx22*qtt2+sal*(nx3+nx4*qtt+nx5*qtt2)
        zenx=1.e-9*exp(znx)
       enddo
      enddo  ! End of time loop 2 

!  4) compute N2O flux ---------------------------------------------------------------------------------------------------------------------------

      do i=1,LT ! Begin of time loop 3
       do n=1,NOC_CBR
        flun2o=kgn2o(i,n)*ON2O(i,1,n)*MGT(i,1,n)
        fldn2o=kgn2o(i,n)*n2oatm*zenx*MGT(i,1,n)
        ON2O(i,1,n)=ON2O(i,1,n)+(fldn2o-flun2o)/zz(1)
       enddo
      enddo  ! End of time loop 3 

!  5) compute the annual N2O flux ------------------------------------------------------------------------------------------------------------
!vm : Est-ce que ce programme est lu a chaque pas de temps ?
!vm : il y en a besoin pour le calcul annuel ; equivalent LOVECLIM nspyr ?

       END SUBROUTINE Compute_oxnitrous_flux
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       SUBROUTINE Compute_pcO2_flux( TM, FRICE, oxpCO2, osCO2, oxCO2, oxHCO3, oxCO3,  &
                                                                  ODIC, FOC14, FODIC, FOC13, OC14, OC13, WS_OC )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Cette routine sert a appliquer les flux de surface Océan-Atmosphère de la pCO2 
!      Peut fonctionner avec les flags : KC14 == 0 ou 1, WINDINCC == 0 ou == 1
!
!      Auteur : J. Bendtsen
!      Date   : 15 juillet 2010
!      Derniere modification : Formation subroutine le 30/10/2024 (dmr, ecl)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       use marine_bio_mod, only: FDIC, FOAC14, oc13fair, oc13frac
       use loveclim_transfer_mod, only: ZZ, MGT, SQRO2
       use mod_sync_time, only: KENDY, TSTOC
       use C_res_mod, only: c13atm
       use carbone_co2, only: C14ATM, C14DEC, PA_C
       use carbonate_speciation_mod, only: incche

       IMPLICIT NONE 

       ! dmr Variables from the physics: temperature, salinity, sea-ice
       real(kind=dblp), dimension(LT,JT,NOC_CBR), intent(in) :: TM        ! from loveclim_transfer_mod
       real(kind=dblp), dimension(LT, NOC_CBR), intent(in)   :: FRICE     ! from loveclim_transfer_mod
      
       real, dimension(LT, NOC_CBR), intent(in):: oxpCO2 
       real, dimension(LT, NOC_CBR), intent(in):: osCO2             ! from marine_bio_mod
       real, dimension(LT, NOC_CBR), intent(in):: oxCO2             ! from marine_bio_mod
       real, dimension(LT, NOC_CBR), intent(in):: oxHCO3           ! from marine_bio_mod
       real, dimension(LT, NOC_CBR), intent(in):: oxCO3             ! from marine_bio_mod
       
       real(kind=8), dimension(LT, JT, NOC_CBR), intent(inout) ::    ODIC ! from marine_bio_mod

      !dmr --- Surface fluxes
      real, dimension(LT, NOC_CBR), intent(inout):: FOC14                ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(inout):: FOC13                ! from marine_bio_mod
      real, dimension(LT, NOC_CBR), intent(inout):: FODIC                ! from marine_bio_mod

      real, dimension(LT, JT, NOC_CBR), intent(inout):: OC14             ! from marine_bio_mod
      real, dimension(LT, JT, NOC_CBR), intent(inout):: OC13             ! from marine_bio_mod
      
      !dmr --- Optional arguments ...
      real, dimension(LT, NOC_CBR), optional, intent(in):: WS_OC             ! from marine_bio_mod

      ! dmr --- These are the local variables
      INTEGER(kind=ip) :: I ,N
      REAL(kind=dblp)  :: SCO2, FRAC1, FRAC2, FRAC3, FRACTOPEN, AFACTOR
      REAL(kind=dblp)  :: xpCO2, wind_fac

! Calculating the pCO2 at z=0 -------------------------------------------------------------------------------------------------------------------

!vm Initialisation du flux de C12 ocn -> atm (si KC14 == 1)
#if ( KC14 == 1 )
      FDIC = 0.d0
      FOAC14 = 0.d0
      FOC14 = 0.d0
#endif

      do n=1,NOC_CBR
       do i=1,LT

         if (MGT(i,1,n).eq.1) then
          ! pressure in bar=rho * g * z = Kg/m3 * N/kg * m = N/m2 = Pa
          ! 1 bar = 1e5 Pa
          ! rho = b * rho0 /g ?
          ! en attendant use rho0
          ! beware mid_level is negative...
          ! write(*,*) 'computation of p_bar', rho0, gpes,mid_level(j+1)
          ! p_bar= rho0 * gpes * (-1) * mid_level(j+1)*1e-5
          ! write(*,*) ' p_bar', p_bar

!nb     ! call incche(TM(i,j,n),SM(i,j,n),ODIC(i,j,n), 
          ! call incche(TM(i,j,n)+273.15,SM(i,j,n),p_bar,ODIC(i,j,n), OALK(i,j,n),sCO2,xpCO2,xCO2,xHCO3,xCO3,z_h)

          ! oxCO2(i,n)=xCO2

           xpCO2=oxpCO2(i,n)
           sCO2=osCO2(i,n)

           frac1=oxCO2(i,n)/ODIC(i,1,n) !(oxCO2(i,n)+oxHCO3(i,n)+oxCO3(i,n))
           frac2=oxHCO3(i,n)/ODIC(i,1,n) !(oxCO2(i,n)+oxHCO3(i,n)+oxCO3(i,n))
           frac3=oxCO3(i,n)/ODIC(i,1,n) !(oxCO2(i,n)+oxHCO3(i,n)+oxCO3(i,n))

!#if ( CORAL == 1 )
!          OCO3(i,j,n)=xCO3
!#endif

           FRACTOPEN= (1.-FRICE(i,n))
           if (FRACTOPEN.lt.0.2) FRACTOPEN=0.     !nb --- modif d'apres OCMIP

! dmr&nb pas LCM-CC         if (FRACTOPEN.lt.0.5) FRACTOPEN=0.
! dmr&nb pas LCM-CC        if (i.gt.64) FRACTOPEN=0.000
           afactor=1.
! dmr&nb pas LCM-CC        if (i.gt.12.and.i.le.24) afactor=1.5
! dmr&nb pas LCM-CC         if (i.gt.48.and.i.le.60) afactor=1.5

#if ( WINDINCC == 0 )
           wind_fac = fco2ex
#else
           !dmr&vm real wind ...
           wind_fac = WS_OC(i,n)
#endif

! ec : Value of FODIC && FOC14 (with or without wind) -- Put in a single loop rather than 2
           FODIC(i,n)=afactor*FRACTOPEN*fco2ex*sCO2*(PA_C*1.0e-6-xpCO2)* SQRO2(i,n)
           FOC14(i,n)=afactor*FRACTOPEN*fco2ex*sCO2*(1.0e-6*C14ATM-xpCO2*OC14(i,1,n)/ODIC(i,1,n))*SQRO2(i,n)

!vm : diagnostique du flux de C14 (FOAC14) et C12 (FDIC) de l'ocn vers l'atm  
#if ( KC14 == 1 )
           FOAC14 = FOAC14 - FOC14(i,n) /(SQRO2(i,n)*ZZ(1))*TSTOC
           FDIC = FDIC - FODIC(i,n) /(SQRO2(i,n)*ZZ(1))*TSTOC
#endif

! C13 fractionation from Marchal 98

           oc13fair=0.99912*(0.99869+4.9*TM(i,1,n)/1000./1000.)

           frac1=frac1*(0.99869+4.9*TM(i,1,n)*1.E-06)
           frac2=frac2*(0.99869+4.9*TM(i,1,n)*1.E-06)/(1.01078-114*TM(i,1,n)*1.E-06)
           frac3=frac3*(0.99869+4.9*TM(i,1,n)*1.E-06)/(1.00722-52*TM(i,1,n)*1.E-06)

           oc13frac=0.99912*(frac1+frac2+frac3)



! end of C13 calcul according to the wind          
#if ( WINDINCC == 0 )
         FOC13(i,n)=afactor*FRACTOPEN* &
          fco2ex*sCO2*(PA_C*1.0e-6*c13atm*oc13fair-xpCO2 & 
          *(1000+OC13(i,1,n)/ODIC(i,1,n))*oc13frac)*SQRO2(i,n)
#else
         !dmr&vm real wind ...
         FOC13(i,n)=afactor*FRACTOPEN* &
           WS_OC(i,n)*sCO2*(PA_C*1.0e-6*c13atm*oc13fair- xpCO2 & 
           *(1000+OC13(i,1,n)/ODIC(i,1,n))*oc13frac)*SQRO2(i,n)
#endif

         FOC13(i,n)=FOC13(i,n)-1000*FODIC(i,n)

! ec : Décroissance radioactive non utilisé ? 
!         IF(KENDY.EQ.1) then
!          do ja=1,JT
!           if (MGT(i,ja,n).eq.1) OC14(i,ja,n)=(1.-C14dec)*OC14(i,ja,n)
!          enddo
!         endif

         endif !ec : End of if MGT == 1
       enddo 
      enddo ! ec : End of spatial loop

! dmr&vm --- Decroissance radioactive !!
       IF(KENDY.EQ.1) then ! if end_of_year
       WHERE(MGT.eq.1)
         OC14 = (1.-C14dec) * OC14
       ENDWHERE
         WRITE(*,*) 'ocn_bio.f: ocn C14 decay'
      ENDIF

!nb --- carbonate
!     call CARBONATE

       END SUBROUTINE Compute_pcO2_flux
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       SUBROUTINE CC_surface_flux(flux,champs,ni,nj,nk)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Cette routine ade sert a appliquer les flux de surface de CYCC dans la version couplee a LOVECLIM
!      (realisee dans ADE sinon) (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche (pour le portage)
!      Date   : 15 juillet 2010
!      Derniere modification : 15 juillet 2010
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        USE mod_sync_time,         ONLY: TSTOC
        USE loveclim_transfer_mod, ONLY: SQRO2, ZZ

        IMPLICIT NONE

        INTEGER :: ni, nj, nk

        REAL, dimension(ni,nj,nk) :: champs
        REAL, dimension(ni,nk) :: flux

        champs(:,1,:) = champs(:,1,:) + flux(:,:)/(SQRO2(:,:)*ZZ(1))*TSTOC

      END SUBROUTINE CC_surface_flux
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+--


      END MODULE mod_bgc_fluxes_OCNATM
