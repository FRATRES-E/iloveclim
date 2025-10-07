!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of [CLIO-component, in following isotopic model of oceanic oxygen: iso_dioxygen_mod]
!!      Iso_dioxygen_mod is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
!!      as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!!
!!      Iso_dioxygen_mod is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!!      of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!!
!!      You should have received a copy of the GNU General Public License along with Foobar.
!!      If not, see <http://www.gnu.org/licenses/>.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [iso_dioxygen_mod]
!
!>     @author : Emeline Clermont (ecl)
!>     @author : Didier M. Roche (dmr)
!>     @author : Ji-Woong Yang (jwy)
!
!>     @brief This module [iso_dioxygen_mod] calculates oxygen isotopes in the ocean.
!
!>     @date Creation date: October, 22, 2024
!>     @date Last modification: February, 11, 2025
!>     @author Last modified by : ecl
!
!>     @version This is svn version: $LastChangedRevision$
!
!      DESCRIPTION : In this file, there are subroutines and functions that are
!      called to calculate the oxygen isotopes in the ocean.
! 
!      - Subroutine compute_ISOO2_mbiodyn: Calculation of the net O2 flux after
!      production and respiration of oxygen isotopes (called in the mbiota_mod
!      file, in the subroutine mbiodyn)
!
!      - Subroutine compute_ISOO2_maphot: Calculation of the flow of O2 which is
!      remineralised (called in the mod_aphotic_zone file (formerly maphot)  
!
!      [NOTA] These subroutines are calculated differently depending on
!      whether the Rayleigh equation is used. Please choose whether
!      Rayleigh = 1 (activate) or 0 (in choixcomposantes file).  
!
!      REVISION HISTORY:
!      2024-10-22 - Initial Version
!      TODO_YYYY-MM-DD - TODO_describe_appropriate_changes to be done or discussed - TODO_name
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       module iso_dioxygen_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       use global_constants_mod, only: dblp=>dp, ip

       implicit none

       PRIVATE

       ! jwy -- Nombre d'isotopes d'oxygene d'O2 :
       integer(kind=ip), PARAMETER, PUBLIC :: iair    = 1     ! L'O2 totale

       ! ecl -- Ototal minimum 
       real(kind=dblp), PARAMETER, PUBLIC :: Ototal_min = 30 ! Pas en dessous de 29 umolO2/kg 

#if ( OOISO == 1 )
       integer(kind=ip), PARAMETER, PUBLIC :: iair16  = 2     ! Molecule 16O-16O
       integer(kind=ip), PARAMETER, PUBLIC :: iair17  = 3     ! Molecule 17O-16O
       integer(kind=ip), PARAMETER, PUBLIC :: iair18  = 4     ! Molecule 18O-16O
       integer(kind=ip), PARAMETER, PUBLIC :: nairiso = 4     ! Number of isotopes in O2
       
       ! jwy -- d18O of ambient air O2 against SMOW (Barkan & Luz, 2005)
       ! jwy -- To convert delta_O2_smow into delta_O2_air
       real(kind=dblp), PARAMETER, PUBLIC  :: r18air = (1+23.88D-3)*(2005.2D-6)
       real(kind=dblp), PARAMETER, PUBLIC  :: r17air = (1+12.08D-3)*(379.9D-6)
       real(kind=dblp), PARAMETER          :: rair(nairiso) = [1.d0, 1.d0, r17air, r18air]

       ! ecl -- Standard Mean Ocean Water : Standard Ratios
       real(kind=dblp), PARAMETER, PUBLIC             :: r18smow = 2005.2E-6_dblp !(no unit)
       real(kind=dblp), PARAMETER, PUBLIC             :: r17smow = 379.9E-6_dblp
       real(kind=dblp), PARAMETER, dimension(nairiso) :: rsmow = [1.0_dblp,1.0_dblp,r17smow,r18smow]

       ! ecl -- GPP to NPP factor
       real(kind=dblp) :: factor_GNPP_to_NPP = 2._dblp

       ! ecl -- Parameter for fractionnement during respiration
       integer, PARAMETER :: alphaR18 = 1, alphaR17 = 2
       integer, PARAMETER :: alphaR = 2

       ! ecl -- sum of relative isotope abundances according to rsmow 
       real(kind=dblp), PARAMETER :: Rsum_iso_smow = 1.0d0 + r17smow + r18smow

       ! ecl -- Public subroutines or functions 
       PUBLIC :: compute_ISOO2_mbiodyn, compute_ISOO2_maphot
       PUBLIC :: GPPO2_func, O2_transfer_velocity_iso, OO2_saturation
       PUBLIC :: Ray_respO2, Ray_reminO2

      ! ecl - Parameters of fractionation factor during respiration
      real(kind=dblp), PARAMETER :: alpha18_respiration = 0.980_dblp !0.980_dblp !0.988_dblp
      real(kind=dblp), PARAMETER :: theta_17 = 0.518_dblp !0.518_dblp !0.520_dblp

      ! ecl - Parameters of fractionation factor during photosynthesis
      real(kind=dblp), PARAMETER :: fphoto18 = 1.0_dblp !1.004 !Not defined =1.0
      real(kind=dblp), PARAMETER :: fphoto17 = 1.0_dblp !1.002078 ! Not defined = 1.0
      contains


! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SUBROUTINE PART >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! **********************************************************************************************************************************
      SUBROUTINE compute_ISOO2_mbiodyn(NPPO2,NCP,OO2_cell,OO2_netflux,waterTemp)
! **********************************************************************************************************************************

!     DESCRIPTION : Subroutine to calculate the net isotopic oxygen flux after production and respiration for 16O, 17O and 18O.
!     Call in mbiota_mod file. 

       real(kind=dblp), intent(in)                      :: NPPO2        ! [umolO2/kg]
       real(kind=dblp), intent(in)                      :: NCP          ! [umolO2/kg]
       real(kind=dblp), DIMENSION(nairiso), intent(in)  :: OO2_cell     ! oxygen content of the ocean cell
       real(kind=dblp), DIMENSION(nairiso), intent(out) :: OO2_netflux  ! net isotopic oxgen flux
       real(kind=dblp), intent(in), optional            :: waterTemp    ! water temperature of the ocean cell

       ! Local variables
       real(kind=dblp), dimension(nairiso) :: GrossPrimProdO2, consum_oxygen
       real(kind=dblp), dimension(alphaR)  :: aresp

       if (present(waterTemp)) then 
         aresp = alpha_resp(waterTemp)
       else 
         aresp = alpha_resp()
       endif 

       GrossPrimProdO2 = GPPO2_func(NPPO2)
       consum_oxygen = consumption_O2(NCP, GrossPrimProdO2(iair), aresp, OO2_cell)
       OO2_netflux   = OO2_diff(GrossPrimProdO2, consum_oxygen)

      END SUBROUTINE compute_ISOO2_mbiodyn


! **********************************************************************************************************************************
      SUBROUTINE compute_ISOO2_maphot(OO2_dif,OO2_cell,OO2_flux_remin,waterTemp)
! **********************************************************************************************************************************

!     DESCRIPTION : Subroutine to calculate the net isotopic oxygen flux after remineralisation for 16O, 17O and 18O. 
!     Call in mod_aphotic_zone (maphot in old version) file. 

       real(kind=dblp), intent(in)                      :: OO2_dif        ! [umolO2/kg]
       real(kind=dblp), DIMENSION(nairiso), intent(in)  :: OO2_cell       ! oxygen content of the ocean cell
       real(kind=dblp), DIMENSION(nairiso), intent(out) :: OO2_flux_remin
       real(kind=dblp), intent(in),optional             :: waterTemp      ! water temperature of the ocean cell

       ! Local variables
       real(kind=dblp), dimension(alphaR) :: aresp

       if (present(waterTemp)) then 
         aresp = alpha_resp(waterTemp)
       else 
         aresp = alpha_resp()
       endif 

       OO2_flux_remin = remin_O2(OO2_dif, aresp, OO2_cell)

      END SUBROUTINE compute_ISOO2_maphot
! **********************************************************************************************************************************



! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  FUNCTIONS PART >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


! Functions required for the compute_ISOO2_mbiodyn subroutine (Without Rayleigh) ---------------------------------------------------
!===================================================================================================================================
      function GPPO2_func(NPPO2) result(GPP_O2)
!===================================================================================================================================

!     DESCRIPTION : GPPO2_func describes the gross primary production of oxygen.
!     Input variable : - NPPO2: the net primary production of oxygen.
!                      - factor: to be adjust according to oxygen production.
!                      - fphoto : fractionation during photosynthesis.
!                      - r18smow and r17smow: refer to the smow reference std for oxygen 17 and 18 conversaly.
!                      - Rsum_iso_smow : sum of relative isotope abundances according to rsmow. 
!     Output variable : GPP_O2.

       real(kind=dblp), intent(in)         :: NPPO2           ! [umolO2/kg]
       real(kind=dblp), dimension(nairiso) :: GPP_O2
       real(kind=dblp)                     :: ratio

       GPP_O2(iair) = factor_GNPP_to_NPP*NPPO2

       GPP_O2(iair17) = fphoto17 * rsmow(iair17) * GPP_O2(iair) / Rsum_iso_smow

       GPP_O2(iair18) = fphoto18 * rsmow(iair18) * GPP_O2(iair) / Rsum_iso_smow

       GPP_O2(iair16) = GPP_O2(iair) - GPP_O2(iair17) - GPP_O2(iair18)
         
      end function GPPO2_func


!===================================================================================================================================
      function alpha_resp(watTemp) result(aresp)
!===================================================================================================================================

!     DESCRIPTION : alpha_resp describes the fractionation during respiration of 18O in relation
!     to 16O, and same for the 17O. It can be calculated according to temperature or fixed.
!     Input variable : optional -> watTemp: Temperature function.
!     Output variable : - aresp(1) = alpha_18
!                       - aresp(2) = alpha_17

       real(kind=dblp), intent(in), optional :: watTemp
       real(kind=dblp), dimension(alphaR)    :: aresp
       real(kind=dblp) :: theta_O2 ! diff integration between 17O et 18O


      if (present(watTemp)) then 
       aresp(alphaR18) = (1000.0_dblp-(0.11_dblp*watTemp+12.7_dblp))/1000.0_dblp
       theta_O2 = (2.1*10**(-4))*watTemp+0.5054_dblp
 
      else 
       aresp(alphaR18) = alpha18_respiration
       theta_O2 = theta_17
      endif 
       
      aresp(alphaR17) = aresp(alphaR18)**theta_O2

      end function alpha_resp


!===================================================================================================================================
      function consumption_O2(NCP, GPPO2, arespi, OO2_cell) result(consum_O2)
!===================================================================================================================================

!     DESCRIPTION : consumption_O2 describes the net consumption of oxygen, including respiration 
!     (autotrophe and heterotrophe) and remineralisation part. 
!     Input variable : - NCP: the net oxygen balance.
!                      - GPPO2: Oxygen production.
!                      - aresp18: fractionantion between 18 and 16O during respiration.
!                      - thetaO2: fractionation between 17 and 18O during respiration.
!                      - OO2_cell: OO2 in the cell.
!     Output variable : consum_O2.

       real(kind=dblp), intent(in)                     :: NCP, GPPO2        ! [umolO2/kg]
       real(kind=dblp), dimension(alphaR), intent(in)  :: arespi
       real(kind=dblp), dimension(nairiso), intent(in) :: OO2_cell          ! [umolO2/kg]
       real(kind=dblp), dimension(nairiso)             :: consum_O2         ! [umolO2/kg]
       real(kind=dblp)                                 :: Rsum_iso, r17_cell, r18_cell

       r17_cell = OO2_cell(iair17)/OO2_cell(iair16) 
       r18_cell = OO2_cell(iair18)/OO2_cell(iair16)  
       Rsum_iso = 1.0d0 + r17_cell + r18_cell 

       consum_O2(iair) = GPPO2 - NCP
       consum_O2(iair17) = arespi(alphaR17)*r17_cell*consum_O2(iair)/Rsum_iso
       consum_O2(iair18) = arespi(alphaR18)*r18_cell*consum_O2(iair)/Rsum_iso
       consum_O2(iair16) = consum_O2(iair) - consum_O2(iair17) - consum_O2(iair18)

      end function consumption_O2


!===================================================================================================================================
      function OO2_diff(GPPO2, consum_O2) result(OO2_ddiff)
!===================================================================================================================================

!     DESCRIPTION : OO2_diff is calculated after the production and the consumption of oxygen.
!     Input variable : - GPPO2: Oxygen production.
!                      - consum_O2: Oxygen consumption.
!     Output variable : OO2_ddiff.

       real(kind=dblp), dimension(nairiso), intent(in) :: GPPO2      ! [umolO2/kg]
       real(kind=dblp), dimension(nairiso), intent(in) :: consum_O2
       real(kind=dblp), dimension(nairiso)             :: OO2_ddiff

       OO2_ddiff(iair:nairiso) = GPPO2(iair:nairiso) - consum_O2(iair:nairiso)

      end function OO2_diff


!===================================================================================================================================
      function OO2_saturation(OO2_cell,OO2_sat) result(OO2_sati)
!===================================================================================================================================

!     DESCRIPTION : OO2_saturation is calculated after the saturation of oxygen total.
!     If oxygen exceeds saturation, then it will take on the value of oxygen at saturation 
!     in the ocean. If this is the case, we need to reduce the isotopes while maintaining 
!     their proportions.
!     Input variable : - OO2_cell: Oxygen in the cell.
!                      - OO2_sat : Oxygen sat
!     Output variable : OO2_sati.

       real(kind=dblp), dimension(nairiso), intent(in)  :: OO2_cell
       real(kind=dblp), intent(in)                      :: OO2_sat
       real(kind=dblp), dimension(nairiso)              :: ratio
       real(kind=dblp), dimension(nairiso)              :: OO2_sati

      ! Calculation of isotopic ratios before adjusting
      ratio(iair16:nairiso) = OO2_cell(iair16:nairiso) / OO2_cell(iair)

      OO2_sati(iair) = min(OO2_sat, OO2_cell(iair))
      OO2_sati(iair16:nairiso) = ratio(iair16:nairiso) * OO2_sati(iair)

      end function OO2_saturation


!===================================================================================================================================
      FUNCTION Isotopic_Ratio_R18(OO2_cell) result(R18_cell)
!===================================================================================================================================

!     DESCRIPTION : R18 describes the isotopic ratio (18O/16O). 
!     Input variable : OO2_cell: OO2 in the cell.
!     Output variable : R18_cell.

       real(kind=dblp), dimension(nairiso), intent(in) :: OO2_cell
       real(kind=dblp)                                 :: R18_cell
       
       R18_cell  = OO2_cell(iair18) / OO2_cell(iair16) 

      END FUNCTION Isotopic_Ratio_R18


!===================================================================================================================================
      FUNCTION Isotopic_Ratio_R17(OO2_cell) result(R17_cell)
!===================================================================================================================================

!     DESCRIPTION : R17 describes the isotopic ratio (17O/16O). 
!     Input variable : OO2_cell: OO2 in the cell.
!     Output variable : R17_cell.

       real(kind=dblp), dimension(nairiso), intent(in) :: OO2_cell
       real(kind=dblp)                                 :: R17_cell
       
       R17_cell  = OO2_cell(iair17) / OO2_cell(iair16) 

      END FUNCTION Isotopic_Ratio_R17



! Functions required for the compute_ISOO2_maphot subroutine ----------------------------------------------------------------------
!===================================================================================================================================
      function remin_O2(OO2_dif, arespi, OO2_cell) result(remin_oxy)
!===================================================================================================================================

!     DESCRIPTION : remin_O2 describes the remineralisation part.
!     To manage Oxygen Minimum Zone (OMZ), we can introduced the Michaelis-Menton 
!     equation [Canfield, 2019], with km the half-saturation constant and 
!     ratemax the remineralisation rate. [mod_aphotic_zone]
!     Input variable : - OO2_dif : dissolved oxygen balance (calculated in maphot).
!                      - aresp18: fractionantion between 18 and 16O during respiration.
!                      - thetaO2 : fractionation between 17 and 18O during respiration.
!                      - OO2 : Oxygen Isotopes
!     Output variable : remin_oxy.

       real(kind=dblp), intent(in)                     :: OO2_dif
       real(kind=dblp), dimension(alphaR), intent(in)  :: arespi
       real(kind=dblp), dimension(nairiso), intent(in) :: OO2_cell
       real(kind=dblp), dimension(nairiso)             :: remin_oxy
       real(kind=dblp)                                 :: Rsum_iso, r17_cell, r18_cell

       r17_cell = OO2_cell(iair17)/OO2_cell(iair16) 
       r18_cell = OO2_cell(iair18)/OO2_cell(iair16)  
       Rsum_iso = 1.0d0 + r17_cell + r18_cell 

       remin_oxy(iair) = OO2_dif

       if ( remin_oxy(iair) == 0.0 ) then 
          remin_oxy = 0.0
       else 
          remin_oxy(iair17) = arespi(alphaR17)*r17_cell*remin_oxy(iair)/Rsum_iso

          remin_oxy(iair18) = arespi(alphaR18)*r18_cell*remin_oxy(iair)/Rsum_iso

          remin_oxy(iair16) = remin_oxy(iair) - remin_oxy(iair17) - remin_oxy(iair18)
       endif 

      end function remin_O2



! Rayleigh's fonctions ------------------------------------------------------------------------------------------------------------
!===================================================================================================================================
      FUNCTION Ray_respO2(GPPO2,NCP,OO2_before,OO2_cell,waterTemp) result(consum_O2)
!===================================================================================================================================

!     DESCRIPTION : resp_O2 describes the consumption of oxygen in ocean surface.
!     Here, Rayleigh fractionationation describes the evolution of a system in which one phase 
!     (respiration) is continuously removed from the system by fractional distilation. 
!     Call in mbiota_mod file. 

!     Why apply it? In the model, fractionation during respiration is applied without taking into account the 
!     quantity of oxygen present. When there is a little oxygen, marine organisms will not fractionate in the same way. 
!     [Warming] Use after production and oxygen update !    

!     Input variable : - GPPO2: the gross pirmary production of oxygen.
!                      - NCP: the net oxygen balance.
!                      - R18_cell and R17_cell: Isotopic ratios *O/16O
!                      - OO2_before: OO2 in the cell at time t-1.
!                      - OO2_cell: OO2 in the cell at time t. 
!     Output variable : consum_O2.


       real(kind=dblp), intent(in)                     :: NCP, GPPO2    
       real(kind=dblp), intent(in)                     :: OO2_before    
       real(kind=dblp), intent(in), dimension(nairiso) :: OO2_cell
       real(kind=dblp), intent(in),optional            :: waterTemp     ! water temperature of the ocean cell
       real(kind=dblp), dimension(nairiso)             :: consum_O2     

       REAL(kind=dblp), dimension(alphaR) ::  a_resp
       REAL(kind=dblp) ::  OO2_interm
       REAL(kind=dblp) ::  f_oxy, R18_interm, R17_interm
       REAL(kind=dblp) ::  R18_out, R17_out
       REAL(kind=dblp) ::  O2t, O16t, O17t, O18t

       ! Time (t-1) = OO2_before ------------------------------------------------------------

       ! Time (t-1/2) = OO2_before + GPPO2 = OO2_cell ---------------------------------------
       R18_interm = Isotopic_Ratio_R18(OO2_cell)
       R17_interm = Isotopic_Ratio_R17(OO2_cell)
       OO2_interm = OO2_before + GPPO2 

       ! Time (t) = Rayleigh -----------------------------------------------------------------
    
       ! 1) Calculation of f after t-1 : f_oxy = O2 restant / O2 init
       ! with  OO2_cell: Oxy tot in the cell after prod + resp
       consum_O2(iair) = GPPO2 - NCP
       O2t = OO2_interm - consum_O2(iair)
       f_oxy = O2t / OO2_interm

       if ( f_oxy > 1.0) then 
         WRITE(*,*) "RAYLEIGH FAIT N IMPORTE QUOI !!!"
!         WRITE(*,*) "R intermediaire", R18_interm, R17_interm
!         WRITE(*,*) "f_oxy", f_oxy
!         WRITE(*,*) "R out", R18_out, R17_out
!         WRITE(*,*) "O estimé", O16t, O17t, O18t
!         WRITE(*,*) "conso", consum_O2(iair16), consum_O2(iair17), consum_O2(iair18)
!         WRITE(*,*) "GPP,NCP", GPPO2, NCP
         f_oxy = 1.0_dblp
       endif 

       ! 2) Isotopic ratio calculated according to f_oxy. 
       if (present(waterTemp)) then 
         a_resp = alpha_resp(watTemp=waterTemp)
       else 
         a_resp = alpha_resp()
       endif 

       R18_out  = R18_interm * f_oxy**(a_resp(alphaR18)-1.0)
       R17_out = R17_interm * f_oxy**(a_resp(alphaR17)-1.0)

       ! 3) Concentration final : OO2 avec -> O2_restant = 16O+17O+18O = 16O * (1+R17 + R18)
       !O16t = OO2_cell(iair) / ( 1 + R17_out + R18_out)
       O16t = O2t / ( 1 + R17_out + R18_out)
       O18t = R18_out * O16t
       O17t = R17_out * O16t

       ! 4) Respiration flux
       consum_O2(iair17) = ABS(OO2_cell(iair17) - O17t) 
       consum_O2(iair18) = ABS(OO2_cell(iair18) - O18t) 
       consum_O2(iair16) = consum_O2(iair) - consum_O2(iair17) - consum_O2(iair18)


      END FUNCTION Ray_respO2


!===================================================================================================================================
      FUNCTION Ray_reminO2(OO2_dif, OO2_before_remin,OO2_cell,waterTemp) result(remin_oxy)
!===================================================================================================================================

!     DESCRIPTION : Subroutine to calculate the net isotopic oxygen flux after remineralisation for 16O, 17O 
!     and 18O, according to Rayleigh. 
!     Here, Rayleigh fractionationation describes the evolution of a system in which one phase 
!     (respiration) is continuously removed from the system by fractional distilation. 
!     Rayleigh fractionation describes the evolution of a system in which one phase is continuously removed from 
!     the system by fractional distillation: R = R_init*f**(alpha-1) 
!     Call in mod_aphotic_zone (maphot in old version) file. 

!     Input variable : - GPPO2: the gross pirmary production of oxygen.
!                      - NCP: the net oxygen balance.
!                      - R18_cell and R17_cell: Isotopic ratios *O/16O
!                      - OO2_before: OO2 in the cell at time t-1.
!                      - OO2_cell: OO2 in the cell at time t. 
!     Output variable : consum_O2.


       real(kind=dblp), intent(in)                     :: OO2_dif
       real(kind=dblp), intent(in)                     :: OO2_before_remin 
       real(kind=dblp), intent(in), dimension(nairiso) :: OO2_cell
       real(kind=dblp), intent(in), optional           :: waterTemp     ! water temperature of the ocean cell
       real(kind=dblp), dimension(nairiso)             :: remin_oxy     

       REAL(kind=dblp), dimension(alphaR) ::  a_resp
       REAL(kind=dblp) ::  OO2_before_resp
       REAL(kind=dblp) ::  f_remin, R18_init, R17_init
       REAL(kind=dblp) ::  R18_out, R17_out
       REAL(kind=dblp) ::  O16t, O17t, O18t

       ! Time (t-1) = OO2_before ------------------------------------------------------------
       R18_init = Isotopic_Ratio_R18(OO2_cell)
       R17_init = Isotopic_Ratio_R17(OO2_cell)

       ! Time (t) = Rayleigh on remineralisaton ---------------------------------------------
    
       ! 1) Calculation of f : f_oxy = O2 restant / O2 init 
       remin_oxy(iair) = OO2_dif
       f_remin = OO2_cell(iair) / OO2_before_remin

       ! 2) Isotopic ratio calculated according to f_oxy. 
       if (present(waterTemp)) then 
         a_resp = alpha_resp(waterTemp)
       else 
         a_resp = alpha_resp()
       endif 

       R18_out  = R18_init * f_remin**(a_resp(alphaR18)-1.0)
       R17_out = R17_init * f_remin**(a_resp(alphaR17)-1.0)

       ! 3) Final Concentration: OO2 avec -> OO2t = 16O+17O+18O = 16O * (1+R17 + R18)
       O16t = OO2_cell(iair) / ( 1 + R17_out + R18_out)
       O18t = R18_out * O16t
       O17t = R17_out * O16t

       ! 4) Respiration flux
       remin_oxy(iair17) = - (ABS(OO2_cell(iair17) - O17t)) !Keep remin negative 
       remin_oxy(iair18) = - (ABS(OO2_cell(iair18) - O18t)) 
       remin_oxy(iair16) = remin_oxy(iair) - remin_oxy(iair17) - remin_oxy(iair18)


       if ( f_remin > 1.0 ) then 
         WRITE(*,*) "RAYREMIN FAIT N IMPORTE QUOI !!!"
!         WRITE(*,*) "R init maphot", R18_init, R17_init
!         WRITE(*,*) "f_remin maphot", f_remin
!         WRITE(*,*) "R out maphot", R18_out, R17_out
!         WRITE(*,*) "O estimé maphot", O16t, O17t, O18t
!         WRITE(*,*) "remin", remin_oxy(iair16), remin_oxy(iair17), remin_oxy(iair18)
         f_remin = 1.0_dblp
       endif 

      END FUNCTION Ray_reminO2



! Function for mod_bgc_OCN file ----------------------------------------------------------------------------------------------------
!===================================================================================================================================
      FUNCTION O2_transfer_velocity_iso(tm_cell,kg_O2,O2_sat,OO2_cell) result(kg_times_O2dif)
!===================================================================================================================================

!     DESCRIPTION : Function for gaz transfer velocity * O2_diff (kg_times_O2dif). 
!     -> Formulation following equations of Wanninkhof (1992),  Najjar & al (2007), Nicholson, 2012 (Eq.7)
!     -> Avec 18alpha_gek = 0.9972 et pslp = 1.0
!     -> Avec : pslp was 1013./1013, indeed 1 

!     [NOTA]dmr&ec : This is written assuming that the flux associated can be computed as:
!                          FOO2 = kgO2*(1-FRICE(i,n))*O2_dif*SQRO2(i,n)

!     [NOTA] ec : if ( OO2 > O2_sat) then flux ocean to atm  
!                 if (OO2 < O2_sat) donc flux atm to océan

!     Input variable : - tm_cell: temperature of the cell
!                      - kg_O2: Gaz transfer velocity (function of wind and Schmidt number)
!                      - O2_sat: Saturation of oxygen in the surface layer
!                      - OO2_cell: Isotopes of oxygen in the cell 
!     Output variable : - kg_times_O2dif: kg_O2 * O2_dif

        REAL(kind=dblp),                     INTENT(in) :: tm_cell, kg_O2, O2_sat
        REAL(kind=dblp), DIMENSION(nairiso), INTENT(in) :: OO2_cell
        REAL(kind=dblp), DIMENSION(iair16:nairiso)      :: kg_times_O2dif ! dmr&ec [WARNING] dimension is (NISOO2-1)

        REAL(kind=dblp) ::  alpha_eq18, alpha_eq17
        REAL(kind=dblp) ::  kg_18, kg_17
        REAL(kind=dblp) ::  Rsum_iso, C_sat18, C_sat17


        ! Fractionnement à l'équilibre
        alpha_eq18 = 1+((-0.73+427/(tm_cell+273.15))/1000.)
        alpha_eq17 = (alpha_eq18**0.518d0)*exp((0.6d0*tm_cell+1.8d0)*1E-6)

        ! Piston velocities
        kg_18 = kg_O2 * 0.9972d0
        kg_17 = kg_O2 * (0.9972d0**0.518d0)

        ! Saturation
        Rsum_iso = 1.0d0 + r17air + r18air
        C_sat18 = r18air*O2_sat/Rsum_iso 
        C_sat17 = r17air*O2_sat/Rsum_iso

       ! FLux isotopes 
        kg_times_O2dif(iair18) = (-1) * kg_18 * (OO2_cell(iair18) - alpha_eq18 * C_sat18)
        kg_times_O2dif(iair17) = (-1) * kg_17 * (OO2_cell(iair17) - alpha_eq17 * C_sat17)
        kg_times_O2dif(iair16)= ((-1)*kg_O2*(OO2_cell(iair)-O2_sat)) - (kg_times_O2dif(iair18) + kg_times_O2dif(iair17))


! dmr&ec : Oxygen isotopes fluxes were calculated as follows, kg_times_O2diff represents part of the calculation 
!        flu18(i,n)=(-1)*0.9972*(1-FRICE(i,n))*O2_factor(i,n)*((OO2(i,j,n,4)/OO2(i,j,n,1))*OO2(i,j,n,1) &
!                              -pslp(i,n)*alpha_eq(ttc,4)*r18air*O2_sat_thistime(i,j,n))
!
!        flu17(i,n)=(-1)*(1+0.518*(0.9972-1))*(1-FRICE(i,n))*O2_factor(i,n)*((OO2(i,j,n,3)/OO2(i,j,n,1))*OO2(i,j,n,1) &
!                              -pslp(i,n)*alpha_eq(ttc,3)*r17air*O2_sat_thistime(i,j,n))
!                            
!        flu16(i,n)=(1)*(1-FRICE(i,n))*O2_factor(i,n)*(pslp(i,n)*O2_sat_thistime(i,j,n)-OO2(i,j,n,1))

      END FUNCTION O2_transfer_velocity_iso


#endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module Iso_dioxygen_mod here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      end module iso_dioxygen_mod

