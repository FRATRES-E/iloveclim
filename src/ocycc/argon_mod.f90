!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of [CLIO-component, in following isotopic model of oceanic oxygen: argon_mod]
!!      Argon_mod is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
!!      as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!!
!!      Argon_mod is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!!      of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!!
!!      You should have received a copy of the GNU General Public License along with Foobar.
!!      If not, see <http://www.gnu.org/licenses/>.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [argon_mod]
!
!>     @author : Emeline Clermont (ecl)
!
!>     @brief This module [argon_mod] calculates argon in the ocean.
!
!>     @date Creation date: February, 12, 2026
!>     @date Last modification: --
!>     @author Last modified by : --
!
!>     @version This is svn version: $LastChangedRevision$
!
!      DESCRIPTION : In this file, there are subroutines and functions that are
!      called to calculate the argon in the ocean.
! 
!      - Subroutine ? : Calculation of the ??? (called in the ?? mbiota_mod
!      file, in the subroutine mbiodyn)
!
!      REVISION HISTORY:
!      2024-10-22 - Initial Version
!      TODO_YYYY-MM-DD - TODO_describe_appropriate_changes to be done or discussed - TODO_name
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       module argon_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       use global_constants_mod, only: dblp=>dp, ip

       implicit none

       PRIVATE

#if ( ARGON == 1 )

       PUBLIC :: calcul_Argon_fluxes

      contains


! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SUBROUTINE PART >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! **********************************************************************************************************************************
      SUBROUTINE calcul_Argon_fluxes( waterTemp, waterSalt, surfWind, Ar_cell, Ar_fluxes )
! **********************************************************************************************************************************

!     DESCRIPTION : Subroutine to calculate the argon fluxes.
!     Call in mod_bgc_OCNATM file. 

       real(kind=dblp), intent(in)  :: waterTemp    ! [degreeC]
       real(kind=dblp), intent(in)  :: waterSalt    ! [Na]
       real(kind=dblp), intent(in)  :: surfWind     ! [m/s]
       real(kind=dblp), intent(in)  :: Ar_cell      ! [umol/kg]
       real(kind=dblp), intent(out) :: Ar_fluxes    ! net argon flux

       ! Local variables
       real(kind=dblp) :: Argon_Sat
       real(kind=dblp) :: kg_Ar

       Argon_Sat = Ar_saturation(waterTemp, waterSalt)
       kg_Ar     = Ar_transfer_velocity(waterTemp, surfWind)
       Ar_fluxes = (-1) * kg_Ar * (Ar_cell - Argon_Sat)


      END SUBROUTINE calcul_Argon_fluxes


! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  FUNCTIONS PART >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!===================================================================================================================================
      function Ar_saturation(tm_cell, sm_cell) result(argon_sat)
!===================================================================================================================================

!     DESCRIPTION : Argon content at equilibrium, dependent on temperature and salinity.
!     REFERENCE : Hamme and Emerson, 2004
!     Input variable : tm_cell: Temperature in the cell.
!     Output variable : argon_sat.

       real(kind=dblp), intent(in)  :: tm_cell
       real(kind=dblp), intent(in)  :: sm_cell
       real(kind=dblp)              :: argon_sat !umol/kg

       ! local variables
       real(kind=dblp)              :: temp_K, Ts                 
       real(kind=dblp)              :: log_argon_sat                
       real(kind=dblp), parameter   :: A0 = 2.79150_dblp           
       real(kind=dblp), parameter   :: A1 = 3.17609_dblp             
       real(kind=dblp), parameter   :: A2 = 4.13116_dblp                
       real(kind=dblp), parameter   :: A3 = 4.90379_dblp            
       real(kind=dblp), parameter   :: B0 = -6.96233e-3_dblp            
       real(kind=dblp), parameter   :: B1 = -7.66670e-3_dblp      
       real(kind=dblp), parameter   :: B2 = -1.16888e-2_dblp       


       ! 1) Temperature conversion  
       temp_K = tm_cell + 273.15_dblp


       ! 2) Normalisation logarithmic 
       Ts =  lOG((298.5_dblp - tm_cell) / temp_K)

       ! 3) Polynomial calculation
       log_argon_sat = A0 + A1*Ts + A2*Ts**2 + A3*Ts**3   &
                     + sm_cell*(B0 + B1*Ts + B2*Ts**2) 

       ! 4) Conversion in concentration 
       argon_sat = EXP(log_argon_sat)


      end function Ar_saturation


!===================================================================================================================================
      FUNCTION Ar_transfer_velocity(tm_cell,surf_wind) result(kg_Ar)
!===================================================================================================================================

!     DESCRIPTION : Function for gaz transfer velocity. 
!     X_conv is a constant factor to convert the piston velocity from [cm/hr] to [m/s].
!     Schmidt number is calculated according Wanninkhof, 2014 
!     Kg_Ar is calculated according Wanninkhof, 1992. 

!     Input variable : - tm_cell : Temperature in the cell.
!                      - surf_wind : Wind at the surface ocean. 
!     Output variable : kg_Ar.


        REAL(kind=dblp), INTENT(in) :: tm_cell, surf_wind
        REAL(kind=dblp)             :: kg_Ar

        ! Local variables
        REAL(kind=dblp), PARAMETER  :: Xconv=1._dblp/3.6e+05 
        REAL(kind=dblp)             :: schmitt_Ar
        REAL(kind=dblp), PARAMETER  :: piston_vel = 0.251_dblp

       
       ! 1) Schimidt number for Argon : Wanninkhof, 2014 -> A + BT + CT**2 + dT**3 + ET**4
        schmitt_Ar = 2078.1_dblp - 146.74_dblp*tm_cell                  &
                   + 5.6403_dblp*tm_cell**2 - 0.11838_dblp*tm_cell**3   &
                   + 0.0010148_dblp*tm_cell**4


       ! 2) Piston velocity for Argon : Wannninkhof, 1992
        kg_Ar = Xconv * piston_vel * surf_wind**2 * (schmitt_Ar/660._dblp)**(-0.5)


      END FUNCTION Ar_transfer_velocity


#endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module argon_mod here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      end module argon_mod

