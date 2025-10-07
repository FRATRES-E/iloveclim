!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of iLOVECLIM/ECBILT
!!      iLOVECLIM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
!!      as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!!
!!      Foobar is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!!      of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!!
!!      You should have received a copy of the GNU General Public License along with Foobar.
!!      If not, see <http://www.gnu.org/licenses/>.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"
!dmr -- Adding the choice of components through the pre-processing options

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: atmos_composition_mod
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module atmos_composition_mod is a surrogate to centralize all variables for atmospheric composition
!
!>     @date Creation date: October, 07th, 2019
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module atmos_composition_mod

       use global_constants_mod, only: dblp=>dp

       implicit none

       public :: get_PGA_CO2, set_PGA_CO2, get_PCO2_REF, atmos_carbon_update, get_PCO2_VEG

       private


! dmr PO2ref is the reference di-oxygen concentration in the atmosphere
       real(kind=dblp) :: PO2ref = 200000.0_dblp
       real(kind=dblp) :: PCO2ref = 277.147_dblp
#if ( OOISO == 1 )
       real(kind=dblp) :: PA0_O
#endif
! dmr Originally initialized as
!      PGACO2 = PCO2ref
!     In emic.f
!     Let's see how this should be replaced ...

       real(kind=dblp) :: PGACO2

! dmr patmCO2 is used to provide the varying (or not) CO2 to VECODE
! dmr Originally initialized as
!      patmCO2 = PGACO2
!     In emic.f
!     Let's see how this should be replaced ...

       real(kind=dblp) :: patmCO2


      ! NOTE_avoid_public_variables_if_possible

      contains

      function get_PGA_CO2() result(returnValue)

        real(kind=dblp) :: returnValue

        returnValue = PGACO2

        return
      end function get_PGA_CO2

      function set_PGA_CO2(new_co2_value) result(returnValue)

        real(kind=dblp), intent(in) :: new_co2_value

        logical :: returnValue

        PGACO2 = new_co2_value

        returnValue = .TRUE.

        return
      end function set_PGA_CO2

      function get_PCO2_REF() result(returnValue)

        real(kind=dblp) :: returnValue

        returnValue = PCO2ref

        return
      end function get_PCO2_REF

      function get_PCO2_VEG() result(returnValue)

        real(kind=dblp) :: returnValue

        returnValue = patmCO2

        return
      end function get_PCO2_VEG

      function set_PCO2_VEG(new_co2_value) result(returnValue)

        real(kind=dblp), intent(in) :: new_co2_value

        logical :: returnValue

        patmCO2 = new_co2_value

        returnValue = .TRUE.

        return
      end function set_PCO2_VEG


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: atmos_carbon_update
!
!>     @brief This function is updating the different pCO2s at every atmospheric timestep ...
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function atmos_carbon_update() result(returnValue)


        use comemic_mod, only: lferco2, lradco2

#if ( CYCC >= 2 )
        use C_res_mod,   only: c13atm
        use carbone_co2, only: C14ATM0, c14rstd, new_run_c, PA0_C
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  void
!>    @param[out] void
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical             :: returnValue

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!     Fertilization:
         if(lferCO2) then
            patmCO2=PGACO2
         else
            patmCO2=PCO2ref
         endif            ! lferCO2
!
!     Radiative:
         if(.NOT.lradCO2) then
             PGACO2 = PCO2ref
         endif

#if ( CYCC >= 2 )
!cnb new_run_c=1 -> fresh start, otherwise restart_mb and restart_cc are
!used
      IF (new_run_c.EQ.1) THEN
         PA0_C = PCO2ref
         c13atm = 1000.-6.5
! Init of 14C (VM version) in case of a new run
         C14ATM0 = PA0_C * c14rstd
#if ( O2ATM == 1 )
         PA0_O = PO2ref
#endif
      ENDIF
#endif /* On CYCC */

        returnValue=.TRUE.
        return
      end function atmos_carbon_update

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module atmos_composition_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
