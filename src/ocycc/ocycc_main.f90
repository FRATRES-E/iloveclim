!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of iLOVECLIM/OCYCC
!!      iLOVECLIM/OCYCC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
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
!      MODULE: ocycc_main
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module ocycc_main is handling the setup & step of the oceanic carbon and PaTh model.
!
!>     @date Creation date: October, 11th, 2019
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

      module ocycc_main

#if ( OCYCC == 1 )

       implicit none
       private

       public :: ocycc_step, ocycc_ini

      ! NOTE_avoid_public_variables_if_possible

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ocycc_ini
!
!>     @brief This function is handling a daily step in the oceanic carbon cycle and associated affairs (e.g. PaTh)
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ocycc_ini(day,month) result(returnValue)

       use MOD_SYNC_TIME, only: fait_pointer_CC_TIME, sync_timer
       use fait_pointer, only: fait_pointer_CC_OCN
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  day, month
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer, intent(in)    :: day, month

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical             :: returnValue

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       call fait_pointer_CC_OCN(0)
       call sync_lcm_ocycc(1) ! CLIO-> OCYCC to get temp and sal
       call fait_pointer_CC_TIME
       call sync_timer(day,month)
       call INIT_OCC
       call sync_lcm_ocycc(2) ! OCYCC-> CLIO to put carbon vars in scal

       returnValue = .true.

      end function ocycc_ini

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ocycc_step
!
!>     @brief This function is handling a daily step in the oceanic carbon cycle and associated affairs (e.g. PaTh)
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ocycc_step(day,month) result(returnValue)

      use mod_photic_zone, only: mbiota
      use MOD_SYNC_TIME, only: sync_timer
      use fait_pointer, only: fait_pointer_CC_OCN
      use mod_bgc_fluxes_OCNATM, only: OCN_BIO_FLUXES

#if ( PATH >= 1 )
       use path_mod, only: path_step
#endif

#if ( NEOD >= 1 )
       use neodymium_mod, only: neodymium_step
#endif

#if ( MEDUSA == 1 )
       use flux_from_sediments_mod, only: flux_from_sediments
#endif
       use DICspeciation_mod, only: DICspeciation_surface
       use O2sat_mod, only: O2sat_in_photic_zone

       ! gm: added the following USE clauses to access the dummy
       !     variables to provide to mbiota and ocn_bio below
       use loveclim_transfer_mod, ONLY: TM, SM, FRICE, MGT, DVOL, SQRO2
       use O2SAT_mod, ONLY: O2_sat_thistime
       use marine_bio_mod, only: OO2, phyto_prod, &
                                 oxpco2, osco2, oxco2, oxHCO3, oxCO3, &
                                 ODIC, FOC14, FOALK, FODOC, FODOC13, &
                                 FODOCS, FODOCS13, FONO3, FOO2, FOPO4, &
                                 FODIC, FOC13, OC14, OALK, ODOC, ODOC13, &
                                 ODOCS, ODOCS13, ONO3, OPO4, OC13, OC14
       use marine_bio_mod, only: OetaC_POMoxid, OetaN_POMoxid, OetaO2_POMoxid, &
                                 OetaC_DOMoxid_1D, OetaN_DOMoxid_1D
       use mbiota_mod, only: oc_bottom_cell
       use mbiota_mod, only: PHYTO_M, ZOO_M, PHYTO_M13, ZOO_M13, &
                             TPP_ma, caco3_ma



       implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  day, month
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer, intent(in)    :: day, month

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical             :: returnValue

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     dmr   Recopie des tableaux ad hoc pour le cycle du carbone
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

               call fait_pointer_CC_OCN(1) !TM, SM
!       print*, 'NATH sync_lcm dans step'
               CALL sync_lcm_ocycc(1) ! CLIO-> OCYCC

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     dmr   Integration du cycle du carbone ocean avant appel advection
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
               call sync_timer(day,month)

               call O2sat_in_photic_zone()
               call DICspeciation_surface()  !incche

               call mbiota(tm, frice, mgt, dvol, oc_bottom_cell, &
                      phyto_m, zoo_m, phyto_m13, zoo_m13, phyto_prod, &
                      oo2, opo4, ono3, oalk, odic, oc13, oc14, &
                      odoc, odoc13, odocs, odocs13, oxCO2, &
                      OetaC_POMoxid, OetaN_POMoxid, OetaO2_POMoxid, &
                      OetaC_DOMoxid_1D, OetaN_DOMoxid_1D, &
                      tpp_ma, caco3_ma)

#if ( OXNITREUX == 0 )
               call OCN_BIO_FLUXES(tm, sm, frice, oo2, o2_sat_thistime, &
                      oxpco2, osco2, oxco2, &
                      oxhco3, oxco3, odic, foc14, foalk, fodoc, &
                      fodoc13, fodocs, fodocs13, fono3, foo2,fopo4, fodic, foc13, &
                      oc14, oalk, odoc, odoc13, odocs, odocs13, ono3, opo4, oc13)
#else
               call OCN_BIO_FLUXES(tm, sm, frice, oo2, o2_sat_thistime, &
                      oxpco2, osco2, oxco2, &
                      oxhco3, oxco3, odic, foc14, foalk, fodoc, &
                      fodoc13, fodocs, fodocs13, fono3, foo2,fopo4, fodic, foc13, &
                      oc14, oalk, odoc, odoc13, odocs, odocs13, ono3, opo4, oc13, ONO2=ONO2)
#endif              


#if ( PATH >= 1 )
! --- lim&dmr recuperation des champs de particules
               call sync_path_ocycc(2)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      lim&dmr Appel au "step" PaTh avant advection ocean
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
               call path_step
#endif

#if ( NEOD >= 1 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      tva Appel au "step" neodymium avant advection ocean
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
               call neodymium_step
#endif

#if ( MEDUSA == 1 )
               CALL flux_from_sediments()
                                    ! From FLUX_FROM_SEDIMENT_MOD
                                    !  - steps the ODIC, OALK, etc. variables
                                    !    for one day
#endif



               CALL sync_lcm_ocycc(2) ! OCYCC-> CLIO

       returnValue = .true.

      end function ocycc_step

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#endif /* on OCYCC == 1 */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module ocycc_main

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
