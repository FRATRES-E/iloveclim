!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of [insert sub-component name here, in following Foobar]
!!      Foobar is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
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
!      MODULE: comrunlabel_mod
!
!>     @author  Pepijn Bakker (PB)
!
!>     @brief This module newunit_mod is a module to initialize the variable that is used to read files.
!
!>     @date Creation date: May, 07th, 2021
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : PB
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       module newunit_mod


       use global_constants_mod, only: ip

       implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** Contents: real start year of the run
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** PARAMETERS
! *** COMMON  /ilabel/  irunlabelf
!     irunlabelf: real date at the starting point of this run
!                 (irunlabelf=irunlabel-initdate)
      integer(ip) :: newunit_id, wisocpl_restart_id, error_id, ocbasin_id, ocheattr_id, sum_dat_id, win_dat_id, berg_dat_id    &
                   , coef_dat_id, lwrref_dat_id, lwrcoef_dat_id, swrref_dat_id, swrcoef_dat_id, GHG_dat_id, TSI_RM_dat_id      &
                   ,VOLC_dat_id, SUL_dat_id, OZONE_dat_id, scenario2Xco2_dat_id, info_id, mbcs2_cor_id, book_id, ipcc_id       &
                   , runoff_id, namelistecbilt_id, gauss_asc_id, parameterschk_id, carbon_emission_dat_id

!~       common /ilabel/  irunlabelf
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       end module newunit_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
