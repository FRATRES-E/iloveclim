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
!      MODULE: comoutlocal_mod
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module comloutlocal_mod is a fortran90 module port of the original comoutlocal.h f77 include.
!
!>     @date Creation date: June, 07th, 2018
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

       module comoutlocal_mod


       use global_constants_mod, only: dblp=>dp, silp=>sp, ip

       use comatm, only: nlat, nlon, nvl

       implicit none

       private:: nlat, nlon, nvl

       integer(ip), dimension(nvl)      ::  ivlevel(nvl)
       integer(ip), dimension(0:nvl)    ::  itlevel(0:nvl)

       real(dblp), dimension(nlat,nlon) ::  tsurf1, temp4g1, temp2g1, tempsg1, albep, temp0g1, dyrain1, corain1, torain1, evap1    &
                                           ,eminp1, hesw, nlrads, runofl1, runofo1, winstu1, winstv1, snow1, swrs1, lwrs1

#if ( ISOATM >= 1 )
       real(dblp), dimension(nlat,nlon) :: iso18torain, iso18tosnow, iso17torain, iso17tosnow, isodtorain, isodtosnow
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       end module comoutlocal_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
