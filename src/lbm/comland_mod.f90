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
!      MODULE: comland_mod
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module comland_mod is a fortran90 module port of the original comland.h f77 include.
!
!>     @date Creation date: January, 29th, 2016
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

      module comland_mod

       use comatm, only: nlat, nlon, nwisos
       use global_constants_mod, only: dblp=>dp, ip, radius=>radius_earth


       implicit none

       ! dmr so far this module is a data holder for sharing, needs public variables sadly

       ! private
       private :: nlat, nlon, nwisos, dblp, ip

       public

       ! NOTE_avoid_public_variables_if_possible

       integer(kind=ip)           ::   nbasins

       real(kind=dblp), parameter :: epsl    =1.d-10
       real(kind=dblp), parameter :: pi      = 4_dblp*datan(1.d0)     &
!~                                    , radius  = 6.37e+6                &
                                   , tzero   = 273.15_dblp            &
                                   , rowat   = 1000._dblp             &
                                   , rlatvap = 2.5e+06_dblp           &
                                   , rlatsub = 2.8e+06_dblp           &
                                   , rlatfus = 0.3e+06_dblp

       real(kind=dblp), parameter :: betam = 1._dblp/(rlatfus*rowat)

       ! [TODO] There is no reason to run the landmodel only at ECBilt resolution. nlat and nlon should be changeable.
       integer, parameter          :: mbasins=26 !, nlat=32, nlon=64

       real(kind=dblp), dimension(nlat,nlon,nwisos)   :: bmoisg
       real(kind=dblp), dimension(nlat,nlon)   :: bmoism
       real(kind=dblp)                         :: lhcap, rlhcap, dtland, rdtland

       real(kind=dblp)                         :: tareas
       real(kind=dblp), save                   :: bmoismfix
       real(kind=dblp), dimension(nlat,nlon)   :: nethfxland,landheat = 0.0d0

       real(kind=dblp), dimension(nlat,nlon,nwisos)   :: dsnow, evapl

       real(kind=dblp), dimension(nlat,nlon)   :: meltheat, evapoc,tland
       real(kind=dblp)                         :: dsnm
       integer, dimension(nlat,nlon)            :: ilabas,iocbas
       integer                                  :: iscenland,islndstrt

       real(kind=dblp), dimension(nlat,nlon,nwisos)   :: runofl = 0.0d0, runofo = 0.0d0


!~        real(kind=dblp)                          :: rlatfus, rlatsub, rlatvap, rowat
       real(kind=dblp), dimension(nlat)         :: dareas, albsnow
       real(kind=dblp), dimension(mbasins)      :: arocbas, runo_yearly
       real(kind=dblp), dimension(nlat,nlon)    :: forestfr, fractl, alblbm
       real(kind=dblp), dimension(nlat,nlon,4)  :: albland

#if ( SMB_TYP >= 1 )
        real(kind=dblp), dimension(nlat,nlon,4)   :: albkeep
#endif

       real(kind=dblp), dimension(nlat,nlon,nwisos)    :: rainf, snowf

#if (ICEBERG == 2 && ISM != 2 )
       real(kind=dblp), dimension(nlat,nlon)   :: sntoicebl, sntoicebo ! arrays that contains the volume of excess snow that will be converted to icebergs over land or ocean
       real(kind=dblp)                         :: iceb_totl, iceb_toto
#endif
       real(kind=dblp)                         :: heatsnown,heatsnows

#if ( F_PALAEO_FWF == 2 )
!nb water flux from ice sheet
      real(kind=dblp), dimension(nlat,nlon)    :: wf_ice_sheet
#endif

! open files
      integer(kind=ip):: alb_dat_id, forfr_dat_id


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module comland_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
