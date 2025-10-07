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
!      MODULE: comsurf_mod
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module comsurf_mod is Fortran90 replacement of comsurf.h
!
!>     @date Creation date: November, 18th, 2015
!>     @date Last modification: SLastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module comsurf_mod

       use comatm, only: nlat,nlon, nwisos
       use global_constants_mod, only: dblp=>dp, ip

       implicit none

       private :: nlat, nlon, nwisos, dblp, ip

       public

! *** Contents: Common declarations for surface dependent variables

      integer(kind=ip), parameter :: noc = 1 ,nse = 2, nld = 3,ntyps=3
      integer(kind=ip) :: iclimflux

      real(kind=dblp), dimension(nlat,nlon)       ::  tsurf
      real(kind=dblp), dimension(nlat,nlon,ntyps) ::  tsurfn

      real(kind=dblp), dimension(nlat,nlon)       :: fractoc
      real(kind=dblp), dimension(nlat,nlon,ntyps) :: fractn

      real(kind=dblp), dimension(nlat,nlon)       :: albes
      real(kind=dblp), dimension(nlat,nlon,ntyps) :: albesn

      real(kind=dblp), dimension(nlat,nlon)       :: alb2es
      real(kind=dblp), dimension(nlat,nlon,ntyps) :: alb2esn

      real(kind=dblp)                             :: albocef
      real(kind=dblp), dimension(ntyps)           :: emisn

      real(kind=dblp), dimension(nlat,nlon)       :: hesws
      real(kind=dblp), dimension(nlat,nlon,ntyps) :: heswsn

      real(kind=dblp), dimension(nlat,nlon)       :: hesw0
      real(kind=dblp), dimension(nlat,nlon,ntyps) :: hesw0n

      real(kind=dblp), dimension(nlat,nlon)       :: hesw1
      real(kind=dblp), dimension(nlat,nlon,ntyps) :: hesw1n

      real(kind=dblp), dimension(nlat,nlon)       :: hesw2
      real(kind=dblp), dimension(nlat,nlon,ntyps) :: hesw2n

      real(kind=dblp), dimension(nlat,nlon)       :: ulrads
      real(kind=dblp), dimension(nlat,nlon,ntyps) :: ulradsn

      real(kind=dblp), dimension(nlat,nlon)       :: ulrad0
      real(kind=dblp), dimension(nlat,nlon,ntyps) :: ulrad0n

      real(kind=dblp), dimension(nlat,nlon)       :: ulrad1
      real(kind=dblp), dimension(nlat,nlon,ntyps) :: ulrad1n

      real(kind=dblp), dimension(nlat,nlon)       :: ulrad2
      real(kind=dblp), dimension(nlat,nlon,ntyps) :: ulrad2n

      real(kind=dblp), dimension(nlat,nlon)       :: dlrads
      real(kind=dblp), dimension(nlat,nlon,ntyps) :: dlradsn


      real(kind=dblp), dimension(nlat,nlon,nwisos)        :: evap
      real(kind=dblp), dimension(nlat,nlon,ntyps,nwisos)  :: evapn=0.0d0

      real(kind=dblp), dimension(nlat,nlon)        :: eflux
      real(kind=dblp), dimension(nlat,nlon,ntyps)  :: efluxn = 0.0d0

      real(kind=dblp), dimension(nlat,nlon)        :: hflux
      real(kind=dblp), dimension(nlat,nlon,ntyps)  :: hfluxn

      real(kind=dblp), dimension(nlat,nlon)        :: cdragv
      real(kind=dblp), dimension(nlat,nlon,ntyps)  :: cdragvn = 0.0d0 ! [NOTA] dummy init

      real(kind=dblp), dimension(nlat,nlon)        :: q10
      real(kind=dblp), dimension(nlat,nlon,ntyps)  :: q10n

      real(kind=dblp), dimension(nlat,nlon)        :: qsurf
      real(kind=dblp), dimension(nlat,nlon,ntyps)  :: qsurfn = 0.0d0 ! [NOTA] dummy init

      real(kind=dblp), dimension(nlat,nlon)        :: tempsg
      real(kind=dblp), dimension(nlat,nlon,ntyps)  :: tempsgn

      real(kind=dblp), dimension(nlat,nlon)        :: pground
      real(kind=dblp), dimension(nlat,nlon,ntyps)  :: pgroundn

      real(kind=dblp), dimension(nlat,nlon)        :: hficof, lwrmois, winstu, winstv, abmoism

      real(kind=dblp), dimension(nlat,nlon,ntyps)  :: rmountn, evfacan

      real(kind=dblp), parameter                   :: epss = 1.E-10_8

      real(kind=dblp), dimension(nlat,nlon,nwisos)        ::  adsnow,abmoisg

      real(kind=dblp), dimension(nlat,nlon)        ::  arunofl=0.0d0,arunofo=0.0d0, ahic = 0.0d0 ! [NOTA] dummy init

! *** adding one more variable for shorwave radiation. !mohr
      real(kind=dblp), dimension(nlat,nlon)        :: swrad, lwrad

#if ( F_PALAEO_FWF == 2 )
!nb ice sheet thicknes change
      real(kind=dblp), dimension(nlat,nlon)         ::  thi_chge
#endif


end module comsurf_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
