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

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [ec_ca2lbm]
!
!>     @author  Thomas Extier (tex)
!>     @author  Didier M. Roche (dmr)

!>     @brief This module communicates coupler data to CARAIB model
!
!>     @date Creation date: May, 30th, 2018
!>     @date Last modification: May, 30th, 2018
!>     @author Last modified by : tex&dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

module ec_ca2lbm

#if (CARAIB > 0)

  use global_constants_mod, only : dp, ip, nbdays => days_year365d_i, ndday => days_year360d_i
  use taillesGrilles,       only : nlat => iEcb, nlon => jEcb
  use ec_co2ca,             only : shift_indexes, n_pix_caraib

  implicit none

  real(kind=dp), dimension(nlat,nlon,nbdays) :: alb_transit
#if ( CYCC > 1 )
  !real(kind=dp), dimension(nlat,nlon) :: stock_carbon_caraib
  real(kind=dp)                              :: stock_carbon_caraib
  real, dimension(n_pix_caraib)              :: pixarea
#endif
  real(kind=dp), dimension(nlat,nlon,ndday)  :: alb_for_lbm
  real(kind=dp), dimension(nlat,nlon)        :: tree_frac
  real(kind=dp), dimension(nlat,nlon)        :: grass_frac
  real(kind=dp), dimension(nlat,nlon)        :: desert_frac
  real(kind=dp), dimension(nlat,nlon)        :: veget_frac
  logical                                    :: est_calculee = .false.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


  contains

  subroutine shift_albedo

  integer(ip) :: i,r_index,timer_ca2lbm

  r_index = 1
  timer_ca2lbm = 0

  do i=1,ndday
     timer_ca2lbm = timer_ca2lbm + 1

     if ( i == shift_indexes(r_index) ) then
        timer_ca2lbm = timer_ca2lbm + 1
        r_index = r_index + 1
     endif

     alb_for_lbm(:,:,i) = alb_transit(:,:,timer_ca2lbm)

!~ 	 write(*,*) "alb_for_lbm ==", alb_for_lbm(:,:,i)
!~      write(*,*) "Time indexes == ", i, timer_ca2lbm
!~      write(*,*) "in ec_ca2lbm"

  enddo

  est_calculee = .true.

  end subroutine shift_albedo

#endif

end module ec_ca2lbm

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
