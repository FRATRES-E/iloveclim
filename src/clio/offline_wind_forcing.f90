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
#include 'choixcomposantes.h'
!dmr -- Adding the choice of components through the pre-processing options

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: offline_wind_forcing
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module, offline_wind_forcing, allows to replace the winds in the coupler by fields read in external files 
!
!>     @date Creation date: October, 27th, 2022
!>     @date Last modification: SLastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!      DESCRIPTION
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!>      Compute formula: \f$ \frac{d\lambda}{dt} , \frac{d\phi}{dt},  \frac{dz}{dt} \f$
!!
!!      REFERENCES: papers or other documents to be cited...(including link when possible)    
!
!       REVISION HISTORY:
!       2022_10_27 - Initial Version
!       TODO_YYYY-MM-DD - TODO_describe_appropriate_changes to be done or discussed - TODO_name
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module offline_wind_forcing

       use global_constants_mod, only: str_len, dblp=>dp, months_year_i, days_year360d_i
       use para0_mod,            only: imax, jmax

       implicit none
       private
       public :: update_winds_off, get_winds

         logical                                               :: is_initialized = .false.
         character(len=str_len), parameter                     :: filename_NC_winds="inputdata/forcing/monthly_climatology-CLIOT.nc"
         real(kind=dblp), dimension(imax,jmax,months_year_i)   :: monthly_winds_CLIO_u, monthly_winds_CLIO_v
         real(kind=dblp), dimension(imax,jmax,days_year360d_i) :: daily_winds_CLIO_u, daily_winds_CLIO_v

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: read_winds_from_file
!
!>     @brief This subroutine use the given path to read winds in an external netCDF
!
!      DESCRIPTION:
!
!       REVISION HISTORY:
!       2022-10-28 - Initial Version
!       TODO_YYYY-MM-DD - TODO_describe_appropriate_changes to be done or discussed - TODO_name
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function read_winds_from_file() result(returnValue)

       use ncio,      only: nc_read

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables [NONE]
!>    @return returnValue
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical             :: returnValue   !> @var logical to test success/failure of the function
       logical             :: f_exists

       integer, parameter  :: imaxm2 = imax-2, imaxm1 = imax-1
       real(kind=dblp), dimension(imaxm2,jmax,months_year_i) :: monthly_winds_to_read

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         inquire(file=filename_NC_winds,exist=f_exists)

         if (f_exists) then
               ! proceed with opening and reading

               call nc_read(filename_NC_winds,"10u",monthly_winds_to_read)

               monthly_winds_CLIO_u(2:imaxm1,:,:) = monthly_winds_to_read(:,:,:)
               monthly_winds_CLIO_u(1,:,:) = monthly_winds_to_read(imaxm2,:,:)
               monthly_winds_CLIO_u(imax,:,:) = monthly_winds_to_read(1,:,:)

               call nc_read(filename_NC_winds,"10v",monthly_winds_to_read)

               monthly_winds_CLIO_v(2:imaxm1,:,:) = monthly_winds_to_read(:,:,:)
               monthly_winds_CLIO_v(1,:,:) = monthly_winds_to_read(imaxm2,:,:)
               monthly_winds_CLIO_v(imax,:,:) = monthly_winds_to_read(1,:,:)

         else
               WRITE(*,*) "File "//TRIM(filename_NC_winds)//" does not exist! [STOP]"
               stop
         endif


      end function read_winds_from_file

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function read_winds_from_file
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: update_winds_off
!
!>     @brief This subroutine use the given path to read winds in an external netCDF
!
!      DESCRIPTION:
!
!       REVISION HISTORY:
!       2022-10-28 - Initial Version
!       TODO_YYYY-MM-DD - TODO_describe_appropriate_changes to be done or discussed - TODO_name
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function update_winds_off() result(returnValue)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables [NONE]
!>    @return returnValue
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical             :: returnValue   !> @var logical to test success/failure of the function

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         if (.not. is_initialized) then
               returnValue = read_winds_from_file() ! I am assuming variables are read with monthly climatological mean for now ...
               returnValue = interpolate_monthly_to_daily(monthly_winds_CLIO_u, daily_winds_CLIO_u)
               returnValue = interpolate_monthly_to_daily(monthly_winds_CLIO_v, daily_winds_CLIO_v)
               is_initialized = .true.
         endif


      end function update_winds_off

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function read_winds_from_file
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      function get_winds(direction,iday) result(twoD_winds)

              use global_constants_mod, only: sp

              integer(sp), intent(in) :: direction, iday

              real(kind=dblp), dimension(imax,jmax) :: twoD_winds

              if (direction .eq. 0 ) then
                      twoD_winds(:,:) = daily_winds_CLIO_u(:,:,iday)
              else
                      twoD_winds(:,:) = daily_winds_CLIO_v(:,:,iday)
              endif
      end function get_winds




!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: interpolate_monthly_to_daily
!
!>     @brief This subroutine interpolate monthly data to daily data using a spline
!
!      DESCRIPTION:
!
!       REVISION HISTORY:
!       2022-10-28 - Initial Version
!       TODO_YYYY-MM-DD - TODO_describe_appropriate_changes to be done or discussed - TODO_name
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function interpolate_monthly_to_daily(array_to_interpol, array_interpolatd) result(returnValue)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables [NONE]
!>    @return returnValue
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
         use global_constants_mod, only: days_month360d_i, sp
         use mod_splines, only: spline

         real(kind=dblp), dimension(:,:,:), intent(in)  :: array_to_interpol
         real(kind=dblp), dimension(:,:,:), intent(out) :: array_interpolatd

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical             :: returnValue   !> @var logical to test success/failure of the function
       real(kind=dblp), dimension(:), allocatable    :: index_daily, index_monthly
       integer(kind=sp)    :: i,j, dimxl, dimxr, dimyl, dimyr
       real(kind=dblp), dimension(:,:,:), allocatable :: local_array_copy
       type(spline) :: local_spline

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       dimxl = LBOUND(array_to_interpol,dim=1)
       dimxr = UBOUND(array_to_interpol,dim=1)
       dimyl = LBOUND(array_to_interpol,dim=2)
       dimyr = UBOUND(array_to_interpol,dim=2)

       allocate(local_array_copy(dimxl:dimxr,dimyl:dimyr,0:months_year_i))

       local_array_copy(:,:,1:months_year_i) = array_to_interpol(:,:,:)
       local_array_copy(:,:,0) = array_to_interpol(:,:,months_year_i)

       index_monthly = [(real((i-1)*days_month360d_i+days_month360d_i/2,kind=dblp),i=0,months_year_i,1)]
       index_daily   = [(real(j,kind=dblp),j=1,days_year360d_i,1)]

       do i=dimxl,dimxr
       do j=dimyl,dimyr
         local_spline = spline(index_monthly(:),local_array_copy(i,j,:)) 
         array_interpolatd(i,j,:) = local_spline%value(index_daily)
       enddo
       enddo

!       do i=2,UBOUND(index_monthly,dim=1)
!          write(*,*) "before ::", index_monthly(i), array_to_interpol(12,12,i-1)
!       enddo
!       do i=1,UBOUND(index_daily,dim=1)
!          write(*,*) "after ::", index_daily(i), array_interpolatd(12,12,i)
!       enddo
      end function interpolate_monthly_to_daily

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function read_winds_from_file
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module Foo_mod here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module offline_wind_forcing

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
