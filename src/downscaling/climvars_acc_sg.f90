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

#if ( CARAIB > 0 || CARAIB_FORC_W > 0 )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [climvars_acc]
!
!>     @brief This subroutine calculates the climatic variables over a daily step
!
!      DESCRIPTION:
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  subroutine climvars_acc_sg(num_days,steps_inday,curr_step,shift_indexes,current_shift,max_shift)

    use global_constants_mod, only  : tK_zero_C, seconds_solar_day_i, i_days_yr => days_year365d_i

    use taillesGrilles, only       : nlat => iEcb, nlon => jEcb !utile?

    use ec_co2ca, only             : frac_land, temp_air_sg2ca, pprc_sg2ca, shr_sg2ca, rhum_sg2ca, wndspeed_sg2ca
    
    use global_constants_mod, only: ip

#if ( CARAIB_FORC_W > 0 )
    use file_libs, only           : fileDescriptor, open_f, close_f
    use global_constants_mod, only: str_len, nbdays => days_year360d_i
#endif

    use comemic_mod, only: iyear, iatm
    use comcoup_mod, only: couprf, coupsf, couptcc
    use comphys,     only: relhum, uv10
    use comsurf_mod, only: nld, tempsgn, fractn
    use comatm,      only: iwater

    use input_subgrid2L, only: lon_sg,lat_sg,nbpointssg,sub_grid_notflat,max_nb_points,torain_sg,tosnow_sg,tsurf_sg,relhum_sg    &
                             , topo_sg
    
    implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

    integer(ip), intent(in) ::  num_days,steps_inday,curr_step,current_shift,max_shift
    integer(ip), dimension(max_shift), intent(in)  :: shift_indexes

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

    integer(ip)             :: ij,i,j,n_point

    integer(ip)             :: current_shift_loc, shift_num_days_loc

    integer(ip), parameter  :: j_antarctique = 7
    
#if ( CARAIB_FORC_W > 0 )

    type(fileDescriptor)                :: tempair_f,pprc_f,shr_f,rhum_f,wndspeed_f,dhu_f,dom_f


! dmr --- File name variable (nm)
    character(len=str_len)              ::                                               &
                     tempair_nm   = "outputdata/coupler/temsg",                          &
                     pprc_nm      = "outputdata/coupler/prcsg",                          &
                     shr_nm       = "outputdata/coupler/shrsg",                          &
                     rhum_nm      = "outputdata/coupler/rhusg",                          &
                     wndspeed_nm  = "outputdata/coupler/wndsg",                          &
                     dhu_nm       = "outputdata/coupler/dhusg",                          &
                     dom_nm       = "outputdata/coupler/domsg"

! dmr --- Formatted flag: .false. is unformatted file
    logical, parameter                  ::                                              &
                     tempair_fm   = .true. ,                                            &
                     pprc_fm      = .true. ,                                            &
                     shr_fm       = .true. ,                                            &
                     rhum_fm      = .true. ,                                            &
                     wndspeed_fm  = .true.,                                             &
                     dhu_fm       = .true.,                                             &
                     dom_fm       = .true.

    character(len=10)                   :: cyear
    integer(ip)                         :: nb_hours = 0

#endif

    shift_num_days_loc = num_days + current_shift

    ij = 0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! afq -- We do not have restart info for the downscaling, and the first climate fields computed here are from the restart
!        Here we use the coarse grid info for the first year
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
    if ( iyear .gt. 1 ) then
       do j=1,nlat
          do i=1,nlon
             !if(fractn(j,i,nld) .GE. frac_land .AND. j .GE. j_antarctique .AND. sub_grid_notflat(j,i).GT.0 ) then
             if(fractn(j,i,nld) .GE. 0d0 .AND. j .GE. j_antarctique .AND. sub_grid_notflat(j,i).GT.0 ) then
                do n_point=1,nbpointssg(j,i)
                   if (topo_sg(j,i,n_point) .GT. 0) then
                     ij = ij + 1
                     temp_air_sg2ca(shift_num_days_loc,ij) = tsurf_sg(j,i,n_point)-tK_zero_C
                     pprc_sg2ca(shift_num_days_loc,ij) = (torain_sg(j,i,n_point)+tosnow_sg(j,i,n_point))*1000.*seconds_solar_day_i
                     shr_sg2ca(shift_num_days_loc,ij) = (1-couptcc(j,i))*100.
                     rhum_sg2ca(shift_num_days_loc,ij) = relhum_sg(j,i,n_point)*100.
                     wndspeed_sg2ca(shift_num_days_loc,ij) = uv10(j,i)
                   endif
                enddo
             end if
          end do
       end do
    else
       do j=1,nlat
          do i=1,nlon
             !if(fractn(j,i,nld) .GE. frac_land .AND. j .GE. j_antarctique .AND. sub_grid_notflat(j,i).GT.0 ) then
             if(fractn(j,i,nld) .GE. 0d0 .AND. j .GE. j_antarctique .AND. sub_grid_notflat(j,i).GT.0 ) then
                do n_point=1,nbpointssg(j,i)
                   if (topo_sg(j,i,n_point) .GT. 0) then
                     ij = ij + 1
                     temp_air_sg2ca(shift_num_days_loc,ij) = tempsgn(j,i,nld)-tK_zero_C
                     pprc_sg2ca(shift_num_days_loc,ij) = (couprf(j,i,iwater)+coupsf(j,i,iwater))*1000.*seconds_solar_day_i
                     shr_sg2ca(shift_num_days_loc,ij) =  ((1-couptcc(j,i))/steps_inday)*100.
                     rhum_sg2ca(shift_num_days_loc,ij) = (relhum(j,i)*100.)/steps_inday
                     wndspeed_sg2ca(shift_num_days_loc,ij) = uv10(j,i)
                   endif
                enddo
             end if
             !write(*,*) "nbpix", ij !use this to get the number of pixels
          end do
       end do
    endif
   
    if ( curr_step  == steps_inday ) then
       if (ANY( shift_indexes == shift_num_days_loc) ) then

          current_shift_loc = current_shift + 1
          shift_num_days_loc = num_days + current_shift_loc

          temp_air_sg2ca(shift_num_days_loc,:) = temp_air_sg2ca(shift_num_days_loc-1,:)
          pprc_sg2ca(shift_num_days_loc,:) = pprc_sg2ca(shift_num_days_loc-1,:) 
          shr_sg2ca(shift_num_days_loc,:) = shr_sg2ca(shift_num_days_loc-1,:)
          rhum_sg2ca(shift_num_days_loc,:) = rhum_sg2ca(shift_num_days_loc-1,:)
          wndspeed_sg2ca(shift_num_days_loc,:) = wndspeed_sg2ca(shift_num_days_loc-1,:)

          write(*,*) "NUM_day in climvars_acc_sg ##", shift_num_days_loc, curr_step, current_shift_loc

       endif
    endif

#if ( CARAIB_FORC_W > 0 )

    nb_hours = nb_hours + 1
    if(nb_hours.ge.iatm) then
       nb_hours = 0_ip
    end if

    if (mod(num_days,nbdays).eq.0) then

       if(nb_hours.eq.0) then

          write(cyear,'(i10)'),iyear
          print*, "cyear", cyear

          tempair_f%isFormatted = tempair_fm
          call open_f(tempair_f, trim(tempair_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
          pprc_f%isFormatted = pprc_fm
          call open_f(pprc_f, trim(pprc_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
          shr_f%isFormatted = shr_fm
          call open_f(shr_f, trim(shr_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
          rhum_f%isFormatted = rhum_fm
          call open_f(rhum_f, trim(rhum_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
          wndspeed_f%isFormatted = wndspeed_fm
          call open_f(wndspeed_f, trim(wndspeed_nm)//trim(adjustl(cyear))//".dat", o_stat="new")

          ij = 0

          do j=1,nlat
             do i=1,nlon
                !if(fractn(j,i,nld) .GE. frac_land .AND. j .GE. j_antarctique .AND. sub_grid_notflat(j,i).GT.0 ) then
                if(fractn(j,i,nld) .GE. 0d0 .AND. j .GE. j_antarctique .AND. sub_grid_notflat(j,i).GT.0 ) then

                   do n_point=1,nbpointssg(j,i)

                      if (topo_sg(j,i,n_point) .GT. 0) then

                        ij = ij + 1

                        write(tempair_f%id,'(370f10.4)') lon_sg(j,i,n_point),lat_sg(j,i,n_point),temp_air_sg2ca(1:i_days_yr,ij)
                        write(pprc_f%id,'(370f10.4)') lon_sg(j,i,n_point),lat_sg(j,i,n_point),pprc_sg2ca(1:i_days_yr,ij)
                        write(shr_f%id,'(370f10.4)') lon_sg(j,i,n_point),lat_sg(j,i,n_point),shr_sg2ca(1:i_days_yr,ij)
                        write(rhum_f%id,'(370f10.4)') lon_sg(j,i,n_point),lat_sg(j,i,n_point),rhum_sg2ca(1:i_days_yr,ij)
                        write(wndspeed_f%id,'(370f10.4)') lon_sg(j,i,n_point),lat_sg(j,i,n_point),wndspeed_sg2ca(1:i_days_yr,ij)

                      endif

                   enddo

                end if
             end do
          end do

          call close_f(tempair_f)
          call close_f(pprc_f)
          call close_f(shr_f)
          call close_f(rhum_f)
          call close_f(wndspeed_f)

       endif
    endif

#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  end subroutine climvars_acc_sg


#endif
