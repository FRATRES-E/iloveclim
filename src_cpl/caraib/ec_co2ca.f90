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
!      MODULE: [ec_co2ca]
!
!>     @author  Thomas Extier (tex)
!>     @author  Didier M. Roche (dmr)

!>     @brief This module communicates coupler data to CARAIB model
!
!>     @date Creation date: October, 16th, 2017
!>     @date Last modification: November, 27th, 2019
!>     @author Last modified by : tex&dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

module ec_co2ca



  use global_constants_mod, only : dp, sp, ip, nmth => months_year_i, nbdays => days_year360d_i, i_days_yr => days_year365d_i     &
                                 , str_len
  use taillesGrilles, only       : nlat => iEcb, nlon => jEcb

  implicit none

  public :: reset_climvars
  public :: climvars_acc
  public :: orbitals_for_CARAIB
  public :: co2_for_CARAIB
#if (CARAIB > 0)
  public :: sync_coupler_caraib
#endif

#if ( CARAIB > 0 || CARAIB_FORC_W > 0 )
  integer(ip), parameter                      :: n_pix_caraib = 541          !622 for LGM
#if ( DOWNSCALING == 2 )
  integer(ip), parameter                      :: n_pix_caraib_nh40 = 27691
  integer(ip), parameter                      :: n_pix_caraib_africa = 39999 !todo
  integer(ip), parameter                      :: n_pix_caraib_europe = 22060 !24224
#endif
#endif

  integer(ip), parameter                      :: j_antarctique = 7
#if ( HOURLY_RAD > 0 )
  integer(ip), parameter                      :: nb_hour_steps_perday = 6
#endif
  integer(ip), parameter                      :: npft = 26
  real(dp),    parameter                      :: frac_land = 0.3_dp


! --- dmr&tex Listes des tableaux de transmission coupleur => CARAIB

  real(sp), dimension(i_days_yr,n_pix_caraib)              :: temp_air     ! daily air temperature at 10 meter high for each surface type in °C
  real(sp), dimension(i_days_yr,n_pix_caraib)              :: pprc         ! daily average precipitation in mm/mth
  real(sp), dimension(i_days_yr,n_pix_caraib)              :: shr          ! daily percentage of sunshine hours
  real(sp), dimension(i_days_yr,n_pix_caraib)              :: rhum         ! daily air relative humidity in %
  real(sp), dimension(i_days_yr,n_pix_caraib)              :: wndspeed     ! daily values of wind speed at the surface in m/s
  real(sp), dimension(i_days_yr,n_pix_caraib)              :: dom          ! groundwater d18O isotopes
  real(sp), dimension(i_days_yr,n_pix_caraib)              :: dom17        ! groundwater d17O isotopes
  real(sp), dimension(i_days_yr,n_pix_caraib)              :: dom_D        ! groundwater dD isotopes
  real(sp), dimension(i_days_yr,n_pix_caraib)              :: dhu          ! water vapor d18O isotopes
  real(sp), dimension(i_days_yr,n_pix_caraib)              :: dhu17O       ! water vapor d17O isotopes
  real(sp), dimension(i_days_yr,n_pix_caraib)              :: dhu_D        ! water vapor dD isotopes
  real(sp), dimension(i_days_yr,n_pix_caraib)              :: dp_O         ! precipitation dO isotopes
  real(sp), dimension(i_days_yr,n_pix_caraib)              :: dp_D         ! precipitation dD isotopes

#if ( HOURLY_RAD > 0 )
  real(sp), dimension(i_days_yr,n_pix_caraib)                       :: atm_transmissivity ! atmosphere transmissivity, TOA to surface
  real(sp), dimension(i_days_yr,nb_hour_steps_perday,n_pix_caraib)  :: solar_flux_hourly  ! hourly TOA solar flux in W.m-2
#endif

  real(sp), dimension(nlat,nlon)              :: d_vap
  real(sp), dimension(nlat,nlon)              :: d_vap_17O
  real(sp), dimension(nlat,nlon)              :: d_vap_D
  real(sp), dimension(nlat,nlon)              :: d_gw
  real(sp), dimension(nlat,nlon)              :: d_gw17
  real(sp), dimension(nlat,nlon)              :: d_gw_D
  real(sp), dimension(nlat,nlon)              :: d_precip_D
  real(sp), dimension(nlat,nlon)              :: d_precip

  real(sp)                                    :: stdCO2

  real(sp)                                    :: orb_ecc
  real(sp)                                    :: orb_obl ! for CARAIB, needs to be in rad.
  real(sp)                                    :: orb_per ! for CARAIB, needs to be in rad.

!~   integer(ip), parameter                      :: nbdays_mnth = nbdays / nmth     ! number of days per month

  integer(ip), parameter                      :: max_shift = 5
  integer(ip), dimension(max_shift)           :: shift_indexes
  integer(ip)                                 :: current_shift, shift_num_days

#if ( CARAIB > 0 || CARAIB_FORC_W > 0 )
#if ( DOWNSCALING == 2 )
    real(sp), dimension(i_days_yr,n_pix_caraib_europe) :: temp_air_sg2ca, pprc_sg2ca, shr_sg2ca, rhum_sg2ca, wndspeed_sg2ca
#endif
#endif

  contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [climvars_acc]
!
!>     @brief This subroutine calculates the climatic variables over a daily step
!
!      DESCRIPTION:
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  subroutine climvars_acc(num_days,steps_inday,curr_step)

    use global_constants_mod, only  : tK_zero_C,seconds_solar_day_i
    use comatm, only                : nsh2, nvl, ntl, nsh, nm

#if ( CARAIB_FORC_W > 0 )
    use file_libs, only             : fileDescriptor, open_f, close_f
#endif

#if ( ISOATM >= 1 )
    use iso_param_mod, only         : ieau, ieau18, ieau17, ieaud, rsmow
#endif

#if ( COMATM == 1 )
    use comemic_mod, only: iyear, iatm
    use comcoup_mod, only: couprf, coupsf, couptcc
    use comphys,     only: relhum, uv10, rmoisg, fswdsfc, torain, tosnow
    use comsurf_mod, only: nld, tempsgn, fractn
    use comland_mod, only: bmoisg, bmoismfix

#if ( HOURLY_RAD > 0 )
    use comphys,              only: solflux_hourly, fswdsfc
    use commons_mod,          only: fswdtoa0
    use global_constants_mod, only: eps_sp
#endif

#endif

    implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! tex&dmr --- num_days is the day in the year
    integer(ip), intent(in) :: num_days,steps_inday,curr_step

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

    integer(ip)             :: ij,i,j

#if ( CARAIB_FORC_W > 0 )

    type(fileDescriptor)    :: tempair_f,pprc_f,shr_f,rhum_f,wndspeed_f,dhu_f,dhu17O_f,dhuD_f,dom_f,dom17_f,domD_f,dpD_f,dp_O_f


! dmr --- File name variable (nm)
    character(len=str_len)              ::                                               &
                     tempair_nm   = "outputdata/coupler/tem",                            &
                     pprc_nm      = "outputdata/coupler/prc",                            &
                     shr_nm       = "outputdata/coupler/shr",                            &
                     rhum_nm      = "outputdata/coupler/rhu",                            &
                     wndspeed_nm  = "outputdata/coupler/wnd",                            &
                     dhu_nm       = "outputdata/coupler/dhu",                            &
                     dhu17O_nm    = "outputdata/coupler/dhu17O",                         &
                     dhuD_nm      = "outputdata/coupler/dhuD",                           &
                     dom_nm       = "outputdata/coupler/dom",                            &
                     dom17_nm      = "outputdata/coupler/dom17",                         &
                     domD_nm      = "outputdata/coupler/domD",                           &
                     dpD_nm       = "outputdata/coupler/dpD",                            &
                     dp_O_nm      = "outputdata/coupler/dpO_"

! dmr --- Formatted flag: .false. is unformatted file
    logical, parameter                  ::                                               &
                     tempair_fm   = .true.,                                              &
                     pprc_fm      = .true.,                                              &
                     shr_fm       = .true.,                                              &
                     rhum_fm      = .true.,                                              &
                     wndspeed_fm  = .true.,                                              &
                     dhu_fm       = .true.,                                              &
                     dhu17O_fm    = .true.,                                              &
                     dhuD_fm      = .true.,                                              &
                     dom_fm       = .true.,                                              &
                     dom17_fm     = .true.,                                              &
                     domD_fm      = .true.,                                              &
                     dpD_fm       = .true.,                                              &
                     dp_O_fm      = .true.

    real(dp), allocatable, dimension(:) :: templon,templat                              ! Value of the latitude and longitude at 5.625°

    real(dp)                            :: phi(nlat),rj,rw,dumwei,ari(nlat/2),land_frac,lon1,lat1,oro_g
    real(dp), parameter                 :: pi = 2*acos(0.)
    real(dp), parameter                 :: radian_to_degree = 180_dp/pi
    integer(ip)                         :: k,l,ilat,ird
    character(len=10)                   :: cyear
    integer(ip)                         :: nb_hours = 0, gauss_asc_id

#endif

    shift_num_days = num_days + current_shift

    ij = 0
    bmoismfix=0.15

#if (ISOATM >= 2)

    do j=1,nlat
      do i=1,nlon
         if(fractn(j,i,nld) .GE. frac_land .AND. j .GE. j_antarctique) then
            if (rmoisg(j,i,ieau).gt.0.0) then                                                         ! Water vapor
              d_vap(j,i) = (rmoisg(j,i,ieau18) / rmoisg(j,i,ieau) / rsmow(ieau18) - 1.0d0) * 1000
              d_vap_17O(j,i) = (rmoisg(j,i,ieau17) / rmoisg(j,i,ieau) / rsmow(ieau17) - 1.0d0) * 1000
              d_vap_D(j,i) = (rmoisg(j,i,ieaud) / rmoisg(j,i,ieau) / rsmow(ieaud) - 1.0d0) * 1000
            else
              d_vap(j,i) = 0.0
              d_vap_17O(j,i) = 0.0
              d_vap_D(j,i) = 0.0
            endif
            if (bmoisg(j,i,ieau).gt.0.0) then                                                         ! Groundwater
              d_gw(j,i) = (bmoisg(j,i,ieau18) / bmoisg(j,i,ieau) / rsmow(ieau18) - 1.0d0) *1000
              d_gw17(j,i) = (bmoisg(j,i,ieau17) / bmoisg(j,i,ieau) / rsmow(ieau17) - 1.0d0) *1000
              d_gw_D(j,i) = (bmoisg(j,i,ieaud) / bmoisg(j,i,ieau) / rsmow(ieaud) - 1.0d0) *1000
            else
              d_gw(j,i) = 0.0
              d_gw17(j,i) = 0.0
              d_gw_D(j,i) = 0.0
            endif
            if (torain(j,i,ieau).ge.0.0 .AND. tosnow(j,i,ieau).ge.0.0) then                           ! Precipitation
              d_precip(j,i) = ((torain(j,i,ieau18)+tosnow(j,i,ieau18)) / (torain(j,i,ieau)+tosnow(j,i,ieau)) / rsmow(ieau18) - 1.0d0) *1000
              d_precip_D(j,i) = ((torain(j,i,ieaud)+tosnow(j,i,ieaud)) / (torain(j,i,ieau)+tosnow(j,i,ieau)) / rsmow(ieaud) - 1.0d0) *1000
            else
              d_precip(j,i) = 0.0
              d_precip_D(j,i) = 0.0
            endif
         endif
      enddo
    enddo

#endif

#if ( DOWNSCALING == 2 )
    call climvars_acc_sg(num_days,steps_inday,curr_step,shift_indexes,current_shift,max_shift) ! have to be before the standard block below
#endif

    do j=1,nlat
      do i=1,nlon
         if(fractn(j,i,nld) .GE. frac_land .AND. j .GE. j_antarctique) then

            ij = ij + 1
            temp_air(shift_num_days,ij) = tempsgn(j,i,nld)-tK_zero_C
            rhum(shift_num_days,ij) = relhum(j,i)*100
            wndspeed(shift_num_days,ij) = uv10(j,i)

#if ( ISOATM >= 2 )
            pprc(shift_num_days,ij) = (couprf(j,i,ieau)+coupsf(j,i,ieau))*1000*seconds_solar_day_i
#else
            pprc(shift_num_days,ij) = (couprf(j,i,1)+coupsf(j,i,1))*1000*seconds_solar_day_i
#endif

#if ( HOURLY_RAD == 0 )
! --- dmr  COUPLING through cloud cover ...
            shr(shift_num_days,ij) = ((1-couptcc(j,i)))*100
#endif

#if ( HOURLY_RAD > 0 ) /* Coupling through solar flux and transmissivity */
            solar_flux_hourly(shift_num_days,:,ij) = solflux_hourly(:,j)
            if ( fswdtoa0(j,i).GT.eps_sp ) then
              atm_transmissivity(shift_num_days,ij) = fswdsfc(j,i)/fswdtoa0(j,i)
            else
              atm_transmissivity(shift_num_days,ij) = 0.0_dp
            endif
#endif

!~            if ((j.eq.28).and.(i.eq.12)) then
!~              WRITE(*,'(A)') "day transmitted"
!~              write(*,*) "TRANSMIT = ", atm_transmissivity(shift_num_days,ij), ij
!~              write(*,*) "TOA ==", fswdtoa0(j,i), SUM(solflux_hourly(:,j),dim=1)/6.
!~            endif

#if (ISOATM >= 2)
            dhu(shift_num_days,ij) = dhu(shift_num_days,ij) + d_vap(j,i)
            dhu17O(shift_num_days,ij) = dhu17O(shift_num_days,ij) + d_vap_17O(j,i)
            dhu_D(shift_num_days,ij) = dhu_D(shift_num_days,ij) + d_vap_D(j,i)
            dom(shift_num_days,ij) = dom(shift_num_days,ij) + d_gw(j,i)
            dom17(shift_num_days,ij) = dom17(shift_num_days,ij) + d_gw17(j,i)
            dom_D(shift_num_days,ij) = dom_D(shift_num_days,ij) + d_gw_D(j,i)
            dp_O(shift_num_days,ij) = dp_O(shift_num_days,ij) + d_precip(j,i)
            dp_D(shift_num_days,ij) = dp_D(shift_num_days,ij) + d_precip_D(j,i)
#endif

         end if
      end do
    end do


    if ( curr_step  == steps_inday ) then
       if (ANY( shift_indexes == shift_num_days) ) then

          current_shift = current_shift + 1
          shift_num_days = num_days + current_shift

          temp_air(shift_num_days,:) = temp_air(shift_num_days-1,:)
          pprc(shift_num_days,:) = pprc(shift_num_days-1,:)
          shr(shift_num_days,:) = shr(shift_num_days-1,:)
          rhum(shift_num_days,:) = rhum(shift_num_days-1,:)
          wndspeed(shift_num_days,:) = wndspeed(shift_num_days-1,:)

#if ( HOURLY_RAD > 0 )
          solar_flux_hourly(shift_num_days,:,:) = solar_flux_hourly(shift_num_days-1,:,:)
#endif

#if (ISOATM >= 2)
          dhu(shift_num_days,:) = dhu(shift_num_days-1,:)
          dhu17O(shift_num_days,:) = dhu17O(shift_num_days-1,:)
          dhu_D(shift_num_days,:) = dhu_D(shift_num_days-1,:)
          dom(shift_num_days,:) = dom(shift_num_days-1,:)
          dom17(shift_num_days,:) = dom17(shift_num_days-1,:)
          dom_D(shift_num_days,:) = dom_D(shift_num_days-1,:)
          dp_O(shift_num_days,:) = dp_O(shift_num_days-1,:)
          dp_D(shift_num_days,:) = dp_D(shift_num_days-1,:)
#endif

          write(*,*) "NUM_day in ec_co2ca ##", shift_num_days, curr_step, current_shift

       endif
    endif

#if ( CARAIB_FORC_W > 0 )

        open(newunit=gauss_asc_id, file='inputdata/gauss.asc',status='old',form='formatted')
        rewind(gauss_asc_id)
        ilat=nlat/2
	10 continue
        read(gauss_asc_id,220,end=15) rj,rw
        ird=int(rj)
        if (ird.eq.ilat) then
           do k=1,ird
              read(gauss_asc_id,220) ari(k),dumwei
           enddo
           goto 20
        else
           goto 10
        endif
	15 continue
	20 continue

        do j=1,ilat
            phi(j)=-ari(ilat+1-j)
            phi(ilat+j)=ari(j)
        enddo
        do j=1,nlat
            phi(j)=asin(phi(j))
        enddo

220  format(f18.10,f17.10)
        close(gauss_asc_id)

        allocate(templon(nlon))
        templon = (/ ((360.0d0*l)/nlon,l=0,nlon-1) /)
        allocate(templat(nlat))
        templat = phi(1:nlat)*radian_to_degree

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
#if (ISOATM >= 2)
             dhu_f%isFormatted = dhu_fm
             call open_f(dhu_f, trim(dhu_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
             dhu17O_f%isFormatted = dhu17O_fm
             call open_f(dhu17O_f, trim(dhu17O_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
             dhuD_f%isFormatted = dhuD_fm
             call open_f(dhuD_f, trim(dhuD_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
             dom_f%isFormatted = dom_fm
             call open_f(dom_f, trim(dom_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
             dom17_f%isFormatted = dom17_fm
             call open_f(dom17_f, trim(dom17_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
             domD_f%isFormatted = domD_fm
             call open_f(domD_f, trim(domD_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
             dpD_f%isFormatted = dpD_fm
             call open_f(dpD_f, trim(dpD_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
             dp_O_f%isFormatted = dp_O_fm
             call open_f(dp_O_f, trim(dp_O_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
#endif

             ij = 0

             do j=1,nlat
               do i=1,nlon
                 if(fractn(j,i,nld) .GE. frac_land .AND. j .GE. j_antarctique ) then

                     ij = ij + 1

                     write(tempair_f%id,'(370f10.4)') templon(i),templat(j),temp_air(1:i_days_yr,ij)
                     write(pprc_f%id,'(370f10.4)') templon(i),templat(j),pprc(1:i_days_yr,ij)
                     write(shr_f%id,'(370f10.4)') templon(i),templat(j),shr(1:i_days_yr,ij)
                     write(rhum_f%id,'(370f10.4)') templon(i),templat(j),rhum(1:i_days_yr,ij)
                     write(wndspeed_f%id,'(370f10.4)') templon(i),templat(j),wndspeed(1:i_days_yr,ij)

#if ( HOURLY_RAD > 0 )
!dmr NEED TO ADD HERE THE NECESSARY CREATION OF SOLAR_FLUX_HOURLY file ...
!dmr NEED TO ADD HERE THE NECESSARY CREATION OF TRANSMISSIVITY file ...
#endif

#if (ISOATM >= 2)
                     write(dhu_f%id,'(370f10.4)') templon(i),templat(j),dhu(1:i_days_yr,ij)
                     write(dhu17O_f%id,'(370f10.4)') templon(i),templat(j),dhu17O(1:i_days_yr,ij)
                     write(dhuD_f%id,'(370f10.4)') templon(i),templat(j),dhu_D(1:i_days_yr,ij)
                     write(dom_f%id,'(370f10.4)') templon(i),templat(j),dom(1:i_days_yr,ij)
                     write(dom17_f%id,'(370f10.4)') templon(i),templat(j),dom17(1:i_days_yr,ij)
                     write(domD_f%id,'(370f10.4)') templon(i),templat(j),dom_D(1:i_days_yr,ij)
                     write(dpD_f%id,'(370f10.4)') templon(i),templat(j),dp_D(1:i_days_yr,ij)
                     write(dp_O_f%id,'(370f10.4)') templon(i),templat(j),dp_O(1:i_days_yr,ij)
#endif

                 end if
               end do
             end do

             call close_f(tempair_f)
             call close_f(pprc_f)
             call close_f(shr_f)
             call close_f(rhum_f)
             call close_f(wndspeed_f)
#if (ISOATM >= 2)
             call close_f(dhu_f)
             call close_f(dhu17O_f)
             call close_f(dhuD_f)
             call close_f(dom_f)
             call close_f(dom17_f)
             call close_f(domD_f)
             call close_f(dp_O_f)
             call close_f(dpD_f)
#endif

          endif
        endif

        deallocate(templon)
        deallocate(templat)

#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  end subroutine climvars_acc


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [co2_for_transient]
!
!>     @brief This subroutine define pointer for CO2

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  subroutine co2_for_CARAIB

    use comphys, only: ghg

    implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        stdCO2 = ghg(1)
        write(*,*) "PCO2 IN CO2CA ===", stdCO2

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  end subroutine co2_for_CARAIB


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [orb_for_transient]
!
!>     @brief This subroutine define pointers for orbital parameters

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

   subroutine orbitals_for_CARAIB

     use global_constants_mod, only  : pi_dp, rad_to_deg
     use comphys, only: ecc, obl, omweb

     implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         orb_ecc = ecc
         orb_obl = obl
         orb_per = omweb
         write(*,*) "ORBITAL PARAMETERS IN CO2CA ===", orb_ecc, orb_obl, orb_per

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

   end subroutine orbitals_for_CARAIB


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [reset_climvars]
!
!>     @brief This subroutine reset the values of the climatic variables and use random numbers to add 5 days to LOVECLIM
!>     simulation
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  subroutine reset_climvars()

    use randoms, only: get_random_array

    implicit none

    integer         :: i

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

    temp_air(:,:) = 0.0_sp
    pprc(:,:)     = 0.0_sp
    shr(:,:)      = 0.0_sp
    rhum(:,:)     = 0.0_sp
    wndspeed(:,:) = 0.0_sp
#if ( HOURLY_RAD > 0 )
    solar_flux_hourly(:,:,:) = 0.0_sp
#endif

#if ( DOWNSCALING == 2 )
    temp_air_sg2ca(:,:) = 0.0_dp
    pprc_sg2ca(:,:)     = 0.0_dp
    shr_sg2ca(:,:)      = 0.0_dp
    rhum_sg2ca(:,:)     = 0.0_dp
    wndspeed_sg2ca(:,:) = 0.0_dp
    ! swrs_sg2ca(:,:)     = 0.0_sp
    ! lwrs_sg2ca(:,:)     = 0.0_sp
#endif

#if (ISOATM >=2 )
    dhu(:,:)      = 0.0_sp
    dhu17O(:,:)   = 0.0_sp
    dhu_D(:,:)    = 0.0_sp
    dom(:,:)      = 0.0_sp
    dom17(:,:)    = 0.0_sp
    dom_D(:,:)    = 0.0_sp
    dp_O(:,:)     = 0.0_sp
    dp_D(:,:)     = 0.0_sp
#endif

    shift_indexes = 0_ip

    do while (MINVAL(abs(shift_indexes-cshift(shift_indexes,1))).lt.2)
       call get_random_array(shift_indexes,sort=.true.)
    enddo

    write(*,*) "shift_indexes == "
    write(*,*) shift_indexes
!~       write(*,*) cshift(shift_indexes,1)
!~       write(*,*) abs(shift_indexes-cshift(shift_indexes,1))

    current_shift = 0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  end subroutine reset_climvars

#if (CARAIB > 0)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [fait_pointer]
!
!>     @brief This subroutine creates pointers for CARAIB targets

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  subroutine sync_coupler_caraib(tcel1year,prc1year,sunhour1year,rhu1year,win1year, orbit_ecc, orbit_obl, orbit_per, CO2year     &
                                ,dhu1year,D17O1year,dhuD1year,dom1year,domD1year,dom1year17O,dpD1year,trm1year,sfx1year)

       real(sp), dimension(:,:), intent(out) :: tcel1year,prc1year,sunhour1year,rhu1year,win1year
       real(sp)                , intent(out) :: orbit_ecc, orbit_obl, orbit_per, CO2year
       real(sp), dimension(:,:), intent(out), optional :: dhu1year,D17O1year,dhuD1year,dom1year,domD1year,dom1year17O,dpD1year
       real(sp), dimension(:,:), intent(out), optional :: trm1year
       real(sp), dimension(:,:,:),intent(out),optional :: sfx1year

!      use MOD_NETCDFCARAIB, only: tcel1year,prc1year,sunhour1year,rhu1year,win1year,dhu1year,dom1year ! --- dmr modif temporaire ,pco2_rd,exc,obl,xlsper

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       tcel1year(:,:)    = temp_air(:,:)
       prc1year(:,:)     = pprc(:,:)
       sunhour1year(:,:) = shr(:,:)
       rhu1year(:,:)     = rhum(:,:)
       win1year(:,:)     = wndspeed(:,:)

       orbit_ecc         = orb_ecc
       orbit_obl         = orb_obl
       orbit_per         = orb_per
       CO2year           = stdCO2

#if ( HOURLY_RAD > 0 )
! Here add transmission of solar flux TOA hourly dimension 12,365,npix
       sfx1year(:,:,:)   = solar_flux_hourly(:,:,:)
! Here add transmission of atmospheric transmissivity values ...
       trm1year(:,:)     = atm_transmissivity(:,:)
#endif

#if ( ISOATM >=2 )
       if (present(dom1year)) then
         dom1year(:,:)     = dom(:,:)
       endif
       if (present(domD1year)) then
         domD1year(:,:)    = dom_D(:,:)
       endif
       if (present(dom1year17O)) then
         dom1year17O(:,:)  = dom17(:,:)
       endif
       if (present(dhu1year)) then
         dhu1year(:,:)     = dhu(:,:)
       endif
       if (present(D17O1year)) then
         D17O1year(:,:)     = dhu17O(:,:)
       endif
       if (present(dhuD1year)) then
         dhuD1year(:,:)    = dhu_D(:,:)
       endif
       if (present(dpD1year)) then
         dpD1year(:,:)     = dp_D(:,:)
       endif
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

  end subroutine sync_coupler_caraib

#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module ec_co2ca

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
