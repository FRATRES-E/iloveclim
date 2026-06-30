#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: Vegetation_Output
!
!>     @author  !mohr, BdB (05-2019), refactoring by <author>
!
!>     @brief   NetCDF output module for the VECODE vegetation model.
!
!>     @details Provides open / close / write / output to write annual vegetation
!>              fields to a NetCDF file.  Internal state (file ID, initialisation
!>              flag) is module-level save.  All grid metadata (phi, numvegvar, …)
!>              is read from veget_mod; time metadata comes from comemic_mod.
!
!>     @date    Original : !mohr
!>     @date    Time axis fix : BdB, 05-2019
!>     @date    Refactored to clean F90 : <date>
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module Vegetation_Output

        use comemic_mod,only: iyear, imonth, fini, irunlabel, nwrskip,  &
                              new_year_veg, time_in_years, current_int_veg
        use comatm,     only: nlat, nlon
        use veget_mod,  only: phi, numvegvar, newvegvar, namevegvar,     &
                              veg_fill_value, veg_missing_value
        use error0_mod, only: ec_error

        implicit none

        include 'netcdf.inc'

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Public interface
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        private

        public :: open, close, write, output, outputStdDev
        public :: ncid, status
        public :: Monthly_Means, Yearly_Means
        public :: vegetzav_id, vegetoutp_id

        ! field index constants (used in veget_wr to identify output variables)
        public :: tree_fraction, grass_fraction, desert_fraction
        public :: needle_tree_fraction
        public :: tree_leaf_area_index, grass_leaf_area_index
        public :: vegetation_albedo, surface_temperature
        public :: Annual_gdd0_index, precipitation, reduced_precipitaion
        public :: Biomass_carbon_uptake, Biomass_carbon_stock
        public :: leaves_biomass_carbon_stock, stems_root_biomass_stock
        public :: litter_carbon_stock
        public :: mortmass_and_soil_organic_matter_carbon_stock
        public :: net_primary_productivity, ice_sheet_fraction

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Derived types for NetCDF metadata
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        type :: Info
          character(len=8)  :: short
          character(len=56) :: std
          character(len=56) :: long
          character(len=64) :: unit
          character(len=8)  :: axis
          character(len=56) :: more
        end type Info

        type :: Infovar
          character(len=8)  :: short
          character(len=56) :: std
          character(len=56) :: long
          character(len=64) :: unit
          character(len=56) :: fill
          character(len=56) :: miss
        end type Infovar

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Module parameters
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! output file type codes
        integer, parameter :: Instantaneous_Data    = 1
        integer, parameter :: Monthly_Means         = 2
        integer, parameter :: Seasonal_Means        = 3
        integer, parameter :: Total_Monthly_Means   = 4
        integer, parameter :: Total_Seasonal_Means  = 5
        integer, parameter :: Yearly_Means          = 6

        ! time computation mode codes
        integer, parameter :: Compute_Time_in_Months_Only = 1
        integer, parameter :: Compute_Time_in_Years_Only  = 2

        ! coordinate and variable metadata templates
        type(Info), parameter :: &
          I_info           = Info("name", "standard_name", "long_name",     &
                                  "units",                       "axis", ""),&
          I_longitude      = Info("lon",  "longitude",     "Longitude",     &
                                  "degrees_east",               "X",    ""),&
          I_latitude       = Info("lat",  "latitude",      "Latitude",      &
                                  "degrees_north",              "Y",    ""),&
          I_time_in_months = Info("time", "time",          "Model Time",    &
                                  "months since 1850-1-1","T", "360_day"),  &
          I_time_in_years  = Info("time", "time",          "Model Time",    &
                                  "years since 1850-1-1", "T", "360_day")

        type(Infovar), parameter :: &
          I_infovar = Infovar("name", "standard_name", "long_name",         &
                              "units", "_FillValue", "missing_value")

        ! 2-D field index constants
        integer, parameter :: tree_fraction               =  1
        integer, parameter :: grass_fraction              =  2
        integer, parameter :: desert_fraction             =  3
        integer, parameter :: needle_tree_fraction        =  4
        integer, parameter :: tree_leaf_area_index        =  5
        integer, parameter :: grass_leaf_area_index       =  6
        integer, parameter :: vegetation_albedo           =  7
        integer, parameter :: surface_temperature         =  8
        integer, parameter :: Annual_gdd0_index           =  9
        integer, parameter :: precipitation               = 10
        ! note: spelling preserved from original to avoid breaking existing NetCDF files
        integer, parameter :: reduced_precipitaion        = 11
        integer, parameter :: Biomass_carbon_uptake       = 12
        integer, parameter :: Biomass_carbon_stock        = 13
        integer, parameter :: leaves_biomass_carbon_stock = 14
        integer, parameter :: stems_root_biomass_stock    = 15
        integer, parameter :: litter_carbon_stock         = 16
        integer, parameter :: &
          mortmass_and_soil_organic_matter_carbon_stock   = 17
        integer, parameter :: net_primary_productivity    = 18
        integer, parameter :: ice_sheet_fraction          = 19

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Module variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        logical, save :: initialize_module = .true.
        integer, save :: how_to_compute_time
        integer, save :: ncid
        integer       :: status

        ! Private log unit for NetCDF error messages — opened on first use via newunit
        integer, save :: nc_err_id = -1

        ! ASCII file unit IDs (set by veget_wr, read by this module for error reporting)
        integer :: vegetzav_id, vegetoutp_id

        ! coordinate bound arrays (allocated in initializeModule, freed after use)
        real(kind=8), allocatable, dimension(:) :: lon_bounds, lat_bounds, time_bounds

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: write_nc_err  (private)
!>     @brief  Write a NetCDF error message to the vegetation log file.
!>             The log file is opened on the first call via newunit so that no
!>             fixed unit number is required (iveg+29 no longer exists).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine write_nc_err(msg)
        character(len=*), intent(in) :: msg
        logical :: opened
        if (nc_err_id < 0) then
          inquire(file='outputdata/vegetation/veget_nc_errors.log', opened=opened)
          if (.not. opened) then
            open(newunit=nc_err_id, file='outputdata/vegetation/veget_nc_errors.log', &
                 status='unknown', position='append')
          end if
        end if
        write(nc_err_id, '(A)') trim(msg)
      end subroutine write_nc_err

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: open
!>     @brief  Open (or create) the NetCDF output file for the current restart interval.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine open(file)

        integer, intent(in) :: file

        character(len=256) :: output_filename
        logical            :: existe
        character(len=6), save :: finichf

1       format(i6.6)

        ! update filename label at run start and at each restart boundary
        if (iyear == 1) finichf = fini
        if (mod(iyear-1, nwrskip) == 0) write(finichf, 1) irunlabel + iyear - 1

        select case (file)
        case (Monthly_Means)
          output_filename     = 'outputdata/vegetation/vegmmyl' // finichf // '.nc'
          how_to_compute_time = Compute_Time_in_Months_Only
        case (Yearly_Means)
          output_filename     = 'outputdata/vegetation/vegym' // finichf // '.nc'
          how_to_compute_time = Compute_Time_in_Years_Only
        case default
          call ec_error(123)
        end select

        inquire(file=output_filename, exist=existe)
        initialize_module = .not. existe

        if (existe) then
          status = nf_open(output_filename, NF_WRITE, ncid)
        else
          status = nf_create(output_filename, nf_clobber, ncid)
        end if
        if (status /= nf_noerr) then
          call write_nc_err(nf_strerror(status))
          call ec_error(123)
        end if

        if (initialize_module) call initializeModule()

      end subroutine open

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: close
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine close()

        if (initialize_module) call initializeModule()
        status = NF_CLOSE(ncid)
        if (status /= nf_noerr) then
          call write_nc_err(nf_strerror(status))
          call ec_error(123)
        end if

      end subroutine close

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: incretime  (private)
!>     @brief  Write time coordinate values and bounds for N_num time steps.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine incretime(N_num)

        integer, intent(in) :: N_num

        integer :: varid, i
        integer :: istart(2), cnt(2)
        real(kind=8) :: D_data(N_num)
        real(kind=8) :: D_data2(2,N_num)

        istart(:) = 1
        cnt(:)    = N_num
        status = NF_INQ_VARID(ncid, 'time', varid)

        do i = 1, N_num
          select case (how_to_compute_time)
          case (Compute_Time_in_Months_Only)
            D_data(i) = time_in_years(current_int_veg) + real(i, 8) / 12.0d0
          case (Compute_Time_in_Years_Only)
            D_data(i) = time_in_years(current_int_veg + i)
          end select
        end do
        status = NF_PUT_VARA_DOUBLE(ncid, varid, istart, cnt, D_data)

        status = NF_INQ_VARID(ncid, 'time_bounds', varid)
        cnt(1) = 2
        do i = 1, N_num
          select case (how_to_compute_time)
          case (Compute_Time_in_Months_Only)
            D_data2(1,i) = time_in_years(current_int_veg) + (real(i,8) - 0.5d0) / 12.0d0
            D_data2(2,i) = time_in_years(current_int_veg) + (real(i,8) + 0.5d0) / 12.0d0
          case (Compute_Time_in_Years_Only)
            D_data2(1,i) = time_in_years(current_int_veg + i) - 0.5d0
            D_data2(2,i) = time_in_years(current_int_veg + i) + 0.5d0
          end select
        end do
        status = NF_PUT_VARA_DOUBLE(ncid, varid, istart, cnt, D_data2)
        if (status /= nf_noerr) then
          call write_nc_err(nf_strerror(status))
          call ec_error(123)
        end if

      end subroutine incretime

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: write
!>     @brief  Write a 2-D field slice (nlat×nlon) at the current time step.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine write(numid, D_data)

        integer,                       intent(in) :: numid
        real(kind=8), dimension(:,:),  intent(in) :: D_data

        integer :: tab_dim(5), istart(5)
        integer :: l, varid, dimid
        logical :: write_test
        real(kind=8), allocatable :: tmp(:,:,:)
        integer :: t

        write_test = .false.
        istart(:)  = 1
        tab_dim(:) = 1

        select case (how_to_compute_time)
        case (Compute_Time_in_Months_Only)
          write_test = (newvegvar(numid,3) == 1.0d0)
        case (Compute_Time_in_Years_Only)
          write_test = (newvegvar(numid,4) == 1.0d0)
        end select

        if (.not. write_test) return

        status = NF_INQ_VARID(ncid, namevegvar(numid,2), varid)
        tab_dim(1) = nlon
        tab_dim(2) = nlat

        status = NF_INQ_DIMID(ncid, 'time', dimid)
        status = NF_INQ_DIMLEN(ncid, dimid, l)
        t = computeTime()
        istart(3) = t
        if (l /= t) call incretime(t)

        allocate(tmp(nlon, nlat, t))
        do l = 1, nlat
          tmp(:,l,1) = D_data(l,:)
        end do
        status = NF_PUT_VARA_DOUBLE(ncid, varid, istart, tab_dim, tmp)
        deallocate(tmp)

        if (status /= nf_noerr) then
          call write_nc_err(nf_strerror(status))
          call ec_error(123)
        end if

      end subroutine write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      FUNCTION: output
!>     @brief  Return .true. if a field's output flag is non-zero.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function output(flag) result(res)
        real(kind=8), intent(in) :: flag
        logical :: res
        res = (flag /= 0.0d0)
      end function output

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      FUNCTION: outputStdDev
!>     @brief  Return .true. if a field's flag requests standard deviation output.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function outputStdDev(flag) result(res)
        integer, intent(in) :: flag
        logical :: res
        res = (flag == 2)
      end function outputStdDev

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: initializeModule  (private)
!>     @brief  Define NetCDF dimensions, coordinate variables, bounds, and field variables.
!>             Called once per output file (on creation).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine initializeModule()

        real(kind=8), parameter :: pi               = 2.0d0 * acos(0.0d0)
        real(kind=8), parameter :: radian_to_degree = 180.0d0 / pi

        real(kind=8), allocatable :: templon(:), templat(:), temp(:,:)
        integer :: tab_dimid(5), tab_dim(5)
        integer :: i, l, dimid, varid
        logical :: write_test
        character(len=120) :: longdata
        integer :: cvalues(8)

999     format(i4.4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',i2.2,'Z')

        if (.not. initialize_module) return

        !--- 1. longitude coordinate ---
        allocate(templon(nlon))
        templon = [(360.0d0 * l / nlon, l=0, nlon-1)]

        allocate(lon_bounds(nlon+1))
        lon_bounds(1) = templon(1) - (templon(2) - templon(1)) / 2.0d0
        do i = 2, nlon
          lon_bounds(i) = templon(i-1) - (templon(i-1) - templon(i)) / 2.0d0
        end do
        lon_bounds(nlon+1) = templon(nlon) + &
                             (lon_bounds(nlon) - lon_bounds(nlon-1)) / 2.0d0

        call initdim(I_longitude, nlon, .false., .false., templon)
        status = NF_DEF_DIM(ncid, 'bnds', 2, dimid)
        tab_dimid(1) = dimid

        allocate(temp(2,nlon))
        do i = 1, nlon
          temp(1,i) = lon_bounds(i)
          temp(2,i) = lon_bounds(i+1)
        end do
        status = NF_INQ_DIMID(ncid, I_longitude%short, dimid)
        tab_dimid(2) = dimid
        longdata = 'lon_bounds'
        call initbnds(longdata, nlon, tab_dimid, temp)
        deallocate(temp, lon_bounds)

        !--- 2. latitude coordinate ---
        allocate(templat(nlat))
        templat = phi(1:nlat) * radian_to_degree

        allocate(lat_bounds(nlat+1))
        lat_bounds(1) = templat(1) - (lat_bounds(3) - lat_bounds(2)) / 2.0d0
        do i = 2, nlat
          lat_bounds(i) = templat(i-1) - (templat(i-1) - templat(i)) / 2.0d0
        end do
        lat_bounds(nlat+1) = templat(nlat) + &
                             (lat_bounds(nlat) - lat_bounds(nlat-1)) / 2.0d0

        call initdim(I_latitude, nlat, .false., .false., templat)

        allocate(temp(2,nlat))
        do i = 1, nlat
          temp(1,i) = lat_bounds(i)
          temp(2,i) = lat_bounds(i+1)
        end do
        status = NF_INQ_DIMID(ncid, I_latitude%short, dimid)
        tab_dimid(2) = dimid
        longdata = 'lat_bounds'
        call initbnds(longdata, nlat, tab_dimid, temp)
        deallocate(temp, lat_bounds)

        !--- 3. time coordinate ---
        select case (how_to_compute_time)
        case (Compute_Time_in_Months_Only)
          call initdim(I_time_in_months, NF_UNLIMITED, .true., .false., [1.0d0])
        case (Compute_Time_in_Years_Only)
          call initdim(I_time_in_years,  NF_UNLIMITED, .true., .false., [1.0d0])
        end select

        allocate(time_bounds(2))
        time_bounds = [0.5d0, 1.5d0]
        allocate(temp(2,1))
        temp(1,1) = time_bounds(1)
        temp(2,1) = time_bounds(2)
        status = NF_INQ_DIMID(ncid, I_time_in_years%short, dimid)
        tab_dimid(2) = dimid
        longdata = 'time_bounds'
        call initbnds(longdata, 1, tab_dimid, temp)
        deallocate(temp, time_bounds, templon, templat)
        status = NF_ENDDEF(ncid)

        !--- 4. field variables ---
        write_test = .false.
        do i = 1, numvegvar
          select case (how_to_compute_time)
          case (Compute_Time_in_Months_Only)
            write_test = (newvegvar(i,3) == 1.0d0)
          case (Compute_Time_in_Years_Only)
            write_test = (newvegvar(i,4) == 1.0d0)
          end select

          if (write_test) then
            status = NF_REDEF(ncid)
            status = NF_INQ_DIMID(ncid, I_longitude%short, dimid)
            tab_dimid(1) = dimid ;  tab_dim(1) = nlon
            status = NF_INQ_DIMID(ncid, I_latitude%short, dimid)
            tab_dimid(2) = dimid ;  tab_dim(2) = nlat
            status = NF_INQ_DIMID(ncid, I_time_in_years%short, dimid)
            tab_dimid(3) = dimid ;  tab_dim(3) = 1
            call initvar(3, namevegvar(i,:), tab_dimid)
            status = NF_ENDDEF(ncid)
          end if
          write_test = .false.
        end do

        !--- 5. global attributes (creation date) ---
        status = NF_REDEF(ncid)
        call date_and_time(values=cvalues)
        longdata = ''
        write(longdata, 999) cvalues(1), cvalues(2), cvalues(3), &
                             cvalues(5), cvalues(6), cvalues(7)
        status = NF_ENDDEF(ncid)

        initialize_module = .false.

      contains

        subroutine initdim(I_num, N_num, Is_time, Is_pressure, D_data)
          type(Info), intent(in) :: I_num
          integer,    intent(in) :: N_num
          logical,    intent(in) :: Is_time, Is_pressure
          real(kind=8), intent(in) :: D_data(N_num)
          integer :: ldimid, lvarid, lstart(1), lcnt(1)

          status = NF_DEF_DIM(ncid, I_num%short, N_num, ldimid)
          lstart(1) = ldimid
          status = NF_DEF_VAR(ncid, I_num%short, NF_FLOAT, 1, lstart, lvarid)
          status = NF_PUT_ATT_TEXT(ncid, lvarid, I_info%std,  &
                   len_trim(I_num%std),  trim(I_num%std))
          status = NF_PUT_ATT_TEXT(ncid, lvarid, I_info%long, &
                   len_trim(I_num%long), trim(I_num%long))
          status = NF_PUT_ATT_TEXT(ncid, lvarid, I_info%axis, &
                   len_trim(I_num%axis), trim(I_num%axis))
          if (Is_time) then
            status = NF_PUT_ATT_TEXT(ncid, lvarid, I_info%unit, &
                     len_trim(I_num%unit), trim(I_num%unit))
            status = NF_PUT_ATT_TEXT(ncid, lvarid, 'calendar', &
                     len_trim(I_num%more), trim(I_num%more))
          else
            status = NF_PUT_ATT_TEXT(ncid, lvarid, I_info%unit, &
                     len_trim(I_num%unit), trim(I_num%unit))
          end if
          if (Is_pressure) status = NF_PUT_ATT_TEXT(ncid, lvarid, 'positive', &
                                     len_trim(I_num%more), trim(I_num%more))
          status = NF_ENDDEF(ncid)
          lstart(1) = 1 ;  lcnt(1) = N_num
          if (Is_time) lcnt(1) = 1
          status = NF_PUT_VARA_DOUBLE(ncid, lvarid, lstart, lcnt, D_data)
          status = NF_REDEF(ncid)
          if (status /= nf_noerr) then
            call write_nc_err(nf_strerror(status))
            call ec_error(123)
          end if
        end subroutine initdim

        subroutine initbnds(name, I_num, tab_dimid_in, D_data)
          character(len=*),         intent(in) :: name
          integer,                  intent(in) :: I_num, tab_dimid_in(2)
          real(kind=8), dimension(:,:), intent(in) :: D_data
          integer :: lvarid, lstart(2), lcnt(2)
          status = NF_DEF_VAR(ncid, name, NF_DOUBLE, 2, tab_dimid_in, lvarid)
          status = NF_ENDDEF(ncid)
          lstart(:) = 1 ;  lcnt(1) = 2 ;  lcnt(2) = I_num
          status = NF_PUT_VARA_DOUBLE(ncid, lvarid, lstart, lcnt, D_data)
          status = NF_REDEF(ncid)
          if (status /= nf_noerr) then
            call write_nc_err(nf_strerror(status))
            call ec_error(123)
          end if
        end subroutine initbnds

        subroutine initvar(N_num, tab_c, tab_dimid_in)
          integer,       intent(in) :: N_num, tab_dimid_in(5)
          character(len=60), intent(in) :: tab_c(5)
          integer :: lvarid
          real(kind=8) :: tmp(1)
          status = NF_DEF_VAR(ncid, tab_c(2), NF_DOUBLE, N_num, tab_dimid_in, lvarid)
          status = NF_PUT_ATT_TEXT(ncid, lvarid, I_infovar%std,  &
                   len_trim(tab_c(3)), trim(tab_c(3)))
          status = NF_PUT_ATT_TEXT(ncid, lvarid, I_infovar%long, &
                   len_trim(tab_c(1)), trim(tab_c(1)))
          status = NF_PUT_ATT_TEXT(ncid, lvarid, I_infovar%unit, &
                   len_trim(tab_c(4)), trim(tab_c(4)))
          tmp(1) = veg_fill_value
          status = NF_PUT_ATT_DOUBLE(ncid, lvarid, I_infovar%fill, NF_DOUBLE, 1, tmp)
          tmp(1) = veg_missing_value
          status = NF_PUT_ATT_DOUBLE(ncid, lvarid, I_infovar%miss, NF_DOUBLE, 1, tmp)
          if (status /= nf_noerr) then
            call write_nc_err(nf_strerror(status))
            call ec_error(123)
          end if
        end subroutine initvar

      end subroutine initializeModule

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      FUNCTION: computeTime  (private)
!>     @brief  Return the current time index within the output file.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function computeTime() result(res)
        integer :: res
        select case (how_to_compute_time)
        case (Compute_Time_in_Months_Only)
          res = new_year_veg * 12 + imonth
        case (Compute_Time_in_Years_Only)
          res = new_year_veg
        end select
      end function computeTime

      end module Vegetation_Output

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
