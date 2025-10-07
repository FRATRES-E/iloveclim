!dmr -- Added optional components choice - Fri Oct 28 10:57:45 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Fri Oct 28 10:57:45 CEST 2011
      MODULE Atmosphere_Output

!      USE uuid_fort_wrap, ONLY: UUID_V4, uuid_size
      use uuid_module, only: generate_uuid
      use global_constants_mod, only: uuid_size

#if ( COMATM == 1 )
      USE comatm
      USE comdiag
      use comemic_mod, only: iyear, imonth, iday, fini, irunlabel, irunlabeld
     &               , globalatt
      use comunit
#endif

! --- BdB 05-2019: added variables for writing time
      use comemic_mod, only: nwrskip, new_year_atm
     &               , time_in_years, current_int_atm

      use newunit_mod, only: error_id
      
      IMPLICIT NONE
! TODO replace netcdf.inc with "use netcdf"
#include "netcdf.inc"

#if ( COMATM == 0 )
#include "comatm.h"
#include "comdiag.h"
#include "comemic.h"
#include "comunit.h"
#endif


!     ----------------------------------------------------------------------------
!     Module Interface Specification

      PRIVATE
      PUBLIC :: open
      PUBLIC :: close
      PUBLIC :: write
      PUBLIC :: ncid
      PUBLIC :: status
      PUBLIC :: output

!     Codes for output file types
      PUBLIC :: Instantaneous_Data
      PUBLIC :: Monthly_Means
      PUBLIC :: Seasonal_Means
      PUBLIC :: Total_Monthly_Means
      PUBLIC :: Total_Seasonal_Means
      PUBLIC :: Yearly_Means

!     Codes for one dimensional fields
      REAL*8, ALLOCATABLE, DIMENSION(:) :: lon_bounds, lat_bounds, time_bounds

!     Codes for one dimensional fields
      PUBLIC :: lon
      PUBLIC :: lat
      PUBLIC :: P_T2
      PUBLIC :: P_T3
      PUBLIC :: P_T4
      PUBLIC :: P_U3
      PUBLIC :: time

!     Codes for two dimensional fields
      PUBLIC :: Specific_Humidity
      PUBLIC :: Surface_Temperature
      PUBLIC :: Stratospheric_Temperature
      PUBLIC :: Two_Meter_Temperature
      PUBLIC :: Bottom_Moisture
      PUBLIC :: Land_Snow_Depth
      PUBLIC :: Large_Scale_Precipitation
      PUBLIC :: Convective_Precipitation
      PUBLIC :: Surface_Sensible_Heat_Flux
      PUBLIC :: Surface_Latent_Heat_Flux
      PUBLIC :: Relative_Humidity
      PUBLIC :: Wind_at_10_Meter
      PUBLIC :: Ocean_Surface_Runoff
      PUBLIC :: Land_Surface_Runoff
      PUBLIC :: Total_Cloud_Cover
      PUBLIC :: Planetary_Albedo
      PUBLIC :: Surface_Albedo
      PUBLIC :: Surface_Solar_Radiation
      PUBLIC :: Top_Solar_Radiation
      PUBLIC :: Surface_Thermal_Radiation
      PUBLIC :: Top_Thermal_Radiation
      PUBLIC :: U_Stress
      PUBLIC :: V_Stress
      PUBLIC :: Surface_Evaporation
      PUBLIC :: Total_Precipitation
      PUBLIC :: Surface_Pressure
      PUBLIC :: Evap_Minus_Precip
      PUBLIC :: Sea_Ice_Thickness
      PUBLIC :: Drag_Coefficient_W
      PUBLIC :: Drag_Coefficient_V
      PUBLIC :: Total_Snow_Fall
      PUBLIC :: Surface_Shortwave_Radiation !mohr
      PUBLIC :: Surface_Longwave_Radiation  !mohr

#if ( ISOATM >= 1 )
      PUBLIC :: Total_Precip18
      PUBLIC :: Total_Snow18
      PUBLIC :: Total_Precip17
      PUBLIC :: Total_Snow17
      PUBLIC :: Total_Precipd
      PUBLIC :: Total_Snowd
      PUBLIC :: Bottom_Mois18
      PUBLIC :: Bottom_Mois17
      PUBLIC :: Bottom_Moisd
      PUBLIC :: Land_Snow_D18
      PUBLIC :: Land_Snow_D17
      PUBLIC :: Land_Snow_Dd
      PUBLIC :: Spec_Hum_d18
      PUBLIC :: Spec_Hum_d17
      PUBLIC :: Spec_Hum_dd
#endif

!     Codes for three dimensional fields
      PUBLIC :: Temperature
      PUBLIC :: Wind_U
      PUBLIC :: Wind_V
      PUBLIC :: Vertical_Pressure_Wind
      PUBLIC :: Stream_Function
      PUBLIC :: Velocity_Potential
      PUBLIC :: Heating_Force
      PUBLIC :: Potential_Vorticity_Forcing
      PUBLIC :: Ageostrophic_Wind_U
      PUBLIC :: Ageostrophic_Wind_V
      PUBLIC :: QG_Potential_Vorticity
      PUBLIC :: Geopotential_Height
      PUBLIC :: User_Assigned_T1
      PUBLIC :: User_Assigned_T2
      PUBLIC :: User_Assigned_U1
      PUBLIC :: User_Assigned_U2

!      INTERFACE write
!      MODULE PROCEDURE writeField2D
!      MODULE PROCEDURE writeField3D
!      END INTERFACE

!     ----------------------------------------------------------------------------
!     Type Declarations

      TYPE Info
      CHARACTER(len=8)  :: short
      CHARACTER(len=56) :: std
      CHARACTER(len=56) :: long
      CHARACTER(len=64) :: unit
      CHARACTER(len=8)  :: axis
      CHARACTER(len=56) :: more
      END TYPE Info

      TYPE Infovar
      CHARACTER(len=8)  :: short
      CHARACTER(len=56) :: std
      CHARACTER(len=56) :: long
      CHARACTER(len=64) :: unit
      CHARACTER(len=56) :: fill
      CHARACTER(len=56) :: miss
      END TYPE Infovar

!     ----------------------------------------------------------------------------
!     Module Parameters

!     Output file types
      INTEGER, PARAMETER :: Instantaneous_Data        =  1
      INTEGER, PARAMETER :: Monthly_Means             =  2
      INTEGER, PARAMETER :: Seasonal_Means            =  3
      INTEGER, PARAMETER :: Total_Monthly_Means       =  4
      INTEGER, PARAMETER :: Total_Seasonal_Means      =  5
      INTEGER, PARAMETER :: Yearly_Means              =  6

!     Grid coordinates
      TYPE(Info), PARAMETER ::
     & I_info           = Info("name",  "standard_name", "long_name",
     &          "units",                            "axis", ""),
     & I_longitude      = Info("lon",   "longitude",     "Longitude",
     &          "degrees_east",                     "X",    ""),
     & I_latitude       = Info("lat",   "latitude",      "Latitude",
     &          "degrees_north",                    "Y",    ""),
     & I_pressure_T2    = Info("P_T2",  "pressure",      "Air Pressure",
     &          "Pa",                               "Z",    "down"),
     & I_pressure_T3    = Info("P_T3",  "pressure",      "Air Pressure",
     &          "Pa",                               "Z",    "down"),
     & I_pressure_T4    = Info("P_T4",  "pressure",      "Air Pressure",
     &          "Pa",                               "Z",    "down"),
     & I_pressure_U3    = Info("P_U3",  "pressure",      "Air Pressure",
     &          "Pa",                               "Z",    "down"),
     & I_time_in_days   = Info("time",  "time",          "Model Time",
     &          "days",   "T",    "360_day"),
     & I_time_in_months = Info("time",  "time",          "Model Time",
     &          "months", "T",    "360_day"),
     & I_time_in_years  = Info("time",  "time",          "Model Time",
     &          "years",  "T",    "360_day")
!     Grid coordinates
      TYPE(Infovar), PARAMETER ::
     & I_infovar        = Infovar("name",  "standard_name", "long_name",
     &          "units",                 "_FillValue", "missing_value")

!     One dimensional fields
      INTEGER, PARAMETER :: lon  = 1
      INTEGER, PARAMETER :: lat  = 2
      INTEGER, PARAMETER :: P_T2 = 3
      INTEGER, PARAMETER :: P_T3 = 4
      INTEGER, PARAMETER :: P_T4 = 5
      INTEGER, PARAMETER :: P_U3 = 6
      INTEGER, PARAMETER :: time = 7

!     Two dimensional fields
      INTEGER, PARAMETER :: Specific_Humidity           =  1
      INTEGER, PARAMETER :: Surface_Temperature         =  2
      INTEGER, PARAMETER :: Stratospheric_Temperature   =  3
      INTEGER, PARAMETER :: Two_Meter_Temperature       =  4
      INTEGER, PARAMETER :: Bottom_Moisture             =  5
      INTEGER, PARAMETER :: Land_Snow_Depth             =  6
      INTEGER, PARAMETER :: Large_Scale_Precipitation   =  7
      INTEGER, PARAMETER :: Convective_Precipitation    =  8
      INTEGER, PARAMETER :: Surface_Sensible_Heat_Flux  =  9
      INTEGER, PARAMETER :: Surface_Latent_Heat_Flux    = 10
      INTEGER, PARAMETER :: Relative_Humidity           = 11
      INTEGER, PARAMETER :: Wind_at_10_Meter            = 12
      INTEGER, PARAMETER :: Ocean_Surface_Runoff        = 13
      INTEGER, PARAMETER :: Land_Surface_Runoff         = 14
      INTEGER, PARAMETER :: Total_Cloud_Cover           = 15
      INTEGER, PARAMETER :: Planetary_Albedo            = 16
      INTEGER, PARAMETER :: Surface_Albedo              = 17
      INTEGER, PARAMETER :: Surface_Solar_Radiation     = 18
      INTEGER, PARAMETER :: Top_Solar_Radiation         = 19
      INTEGER, PARAMETER :: Surface_Thermal_Radiation   = 20
      INTEGER, PARAMETER :: Top_Thermal_Radiation       = 21
      INTEGER, PARAMETER :: U_Stress                    = 22
      INTEGER, PARAMETER :: V_Stress                    = 23
      INTEGER, PARAMETER :: Surface_Evaporation         = 24
      INTEGER, PARAMETER :: Total_Precipitation         = 25
      INTEGER, PARAMETER :: Surface_Pressure            = 26
      INTEGER, PARAMETER :: Evap_Minus_Precip           = 27
      INTEGER, PARAMETER :: Sea_Ice_Thickness           = 28
      INTEGER, PARAMETER :: Drag_Coefficient_W          = 29
      INTEGER, PARAMETER :: Drag_Coefficient_V          = 30
      INTEGER, PARAMETER :: Total_Snow_Fall             = 31
      INTEGER, PARAMETER :: Surface_Shortwave_Radiation = 32 !mohr, the value must be 32 in case ISOATM < 1
      INTEGER, PARAMETER :: Surface_Longwave_Radiation  = 33 !mohr
      INTEGER, PARAMETER :: num_noiso_var               = 33
#if ( ISOATM >= 1 )
      INTEGER, PARAMETER :: Total_Precip18          = 34
      INTEGER, PARAMETER :: Total_Snow18            = 35
      INTEGER, PARAMETER :: Total_Precip17          = 36
      INTEGER, PARAMETER :: Total_Snow17            = 37
      INTEGER, PARAMETER :: Total_Precipd           = 38
      INTEGER, PARAMETER :: Total_Snowd             = 39
      INTEGER, PARAMETER :: Bottom_Mois18           = 40
      INTEGER, PARAMETER :: Bottom_Mois17           = 41
      INTEGER, PARAMETER :: Bottom_Moisd            = 42
      INTEGER, PARAMETER :: Land_Snow_D18           = 43
      INTEGER, PARAMETER :: Land_Snow_D17           = 44
      INTEGER, PARAMETER :: Land_Snow_Dd            = 45
      INTEGER, PARAMETER :: Spec_Hum_d18            = 46
      INTEGER, PARAMETER :: Spec_Hum_d17            = 47
      INTEGER, PARAMETER :: Spec_Hum_dd             = 48
#endif

#if ( ISOATM >= 1 )
      INTEGER, PARAMETER :: Number_of_2D_Fields     = 48
#else
      INTEGER, PARAMETER :: Number_of_2D_Fields     = 33
#endif

!dmr ###     INTEGER, PARAMETER :: Number_of_2D_Fields       = 31

#if ( ISOATM >= 1 )
!     Three dimensional fields
      INTEGER, PARAMETER :: Temperature                 =  49
      INTEGER, PARAMETER :: Wind_U                      =  50
      INTEGER, PARAMETER :: Wind_V                      =  51
      INTEGER, PARAMETER :: Vertical_Pressure_Wind      =  52
      INTEGER, PARAMETER :: Stream_Function             =  53
      INTEGER, PARAMETER :: Velocity_Potential          =  54
      INTEGER, PARAMETER :: Heating_Force               =  55
      INTEGER, PARAMETER :: Potential_Vorticity_Forcing =  56
      INTEGER, PARAMETER :: Ageostrophic_Wind_U         =  57
      INTEGER, PARAMETER :: Ageostrophic_Wind_V         = 58
      INTEGER, PARAMETER :: QG_Potential_Vorticity      = 59
      INTEGER, PARAMETER :: Geopotential_Height         = 60
      INTEGER, PARAMETER :: User_Assigned_T1            = 61
      INTEGER, PARAMETER :: User_Assigned_T2            = 62
      INTEGER, PARAMETER :: User_Assigned_U1            = 63
      INTEGER, PARAMETER :: User_Assigned_U2            = 64
#else
!     Three dimensional fields
      INTEGER, PARAMETER :: Temperature                 =  34
      INTEGER, PARAMETER :: Wind_U                      =  35
      INTEGER, PARAMETER :: Wind_V                      =  36
      INTEGER, PARAMETER :: Vertical_Pressure_Wind      =  37
      INTEGER, PARAMETER :: Stream_Function             =  38
      INTEGER, PARAMETER :: Velocity_Potential          =  39
      INTEGER, PARAMETER :: Heating_Force               =  40
      INTEGER, PARAMETER :: Potential_Vorticity_Forcing =  41
      INTEGER, PARAMETER :: Ageostrophic_Wind_U         =  42
      INTEGER, PARAMETER :: Ageostrophic_Wind_V         = 43
      INTEGER, PARAMETER :: QG_Potential_Vorticity      = 44
      INTEGER, PARAMETER :: Geopotential_Height         = 45
      INTEGER, PARAMETER :: User_Assigned_T1            = 46
      INTEGER, PARAMETER :: User_Assigned_T2            = 47
      INTEGER, PARAMETER :: User_Assigned_U1            = 48
      INTEGER, PARAMETER :: User_Assigned_U2            = 49

#endif

      INTEGER, PARAMETER :: Number_of_3D_Fields = 16

      INTEGER, PARAMETER :: Compute_Time_in_Days         = 1
      INTEGER, PARAMETER :: Compute_Time_in_Years_Months = 3
      INTEGER, PARAMETER :: Compute_Time_in_Months_Only  = 2
      INTEGER, PARAMETER :: Compute_Time_in_Years_Only   = 4

!     ----------------------------------------------------------------------------
!     Module Variables

      LOGICAL, SAVE :: initialize_module = .TRUE.
      INTEGER, SAVE :: how_to_compute_time, ncid
      INTEGER       :: status

      CONTAINS

!     ----------------------------------------------------------------------------
!     PUBLIC FUNCTIONS AND SUBROUTINES
!     ----------------------------------------------------------------------------

!     ----------------------------------------------------------------------------
!     open

      SUBROUTINE open(file)
      INTEGER, INTENT(in) :: file

      CHARACTER(len=256) :: output_filename
      LOGICAL            :: existe
      CHARACTER*6, SAVE  :: finichf


! --- BdB 05-2019: update name every nwrskip (the restart output interval)
      ! initialise 
      if (iyear==1) finichf=fini
      ! update the name
      if (mod(iyear-1,nwrskip).eq.0) then
       write(finichf,1) irunlabel+iyear-1
 1     format(i6.6)
      end if

      SELECT CASE (file)
      CASE (Instantaneous_Data)
         output_filename = "outputdata/atmos/atminst"//fini//".nc"
         how_to_compute_time = Compute_Time_in_Days
      CASE (Monthly_Means)
! --- BdB 05-2019: update output name
         output_filename = "outputdata/atmos/atmmmyl"//finichf//".nc"
         how_to_compute_time = Compute_Time_in_Months_Only
      CASE (Seasonal_Means)
         output_filename = "outputdata/atmos/atmsmyl"//fini//".nc"
         how_to_compute_time = Compute_Time_in_Months_Only
      CASE (Total_Monthly_Means)
         output_filename = "outputdata/atmos/atmmmwp"//fini//".nc"
         how_to_compute_time = Compute_Time_in_Years_Months
      CASE (Total_Seasonal_Means)
         output_filename = "outputdata/atmos/atmsmwp"//fini//".nc"
         how_to_compute_time = Compute_Time_in_Years_Months
      CASE (Yearly_Means)
! --- BdB 05-2019: update output name
         output_filename = "outputdata/atmos/atmym"//finichf//".nc"
         how_to_compute_time = Compute_Time_in_Years_Only
      CASE DEFAULT
         call ec_error(123)
      END SELECT

      INQUIRE(file=output_filename,exist=existe)
      IF (existe) THEN
         initialize_module = .FALSE.
      ELSE
         initialize_module = .TRUE.
      END IF

      IF (existe) THEN
         status=nf_open(output_filename,NF_WRITE,ncid)
      ELSE
         status=nf_create(output_filename,nf_clobber,ncid)
      END IF
      IF (status/=nf_noerr) THEN
         write(error_id,*) nf_strerror(status)
         call ec_error(123)
      END IF

      IF (initialize_module) CALL initializeModule

      RETURN
      END SUBROUTINE open

!     ----------------------------------------------------------------------------
!     close

      SUBROUTINE close()

      IF (initialize_module) CALL initializeModule

      status=NF_CLOSE(ncid)
      IF (status/=nf_noerr) THEN
         write(error_id,*) nf_strerror(status)
         call ec_error(123)
      END IF

      RETURN
      END SUBROUTINE close

      SUBROUTINE incretime(N_num)

      INTEGER             :: varid,i
      INTEGER             :: start(2)
      INTEGER             :: count(2)
      INTEGER, INTENT(in) :: N_num
      REAL*8              :: D_data(N_num)
      REAL*8              :: D_data2(2,N_num)
    
!     Output data for Variable
      start(:) = 1
      count(:) = N_num
      status=NF_INQ_VARID(ncid,'time',varid)

      do i=1,N_num
         !D_data(i) = i
! ---    BdB 05-2019: update time relative to initial irunlabel
         SELECT CASE (how_to_compute_time)
! ---    dmr 2020-01-25         
         CASE (Compute_Time_in_Days)
          ! determine day
          D_data(i) = i
         CASE (Compute_Time_in_Months_Only)
          ! determine year
          D_data(i) = time_in_years(current_int_atm) + real(i) / 12.
         CASE (Compute_Time_in_Years_Only)
          D_data(i) = time_in_years(current_int_atm+i)
         END SELECT
      enddo
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,D_data)

      status=NF_INQ_VARID(ncid,"time_bounds",varid)
      count(1) = 2

      do i=1,N_num
         !D_data2(1,i) = i-0.5
         !D_data2(2,i) = i+0.5
! ---    BdB 05-2019: update time relative to initial irunlabel
         SELECT CASE (how_to_compute_time)
! ---    dmr 2020-01-25                  
         CASE (Compute_Time_in_Days)
           ! determine day
           D_data2(1,i) = time_in_years(current_int_atm) + (real(i) - 0.5) / 12.
           D_data2(2,i) = time_in_years(current_int_atm) + (real(i) + 0.5) / 12.
         CASE (Compute_Time_in_Months_Only)
           D_data2(1,i) = time_in_years(current_int_atm) + (real(i) - 0.5) / 12.
           D_data2(2,i) = time_in_years(current_int_atm) + (real(i) + 0.5) / 12.
         CASE (Compute_Time_in_Years_Only)
           D_data2(1,i) = time_in_years(current_int_atm+i) - 0.5
           D_data2(2,i) = time_in_years(current_int_atm+i) + 0.5
         END SELECT
      enddo
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,D_data2)

      IF (status/=nf_noerr) THEN
         write(error_id,*) nf_strerror(status)
         call ec_error(123)
      END IF

      RETURN
      END SUBROUTINE incretime

      SUBROUTINE write(numid,D_data)

      INTEGER, INTENT(in) :: numid
      INTEGER             :: tab_dim(5)
      INTEGER             :: start(5)
      INTEGER             :: l,k,varid,dimid
      LOGICAL             :: write_test

      REAL*8, DIMENSION(:,:,:,:), allocatable :: tmp
      REAL*8, DIMENSION(:,:,:), INTENT(in)    :: D_data

      write_test = .FALSE.
      start(:)   = 1
      tab_dim(:) = 1

      SELECT CASE (how_to_compute_time)
      CASE (Compute_Time_in_Days)
         if ( newtotvar(numid,1)==1 ) write_test = .TRUE.
      CASE (Compute_Time_in_Years_Months)
         if ( newtotvar(numid,3)==1 ) write_test = .TRUE.
      CASE (Compute_Time_in_Months_Only)
         if ( newtotvar(numid,2)==1 ) write_test = .TRUE.
      CASE (Compute_Time_in_Years_Only)
         if ( newtotvar(numid,4)==1 ) write_test = .TRUE.
      END SELECT

      if ( write_test ) then

      status=NF_INQ_VARID(ncid,nametotvar(numid,2),varid)
      tab_dim(1)=nlon
      tab_dim(2)=nlat

      if ( nametotvar(numid,5)/="N" ) then
         SELECT CASE (nametotvar(numid,5))
         CASE ("T2")
           status=NF_INQ_DIMID(ncid,I_pressure_T2%short,dimid)
         CASE ("T3")
           status=NF_INQ_DIMID(ncid,I_pressure_T3%short,dimid)
         CASE ("T4")
           status=NF_INQ_DIMID(ncid,I_pressure_T4%short,dimid)
         CASE ("U3")
           status=NF_INQ_DIMID(ncid,I_pressure_U3%short,dimid)
         END SELECT
         status=NF_INQ_DIMLEN(ncid,dimid,l)
         tab_dim(3) = l
         allocate(tmp (nlon,nlat,l,computeTime()))
         status=NF_INQ_DIMID(ncid,'time',dimid)
         status=NF_INQ_DIMLEN(ncid,dimid,l)
         start(4) = computeTime()
         if ( l/=computeTime() ) call incretime(computeTime())
!         write(*,*) "nametotvar", nametotvar(numid,:)
!         write(*,*) "I am dealing with ===", numid
!         write(*,*) "U_BOUND", UBOUND(D_data,DIM=3), UBOUND(tmp,DIM=3)
!         write(*,*) "==", l, computeTime(), Temperature, Number_of_2D_Fields
         do l=1,nlat
            do k=1,nlon
               tmp(k,l,:,1)=D_data(l,k,:)
            enddo
         enddo
      else
         status=NF_INQ_DIMID(ncid,'time',dimid)
         status=NF_INQ_DIMLEN(ncid,dimid,l)
         start(3) = computeTime()
         if ( l/=computeTime() ) call incretime(computeTime())
         allocate(tmp (nlon,nlat,computeTime(),1))
         do l=1,nlat
            tmp(:,l,1,1)=D_data(l,:,1)
         enddo
      end if

!     Output data for Variable
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,tab_dim,tmp)
      deallocate(tmp)

      IF (status/=nf_noerr) THEN
         write(error_id,*) nf_strerror(status)
         call ec_error(123)
      END IF

      end if

      RETURN
      END SUBROUTINE write

!     ----------------------------------------------------------------------------
!     output
!
!     This function checks whether the output flags of a field indicate that
!     the field should be written to any of the output files.

      FUNCTION output(flag) RESULT(result)

      REAL*8, INTENT(in) :: flag
      LOGICAL            :: result

      result = ( flag /= 0 )

      RETURN
      END FUNCTION output

!     ----------------------------------------------------------------------------
!     outputStdDev
!
!     This function checks whether the output flags of a field indicate that
!     the standard deviation should be written to any of the output files.

      FUNCTION outputStdDev(flag) RESULT(result)
      LOGICAL :: result
      INTEGER, INTENT(in) :: flag

      result = ( flag == 2 )

      RETURN
      END FUNCTION outputStdDev


!     ----------------------------------------------------------------------------
!     PRIVATE FUNCTIONS AND SUBROUTINES
!     ----------------------------------------------------------------------------

!     ----------------------------------------------------------------------------
!     initialize

      SUBROUTINE initializeModule

      use global_constants_mod, only: dblp=>dp, ip

      REAL*8, PARAMETER :: pi               = 2*acos(0.)
      REAL*8, PARAMETER :: radian_to_degree = 180.0/pi

      REAL*8, ALLOCATABLE, DIMENSION(:)   :: templon, templat
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: temp

      integer tab_dimid(5)
      integer tab_dim(5)

      INTEGER :: i,l,dimid,varid
      LOGICAL :: write_test = .FALSE.

      character(len=120):: longdata
      CHARACTER(len=uuid_size) :: uuid_obt

      integer,dimension(8) :: cvalues
      REAL*8 :: realmonth

      integer(kind=ip):: phi_lat_dat_id

      realmonth=1.0d0
      if ( irunlabeld /= 360 ) realmonth=(irunlabeld/30)+1.0d0

      IF (initialize_module) THEN
!     -----------------------------------------------------------------------
!     1. Setup two- and three-dimensional grids.

	     allocate(templon(nlon))
	     templon = (/ ((360.0d0*l)/nlon,l=0,nlon-1) /)
	     allocate(templat(nlat))
#if ( PMIP == 1 )
!mab
         open(newunit=phi_lat_dat_id,
     &    file='phi_lat.dat',form='unformatted',status='old')
         read(phi_lat_dat_id) (phi(i),i=1,nlat)
         close(phi_lat_dat_id)
!         print*, 'nach einlesen!!! phi',phi
#endif
	     templat = phi(1:nlat)*radian_to_degree

         allocate(lat_bounds(nlat+1))
         DO i = 2, nlat
           lat_bounds(i) = templat(i-1) - (templat(i-1) - templat(i))/2
         enddo
         lat_bounds(1) = templat(1) - (lat_bounds(3) -lat_bounds(2))/2.
         lat_bounds(nlat+1) = templat(nlat) +
     &                         (lat_bounds(nlat)-lat_bounds(nlat-1))/2.
         allocate(lon_bounds(nlon+1))
         DO i = 2, nlon
           lon_bounds(i) = templon(i-1) - (templon(i-1) - templon(i))/2
         enddo
         lon_bounds(1) = templon(1) - (lon_bounds(3) -lon_bounds(2))/2.
         lon_bounds(nlon+1) = templon(nlon) +
     &                         (lon_bounds(nlon)-lon_bounds(nlon-1))/2.
         allocate(time_bounds(2))
         time_bounds(1) = 0.5
         time_bounds(2) = 1.5

         CALL initdim(I_longitude,nlon,.FALSE.,.FALSE.,templon)

	     status=NF_DEF_DIM(ncid,"bnds",2,dimid)
	     tab_dimid(1)=dimid

	     allocate(temp(2,nlon))
	     DO i = 1, nlon
	        temp(1,i) = lon_bounds(i)
	        temp(2,i) = lon_bounds(i+1)
	     ENDDO
	     status=NF_INQ_DIMID(ncid,I_longitude%short,dimid)
         tab_dimid(2)=dimid
         longdata="lon_bounds"
	     call initbnds(longdata,nlon,tab_dimid,temp)
         deallocate(temp)
         deallocate(lon_bounds)

         CALL initdim(I_latitude,nlat,.FALSE.,.FALSE.,templat)

         allocate(temp(2,nlat))
         DO i = 1, nlat
            temp(1,i) = lat_bounds(i)
            temp(2,i) = lat_bounds(i+1)
         ENDDO
         status=NF_INQ_DIMID(ncid,I_latitude%short,dimid)
         tab_dimid(2)=dimid
         longdata="lat_bounds"
         call initbnds(longdata,nlat,tab_dimid,temp)
         deallocate(temp)
         deallocate(lat_bounds)

         IF ( thirdd(1)==1 )
     &   CALL initdim(I_pressure_T2,2,.FALSE.,.TRUE.,
     &                (/350.0d2, 650.0d2/))
         IF ( thirdd(2)==1 )
     &   CALL initdim(I_pressure_T3,3,.FALSE.,.TRUE.,
     &                (/100.0d2, 350.0d2, 650.0d2/))
         IF ( thirdd(3)==1 )
     &   CALL initdim(I_pressure_T4,4,.FALSE.,.TRUE.,
     &                (/100.0d2, 350.0d2, 650.0d2, 1000.0d2/))
         IF ( thirdd(4)==1 )
     &   CALL initdim(I_pressure_U3,3,.FALSE.,.TRUE.,
     &                (/200.0d2, 500.0d2, 800.0d2/))

!     -----------------------------------------------------------------------
!     2. Time coordinate axes

         SELECT CASE (how_to_compute_time)
         CASE (Compute_Time_in_Days)
            CALL initdim(I_time_in_days,NF_UNLIMITED,.TRUE.,.FALSE.,
     &                (/1.0d0/))
         CASE (Compute_Time_in_Years_Months)
            CALL initdim(I_time_in_months,NF_UNLIMITED,.TRUE.,.FALSE.,
     &                (/1.0d0/))
         CASE (Compute_Time_in_Months_Only)
            CALL initdim(I_time_in_months,NF_UNLIMITED,.TRUE.,.FALSE.,
     &                (/realmonth/))
         CASE (Compute_Time_in_Years_Only)
            CALL initdim(I_time_in_years,NF_UNLIMITED,.TRUE.,.FALSE.,
     &                (/1.0d0/))
         END SELECT

         allocate(temp(2,1))
         temp(1,1) = time_bounds(1)
         temp(2,1) = time_bounds(2)
         status=NF_INQ_DIMID(ncid,I_time_in_years%short,dimid)
         tab_dimid(2)=dimid
         longdata="time_bounds"
         call initbnds(longdata,1,tab_dimid,temp)
         deallocate(temp)
         deallocate(time_bounds)

         deallocate(templon)
         deallocate(templat)
         status=NF_ENDDEF(ncid)

!     -----------------------------------------------------------------------
!     3. Create fields.

         do i=1,numtotvar

         SELECT CASE (how_to_compute_time)
         CASE (Compute_Time_in_Days)
            if ( newtotvar(i,1)==1 ) write_test = .TRUE.
         CASE (Compute_Time_in_Years_Months)
            if ( newtotvar(i,3)==1 ) write_test = .TRUE.
         CASE (Compute_Time_in_Months_Only)
            if ( newtotvar(i,2)==1 ) write_test = .TRUE.
         CASE (Compute_Time_in_Years_Only)
            if ( newtotvar(i,4)==1 ) write_test = .TRUE.
         END SELECT

         if ( write_test ) then

         status=NF_REDEF(ncid)

         status=NF_INQ_DIMID(ncid,I_longitude%short,dimid)
         tab_dimid(1)=dimid
         tab_dim(1)=nlon
         status=NF_INQ_DIMID(ncid,I_latitude%short,dimid)
         tab_dimid(2)=dimid
         tab_dim(2)=nlat

         if ( nametotvar(i,5)/="N" ) then
            SELECT CASE (nametotvar(i,5))
            CASE ("T2")
               status=NF_INQ_DIMID(ncid,I_pressure_T2%short,dimid)
            CASE ("T3")
               status=NF_INQ_DIMID(ncid,I_pressure_T3%short,dimid)
            CASE ("T4")
               status=NF_INQ_DIMID(ncid,I_pressure_T4%short,dimid)
            CASE ("U3")
               status=NF_INQ_DIMID(ncid,I_pressure_U3%short,dimid)
            END SELECT
            status=NF_INQ_DIMLEN(ncid,dimid,l)
            tab_dimid(3)=dimid
            tab_dim(3)=l
            status=NF_INQ_DIMID(ncid,I_time_in_years%short,dimid)
            tab_dimid(4)=dimid
            tab_dim(4)=1
            CALL initvar(4,nametotvar(i,:),tab_dimid)
            status=NF_INQ_VARID(ncid,nametotvar(i,2),varid)
         else
            status=NF_INQ_DIMID(ncid,I_time_in_years%short,dimid)
            tab_dimid(3)=dimid
            tab_dim(3)=1
            CALL initvar(3,nametotvar(i,:),tab_dimid)
            status=NF_INQ_VARID(ncid,nametotvar(i,2),varid)
         end if

         status=NF_ENDDEF(ncid)
         end if

         write_test = .FALSE.

         enddo

!     Write the global attributes
         status=NF_REDEF(ncid)
         globalatt(6,1)="table_id"
         globalatt(6,2)="Table Amon ??? ref to be updated"
         globalatt(22,1)="tracking_id"
!~ !         call system('uuidgen > uuidgen_file')
!~          open(iuo+84,file = 'uuidgen_file',status='old',form='formatted')
!~          read(iuo+84,*) longdata
!~          close(iuo+84)
!~ !         call system('rm -rf uuidgen_file')
         uuid_obt = generate_uuid(4)
         globalatt(22,2)=trim(uuid_obt)
c~          globalatt(22,2)=trim(longdata)
         globalatt(24,1)="frequency"
         SELECT CASE (how_to_compute_time)
	      CASE (Compute_Time_in_Days)
	         globalatt(24,2)="day"
	      CASE (Compute_Time_in_Years_Months)
	         globalatt(24,2)="mon"
	      CASE (Compute_Time_in_Months_Only)
	         globalatt(24,2)="mon"
	      CASE (Compute_Time_in_Years_Only)
	         globalatt(24,2)="yr"
	     END SELECT
	     call date_and_time(VALUES=cvalues)
	     longdata=''
             write(longdata,999) cvalues(1), cvalues(2), cvalues(3), cvalues(5), 
     &                           cvalues(6), cvalues(7)
 999     format(i4.4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',i2.2,'Z')
         globalatt(25,1)="creation_date"
         globalatt(25,2)=trim(longdata)
         globalatt(26,1)="modeling_realm"
         globalatt(26,2)="atmos"
	     do i=1,26
                status=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,globalatt(i,1),
     &                  len_trim(globalatt(i,2)),trim(globalatt(i,2)))
	     enddo
         status=NF_ENDDEF(ncid)

!     -----------------------------------------------------------------------
!     6. Indicate that initialization has been done.

         initialize_module = .FALSE.

      END IF
      RETURN

      CONTAINS

      SUBROUTINE initdim(I_num,N_num,Is_time,Is_pressure,D_data)

      TYPE(info), INTENT(in)  :: I_num
      INTEGER,    INTENT(in)  :: N_num
      LOGICAL,    INTENT(in)  :: Is_time
      LOGICAL,    INTENT(in)  :: Is_pressure
      INTEGER                 :: dimid,varid
      INTEGER                 :: start(1)
      INTEGER                 :: count(1)
      REAL*8,     INTENT(in)  :: D_data(N_num)
      CHARACTER(len=256)      :: string,tmpstri,tmpstr
      INTEGER                 :: realmonth

      realmonth=1
      if ( irunlabeld /= 360 ) realmonth=(irunlabeld/30)+1

!     Definition
      status=NF_DEF_DIM(ncid,I_num%short,N_num,dimid)
      start(1)=dimid
      status=NF_DEF_VAR(ncid,I_num%short,NF_FLOAT,1,start,varid)

      IF (I_num%std=="time") THEN
         write(tmpstr,*) realmonth
         write(tmpstri,*) irunlabel
         SELECT CASE (I_num%unit)
         CASE ("years")
           string=              
     &  "years since 1-"//trim(adjustl(tmpstr))//"-"//trim(adjustl(tmpstri))
         CASE ("months")
           string=       
     &  "months since 1-"//trim(adjustl(tmpstr))//"-"//trim(adjustl(tmpstri))
         CASE ("days")
           string=       
     &  "days since 1-"//trim(adjustl(tmpstr))//"-"//trim(adjustl(tmpstri))     
         END SELECT
      END IF

!     Define the attributes
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_info%std,len_trim(I_num%std),trim(I_num%std))
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_info%long,len_trim(I_num%long),trim(I_num%long))
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_info%axis,len_trim(I_num%axis),trim(I_num%axis))
      IF (I_num%std=="time") THEN
         status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_info%unit,len_trim(string),trim(string))
      ELSE
         status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_info%unit,len_trim(I_num%unit),trim(I_num%unit))
      END IF

      IF (Is_pressure) status=NF_PUT_ATT_TEXT(ncid,varid,
     &      'positive',len_trim(I_num%more),trim(I_num%more))
      IF (Is_time) status=NF_PUT_ATT_TEXT(ncid,varid,
     &      'calendar',len_trim(I_num%more),trim(I_num%more))

      status=NF_ENDDEF(ncid)
!     Output data for Variable
      start(1) = 1
      count(1)   = N_num
c~       IF (Is_time) count(1)   = 1
c~  When using a time axis, N_num is NF_UNLIMITED = 0, hence D_data does not exist!
      IF (I_num%std=="time") then
        count(1)   = 1
        status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,(/1.0/))
      else
        status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,D_data)
      ENDIF
      status=NF_REDEF(ncid)

      IF (status/=nf_noerr) THEN
         write(error_id,*) nf_strerror(status)
         call ec_error(123)
      END IF

      RETURN
      END SUBROUTINE initdim

      SUBROUTINE initbnds(name,I_num,tab_dimid,D_data)

      character*120, INTENT(in) :: name
      INTEGER                   :: varid, start(2), count(2)
      INTEGER,      INTENT(in)  :: I_num, tab_dimid(2)

      REAL*8, DIMENSION(:,:), INTENT(in)  :: D_data

!     Definition
      status=NF_DEF_VAR(ncid,name,NF_DOUBLE,2,tab_dimid,varid)
      status=NF_ENDDEF(ncid)
!     Output data for Variable
      start(:) = 1
      count(1) = 2
      count(2) = I_num
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,D_data)
      status=NF_REDEF(ncid)

      IF (status/=nf_noerr) THEN
         write(error_id,*) nf_strerror(status)
         call ec_error(123)
      END IF

      RETURN
      END SUBROUTINE initbnds

      SUBROUTINE initvar(N_num,tab_c,tab_dimid)

      character*60, INTENT(in)  :: tab_c(5)
      INTEGER                   :: varid
      INTEGER,      INTENT(in)  :: tab_dimid(5)
      INTEGER,      INTENT(in)  :: N_num
      REAL*8                    :: tmp(1)

!     Definition
      status=NF_DEF_VAR(ncid,tab_c(2),NF_DOUBLE,N_num,tab_dimid,varid)

!     Define the attributes
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_infovar%std,len_trim(tab_c(3)),trim(tab_c(3)))
!      write(*,*) "name var to write: ", I_infovar%std
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_infovar%long,len_trim(tab_c(1)),trim(tab_c(1)))
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_infovar%unit,len_trim(tab_c(4)),trim(tab_c(4)))

      tmp(1)=fill_value
      status=NF_PUT_ATT_DOUBLE(ncid,varid,
     &      I_infovar%fill,NF_DOUBLE,1,tmp)
      tmp(1)=missing_value
      status=NF_PUT_ATT_DOUBLE(ncid,varid,
     &      I_infovar%miss,NF_DOUBLE,1,tmp)

      IF (status/=nf_noerr) THEN
         write(error_id,*) nf_strerror(status)
         call ec_error(123)
      END IF

      RETURN
      END SUBROUTINE initvar

      END SUBROUTINE initializeModule

!     ----------------------------------------------------------------------------
!     computeTime
!
!     This function computes the current time in days or in months depending
!     on the type of data that is being written to the output file.

      FUNCTION computeTime() RESULT(result)
      INTEGER :: result
      INTEGER :: realyear,realmonth

      realyear=iyear
      realmonth=0
      if ( irunlabeld == 360 ) realyear=iyear-1
      if ( irunlabeld /= 360 ) realmonth=(irunlabeld/30)

      SELECT CASE (how_to_compute_time)
      CASE (Compute_Time_in_Days)
         result = realyear*360.0 + (imonth-1)*30.0 + iday
      CASE (Compute_Time_in_Months_Only)
         !result = realyear*12.0 + imonth - realmonth
! --- BdB 05-2019: adjust time, using the resetted year counter
         result = new_year_atm*12.0 + imonth - realmonth
      CASE (Compute_Time_in_Years_Months)
         result = imonth
      CASE (Compute_Time_in_Years_Only)
         !result = iyear
! --- BdB 05-2019: adjust time, using the resetted year counter
         result=new_year_atm
      END SELECT

      RETURN
      END FUNCTION computeTime

      END MODULE Atmosphere_Output
