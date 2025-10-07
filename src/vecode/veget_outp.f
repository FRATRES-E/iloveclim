!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:34 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:34 CET 2009

      MODULE Vegetation_Output

#if ( COMATM == 1 )
      use comdiag, only:
      use comemic_mod, only: iyear, imonth, fini
      use comatm, only: nlat, nlon
      use veget_mod
#endif
! --- BdB 05-2019: added variables for writing time
      use comemic_mod, only: irunlabel, nwrskip, new_year_veg
     &               , time_in_years, current_int_veg


      IMPLICIT NONE

      include 'netcdf.inc'

!     include 'comemic.h'  !mohr commenting these lines in original
!     include 'comdiag.h'  !to avoid double definition of variables

#if ( COMATM == 0 )
#include "comdiag.h"
#include "comemic.h"
#include "veget.h"
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
      PUBLIC :: Monthly_Means
      PUBLIC :: Yearly_Means

!     Codes for one dimensional fields
      REAL*8, ALLOCATABLE, DIMENSION(:) :: lon_bounds, lat_bounds, time_bounds

! --- Added global variables for newunit construction
      integer, public :: vegetzav_id, vegetoutp_id

!     Codes for two dimensional fields
      PUBLIC :: tree_fraction
      PUBLIC :: grass_fraction
      PUBLIC :: desert_fraction
      PUBLIC :: needle_tree_fraction
      PUBLIC :: tree_leaf_area_index
      PUBLIC :: grass_leaf_area_index
      PUBLIC :: vegetation_albedo
      PUBLIC :: surface_temperature
      PUBLIC :: Annual_gdd0_index
      PUBLIC :: precipitation
      PUBLIC :: reduced_precipitaion
      PUBLIC :: Biomass_carbon_uptake
      PUBLIC :: Biomass_carbon_stock
      PUBLIC :: leaves_biomass_carbon_stock
      PUBLIC :: stems_root_biomass_stock
      PUBLIC :: litter_carbon_stock
      PUBLIC :: mortmass_and_soil_organic_matter_carbon_stock
      PUBLIC :: net_primary_productivity
      PUBLIC :: ice_sheet_fraction

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
     & I_time_in_months = Info("time",  "time",          "Model Time",
     &          "months since 1850-1-1", "T",    "360_day"),
     & I_time_in_years  = Info("time",  "time",          "Model Time",
     &          "years since 1850-1-1",  "T",    "360_day")
!     Grid coordinates
      TYPE(Infovar), PARAMETER ::
     & I_infovar        = Infovar("name",  "standard_name", "long_name",
     &          "units",                 "_FillValue", "missing_value")

!     One dimensional fields
!      INTEGER, PARAMETER :: lon  = 1
!      INTEGER, PARAMETER :: lat  = 2
!      INTEGER, PARAMETER :: time = 7

!     Two dimensional fields
      INTEGER, PARAMETER :: tree_fraction               =  1
      INTEGER, PARAMETER :: grass_fraction              =  2
      INTEGER, PARAMETER :: desert_fraction             =  3
      INTEGER, PARAMETER :: needle_tree_fraction        =  4
      INTEGER, PARAMETER :: tree_leaf_area_index        =  5
      INTEGER, PARAMETER :: grass_leaf_area_index       =  6
      INTEGER, PARAMETER :: vegetation_albedo           =  7
      INTEGER, PARAMETER :: surface_temperature         =  8
      INTEGER, PARAMETER :: Annual_gdd0_index           =  9
      INTEGER, PARAMETER :: precipitation               =  10
      INTEGER, PARAMETER :: reduced_precipitaion        =  11
      INTEGER, PARAMETER :: Biomass_carbon_uptake       =  12
      INTEGER, PARAMETER :: Biomass_carbon_stock        =  13
      INTEGER, PARAMETER :: leaves_biomass_carbon_stock =  14
      INTEGER, PARAMETER :: stems_root_biomass_stock    =  15
      INTEGER, PARAMETER :: litter_carbon_stock         =  16
      INTEGER, PARAMETER ::
     &    mortmass_and_soil_organic_matter_carbon_stock =  17
      INTEGER, PARAMETER :: net_primary_productivity    =  18
      INTEGER, PARAMETER :: ice_sheet_fraction          =  19

      INTEGER, PARAMETER :: Compute_Time_in_Months_Only  = 1
      INTEGER, PARAMETER :: Compute_Time_in_Years_Only   = 2

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

! --- BdB 05-2019: update name every nwrskip (the restart output interval)
      CHARACTER*6, SAVE  :: finichf

      ! initialise 
      if (iyear==1) finichf=fini
      ! update
      if (mod(iyear-1,nwrskip).eq.0) then
       write(finichf,1) irunlabel+iyear-1
 1     format(i6.6)
      end if

! --- BdB 05-2019: update output names
      SELECT CASE (file)
      CASE (Monthly_Means)
         output_filename = "outputdata/vegetation/vegmmyl"//finichf//".nc"
         how_to_compute_time = Compute_Time_in_Months_Only
      CASE (Yearly_Means)
         output_filename = "outputdata/vegetation/vegym"//finichf//".nc"
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
         write(iveg+29,*) nf_strerror(status)
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
         write(iveg+29,*) nf_strerror(status)
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
         CASE (Compute_Time_in_Months_Only)
          ! determine year
          D_data(i) = time_in_years(current_int_veg) + real(i) / 12.
         CASE (Compute_Time_in_Years_Only)
          D_data(i) = time_in_years(current_int_veg+i)
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
         CASE (Compute_Time_in_Months_Only)
          ! determine year
          D_data2(1,i) = time_in_years(current_int_veg) + (real(i) - 0.5) / 12.
          D_data2(2,i) = time_in_years(current_int_veg) + (real(i) + 0.5) / 12.
         CASE (Compute_Time_in_Years_Only)
          D_data2(1,i) = time_in_years(current_int_veg+i) - 0.5
          D_data2(2,i) = time_in_years(current_int_veg+i) + 0.5
         END SELECT
      enddo
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,D_data2)

      IF (status/=nf_noerr) THEN
         write(iveg+29,*) nf_strerror(status)
         call ec_error(123)
      END IF

      RETURN
      END SUBROUTINE incretime

      SUBROUTINE write(numid,D_data)

      INTEGER, INTENT(in) :: numid
      INTEGER             :: tab_dim(5)
      INTEGER             :: start(5)
      INTEGER             :: l,varid,dimid
      LOGICAL             :: write_test

      REAL*8, DIMENSION(:,:,:), allocatable :: tmp
      REAL*8, DIMENSION(:,:), INTENT(in)    :: D_data

      write_test = .FALSE.
      start(:)   = 1
      tab_dim(:) = 1
      SELECT CASE (how_to_compute_time)
      CASE (Compute_Time_in_Months_Only)
         if ( newvegvar(numid,3)==1 ) write_test = .TRUE.
      CASE (Compute_Time_in_Years_Only)
         if ( newvegvar(numid,4)==1 ) write_test = .TRUE.
      END SELECT

      if ( write_test ) then

      status=NF_INQ_VARID(ncid,namevegvar(numid,2),varid)
      tab_dim(1)=nlon
      tab_dim(2)=nlat

      status=NF_INQ_DIMID(ncid,'time',dimid)
      status=NF_INQ_DIMLEN(ncid,dimid,l)
      start(3) = computeTime()
      if ( l/=computeTime() ) call incretime(computeTime())
      allocate(tmp (nlon,nlat,computeTime()))
      do l=1,nlat
         tmp(:,l,1)=D_data(l,:)
      enddo

!     Output data for Variable
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,tab_dim,tmp)
      deallocate(tmp)

      IF (status/=nf_noerr) THEN
         write(iveg+29,*) nf_strerror(status)
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

      REAL*8, PARAMETER :: pi               = 2*acos(0.)
      REAL*8, PARAMETER :: radian_to_degree = 180.0/pi

      REAL*8, ALLOCATABLE, DIMENSION(:)   :: templon, templat
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: temp

      integer tab_dimid(5)
      integer tab_dim(5)

      INTEGER :: i,l,dimid,varid
      LOGICAL :: write_test = .FALSE.

      character(len=120):: longdata
      integer,dimension(8) :: cvalues

      IF (initialize_module) THEN
!     -----------------------------------------------------------------------
!     1. Setup two- and three-dimensional grids.

	     allocate(templon(nlon))
	     templon = (/ ((360.0d0*l)/nlon,l=0,nlon-1) /)
	     allocate(templat(nlat))
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

!     -----------------------------------------------------------------------
!     2. Time coordinate axes

         SELECT CASE (how_to_compute_time)
         CASE (Compute_Time_in_Months_Only)
            CALL initdim(I_time_in_months,NF_UNLIMITED,.TRUE.,.FALSE.,
     &                (/1.0d0/))
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

         do i=1,numvegvar

         SELECT CASE (how_to_compute_time)
         CASE (Compute_Time_in_Months_Only)
            if ( newvegvar(i,3)==1 ) write_test = .TRUE.
         CASE (Compute_Time_in_Years_Only)
            if ( newvegvar(i,4)==1 ) write_test = .TRUE.
         END SELECT

         if ( write_test ) then

         status=NF_REDEF(ncid)

         status=NF_INQ_DIMID(ncid,I_longitude%short,dimid)
         tab_dimid(1)=dimid
         tab_dim(1)=nlon
         status=NF_INQ_DIMID(ncid,I_latitude%short,dimid)
         tab_dimid(2)=dimid
         tab_dim(2)=nlat

         status=NF_INQ_DIMID(ncid,I_time_in_years%short,dimid)
         tab_dimid(3)=dimid
         tab_dim(3)=1
         CALL initvar(3,namevegvar(i,:),tab_dimid)
         status=NF_INQ_VARID(ncid,namevegvar(i,2),varid)

         status=NF_ENDDEF(ncid)
         end if

         write_test = .FALSE.

         enddo

!     Write the global attributes, here is where the code crashes!!! !mohr
         status=NF_REDEF(ncid)
c~          globalatt(6,1)="table_id"
c~          globalatt(6,2)="Table Lmon (13 November 2004)"
c~          globalatt(22,1)="tracking_id"
c~          call system('uuidgen > uuidgen_file')
c~          open(iveg+84,file = 'uuidgen_file',status='old',form='formatted')
c~          read(iveg+84,*) longdata
c~          close(iveg+84)
c~          call system('rm -rf uuidgen_file')
c~          globalatt(22,2)=trim(longdata)
c~          globalatt(24,1)="frequency"
c~          globalatt(24,2)="yr"
	     call date_and_time(VALUES=cvalues)
	     longdata=''
	     write(longdata,999) cvalues(1), cvalues(2), cvalues(3), cvalues(5),
     &    cvalues(6), cvalues(7)
 999     format(i4.4,'-',i2.2,'-',i2.2,'T',i2.2,':',i2.2,':',i2.2,'Z')
c~          globalatt(25,1)="creation_date"
c~          globalatt(25,2)=trim(longdata)
c~          globalatt(26,1)="modeling_realm"
c~          globalatt(26,2)="land"
	     do i=1,26
	       !status=NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,globalatt(i,1),len_trim(globalatt(i,2)),trim(globalatt(i,2)))
               !mohr: I had to comment this line, other way the writing in the netCDF was crashing.
               !mohr: I should check why but, by the moment, it works.
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
      CHARACTER(len=256)      :: string

!     Definition
      status=NF_DEF_DIM(ncid,I_num%short,N_num,dimid)
      start(1)=dimid
      status=NF_DEF_VAR(ncid,I_num%short,NF_FLOAT,1,start,varid)

!      IF (I_info%std=="time") THEN
!         string="months since "//globalatt(19,2)//"-1-1"
!      END IF

!     Define the attributes
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_info%std,len_trim(I_num%std),trim(I_num%std))
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_info%long,len_trim(I_num%long),trim(I_num%long))
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_info%axis,len_trim(I_num%axis),trim(I_num%axis))
      IF (I_info%std=="time") THEN
         status=NF_PUT_ATT_TEXT(ncid,varid,
     &      string,len_trim(I_num%unit),trim(I_num%unit))
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
      IF (Is_time) count(1)   = 1
      status=NF_PUT_VARA_DOUBLE(ncid,varid,start,count,D_data)
      status=NF_REDEF(ncid)

      IF (status/=nf_noerr) THEN
         write(iveg+29,*) nf_strerror(status)
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
         write(iveg+29,*) nf_strerror(status)
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
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_infovar%long,len_trim(tab_c(1)),trim(tab_c(1)))
      status=NF_PUT_ATT_TEXT(ncid,varid,
     &      I_infovar%unit,len_trim(tab_c(4)),trim(tab_c(4)))

      tmp(1)=veg_fill_value
      status=NF_PUT_ATT_DOUBLE(ncid,varid,
     &      I_infovar%fill,NF_DOUBLE,1,tmp)
      tmp(1)=veg_missing_value
      status=NF_PUT_ATT_DOUBLE(ncid,varid,
     &      I_infovar%miss,NF_DOUBLE,1,tmp)

      IF (status/=nf_noerr) THEN
         write(iveg+29,*) nf_strerror(status)
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

      SELECT CASE (how_to_compute_time)
      CASE (Compute_Time_in_Months_Only)
         !result = imonth
! ---    BdB 05-2019: adjust time, using the resetted year counter
         result = new_year_veg*12.0 + imonth
      CASE (Compute_Time_in_Years_Only)
         !result = iyear
! ---    BdB 05-2019: adjust time, using the resetted year counter
         result=new_year_veg
      END SELECT

      RETURN
      END FUNCTION computeTime

      END MODULE Vegetation_Output
