!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2020-2021 Didier M. Roche (a.k.a. dmr)

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#include "choixcomposantes.h"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      MODULE IO_NC_MOD

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~         USE uuid_fort_wrap, ONLY: uuid_size
        USE global_constants_mod, ONLY: dblp=>dp, silp=>sp, ip, str_len, uuid_size

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        IMPLICIT NONE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   History
! dmr           Change from 0.0.0: created the actual grid object and set the initialization
! dmr           Change from 0.1.0: added objects axes and file
! dmr           Change from 0.1.5: created the file initialization function
! dmr           Change from 0.2.0: added support for 3D file without time, some cleaning
! dmr           Change from 0.3.0: added support for 3D file with    time, some cleaning
! dmr           Change from 0.4.0: always compiled (no longer CLIO_OUT_NEWGEN only), geolocation lookup removed,
! dmr                              filename no longer limited to 18 characters, added exact (raw netCDF) readers:
! dmr                              IO_NC_FILE%open (existing file), IO_NC_FILE%nrec, IO_GRID_VAR%read (2D + time)
! dmr&clo       Change from 0.5.0: CF-1.8 metadata. Global attributes title, institution, source, author, history,
! dmr&clo                          Conventions (lower case, CF), institution/author/source settable via io_nc_set_metadata.
! dmr&clo                          Variables: long_name, standard_name, units, optional cell_methods, no dangling
! dmr&clo                          grid_mapping. Axes: long_name, standard_name, axis, optional bounds.
! dmr&clo                          Time writers take an optional explicit time value (default: record index).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        CHARACTER(LEN=5), PARAMETER :: version_mod ="0.6.0"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER,          PARAMETER :: axisNmSize = 6
        CHARACTER(LEN=7), PARAMETER :: undef_str  = "NotSet"
        REAL(KIND=dblp),  PARAMETER :: undef_dblp = (-1.0_dblp)*HUGE(1.0_silp)
        REAL(KIND=silp),  PARAMETER :: undef_silp = (-1.0_silp)*HUGE(1.0_silp)
        INTEGER(KIND=ip), PARAMETER :: undef_ip   = (-1_ip)*HUGE(1_ip)
        INTEGER(KIND=ip), PARAMETER :: time_chunk = 10

! dmr&clo   CF version written in the global attribute Conventions.
        CHARACTER(LEN=*), PARAMETER :: CF_CONVENTIONS = "CF-1.8"

! dmr&clo   Run-wide global attributes, set once by io_nc_set_metadata (from the global namelist group ncmeta).
        CHARACTER(LEN=str_len), SAVE :: md_institution = undef_str
        CHARACTER(LEN=str_len), SAVE :: md_author      = undef_str
        CHARACTER(LEN=str_len), SAVE :: md_source      = undef_str

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  dmr   Holder type for the axes variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        TYPE :: IO_NC_AXIS

          PRIVATE
          LOGICAL                   :: initialized =.false.
          CHARACTER(LEN=uuid_size)  :: uuid_axis
          CHARACTER(LEN=axisNmSize) :: axis_name   = undef_str
          CHARACTER(LEN=str_len)    :: axis_unit   = undef_str
          INTEGER(kind=ip)          :: axis_size
          LOGICAL                   :: is_time     = .false.
          CHARACTER(LEN=str_len)    :: calendar    = undef_str
          CHARACTER(LEN=str_len)    :: long_name   = undef_str
          CHARACTER(LEN=str_len)    :: std_name    = undef_str
          CHARACTER(LEN=str_len)    :: axis_attr   = undef_str          ! CF axis attribute: X, Y, Z or T
          REAL(kind=dblp), dimension(:,:), allocatable :: axis_bounds   ! (2, axis_size), written as <axis>_bnds
          REAL(kind=dblp), dimension(:)  , allocatable :: axis_array_rdblp
          REAL(kind=silp), dimension(:)  , allocatable :: axis_array_rsilp
          INTEGER(kind=ip), dimension(:) , allocatable :: axis_array_isilp


          CONTAINS
            PROCEDURE, PUBLIC  :: init => init_io_nc_axis
            PROCEDURE, PUBLIC  :: show => show_io_nc_axis
            PROCEDURE, PUBLIC  :: wrte => wrte_io_nc_axis

        END TYPE IO_NC_AXIS

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  dmr   Holder type for the actual netCDF files
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        TYPE :: IO_NC_FILE

          PRIVATE
          LOGICAL                 :: initialized       =.false.
          CHARACTER(LEN=uuid_size):: uuid_file
          CHARACTER(LEN=str_len) :: filename
          CHARACTER(LEN=str_len)  :: institution_name  = undef_str
          CHARACTER(LEN=str_len)  :: short_author_name = undef_str
          CHARACTER(LEN=str_len)  :: title_file        = "Output generated by io_nc v."//version_mod
          LOGICAL                 :: overwrte_s        = .true.
          LOGICAL                 :: netcdf4f_s        = .true.

          CONTAINS
            PROCEDURE, PUBLIC :: init => init_io_nc_file
            PROCEDURE, PUBLIC :: open => open_io_nc_file
            PROCEDURE, PUBLIC :: nrec => nrec_io_nc_file

        END TYPE IO_NC_FILE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  dmr   Holder type for the metada and control data
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        TYPE :: IO_GRID_VAR

          PRIVATE
          LOGICAL                  :: initialized=.false.
          CHARACTER(LEN=uuid_size) :: uuid_var
          CHARACTER(LEN=str_len)   :: VarName
!~        VarSource
!~        MODSource
          CHARACTER(LEN=str_len)   :: Long_Name= undef_str
          CHARACTER(LEN=str_len)   :: STD_Name = undef_str
          CHARACTER(LEN=str_len)   :: Units    = undef_str
          CHARACTER(LEN=str_len)   :: Cell_Methods = undef_str

! dmr The size of this variable should be calculated on the longest name requested automatically
          CHARACTER(LEN=7), DIMENSION(:), ALLOCATABLE :: Axes_List
          INTEGER                  :: nbaxes=0
          TYPE(IO_NC_FILE), POINTER:: nc_file

          CONTAINS

            PROCEDURE, PUBLIC  :: init => init_io_grid_var
            PROCEDURE, PUBLIC  :: show => show_io_grid_var
            PROCEDURE, PUBLIC  :: setf => set_NCfilename_io_grid_var
            GENERIC  , PUBLIC  :: wrte => write_io_grid_2Dvar_TIME, write_io_grid_3Dvar_noTIME, write_io_grid_3Dvar_TIME
            PROCEDURE, PRIVATE :: write_io_grid_2Dvar_TIME
            PROCEDURE, PRIVATE :: write_io_grid_3Dvar_noTIME
            PROCEDURE, PRIVATE :: write_io_grid_3Dvar_TIME
            GENERIC  , PUBLIC  :: read => read_io_grid_2Dvar_TIME
            PROCEDURE, PRIVATE :: read_io_grid_2Dvar_TIME

        END TYPE IO_GRID_VAR


       CONTAINS


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE init_io_nc_axis(this, AxisName, OPT_vals_to_write1D, OPT_AxisUnit, OPT_isTime, OPT_calendar,            &
                                  OPT_LongName, OPT_StdName, OPT_Axis, OPT_bounds)

!~        USE uuid_fort_wrap, ONLY: UUID_V4
       use uuid_module, only: generate_uuid

       CLASS(IO_NC_AXIS)               , INTENT(OUT):: this
       CHARACTER(LEN=*)                , INTENT(IN) :: AxisName
       CHARACTER(LEN=*),       OPTIONAL, INTENT(IN) :: OPT_AxisUnit
       LOGICAL,                OPTIONAL, INTENT(IN) :: OPT_isTime
       CHARACTER(LEN=*),       OPTIONAL, INTENT(IN) :: OPT_calendar
       CHARACTER(LEN=*),       OPTIONAL, INTENT(IN) :: OPT_LongName
       CHARACTER(LEN=*),       OPTIONAL, INTENT(IN) :: OPT_StdName
       CHARACTER(LEN=*),       OPTIONAL, INTENT(IN) :: OPT_Axis
       REAL(kind=dblp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: OPT_bounds   ! (2, axis size)

       CLASS(*), DIMENSION(:), OPTIONAL, intent(in)     :: OPT_vals_to_write1D


         IF ( .NOT. this%initialized ) then

!~            CALL UUID_V4(this%uuid_axis)
           this%uuid_axis = generate_uuid(4)
           
           this%axis_name = AxisName

           IF (PRESENT(OPT_isTime)) THEN
              this%is_time = OPT_isTime
           ENDIF

           IF (this%is_time) THEN
                ! DO THE TIME AXIS
                this%axis_size = time_chunk
                ALLOCATE(this%axis_array_rsilp(this%axis_size))
                this%axis_array_rsilp(:) = undef_silp
           ELSE
             this%axis_size = SIZE(OPT_vals_to_write1D,DIM=1)
             SELECT TYPE(OPT_vals_to_write1D)

               TYPE IS (REAL(kind=dblp))
                  ALLOCATE(this%axis_array_rdblp(this%axis_size))
                  this%axis_array_rdblp(:) = OPT_vals_to_write1D(:)

               TYPE IS (REAL(kind=silp))
                  ALLOCATE(this%axis_array_rsilp(this%axis_size))
                  this%axis_array_rsilp= OPT_vals_to_write1D(:)

               TYPE IS (INTEGER(kind=silp))
                  ALLOCATE(this%axis_array_isilp(this%axis_size))
                  this%axis_array_isilp= OPT_vals_to_write1D(:)

               CLASS DEFAULT
                 WRITE(*,*) "UNKNOWN TYPE IN init_io_nc_axis"

                END SELECT

           ENDIF

           IF ( PRESENT(OPT_AxisUnit) ) THEN
              this%axis_unit = OPT_AxisUnit
           ENDIF

           IF ( PRESENT(OPT_calendar) ) THEN
              this%calendar = OPT_calendar
           ENDIF

           IF ( PRESENT(OPT_LongName) ) this%long_name = OPT_LongName
           IF ( PRESENT(OPT_StdName) )  this%std_name  = OPT_StdName
           IF ( PRESENT(OPT_Axis) )     this%axis_attr = OPT_Axis
           IF ( PRESENT(OPT_bounds) ) THEN
              ALLOCATE(this%axis_bounds(2,SIZE(OPT_bounds,DIM=2)))
              this%axis_bounds(:,:) = OPT_bounds(:,:)
           ENDIF

         ELSE ! I am called to initialize a nc_file that already exist ... bizarre!

            WRITE(*,*) "un-implemented [ABORT]"
            STOP 1
        ENDIF

        this%initialized = .TRUE.

       END SUBROUTINE init_io_nc_axis

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE init_io_nc_file(this, FileName, OPT_TitleFile, OPT_Institution, OPT_overwrite, OPT_netCDF4)

!~        USE uuid_fort_wrap, ONLY: UUID_V4
       use uuid_module, only: generate_uuid
       
       USE ncio, ONLY: nc_create, nc_write_attr, nc_write_dim, nc_write

       CLASS(IO_NC_FILE)         , INTENT(OUT):: this
       CHARACTER(LEN=*)          , INTENT(IN) :: FileName
       CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: OPT_TitleFile
       CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: OPT_Institution
       LOGICAL,          OPTIONAL, INTENT(IN) :: OPT_overwrite
       LOGICAL,          OPTIONAL, INTENT(IN) :: OPT_netCDF4


       INTEGER(kind=ip)       :: length, rc
       CHARACTER(LEN=str_len) :: logname
       CHARACTER(LEN=8)       :: cdate
       CHARACTER(LEN=10)      :: ctime


         IF ( .NOT. this%initialized ) then

!~            CALL UUID_V4(this%uuid_file)
           this%uuid_file = generate_uuid(4)
           this%filename = FileName


           IF (md_author /= undef_str) THEN
              this%short_author_name = md_author
           ELSE
              CALL GET_ENVIRONMENT_VARIABLE('LOGNAME', logname, length, rc)
              IF (rc == 0) THEN
                 this%short_author_name = TRIM(logname)
              ENDIF
           ENDIF

           this%institution_name = md_institution

           IF ( PRESENT(OPT_TitleFile) ) THEN
              this%title_file = OPT_TitleFile
           ENDIF

           IF ( PRESENT(OPT_Institution) ) THEN
              this%institution_name = OPT_Institution
           ENDIF

           IF ( PRESENT(OPT_overwrite) ) THEN
              this%overwrte_s = OPT_overwrite
           ENDIF

           IF ( PRESENT(OPT_netCDF4) ) THEN
              this%netcdf4f_s = OPT_netCDF4
           ENDIF

           CALL DATE_AND_TIME(DATE=cdate, TIME=ctime)

           CALL nc_create(this%filename,overwrite=this%overwrte_s,netcdf4=this%netcdf4f_s,                                &
                          author=TRIM(this%short_author_name), institution=TRIM(this%institution_name))
           CALL nc_write_attr(this%filename,"title",TRIM(this%title_file))
           CALL nc_write_attr(this%filename,"source",TRIM(md_source))
           CALL nc_write_attr(this%filename,"Conventions",CF_CONVENTIONS)
           CALL nc_write_attr(this%filename,"history",cdate(1:4)//"-"//cdate(5:6)//"-"//cdate(7:8)//"T"//ctime(1:2)//":"//  &
                              ctime(3:4)//":"//ctime(5:6)//" created by iLOVECLIM (io_nc v"//version_mod//")")

         ELSE ! I am called to initialize a nc_file that already exist ... bizarre!

            WRITE(*,*) "un-implemented [ABORT]"
            STOP 1
        ENDIF

        this%initialized = .TRUE.

       END SUBROUTINE init_io_nc_file

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      SUBROUTINE set_NCfilename_io_grid_var(this,nc_file)

        CLASS(IO_GRID_VAR),         INTENT(INOUT) :: this
        TYPE(IO_NC_FILE), TARGET,   INTENT(IN)  :: nc_file


            ! nullify(this%nc_file)
            this%nc_file => nc_file

      END SUBROUTINE set_NCfilename_io_grid_var
! ---
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      SUBROUTINE init_io_grid_var(this,varname,nc_file,axeslist,OPT_longname,OPT_stdname,OPT_units,OPT_cellmethods)

!~         USE uuid_fort_wrap, ONLY: UUID_V4
        use uuid_module, only: generate_uuid
       
        CLASS(IO_GRID_VAR),         INTENT(OUT) :: this
        CHARACTER(LEN=*),           INTENT(IN)  :: varname
        TYPE(IO_NC_FILE), TARGET,   INTENT(IN)  :: nc_file
        CHARACTER(LEN=*),           INTENT(IN)  :: axeslist
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: OPT_longname
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: OPT_stdname
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: OPT_units
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: OPT_cellmethods

        INTEGER :: i

        IF ( .NOT. this%initialized ) then

            this%VarName = varname

            WRITE(*,*) "DEBUG ===", this%nbaxes, (axeslist(i:i), i=1, len_trim(axeslist)-1)
            this%nbaxes = count( (/ (axeslist(i:i), i=1, len_trim(axeslist)) /) == " ")+1

            ALLOCATE(this%Axes_List(this%nbaxes))
            read(axeslist,fmt=*) (this%Axes_List(i),i=1,this%nbaxes)

!~             CALL UUID_V4(this%uuid_var)
            this%uuid_var = generate_uuid(4)
            this%nc_file => nc_file

            IF ( PRESENT(OPT_longname) ) then
               this%long_name = OPT_longname
            ENDIF

            IF ( PRESENT(OPT_stdname) ) then
               this%STD_Name = OPT_stdname
            ENDIF

            IF ( PRESENT(OPT_units) ) then
               this%Units = OPT_units
            ENDIF

            IF ( PRESENT(OPT_cellmethods) ) then
               this%Cell_Methods = OPT_cellmethods
            ENDIF


        ELSE ! I am called to initialize a variable that already exist ... bizarre!
            WRITE(*,*) "un-implemented [ABORT]"
            STOP 1
        ENDIF

        this%initialized = .TRUE.

      END SUBROUTINE init_io_grid_var

! ---

      SUBROUTINE show_io_grid_var(this)

        CLASS(IO_GRID_VAR),       INTENT(IN) :: this

        WRITE(*,*)
        WRITE(*,*) "========================================"
        WRITE(*,*) "INFORMATION FOR /clio_grid_var/ variable"
        WRITE(*,*) "VarName = ", trim(this%VarName)
        WRITE(*,*) "uuidVar = ", this%uuid_var
        WRITE(*,*) "nbaxes  = ", this%nbaxes
        WRITE(*,*) "AxesList= ", this%Axes_List(:)
        WRITE(*,*) "LongName= ", trim(this%Long_Name)
        WRITE(*,*) "STD_Name= ", trim(this%STD_Name)
        WRITE(*,*) "Units   = ", trim(this%Units)
        WRITE(*,*) "FileNC  = ", trim(this%nc_file%filename)
        WRITE(*,*) "========================================"

      END SUBROUTINE show_io_grid_var

! ---

      SUBROUTINE show_io_nc_axis(this)

        CLASS(IO_NC_AXIS),       INTENT(IN) :: this

        WRITE(*,*)
        WRITE(*,*) "========================================"
        WRITE(*,*) "INFORMATION FOR /io_nc_axis/ Axis"
        WRITE(*,*) "AxisName = ", trim(this%axis_name)
        WRITE(*,*) "uuidAxis = ", this%uuid_axis
        WRITE(*,*) "AxesSize = ", this%axis_size
        WRITE(*,*) "Units   = ", trim(this%axis_unit)
        WRITE(*,*) "========================================"

      END SUBROUTINE show_io_nc_axis

! ---

      SUBROUTINE write_io_grid_2Dvar_TIME(this, values_to_write, pass_nc, OPT_TimeValue)

        USE ncio, only: nc_write, nc_write_attr

        CLASS(IO_GRID_VAR),       INTENT(IN) :: this

        INTEGER, INTENT(IN)                  :: pass_nc
        REAL(kind=dblp), OPTIONAL, INTENT(IN):: OPT_TimeValue   ! value of the time coordinate (default: pass_nc)

        CLASS(*), DIMENSION(:,:), intent(in) :: values_to_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (present(OPT_TimeValue)) then
           call nc_write(this%nc_file%filename,"time",OPT_TimeValue,dim1="time",start=[pass_nc],count=[1])
        else
           call nc_write(this%nc_file%filename,"time",pass_nc,dim1="time",start=[pass_nc],count=[1])
        endif

        select type(values_to_write)

           type is (REAL(kind=dblp))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:)                                        &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3)                           &
                          ,start=[1,1,pass_nc],count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2),1],      &
                            long_name=att(this%Long_Name), standard_name=att(this%STD_Name), units=att(this%Units),  &
                            grid_mapping="", missing_value=undef_dblp)

           type is (INTEGER(kind=ip))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:)                                        &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3)                           &
                          ,start=[1,1,pass_nc],count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2),1],      &
                            long_name=att(this%Long_Name), standard_name=att(this%STD_Name), units=att(this%Units),  &
                            grid_mapping="", missing_value=undef_ip)

           type is (REAL(KIND=silp))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:)                                        &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3)                           &
                          ,start=[1,1,pass_nc],count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2),1],      &
                            long_name=att(this%Long_Name), standard_name=att(this%STD_Name), units=att(this%Units),  &
                            grid_mapping="", missing_value=undef_silp)
            CLASS DEFAULT
              WRITE(*,*) "UNKNOWN TYPE IN io_nc_mod / nc_write"

         END SELECT

        if (this%Cell_Methods /= undef_str)                                                                                 &
           call nc_write_attr(this%nc_file%filename, this%VarName, "cell_methods", trim(this%Cell_Methods))


      END SUBROUTINE write_io_grid_2Dvar_TIME

! ---

      SUBROUTINE write_io_grid_3Dvar_noTIME(this, values_to_write)

        USE ncio, only: nc_write, nc_write_attr

        CLASS(IO_GRID_VAR),         INTENT(IN) :: this

        CLASS(*), DIMENSION(:,:,:), intent(in) :: values_to_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


        select type(values_to_write)

           type is (REAL(kind=dblp))
!~               write(*,*) "wrte tbl_wrte ....", __LINE__, __FILE__ ,                                                  &
!~                           UBOUND(values_to_write,dim=1), UBOUND(values_to_write,dim=2),UBOUND(values_to_write,dim=3)
!~               write(*,*) this%nc_file%filename,this%VarName, this%Axes_List(1), this%Axes_List(2), this%Axes_List(3)

              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:,:)                                      &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3)                           &
                          , long_name=att(this%Long_Name), standard_name=att(this%STD_Name), units=att(this%Units),  &
                            grid_mapping="", missing_value=undef_dblp)

           type is (INTEGER(kind=ip))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:,:)                                      &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3)                           &
                          , start=[1,1,1],count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2)               &
                          ,UBOUND(values_to_write,dim=3)],                                                                 &
                            long_name=att(this%Long_Name), standard_name=att(this%STD_Name), units=att(this%Units),  &
                            grid_mapping="", missing_value=undef_ip)

           type is (REAL(KIND=silp))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:,:)                                      &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3)                           &
                          , start=[1,1,1],count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2)               &
                          ,UBOUND(values_to_write,dim=3)],                                                                 &
                            long_name=att(this%Long_Name), standard_name=att(this%STD_Name), units=att(this%Units),  &
                            grid_mapping="", missing_value=undef_silp)
            CLASS DEFAULT
              WRITE(*,*) "UNKNOWN TYPE IN io_nc_mod / nc_write"

         END SELECT

        if (this%Cell_Methods /= undef_str)                                                                                 &
           call nc_write_attr(this%nc_file%filename, this%VarName, "cell_methods", trim(this%Cell_Methods))

!~                write(*,*) "wrte tbl_wrte ....", __LINE__, __FILE__
      END SUBROUTINE write_io_grid_3Dvar_noTIME

! ---

      SUBROUTINE write_io_grid_3Dvar_TIME(this, values_to_write, pass_nc, OPT_TimeValue)

        USE ncio, only: nc_write, nc_write_attr

        CLASS(IO_GRID_VAR),       INTENT(IN)   :: this

        INTEGER, INTENT(IN)                    :: pass_nc
        REAL(kind=dblp), OPTIONAL, INTENT(IN)  :: OPT_TimeValue   ! value of the time coordinate (default: pass_nc)

        CLASS(*), DIMENSION(:,:,:), intent(in) :: values_to_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (present(OPT_TimeValue)) then
           call nc_write(this%nc_file%filename,"time",OPT_TimeValue,dim1="time",start=[pass_nc],count=[1])
        else
           call nc_write(this%nc_file%filename,"time",pass_nc,dim1="time",start=[pass_nc],count=[1])
        endif

        select type(values_to_write)

           type is (REAL(kind=dblp))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:,:)                                           &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3),dim4=this%Axes_List(4)         &
                          ,start=[1,1,1,pass_nc]                                                                                &
                          ,count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2),UBOUND(values_to_write,dim=3),1], &
                            long_name=att(this%Long_Name), standard_name=att(this%STD_Name), units=att(this%Units),  &
                            grid_mapping="", missing_value=undef_dblp)

           type is (INTEGER(kind=ip))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:,:)                                           &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3),dim4=this%Axes_List(4)         &
                          ,start=[1,1,1,pass_nc]                                                                                &
                          ,count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2),UBOUND(values_to_write,dim=3),1], &
                            long_name=att(this%Long_Name), standard_name=att(this%STD_Name), units=att(this%Units),  &
                            grid_mapping="", missing_value=undef_ip)

           type is (REAL(KIND=silp))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:,:)                                           &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3),dim4=this%Axes_List(4)         &
                          ,start=[1,1,1,pass_nc]                                                                                &
                          ,count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2),UBOUND(values_to_write,dim=3),1], &
                            long_name=att(this%Long_Name), standard_name=att(this%STD_Name), units=att(this%Units),  &
                            grid_mapping="", missing_value=undef_silp)
            CLASS DEFAULT
              WRITE(*,*) "UNKNOWN TYPE IN io_nc_mod / nc_write"

         END SELECT

        if (this%Cell_Methods /= undef_str)                                                                                 &
           call nc_write_attr(this%nc_file%filename, this%VarName, "cell_methods", trim(this%Cell_Methods))


      END SUBROUTINE write_io_grid_3Dvar_TIME

! ---

      SUBROUTINE wrte_io_nc_axis(this, filename)

        USE ncio, only: nc_write_dim, nc_write, nc_write_attr

        CLASS(IO_NC_AXIS),       INTENT(IN) :: this
        CHARACTER(LEN=*) ,       INTENT(IN) :: filename

        call this%show()
        if (this%is_time) then

          CALL nc_write_dim(filename,trim(this%axis_name),x=1.0, units=trim(this%axis_unit),calendar=trim(this%calendar)   &
                          , unlimited=.TRUE., long_name=att(this%long_name), standard_name=att(this%std_name)            &
                          , axis=att(this%axis_attr))
        else

          if ( allocated(this%axis_array_rsilp) ) then
             CALL nc_write_dim(filename,trim(this%axis_name),x=this%axis_array_rsilp(:),units=trim(this%axis_unit)        &
                             , long_name=att(this%long_name), standard_name=att(this%std_name), axis=att(this%axis_attr))
          elseif ( allocated(this%axis_array_rdblp) ) then
             CALL nc_write_dim(filename,trim(this%axis_name),x=this%axis_array_rdblp(:),units=trim(this%axis_unit)        &
                             , long_name=att(this%long_name), standard_name=att(this%std_name), axis=att(this%axis_attr))
          elseif ( allocated(this%axis_array_isilp) ) then
             CALL nc_write_dim(filename,trim(this%axis_name),x=this%axis_array_isilp(:),units=trim(this%axis_unit)        &
                             , long_name=att(this%long_name), standard_name=att(this%std_name), axis=att(this%axis_attr))
          else
             WRITE(*,*) "No values allocated for given axis, cannot write it in"
             WRITE(*,*) "TROUBLESOME AXIS == ", trim(this%axis_name)
          endif

          ! dmr&clo   CF cell bounds: dimension "bnds" (created once per file), variable <axis>_bnds, attribute bounds
          if ( allocated(this%axis_bounds) ) then
             if (.not. nc_has_dim(filename, "bnds")) CALL nc_write_dim(filename, "bnds", x=1, nx=2)
             CALL nc_write(filename, trim(this%axis_name)//"_bnds", this%axis_bounds(:,:), dim1="bnds",                  &
                           dim2=trim(this%axis_name), grid_mapping="")
             CALL nc_write_attr(filename, trim(this%axis_name), "bounds", trim(this%axis_name)//"_bnds")
          endif

        endif


      END SUBROUTINE wrte_io_nc_axis




!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Run-wide metadata setter and small helpers.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE io_nc_set_metadata(OPT_institution, OPT_author, OPT_source)

       CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: OPT_institution   ! global attribute institution
       CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: OPT_author        ! global attribute author (default: $LOGNAME)
       CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: OPT_source        ! global attribute source (model and version)

         IF (PRESENT(OPT_institution)) THEN
            IF (LEN_TRIM(OPT_institution) > 0) md_institution = OPT_institution
         ENDIF
         IF (PRESENT(OPT_author)) THEN
            IF (LEN_TRIM(OPT_author) > 0) md_author = OPT_author
         ENDIF
         IF (PRESENT(OPT_source)) THEN
            IF (LEN_TRIM(OPT_source) > 0) md_source = OPT_source
         ENDIF

       END SUBROUTINE io_nc_set_metadata

! ---

       SUBROUTINE io_nc_read_ncmeta(nml_unit)

       ! Read the optional group ncmeta (institution, author, source) from an open namelist file and set the run-wide
       ! global attributes. The file is rewound first, so the group can be anywhere. Absent group: defaults kept.

       USE, INTRINSIC :: iso_fortran_env, ONLY: iostat_end

       INTEGER, INTENT(IN) :: nml_unit

       CHARACTER(LEN=str_len) :: institution, author, source
       INTEGER                :: ios

       NAMELIST /ncmeta/ institution, author, source

         institution = ""
         author      = ""
         source      = ""
         REWIND(nml_unit)
         READ(nml_unit, NML=ncmeta, IOSTAT=ios)
         IF (ios == iostat_end) RETURN
         IF (ios /= 0) THEN
            WRITE(*,*) "io_nc: error reading namelist group ncmeta, iostat = ", ios
            STOP 1
         ENDIF
         CALL io_nc_set_metadata(OPT_institution=institution, OPT_author=author, OPT_source=source)

       END SUBROUTINE io_nc_read_ncmeta

! ---

       PURE FUNCTION att(str) result(attval)
       ! attribute value to hand to ncio: blank (= attribute not written) when not set
       CHARACTER(LEN=*), INTENT(IN)  :: str
       CHARACTER(LEN=:), ALLOCATABLE :: attval

         IF (TRIM(str) == undef_str) THEN
            attval = ""
         ELSE
            attval = TRIM(str)
         ENDIF

       END FUNCTION att

! ---

       FUNCTION nc_has_dim(filename, dimname) result(has_dim)

       USE netcdf, ONLY: nf90_open, nf90_close, nf90_inq_dimid, NF90_NOWRITE, NF90_NOERR

       CHARACTER(LEN=*), INTENT(IN) :: filename, dimname
       LOGICAL                      :: has_dim

       INTEGER :: ncid, dimid

         CALL check_nc(nf90_open(TRIM(filename), NF90_NOWRITE, ncid), "nf90_open", filename)
         has_dim = (nf90_inq_dimid(ncid, TRIM(dimname), dimid) == NF90_NOERR)
         CALL check_nc(nf90_close(ncid), "nf90_close", filename)

       END FUNCTION nc_has_dim

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  dmr   Reading side. Uses netcdf-fortran directly on purpose: ncio's readers post-process the values
!        (|x| < 1e-30 set to 0, NaN replaced, scale/offset, missing values) which prevents an exact round trip.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE open_io_nc_file(this, FileName)

       CLASS(IO_NC_FILE), INTENT(OUT):: this
       CHARACTER(LEN=*) , INTENT(IN) :: FileName

       LOGICAL :: exists

         INQUIRE(FILE=TRIM(FileName), EXIST=exists)
         IF (.NOT. exists) THEN
            WRITE(*,*) "io_nc: cannot open for reading, no such file: ", TRIM(FileName)
            STOP 1
         ENDIF
         this%filename    = FileName
         this%initialized = .TRUE.

       END SUBROUTINE open_io_nc_file

! ---

       FUNCTION nrec_io_nc_file(this, OPT_TimeName) result(nrec)

       USE netcdf, ONLY: nf90_open, nf90_close, nf90_inq_dimid, nf90_inquire_dimension, NF90_NOWRITE

       CLASS(IO_NC_FILE),          INTENT(IN) :: this
       CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: OPT_TimeName
       INTEGER(kind=ip)                       :: nrec

       INTEGER :: ncid, dimid, dimlen

         CALL check_nc(nf90_open(TRIM(this%filename), NF90_NOWRITE, ncid), "nf90_open", this%filename)
         IF (PRESENT(OPT_TimeName)) THEN
           CALL check_nc(nf90_inq_dimid(ncid, TRIM(OPT_TimeName), dimid), "nf90_inq_dimid", this%filename)
         ELSE
           CALL check_nc(nf90_inq_dimid(ncid, "time", dimid), "nf90_inq_dimid", this%filename)
         ENDIF
         CALL check_nc(nf90_inquire_dimension(ncid, dimid, len=dimlen), "nf90_inquire_dimension", this%filename)
         CALL check_nc(nf90_close(ncid), "nf90_close", this%filename)
         nrec = dimlen

       END FUNCTION nrec_io_nc_file

! ---

      SUBROUTINE read_io_grid_2Dvar_TIME(this, values_to_read, pass_nc)

        USE netcdf, ONLY: nf90_open, nf90_close, nf90_inq_varid, nf90_get_var, NF90_NOWRITE

        CLASS(IO_GRID_VAR),              INTENT(IN)  :: this
        REAL(kind=dblp), DIMENSION(:,:), INTENT(OUT) :: values_to_read
        INTEGER,                         INTENT(IN)  :: pass_nc

        INTEGER :: ncid, varid

         CALL check_nc(nf90_open(TRIM(this%nc_file%filename), NF90_NOWRITE, ncid), "nf90_open", this%nc_file%filename)
         CALL check_nc(nf90_inq_varid(ncid, TRIM(this%VarName), varid), "nf90_inq_varid "//TRIM(this%VarName),       &
                       this%nc_file%filename)
         CALL check_nc(nf90_get_var(ncid, varid, values_to_read, start=[1,1,pass_nc],                                  &
                       count=[SIZE(values_to_read,dim=1),SIZE(values_to_read,dim=2),1]),                               &
                       "nf90_get_var "//TRIM(this%VarName), this%nc_file%filename)
         CALL check_nc(nf90_close(ncid), "nf90_close", this%nc_file%filename)

      END SUBROUTINE read_io_grid_2Dvar_TIME

! ---

      SUBROUTINE check_nc(status, what, filename)

        USE netcdf, ONLY: NF90_NOERR, nf90_strerror

        INTEGER,          INTENT(IN) :: status
        CHARACTER(LEN=*), INTENT(IN) :: what, filename

         IF (status /= NF90_NOERR) THEN
            WRITE(*,*) "io_nc: ", TRIM(what), " failed on ", TRIM(filename), ": ", TRIM(nf90_strerror(status))
            STOP 1
         ENDIF

      END SUBROUTINE check_nc

! ---

      FUNCTION int_to_str(i) result(res)

      CHARACTER(:),    ALLOCATABLE:: res
      INTEGER(kind=ip),intent(in) :: i

      CHARACTER(RANGE(i)+2) :: tmp

      WRITE(tmp,'(i0)') i

      res = TRIM(tmp)

      END FUNCTION int_to_str

      END MODULE IO_NC_MOD



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
