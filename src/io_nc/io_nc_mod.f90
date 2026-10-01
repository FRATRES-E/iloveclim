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
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        CHARACTER(LEN=5), PARAMETER :: version_mod ="0.5.0"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER,          PARAMETER :: axisNmSize = 6
        CHARACTER(LEN=7), PARAMETER :: undef_str  = "NotSet"
        REAL(KIND=dblp),  PARAMETER :: undef_dblp = (-1.0_dblp)*HUGE(1.0_silp)
        REAL(KIND=silp),  PARAMETER :: undef_silp = (-1.0_silp)*HUGE(1.0_silp)
        INTEGER(KIND=ip), PARAMETER :: undef_ip   = (-1_ip)*HUGE(1_ip)
        INTEGER(KIND=ip), PARAMETER :: time_chunk = 10

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

       SUBROUTINE init_io_nc_axis(this, AxisName, OPT_vals_to_write1D, OPT_AxisUnit, OPT_isTime, OPT_calendar)

!~        USE uuid_fort_wrap, ONLY: UUID_V4
       use uuid_module, only: generate_uuid

       CLASS(IO_NC_AXIS)               , INTENT(OUT):: this
       CHARACTER(LEN=*)                , INTENT(IN) :: AxisName
       CHARACTER(LEN=*),       OPTIONAL, INTENT(IN) :: OPT_AxisUnit
       LOGICAL,                OPTIONAL, INTENT(IN) :: OPT_isTime
       CHARACTER(LEN=*),       OPTIONAL, INTENT(IN) :: OPT_calendar

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


         IF ( .NOT. this%initialized ) then

!~            CALL UUID_V4(this%uuid_file)
           this%uuid_file = generate_uuid(4)
           this%filename = FileName


           CALL GET_ENVIRONMENT_VARIABLE('LOGNAME', logname, length, rc)
           IF (rc == 0) THEN
              this%short_author_name = TRIM(logname)
           ENDIF

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

           CALL nc_create(this%filename,overwrite=this%overwrte_s,netcdf4=this%netcdf4f_s, author=TRIM(this%short_author_name))
           CALL nc_write_attr(this%filename,"Title",TRIM(this%title_file))
           CALL nc_write_attr(this%filename,"Institution", TRIM(this%institution_name))

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

      SUBROUTINE init_io_grid_var(this,varname,nc_file,axeslist,OPT_longname,OPT_stdname,OPT_units)

!~         USE uuid_fort_wrap, ONLY: UUID_V4
        use uuid_module, only: generate_uuid
       
        CLASS(IO_GRID_VAR),         INTENT(OUT) :: this
        CHARACTER(LEN=*),           INTENT(IN)  :: varname
        TYPE(IO_NC_FILE), TARGET,   INTENT(IN)  :: nc_file
        CHARACTER(LEN=*),           INTENT(IN)  :: axeslist
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: OPT_longname
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: OPT_stdname
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: OPT_units

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

      SUBROUTINE write_io_grid_2Dvar_TIME(this, values_to_write, pass_nc)

        USE ncio, only: nc_write

        CLASS(IO_GRID_VAR),       INTENT(IN) :: this

        INTEGER, INTENT(IN)                  :: pass_nc

        CLASS(*), DIMENSION(:,:), intent(in) :: values_to_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        call nc_write(this%nc_file%filename,"time",pass_nc,dim1="time",start=[pass_nc],count=[1])

        select type(values_to_write)

           type is (REAL(kind=dblp))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:)                                        &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3)                           &
                          ,start=[1,1,pass_nc],count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2),1],      &
                            missing_value=undef_dblp)

           type is (INTEGER(kind=ip))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:)                                        &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3)                           &
                          ,start=[1,1,pass_nc],count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2),1],      &
                            missing_value=undef_ip)

           type is (REAL(KIND=silp))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:)                                        &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3)                           &
                          ,start=[1,1,pass_nc],count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2),1],      &
                            missing_value=undef_silp)
            CLASS DEFAULT
              WRITE(*,*) "UNKNOWN TYPE IN io_nc_mod / nc_write"

         END SELECT


      END SUBROUTINE write_io_grid_2Dvar_TIME

! ---

      SUBROUTINE write_io_grid_3Dvar_noTIME(this, values_to_write)

        USE ncio, only: nc_write

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
                          , missing_value=undef_dblp)

           type is (INTEGER(kind=ip))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:,:)                                      &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3)                           &
                          , start=[1,1,1],count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2)               &
                          ,UBOUND(values_to_write,dim=3)],                                                                 &
                            missing_value=undef_ip)

           type is (REAL(KIND=silp))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:,:)                                      &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3)                           &
                          , start=[1,1,1],count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2)               &
                          ,UBOUND(values_to_write,dim=3)],                                                                 &
                            missing_value=undef_silp)
            CLASS DEFAULT
              WRITE(*,*) "UNKNOWN TYPE IN io_nc_mod / nc_write"

         END SELECT

!~                write(*,*) "wrte tbl_wrte ....", __LINE__, __FILE__
      END SUBROUTINE write_io_grid_3Dvar_noTIME

! ---

      SUBROUTINE write_io_grid_3Dvar_TIME(this, values_to_write, pass_nc)

        USE ncio, only: nc_write

        CLASS(IO_GRID_VAR),       INTENT(IN)   :: this

        INTEGER, INTENT(IN)                    :: pass_nc

        CLASS(*), DIMENSION(:,:,:), intent(in) :: values_to_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        call nc_write(this%nc_file%filename,"time",pass_nc,dim1="time",start=[pass_nc],count=[1])

        select type(values_to_write)

           type is (REAL(kind=dblp))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:,:)                                           &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3),dim4=this%Axes_List(4)         &
                          ,start=[1,1,1,pass_nc]                                                                                &
                          ,count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2),UBOUND(values_to_write,dim=3),1], &
                            missing_value=undef_dblp)

           type is (INTEGER(kind=ip))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:,:)                                           &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3),dim4=this%Axes_List(4)         &
                          ,start=[1,1,1,pass_nc]                                                                                &
                          ,count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2),UBOUND(values_to_write,dim=3),1], &
                            missing_value=undef_ip)

           type is (REAL(KIND=silp))
              call nc_write(this%nc_file%filename,this%VarName,values_to_write(:,:,:)                                           &
                          , dim1=this%Axes_List(1),dim2=this%Axes_List(2),dim3=this%Axes_List(3),dim4=this%Axes_List(4)         &
                          ,start=[1,1,1,pass_nc]                                                                                &
                          ,count=[UBOUND(values_to_write,dim=1),UBOUND(values_to_write,dim=2),UBOUND(values_to_write,dim=3),1], &
                            missing_value=undef_silp)
            CLASS DEFAULT
              WRITE(*,*) "UNKNOWN TYPE IN io_nc_mod / nc_write"

         END SELECT


      END SUBROUTINE write_io_grid_3Dvar_TIME

! ---

      SUBROUTINE wrte_io_nc_axis(this, filename)

        USE ncio, only: nc_write_dim

        CLASS(IO_NC_AXIS),       INTENT(IN) :: this
        CHARACTER(LEN=*) ,       INTENT(IN) :: filename

        call this%show()
        if (this%is_time) then

          CALL nc_write_dim(filename,trim(this%axis_name),x=1.0, units=trim(this%axis_unit),calendar=trim(this%calendar)   &
                          , unlimited=.TRUE.)
        else

          if ( allocated(this%axis_array_rsilp) ) then
             CALL nc_write_dim(filename,trim(this%axis_name),x=this%axis_array_rsilp(:),units=trim(this%axis_unit))
          elseif ( allocated(this%axis_array_rdblp) ) then
             CALL nc_write_dim(filename,trim(this%axis_name),x=this%axis_array_rdblp(:),units=trim(this%axis_unit))
          elseif ( allocated(this%axis_array_isilp) ) then
             CALL nc_write_dim(filename,trim(this%axis_name),x=this%axis_array_isilp(:),units=trim(this%axis_unit))
          else
             WRITE(*,*) "No values allocated for given axis, cannot write it in"
             WRITE(*,*) "TROUBLESOME AXIS == ", trim(this%axis_name)
          endif

        endif


      END SUBROUTINE wrte_io_nc_axis




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
