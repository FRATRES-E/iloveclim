!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009

#if ( CLIO_OUT_NEWGEN == 1 )
      MODULE CLIOGRID_OUTPUTGEN
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      

      USE uuid_fort_wrap, ONLY: uuid_size

      
      IMPLICIT NONE

      PRIVATE
      
      PUBLIC :: CLIOGRID_INIT, clio_grid_var, tracer_grid_Ys, tracer_grid_Xs
      
      CHARACTER(LEN=15) :: filename = "CLIO3_NewGen.nc"
      INTEGER, PARAMETER :: size_x = 120, size_y = 65      

      !--- tracer grid values for a Staggered grid
      DOUBLE PRECISION, DIMENSION(size_x,size_y) :: tracer_grid_Ys, tracer_grid_Xs
      
      INTEGER, PARAMETER :: strlen=256
      CHARACTER(LEN=7), PARAMETER :: undef_str="Not_Set"
      
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      
!  dmr   Holder type for the metada and control data
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      
      
      TYPE :: clio_grid_var

        PRIVATE
        LOGICAL                 :: initialized=.false.
        CHARACTER(LEN=uuid_size):: uuid_var
		CHARACTER(LEN=strlen)   :: VarName
!~      VarSource
!~      MODSource
        CHARACTER(LEN=strlen)   :: Long_Name= undef_str
        CHARACTER(LEN=strlen)   :: STD_Name = undef_str
        CHARACTER(LEN=strlen)   :: Units    = undef_str
        CHARACTER(LEN=strlen)   :: Calendar = "360_day"
                
        CHARACTER(LEN=6), DIMENSION(:), ALLOCATABLE :: Axes_List
        INTEGER                 :: nbaxes=0
  

        CONTAINS
        
          PROCEDURE, PUBLIC  :: init => init_clio_grid_var
          PROCEDURE, PUBLIC  :: show => show_clio_grid_var
          
          PROCEDURE, PRIVATE :: wrte_int => write_clio_grid_var_dbl
          PROCEDURE, PRIVATE :: wrte_dbl => write_clio_grid_var_int
          
          GENERIC, PUBLIC    :: wrte => wrte_int, wrte_dbl
          
      END TYPE clio_grid_var

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|            
      
       CONTAINS

! ---

      SUBROUTINE init_clio_grid_var(this,varname,axeslist,OPT_longname,OPT_stdname,OPT_units,OPT_calendar)
      
        USE uuid_fort_wrap, ONLY: UUID_V4
      
        CLASS(clio_grid_var),       INTENT(OUT) :: this
        CHARACTER(LEN=*),           INTENT(IN)  :: varname
        CHARACTER(LEN=*),           INTENT(IN)  :: axeslist
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: OPT_longname
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: OPT_stdname
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: OPT_units
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: OPT_calendar 
      
        INTEGER :: i  
          
        IF ( .NOT. this%initialized ) then
        
            this%VarName = varname
            
            this%nbaxes = count( (/ (axeslist(i:i), i=1, len_trim(axeslist)) /) == " ")+1
            ALLOCATE(this%Axes_List(this%nbaxes))
            read(axeslist,fmt=*) (this%Axes_List(i),i=1,this%nbaxes)
        
            CALL UUID_V4(this%uuid_var)

            IF ( PRESENT(OPT_longname) ) then
               this%long_name = OPT_longname
            ENDIF
            
            IF ( PRESENT(OPT_stdname) ) then
               this%STD_Name = OPT_stdname
            ENDIF
            
            IF ( PRESENT(OPT_units) ) then
               this%Units = OPT_units
            ENDIF

            IF ( PRESENT(OPT_calendar) ) then
               this%Calendar = OPT_calendar
            ENDIF            
            
        ELSE ! I am called to initialize a variable that already exist ... bizarre!
            WRITE(*,*) "un-implemented [ABORT]"
            STOP 1
        ENDIF
        
        this%initialized = .TRUE.      
      
      END SUBROUTINE init_clio_grid_var

      SUBROUTINE show_clio_grid_var(this)
      
        CLASS(clio_grid_var),       INTENT(IN) :: this

        WRITE(*,*) "INFORMATION FOR /clio_grid_var/ variable"
        WRITE(*,*) "VarName = ", trim(this%VarName)
        WRITE(*,*) "uuidVar = ", this%uuid_var
        WRITE(*,*) "nbaxes  = ", this%nbaxes
        WRITE(*,*) "AxesList= ", this%Axes_List(:)
        WRITE(*,*) "LongName= ", trim(this%Long_Name)
        WRITE(*,*) "STD_Name= ", trim(this%STD_Name)
        WRITE(*,*) "Units   = ", trim(this%Units)
        WRITE(*,*) "Calendar= ", trim(this%Calendar)
        WRITE(*,*) "========================================"

      END SUBROUTINE show_clio_grid_var


      SUBROUTINE CLIOGRID_INIT()  
      
        USE ncio, ONLY: nc_create, nc_write_attr, nc_write_dim, nc_write
        
        REAL, DIMENSION(size_y) :: y_axis_values =                                                                                &
                 (/  -79.5, -76.5, -73.5, -70.5, -67.5, -64.5, -61.5, -58.5, -55.5, -52.5, -49.5, -46.5, -43.5, -40.5, -37.5,     &
                     -34.5, -31.5, -28.5, -25.5, -22.5, -19.5, -16.5, -13.5, -10.5, -7.5, -4.5, -1.5, 1.5, 4.5, 7.5, 10.5, 13.5,  &
                      16.5, 19.5, 22.5, 25.5, 28.5, 31.5, 34.5, 37.5, 40.5, 43.5, 46.5, 49.5, 52.5, 55.5, 58.5, 61.5, 64.5, 67.5, &
                      70.5, 73.5, 76.5, 79.5, 82.5, 85.5, 88.5, 91.5, 94.5, 97.5, 100.5, 103.5, 106.5, 109.5, 112.5 /)

        REAL, DIMENSION(size_x) :: x_axis_values =                                                                                &
                   (/ 28.5, 31.5, 34.5, 37.5, 40.5, 43.5, 46.5, 49.5, 52.5, 55.5, 58.5, 61.5, 64.5, 67.5, 70.5, 73.5, 76.5, 79.5, &
                      82.5, 85.5, 88.5, 91.5, 94.5, 97.5, 100.5, 103.5, 106.5, 109.5, 112.5, 115.5, 118.5, 121.5, 124.5, 127.5,   &
                      130.5, 133.5, 136.5, 139.5, 142.5, 145.5, 148.5, 151.5, 154.5, 157.5, 160.5, 163.5, 166.5, 169.5, 172.5,    &
                      175.5, 178.5, 181.5, 184.5, 187.5, 190.5, 193.5, 196.5, 199.5, 202.5, 205.5, 208.5, 211.5, 214.5, 217.5,    &
                      220.5, 223.5, 226.5, 229.5, 232.5, 235.5, 238.5, 241.5, 244.5, 247.5, 250.5, 253.5, 256.5, 259.5, 262.5,    &
                      265.5, 268.5, 271.5, 274.5, 277.5, 280.5, 283.5, 286.5, 289.5, 292.5, 295.5, 298.5, 301.5, 304.5, 307.5,    &
                      310.5, 313.5, 316.5, 319.5, 322.5, 325.5, 328.5, 331.5, 334.5, 337.5, 340.5, 343.5, 346.5, 349.5, 352.5,    &
                      355.5, 358.5, 361.5, 364.5, 367.5, 370.5, 373.5, 376.5, 379.5, 382.5, 385.5 /)
        
        CHARACTER(LEN=256) :: logname, line
        CHARACTER(LEN=256) :: institution_name
        CHARACTER(LEN=256) :: command_line
        
        integer :: length, rc
        
        command_line='curl -s https://json.geoiplookup.io/$(curl -s https://ipinfo.io/ip)'
        command_line=TRIM(command_line)//'|grep country_code | cut --delimiter=: -f2 > country_code.info'
        
        CALL EXECUTE_COMMAND_LINE(""//TRIM(command_line))

        OPEN(10540,file="country_code.info",form='formatted')
        READ(10540,*) line
        CLOSE(10540)
        
        IF (TRIM(line)=="FR") THEN
           institution_name="Laboratoire de Sciences du Climat et de l'Environnement"
        ELSE IF(TRIM(line)=="NL") THEN
           institution_name="Vrije Universiteit Amsterdam"
        ELSE
           WRITE(*,*) "Institution == ", TRIM(line)
           institution_name="Unknown"
        ENDIF
                         
        CALL GET_ENVIRONMENT_VARIABLE('LOGNAME', logname, length, rc)
        
        IF (rc == 0) THEN
            CALL nc_create(filename,overwrite=.TRUE.,netcdf4=.TRUE., author=TRIM(logname)) 
        ELSE
            CALL nc_create(filename,overwrite=.TRUE.,netcdf4=.TRUE., author="foo") 
        ENDIF
    
    
        CALL nc_write_attr(filename,"Title","CLIO_grid NewGen output")
        CALL nc_write_attr(filename,"Institution", TRIM(institution_name))

        CALL nc_write_dim(filename,"ptlon",x=x_axis_values,units="degrees_east")
        CALL nc_write_dim(filename,"ptlat",x=y_axis_values,units="degrees_north")
        CALL nc_write_dim(filename,"time",x=1.0, &
                      units="years",calendar="360_day", unlimited=.TRUE.)
      
      
        CALL nc_write(filename,"tlat",tracer_grid_Ys(:,:),dim1="ptlon",dim2="ptlat",start=[1,1],count=[size_x,size_y])
        CALL nc_write(filename,"tlon",tracer_grid_Xs(:,:),dim1="ptlon",dim2="ptlat",start=[1,1],count=[size_x,size_y])             
                     
      END SUBROUTINE CLIOGRID_INIT    
      
! ---

      SUBROUTINE write_clio_grid_var_dbl(this, values_to_write, pass_nc)  
      
        USE ncio, only: nc_write
      
        CLASS(clio_grid_var),       INTENT(IN) :: this
      
        INTEGER, INTENT(IN) :: pass_nc

        DOUBLE PRECISION, DIMENSION(size_x,size_y), intent(in) :: values_to_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
        call nc_write(filename,"time",pass_nc,dim1="time",start=[pass_nc],count=[1])
     
        call nc_write(filename,this%VarName,values_to_write(:,:), dim1=this%Axes_List(1),dim2=this%Axes_List(2)      &
                    , dim3=this%Axes_List(3),start=[1,1,pass_nc],count=[size_x,size_y,1])
      
      END SUBROUTINE write_clio_grid_var_dbl

! ---

      SUBROUTINE write_clio_grid_var_int(this, values_to_write, pass_nc)  
      
        USE ncio, only: nc_write
      
        CLASS(clio_grid_var),       INTENT(IN) :: this
      
        INTEGER, INTENT(IN) :: pass_nc

        INTEGER, DIMENSION(size_x,size_y), intent(in) :: values_to_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
        call nc_write(filename,"time",pass_nc,dim1="time",start=[pass_nc],count=[1])
     
        call nc_write(filename,this%VarName,values_to_write(:,:), dim1=this%Axes_List(1),dim2=this%Axes_List(2)      &
                    , dim3=this%Axes_List(3),start=[1,1,pass_nc],count=[size_x,size_y,1])
      
      END SUBROUTINE write_clio_grid_var_int
      
! ---
      
      END MODULE CLIOGRID_OUTPUTGEN
#endif
