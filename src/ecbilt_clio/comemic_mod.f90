!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module est un portage en FORTRAN 90 du module coupler initial
!       cree en FORTRAN 77.
!      (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche
!      Date   : 19 Aout 2014
!      Derniere modification : 10 decembre 2015
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       module comemic_mod


       use global_constants_mod, only: dblp=>dp, ip, str_len, stdout
       use comatm, only: nlat, nlon

       implicit none

       private :: nlat, nlon, dblp, ip, str_len, stdout

!-----------------------------------------------------------------------
! *** File:     comemic.h
! *** Contents: Common declarations global control variables of EMIC
!-----------------------------------------------------------------------
! *** COMMON /ec_timectl/ day,ntstep,nyears,nstpyear,nocstpyear,iatm,
!                      iyear,imonth,iday
!     day:       keeps track of time in days during the integration
!     ntstep:    total number of timesteps of the integration
!     nyears:    total integration period in years
!     nstpyear:  number of atmospheric timesteps that fit in one year
!     nocstpyear:number of ocean timesteps that fit in one year
!     iatm:      number of atmospheric timesteps that fit in one day
!     nwrskip    number of years between writing the model state to disk
!     iyear:     counts the years during the integration
!     imonth:    counts the months during the integration
!     iday:      counts the days during the integration
!     iseason:   counts the seasons during the integration
!
! *** COMMON /ec_startctl/ irunlabel,fini,fend
!     irunlabel:index of startfile containing initial conditions
!               if zero, then initial conditions are prescribed
!               read in iatmdyn,iatmphys, iniland,
!               iniocean, iniseaice
!     fini:     character string which contains irunlabel
!     fend:     character string which contains irunlabel+nyears
!
! *** COMMON /ec_coupctl/ idtbclin,idtbtrop,ibclins,nbclins,ibtrops,nbtrops
!     idtbclin: baroclinic timestep of ocean in days
!     idtbtrop: barotropic timestep of ocean in days
!     ibclins:  counts the number of atmospheric steps in one baroclinic
!               ocean step
!     nbclins:  total number of atmospheric timesteps in one barocinic
!               ocean step
!     ibtrops:  counts the number of atmospheric steps in one barotropic
!               ocean step
!     nbtrops:  total number of atmospheric timesteps in one barotropic
!               ocean step
!
!     lferCO2: fertilization effect of CO2 if T
!     lradCO2: radiative forcing of CO2 if T
!
!-----------------------------------------------------------------------
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8)                       :: day,undef
      real(kind=8), dimension(nlat)      :: dareafac
      real(kind=8), dimension(nlat,nlon) :: fracto

      integer                            :: nyears,ndays,irunlabel,irunlabeld,iatm,ilan,iice,iobtrop,iobclin,nwrskip,nwrskip_days &
             ,nbclins,nbtrops,end_year,end_day


      integer                            :: ntstep,nstpyear,isatfor,ntotday
      integer                            :: iseason
      integer, target                    :: imonth,iday
      integer, target                    :: nocstpyear
      integer, target                    :: iyear
      integer                            :: is_ism,if_ism,kism,tstartism
      character(len=6)                   :: fini
      logical                            :: flgveg,flgism,flgloch,flgisma,flgismg

!dmr @-@ iceb0
      logical                            :: flag_snow_to_flux
!dmr @-@ iceb0

      logical                            :: flgtsi,flgvol,flgghg,flgsul
      logical                            :: lferCO2,lradCO2

!~ [UNUSED???]
!~       logical                            :: initialization

      character(len=120), dimension(26,2):: globalatt

! --- BdB 05-2019: new variables for writing output
      integer                            :: new_year_atm        ! restarting counting the years every nwrskip interval
      integer                            :: current_int_atm     ! when time interval is resetted, add nwrskip for writing time
      integer                            :: new_year_veg        ! restarting counting the years every nwrskip interval
      integer                            :: current_int_veg     ! when time interval is resetted, add nwrskip for writing time
      integer                            :: init_year_clio      ! first year of clio timer
      real(kind=8), dimension(:), allocatable :: time_in_years  ! time in years for ecbilt and clio, folliwing the runlabel
! --- BdB - end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "albedoclio"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real(dblp), dimension(nlat,4) :: albsea


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Pretty print of execution time in emic ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!     ! ADDING DATE and TIME !mohr
      character(len=8)     :: date
      character(len=10)    :: time
      character(len=5)     :: zone
      integer,dimension(8) :: values

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       contains

       subroutine pretty_print_exec_time(str_in)

         character(len=*), intent(in) :: str_in

!! ADDING DATE AND TIME TO CHECK CORRECT EXECUTION... !mohr
         write(stdout,*) "      "
         write(stdout,*) "      "
         call date_and_time(date,time,zone,values) !mohr
         write(stdout,*) "iLOVECLIM@: "//trim(str_in)
         write(stdout,*) "DT%% "//date(1:4)//"-"//date(5:6)//"-"//date(7:8)," ", time(1:2)//":"//time(3:4)//":"//time(5:6)
         write(stdout,*) "      "
         write(stdout,*) "      "

       end subroutine pretty_print_exec_time

      end module comemic_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
