!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of iLOVECLIM/COUPLER
!!      iLOVECLIM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
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
!dmr -- Adding the choice of components through the pre-processing options

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [initcoupledmodel_mod]
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module [initcoupledmodel_mod] is handling the initialization of the coupled mode iLOVECLIM
!
!>     @date Creation date: October, 04th, 2019
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     The aim is to progressively move all initialization statements within this module
!!      v0.1 : added ec_initemic
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module initcoupledmodel_mod

       use file_libs, only: fileDescriptor ! file type for reading
       use global_constants_mod, only: str_len, dblp=>dp, ip

       implicit none

       private

       public :: ec_initemic, init_coupled_components

      ! NOTE_avoid_public_variables_if_possible


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
! dmr  First define input files path and format ...
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      ! in the old version:
      ! iuo+45 == mozaic
      ! iuo+46 == namelist
      ! iuo+48 == fractoc
      ! iuo+50 == darea

! dmr --- File variable
       type(fileDescriptor), public,save :: mozaic_f, namelist_f, fractoc_f, darea_f

! dmr [2024-01-25] mozaic_nm is not used here! REMOVE?

! dmr --- File name variable (nm)
       character(len=str_len), private ::                                               &
                                         mozaic_nm   = "inputdata/clio/mozaic.w",       &
                                         namelist_nm = "namelist",                      &
                                         fractoc_nm  = "inputdata/fractoc.dat",         &
                                         darea_nm    = "inputdata/darea.dat"

! dmr --- Formatted flag: .false. is unformatted file
       logical, parameter, private      ::                                              &
                                         mozaic_fm   = .false. ,                        &
                                         namelist_fm = .true. ,                         &
                                         fractoc_fm  = .true. ,                         &
                                         darea_fm    = .true.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
! dmr  Second define output files path and format ...
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr --- File variable
       type(fileDescriptor), public    :: parlist_f

! dmr --- File name variable (nm)
       character(len=str_len), private :: parlist_nm = "outputdata/globals/parlist"

! dmr --- Formatted flag: .false. is unformatted file
       logical, parameter, private      :: parlist_fm = .true.




      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [ec_initemic]
!
!>     @brief This function is a module f90 port of the initial subroutine with identical name
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ec_initemic() result(returnValue)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Module imports
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( COMATM == 1 )
      use comatm,      only: nlat, nlon
      use comemic_mod, only: fini, flgism, flgisma, flgismg, flgloch, flgveg, iatm, if_ism, iice, ilan, iobclin, iobtrop          &
     &               , irunlabel, irunlabeld, is_ism, kism, lferco2, lradco2, nbclins, nbtrops, ndays, nocstpyear, nstpyear       &
     &               , ntotday, ntstep, nwrskip, nwrskip_days, nyears, undef, fracto, dareafac

      ! BdB 05-2019: added time_in_years for allocation
      use comemic_mod, only: globalatt, iyear, time_in_years
      use comunit,     only: iuo
#endif

#if ( BATHY >= 2 || NC_BERG == 2 )
      use update_clio_bathy_tools, only: la_date
      use palaeo_timer_mod, only: palaeo_year
#endif

#if ( BRINES >= 3 )
      use brines_mod, only: la_date_brines
#endif

      use file_libs, ONLY: fileDescriptor, open_f, close_f

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  void
!>    @param[out] void
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical     :: returnValue

       integer            :: ija,i,j,k

      TYPE(fileDescriptor) emicparam

       integer(ip), parameter :: ijatm=nlat*nlon, ismfile = 400
       real(kind=dblp)    fractocn(ijatm)

       character(len=6)   :: num_startyear
       character(len=3)   :: num_startday

!nb
#if (BATHY >= 2 )
       character*30 name_file
       character(len=5) :: charI
#endif

       NAMELIST /tstepctl/nyears,irunlabel,iatm,ilan,iice,iobtrop,iobclin,nwrskip


!dmr ### lln 1.22      NAMELIST /tstepctl/nyears,ndays,irunlabel,irunlabeld,iatm,ilan,iice,iobtrop,
!dmr ### lln 1.22     &                   iobclin,nwrskip,nwrskip_days
!dmr --> no init of ndays, irunlabeld, nwrskip_days

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!     *** open emic.param

      call open_f(emicparam, 'emic.param')

      read(emicparam%id,*)
      read(emicparam%id,*)
      read(emicparam%id,'(L4)') flgveg
!dmr [DEPRECATED]      read(emicparam%id,'(L4)') flgicb
!dmr following just passes the line with flgicb
      read(emicparam%id,*)
      read(emicparam%id,'(L4)') flgismg
      read(emicparam%id,'(L4)') flgisma
      read(emicparam%id,'(L4)') flgloch
      read(emicparam%id,*)
      read(emicparam%id,*)
      read(emicparam%id,*)
      read(emicparam%id,'(L4)') lradCO2
      read(emicparam%id,'(L4)') lferCO2
      read(emicparam%id,*)
      read(emicparam%id,*)

      do i=1,26
       if (i.ne.6.AND.i.ne.19.AND.i.ne.22.AND.i.ne.24.AND.i.ne.25.AND.i.ne.26) then
            read(emicparam%id,*)
            read(emicparam%id,'(A)') globalatt(i,1)
            read(emicparam%id,'(A)') globalatt(i,2)
         end if
      enddo

      call close_f(emicparam)


      write(*,'(A,L4)') 'Vecode:',flgveg
! dmr deprecated
!      write(*,'(A,L4)') 'GISM:',flgismg
!      write(*,'(A,L4)') 'AISM:',flgisma
!      write(*,'(A,L4)') 'Loch:',flgloch
! dmr deprecated
      write(*,*)
      write(*,'(A,L4)') 'CO2 radiative forcing:',lradCO2
      write(*,'(A,L4)') 'CO2 fertilization:',lferCO2
      write(*,*)

!     *** open namelist

!-----------------------------------------------------------------------
! *** open & read statements of global input files:
!-----------------------------------------------------------------------

      nyears=10
      ndays=0
      irunlabel=000000
      irunlabeld=0
      iatm=6
      ilan=1
      iice=3
      iobtrop=1
      iobclin=1
      nwrskip=50
      nwrskip_days=0

! dmr --- Set the formatted type and then open ...
      namelist_f%isFormatted = namelist_fm
      call open_f(namelist_f, namelist_nm)

! dmr --- iuo+46 == "namelist"
! dmr --- read the namelist data, tstepctl defined above
      read(namelist_f%id, NML = tstepctl)

      call close_f(namelist_f)

      write(fini,1) irunlabel
 1    format(i6.6)

      irunlabeld=(iyear+1)*360
      if(irunlabeld.lt.360) then
        write(num_startyear,'(i6.6)')irunlabel
        write(num_startday,'(i3.3)')irunlabeld+1
      else
        write(num_startyear,'(i6.6)')irunlabel+1
        write(num_startday,'(i3.3)')1
      endif
      globalatt(19,1)="branch_time"
      globalatt(19,2)=""//fini

      kism=1

      is_ism=0
      if_ism=33

      if(flgismg.or.flgisma) then
         flgism=.TRUE.
      else
         flgism=.FALSE.
      endif

      if(.not.flgismg) then
         if_ism=15
      endif
      if(.not.flgisma) then
         is_ism=15
      endif


! dmr --- Set the formatted type and then open ...
      parlist_f%isFormatted = parlist_fm
      call open_f(parlist_f, parlist_nm)

! dmr --- Write the ncessary data
      write(parlist_f%id, 900) 'nyears   =', nyears
      write(parlist_f%id, 900) 'irunlabel=', irunlabel
      write(parlist_f%id, 900) 'iatm     =', iatm
      write(parlist_f%id, 900) 'ilan     =', ilan
      write(parlist_f%id, 900) 'iice     =', iice
      write(parlist_f%id, 900) 'iobtrop  =', iobtrop
      write(parlist_f%id, 900) 'iobclin  =', iobclin
      write(parlist_f%id, 900) 'nwrskip  =', nwrskip
      write(parlist_f%id, 900) 'flgveg   =', abs(transfer(flgveg,nwrskip))
!dmr [DEPRECATED]      write(parlist_f%id, 900) 'flgicb   =', abs(transfer(flgicb,nwrskip))
      write(parlist_f%id, 900) 'flgismg  =', abs(transfer(flgismg,nwrskip))
      write(parlist_f%id, 900) 'flgisma  =', abs(transfer(flgisma,nwrskip))
      write(parlist_f%id, 900) 'flgloch  =', abs(transfer(flgloch,nwrskip))

      call close_f(parlist_f)

      undef = 9.99E10

!     *** nstpyear is number of atmospheric timesteps per year
!     *** nocstpyear is number of ocean timesteps per year
!     *** ntstep is total number of timesteps
!     *** nbclins is number of atmospheric time steps per baroclinic ocean
!     *** timestep
!     *** nbtrops is number of atmospheric time steps per barotropic ocean
!     *** timestep


      nstpyear   = iatm*360
      nocstpyear = 360/iobclin
      ntstep     = nstpyear*nyears
      ntotday    = nyears*360
      nbclins    = iatm*iobclin
      nbtrops    = iatm*iobtrop

! --- BdB 05-2019: Create array of output year.
      allocate(time_in_years(0:nyears))
      do j = 0,nyears
        time_in_years(j) = real(irunlabel+j)
      end do

! dmr --- Set the formatted type
      fractoc_f%isFormatted = fractoc_fm

#if ( BATHY <= 1 )
! dmr --- then open ...
      call open_f(fractoc_f, fractoc_nm)

#if ( NC_BERG == 2 )
      la_date=palaeo_year-iyear

      write(*,*) 'The date for ice sheet update is ', la_date, palaeo_year, iyear
#endif

#else
! dmr --- transient bathymetry, not the same file name

      la_date=palaeo_year-iyear

      write(*,*) 'The date for bathy update is ', la_date

      write(charI,'(I5.5)'), la_date
      name_file ='inputdata/fractoc_'//trim(charI)//'.dat'

! dmr --- then open ...

      call open_f(fractoc_f, name_file)

#endif

#if ( BRINES >= 3 )
      la_date_brines=palaeo_year-iyear

#endif

      read(fractoc_f%id,*)
      read(fractoc_f%id,*) (fractocn(ija),ija=1,ijatm)

! dmr
!>    @bug Since there is a rewind here, this means that it is used elsewhere! Coherency !!
      rewind(fractoc_f%id)
! dmr Replace this by a close_f for coherency, but bad idea in principle!
! dmr
      call close_f(fractoc_f)
      do ija=1,ijatm
         j=int((ija-1)/nlon)+1
         i=ija-(j-1)*nlon
         fracto(j,i)=fractocn(ija)
         if (fracto(j,i).gt.0.990) fracto(j,i)=1.0d0
      enddo

! dmr --- Set the formatted type and then open ...
      darea_f%isFormatted = darea_fm
      call open_f(darea_f, darea_nm)


      do i=1,nlat
         read(darea_f%id,*) dareafac(i)
      enddo

      call close_f(darea_f)

 900  format(a12,1x,i6)

      returnValue = .true.

      return

      end function ec_initemic

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: init_coupled_components
!
!>     @brief This function is a module f90 port of the initial subroutine with identical name
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function init_coupled_components(day,month) result(returnValue)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Module imports
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


        use comemic_mod, only: irunlabel

#if ( OCYCC == 1 )
        use ocycc_main, only: ocycc_ini
#endif

#if ( ICEBERG > 0 )
        use iceberg_mod, only: read_positb, init_iceberg, read_iceberg_restartdata
#endif

        use landmodel_mod, only: ec_initlbm

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  day, month
!>    @param[out] void
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer(ip), intent(in) :: day, month

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical     :: returnValue

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       call ec_initecbilt
       call init_clio(irunlabel)

#if ( PERTOCEAN == 1 )
       !*** PERTURBATION OF THE OCEAN VIA scal
       call ec_pertscal(1)
#endif
#if ( CFC == 1 )
       call lire_cfc
#endif

       call ec_initlbm
       call ec_initcoup

#if ( OCYCC == 1 )
       returnValue = ocycc_ini(day,month)
#endif

#if ( ICEBERG > 0 )
! Initialize iceberg model
      call init_iceberg

! read positb.init
      call read_positb

! Read iceberg_model restartdata
      call read_iceberg_restartdata

#endif

       returnValue = .true.

       return

      end function init_coupled_components

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module initcoupledmodel_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
