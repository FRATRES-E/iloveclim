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
!dmr -- Adding the choice of components through the pre-processing options

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [Module_Name here]
!
!>     @author  Aurélien Quiquet (afq)
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module [Foo_mod] is handling the creation of a particuliar type of fried noodles ...
!
!>     @date Creation date: November, 41st, 2999
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : tfa, tsa
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module zonesOrmenfwf_mod

       use global_constants_mod, only: dp, ip


       implicit none
       private

!-----------------------------------------------------------------------
! *** file:     zonesormen.h
!dmr --- ---------------------------------------------------------------
! *** contents: variables needed for an addional hyperpycnal forcing
!-----------------------------------------------------------------------
!dmr     --- created : Didier M. Roche
!dmr     --- last updated : september, 26th, 2007
!dmr&afq --- portage .f90 : april, 2017
!dmr --- ---------------------------------------------------------------

! *** variables
       ! integer, parameter :: nbannees = 25, nb_zonesmax = 20
       ! dmr --- a priori nbannees n est pas utilisee ...
       integer, parameter :: nb_zonesmax = 20
       integer, parameter :: size_lon = 122, size_lat = 65

! parametres de tab_zonesormen : nbzones, nbmaxcases, (info, i, j)
       integer, dimension(:,:,:), allocatable, save :: tab_zonesormen     ! toutes les zones possibles (definitions geogra.)
       integer, dimension(:)    , allocatable, save :: datefwf_zonesormen !
       integer, dimension(:)    , allocatable, save :: ID_zonesormen      !
       real(dp), dimension(:)   , allocatable, save :: aires_zonesormen   ! in ?, pour toutes les zones possibles
       real(dp), dimension(:,:),  allocatable, save :: fwf_zonesormen     ! in Sv, seulement pour les zones voulues

! variable interne pour maj temporelle
       real(dp), dimension(:)   , allocatable       :: fwfzones ! in Sv


! variable globale pour la definition du temps

       integer, save :: tps_deb, tps_fin

       public :: init_zonesOrmenfwf_mod
       public :: update_fwf
       public :: set_datesBP

      ! NOTE_avoid_public_variables_if_possible

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [Routine or Function Name Here]
!
!>     @brief This subroutine / function is adding a wonderful knew unused functionality
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine init_zonesOrmenfwf_mod(area_ocean)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real, dimension(:,:), intent(in) :: area_ocean

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer :: nbzones, i, j, k

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       call zonesormen_init()
       call aires_init(area_ocean)
       call lect_fwf_init(tps_deb,tps_fin)

       write(*,*) "DATES SET in ORMEN ==", tps_deb, tps_fin

       if (allocated(fwfzones)) then
         deallocate(fwfzones)
       endif
       allocate(fwfzones(ubound(fwf_zonesormen,dim=2)))

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       end subroutine init_zonesOrmenfwf_mod ! of subroutine

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [Routine or Function Name Here]
!
!>     @brief This subroutine / function is adding a wonderful knew unused functionality
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine zonesormen_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr&afq <void>
      use global_constants_mod, only: dblp=>dp, ip


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer j, IOstatus, indx_zone, nb_cazes_vraies
       integer, parameter :: iuo = 1234, nb_cazesmax = 100, nb_geo = 2

       integer, dimension(nb_zonesmax,0:nb_cazesmax,0:nb_geo) :: locale_zone

       integer(kind=ip) :: zonesOrmenfwf_data_dat_id
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!>    @bug Description of the stupid sticky bug that we know exist there but is not corrected yet!

       print*, "trying to read the imposed zones for ORMEN"

       open(newunit=zonesOrmenfwf_data_dat_id,  file='inputdata/clio/zonesOrmenfwf-data.dat',iostat=IOstatus)

       read(zonesOrmenfwf_data_dat_id,*)

       indx_zone = 0
       nb_cazes_vraies = 0

       do while (IOstatus.NE.-1)

            indx_zone = indx_zone + 1

            read(zonesOrmenfwf_data_dat_id,*) locale_zone(indx_zone,0,0)
            nb_cazes_vraies = max(nb_cazes_vraies,locale_zone(indx_zone,0,0))

            do j=1,locale_zone(indx_zone,0,0)
               read(zonesOrmenfwf_data_dat_id,*) locale_zone(indx_zone,j,1),locale_zone(indx_zone,j,2)
            enddo

            read(zonesOrmenfwf_data_dat_id,*,iostat=IOstatus)
       enddo

       if (allocated(tab_zonesormen)) then
          deallocate(tab_zonesormen)
       endif

       allocate(tab_zonesormen(1:indx_zone,0:nb_cazes_vraies,0:nb_geo))

       tab_zonesormen(:,:,:) = locale_zone(1:indx_zone,0:nb_cazes_vraies,0:nb_geo)

       close(zonesOrmenfwf_data_dat_id)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       end subroutine zonesormen_init! of subroutine

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [Routine or Function Name Here]
!
!>     @brief This subroutine / function is adding a wonderful knew unused functionality
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine aires_init(aires_ocean)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real, dimension(:,:), intent(in) :: aires_ocean

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer :: nbzones, i, j, k

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       nbzones = UBOUND(tab_zonesormen,dim=1)

       if (allocated(aires_zonesormen)) then
         deallocate(aires_zonesormen)
       endif

       allocate(aires_zonesormen(nbzones))

       aires_zonesormen(:) = 0.0

       do i=1,nbzones  ! nombre de zones forcees

          do j=1,tab_zonesormen(i,0,0) ! nombre cazes dans la zone

            if (tab_zonesormen(i,j,1).eq.0) then ! if longitude is 0 in model index, then take all longitudes !!

              do k=1,size_lon
                aires_zonesormen(i) = aires_zonesormen(i) + aires_ocean(k,tab_zonesormen(i,j,2))
              enddo
            else
              aires_zonesormen(i)=aires_zonesormen(i) + aires_ocean(tab_zonesormen(i,j,1),tab_zonesormen(i,j,2))
            endif
          enddo ! sur j

       enddo ! sur i

       print*, "... fini !"

       end subroutine aires_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [Routine or Function Name Here]
!
!>     @brief This subroutine / function is adding a wonderful knew unused functionality
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine lect_fwf_init(t_debut,t_fin)


        integer, intent(in) :: t_debut, t_fin
        integer, parameter ::  tim_max = 10000

        integer :: IOstatus, IOstatline, eff_nbzones, indx_tim, i, indx_debut, indx_fin, taille_tim
        integer, dimension(nb_zonesmax) :: zones_nb
        integer, dimension(tim_max)     :: curr_date
        character(len=256) :: ligne

        real(dp), dimension(tim_max,nb_zonesmax) :: local_fwf

        integer(kind=ip):: fwf_input_dat_id


        open(newunit=fwf_input_dat_id,file="inputdata/clio/fwf_input.dat")

        zones_nb(:) = 0

        read(fwf_input_dat_id,*)
        read(fwf_input_dat_id,'(a)') ligne
        read(ligne,*,iostat=IOstatline) zones_nb(:)

        eff_nbzones = minloc(zones_nb,dim=1)-1

        read(fwf_input_dat_id,*)

        indx_tim = 0
        do while (IOstatus.NE.-1)
          indx_tim = indx_tim + 1
          read(fwf_input_dat_id,*,iostat=IOstatus) curr_date(indx_tim), local_fwf(indx_tim,1:eff_nbzones)
        enddo

        indx_tim = indx_tim - 1

        indx_debut = 1
        indx_fin   = indx_tim

        do i=1,indx_tim
          if (t_debut.ge.curr_date(i)) then
             indx_debut = i
             exit
          endif
        enddo

        if ((t_debut.ne.curr_date(indx_debut)).and.(indx_debut.gt.1)) then
          indx_debut = indx_debut -1
        endif

        write(*,*) "t_debut == ", t_debut, curr_date(indx_debut), indx_debut

        do i=1,indx_tim
          if (t_fin.ge.curr_date(i)) then
             indx_fin = i
             exit
          endif
        enddo

        write(*,*) "t_fin == ", t_fin, curr_date(indx_fin), indx_fin

        taille_tim = indx_fin - indx_debut + 1

        if (allocated(fwf_zonesormen)) then
          deallocate(fwf_zonesormen)
        endif

        allocate(fwf_zonesormen(1:taille_tim,1:eff_nbzones))

        fwf_zonesormen(:,:) = local_fwf(indx_debut:indx_fin,1:eff_nbzones)

        if (allocated(ID_zonesormen)) then
          deallocate(ID_zonesormen)
        endif

        allocate(ID_zonesormen(1:eff_nbzones))
        ID_zonesormen(:) = zones_nb(1:eff_nbzones)


        if (allocated(datefwf_zonesormen)) then
          deallocate(datefwf_zonesormen)
        endif

        allocate(datefwf_zonesormen(1:taille_tim))

        datefwf_zonesormen(:) = curr_date(indx_debut:indx_fin)

        close(fwf_input_dat_id)

       end subroutine lect_fwf_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [Routine or Function Name Here]
!
!>     @brief This subroutine / function is adding a wonderful knew unused functionality
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine update_fwf(dateBP,fwf_clio)

        integer, intent(in)                 :: dateBP
        real(dp),dimension(:,:),intent(out) :: fwf_clio ! en m.s-1

         call getfwfzones(dateBP)
         call fwf_zones2clio(fwf_clio)

       end subroutine update_fwf



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [Routine or Function Name Here]
!
!>     @brief This subroutine / function is adding a wonderful knew unused functionality
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine getfwfzones(dateloc)

        integer, intent(in)                 :: dateloc

        integer :: timemax, indx_loc

        timemax=ubound(datefwf_zonesormen,dim=1)

!~         write(*,*) "Dateloc in getfwfzones", dateloc

        if (dateloc.gt.datefwf_zonesormen(1)) then
           indx_loc = 1
           fwfzones(:)=fwf_zonesormen(indx_loc,:)

        elseif (dateloc.lt.datefwf_zonesormen(timemax)) then
           indx_loc = timemax
           fwfzones(:)=fwf_zonesormen(indx_loc,:)

        else
!afq&dmr --- On cherche la date la plus proche dans le fichier de forcage
           indx_loc=minloc(abs(datefwf_zonesormen(:)-dateloc),dim=1)
           if ((datefwf_zonesormen(indx_loc)-dateloc)*(datefwf_zonesormen(indx_loc+1)-dateloc).GT.0) then
              indx_loc = indx_loc-1
           endif

           fwfzones(:)=                                                              &
            (fwf_zonesormen(indx_loc,:)-fwf_zonesormen(indx_loc+1,:))/(datefwf_zonesormen(indx_loc)-datefwf_zonesormen(indx_loc+1))&
           * dateloc + fwf_zonesormen(indx_loc,:)                                                                                  &
           -(fwf_zonesormen(indx_loc,:)-fwf_zonesormen(indx_loc+1,:))/(datefwf_zonesormen(indx_loc)-datefwf_zonesormen(indx_loc+1))&
           * datefwf_zonesormen(indx_loc)

        endif

!~         write(*,*) "dans update", fwf_zonesormen(indx_loc,1), fwf_zonesormen(indx_loc+1,1),fwfzones(1),  &
!~          datefwf_zonesormen(indx_loc), datefwf_zonesormen(indx_loc+1),dateloc



       end subroutine getfwfzones


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [Routine or Function Name Here]
!
!>     @brief This subroutine / function is adding a wonderful knew unused functionality
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine fwf_zones2clio(cliofwf)

!dmr --- Notre flux additionel : variable cliofwfsv, [sv] = Sverdrups
        real(dp),dimension(:,:),intent(out) :: cliofwf

        integer :: nb_zoneseff, ii, jj, kk, n_zone, nb_points

        nb_zoneseff = ubound(fwf_zonesormen,dim=2)

        cliofwf(:,:) = 0.0_dp

        do n_zone = 1, nb_zoneseff

          do nb_points=1,tab_zonesormen(ID_zonesormen(n_zone),0,0)

              ii=tab_zonesormen(ID_zonesormen(n_zone),nb_points,1)
              jj=tab_zonesormen(ID_zonesormen(n_zone),nb_points,2)

              if (ii.eq.0) then
                do kk=1,size_lon
                  cliofwf(kk,jj)=cliofwf(kk,jj)+fwfzones(n_zone)/aires_zonesormen(ID_zonesormen(n_zone))
                enddo
              else
                cliofwf(ii,jj)=cliofwf(ii,jj)+fwfzones(n_zone)/aires_zonesormen(ID_zonesormen(n_zone))
              endif

          enddo
        enddo

        cliofwf(:,:) = cliofwf(:,:) * 1.E6_dp ! conversion factor from Sverdrups.m-2 to m3.s-1.m-2 = m.s-1

       end subroutine fwf_zones2clio


       subroutine set_datesBP(t_deb, t_fin)

       integer, intent(in) :: t_deb, t_fin

       tps_deb = t_deb
       tps_fin = t_fin

       end subroutine set_datesBP

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module zonesOrmenfwf_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
