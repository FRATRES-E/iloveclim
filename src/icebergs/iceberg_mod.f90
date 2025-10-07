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
!      MODULE: iceberg_mod
!
!>     @author  Didier M. Roche (dmr), Pepijn Bakker (PB)
!
!>     @brief This module iceberg_mod is the FORTRAN 90 transcription of the previous include
!
!>     @date Creation date: June, 19th, 2018
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : dmr, PB
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module iceberg_mod

       use global_constants_mod, only: dp, ip, vol_mass_dens_ice &
                                       , days_year360d_i, days_year360d
                                       
       use para0_mod, only: imax, jmax, kmax

       implicit none
      
       public :: write_iceberg_restartdata
       public :: introduce_icebergs
       public :: iceb_list2clio         
       public :: rearrange_lmx_arrays
       public :: iceberg_forcing         
       public :: convert_flux_to_icebergs
       public :: read_positb
       public :: init_iceberg
       public :: read_iceberg_restartdata
       public :: iceb_circbounc
       public :: iceb_findocean
       public :: iceb_writetracks
       public :: vol_icb, dVol_icb, heat_icb, hiceb_class &
                 , wiceb_class, pond_icb_class, uiceb_class, viceb_class
       public :: it_icb, lmx, kiceb, iceberg_info_out_id &
                 ,forcing_icb, iceflux_in_snow, iceflux_in_calv &
                 ,wiceb, hiceb, uiceb, viceb, pond_icb, xn, yn &
                 ,icbmax, fw_icb, numclass, uiceb_reg, viceb_reg
                
                                 
       private


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! --- Set of constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer(kind=ip), parameter, public:: dticeb = 86400    ! iceberg timestep (seconds / day)
      real(kind=dp), parameter, public :: repuls = 0.03    &! repulsion speed (m/s) of stranded icebergs
                                       , remsim = 0.1494  &! formulate in terms of density differences??
                                       , cdragw = 0.9     &! drag coefficient water
                                       , cdraga = 1.3     &! drag coefficient air      
                                       , Awater = 1.0     &! effective area for waveradiation
                                       , Aair = 1.768     &! effective area for waveradiation    
                                       , d_edge = 0.000001 &! Repulsed bergs are put at fron, d_edge degrees from the grid edge
                                       , d_cotes = 0.7    &! minimal distance from the coast
                                       , ti = 269.15      &!-4.+tkc; ti is iceberg temperature [NOTA] -4C, is this not minimal oceanic temperature ?
                                       , fbcf = 0.58/86400&! basal turbulent melting rate coefficient
                                       , flcf1 = 7.62*0.001/86400 &! buoyant convection melting rate coefficient   
                                       , flcf2 = 1.29*0.001/86400 &! buoyant convection melting rate coefficient
                                       , flcf3 = 0.5/86400 &! wave erosion melting rate
                                       , cfst1 = 0.92     &! coefficients to define if berg can roll over
                                       , cfst2 = 58.32    &! coefficients to define if berg can roll over
                                       , cdragi= 0.9      & ! sea-ice drag coefficient   
                                       , Cfusheaticb = 334000.0 ! constant to convert ice volume to heat flux (thersf.f)
                                       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer(kind=ip), parameter     :: icbmax = 500000  &! Maximum number of icebergs ... could be moved to an allocation procedure to avoid allocating monster-arrays
                                       , linitmax = 10000  &! initial maximum ? Used to define a few arrays below ... initial set of variables? Not sure it is useful on the long run
                                       , numclass = 10      & ! number of iceberg classes parameter in convert_flux_to_icebergs subroutine. Note: if different from default value of 10, the percentage classes need to be defined differently
                                       , skew_prct_class = 1 ! parameter in convert_flux_to_icebergs subroutine to skew distribution towards bigger or smaller bergs:
! 0 only small bergs (classes 1-3)
! 1 default
! 2 mainly large bergs (classes 8-10)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! --- Geographical-related variables (linked to the CLIO spatial grid), should be kept at minimum
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=dp), dimension(imax,jmax), public :: urep, vrep  &  ! repulsion speed, computed from repuls above
                                                  , wind10_u, wind10_v   ! 10 meters wind speed, gathered from ec_co2oc

      real(kind=dp), dimension(imax,jmax,kmax):: dVol_icb = 0.0_dp        ! dVol_icb is the primary output variable of the iceberg model

      real(kind=dp), dimension(imax,jmax)   :: iceflux_in_snow &         ! ice volume flux input on CLIO grid (m3) from excess snow; note this flux is water! used in convert_flux_to_icebergs subroutine
                                            ,  iceflux_in_calv           ! ice volume flux input on CLIO grid (m3) from calving in grisli; note this flux is ice!
                                    
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! --- Variables to define initial set of icebergs when using iceberg armada
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer(kind=ip), dimension(linitmax) :: pond_icb0
      real(kind=dp),   dimension(linitmax)  :: xn0, yn0, hiceb0, wiceb0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! --- Variables to define seasonal cycle of icebergs when using iceberg production from excess snow or from ISM
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer(kind=ip), dimension(linitmax):: pond_icb_seas
      real(kind=dp), dimension(linitmax)   :: xn_seas, yn_seas, hiceb_seas, wiceb_seas
      integer(kind=ip)                     :: num_icb_seas ! number of icebergs for this year
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! --- Arrays on the number of icebergs (individual characteristics of bergs)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=dp), dimension(icbmax)      :: uiceb, viceb &
                                            ,  xn, yn, xn_old, yn_old, hiceb, wiceb &! hiceb is submarine heigth of iceberg
                                            ,  vol_orig, vol_melt &
                                            ,  uiceb_reg, viceb_reg
      integer(kind=ip), dimension(icbmax)   :: kiceb, pond_icb                                                     


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! --- Restart and armada forcing
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      

      integer(kind=ip)                      :: linit  &                   ! Number of icebergs read from positb.init 
                                            ,  restart_icb, forcing_icb   ! Flags in positb.init to state if cold restart is to be used and if 'armada of icebergs' forcing is to be used
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! --- convert_flux_to_icebergs
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      

      real(kind=dp), dimension(numclass)    :: prct_class & ! percentage of iceberg mass per class
                                            ,  hiceb_c    & ! iceberg height per class (m)
                                            ,  wiceb_c    & ! iceberg width per class (m)
                                            ,  pond_icb_c   ! number of icebergs per class (-)
      real(kind=dp)                            hiceb_c_orig ! keep a copy of the original hiceb_c(1)
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! --- Clio-shaped arrays for outputing to netcdf
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      

      real(kind=dp), dimension(imax,jmax,kmax):: fw_icb, heat_icb  ! Freshwater flux [m/time] and heat flux [K/time] due to iceberg melting
      
      real(kind=dp), dimension(imax,jmax):: vol_icb      ! Total volume in gridcell
                                                 
! Icebergs on CLIO grid, averages per width-class (10 in total)
      real(kind=dp), dimension(imax,jmax,numclass):: hiceb_class &    ! iceberg height per width-class
                                            , wiceb_class &    ! iceberg width  per width-class
                                            , pond_icb_class & ! number of icebergs  per width-class
                                            , uiceb_class,viceb_class ! zonal & meriodonal iceberg velocity  per width-class
 

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! --- id to open files
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      
      integer(kind=ip) :: iceberg_info_out_id
      integer(kind=ip) :: iceberg_tracksoutput_id
                                              
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! --- Other
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      
 
      integer(kind=ip)                      :: it_icb  &                 ! Count iterations of iceberg module                 
                                            ,  lmx                       ! Number (label) of last iceberg in table
        
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! 
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! 
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      





      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: init_iceberg
!
!      DESCRIPTION:
!
!>     Initialization of various iceberg module variables.
!      Also the repulsion velocity is calculated here.
!      Initialize the arrays used in convert_flux_to_icebergs subroutine: percentage of ice volume per iceberg size class, and default widht and height of icebergs per class
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine init_iceberg

        use bloc0_mod, only: ks2, tms
        use united_Grid_mod, only: united_Grid_init => init_module

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(kind=ip) :: i,j,l
        logical          :: itexists
                                              
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Open new iceberg module output text files.                                          |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        inquire(file='outputdata/icebergs/iceberg_info.out', EXIST=itexists) 
        if ( itexists ) then
          open(newunit=iceberg_info_out_id,file='outputdata/icebergs/iceberg_info.out',status='replace')
        else
          open(newunit=iceberg_info_out_id,file='outputdata/icebergs/iceberg_info.out',status='new')
        endif
        
        inquire(file='outputdata/icebergs/iceberg_tracks.txt', EXIST=itexists) 
        if ( itexists ) then
          open(newunit=iceberg_tracksoutput_id,file='outputdata/icebergs/iceberg_tracks.txt',status='replace')
        else
          open(newunit=iceberg_tracksoutput_id,file='outputdata/icebergs/iceberg_tracks.txt',status='new')
        endif

        inquire(file = "outputdata/icebergs/iceberg_info.out", EXIST=itexists)
        if ( itexists ) then
          open(newunit=iceberg_info_out_id,file='outputdata/icebergs/iceberg_info.out',status='replace')
        else
          open(newunit=iceberg_info_out_id,file='outputdata/icebergs/iceberg_info.out',status='new')
        end if
                       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Some Initialisations of the run.                                          |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        write(iceberg_info_out_id,*) '==================================='
        write(iceberg_info_out_id,*) '===Initializing iceberg module====='
#if ( ICEBERG == 1 )
        write(iceberg_info_out_id,*) 'ICEBERG == 1 -> no iceberg fluxes from land or icesheet model'
#elif ( ICEBERG == 2 && ISM != 2 )
        write(iceberg_info_out_id,*) 'ICEBERG == 2 && ISM != 2 -> Icebergs created from excess snow'
#elif ( ICEBERG == 2 && ISM == 2)
        write(iceberg_info_out_id,*) 'ICEBERG == 2 && ISM == 2 -> Icebergs created from the icesheet model ice flux'
#endif        
        it_icb=0                ! seems to define the time = number of days. Could be deprecated in favor of a more generic timer

        num_icb_seas = 0 ! number of icebergs produced for upcoming year.      

        dVol_icb = 0 ! melt volume that will be pased to the ocean every year
        
        do l=1,icbmax
          uiceb(l)=0.
          viceb(l)=0.
          uiceb_reg(l)=0.
          viceb_reg(l)=0.          
          xn(l)=0.
          yn(l)=0.
          kiceb(l)=0
          hiceb(l)=0.
          wiceb(l)=0.
          pond_icb(l)=0
        end do

! initialization of 'old' position to zero
        xn_old=0.
        yn_old=0.
       
! initialization of repulsion velocity                                         |

        urep=0
        vrep=0
        do j=2,jmax
          do i=2,imax
            urep(i,j)=repuls*(tms(i,j,ks2) -tms(i-1,j-1,ks2) &
                      +tms(i,j-1,ks2)-tms(i-1,j,ks2))
            vrep(i,j)=repuls*(tms(i,j,ks2) -tms(i-1,j-1,ks2) &
                      +tms(i-1,j,ks2)-tms(i,j-1,ks2))
          end do
        end do

!-------------------------------------------------------
! Initialize arrays used in convert_flux_to_icebergs subroutine
! Percentage of mass depending on the size class of the berg
      if (skew_prct_class.eq.0) then
        prct_class(1) = 0.333
        prct_class(2) = 0.333
        prct_class(3) = 0.333
        prct_class(4) = 0.
        prct_class(5) = 0.
        prct_class(6) = 0.
        prct_class(7) = 0.
        prct_class(8) = 0.
        prct_class(9) = 0.
        prct_class(10) =0.
        write(iceberg_info_out_id,*)'Using only small icebergs: ',skew_prct_class
      else if (skew_prct_class.eq.1) then
        prct_class(1) = 0.15
        prct_class(2) = 0.15
        prct_class(3) = 0.2
        prct_class(4) = 0.15
        prct_class(5) = 0.08
        prct_class(6) = 0.07
        prct_class(7) = 0.05
        prct_class(8) = 0.05
        prct_class(9) = 0.05
        prct_class(10) =0.05 
        write(iceberg_info_out_id,*)'Using default iceberg distribution: ',skew_prct_class              
      else if (skew_prct_class.eq.2) then
        prct_class(1) = 0.
        prct_class(2) = 0.
        prct_class(3) = 0.
        prct_class(4) = 0.
        prct_class(5) = 0.
        prct_class(6) = 0.
        prct_class(7) = 0.
        prct_class(8) = 0.333
        prct_class(9) = 0.333
        prct_class(10) =0.333     
        write(iceberg_info_out_id,*)'Using only large icebergs: ',skew_prct_class
      end if
           
! height of bergs depending on classes (Note: hiceb and hiceb_c are submarine heigth of icebergs)
      hiceb_c(1) = 67
      hiceb_c(2) = 133
      hiceb_c(3) = 200
      hiceb_c(4) = 267
      hiceb_c(5) = 300
      hiceb_c(6) = 300
      hiceb_c(7) = 300
      hiceb_c(8) = 300
      hiceb_c(9) = 300
      hiceb_c(10)= 300

! Keep a copy of the original hiceb_c(1) because it will be adjusted in convert_flux_to_icebergs
      hiceb_c_orig = hiceb_c(1)
      
! width of bergs depending on classes (also used to create CLIO output of mean values in iceberg classes)
      wiceb_c(1) = 67
      wiceb_c(2) = 133
      wiceb_c(3) = 200
      wiceb_c(4) = 267
      wiceb_c(5) = 333
      wiceb_c(6) = 400
      wiceb_c(7) = 500
      wiceb_c(8) = 600
      wiceb_c(9) = 800
      wiceb_c(10)= 1000
      
!-------------------------------------------------------
! Check if sum of percentages given for input classes equals 1
      if ( .not. sum(prct_class).gt.0.99) then
        write(iceberg_info_out_id,*)'Sum of iceberg class percentages .ne. 1!',sum(prct_class)
      end if     


      call united_Grid_init()

        write(iceberg_info_out_id,*) '==================================='

      return
      end subroutine init_iceberg

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: read_iceberg_restartdata
!
!      DESCRIPTION:
!
!>     Read iceberg restart data from resicb.om file
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine read_iceberg_restartdata
     
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(kind=ip)          :: l, resicb_om_id
        logical                   :: rexist=.FALSE.
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! Get restart data, depending on restart_icb and on CALVFLUX flag
        if (restart_icb == 0) then ! use restart file

        write(iceberg_info_out_id,*)   '==================================='
          write(iceberg_info_out_id,*) '========Reading restart ==========='
          
          inquire(file='startdata/resicb.om',exist=rexist)
          if (rexist) then
            write(iceberg_info_out_id,*) 'Using resicb.om restart file for icebergs'
            open(newunit=resicb_om_id,file='startdata/resicb.om',status='old' &
                        ,form='unformatted')
            ! Read information on current icebergs
            read (resicb_om_id) lmx
            do l=1,lmx
              read(resicb_om_id) xn(l),yn(l),kiceb(l),uiceb(l) &
                       ,viceb(l),hiceb(l),wiceb(l),pond_icb(l)
            end do
            ! Read information on iceberg seasonal cycle for next year
            read (resicb_om_id) num_icb_seas
            do l=1,num_icb_seas
            read(resicb_om_id) pond_icb_seas(l),xn_seas(l), yn_seas(l) &
                       , hiceb_seas(l), wiceb_seas(l)
            end do        
            close(resicb_om_id)  
                              
          else ! When restart is called for but not found 
            write(iceberg_info_out_id,*) 'Iceberg module: resicb.om not found'
            lmx = 0
          endif
          
        else ! restart_icb == 1 (cold restart)
          write(iceberg_info_out_id,*) 'Cold restart for icebergs'
          lmx = 0
          num_icb_seas = 0
        endif
        
        write(iceberg_info_out_id,*) '==================================='

        return
      end subroutine read_iceberg_restartdata

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: write_iceberg_restartdata
!
!>     Write iceberg restart data to resicb.om file
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine write_iceberg_restartdata

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(kind=ip)          :: l, nbwires, resicb_om_id

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        write(iceberg_info_out_id,*) '==================================='
        write(iceberg_info_out_id,*) '===Writing iceberg restart data===='
        
!Write restart file resicb.om
        nbwires = 0
        write(iceberg_info_out_id,*)'Writing restart file resicb.om'        
        open(newunit=resicb_om_id,file='restartdata/resicb.om',status='replace' &
                      ,form='unformatted')
!PB How to get the restart file in the proper place?     
        do l=1,lmx
          if (kiceb(l) .ge. -1) nbwires=nbwires+1 ! Count number of icebergs to write to restart file
        end do
        write(resicb_om_id) nbwires
!        write(iceberg_info_out_id,*) 'test write nbwires',nbwires
          
! PB Write current iceberg information          
        do l=1,lmx                 
          if(kiceb(l) .ge. -1) then
            write(resicb_om_id) xn(l),yn(l),kiceb(l),uiceb(l) &
                      ,viceb(l),hiceb(l),wiceb(l),pond_icb(l)
          endif
        end do

! Write information on iceberg seasonal cycle for next year
        write (resicb_om_id) num_icb_seas
        do l=1,num_icb_seas
          write(resicb_om_id) pond_icb_seas(l),xn_seas(l), yn_seas(l) &
                       , hiceb_seas(l), wiceb_seas(l)
        end do 
        close(resicb_om_id)

        write(iceberg_info_out_id,*) '==================================='
        return
      end subroutine write_iceberg_restartdata

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: read_positb
!
!>     Read the iceberg parameter file positb.init.
!      Including restart_icb, forcing_icb, linit
!      If iceberg model is forced by Iceberg Armada (forcing_icb == 1) then these data are also read
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine read_positb

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        logical                   :: rexist=.FALSE.
        integer(kind=ip)          :: l         
        integer(kind=ip)          :: positb_init_id

         
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        write(iceberg_info_out_id,*) '==================================='
        write(iceberg_info_out_id,*) '========Reading positb.init========'
        
        inquire(file='inputdata/positb.init',exist=rexist)
        if (rexist) then
          open (newunit=positb_init_id, file="inputdata/positb.init", status='old', form='formatted')
          write(iceberg_info_out_id,*) 'Read positb'
          read(positb_init_id,*)
          read(positb_init_id,*) restart_icb    
          read(positb_init_id,*) 
          read (positb_init_id,*) forcing_icb
          read (positb_init_id,*)
          read (positb_init_id,*) linit
          read (positb_init_id,*)
        else 
          write(iceberg_info_out_id,*) 'positb.init file is missing'
        endif

!-------------------------------------------------------
! Read positb.init forcing of iceberg Armada if called for
        if (forcing_icb == 1) then

! JONO APW instead of reading iceberg variables xn(l),..,... from positb,
! we put these vars into basic production vars xn0(1-linit),..0,..
! READ 30=positb.init icebergs
! (pond_icb is a factor to simulate a number of bergs as one particle)

! linit: number of bergs in positb.init

          write(iceberg_info_out_id,*)'Iceberg Armada read from positb.init'

          do l=1,linit
            read(positb_init_id,*) xn0(l),yn0(l),hiceb0(l),wiceb0(l),pond_icb0(l)
            write(iceberg_info_out_id,*) 'Iceberg forcing',&
            xn0(l), yn0(l), hiceb0(l),wiceb0(l),pond_icb0(l)
          end do

        endif !(forcing_icb == 1)

        write(iceberg_info_out_id,*) '==================================='

        return
      end subroutine read_positb

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: iceberg_forcing
!
!>     Impose the iceberg forcing that is read in iceberg parameter file positb.init.
!>     Check if new iceberg is on land or of the grid, if so it is ignored
!>     The new iceberg are positioned at the first empty spot in the iceberg table
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine iceberg_forcing
      
         use bloc0_mod, only: tms, zw, ks2
         use to_and_from_CLIO, only: get_indexes_C, get_lonlat_C

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(kind=ip)          :: kk,l,ll, empty_spot, indexes_C(2)      
          
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        write(iceberg_info_out_id,*) '==================================='
        write(iceberg_info_out_id,*) '====Imposing Armade of icebergs===='
        write(iceberg_info_out_id,*) 'Total of ',linit,' iceberg groups'

! Add new icebergs to empty spots in iceberg table
        do l=1,linit  ! loop over new icebergs
          if (pond_icb0(l).ne.0) then
            indexes_C = get_indexes_C(xn0(l), yn0(l))
            if (tms(indexes_C(1),indexes_C(2),ks2).eq.0) then
              write(iceberg_info_out_id,*) 'Iceberg from positb.init is on land, removed', &
                                           l,xn0(l),yn0(l), &
                                           tms(indexes_C(1),indexes_C(2),ks2)
            else
! Find empty spot in iceberg table
              empty_spot=0
              do ll=1,lmx
                if (kiceb(ll).eq.(-2)) then
                  empty_spot=ll
                  goto 202
                end if
              end do
 202   continue
 
! If no empty spot, add icebergs to the end of the tabel
              if (empty_spot.eq.0) then
                lmx=lmx+1
                empty_spot=lmx
              end if
         
! Fill empty spot
              xn(empty_spot)=xn0(l)
              yn(empty_spot)=yn0(l)
              wiceb(empty_spot)=wiceb0(l)
              pond_icb(empty_spot)=pond_icb0(l)
              uiceb(empty_spot)=0.
              viceb(empty_spot)=0.
              hiceb(empty_spot)=hiceb0(l)              
                       
!----------------------------------------------------------------------
! Calculate water-depth of the iceberg in k-levels
              kiceb(empty_spot)=ks2
              do kk=ks2,2,-1
                if (-zw(kiceb(empty_spot)).le.hiceb(empty_spot)) then
                  kiceb(empty_spot)=kk-1
                else
                  goto 302
                end if
              end do
 302    continue
 
            if (kiceb(empty_spot).le.-1) kiceb(empty_spot)=ks2
                 
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
            end if ! if (tms(i,j,ks2).eq.0)
          end if ! if (pond_icb0(l).ne.0)
        end do ! do l=1,linit

        write(iceberg_info_out_id,*) '==================================='

        return
      end subroutine iceberg_forcing

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|




!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: convert_flux_to_icebergs
!
!>     @Routine that describes how a volume flux of ice should be converted into iceberg with certain characteristics
!
!      DESCRIPTION:
!
!      With a double loop, loop over icevolume input file.
!      First put ice volume in iceberg according to percentages and default height and width.
!      If ice remains, create additional icebergs in the classes that already have some icebergs, starting with the second to largest class.
!      Finally, still some ice could remain and the budget is closed by adjusting the height of the smallest iceberg class. Because the number of bergs in this class is generally large, the height adjustment is rather small
!      Then the new icebergs are added to the iceberg table
!      The results are all icebergs for the whole year put in arrays ending with *_seas and those can be used later to actually create icebergs per ocean timestep
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine convert_flux_to_icebergs(iceflux_in)

      use ice_mod, only: imax, jmax
      use to_and_from_CLIO, only: get_indexes_C, get_lonlat_C
      use bloc0_mod, only: tms, ks2 !zw

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(kind=ip)                       :: i,j,k,l,ll &
                                               ,  empty_spot
!                                              ,  indexes_C(2)   
        real(kind=dp)                          :: vol_singleberg_k &
                                               ,  vol_tot,vol_remain &
                                               ,  hiceb_adjusted &
                                               ,  pond_icb_temp &
                                               ,  iceb_toti
        real(kind=dp), dimension(2)            :: lon_lat_C 
        real(kind=dp), dimension(imax,jmax)    :: iceflux_in
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!    Initialize arrays
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        iceb_toti = 0.
        pond_icb_seas(:) = 0.
        xn_seas(:) = 0.
        yn_seas(:) = 0.
        hiceb_seas(:) = 0.
        wiceb_seas(:) = 0.
        empty_spot = 1
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!    Main loop
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! Loop over iceflux_in grid to find locations of new icebergs
        do j=2,jmax-1
          do i=2,imax-1
            hiceb_c(1) = hiceb_c_orig ! set hiceb_c(1) back to original value
            if (iceflux_in(i,j).gt.0) then ! new icebergs found

              iceb_toti = iceb_toti + iceflux_in(i,j)
!----------------------------------------------------------------------
! Put ice volume in all classes according to percentage in that class and the height and width of that class.
              vol_tot=0.
              pond_icb_c(:)=0.
              do k=1,numclass ! loop over number of iceberg classes
                if (prct_class(k).ne.0) then
                  vol_singleberg_k = (1+remsim)*hiceb_c(k)*wiceb_c(k)*1.5*wiceb_c(k)
                  pond_icb_c(k)= floor(iceflux_in(i,j)*prct_class(k) &
                                / vol_singleberg_k)                                
                  if (pond_icb_c(k).le.0) goto 400
                  vol_tot = vol_tot+(vol_singleberg_k*pond_icb_c(k))
                end if
              end do
 400          continue
 
! Put remaining ice volume in the iceberg classes for which pond_icb_k is non-zero, starting from largest
              do ll=k-1,1,-1
                if (prct_class(ll).ne.0) then
                  vol_remain = iceflux_in(i,j)-vol_tot
                  vol_singleberg_k = (1+remsim)*hiceb_c(ll)*wiceb_c(ll) &
                                     *1.5*wiceb_c(ll)
                  pond_icb_temp= floor(vol_remain / vol_singleberg_k)
                  if (pond_icb_temp.le.0) goto 500                
                  vol_tot = vol_tot+(vol_singleberg_k*pond_icb_temp)
                  pond_icb_c(ll)= pond_icb_c(ll)+ pond_icb_temp  
!                  write(iceberg_info_out_id,*) 'test3',vol_remain,vol_singleberg_k  &
!                  ,pond_icb_temp,pond_icb_c(ll),ll

                end if
              end do  
 500          continue    
                      
! Adjust iceberg height in smallest class to close ice volume budget
              vol_singleberg_k = (1+remsim)*hiceb_c(1)*wiceb_c(1)*1.5 &
                                *wiceb_c(1)
              vol_remain = iceflux_in(i,j)-vol_tot+(vol_singleberg_k &
                          *pond_icb_c(1)) ! remaining volume after filling up class 1 with default height
              pond_icb_c(1)= ceiling(vol_remain / vol_singleberg_k)             
              hiceb_adjusted = vol_remain / ((1+remsim)*wiceb_c(1)*1.5*wiceb_c(1) &
                              *pond_icb_c(1))
              hiceb_c(1) = hiceb_adjusted
                  
!----------------------------------------------------------------------
! Check iceberg volume input and output; if difference is larger than 0.1% give error message
              vol_tot=0.
              do k=1,numclass
                vol_tot = vol_tot+((1+remsim)*hiceb_c(k)*wiceb_c(k)*1.5 &
                         *wiceb_c(k)*pond_icb_c(k))
              end do
              if (abs(vol_tot-iceflux_in(i,j))/iceflux_in(i,j).gt.0.001) then
                write(iceberg_info_out_id,*) 'Ice volume input and output differ, something must be going wrong!'
                write(iceberg_info_out_id,*) vol_tot-iceflux_in(i,j),iceflux_in(i,j)
                write(iceberg_info_out_id,*) abs(vol_tot-iceflux_in(i,j))/iceflux_in(i,j)*100,'% difference'
              end if

!----------------------------------------------------------------------
! Add new icebergs to empty spots in seasonal iceberg table
              lon_lat_C = get_lonlat_C(i,j)
              if (tms(i,j,ks2).eq.0) then
                write(iceberg_info_out_id,*) 'Iceberg from excess snow is on land, not used', &
                       l,i,j,lon_lat_C(1),lon_lat_C(2), tms(i,j,ks2)
              else
                do k=1,numclass ! loop over number of iceberg classes
                  if (pond_icb_c(k).ne.0) then
                    xn_seas(empty_spot)=lon_lat_C(1)
                    yn_seas(empty_spot)=lon_lat_C(2)
                    wiceb_seas(empty_spot)=wiceb_c(k)
                    hiceb_seas(empty_spot)=hiceb_c(k)                  
                    pond_icb_seas(empty_spot)=pond_icb_c(k)
                    empty_spot = empty_spot + 1
                  end if
                end do        
              endif    
!----------------------------------------------------------------------
            end if !(iceflux_in(i,j).lt.0)
          end do !j=1,jmax
        end do ! i=1,imax
        
        num_icb_seas = empty_spot-1 ! From empty spot to total number of new icebergs
        
        write(iceberg_info_out_id,*) "Total created new icebergs for this year",num_icb_seas

        return        
      end subroutine convert_flux_to_icebergs

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: introduce_icebergs
!
!      DESCRIPTION:
!
!      Distribute the created yearly sum of icebergs equally over the year
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine introduce_icebergs(istep)

        use bloc0_mod, only: ks2, zw

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(kind=ip) :: istep, day, l, kk, n, empty_spot &
                          , nmin, nmax, num_day_icb, num_rest_icb
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        if (num_icb_seas.gt.0) then
          ! Divide total number of icebergs of this year evenly over the year        
          num_day_icb = max(floor(num_icb_seas/days_year360d),1)
          num_rest_icb = num_icb_seas-(num_day_icb*days_year360d_i)
          day = mod(istep,days_year360d_i)

          nmin = (day-1)*num_day_icb+1
          nmax = day*num_day_icb
          ! Introduce the remaining icebergs at last day of the year
          if (day.eq.0) then
            nmin = num_icb_seas-(num_day_icb+num_rest_icb)+1
            nmax = num_icb_seas
          endif
          
          ! When number of icebergs is smaller than number of days in the year, put 1 iceberg per day until all are gone
          if (num_day_icb.eq.1) then
            if (day.le.num_icb_seas) then
              nmin = 1
              nmax = 1
            else
              nmin = 0
              nmax = 0
            endif
          endif
	         
! Add new icebergs to empty spots in iceberg table
          do n=nmin,nmax ! loop over the iceberg that are to be introduced on this day
            if (pond_icb_seas(n).ne.0) then
              empty_spot=0
! Find empty spot in iceberg table
              do l=1,lmx
                if (kiceb(l).eq.(-2)) then
                  empty_spot=l
                  goto 211
                endif
              end do
211           continue
! If no empty spot, add icebergs to the end of the tabel
              if (empty_spot.eq.0) then
                lmx=lmx+1
                empty_spot=lmx
              endif
         
! Fill empty spot
              xn(empty_spot)=xn_seas(n)
              yn(empty_spot)=yn_seas(n)
              wiceb(empty_spot)=wiceb_seas(n)
              pond_icb(empty_spot)=pond_icb_seas(n)
              uiceb(empty_spot)=0.
              viceb(empty_spot)=0.
              hiceb(empty_spot)=hiceb_seas(n)              
              
!----------------------------------------------------------------------
! Calculate water-depth of the iceberg in k-levels
              kiceb(empty_spot)=ks2
              do kk=ks2,2,-1
                if (-zw(kiceb(empty_spot)).le.hiceb(empty_spot)) then
                  kiceb(empty_spot)=kk-1
                else
                  goto 311
                endif
              end do
311           continue
              if (kiceb(empty_spot).le.-1) kiceb(empty_spot)=ks2
!----------------------------------------------------------------------
            endif !(pond_icb_c(k).ne.0)
          end do !n=nmin,nmax   
        endif !(num_icb_seas.gt.0)

        return        
      END SUBROUTINE introduce_icebergs

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|




!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: rearrange_lmx_arrays
!
!      DESCRIPTION:
!
!      This routine rearranges the table containing the information on the icebergs. 
!      It will remove all empty spots in order for the table not to become too long
!
!      Auteur : Didier M. Roche  & Marianne Buegelmayer
!      Date   : 23 mai 2012
!      Derniere modification :  31 mai 2012
!      30-6-2020 Pepijn Bakker: removed livingicb variable, not used
!      19-10-2021 Pepijn Bakker: added in and jn variables
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine rearrange_lmx_arrays

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer(kind=ip) :: licb, licb1
        
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       licb = 1
       licb1 = 0
       DO WHILE ((kiceb(licb).LE.-1).AND.(licb.LT.lmx))
           licb = licb + 1  
           licb1 = licb1 + 1  
       END DO

       IF (licb1.GE.1) THEN

          xn(1:lmx-licb1) = xn(licb:lmx)
          yn(1:lmx-licb1) = yn(licb:lmx)
          hiceb(1:lmx-licb1) = hiceb(licb:lmx)
          wiceb(1:lmx-licb1) = wiceb(licb:lmx)
          uiceb(1:lmx-licb1) = uiceb(licb:lmx)
          viceb(1:lmx-licb1) = viceb(licb:lmx)
          pond_icb(1:lmx-licb1) = pond_icb(licb:lmx)
          kiceb(1:lmx-licb1) = kiceb(licb:lmx)
          vol_orig(1:lmx-licb1) = vol_orig(licb:lmx)
          vol_melt(1:lmx-licb1) = vol_melt(licb:lmx)
          
          lmx = lmx-licb1

       END IF

! Set lmx to zero if only empty icebergs remain
       if ((lmx.eq.1).and.(kiceb(1).eq.-2)) then
          lmx = 0
          xn(1) = 0
          yn(1) = 0
          hiceb(1) = 0
          wiceb(1) = 0
          uiceb(1) = 0
          viceb(1) = 0
          pond_icb(1) = 0

       end if            

      return        
      END SUBROUTINE rearrange_lmx_arrays


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|




!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: iceb_list2clio
!
!>     @Routine that creates 2D variables on the CLIO grid from the list-type output of the iceberg module
!
!      DESCRIPTION:
!
!      For a limited number of iceberg variables, 2D arrays are created from the list-type iceberg output.
!      To dondence the list-type variables onto a 2D grid, the mean (or sum) is calculated for 10 different 'width-classes'
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine iceb_list2clio

        use para0_mod, only: imax, jmax
! i = call to function of Didier
        use to_and_from_CLIO, only: get_indexes_C

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     Variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        implicit none

        integer(kind=ip) :: j,i,l,c, indexes_C(2)
        integer(kind=ip), dimension(imax,jmax,numclass):: count_icb ! To keep track of number of icebergs per CLIO grid cell
        integer(kind=ip), dimension(1):: wclass ! note: needs to be an array of size 1 because the function 'maxloc' returns an array type variable

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the subroutine starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|      
! Initialize variables
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        vol_icb=0.
        hiceb_class=0.
        wiceb_class=0.
        pond_icb_class=0.
        uiceb_class=0.
        viceb_class=0.

        count_icb = 0
      
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Create iceberg output data      
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Loop over all icebergs in list lmx
        do l=1,lmx
          if (wiceb(l).gt.0.) then
            
!- Position of icebergs on clio grid
            indexes_C = get_indexes_C(xn(l), yn(l))
            i=indexes_C(1)
            j=indexes_C(2)
            
!- Find width-class of the iceberg (classes run from 0-wiceb_c(1); wiceb_c(1)-wiceb_c(2) etc.                              
            wclass = maxloc(wiceb(l)-wiceb_c,mask=(wiceb(l)-wiceb_c < 0))
            
!- Calculate total iceberg volume on CLIO grid
            vol_icb(i,j) = vol_icb(i,j)+1.5*wiceb(l)*wiceb(l) &
                *hiceb(l)*(1+remsim)*pond_icb(l)

!- Calculate sums per width-class on CLIO grid (divided by iceberg counts later to get to the mean)
            hiceb_class(i,j,wclass) = hiceb_class(i,j,wclass)+hiceb(l)
            wiceb_class(i,j,wclass) = wiceb_class(i,j,wclass)+wiceb(l)
            pond_icb_class(i,j,wclass) = pond_icb_class(i,j,wclass)+pond_icb(l) ! output as a sum, not a mean
            uiceb_class(i,j,wclass) = uiceb_class(i,j,wclass)+uiceb_reg(l) ! use iceberg velocities on the regular grid as output in netcdf
            viceb_class(i,j,wclass) = viceb_class(i,j,wclass)+viceb_reg(l) ! use iceberg velocities on the regular grid as output in netcdf
            
!- Count number of icebergs per CLIO grid cell and class
            count_icb(i,j,wclass)=count_icb(i,j,wclass)+1

          endif ! end if (wiceb(l).gt.0) loop
        end do ! end lmx loop
      
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Calculate means per width-class on CLIO grid by dividing sums by number of bergs
        do i=1,imax
          do j=1,jmax
            do c=1,numclass
               if (hiceb_class(i,j,c).ne.0) then
                  hiceb_class(i,j,c)=hiceb_class(i,j,c)/count_icb(i,j,c)
                  wiceb_class(i,j,c)=wiceb_class(i,j,c)/count_icb(i,j,c)
                  uiceb_class(i,j,c)=uiceb_class(i,j,c)/count_icb(i,j,c)
                  viceb_class(i,j,c)=viceb_class(i,j,c)/count_icb(i,j,c)
               endif
            end do
          end do
        end do

      return        
      END SUBROUTINE iceb_list2clio


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|





!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: iceb_writetracks
!
!>     @Routine that writes the iceberg number and position to a text file
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine iceb_writetracks(istep)


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     Variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        implicit none

        integer(kind=ip) :: l,istep
        logical          :: itexist
        
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the subroutine starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Create iceberg track output to text file    
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Loop over all icebergs in list lmx
        inquire(file='outputdata/icebergs/iceberg_tracks.txt', EXIST=itexist) 
        if ( itexist ) then
          open(newunit=iceberg_tracksoutput_id,file='outputdata/icebergs/iceberg_tracks.txt',&
               status='old',position="append")
        end if
       
        do l=1,lmx
          if (wiceb(l).gt.0.) then
            write(iceberg_tracksoutput_id,'(2I4,2F10.2)') istep, l, xn(l), yn(l)
           endif
         enddo
         close(iceberg_tracksoutput_id)
         
      return        
      END SUBROUTINE iceb_writetracks


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|




!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: iceb_circbounc
!
!>     @Routine that applies circular boundary condition
!     
!      DESCRIPTION:
!
!      If an iceberg is found to be off the grid in longitudinal direction, the circular boundary condition can be called and applied.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine iceb_circbounc(x, y)

        use global_constants_mod, only: pi_dp, deg_to_rad, rad_to_deg
        use to_and_from_CLIO, only: get_indexes_C
        
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     Variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        implicit none

        real(KIND=dp),    PARAMETER     :: dblePI       =   2.0_dp * pi_dp
        real(kind=dp),    intent(inout) :: x   
        real(kind=dp),    intent(in)    :: y   
        real(kind=dp)                   :: xnew   
        
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the subroutine starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----| 
!- Apply circular boundary condition if needed
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        xnew = mod(x*deg_to_rad+dblePI,dblePI)*rad_to_deg
!        if (.NOT. (abs(xnew/x - 1) < 1e-4)) then
!          write(iceberg_info_out_id,*) 'circ bound',x,xnew,y
!        end if
        x = xnew
        
        
      return
      end subroutine iceb_circbounc


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


   

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: iceb_findocean
!
!>     @Routine A close-by gridcell is searched for that is on the grid and ocean
!     
!      DESCRIPTION:
!
!      This subroutine should be called in case the iceberg gridcell under consideration is either of the grid in the latitudinal direction or is land.
!      The search routine is fairly simple, it doesn't look really for the closest gridcell, but makes squares around the current gridcell in the clock-wise direction, starting from the bottom-left. If needed, the search radius is increaded until a maximum of 'inc_max' to all sides
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine iceb_findocean(l, i, j, kiceb)

        use para0_mod, only: jmax
        use bloc0_mod, only: tms, ks2  

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     Variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        implicit none

        integer(kind=ip), intent(in) :: l
        integer(kind=ip), intent(inout):: i, j, kiceb

        integer(kind=ip) :: ii, jj, inc
        integer(kind=ip) :: inc_max = 20
  
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the subroutine starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----| 
!- Find location that is both on the grid and is an ocean cell. 
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
          do inc=1,inc_max
            do ii=i-inc,i+inc
              do jj=j-inc,j+inc
                if (jj.le.1.or.jj.ge.imax.or.tms(ii,jj,ks2).eq.0) then
                else
                  i=ii
                  j=jj
                  goto 100
                endif
              end do
            end do
          end do
          write(iceberg_info_out_id,*) &
            'Error, no good ocean cell found, iceberg removed',l
          kiceb = 0

 100  continue 
      
      return
      end subroutine iceb_findocean


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
            
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


end module iceberg_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
