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


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: path_mod
!
!>     @author  Lise Missiaen (lim)
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module [path_mod] is computing for 231Pa and 230Th :
!                                                                      Production of dissolved Pa and Th in the water column
!                                                                      Radioactive decay
!                                                                      Reversible scavenging : exchange of activity between particulate and dissolved phases
!                                                                      Settling, sedimentation of particulate activities
!
!>     @date Creation date: March, 20th, 2017
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : lim
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!     Based on what has been developed for Bern3D model, NEMO-PISCES...
!!     Useful references : Marchal, 2000, Sidall 2005 and 2007, Dutay 2009
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

    module path_mod

       use global_constants_mod, only: dp, str_len

       implicit none
      ! private

       private :: si_calc
       public  :: path_init
       public  :: path_step
       public  :: path_write_rest
       
      ! lim test add subroutine to read fixed particles fields from NEMO
      ! TODO public or private ??? lim : private ??
       private  :: read_particles_fields
       private  ::  transform_particles_fields


       public   :: restart_fnm


       ! NOTE_avoid_public_variables_if_possible

! --- lim&dmr Following variables define the position of path variables wihin the CLIO scal array
       integer, parameter, public :: scalstart_path = 14, scalend_path = 17

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim  Tracers and index
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       
       integer, parameter :: imax_loc = 122, jmax_loc = 65, kmax_loc = 20 ! dimensions CLIO a importer ...!!!! Attention orientation surface-bottom... !!!

       integer, dimension(imax_loc, jmax_loc), public :: k_surf, k_fond ! indexes of the surface of the ocean

       integer, parameter :: nb_trac  = 4 ! (pa, th)*(part,dissous)
       integer, parameter :: npad = 1, npap = 2, nthd = 3, nthp=4 ! pa = 231Pa, th = 230Th, d = dissolved, p = particulate
       integer, parameter :: nb_part  = 4 ! (opal, caco2, poc, poussières)
       integer, parameter :: nopal = 1, ncaco3 = 2, npoc = 3, nlith=4 !! Attention dust = clay ??? 


       real(dp), dimension(:,:,:), pointer, public :: sea_mask ! 1 when grid-cell is ocean 0 when grid-cell is land

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim  Radioactve constants for 231Pa and 230Th
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! lim  Production Constants definition
!     Ref : Henderson & Anderson 2003

        real(dp), parameter :: beta_231Pa = 2.33E-3_dp   ! production rate of 231Pa in the water column (dpm/m3/y)
        real(dp), parameter :: beta_230Th = 2.52E-2_dp   ! production rate of 230Th in the water column (dpm/m3/y)

        real(dp), dimension(nb_trac) :: beta_trac_si ! production rate of (trac) in the water column (dpm/m3/s),
        ! SI = international system units (not opal!!!) we keep dpm to be consistent with observations !
        ! Note: this definition is due to code optimization,
        !       there is no physical sense to define beta_pap and beta_thp the production is made only in the water column and concerns the dissolved phase
        !       thus beta_trac_si(pap)=0 and beta_trac_si(thp)=0

! lim  Decay Constants definition
!     Ref : Audi et al,  2003

        real(dp), parameter :: lambda_231Pa = 2.116E-5_dp   ! Decay constant of 231Pa in y-1
        real(dp), parameter :: lambda_230Th = 9.195E-6_dp   ! Decay constant of 230 in y-1

        real(dp), dimension(nb_trac) :: lambda_trac_si ! Decay constant of (trac) in s-1

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim  Scavenging constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim  definition of the desorption constant for 231Pa and 230Th
        real(dp), parameter :: k_desorp = 2.4_dp ! (desorptions per year) Ref : Rempfer, subm // Marchal, 2000 : 3 desorptions per year
        real(dp) :: k_desorp_si ! k_desorp in s-1

! lim definition of the adsorptions constants for 231Pa and 230Th which are calculated in each grid-cell
        real(dp), dimension(imax_loc,jmax_loc,kmax_loc) :: k_adsorp_pa
        real(dp), dimension(imax_loc,jmax_loc,kmax_loc) :: k_adsorp_th

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim  definition of variables we need may be taken out from clio or somewhere else in the model
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim definition of the settling speed
        real(dp), parameter :: ws= 1000 ! sedimentation speed in m/y
        real(dp)  :: ws_si ! sedimentation speed in SI ie m/s

        
! lim  definition of the number of seconds per year in the model
        integer, parameter :: nb_s_y = 31104000 ! in s ( if in the model we have 360 days per year, 24 h per day, 3600 s per hour)
! lim  Note : to get from global constants ?? depends if we want to convert by real years or model years... 

! lim  definition of the time step
        integer, parameter :: time_step = 86400 ! time step of 1 day in seconds
        
! lim thickness of each water column layer in m,  to take from CLIO
        real(dp), dimension(:), pointer, public :: epais

! lim mask (table that contains 0 or 1) for each grid cell
! grid cells adjacent to continents should contain 1 / other grid cells should get 0         
        integer, dimension(imax_loc, jmax_loc, kmax_loc) :: mask_cont

! lim only used if fixed particles fields from NEMO-PISCES
! lim definition of the table that contains the fields [lon, lat, depth, the 4 particles concentrations (mol/m3) fields]
! lim formated for CLIO grid (repetition of first and last longitude
        real(dp), dimension(imax_loc,jmax_loc,kmax_loc, nb_part) :: particles_fields 

! lim definitions of the tables that will recieve the particles concentrations fields raw format from netcdf
! attention dimensions different because some longitudes are repeated for easier calculations
        real(dp), dimension(120,65,20) :: CaCO3_conc_field ! in mmol/m3
        real(dp), dimension(120,65,20) :: POC_conc_field   ! in mmol/m3
        real(dp), dimension(120,65,20) :: opal_conc_field  ! in mmol/m3

! restart flag for the Pa/Th 
        logical                                   :: restart_path = .true.
        character(len=str_len), parameter         :: restart_fnm  = "path_start_data.bin"        
        
        

! lim definition of the table that contains the fields [lon, lat, the 4 particles fluxes fields]
        real(dp), dimension(imax_loc,jmax_loc,kmax_loc, nb_part) :: particles_fluxes_field ! in mol/m2/s

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim  definition of the table that contains the activities
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! lim definition of the tables that contains the fields [lon, lat, depth, and the 4 tracers A231pa_d, A321pa_p, A230th_d, A230th_p]
        real(dp), dimension(imax_loc,jmax_loc,kmax_loc,nb_trac) :: Aactivity ! in dpm/m3

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim definition of the tend (what happens to 231Pa and 230Th during 1 time step
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        real(dp), dimension(imax_loc,jmax_loc,kmax_loc,nb_trac) :: trends_trac ! units dpm/m3/s (dissolved,particulate)*(Pa, Th)
! Note : trends_trac has to be set to 0 at each new time step performed ! 

      contains
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [read_particles_fields]
!
!>     @brief This subroutine reads particles concentration fields from netcdf file
!
!      DESCRIPTION:
!
!>     lim :  read the particles concentration CaCO3, opal and POC from netcfd file
!!     lim : the fields should be regridded on CLIO grid/ fields from NEMO-PISCES
!      lim be careful, the default file is NEMO-PISCES particles concentrations at modern conditions
!!     
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      subroutine read_particles_fields

      use ncio
        ! filin est le fichier netcdf à lire 
        ! dans le dossier /run

        CHARACTER*200 :: filin

        filin="inputdata/path/Part_conc_NEMO_onCLIO.nc"
        
        call nc_read(filin,"CaCO3", CaCO3_conc_field)
        call nc_read(filin,"POC", POC_conc_field)
        call nc_read(filin,"GSi", opal_conc_field)

      end subroutine read_particles_fields

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [transform_particles_fields]
!
!>     @brief This subroutine transforms the particles fields from netcdf file to usable tables
!
!      DESCRIPTION:
!
!>     lim : first step is to repeat the longitudes i.e. go from a 120 *65 *20 table to 122*65*20
!!     lim : second step to transform in mol/m3 
!      lim 3rd step to transform in mol/m2/s
!!     
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      subroutine transform_particles_fields
      
     ! lim debug loop needs i, j, k
      integer :: i, k, j

      particles_fields (2:121,:,:,nopal)= opal_conc_field(1:120,:,:)
      particles_fields (1,:,:,nopal)=opal_conc_field(120,:,:)
      particles_fields (122,:,:,nopal)=opal_conc_field(1,:,:)


      particles_fields (2:121,:,:,npoc)= POC_conc_field(1:120,:,:)
      particles_fields (1,:,:,npoc)=POC_conc_field(120,:,:)
      particles_fields (122,:,:,npoc)=POC_conc_field(1,:,:)

      particles_fields (2:121,:,:,ncaco3)= CaCO3_conc_field(1:120,:,:)
      particles_fields (1,:,:,ncaco3)=CaCO3_conc_field(120,:,:)
      particles_fields (122,:,:,ncaco3)=CaCO3_conc_field(1,:,:)

  ! lim : we put the particle concentration fields into the particles fluxes fields    
      particles_fluxes_field(:,:,:,nopal) = particles_fields(:,:,:,nopal)
      particles_fluxes_field(:,:,:,npoc) = particles_fields(:,:,:,npoc)
      particles_fluxes_field(:,:,:,ncaco3) = particles_fields(:,:,:,ncaco3)


!lim transformation from mmol/m3 to mol/m3
      particles_fluxes_field(:,:,:,:)= particles_fluxes_field(:,:,:,:)/1000

!lim tranformation from concentration (mol/m3) to flux (mol/m2/s)
!lim the transformation is made considering a constant and uniform sedimentation speed ws_si
      particles_fluxes_field(:,:,:,:)= particles_fluxes_field(:,:,:,:)*ws_si

! lim debug loop for setting to 0 the particles fluxes where mask problems are detected
! necessary otherwise the code wll crash
! TODO when the particles file is fine-erase those lines
      do i = 1,122
         do j = 1,65
            do k = 1,20
               if (particles_fluxes_field(i,j,k,ncaco3)> 1E29)then 
               particles_fluxes_field(i,j,k,ncaco3)= 0 
               particles_fluxes_field(i,j,k,npoc)= 0 
               particles_fluxes_field(i,j,k,nopal)= 0 
               endif
            enddo
         enddo
      enddo

      end subroutine transform_particles_fields


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [si_calc]
!
!>     @brief This subroutine caclulates the constants in SI units
!
!      DESCRIPTION:
!
!>     lim :  Basic caclulations from input data (conversions & cie)
!!     Is it ok to do those calculations here ??? Will we keep the right values for these variables ?
!!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine si_calc

        ! production rates in dpm/m3/s
        beta_trac_si(nthd) = beta_230Th/nb_s_y
        beta_trac_si(nthp) = 0.0_dp

        beta_trac_si(npad) = beta_231Pa/nb_s_y
        beta_trac_si(npap) = 0.0_dp

        ! Decay constants in s-1

        lambda_trac_si(nthd) = lambda_230Th/nb_s_y
        lambda_trac_si(nthp) = lambda_230Th/nb_s_y

        lambda_trac_si(npad) = lambda_231Pa/nb_s_y
        lambda_trac_si(npap) = lambda_231Pa/nb_s_y

        ! Desorption constant in s-1
        k_desorp_si = k_desorp/nb_s_y

        ! Settling speed in m/s
        ws_si = ws/nb_s_y

      end subroutine si_calc


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [k-adsorp_calc]
!
!>     @brief This subroutine caclulates the adsorption coefficients for 231Pa and 230Th in each grid cell
!
!      DESCRIPTION:
!
!>     Method Marchal et al 2000 and Rempfer, 2017
!!
!!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine  k_adsorp_calc

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! lim definition of the scavenging efficiencies (scalar)

        real(dp), parameter :: sigma_0 = 1.0_dp ! Reference scavenging efficiency in (m2/mol) Ref Rempfer, subm// Marchal, 2000 : 0.75 m2/mol

        real(dp), dimension(nb_part) :: sigma_part_th, sigma_part_pa ! in m2/mol

! lim definition of fractionation factors (scalar): f(Th/Pa)(note that there is no unit, it's a ratio)for the different particles types
!     if f = 1 no affinity difference for Pa and Th
!     if f > 1 stronger affinity of Th on the particle type
!     if f < 1 stronger affinity for Pa on the particle type

!     g is the relative affinity factor / describes the relative affinity for Pa for caco3 and opal
!     g > 1 stronger affinity for caco3
!     g < 1 stronger affinity for opal
       integer, parameter :: f_poc=1     ! Ref Rempfer, subm// Note : no units 
       integer, parameter :: f_caco3=10  ! Ref Rempfer, subm// Note : no units 
       integer, parameter :: f_opal=1    ! Ref Rempfer, subm// Note : no units 
       integer, parameter :: f_clay=10   ! Ref Rempfer, subm// Note : no units 
       integer, parameter :: g=1         ! Ref Rempfer, subm// Note : no units 

       real(dp),dimension(imax_loc,jmax_loc,kmax_loc,nb_part) :: vari_interim ! intermediate variable needed for calculations


       integer :: nt ! index for loop

       real(dp) :: scav_coeff           ! Scavenging coefficient for boundary scavenging parametrization

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       ! lim NOTE : is it necessary to set the last dimension of k_adsorp_th and k_adsrop_pa to 0 each time its calculated ??
       ! in order to avoid re-use of wrong figures...

       ! Step 1 : calculation of the sigma for Pa and Th for all particles types considered

#if ( PATH == 1 )
        ! A set of parameters which gives some acceptable Pa/Th results at the sea floor with particle fields from NEMO-PISCES
        ! The opal belt is not well represented, the carrier phase for Pa is opal, the coefficient for POC might not be the best fit
        ! quite bad agreement with the data in the upper part of the water column where POC has its highest concentrations

        ! for Pa
        sigma_part_pa(npoc)   = 1.5464172118167               ! in m2/mol
        sigma_part_pa(nopal)  = 7.62443575121686        ! in m2/mol
        sigma_part_pa(ncaco3) = 1.86698195773014            ! in m2/mol
        sigma_part_pa(nlith)  = sigma_0 * real(f_clay, kind=dp)             ! in m2/mol  !! DUST == CLAY ??????? [TODO]
        ! for Th
        sigma_part_th(npoc) = 5.46884090041493                                      ! in m2/mol
        sigma_part_th(nopal) = 3.76588888674175                 ! in m2/mol
        sigma_part_th(ncaco3) = 76.8301413063891                                    ! in m2/mol
        sigma_part_th(nlith) = sigma_0                                      ! in m2/mol

!endif PATH==1
#endif 

      ! For other cases with coupled particles, please check the tuning, the coefficient strongly rely on particle fields
      ! It might be a good idea to define the sigma parameters independently for Pa and Th (thus remove f_... and g)
#if ( PATH >= 2 )

        ! for Pa
        sigma_part_pa(npoc)   = sigma_0 * real(f_poc,kind=dp)               ! in m2/mol
        sigma_part_pa(nopal)  = sigma_0 * real(f_caco3 * g, kind=dp)        ! in m2/mol
        sigma_part_pa(ncaco3) = sigma_0 * real(f_caco3,kind=dp)             ! in m2/mol
        sigma_part_pa(nlith)  = sigma_0 * real(f_clay, kind=dp)             ! in m2/mol  !! DUST == CLAY ??????? [TODO]
        ! for Th
        sigma_part_th(npoc) = sigma_0                                       ! in m2/mol
        sigma_part_th(nopal) = sigma_0 * f_caco3 * g /f_opal                ! in m2/mol
        sigma_part_th(ncaco3) = sigma_0                                     ! in m2/mol
        sigma_part_th(nlith) = sigma_0                                      ! in m2/mol

!endif PATH>=2
#endif 

        ! Step 2 : calculation of the k_adsorp_pa and k_adsorp_th
        forall (nt=1:nb_part)
           vari_interim(:,:,:,nt) = particles_fluxes_field(:,:,:,nt)* sigma_part_pa(nt)
!             (s-1)                             (mol.m-2.s-1)       *     (m2.mol-1)
        endforall

        k_adsorp_pa(:,:,:)= SUM(vari_interim,DIM=4)
!              (s-1)               (s-1)

        forall (nt=1:nb_part)
           vari_interim(:,:,:,nt) = particles_fluxes_field(:,:,:,nt)* sigma_part_th(nt)
!               (s-1)                          ( mol.m-2.s-1)       *      (m2.mol-1)
        endforall

        k_adsorp_th(:,:,:)= SUM(vari_interim,DIM=4)
!           (s-1)                   (s-1)

! lim option : boundary scavenging
!              in all grid cells adjacent to continent and k_adsorp_pa multiplied by a coefficient

        scav_coeff = 0       ! 0 = no boundary scavenging // 1 value that imitates Rempfer 2017, 
                                                           ! all coeff multiplied by 2 at boundaries
        

        k_adsorp_pa(:,:,:) = k_adsorp_pa(:,:,:) + (k_adsorp_pa*mask_cont(:,:,:)*scav_coeff)
        
      end subroutine k_adsorp_calc

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine  k_adsorp_calc  here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [path_write_rest]
!
!>     @brief This function blaaah
!
!      DESCRIPTION:
!
!>    Dissolved activities initialized at the production ratio
!!    Particulate activities initialized to 0  !! To imporove : initialize at equilibrium ?
!!    How sophisticated the initialization should be ?
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine path_write_rest()
        
        use file_libs, only: fileDescriptor, open_f, close_f
        
! dmr --- File variable
        type(fileDescriptor)              :: restart_fich        
        
        ! dmr --- Formatted flag: .false. is unformatted file
        logical, parameter                :: restart_fm = .false.        


          ! dmr --- Set the formatted type and then open ...
        restart_fich%isFormatted = restart_fm
        call open_f(restart_fich, restart_fnm, o_stat="old")
        write(restart_fich%id) Aactivity
        call close_f(restart_fich)
        
      end subroutine path_write_rest

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [path_init]
!
!>     @brief This function initializes the dissolved and particulate activities in 231Pa and 230Th
!
!      DESCRIPTION:
!
!>    Dissolved activities initialized at the production ratio
!!    Particulate activities initialized to 0  !! To imporove : initialize at equilibrium ?
!!    How sophisticated the initialization should be ?
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine path_init(tracer_CLIO)


        use global_constants_mod, only: str_len
        use file_libs, only: fileDescriptor, open_f, close_f
 

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dp), dimension(:,:,:,:), intent(out) :: tracer_CLIO

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer :: nt


! dmr --- File variable
        type(fileDescriptor)              :: restart_fich        
        ! dmr --- File name variable (nm)
        character(len=str_len)            :: restart_nm   = "startdata/"//restart_fnm
        ! dmr --- Formatted flag: .false. is unformatted file
        logical, parameter                :: restart_fm = .false.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim continental mask for boundary scavenging TODO create a proper mask
          mask_cont(:,:,:) = 1._dp

! --- lim&dmr Call to si_calc necessary to calculate and initialize variables in SI like beta_trac_si
          call si_calc

       if (.not. restart_path) then

! lim Initialization at 0 like in NEMO 
          Aactivity(:,:,:,:) = 0._dp
          tracer_CLIO(:,:,:,:) = Aactivity(:,:,:,:)
! lim debug print
          write(*,*) "Je suis dans path_init / no restart !!", tracer_CLIO(1,1,1,:)

       else ! we read a restart file!!
          
          ! dmr --- Set the formatted type and then open ...
          restart_fich%isFormatted = restart_fm
          call open_f(restart_fich, restart_nm, o_stat="old")
          read(restart_fich%id) Aactivity
          call close_f(restart_fich)
          tracer_CLIO(:,:,:,:) = Aactivity(:,:,:,:)
          
       endif ! on restart_path true or false
       
#if ( PATH == 1 )
! if particles fixed
   ! lim reads the particles netcdf file default NEMO concentrations at modern
   ! Attention penser a traiter l'opale aussi en fixe tant que ca n existe pas dans iloveclim
       call read_particles_fields
   ! lim transforms concentrations in fluxes using constant sedimentation speed 

       call transform_particles_fields
!endif PATH==1
#endif 

#if ( PATH >= 2 )
! if particles interactive from ocycc
       particles_fluxes_field(:,:,:,:)=0._dp
!endif PATH==2
#endif 

!lim initialization of the clay for bottom scavenging
! clay is kept constant
! if no bottom scavenging set to : 0 
! if bottom scavenging is activated clay concentration from 40 to 1650 ug/L TODO convert in flux in mol/m2/s !!
! set to : 
       particles_fluxes_field(:,:,:,nlith)=0._dp


      end subroutine  path_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [prod_decay]
!
!>     @brief This subroutine is updating the trend taking into account the production and decay for particulate and dissolved 231Pa and 230Th
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine prod_decay

!       use AnotherModule_mod, only: some_variable
!       use AnotherModule_mod, only: some_otherfunction

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer :: nt

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        forall (nt=1:nb_trac)
           trends_trac(:,:,:,nt) = trends_trac(:, :,:,nt) + beta_trac_si(nt) - (lambda_trac_si(nt) * Aactivity(:,:,:,nt))
!            (dpm/m3/s)                  (dpm/m3/s)            (dpm/m3/s)            (s-1)              (dpm/m3)
        end forall

      end subroutine  prod_decay

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [scavenging]
!
!>     @brief This subroutine updates the trend taking into account the reversible scavenging
!
!      DESCRIPTION:
!
!>
!!
!!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine scavenging

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! lim trend is in dpm/m3/s
        ! lim NOTE : check units is it ok to multiply activities in dpm/s by adsorption coefficients in s-1 ?


        trends_trac(:,:,:,npad) = trends_trac(:,:,:,npad) + k_desorp_si * Aactivity(:,:,:,npap) &
                                  - k_adsorp_pa(:,:,:) * Aactivity(:,:,:,npad)
!              (dpm/m3/s)            (dpm/m3/s)      +         (s-1)     *     (dpm/m3)
!                                  -     (s-1)       *          (dpm/m3)        
        
        trends_trac(:,:,:,nthd) = trends_trac(:,:,:,nthd) + k_desorp_si * Aactivity(:,:,:,nthp) &
                                  - k_adsorp_th(:,:,:) * Aactivity(:,:,:,nthd)
!              (dpm/m3/s)            (dpm/m3/s)      +         (s-1)     *     (dpm/m3)
!                                  -     (s-1)       *          (dpm/m3)        

        trends_trac(:,:,:,npap) = trends_trac(:,:,:,npap) - k_desorp_si * Aactivity(:,:,:,npap) &
                                  + k_adsorp_pa(:,:,:) * Aactivity(:,:,:,npad)
!              (dpm/m3/s)            (dpm/m3/s)      +         (s-1)     *     (dpm/m3)
!                                  -     (s-1)       *          (dpm/m3)        

        trends_trac(:,:,:,nthp) = trends_trac(:,:,:,nthp) - k_desorp_si * Aactivity(:,:,:,nthp) &
                                  + k_adsorp_th(:,:,:) * Aactivity(:,:,:,nthd)
!              (dpm/m3/s)            (dpm/m3/s)      +         (s-1)     *     (dpm/m3)
!                                  -     (s-1)       *          (dpm/m3)        

      end subroutine  scavenging

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [settling]
!
!>     @brief This subroutine updates the trend taking into account the settling of the particles
!
      !      DESCRIPTION:
      ! The settling has 2 parts : what happens in the water column and in the interface with sediments (bottom layer)
      ! in the water column each grid cell receives dissolved activity from the top and loses some to the bottom
      ! at the interface to the sediments, what goes out of the water cell to the bottom is lost
      ! note :  do we want to pass it to MEDUSA later ?
      ! note : in this version,  there is only 1 particle type and so 1 settling speed
!
!>
!!
!!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine settling

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        real(dp), dimension(imax_loc, jmax_loc, kmax_loc,nb_trac) :: Zflux ! Activity flux for each Z level of the ocean
        ! (particulate Pa npap =2, particulate Th nthp =4)
        ! Note : we actually need only 2 dimensions as we do not settle dissolved activities...
        ! Note : Zflux is in dpm/m2/s
        ! Index problems : fluxes calculated at the borders, activities defined in the center of grid cells ??
        ! Note : local variable

        integer :: i,j,k,n ! local indexes to make loops inside this subroutine
        integer :: inc ! increment for do loop
        
   

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! initialization Zflux to zero as we need the surface row to be zero (there is no flux down from the surface)
        Zflux(:,:,:,:)=0._dp ! dpm/m2/s

        inc=-1
        
! lim : first we calculate the vertical flux for each Z level of the ocean grid        
        DO n=1,nb_trac
           
           DO i=1,imax_loc
           
              DO j=1,jmax_loc
              
                 DO k=k_surf(i,j),k_fond(i,j),inc ! Attention !!! index ok ??
                    
                    ! --- Zflux contient le flux sortant de la case n vers la case n-1
                    !          (indices inversés CLIO, i.e. du haut vers le bas)
                    !     Zflux est défini sur [|k_surf:k_fond|]
                    !     NOTA : Zflux(...,k_surf(i,j)) n'est pas défini
                    Zflux(i,j,k,n) = ws_si *Aactivity(i,j,k,n)  !!! Attention formulae
!                   (dpm/m2/s)       (m/s)  *    (dpm/m3)

                    
                 ENDDO
              
              ENDDO
           
           ENDDO
          
        ENDDO

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! lim then we update the trends with the balance between what comes in or out via sedimentation
! Note : the trend is in Bq/m3/s        

        DO i=1,imax_loc
           
              DO j=1,jmax_loc


! --- lim&dmr Cas spécifique de la surface, uniquement un flux sortant de particules

                 k = k_surf(i,j)
                 trends_trac(i,j,k,npap)=trends_trac(i,j,k,npap) + ((0.0_dp-Zflux(i,j,k,npap))/epais(k))
!                      (dpm/m3/s)               (dpm/m3/s)       +    (dpm/m2/s)                          /  (m)                  
                    
                  trends_trac(i,j,k,nthp)=trends_trac(i,j,k,nthp) + ((0.0_dp-Zflux(i,j,k,nthp))/epais(k))
!                      (dpm/m3/s)               (dpm/m3/s)        +    (dpm/m2/s)                          /  (m)       


! --- lim&dmr cas général pour toutes les cases sauf surface puisque
!             flux entrant est nul
                 DO k=k_surf(i,j)+inc,k_fond(i,j),inc
                    
                    trends_trac(i,j,k,npap)=trends_trac(i,j,k,npap) + ((Zflux(i,j,k-inc,npap)-Zflux(i,j,k,npap))/epais(k))
!                         (dpm/m3/s)               (dpm/m3/s)       +    (dpm/m2/s)                          /  (m)                  
                    
                    trends_trac(i,j,k,nthp)=trends_trac(i,j,k,nthp) + ((Zflux(i,j,k-inc,nthp)-Zflux(i,j,k,nthp))/epais(k))
!                         (dpm/m3/s)               (dpm/m3/s)       +    (dpm/m2/s)                          /  (m)                  

                 ENDDO
              
              ENDDO
           
           ENDDO
          

 ! lim  Burial into the sediments ?? already done in these loops or not ???? I think it's already done but I am not sure
 ! flux_to_burial should be the last index of the local ocean, i.e. Zflux(i,j,k_fond(i,j),:)
        
        return
        
       
      end subroutine  settling

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
     

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [path_step]
!
!>     @brief This subroutine updates the activities A231Pa_p, A231Pa_d, A230Th_p, A230Th_d after CLIO's transport
!
!      DESCRIPTION:
!
!>
!!
!!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine path_step

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer :: nt ! to make a loop on Aactivity
     ! lim debug loop needs i, j, k
      integer :: i, k, j,c,d

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       ! Set the trends to 0 at the beginning of each time steps

        trends_trac(:,:,:,:) = 0.0_dp ! (dpm/m3/s)

        ! Call all the subroutines that modify the trend :

        ! Production and decay
        call prod_decay

        ! Scavenging = adsorption & desorption
        ! lim : calcualtion of the adsorption coefficients
        call k_adsorp_calc
        ! lim k_adsorp caclulation not tested !!! 
        ! lim : updating the trends
        call scavenging

        ! Settling = sedimentation of particles
        call settling

        ! Updating the activities with trends
        Aactivity(:,:,:,:) = Aactivity(:,:,:,:) + trends_trac(:,:,:,:) * time_step
 ! lim debug print Aactivity après update
     
        ! Multiplying by the sea mask 

        forall (nt = 1:nb_trac)
          Aactivity(:,:,:,nt) = Aactivity(:,:,:,nt)*sea_mask(:,:,:)
       endforall 




      end subroutine  path_step

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module path_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
