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
!      MODULE: neodymium_mod
!
!>     @author  Tristan Vadsaria (tva)
!
!>     @brief This module [Nd_mod] is computing Nd143 and Nd144, treated as tracers :
!                                                                      Surface flux from dust and rivers
!                                                                      Sediment source as global flux from seafloor-ocean boundary
!                                                                      Reversible scavenging as adsorption and desorption ont
!                                                                      Inserting 144Nd and 143Nd as prognostic tracers
!
!>     @date Creation date: October, 19th, 2022
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : tva
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!     Based on what has been developed for Bern3D model, NEMO-PISCES, CESM, FAMOUS
!!     Useful references : Arsouze et al. (2009); Rempfer et al. (2011); Gu et al. (2019); Poppelmeier et al. (2020, 2022)
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

    module neodymium_mod
       use para0_mod, only:  imax_loc=>imax, jmax_loc=>jmax, kmax_loc=>kmax
       use bloc0_mod, only: ks2, dz, tms, kfs
       use bloc_mod, only: aire, kfond
       use dynami_mod, only: dxc1, dxc2
       use global_constants_mod, only: dp, str_len
       use bloc0_mod, only: fpoc => fPOC_flx_clio !TmolsC.m-2.timestep-1
       use bloc0_mod, only: fcaco3 => fCAL_flx_clio !TmolsC.m-2.timestep-1
       use ncio
       

       implicit none

       private :: dust_source_calculation
       private :: river_source_calculation
       private :: benthic_flux_calculation
       private :: reversible_scavenging
       private :: downward_flux
       public  :: neodymium_init
       public  :: neodymium_step
       public  :: init_source_netcdf
       public  :: coeff_init ! private? maybe merge with neodymium_init, check with dmr  
       
       


       ! NOTE_avoid_public_variables_if_possible

! --- tva Following variables define the position of nd variables wihin the CLIO scal array
       integer, parameter, public :: scalstart_neodymium = 18, scalend_neodymium = 19

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tva  Tracers and index
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! tva Number of tracer: due to the Nd scheme employed, we only need to simulate explicitly the total concentration of the two isotopes nd143 and nd144 rather than concentration in every phase (dissolved, particulate)
       integer, parameter :: nb_trac  = 2 

! tva indices of the tracers
       integer, parameter :: nd143 = 1, nd144 = 2

! tva in case of prescribed particles
       integer, parameter :: nb_part  = 4 
       integer, parameter :: nopal = 1, ncaco3 = 2, npoc = 3

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tva Definition of the tables that contains the advected [lon, lat, depth, and the 2 tracers: 144Nd and 143Nd (neodymium)
! tva ... for the total concentration of Nd as well as for the particulate and the dissolve phase (neod_diss, neod_part)
! tva ... and for the concentration due to the vertical cycling (z_flux, z_neod)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       real(dp), dimension(imax_loc,jmax_loc,kmax_loc,nb_trac) :: neodymium, neod_diss, neod_part, z_neod, z_flux
!DIAG
       real(dp) :: Nd_inventory, Nddiss_inventory, Ndpart_inventory, z_neod_SUM, z_flux_source, z_neod_SUM2

       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! tva Definitations of CHUR
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       real(dp), parameter :: CHUR = 0.512638
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! tva Definitations for the dust source component
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! tva Definition of the input tables for reading the isotopic ratio (eNd) of dust and the dust flux in kg/m2/s
       real(dp), dimension(120,65) :: ir_dust, dust_flux
! tva Definition of the table that will receive the isotopic ratio (eNd) of dust and the dust flux after longitude transformation (halo point)
       real(dp), dimension(imax_loc,jmax_loc) :: ir_dust_transform, dust_flux_transform, dust_flux_transform_g_m2_s
! tva Definition of the value of the Nd concentration of dust by Goldstein et al. (1984); Grousset et al. (1988, 1998): 20 microgram per gram
       real(dp), parameter :: conc_dust = 0.00002
! tva Definition of the value of the dissolution of dust in seawater, 2 percent, Greaves et al. (1994)
       real(dp), parameter :: beta_dust = 0.02
! tva Definition of the Nd dust field tables
       real(dp), dimension(imax_loc,jmax_loc) :: dust_source_Nd, dust_source_144Nd, dust_source_143Nd 
!DIAG
       real(dp) :: dust_source_Nd_tot

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! tva Definitations for the river source component
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! tva Definition of the input tables for reading the isotopic ratio (eNd) and the Nd concentration of rivers in ppt
       real(dp), dimension(120,65) :: ir_river, conc_river, ir_river_tmp, conc_river_tmp
! tva Same as above after longitude transformation (halo point) and unit conversion
       real(dp), dimension(imax_loc,jmax_loc) :: ir_river_transform, conc_river_transform, conc_river_transform_g
! tva Converting the river runoff from -kg/m2/day to g/m2/s 
       real(dp), dimension(imax_loc,jmax_loc) :: fwruno_g_m2_s  
! tva Dissolution of Nd in river esturaries (REF)       
       real(dp), parameter :: river_remob = 0.3
! tva Definition of the Nd river field tables
       real(dp), dimension(imax_loc,jmax_loc) :: river_source_Nd, river_source_144Nd, river_source_143Nd
!DIAG
       real(dp) :: river_source_Nd_tot, area_tot, area2_tot, fwruno_tot, fwruno_area_tot
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! tva Definitations for the benthic source component
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! tva Definition of the input tables for reading the isotopic ratio (eNd) of the sediment
       real(dp), dimension(120,65) :: ir_benthic
! tva Same as above after longitude transformation (halo point) and unit conversion
       real(dp), dimension(imax_loc,jmax_loc) :: ir_benthic_transform
! tva Definition of the Nd sediment global flux in g/yr (TMP for the implementation: Sediment flux from 3.0*10+9 g/yr (best sim Robinson et al. 2022) to g/s)
! 2be put later in a namelist or whatever for tuning and/or space phase parameterization!!!!!!!!!
       !real(dp), parameter :: benthic_source_Nd_global = 3000000000/(86400*365.25) 
       real(dp), parameter :: benthic_source_Nd_global = 95.06426344208685 ! 3000000000/(86400*365.25) = g/s
! tva Definition of the Nd sediment field tables
       real(dp), dimension(imax_loc,jmax_loc) :: benthic_source_Nd, benthic_source_144Nd, benthic_source_143Nd
!DIAG
       real(dp) :: benthic_source_Nd_tot, seabed_surf_check

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! tva Definitations for reversible scavenging component
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! tva Definition of the equilibrium partition coefficient, first coming from Arsouze et al. 2009 for the implementation part
! 2be put later in a namelist or whatever for tuning and/or space phase parameterization!!!!!!!!!
       real(dp), dimension(imax_loc,jmax_loc,kmax_loc,nb_trac) :: kpoc, kcaco3, kdsi, denum, mult
       real(dp), dimension(imax_loc,jmax_loc,kmax_loc) :: fpoc_g_m2_day, poc_g_m3, poc
       real(dp), dimension(imax_loc,jmax_loc,kmax_loc) :: fcaco3_g_m2_day, caco3_g_m3, caco3  
       real(dp), dimension(imax_loc,jmax_loc,kmax_loc) :: opal      
       real(dp) :: kpoc_143Nd, kpoc_144Nd ,kcaco3_143Nd, kcaco3_144Nd, Zwater, Zcaco3, Zpoc, kdsi_143Nd, kdsi_144Nd, Zdsi
       ! in case of prescribed particles
       real(dp), dimension(imax_loc,jmax_loc,kmax_loc, nb_part) :: particles_prescribed, particles_prescribed_g_g 
       real(dp), dimension(120,65,20) :: CaCO3_conc_field ! in mmol/m3
       real(dp), dimension(120,65,20) :: POC_conc_field   ! in mmol/m3
       real(dp), dimension(120,65,20) :: opal_conc_field  ! in mmol/m3

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! tva Definitations for the downward flux component
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! tva Definition of the vertical sinking of particle in m/s from 1000 m/yr
       real(dp), parameter :: ws=0.0000321502 !1000/(86400*360)
       
      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! INIT PART
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! INIT OF THE NEODYMIUM SCALAR
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine neodymium_init(tracer_CLIO)
        real(dp), dimension(:,:,:,:), intent(out) :: tracer_CLIO
        neodymium(:,:,:,:) = 0._dp
        tracer_CLIO(:,:,:,:) = neodymium(:,:,:,:)
        !write(*,*) "Nd init fond", neodymium(:,:,kfond,nd144)
      end subroutine neodymium_init
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! INIT OF NEODYMIUM SOURCE (NECTDF READING)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine init_source_netcdf
      integer :: i,j
      CHARACTER(len=200) :: filin1, filin2
      CHARACTER(len=200) :: filin3, filin4
      CHARACTER(len=200) :: filin6
      CHARACTER(len=200) :: filin7
      
!DUST

       filin1="inputdata/neodymium/eNd_dust_map_regridded_on_CLIO_remapnn.nc"
       call nc_read(filin1,"eNd_dust", ir_dust)
       
       filin2="inputdata/neodymium/LMDZ_dustdep_1960-2021_regridded_on_CLIO_remapnn.nc"
       call nc_read(filin2,"dustdeptotal_D1cisi", dust_flux)       

       ir_dust_transform (2:121,:)= ir_dust(1:120,:)
       ir_dust_transform (1,:)=ir_dust(120,:)
       ir_dust_transform (122,:)=ir_dust(1,:)

       dust_flux_transform (2:121,:)= dust_flux(1:120,:)
       dust_flux_transform (1,:)=dust_flux(120,:)
       dust_flux_transform (122,:)=dust_flux(1,:)
       
!RIVERS
      
       filin3="inputdata/neodymium/River_eNd_regridded_on_CLIO_remapnn.nc"
       call nc_read(filin3,"river_end", ir_river)

       filin4="inputdata/neodymium/River_Nd_ppt_regridded_on_CLIO_remapnn.nc"
       call nc_read(filin4,"river_nd_ppt", conc_river)
       
       !filin5="inputdata/neodymium/CLIO3_NewGen_fwruno.nc"
       !call nc_read(filin5,"fwruno", fwruno)       

       do i=1,120
          do j=1,65
            if (ir_river(i,j) < -1000) then
              ir_river_tmp(i,j) = 0
              conc_river_tmp(i,j) = 0
            else 
              ir_river_tmp(i,j) = ir_river(i,j)
              conc_river_tmp(i,j) = conc_river(i,j)
            end if
!            write(*,*) conc_river(i,j), conc_river_tmp(i,j)
          end do
       end do
       
       ir_river_transform (2:121,:)= ir_river_tmp(1:120,:)
       ir_river_transform (1,:)=ir_river_tmp(120,:)
       ir_river_transform (122,:)=ir_river_tmp(1,:)

       conc_river_transform (2:121,:)= conc_river_tmp(1:120,:)
       conc_river_transform (1,:)=conc_river_tmp(120,:)
       conc_river_transform (122,:)=conc_river_tmp(1,:)  
       
!BENTHIC FLUX
      
       filin6="inputdata/neodymium/globalSedimentEpsilon_regridded_on_CLIO_remapnn.nc"
       call nc_read(filin6,"globalsedimentepsilonnd", ir_benthic)

       ir_benthic_transform(2:121,:)= ir_benthic(1:120,:)
       ir_benthic_transform(1,:)=ir_benthic(120,:)
       ir_benthic_transform(122,:)=ir_benthic(1,:)     

!READ PARTICLE CONC FIELD IN CASE OF PARTICLE PRESCRIPTION

       filin7="inputdata/path/Part_conc_NEMO_onCLIO.nc"
        
       call nc_read(filin7,"CaCO3", CaCO3_conc_field)
       call nc_read(filin7,"POC", POC_conc_field)
       call nc_read(filin7,"GSi", opal_conc_field)

       particles_prescribed(2:121,:,:,nopal)= opal_conc_field(1:120,:,:)
       particles_prescribed(1,:,:,nopal)=opal_conc_field(120,:,:)
       particles_prescribed(122,:,:,nopal)=opal_conc_field(1,:,:)
 
       particles_prescribed(2:121,:,:,npoc)= POC_conc_field(1:120,:,:)
       particles_prescribed(1,:,:,npoc)=POC_conc_field(120,:,:)
       particles_prescribed(122,:,:,npoc)=POC_conc_field(1,:,:)
 
       particles_prescribed(2:121,:,:,ncaco3)= CaCO3_conc_field(1:120,:,:)
       particles_prescribed(1,:,:,ncaco3)=CaCO3_conc_field(120,:,:)
       particles_prescribed(122,:,:,ncaco3)=CaCO3_conc_field(1,:,:)
       
       ! from mmol/m3 to g/g
       particles_prescribed_g_g(:,:,:,nopal) = ((((particles_prescribed(:,:,:,nopal))/1000)*28)/1025000)*tms
       particles_prescribed_g_g(:,:,:,npoc) = ((((particles_prescribed(:,:,:,npoc))/1000)*12)/1025000)*tms
       particles_prescribed_g_g(:,:,:,ncaco3) = ((((particles_prescribed(:,:,:,ncaco3))/1000)*100.09)/1025000)*tms
       
      end subroutine init_source_netcdf
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! INITIALIZATIONS OF TRACER PARTITION COEFFICIENTS
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      

      subroutine coeff_init
!          conversion factor partition coefficient: micro_mol/l-->g/l (water density)
!          ----------------------------------------------------------------------------
      integer :: i,j,k
      Zwater = 1000
      Zcaco3 = 40. + 12. + 3.*16.
      Zpoc   = 32.7
      !Zgoc  = 32.7
      Zdsi  = 28. + 2.* 16.
      !Zlitho = Zcaco3

      ! The following k values need 2be put later in a namelist or whatever for tuning and/or space phase parameterization!!!!!!!!!      
      kpoc_143Nd = 1.4e+7
      kpoc_144Nd = 1.4e+7
      !kgoc_143Nd = 5.2e+4
      !kgoc_144Nd = 5.2e+4
      kdsi_143Nd = 3.6e+4
      kdsi_144Nd = 3.6e+4
      kcaco3_143Nd = 1.6e+5
      kcaco3_144Nd = 1.6e+5
      !klitho_143Nd = 4.6e+5
      !klitho_144Nd = 4.6e+5

      
       do k=1,kmax_loc
        do j=1,jmax_loc
           do i=1,imax_loc
            kpoc(i,j,k,nd143) = 457800! kpoc_143Nd * Zpoc / Zwater
            kpoc(i,j,k,nd144) = 457800! kpoc_144Nd * Zpoc / Zwater
            !kgoc (ji,jj,jk,1) = kgoc_143 * Zgoc / Zwater
            !kgoc (ji,jj,jk,2) = kgoc_144 * Zgoc / Zwater            
            kdsi (i,j,k,nd143) = 2160 !kdsi_143Nd * Zdsi / Zwater
            kdsi (i,j,k,nd144) = 2160 !kdsi_144Nd * Zdsi / Zwater
            kcaco3(i,j,k,nd143) = 16000!kcaco3_143Nd * Zcaco3 / Zwater
            kcaco3(i,j,k,nd144) = 16000!kcaco3_144Nd * Zcaco3 / Zwater
            !klitho (ji,jj,jk,1) = klitho_143 * Zlitho / Zwater
            !klitho (ji,jj,jk,2) = klitho_144 * Zlitho / Zwater
           end do
         end do
       end do

!       write (*,*) "!!!!!!!!!!!!!!!"
!       write (*,*) "kpoc", kpoc(0,0,0,0)
!       write (*,*) "!!!!!!!!!!!!!!!"
       
!       write (*,*) "!!!!!!!!!!!!!!!"
!       write (*,*) "kcaco3", kcaco3(0,0,0,0)
!       write (*,*) "!!!!!!!!!!!!!!!"       
                        
      end subroutine coeff_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! CALCULATION OF EPSND BASED ON [ND143] and [ND144] (used in grid_io_nc.f08)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      pure elemental function epsilonNd(Nd143, Nd144) result(epsNd)
       real, intent(in) ::Nd143, Nd144
       real ::  epsNd
       epsNd = (((Nd143/Nd144)/CHUR)-1)*10000
      end function epsilonNd

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! SOURCES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! DUST SOURCE CALCULATION
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine dust_source_calculation
#if ( NEOD == 1 )
      integer :: i,j
       
       dust_flux_transform_g_m2_s = dust_flux_transform*1000 !dust flux from kg/m2/s to g/m2/s
       
       write(*,*) "dust source calculation"
       dust_source_Nd = dust_flux_transform_g_m2_s*conc_dust*beta_dust
       dust_source_144Nd = 0.238*dust_source_Nd !the abundance of 144Nd is 23.8 per cent of the total Nd
       dust_source_143Nd = ((ir_dust_transform/10000)+1)*dust_source_144Nd*CHUR ! Using the eNd definition to calculate 143Nd from 144Nd (Arsouze et al. 2009)


!DIAG
       dust_source_Nd_tot =0.
       do i=1,imax_loc
          do j=1,jmax_loc
            dust_source_Nd_tot=dust_source_Nd_tot+dust_source_Nd(i,j)*dxc1(i,j)*dxc2(i,j)*tms(i,j,ks2)*86400*360
          end do
       end do
       write (*,*) "!!!!!!!!!!!!!!!"
       write (*,*) "total Nd dust flux", dust_source_Nd_tot*0.000000001, "(10^9)g/yr"
       write (*,*) "!!!!!!!!!!!!!!!"
#else
       dust_source_144Nd = 0._dp
       dust_source_143Nd = 0._dp
#endif
      end subroutine dust_source_calculation

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! RIVER SOURCE CALCULATION
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine river_source_calculation      
#if ( NEOD == 1 )
      use ice_mod, only: fwruno
      integer :: i,j  
       
       fwruno_g_m2_s = (-1*fwruno*1000)/86400 ! river runoff from -kg/m2/day to g/m2/s
       conc_river_transform_g = conc_river_transform*1.0E-9*1.0E-3 ! converting Nd ppt into g/kg (1.0E-9)... and convert g to kg (1.0E-3)

       
       write(*,*) "river source calculation"
       !write(*,*) "check river source calculation" 
       river_source_Nd =  conc_river_transform_g  * river_remob * fwruno_g_m2_s

       river_source_144Nd = 0.238*river_source_Nd !the abundance of 144Nd is 23.8 per cent of the total Nd
       river_source_143Nd = ((ir_river_transform/10000)+1)*river_source_144Nd*CHUR ! Using the eNd definition to calculate 143Nd from 144Nd (Arsouze et al. 2009) 
 
!       write(*,*) "river source 144Nd calculation", MINVAL(MINVAL(river_source_144Nd,dim=1),dim=1), MAXVAL(MAXVAL(river_source_144Nd,dim=1),dim=1)
       !write(*,*) SUM(SUM(river_source_144Nd,dim=1),dim=1)

!       write(*,*) "river source 143Nd calculation", MINVAL(MINVAL(river_source_143Nd,dim=1),dim=1), MAXVAL(MAXVAL(river_source_143Nd,dim=1),dim=1)
       !write(*,*) SUM(SUM(river_source_143Nd,dim=1),dim=1)       

!DIAG
       river_source_Nd_tot =0.
       area_tot=0.
       area2_tot=0.
       fwruno_tot =0.
       fwruno_area_tot =0.
       do i=1,imax_loc
          do j=1,jmax_loc
            area_tot = area_tot + dxc1(i,j)*dxc2(i,j) ! in m2
            area2_tot = area2_tot + aire(i,j) ! in m2
            fwruno_tot = fwruno_tot + fwruno(i,j) ! in -kg/m2/day
            !fwruno_area_tot = fwruno_area_tot + (fwruno(i,j)*dxc1(i,j)*dxc2(i,j)) ! in -kg/day
            fwruno_area_tot = fwruno_area_tot +(fwruno_g_m2_s(i,j) *dxc1(i,j)*dxc2(i,j))! in g/s
            river_source_Nd_tot=river_source_Nd_tot+river_source_Nd(i,j)*dxc1(i,j)*dxc2(i,j)*tms(i,j,ks2)*86400*360
          end do
       end do
       
       write (*,*) "!!!!!!!!!!!!!!!"
       !write (*,*) "total freshwater to the ocean", fwruno_area_tot*(-1)*360*0.001*0.000000001, "(km3/yr)"
       write (*,*) "total freshwater to the ocean", fwruno_area_tot*360*86400*0.000001*0.000000001, "(km3/yr)"
       write (*,*) "!!!!!!!!!!!!!!!" 

       write (*,*) "!!!!!!!!!!!!!!!"
       write (*,*) "total Nd river flux", river_source_Nd_tot*0.000000001, "(10^9)g/yr"
       write (*,*) "!!!!!!!!!!!!!!!"
       
!       write (*,*) "!!!!!!!!!!!!!!!"
!       write (*,*) "surface of the ocean", area_tot, area2_tot
!       write (*,*) "!!!!!!!!!!!!!!!"      
       
!       write (*,*) "!!!!!!!!!!!!!!!"
!       write (*,*) "total freshwater to the ocean", fwruno_tot, ("-kg/m2/day")
!       write (*,*) "!!!!!!!!!!!!!!!"         
       
!       write (*,*) "!!!!!!!!!!!!!!!"
!       write (*,*) "total freshwater to the ocean", fwruno_area_tot, ("-kg/day")
!       write (*,*) "!!!!!!!!!!!!!!!" 


#else
       river_source_144Nd = 0._dp
       river_source_143Nd = 0._dp
#endif
      end subroutine river_source_calculation

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! BENTHIC SOURCE CALCULATION
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine benthic_flux_calculation
#if ( NEOD == 1 )
      integer :: i,j    
       
       write(*,*) "benthic source calculation"
       !benthic_source_Nd = benthic_source_Nd_global/364012054606928! Sediment flux in g/m2/s by dividing by the surface of seabed (364012054606928. m2 in CLIO)
       benthic_source_Nd = 2.611569101598581E-13 !95.06426344208685/361383969000000
       benthic_source_144Nd = 0.238*benthic_source_Nd !the abundance of 144Nd is 23.8 per cent of the total Nd
       benthic_source_143Nd = ((ir_benthic_transform/10000)+1)*benthic_source_144Nd*CHUR ! Using the eNd definition to calculate 143Nd from 144Nd (Arsouze et al. 2009)
       
       !write(*,*) "benthic source Nd global", benthic_source_Nd_global
       !write(*,*) "benthic source Nd", MINVAL(MINVAL(benthic_source_Nd,dim=1),dim=1), MAXVAL(MAXVAL(benthic_source_Nd,dim=1),dim=1)
       !write(*,*) "benthic source 144Nd", MINVAL(MINVAL(benthic_source_144Nd,dim=1),dim=1), MAXVAL(MAXVAL(benthic_source_144Nd,dim=1),dim=1)       
       !write(*,*) "benthic source 143Nd", MINVAL(MINVAL(benthic_source_143Nd,dim=1),dim=1), MAXVAL(MAXVAL(benthic_source_143Nd,dim=1),dim=1)


!DIAG
       benthic_source_Nd_tot =0.
       seabed_surf_check =0.
       do i=1,imax_loc
          do j=1,jmax_loc
            !write(*,*) "kfs(i,j)", kfs(i,j)
            benthic_source_Nd_tot=benthic_source_Nd_tot+benthic_source_Nd(i,j)*dxc1(i,j)*dxc2(i,j)*tms(i,j,kfs(i,j))*86400*360
            !seabed_surf_check = seabed_surf_check + tms(i,j,kfs(i,j))*dxc1(i,j)*dxc2(i,j)
          end do
       end do

       !write (*,*) "!!!!!!!!!!!!!!!"
       !write (*,*) "seabed surf check", seabed_surf_check
       !write (*,*) "!!!!!!!!!!!!!!!"

       write (*,*) "!!!!!!!!!!!!!!!"
       write (*,*) "total Nd benthic flux", benthic_source_Nd_tot*0.000000001, "(10^9)g/yr"
       write (*,*) "!!!!!!!!!!!!!!!"


#else
       benthic_source_144Nd = 0._dp
       benthic_source_143Nd = 0._dp
#endif
      end subroutine benthic_flux_calculation

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! REVERSIBLE SCAVENGING
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! REVERSIBLE SCAVENGING CALCULATION
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine reversible_scavenging
#if ( NEOD == 1 )
      integer :: i,j,k,n
      real(dp) :: poc_SUM, caco3_SUM
      fpoc_g_m2_day = fpoc*12*1E12 ! TmolsC to g: (1 x 10^12 moles/Tmol) * (12.01 grams/mole) --> flux from TmolsC/m2/day to g/m2/s
      poc_g_m3 = fpoc_g_m2_day/1000  ! To get a concentration in g/m3: division by the settling speed (ws in m/s in ws*360 m/day or 1000/360) 
      !poc = poc_g_m3 /1025000 !  division by the density of seawater to have g/kg

      fcaco3_g_m2_day = fcaco3*100.09*1E12 ! TmolsC to g: (1 x 10^12 moles/Tmol) * (100.09 grams/mole) --> flux from TmolsC/m2/yr to g/m2/yr
      caco3_g_m3 = fcaco3_g_m2_day/1000  ! To get a concentration in g/m3: division by the settling speed (ws in m/s in ws*360 m/day or 1000/360)
      !caco3 = caco3_g_m3 /1025000     !  division by the density of seawater to have g/kg
      
      poc = particles_prescribed_g_g(:,:,:,npoc)
      caco3 = particles_prescribed_g_g(:,:,:,ncaco3)
      opal = particles_prescribed_g_g(:,:,:,nopal)

       !do i=1,imax_loc       
       !   do j=1,jmax_loc
       !      do k=1,kmax_loc   
       !       write(*,*) "fpoc", fpoc(i,j,k)
       !       write(*,*) "fpoc_g_m2_day", fpoc_g_m2_day(i,j,k)
       !       write(*,*) "poc_g_m3", poc_g_m3(i,j,k)
       !       write(*,*) "poc", poc(i,j,k)      
       !       write(*,*) "fcaco3", fcaco3(i,j,k)              
       !       write(*,*) "fcaco3_g_m2_day", fcaco3_g_m2_day(i,j,k)
       !       write(*,*) "caco3_g_m3", caco3_g_m3(i,j,k)
       !       write(*,*) "caco3", caco3(i,j,k)                          
       !       end do
       !   end do
       !end do


      
!  Calculation of [Nd] in the dissolve phase for each isotope and each particle type
       do n=1,nb_trac
          do i=1,imax_loc       
             do j=1,jmax_loc
                 do k=1,kmax_loc   
                   if (neodymium(i,j,k,n)<0) then
                     neodymium(i,j,k,n) = 0._dp
                   endif
                   !denum(i,j,k,n) = (1. +  poc(i,j,k) * kpoc(i,j,k,n) + caco3(i,j,k)  * kcaco3(i,j,k,n))
                   denum(i,j,k,n) = (1. +  poc(i,j,k) * kpoc(i,j,k,n) + caco3(i,j,k)  * kcaco3(i,j,k,n) +  opal(i,j,k)  * kdsi(i,j,k,n))
                   neod_diss(i,j,k,n) = neodymium(i,j,k,n) / denum(i,j,k,n)

                   !if (neod_diss(i,j,k,n)<0) then
                     !neod_diss(i,j,k,n) = 0._dp
                    ! write(*,*) neodymium(i,j,k,n), denum(i,j,k,n)
                   !endif

                   !write(*,*) "1 + K*C",i,j,k, (1. +  poc(i,j,k) * kpoc(i,j,k,n) + caco3(i,j,k)  * kcaco3(i,j,k,n))
                   !write(*,*) "caco3",i,j,k, caco3(i,j,k)
                   !write(*,*) "poc",i,j,k, poc(i,j,k)
                   !write(*,*) "kcaco3*caco3", kcaco3(i,j,k,n)*caco3(i,j,k) 
                   !write(*,*) "kpoc*poc", kpoc(i,j,k,n)*poc(i,j,k) 
                   !write(*,*) "...",i,j,k, (poc(i,j,k) * kpoc(i,j,k,n) + caco3(i,j,k)  * kcaco3(i,j,k,n))
                   !write(*,*) "denum",i,j,k, denum(i,j,k,n)
                   !write(*,*) "NDD, NDT, NDD/NDT",i,j,k, neod_diss(i,j,k,n), neodymium(i,j,k,n), neod_diss(i,j,k,n)/neodymium(i,j,k,n)
                   !write(*,*) "NDD/NDT",i,j,k, neod_diss(i,j,k,n)/neodymium(i,j,k,n)
                 end do
             end do
          end do
       end do
       
!  Calculation of [Nd] in the particulate phase for each isotope and each particle type
       do n=1,nb_trac
          do i=1,imax_loc       
             do j=1,jmax_loc
                 do k=1,kmax_loc
                   if (neodymium(i,j,k,n)<0) then
                     neodymium(i,j,k,n) = 0._dp
                   endif
                   !mult(i,j,k,n) = (poc(i,j,k) * kpoc(i,j,k,n) +  (caco3(i,j,k) * kcaco3(i,j,k,n)))
                   mult(i,j,k,n) = (poc(i,j,k) * kpoc(i,j,k,n) +  (caco3(i,j,k) * kcaco3(i,j,k,n))+  opal(i,j,k)  * kdsi(i,j,k,n))
                   neod_part(i,j,k,n) = neod_diss(i,j,k,n) * mult(i,j,k,n) *tms(i,j,k)

                   !if (neod_part(i,j,k,n)<0) then
                     !neod_diss(i,j,k,n) = 0._dp
                    ! write(*,*) neodymium(i,j,k,n), mult(i,j,k,n)*tms(i,j,k)
                   !endif

                 end do
             end do
          end do
       end do
     
       
       !caco3_SUM = 0
       !poc_SUM = 0
       !do i=1,imax_loc       
       !   do j=1,jmax_loc
       !       k=kmax_loc-7 
             !do k=1,kmax_loc   
       !         caco3_SUM = caco3_SUM + fcaco3(i,j,k)*dxc1(i,j)*dxc2(i,j)*tms(i,j,k)
       !         poc_SUM = poc_SUM + fpoc(i,j,k)*dxc1(i,j)*dxc2(i,j)*tms(i,j,k)
              !end do
       !   end do
       !end do

       !write (*,*) "!!!!!!!!!!!!!!!"
       !write(*,*) "caco3 MIN MAX", MINVAL(MINVAL(MINVAL(caco3,dim=1),dim=1),dim=1), MAXVAL(MAXVAL(MAXVAL(caco3,dim=1),dim=1),dim=1)
       !write (*,*) "fcaco3 SUM", caco3_SUM, "Tmol/j"
       !write (*,*) "!!!!!!!!!!!!!!!"

       !write (*,*) "!!!!!!!!!!!!!!!"
       !write(*,*) "poc MIN MAX", MINVAL(MINVAL(MINVAL(poc,dim=1),dim=1),dim=1), MAXVAL(MAXVAL(MAXVAL(poc,dim=1),dim=1),dim=1)
       !write (*,*) "fpoc SUM", poc_SUM, "Tmol/j"
       !write (*,*) "!!!!!!!!!!!!!!!"


#endif                    

!      write(*,*) "POC val", MINVAL(MINVAL(poc,dim=1),dim=1), MAXVAL(MAXVAL(poc,dim=1),dim=1)
      !write(*,*) "POC surf?", poc(:,:,ks2)

      end subroutine reversible_scavenging
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! DOWNWARD FLUX CALCULATION
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
      subroutine downward_flux
      
#if ( NEOD == 1 )
      integer :: i,j,k,n    
       do n=1,nb_trac
          !do k=2,kmax_loc    
          do k=19,1,-1    
             do j=1,jmax_loc
                 do i=1,imax_loc   
                    z_flux(i,j,k,n) = ws * (neod_part(i,j,k+1,n))*tms(i,j,k+1) ! g/m3 * m/s --> g/m2/s
                    !write(*,*) "diss, part, diss+part, ws, z_flux"
                    !write(*,*) neod_diss(i,j,k,n), neod_part(i,j,k,n), (neod_diss(i,j,k,n) + neod_part(i,j,k,n)), ws, z_flux(i,j,k,n)
                    !write(*,*) "z_flux",i,j,k, z_flux(i,j,k,n), "(g/m2/s)"
                    !write(*,*) neod_diss(i,j,k,n), neod_part(i,j,k,n),  z_flux(i,j,k,n)
                    !if (z_flux(i,j,k,n)< 0) then
                    !if (z_flux(i,j,k,n)>1E-10) then
                    ! write(*,*) i, j, k, n, neod_part(i,j,k+1,n),  z_flux(i,j,k,n)
                    !end if
                 end do
             end do
          end do
       end do
       
      do n=1,nb_trac
          !do k=1,kmax_loc
          do k=20,1,-1              
             do j=1,jmax_loc
                 do i=1,imax_loc   
                    z_neod(i,j,k,n) = ((z_flux(i,j,k,n) - z_flux(i,j,k-1,n)) /dz(k))*tms(i,j,k) ! g/m2/s --> g/m3/s
                    !if (z_neod(i,j,k,n)>1E-10) then
                    ! write(*,*) z_neod(i,j,k,n), z_flux(i,j,k,n), z_flux(i,j,k+1,n)
                    !endif
                    !write(*,*) "terme1-terme2, dz",i,j,k, (z_flux(i,j,k,n) - z_flux(i,j,k+1,n)), dz(k)
                    !write(*,*) "...",i,j,k, (z_flux(i,j,k,n) - z_flux(i,j,k+1,n)) /dz(k)
                    !write(*,*) "z_neod",i,j,k, z_neod(i,j,k,n), "(g/m3/s)"
                    !write(*,*) "terme1, terme2, sous, dz, zneod"
                    !write(*,*) z_flux(i,j,k,n), z_flux(i,j,k+1,n), (z_flux(i,j,k,n) - z_flux(i,j,k+1,n)), dz(k), z_neod(i,j,k,n)
                    !write(*,*) "zneod",i,j,k, z_neod(i,j,k,n), "(g/m3/s)"
                    !if (z_neod(i,j,k,n)>1E-10) then
                    ! write (*,*) "!!!!!!!!!!!!!!!"
                    ! write(*,*) i, j, k, n, z_flux(i,j,k,n), z_flux(i,j,k+1,n), (z_flux(i,j,k,n) - z_flux(i,j,k-1,n)), dz(k), ((z_flux(i,j,k,n) - z_flux(i,j,k-1,n)) /dz(k))*tms(i,j,k)
                    ! write (*,*) "!!!!!!!!!!!!!!!"
                    !end if
                    !write (*,*) z_neod(i,j,k,144)
                    !write (*,*) z_neod(i,j,k,144)*tms(i,j,k)
                    !write (*,*) z_neod(i,j,k,144)*tms(i,j,k)*86400
                    !write (*,*) "!!!!!!!!!!!!!!!"
                    !write(*,*) (z_flux(i,j,k,n) - z_flux(i,j,k+1,n)) /dz(k)
                    !write(*,*) z_neod(i,j,k,n)
                    !write (*,*) "!!!!!!!!!!!!!!!"
              end do
             end do
          end do
       end do
       
       z_neod_SUM = 0.
       !z_flux_source=0.
       !z_flux_sumZ=0.
       do k=1,kmax_loc       
              do j=1,jmax_loc
                     do i=1,imax_loc   
                        !z_flux_source= z_flux_source+ z_flux*dxc1(i,j)*dxc2(i,j)*tms(i,j,ks2)*86400*360
                         !write (*,*) "z_neod(i,j,k,nd144)", z_neod(i,j,k,nd144)  
                         !z_neod_SUM = z_neod_SUM + z_neod(i,j,k,nd144)     
                        z_neod_SUM = z_neod_SUM + ((z_neod(i,j,k,nd144)/0.238))*tms(i,j,k)*dxc1(i,j)*dxc2(i,j)*86400
                        !write(*,*) "z_neod(i,j,k,nd144)", z_neod(i,j,k,nd144)
                        !write(*,*) "z_neod_SUM", z_neod_SUM
                     end do
              end do
       end do
       
       !write (*,*) "!!!!!!!!!!!!!!!"
       !write(*,*) "z_neod_SUM per model time step", z_neod_SUM
       !write (*,*) "!!!!!!!!!!!!!!!"
       

#endif
      end subroutine downward_flux
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! STEP
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine neodymium_step
         integer :: i,j,k,n
         call dust_source_calculation
         call river_source_calculation
         
         write(*,*) "assignation to the tracer variable"
         write(*,*) "surface fluxes"
         
         neodymium(:,:,ks2,nd143) = neodymium(:,:,ks2,nd143) + ((dust_source_143Nd + river_source_143Nd)/dz(ks2))*tms(:,:,ks2)*86400
         neodymium(:,:,ks2,nd144) = neodymium(:,:,ks2,nd144) + ((dust_source_144Nd + river_source_144Nd)/dz(ks2))*tms(:,:,ks2)*86400 

         call benthic_flux_calculation
         write(*,*) "Benthic flux"

         do i=1,imax_loc
              do j=1,jmax_loc
               neodymium(i,j,kfs(i,j),nd143) = neodymium(i,j,kfs(i,j),nd143) + (benthic_source_143Nd(i,j)/dz(kfs(i,j))) &
                                              *tms(i,j,kfs(i,j))*86400
               neodymium(i,j,kfs(i,j),nd144) = neodymium(i,j,kfs(i,j),nd144) + (benthic_source_144Nd(i,j)/dz(kfs(i,j))) &
                                              *tms(i,j,kfs(i,j))*86400
              end do
         end do   
         
         write(*,*) "application of the reversible scavenging"
         call reversible_scavenging
         
         write(*,*) "calculation of the zflux"
         call downward_flux
         neodymium(:,:,:,nd143)= neodymium(:,:,:,nd143) + z_neod(:,:,:,nd143)*tms(:,:,:)*86400
         neodymium(:,:,:,nd144)= neodymium(:,:,:,nd144) + z_neod(:,:,:,nd144)*tms(:,:,:)*86400


          !Negative concentration values are prohibited 
         do n=1,nb_trac
          do k=1,kmax_loc
              do i=1,imax_loc       
                do j=1,jmax_loc
                  if (neodymium(i,j,k,n) < 0) then
                  !write (*,*) neodymium(i,j,k,nd143)
                     neodymium(i,j,k,n) = 0._dp
                  endif 
                end do
             end do
          end do  
         end do   


         Nd_inventory =0.
         Nddiss_inventory =0.
         Ndpart_inventory =0.
         do k=1,kmax_loc
           do i=1,imax_loc       
             do j=1,jmax_loc
               Nd_inventory = Nd_inventory + (neodymium(i,j,k,nd144)/0.238)*tms(i,j,k)*dxc1(i,j)*dxc2(i,j)*dz(k)
               Nddiss_inventory = Nddiss_inventory + (neod_diss(i,j,k,nd144)/0.238)*tms(i,j,k)*dxc1(i,j)*dxc2(i,j)*dz(k)
               Ndpart_inventory = Ndpart_inventory + (neod_part(i,j,k,nd144)/0.238)*tms(i,j,k)*dxc1(i,j)*dxc2(i,j)*dz(k)
             end do
          end do
        end do               
        
       write (*,*) "!!!!!!!!!!!!!!!"
       write (*,*) "TOTAL Nd inventory", Nd_inventory*0.000000000001, "Tg"
       write (*,*) "TOTAL Nd inventory", Nd_inventory*0.000000001, "Gg"
       write (*,*) "!!!!!!!!!!!!!!!"

       write (*,*) "!!!!!!!!!!!!!!!"
       write (*,*) "DISS Nd inventory", Nddiss_inventory*0.000000000001, "Tg"
       write (*,*) "DISS Nd inventory", Nddiss_inventory*0.000000001, "Gg"
       write (*,*) "!!!!!!!!!!!!!!!"      
       
       write (*,*) "!!!!!!!!!!!!!!!"
       write (*,*) "PART Nd inventory", Ndpart_inventory*0.000000000001, "Tg"
       write (*,*) "PART Nd inventory", Ndpart_inventory*0.000000001, "Gg"
       write (*,*) "!!!!!!!!!!!!!!!"          

       write (*,*) "!!!!!!!!!!!!!!!"
       write (*,*) "Nd removal inventory", z_neod_SUM*0.000000000001, "Tg"
       write (*,*) "Nd removal inventory", z_neod_SUM*0.000000001, "Gg"
       write (*,*) "!!!!!!!!!!!!!!!"       
        
!       write (*,*) "!!!!!!!!!!!!!!!"
!       write (*,*) "TOTAL Nd inventory minus removal", (Nd_inventory-z_neod_SUM)*0.000000000001, "Tg"
!       write (*,*) "TOTAL Nd inventory minus removal", (Nd_inventory-z_neod_SUM)z_neod_SUM*0.000000001, "Gg"
!       write (*,*) "!!!!!!!!!!!!!!!"       

      end subroutine neodymium_step

end module neodymium_mod
