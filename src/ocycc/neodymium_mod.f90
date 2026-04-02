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
       use bloc_mod, only: aire
       use loveclim_transfer_mod, only: ZX
       use declars_mod, only: JX
       use dynami_mod, only: dxc1, dxc2
       use global_constants_mod, only: dp, str_len, ip
       use bloc0_mod, only: fpoc => fPOC_flx_clio !TmolsC.m-2.timestep-1
       use bloc0_mod, only: fcaco3 => fCAL_flx_clio !TmolsC.m-2.timestep-1
       use ncio
       

       implicit none

       !private :: dust_source_calculation
       public  :: dust_source_calculation
       private :: river_source_calculation
       !private :: sediment_flux_calculation
       public   :: sediment_flux_calculation
       !private :: boundary_flux_calculation
       public   :: boundary_flux_calculation
       private :: reversible_scavenging
       private :: downward_flux
       public  :: neodymium_init
       public  :: neodymium_step
       public  :: init_source_netcdf
       public  :: coeff_init ! private? maybe merge with neodymium_init, check with dmr

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       integer, parameter :: flag_dust_eNd = 1      ! 0: inactivated, 1:activated
       integer, parameter :: flag_river_eNd = 1     ! 0: inactivated, 1:activated
       integer, parameter :: flag_sediment_eNd = 2  ! 0: inactivated, 1:activated, 2: mixed (sediment map substracted by boundary map)
       integer, parameter :: flag_boundary_eNd = 1  ! 0: inactivated, 1:activated
       integer, parameter :: flag_particule = 2     ! 0: no particle/reversible scavenging, 1: particles prescribed from NEMO, 2: particles from ocycc (no opal), 3: particles from ocycc + NEMO opal
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

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

! tva Definition of the input tables for reading the isotopic ratio (eNd) and the Nd concentration of rivers in ppt on the CLIO grid (OLD)
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
! tva Definitations for the sediment source component
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! tva Definition of the input tables for reading the isotopic ratio (eNd) of the sediment
       real(dp), dimension(120,65) :: ir_sediment
       real(dp), dimension(120,65) :: seafloor_without_marge_mask
! tva Same as above after longitude transformation (halo point) and unit conversion
       real(dp), dimension(imax_loc,jmax_loc) :: ir_sediment_transform
       real(dp), dimension(imax_loc,jmax_loc) :: seafloor_without_marge_mask_transform
       real(dp) ::  seafloor_without_marge_area, seafloor_area
! tva Definition of the Nd sediment global flux in g/yr (TMP for the implementation: Sediment flux from 3.0*10+9 gram of Nd /yr (best sim Robinson et al. 2022) to g/s)
! tva IN CASE OF FULL SEDIMENT SOURCE
! 2be put later in a namelist or whatever for tuning and/or space phase parameterization!!!!!!!!!
       real(dp), parameter :: sediment_source_Nd_global = 3e9!Flux in g per year
! tva Definition of the Nd sediment field tables
       real(dp), dimension(imax_loc,jmax_loc) :: sediment_source_Nd, sediment_source_144Nd, sediment_source_143Nd
!DIAG
       real(dp) :: sediment_source_Nd_tot

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
! tva intermediate variable for the Nd removal       
       real(dp), dimension(imax_loc,jmax_loc,nb_trac) :: z_bottom_store 

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! tva Definitations for the boundary source component
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real(dp), dimension(120,65,20) :: marge, conc_boundary, eps_boundary
       real(dp), dimension(imax_loc,jmax_loc,kmax_loc) :: marge_transform, conc_boundary_transform, eps_boundary_transform
! tva Definition of the sediment boundary source in g/s (Not sure if in gNd/s)
       real(dp), parameter :: sediment_boundary = 350
! tva Definition of the coefficient pour exprimer le dlux en g/m2/s (depend de la grille)
       real(dp), parameter :: coeff_boundary = 17.0922 !80.2546 !17.0922 for NEMO
       real(dp) :: denitide, expide
       real(dp), dimension(imax_loc,jmax_loc,kmax_loc) :: bathy, boundary_source_Nd, boundary_source_144Nd, boundary_source_143Nd
       real(dp) :: boundary_source_Nd_tot


        !dmr&tva --- To write the output text file
        integer(ip) :: eNdinv_id
        !dmr&tva --- To write the output text file



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


        open(newunit=eNdinv_id,file='outputdata/carbon/eNd_inventory.txt',form='formatted')


      end subroutine neodymium_init
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! INIT OF NEODYMIUM SOURCE (NECTDF READING)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine init_source_netcdf
      integer :: i,j,k
      CHARACTER*200 :: filin1, filin2
      CHARACTER*200 :: filin3, filin4, filin5
      CHARACTER*200 :: filin6, filin65, filin66, filin67
      CHARACTER*200 :: filin7
      
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

       !quick check of crazy values
       write(*,*) "min max ir_dust_transform", minval(ir_dust_transform), maxval(ir_dust_transform)
       write(*,*) "min max ir_dust_transform", minval(dust_flux_transform), maxval(dust_flux_transform)
       
!!RIVERS
       ! on CLIO
       filin3="inputdata/neodymium/River_eNd_regridded_on_CLIO_remapnn.nc"
       !filin3="inputdata/neodymium/out-testCLIO3-eps_rivieres-surf.nc"
       call nc_read(filin3,"river_end", ir_river)
       !call nc_read(filin3,"eps_rivieres", ir_river)

       filin4="inputdata/neodymium/River_Nd_ppt_regridded_on_CLIO_remapnn.nc"
       !filin4="inputdata/neodymium/out-testCLIO3-conc_rivieres-surf.nc"
       call nc_read(filin4,"river_nd_ppt", conc_river)
       !call nc_read(filin4,"conc_rivieres", conc_river)
       
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
          end do
       end do
       
       ir_river_transform (2:121,:)= ir_river_tmp(1:120,:)
       ir_river_transform (1,:)=ir_river_tmp(120,:)
       ir_river_transform (122,:)=ir_river_tmp(1,:)

       conc_river_transform (2:121,:)= conc_river_tmp(1:120,:)
       conc_river_transform (1,:)=conc_river_tmp(120,:)
       conc_river_transform (122,:)=conc_river_tmp(1,:)  

       !quick check of crazy values
       write(*,*) "min max ir_river_transform", minval(ir_river_transform), maxval(ir_river_transform)
       write(*,*) "min max conc_river_transform", minval(conc_river_transform), maxval(conc_river_transform)

       

!SEDIMENT FLUX
      
       if (flag_sediment_eNd == 1) then
              filin6="inputdata/neodymium/globalSedimentEpsilon_regridded_on_CLIO_remapnn.nc"
              call nc_read(filin6,"globalsedimentepsilonnd", ir_sediment)
       endif 

       if (flag_sediment_eNd == 2) then
              filin6="inputdata/neodymium/globalSedimentEpsilon_margemasked.nc"
              call nc_read(filin6,"seafloornd_without_margins", ir_sediment)
              call nc_read(filin6,"seafloornd_mask", seafloor_without_marge_mask)

              seafloor_without_marge_mask_transform(2:121,:)= seafloor_without_marge_mask(1:120,:)
              seafloor_without_marge_mask_transform(1,:)=seafloor_without_marge_mask(120,:)
              seafloor_without_marge_mask_transform(122,:)=seafloor_without_marge_mask(1,:)  
       endif 

       ir_sediment_transform(2:121,:)= ir_sediment(1:120,:)
       ir_sediment_transform(1,:)=ir_sediment(120,:)
       ir_sediment_transform(122,:)=ir_sediment(1,:)  
       
       !quick check of crazy values
       write(*,*) "min max ir_sediment_transform 1", minval(ir_sediment_transform), maxval(ir_sediment_transform)
       ! write(*,*) "min max seafloor_without_marge_mask_transform 1", minval(seafloor_without_marge_mask_transform), maxval(seafloor_without_marge_mask_transform)
       
       !quick fix
       where (ir_sediment_transform > 1000._dp)
              ir_sediment_transform = 0._dp
       end where

       where (seafloor_without_marge_mask_transform > 1000._dp)
              seafloor_without_marge_mask_transform = 0._dp
       end where

       write(*,*) "min max ir_sediment_transform 2", minval(ir_sediment_transform), maxval(ir_sediment_transform)
       ! write(*,*) "min max seafloor_without_marge_mask_transform 2", minval(seafloor_without_marge_mask_transform), maxval(seafloor_without_marge_mask_transform)


!BOUNDARY CONTINENTAL MARGIN
       filin65="inputdata/neodymium/marge_clio.nc"
       call nc_read(filin65,"marge", marge)

       marge_transform(2:121,:,:)= marge(1:120,:,:)
       marge_transform(1,:,:)=marge(120,:,:)
       marge_transform(122,:,:)=marge(1,:,:)
       
       filin66="inputdata/neodymium/marge_clio_Nd_epsNd.nc"
       call nc_read(filin66,"Ndconcmarge", conc_boundary)

       conc_boundary_transform(2:121,:,:)= conc_boundary(1:120,:,:)
       conc_boundary_transform(1,:,:)=conc_boundary(120,:,:)
       conc_boundary_transform(122,:,:)=conc_boundary(1,:,:)
       
       filin67="inputdata/neodymium/marge_clio_Nd_epsNd.nc"
       call nc_read(filin67,"epsNdmarge", eps_boundary)

       eps_boundary_transform(2:121,:,:)= eps_boundary(1:120,:,:)
       eps_boundary_transform(1,:,:)=eps_boundary(120,:,:)
       eps_boundary_transform(122,:,:)=eps_boundary(1,:,:)

       do i = 1, imax_loc
              do j = 1, jmax_loc
                  do k = 1, kmax_loc
                      if (isnan(conc_boundary_transform(i, j, k))) then
                            conc_boundary_transform(i, j, k) = 0.0
                      end if
                  end do
              end do
       end do
       
       !quick check of crazy values
       write(*,*) "min max marge 1", minval(marge_transform), maxval(marge_transform)
       write(*,*) "min max conc_boundary_transform 1", minval(conc_boundary_transform), maxval(conc_boundary_transform)
       write(*,*) "min max eps_boundary_transform", minval(eps_boundary_transform), maxval(eps_boundary_transform)

       !quick fix
       where (conc_boundary_transform < 0._dp)
              conc_boundary_transform = 0._dp
       end where

       where (marge_transform > 1000._dp)
              marge_transform = 0._dp
       end where

       write(*,*) "min max marge 2", minval(marge_transform), maxval(marge_transform)
       write(*,*) "min max conc_boundary_transform 2", minval(conc_boundary_transform), maxval(conc_boundary_transform)

!READ AND CONVERT PARTICLE CONC FIELD IN CASE OF PARTICLE PRESCRIPTION

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
       particles_prescribed_g_g(:,:,:,nopal) = ((((particles_prescribed(:,:,:,nopal))/1000)*60)/1025000)*tms
       particles_prescribed_g_g(:,:,:,npoc) = ((((particles_prescribed(:,:,:,npoc))/1000)*32.7)/1025000)*tms
       particles_prescribed_g_g(:,:,:,ncaco3) = ((((particles_prescribed(:,:,:,ncaco3))/1000)*100)/1025000)*tms
 
       if (flag_particule == 1) then
              poc = particles_prescribed_g_g(:,:,:,npoc)
              caco3 = particles_prescribed_g_g(:,:,:,ncaco3)
              opal = particles_prescribed_g_g(:,:,:,nopal)

              WRITE(*,*)'min max poc', minval(poc), maxval(poc)
              WRITE(*,*)'min max caco3', minval(caco3), maxval(caco3)
              WRITE(*,*)'min max opal', minval(opal), maxval(opal)

       end if

       if (flag_particule == 3) then
              opal = particles_prescribed_g_g(:,:,:,nopal)
              WRITE(*,*)'min max opal', minval(opal), maxval(opal)
       end if       
       
      end subroutine init_source_netcdf
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! INITIALIZATIONS OF TRACER PARTITION COEFFICIENTS
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      

      subroutine coeff_init
!          conversion factor partition coefficient: micro_mol/l-->g/l (water density)
!          ----------------------------------------------------------------------------
      integer :: i,j,k,n
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
            !kpoc(i,j,k,nd143) = 457800 * 6! kpoc_143Nd * Zpoc / Zwater
            !kpoc(i,j,k,nd144) = 457800* 6! kpoc_144Nd * Zpoc / Zwater
            kpoc(i,j,k,nd143) = kpoc_143Nd * Zpoc / Zwater
            kpoc(i,j,k,nd144) = kpoc_144Nd * Zpoc / Zwater
            !kgoc (ji,jj,jk,1) = kgoc_143 * Zgoc / Zwater
            !kgoc (ji,jj,jk,2) = kgoc_144 * Zgoc / Zwater            
            !kdsi (i,j,k,nd143) = 2160* 6 !kdsi_143Nd * Zdsi / Zwater
            !kdsi (i,j,k,nd144) = 2160* 6 !kdsi_144Nd * Zdsi / Zwater
            !kcaco3(i,j,k,nd143) = 16000* 6!kcaco3_143Nd * Zcaco3 / Zwater
            !kcaco3(i,j,k,nd144) = 16000* 6!kcaco3_144Nd * Zcaco3 / Zwater
            kdsi(i,j,k,nd143) = kdsi_143Nd * Zdsi / Zwater
            kdsi(i,j,k,nd144) = kdsi_144Nd * Zdsi / Zwater
            kcaco3(i,j,k,nd143) = kcaco3_143Nd * Zcaco3 / Zwater
            kcaco3(i,j,k,nd144) = kcaco3_144Nd * Zcaco3 / Zwater     
            !klitho (ji,jj,jk,1) = klitho_143 * Zlitho / Zwater
            !klitho (ji,jj,jk,2) = klitho_144 * Zlitho / Zwater
           end do
         end do
       end do

       WRITE(*,*) "kpoc(1,1,1,nd143)", kpoc(1,1,1,nd143)
       WRITE(*,*) "kpoc(1,1,1,nd143)", kpoc(1,1,1,nd143)   
       WRITE(*,*) "kdsi(1,1,1,nd143)", kdsi(1,1,1,nd143)
       WRITE(*,*) "kdsi(1,1,1,nd144)", kdsi(1,1,1,nd144)
       WRITE(*,*) "kcaco3(1,1,1,nd143)", kcaco3(1,1,1,nd143)
       WRITE(*,*) "kcaco3(1,1,1,nd144)", kcaco3(1,1,1,nd144)   
                        
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
       
       dust_source_Nd = dust_flux_transform_g_m2_s*conc_dust*beta_dust
       dust_source_144Nd = 0.238*dust_source_Nd !the abundance of 144Nd is 23.8 per cent of the total Nd
       dust_source_143Nd = ((ir_dust_transform/10000)+1)*dust_source_144Nd*CHUR ! Using the eNd definition to calculate 143Nd from 144Nd (Arsouze et al. 2009)


!DIAG
       dust_source_Nd_tot =0._dp
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
       conc_river_transform_g = conc_river_transform*1.0E-9*1.0E-3 ! converting Nd ppt into g/kg (1.0E-9)... and convert kg to g (1.0E-3)

       river_source_Nd =  conc_river_transform_g  * river_remob * fwruno_g_m2_s
       river_source_144Nd = 0.238*river_source_Nd !the abundance of 144Nd is 23.8 per cent of the total Nd
       river_source_143Nd = ((ir_river_transform/10000)+1)*river_source_144Nd*CHUR ! Using the eNd definition to calculate 143Nd from 144Nd (Arsouze et al. 2009)  

!DIAG
       river_source_Nd_tot =0._dp
       area_tot=0._dp
       area2_tot=0._dp
       fwruno_tot =0._dp
       fwruno_area_tot =0._dp
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
       
!       write (*,*) "!!!!!!!!!!!!!!!"
!       !write (*,*) "total freshwater to the ocean", fwruno_area_tot*(-1)*360*0.001*0.000000001, "(km3/yr)"
!       write (*,*) "total freshwater to the ocean", fwruno_area_tot*360*86400*0.000001*0.000000001, "(km3/yr)"
!       write (*,*) "!!!!!!!!!!!!!!!" 

       write (*,*) "!!!!!!!!!!!!!!!"
       write (*,*) "total Nd river flux", river_source_Nd_tot*0.000000001, "(10^9)g/yr"
       write (*,*) "!!!!!!!!!!!!!!!"
       
!       write (*,*) "!!!!!!!!!!!!!!!"
!       write (*,*) "surface of the ocean", area_tot, area2_tot
!       write (*,*) "!!!!!!!!!!!!!!!"      
       
       write (*,*) "!!!!!!!!!!!!!!!"
       write (*,*) "total freshwater to the ocean", fwruno_tot, ("-kg/m2/day")
       write (*,*) "!!!!!!!!!!!!!!!"         
       
!       write (*,*) "!!!!!!!!!!!!!!!"
!       write (*,*) "total freshwater to the ocean", fwruno_area_tot, ("-kg/day")
!       write (*,*) "!!!!!!!!!!!!!!!" 


#else
       river_source_144Nd = 0._dp
       river_source_143Nd = 0._dp
#endif
      end subroutine river_source_calculation

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! SEDIMENT SOURCE CALCULATION
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine sediment_flux_calculation
#if ( NEOD == 1 )
      integer :: i,j    
       

       if (flag_sediment_eNd == 1 ) then
              seafloor_area =0._dp
              do i=1,imax_loc
                     do j=1,jmax_loc
                            seafloor_area = seafloor_area + tms(i,j,kfs(i,j))*dxc1(i,j)*dxc2(i,j)
                     end do
              end do

              sediment_source_Nd = sediment_source_Nd_global/(86400*365.25)! Sediment flux from g/yr to g/s
              !WRITE(*,*) "sediment_source_Nd tot", minval(sediment_source_Nd), maxval(sediment_source_Nd)
              sediment_source_Nd = sediment_source_Nd/seafloor_area! Sediment flux divided by the total surface of seafloor (364012054606928. m2 in CLIO)
              !WRITE(*,*) "sediment_source_Nd tot", minval(sediment_source_Nd), maxval(sediment_source_Nd)
       end if 

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       if (flag_sediment_eNd == 2 ) then 
              seafloor_without_marge_area =0._dp
              do i=1,imax_loc
                     do j=1,jmax_loc
                            seafloor_without_marge_area = seafloor_without_marge_area + tms(i,j,kfs(i,j))& 
                                    *seafloor_without_marge_mask_transform(i,j)*dxc1(i,j)*dxc2(i,j)
                     end do
              end do

              sediment_source_Nd = (sediment_source_Nd_global - boundary_source_Nd_tot)/(86400*360)
              !WRITE(*,*) "sediment_source_Nd_global", sediment_source_Nd_global
              !WRITE(*,*) "boundary_source_Nd_tot", boundary_source_Nd_tot
              !WRITE(*,*) "sediment_source_Nd tot", minval(sediment_source_Nd), maxval(sediment_source_Nd)
              sediment_source_Nd = sediment_source_Nd/(seafloor_without_marge_area)
              !WRITE(*,*) "seafloor_area, seafloor_without_marge_area", seafloor_area, seafloor_without_marge_area
              !WRITE(*,*) "sediment_source_Nd tot", minval(sediment_source_Nd), maxval(sediment_source_Nd)
       end if 

       sediment_source_144Nd = 0.238*sediment_source_Nd !the abundance of 144Nd is 23.8 per cent of the total Nd
       sediment_source_143Nd = ((ir_sediment_transform/10000)+1)*sediment_source_144Nd*CHUR ! Using the eNd definition to calculate 143Nd from 144Nd (Arsouze et al. 2009)

!DIAG
       sediment_source_Nd_tot =0._dp
       do i=1,imax_loc
          do j=1,jmax_loc
            !write(*,*) "kfs(i,j)", kfs(i,j)
            if (flag_sediment_eNd == 1 ) then
              sediment_source_Nd_tot=sediment_source_Nd_tot+sediment_source_Nd(i,j)*dxc1(i,j)*dxc2(i,j)*tms(i,j,kfs(i,j))*86400*360
            end if 
            if (flag_sediment_eNd == 2 ) then
              sediment_source_Nd_tot=sediment_source_Nd_tot+sediment_source_Nd(i,j)*dxc1(i,j)*dxc2(i,j)&
                      *seafloor_without_marge_mask_transform(i,j)*tms(i,j,kfs(i,j))*86400*360
            end if 
          end do
       end do

       !if (flag_sediment_eNd == 1 ) then 
              !write (*,*) "!!!!!!!!!!!!!!!"
              !write (*,*) "seabed surf check", seaflorr_area
              !write (*,*) "!!!!!!!!!!!!!!!"
       !end if

       !if (flag_sediment_eNd == 2 ) then 
       !       write (*,*) "!!!!!!!!!!!!!!!"
       !       write (*,*) "seafloor surf check", seafloor_surf_check
       !       write (*,*) "!!!!!!!!!!!!!!!"              
       !end if

       write (*,*) "!!!!!!!!!!!!!!!"
       write (*,*) "total Nd sediment flux", sediment_source_Nd_tot*0.000000001, "(10^9)g/yr"
       write (*,*) "!!!!!!!!!!!!!!!"


#else
       sediment_source_144Nd = 0._dp
       sediment_source_143Nd = 0._dp
#endif
      end subroutine sediment_flux_calculation

      subroutine boundary_flux_calculation
#if ( NEOD == 1 )
      integer :: i,j,k,k_rev
      integer, parameter :: kmax_plus_one = 21
      real(dp) :: ZX_flipped(kmax_plus_one)
      real(dp), dimension(imax_loc,jmax_loc) :: marge_transform_z_sum
      real(dp), dimension(imax_loc,jmax_loc) :: marge_elsewhere_mask
      real(dp) :: marge_transform_z_sum_sum, marge_elsewhere_mask_sum

       !SWAP ZX
       k_rev = 21
       DO k=1, kmax_plus_one
              ZX_flipped(k) = ZX(k_rev)
              k_rev = k_rev - 1
       END DO

       DO k = 1, kmax_loc
              DO j = 1, jmax_loc
                 DO i = 1, imax_loc
                   expide   = MIN( 8.,( ZX_flipped(k) / 500. )**(-1.5) )
                   denitide = -0.9543 + 0.7662 * LOG( expide ) - 0.235 * LOG( expide )**2
                   bathy(i,j,k) = (marge_transform(i,j,k)* tms(i,j,k)) * MIN( 1., EXP( denitide ) / 0.5 )
                   boundary_source_Nd(i,j,k) = sediment_boundary * coeff_boundary * bathy(i,j,k) * &
                            conc_boundary_transform(i,j,k) * 1e-6 * tms(i,j,k)/(dxc1(i,j)*dxc2(i,j))

                 END DO
              END DO
       END DO   

           boundary_source_144Nd = 0.238*boundary_source_Nd
           boundary_source_143Nd = ((eps_boundary_transform/10000)+1)*boundary_source_144Nd*CHUR

!DIAG
           boundary_source_Nd_tot =0.

           do k=1, kmax_loc
            do i=1,imax_loc
              do j=1,jmax_loc   
                boundary_source_Nd_tot=boundary_source_Nd_tot+boundary_source_Nd(i,j,k)*dxc1(i,j)*dxc2(i,j)*tms(i,j,k)*86400*360
              end do
            end do
           end do
    
           write (*,*) "!!!!!!!!!!!!!!!"
           write (*,*) "total Nd boundary flux", boundary_source_Nd_tot*0.000000001, "(10^9)g/yr"
           write (*,*) "!!!!!!!!!!!!!!!"

#endif
      end subroutine boundary_flux_calculation

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



      if (flag_particule == 2 .or. flag_particule == 3) then
       !CONVERT PARTICLE CONC FIELD FROM OCYCC PARTICLE FLUX
       fpoc_g_m2_day = fpoc*12.01*1E12 ! TmolsC to g: (1 x 10^12 moles/Tmol) * (12.01 grams/mole) --> flux from TmolsC/m2/day to g/m2/day
       poc_g_m3 = fpoc_g_m2_day/2.8  ! To get a concentration in g/m3: division by the settling speed (ws in m/s in ws/360 m/day or 1000/360) 
       poc = poc_g_m3 /1025000 !  division by the density of seawater to have g/g

       fcaco3_g_m2_day = fcaco3*100.09*1E12 ! TmolsC to g: (1 x 10^12 moles/Tmol) * (100.09 grams/mole) --> flux from TmolsC/m2/yr to g/m2/yr
       caco3_g_m3 = fcaco3_g_m2_day/2.8  ! To get a concentration in g/m3: division by the settling speed (ws in m/s in ws*360 m/day or 1000/360)
       caco3 = caco3_g_m3 /1025000     !  division by the density of seawater to have g/g

      end if


      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Calculation of [Nd] in the dissolve phase for each isotope and each particle type
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       do n=1,nb_trac
          do i=1,imax_loc       
             do j=1,jmax_loc
                 do k=1,kmax_loc   
                   if (neodymium(i,j,k,n)<0) then
                     neodymium(i,j,k,n) = 0._dp
                   endif

                   if (flag_particule == 2) then
                     denum(i,j,k,n) = (1. +  poc(i,j,k) * kpoc(i,j,k,n) + caco3(i,j,k)*kcaco3(i,j,k,n))
                   end if
                   if (flag_particule == 1 .or. flag_particule == 3) then                   
                     denum(i,j,k,n) = (1. +  poc(i,j,k) * kpoc(i,j,k,n) + caco3(i,j,k)*kcaco3(i,j,k,n)+opal(i,j,k)*kdsi(i,j,k,n))
                   end if

                   neod_diss(i,j,k,n) = (neodymium(i,j,k,n) / denum(i,j,k,n))*tms(i,j,k)

                 end do
             end do
          end do
       end do
        
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Calculation of [Nd] in the particulate phase for each isotope and each particle type
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       do n=1,nb_trac
          do i=1,imax_loc       
             do j=1,jmax_loc
                 do k=1,kmax_loc
                   if (neodymium(i,j,k,n)<0) then
                     neodymium(i,j,k,n) = 0._dp
                   endif

                   if (flag_particule == 2) then
                     mult(i,j,k,n) = (poc(i,j,k) * kpoc(i,j,k,n) +  (caco3(i,j,k) * kcaco3(i,j,k,n)))
                   end if
                   if (flag_particule == 1 .or. flag_particule == 3) then    
                     mult(i,j,k,n) = (poc(i,j,k) * kpoc(i,j,k,n) +  (caco3(i,j,k)*kcaco3(i,j,k,n))+opal(i,j,k)*kdsi(i,j,k,n))
                   end if

                   neod_part(i,j,k,n) = neod_diss(i,j,k,n) * mult(i,j,k,n) *tms(i,j,k)

                 end do
             end do
          end do
       end do


#endif                    


      end subroutine reversible_scavenging
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! DOWNWARD FLUX CALCULATION
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
      subroutine downward_flux
      
#if ( NEOD == 1 )
      integer :: i,j,k,n    
      z_flux = 0.
       do n=1,nb_trac
          !do k=2,kmax_loc    
          do k=19,1,-1    
             do j=1,jmax_loc
                 do i=1,imax_loc   
                    z_flux(i,j,k,n) = ws * (neod_part(i,j,k+1,n))*tms(i,j,k+1) ! g/m3 * m/s --> g/m2/s
                 end do
             end do
          end do
       end do
       
      do n=1,nb_trac
          !do k=1,kmax_loc
          do k=20,2,-1              
             do j=1,jmax_loc
                 do i=1,imax_loc   
                    z_neod(i,j,k,n) = ((z_flux(i,j,k,n) - z_flux(i,j,k-1,n)) /dz(k))*tms(i,j,k) ! g/m2/s --> g/m3/s
              end do
             end do
          end do
       end do
       
       do n=1,nb_trac        
              do j=1,jmax_loc
                     do i=1,imax_loc   
                        z_neod(i,j,1,n) = ((z_flux(i,j,1,n)) /dz(1))*tms(i,j,1) ! g/m2/s --> g/m3/s
                        z_bottom_store(i,j,n) = z_neod(i,j,1,n)
                        z_neod(i,j,1,n) = 0
                     end do
              end do
        end do
        

       z_neod_SUM = 0._dp
    
       do j=1,jmax_loc
              do i=1,imax_loc   
                     z_neod_SUM = z_neod_SUM + ((z_bottom_store(i,j,nd144)/0.238))*tms(i,j,1)*dxc1(i,j)*dxc2(i,j)*86400
              end do
       end do

       
       write (*,*) "!!!!!!!!!!!!!!!"
       write(*,*) "z_neod_SUM per model time step", z_neod_SUM
       write (*,*) "!!!!!!!!!!!!!!!"
       

#endif
      end subroutine downward_flux
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
! STEP
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine neodymium_step

         !dmr&tva --- To write the output text file
         !dmr&tva --- [NOTA] not clean to add things from an atmospheric module into the ocean part !!!
         use comemic_mod, only: iyear, day, iatm
         use comatm, only: dt
         !dmr&tva --- To write the output text file


         integer :: i,j,k,n

         if (flag_sediment_eNd == 1 .and. flag_boundary_eNd == 1) then
            print *, "FROM NEODYMIUM_MOD:"
            print *, "YOU CANNOT HAVE THE BOUNDARY SOURCE ALONG WITH THE FULL SEDIMENT SOURCE!!!"
            stop
         end if

         if (flag_dust_eNd == 1 .and. flag_river_eNd == 0) then
              !call dust_source_calculation -> called in ini
              neodymium(:,:,ks2,nd143) = neodymium(:,:,ks2,nd143) + ((dust_source_143Nd)/dz(ks2))*tms(:,:,ks2)*86400
              neodymium(:,:,ks2,nd144) = neodymium(:,:,ks2,nd144) + ((dust_source_144Nd)/dz(ks2))*tms(:,:,ks2)*86400 
         end if

         if (flag_dust_eNd == 0 .and. flag_river_eNd == 1) then
              call river_source_calculation
              neodymium(:,:,ks2,nd143) = neodymium(:,:,ks2,nd143) + ((river_source_143Nd)/dz(ks2))*tms(:,:,ks2)*86400
              neodymium(:,:,ks2,nd144) = neodymium(:,:,ks2,nd144) + ((river_source_144Nd)/dz(ks2))*tms(:,:,ks2)*86400 
         end if        

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
         if (flag_dust_eNd == 1 .and. flag_river_eNd == 1) then 
              !call dust_source_calculation -> called in ini
              call river_source_calculation
              neodymium(:,:,ks2,nd143)=neodymium(:,:,ks2,nd143)+((dust_source_143Nd+river_source_143Nd)/dz(ks2))*tms(:,:,ks2)*86400
              neodymium(:,:,ks2,nd144)=neodymium(:,:,ks2,nd144)+((dust_source_144Nd+river_source_144Nd)/dz(ks2))*tms(:,:,ks2)*86400 
         end if

         if (flag_boundary_eNd == 1) then 
              !call boundary_flux_calculation -> called in ini
              do i=1,imax_loc
                     do j=1,jmax_loc
                            do k = 1, kmax_loc
                            neodymium(i,j,k,nd143)=neodymium(i,j,k,nd143)+(boundary_source_143Nd(i,j,k)/dz(k))*tms(i,j,k)*86400
                            neodymium(i,j,k,nd144)=neodymium(i,j,k,nd144)+(boundary_source_144Nd(i,j,k)/dz(k))*tms(i,j,k)*86400 
                            end do
                     end do
              end do   
         end if

         if (flag_sediment_eNd == 1 .or. flag_sediment_eNd == 2) then 
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
              !call sediment_flux_calculation -> called in ini
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
              do i=1,imax_loc
                     do j=1,jmax_loc
                            neodymium(i,j,kfs(i,j),nd143) = neodymium(i,j,kfs(i,j),nd143) & 
                                    + (sediment_source_143Nd(i,j)/dz(kfs(i,j)))*tms(i,j,kfs(i,j))*86400
                            neodymium(i,j,kfs(i,j),nd144) = neodymium(i,j,kfs(i,j),nd144) &
                                    + (sediment_source_144Nd(i,j)/dz(kfs(i,j)))*tms(i,j,kfs(i,j))*86400       
                     end do
              end do   
         end if
        
         call reversible_scavenging
         
         call downward_flux
         neodymium(:,:,:,nd143)= neodymium(:,:,:,nd143) + z_neod(:,:,:,nd143)*tms(:,:,:)*86400
         neodymium(:,:,:,nd144)= neodymium(:,:,:,nd144) + z_neod(:,:,:,nd144)*tms(:,:,:)*86400


          !Negative concentration values are prohibited 
         do n=1,nb_trac
          do k=1,kmax_loc
              do i=1,imax_loc       
                do j=1,jmax_loc
                  if (neodymium(i,j,k,n) < 0._dp) then
                     neodymium(i,j,k,n) = 0._dp
                  endif 
                end do
             end do
          end do  
         end do   


         Nd_inventory =0._dp
         Nddiss_inventory =0._dp
         Ndpart_inventory =0._dp
         do k=1,kmax_loc
           do i=1,imax_loc       
             do j=1,jmax_loc
               Nd_inventory = Nd_inventory + (neodymium(i,j,k,nd144)/0.238)*tms(i,j,k)*dxc1(i,j)*dxc2(i,j)*dz(k)
               Nddiss_inventory = Nddiss_inventory + (neod_diss(i,j,k,nd144)/0.238)*tms(i,j,k)*dxc1(i,j)*dxc2(i,j)*dz(k)
               Ndpart_inventory = Ndpart_inventory + (neod_part(i,j,k,nd144)/0.238)*tms(i,j,k)*dxc1(i,j)*dxc2(i,j)*dz(k)
             end do
          end do
        end do               
        

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
        !dmr&tkv --- To write the output text file
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-
        write(eNdinv_id,110) iyear,int((day+0.5*dt)/(iatm*dt))+1,Nd_inventory, Nddiss_inventory, Ndpart_inventory, z_neod_SUM &
                           , river_source_Nd_tot

  110 format(i8,i8,f17.1,f17.1,f17.1,f17.1,f17.1)

        !dmr&tkv --- To write the output text file

!       write (*,*) "!!!!!!!!!!!!!!!"
!       write (*,*) "TOTAL Nd inventory", Nd_inventory*0.000000000001, "Tg"
!       write (*,*) "TOTAL Nd inventory", Nd_inventory*0.000000001, "Gg"
!       write (*,*) "!!!!!!!!!!!!!!!"

!       write (*,*) "!!!!!!!!!!!!!!!"
!       write (*,*) "DISS Nd inventory", Nddiss_inventory*0.000000000001, "Tg"
!       write (*,*) "DISS Nd inventory", Nddiss_inventory*0.000000001, "Gg"
!       write (*,*) "!!!!!!!!!!!!!!!"      
       
!       write (*,*) "!!!!!!!!!!!!!!!"
!       write (*,*) "PART Nd inventory", Ndpart_inventory*0.000000000001, "Tg"
!       write (*,*) "PART Nd inventory", Ndpart_inventory*0.000000001, "Gg"
!       write (*,*) "!!!!!!!!!!!!!!!"          

       write (*,*) "!!!!!!!!!!!!!!!"
       write (*,*) "Ndpart/Nddiss", Ndpart_inventory/Nddiss_inventory
       write (*,*) "!!!!!!!!!!!!!!!" 

!       write (*,*) "!!!!!!!!!!!!!!!"
!       write (*,*) "Nd removal inventory", z_neod_SUM*0.000000000001, "Tg"
!       write (*,*) "Nd removal inventory", z_neod_SUM*0.000000001, "Gg"
!       write (*,*) "!!!!!!!!!!!!!!!"       
        
!       !write (*,*) "!!!!!!!!!!!!!!!"
!       !write (*,*) "TOTAL Nd inventory minus removal", (Nd_inventory-z_neod_SUM)*0.000000000001, "Tg"
!       !write (*,*) "TOTAL Nd inventory minus removal", (Nd_inventory-z_neod_SUM)z_neod_SUM*0.000000001, "Gg"
!       !write (*,*) "!!!!!!!!!!!!!!!"          

      end subroutine neodymium_step

end module neodymium_mod
