!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of the lcm2ism coupling of GRISLI in iLOVECLIM.
!!      lcm2ism is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
!!      as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!!
!!      lcm2ism is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
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
!      MODULE: transfer_subgrid_ecb_mod
!
!>     @author  Aurelien F. Quiquet (afq) and Didier M. Roche (dmr)
!
!>     @brief This module transf_subgrid_ecb_mod is handling the transfer of topography from a given subgrid
!               to the ECbilt grid. It is a priori grid-specific to the particular grids given.
!               It has been built from the original transfer_grisli_ecb_mod.f90 (rev 953)
!
!>     @date Creation date: February, 04th, 2019
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : afq
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       module transfer_subgrid_ecb_mod

        ! --- afq&dmr, size of the subgrid
        use taillesGrilles, only: sgnx,sgny,sgnxm, sgnym, sgd

        implicit none

      ! NOTE_avoid_public_variables_if_possible

        public  :: subgrid_ecb_wrapper

        ! -- afq, private routine to get the gcm vert. level to use for an sub grid point:
        private  ::compute_weights_virt_levels
        ! -- afq, needs the private routine:
        private :: weights_virt_levels

        private  :: aggreg_one_field
        private  :: aggregSG2ECB

#if ( F_PALAEO == 3 )
        private  :: update_forcedISM
#endif

        private :: topoDiffECBsubgrid
        private :: update_topo_masks
        private :: masks_subgrid_where_flat

        double precision, dimension(sgnxm,sgnym,sgd), public :: weights_low_2d
        integer         , dimension(sgnxm,sgnym,sgd), public ::   index_low_2d

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: subgrid_ecb_wrapper
!
!>     @brief This subroutine allows updating the topography from the subgrid on ECBilt and the pertaining masks
!
!      DESCRIPTION:
!
!>     This subroutine is thought as a replacement of the old Interp_GRISLI2L, not a module and extended with new masks for
!!      the downscaling.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine subgrid_ecb_wrapper()

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Variables through module usage (would be delete in the end if all coming from the same module, this one!)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          use input_subgrid2L, only: topoSG, epaisSG, masqueSGice, areaSG, nbpointsSG, sub_grid_notflat, iEg, jEg
#if (ISM >= 2)
          use input_subgrid2L, only:  topoGRIS, masqueGRISice=>masqueGRIS, epaisglaceGRISnord
#endif
#if ( SMB_TYP == 2)
          use input_subgrid2L, only: annbiasSG
#endif

          use output_ECBilt,  only: masqueECB, topoECB, topdifECB, nb_stat
          use taillesGrilles, only: sgnx, sgny, iEcb, jEcb

#if ( CONSEAU == 1 )
          use varsCONSEAU_mod,      only: smbgrisnord, bmeltgrisnord, calgrisnord,                                &
                                          smbgris, bmeltgris, calgris, flux_gris2ecb, calgrisCLIO, calgrisCLIOms
          use routageEAU_mod,       only: eni, enj, mask_lnd
          use input_subgrid2L,      only: iCLIOg,jCLIOg
          use taillesGrilles,       only: CNX, CNY
          use varsCliotemp_mod,     only: ocn_mask
          use global_constants_mod, only: solar_day, days_year360d
          use comatm,               only: darea, tarea
          use bloc0_mod,            only: tms, ks2
          use dynami_mod,           only: area
          use const_mod,            only: zero
#endif
#if ( F_PALAEO_FWF == 1 )
          use varsCONSEAU_mod,      only: flux_gris2ecb
          use input_subgrid2L,      only: fluxFWF_SG_route
          use comatm,               only: darea
          use routageEAU_mod,       only: eni, enj, mask_lnd
#endif
#if ( DOWNSTS == 1 )
          use comemic_mod, only: fracto !we do not want to do any downscaling when there is only ocean in ECBilt

#if ( DOWNSCALING == 2 )
          use input_subgrid2L, only: max_nb_points, not_flat_subgrid, nbpointssg, area_sg_onecb            &
     &   , cast_ecb_subgrid_in_2d, cluster_var_subgrid_in_ecbilt_dble, cluster_var_subgrid_in_ecbilt_intg  &
     &   , weights_low_sg, index_low_sg, area_sg                                                           &  !clustered variables
     &   , latEcb,lonEcb,latSG,lonSG,nneigh, nexted, index_interpL2G, weights_interpL2G, sumweights_interpL2G & !interpo
     &   , nbpointsSG, topo_sg, difftopo_sg, epaisglace_sg
!#if (ISM >= 2)
!          use input_subgrid2L, only: difftopo_sg, epaisglace_sg
!#endif
#if ( SMB_TYP == 2 )
          use input_subgrid2L, only: annbias_sg
#endif
          use interpolate_mod, only: interpolate_init

          double precision, dimension(iEcb, jEcb, max_nb_points) :: subgrid_not_flat_sg
          integer :: n_dim, nbd
          logical :: success

#endif
          double precision, dimension(iEcb,jEcb)        :: masqueECB_loc,topoECB_loc
#endif

#if ( CONSEAU == 1 )
          integer :: i,j
          double precision, dimension(iEcb, jEcb) :: arouter, smbgrisECB, bmeltgrisECB
#endif
#if ( F_PALAEO_FWF == 1 )
          integer :: i,j
          double precision, dimension(iEcb, jEcb) :: fluxFWFECB_fice

          !real test_FWF_sum3
          real test_FWF_sum4
#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the subroutine starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          call update_topo_masks() ! update the masks

#if ( F_PALAEO == 3 )
!afq --- we hack the topo in order to use the prescribed one
          call update_forcedISM
#if ( F_PALAEO_FWF == 1 )
          call update_forced_FWF
#endif
          masqueSGice=masqueGRISice
          topoSG=topoGRIS
#else
          call update_topo_masks() ! update the masks
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Aggregate topography and masks back to the atmospheric grid ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!#if (ISM >= 2)
!          call aggreg_one_field(masqueGRISice,masqueECB,1)
!#endif
          masqueECB_loc(:,:) = 0.
          topoECB_loc(:,:) = 0.
          call aggreg_one_field(masqueSGice,masqueECB_loc,1)
          ! Beware that topoSG and topoECB are with zeroed-values on the oceans ...
          call aggreg_one_field(topoSG,topoECB_loc,1)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Computes the min and max value of the topography of the subgrid on the Atmos. grid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
          call topoDiffECBsubgrid(topoSG,sgnx,sgny,topdifECB(:,:,nb_stat),iEcb,jEcb,iEG,jEG,topdifECB(:,:,1))

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&afq  Retrieve the GRISLI water fluxes for conservation
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#if ( ISM >= 2 && CONSEAU ==  1 )
          smbgrisECB(:,:)=0.
          bmeltgrisECB(:,:)=0.
          calgrisCLIO(:,:)=0.
          smbgris(:,:,:)=0.
          smbgris(1:sgnx(sgd),1:sgny(sgd),1)=smbgrisnord(:,:)
          bmeltgris(:,:,:)=0.
          bmeltgris(1:sgnx(sgd),1:sgny(sgd),1)=bmeltgrisnord(:,:)
! afq -- grisli provides m3/yr, eventually we want m/s for ec_co2oc.f
          call aggreg_one_field(smbgris,smbgrisECB,0)
          call aggreg_one_field(bmeltgris,bmeltgrisECB,0)
! afq -- calving is transferred to the ocean grid... what about calving which is occuring on land due to inconsistencies 
!        between GRISLI and CLIO in terms of land sea mask?...
          call aggregISM2CLIO(calgrisnord,sgnx(1),sgny(1),calgrisCLIO,CNX,CNY,iCLIOg(:,:,1),jCLIOg(:,:,1),ocn_mask)
          arouter(:,:) = smbgrisECB(:,:) + bmeltgrisECB(:,:)
          flux_gris2ecb(:,:) = 0.
          do j=1,ubound(arouter,dim=2)
             do i=1,ubound(arouter,dim=1)
                if (mask_lnd(i,j).eq.0) then
                   flux_gris2ecb(i,j) = flux_gris2ecb(i,j) + arouter(i,j)                   
                else
                   flux_gris2ecb(eni(i,j),enj(i,j)) = flux_gris2ecb(eni(i,j),enj(i,j)) + arouter(i,j)
                endif
             enddo
          enddo
          
          do j=1,ubound(arouter,dim=2)
             do i=1,ubound(arouter,dim=1)
                flux_gris2ecb(i,j) = flux_gris2ecb(i,j)/(solar_day*days_year360d)/darea(i) 
             enddo
          enddo

          calgrisCLIO = calgrisCLIO/(solar_day*days_year360d) ! -> in m3/s

         !dmr&afq create a new variable for calving that is in m.s-1, coherent
         !        with flux_gris2ecb
          do j=1,CNY
             do i=1,CNX
                if (tms(i,j,ks2).gt.zero) then
                   calgrisCLIOms(i,j) = calgrisCLIO(i,j)/(area(i,j)*tms(i,j,ks2))
                endif
             enddo
          enddo

#endif

#if ( F_PALAEO_FWF == 1 )
! aggreg from grisli grid to ecbilt grid (sum values)
          call aggreg_one_field(fluxFWF_SG_route,fluxFWFECB_fice,0)
          !flux_gris2ecb(:,:) = fluxFWFECB_fice(:,:) !in m/s or m3/s
!nb to have unit change from m3/s -> m/s
          !do i =1, ubound(flux_gris2ecb, dim=1)
          ! flux_gris2ecb(i,:) = fluxFWFECB_fice(i,:)/darea(i)
          !enddo

!          test_FWF_sum3=0.0
!          do j =1, ubound(fluxFWFECB_fice, dim=2)
!           do i =1, ubound(fluxFWFECB_fice, dim=1)
!            test_FWF_sum3=test_FWF_sum3+fluxFWFECB_fice(i,j)
!            !test_FWF_sum3=test_FWF_sum3+flux_gris2ecb(i,j)
!           enddo
!          enddo
!          write(*,*) 'test_FWF_sum3 ', test_FWF_sum3*1e-6

!          do j =1, ubound(flux_gris2ecb, dim=2)
!           do i =1, ubound(flux_gris2ecb, dim=1)
!            if (.not.withini(eni(i,j),1,ubound(flux_gris2ecb,dim=1))) then
!              if (fluxFWFECB_fice(i,j).gt.epsilon(fluxFWFECB_fice(i,j))) then
!                write(*,*) 'GROS PROBLEME i',i,j,eni(i,j), mask_lnd(i,j)
!                write(*,*) 'fluxFWFECB_fice', fluxFWFECB_fice(i,j)
!              endif
!            endif
!            if (.not.withini(enj(i,j),1,ubound(flux_gris2ecb,dim=2))) then
!              if (fluxFWFECB_fice(i,j).gt.epsilon(fluxFWFECB_fice(i,j))) then
!                write(*,*) 'GROS PROBLEME j',i,j,enj(i,j), mask_lnd(i,j)
!                write(*,*) 'fluxFWFECB_fice', fluxFWFECB_fice(i,j)
!              endif
!            endif

!
!           enddo
!          enddo



!nb routage eau sur la grille ecbilt
          flux_gris2ecb(:,:) = 0.0
          do j =1, ubound(flux_gris2ecb, dim=2)
           do i =1, ubound(flux_gris2ecb, dim=1)
             if (mask_lnd(i,j).eq.0) then
           flux_gris2ecb(i,j) = flux_gris2ecb(i,j)                &
          !                                      + fluxFWFECB_fice(i,j)/darea(i)               &
          !                                       *darea(i)
                                                + fluxFWFECB_fice(i,j)
             else
           flux_gris2ecb(eni(i,j),enj(i,j)) = flux_gris2ecb(eni(i,j),enj(i,j))                &
           !                                     + fluxFWFECB_fice(i,j)/darea(i)               &
           !                                      *darea(eni(i,j))
                                                + fluxFWFECB_fice(i,j)
             endif
           enddo
          enddo

          do i =1, ubound(flux_gris2ecb, dim=1)
            flux_gris2ecb(i,:) = flux_gris2ecb(i,:)/darea(i)
          enddo

          !flux_gris2ecb(:,:)=0.0
          test_FWF_sum4=0.0
          do j =1, ubound(flux_gris2ecb, dim=2)
           do i =1, ubound(flux_gris2ecb, dim=1)
            test_FWF_sum4=test_FWF_sum4+flux_gris2ecb(i,j)*darea(i)
            !test_FWF_sum4=test_FWF_sum4+flux_gris2ecb(i,j)
           enddo
          enddo
          write(*,*) 'test_FWF_sum4 in Sv' , test_FWF_sum4*1e-6



!          do i =1, ubound(flux_gris2ecb, dim=1)
!           flux_gris2ecb(i,:) = fluxFWFECB_fice(i,:)/darea(i)
!          enddo


#endif

#if ( DOWNSTS == 1 )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Masks sub-grid in atmos model where the differences in altitude on land in the sub-grid is less than threshold
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
          call masks_subgrid_where_flat(topdifECB,sub_grid_notflat,10.0d0)
          where (fracto.gt.0.99) sub_grid_notflat = 0.
          call compute_weights_virt_levels(topoSG)

#if ( DOWNSCALING == 2 )

!afq -- 11/2022        success = interpolate_init(sgnx,sgny,iEcb,jEcb,latEcb,lonEcb,latSG,lonSG,3,nneigh,nexted,index_interpL2G, &
!afq -- 11/2022     &                             weights_interpL2G, sumweights_interpL2G, sub_grid_notflat)
          do nbd=1,sgd
              success = interpolate_init(sgnx(nbd),sgny(nbd),iEcb,jEcb,latEcb,lonEcb,                                     &
                   latSG(1:sgnx(nbd),1:sgny(nbd),nbd),lonSG(1:sgnx(nbd),1:sgny(nbd),nbd),3,nneigh,nexted,                 &
                   index_interpL2G(:,:,1:sgnx(nbd),1:sgny(nbd),nbd),                                                      &
                   weights_interpL2G(:,1:sgnx(nbd),1:sgny(nbd),nbd), sumweights_interpL2G(1:sgnx(nbd),1:sgny(nbd),nbd),   &
                   sub_grid_notflat)
          enddo

          call cluster_var_subgrid_in_ecbilt_dble(weights_low_2d,weights_low_sg)
          call cluster_var_subgrid_in_ecbilt_intg(index_low_2d,index_low_sg)
          call cluster_var_subgrid_in_ecbilt_dble(areaSG,area_sg)
          area_sg_onecb(:,:) = sum ( area_sg, dim=3 )
          call cluster_var_subgrid_in_ecbilt_dble(topoSG,topo_sg)
!#if (ISM >= 2)
          call cluster_var_subgrid_in_ecbilt_dble(epaisSG,epaisglace_sg)
! afq -- for high-resolution SMB computations, we also compute an elevation difference from the mean elevation
          call topoDiffsubgrid(topo_sg,difftopo_sg)
!#endif
#if ( SMB_TYP == 2 )
          call cluster_var_subgrid_in_ecbilt_dble(annbiasSG,annbias_sg)
#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&afq Create a not flat mask on grisli grid ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
          forall (n_dim=1:max_nb_points)
            subgrid_not_flat_sg(:,:,n_dim) = sub_grid_notflat(:,:)
          endforall

          call cast_ecb_subgrid_in_2d(subgrid_not_flat_sg,not_flat_subgrid)
#endif

          where (nbpointsSG(:,:).gt.0) 
              topoECB(:,:)   = topoECB_loc(:,:)
              masqueECB(:,:) = masqueECB_loc(:,:)
          endwhere
#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Check of what has been done through disk writing out ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if (1)
          CALL output_grd2D(masqueECB,iEcb,jEcb,'outputdata/coupler/checkmasqueECB',1,0.0)
          CALL output_grd2D(topoECB,iEcb,jEcb,  'outputdata/coupler/checktopoECB',1,0.0)
          CALL output_grd2D(topdifECB(:,:,nb_stat),iEcb,jEcb,'outputdata/coupler/topoDiff1',1,0.0)
          CALL output_grd2D(topdifECB(:,:,1),iEcb,jEcb,'outputdata/coupler/topoDiff2',1,0.0)
          CALL output_grd2D(sub_grid_notflat,iEcb,jEcb,'outputdata/coupler/checksub_grid_flat',1,0.0)
#endif

          return
        end subroutine subgrid_ecb_wrapper
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine compute_weights_virt_levels(topo_highres)

         double precision, dimension(:,:,:), intent(in) :: topo_highres

         integer :: i,j,n
        
         do n=1,sgd
            do j=1,sgny(n)
               do i=1,sgnx(n)
                  call weights_virt_levels(topo_highres(i,j,n),index_low_2d(i,j,n),weights_low_2d(i,j,n))
               enddo
            enddo
         enddo

        return
      end subroutine compute_weights_virt_levels


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: weights_virt_levels
!
!>     @brief find index in rmount which corresponds to the altitude just below a given alti_in
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


! dmr [TODO] Cette routine ne devrait être appelée qu'une fois par mise à jour de la topographie
! dmr        Les poids ne changent que si la topo change.
! dmr        Il suffit que les poids soient stockés dans une variable privée à ce module.
! dmr        Ensuite wights_virt_levels doit être appelé dans la routine iLOVECLIM de mis à jour de la topographie (ec_topo ...)
! dmr [TODO]

      subroutine weights_virt_levels(alti_in,ind_low,weight_low)

        use ecbilt_topography, only: nb_levls,rmount_virt

        implicit none

        ! a given altitude:
        double precision, intent(in)          :: alti_in
        ! index in rmount_virt just below alti_in:
        integer, intent(out)                  :: ind_low
        ! distance from remount_virt(ind_low) to alti_in:
        double precision, intent(out)         :: weight_low

        ! loop integer:
        integer :: nb_levl

        ind_low=1
        do nb_levl = 2, nb_levls
           if (rmount_virt(nb_levl).lt.alti_in) then
              ind_low=ind_low+1
           else
              exit
           endif
        end do

        ! linear interpolation between 2 levels in rmount_virt:
        weight_low = 1.
        ! afq if alti_in is lower than rmount_virt(1) weight_low has to be 1.
        if ((ind_low.lt.nb_levls).and.(alti_in.gt.rmount_virt(ind_low))) then
           weight_low = (rmount_virt(ind_low+1)-alti_in)/ &
                (rmount_virt(ind_low+1)-rmount_virt(ind_low))
        else if (ind_low.ge.nb_levls) then
           write (*,*) "Altitude in subgrid greater than max level in atmosphere..."
           ind_low = nb_levls-1
           weight_low = 0.d0
        endif
        
        if (weight_low.gt.1.) then
           write(*,*) "The vertical levels are flat! STOP"
           STOP
        end if
        
        ! at the given alti_in, any variable var can be defined as:
        !   var|alti_in = var|rmount_virt(ind_low)   x weight_low
        !               + var|rmount_virt(ind_low+1) x (1. - weight_low)

        return
      end subroutine weights_virt_levels


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: update_topo_masks
!
!>     @brief This subroutine allows updating the topography from GRISLI on ECBilt and the pertaining masks
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine update_topo_masks()

       use taillesGrilles, only: sgnxm, sgnym, sgd, iEcb, jEcb

!       Module contenant les definitions de la sous grille
       use input_subgrid2L, only: topoSG, masqueSGice, epaisSG
#if ( ISM >= 2 )
       use input_subgrid2L, only: topoGRIS, masqueGRISice=>masqueGRIS, epaisglaceGRISnord, nivo_mer, topoGRIS_sansnivomernord
       use global_constants_mod, only: str_len
#endif 


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer :: ix, jy
        integer :: checkGRISLImasksMKItxt_id
#if ( ISM >= 2 )
        integer :: nbd
        character(len=str_len) :: file_nm
#endif

#if ( ISM >= 2 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      1.0) Computations on the GRISLI or any subgrid grid: topo, masks
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      In GRISLI, the topography is relative to the current sea-level, even if it is lower / higher than present.
!      In ECBilt, the topography is related to constant defined sea-level, so we need to account for this (small) difference
!      nivo_mer is the sealevel variable in GRISLI
!
!      topoGRIS is thus the GRISLI variable corrected for sea-level change
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       topoGRIS(:,:,:) = 0.
       topoGRIS(1:sgnx(1),1:sgny(1),1) = max(topoGRIS_sansnivomernord(:,:) - dble(nivo_mer),0.0E0)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       1.2) Building of an icemask, based on GRISLI ice-sheet thickness
!            Beware that GRISLI assumes 1.0 meters of ice everywhere, hence the limit below
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       masqueGRISice(:,:,:) = 0.
       where(epaisglaceGRISnord(:,:) > 1.0)
         masqueGRISice(1:sgnx(1),1:sgny(1),1)=1.0
       !elsewhere
       !  masqueGRISice(:,:,1)=0.0
       end where

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Test output (optional)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#if (1)
!!        open(1111,file='outputdata/coupler/checkGRISLI_masks_LSM.txt',form='formatted')
!!        DO ix=1, GNX
!!           DO jy=1, GNY
!!             write(1111,'(2I5,12F15.3)') ix,jy, masqueGRISland(ix,jy)
!!           ENDDO
!!        ENDDO
!!        close(1111)

!!        open(1111,file='outputdata/coupler/checkGRISLI_masks_ISF.txt',form='formatted')
!!        DO ix=1, GNX
!!           DO jy=1, GNY
!!             write(1111,'(2I5,12F15.3)') ix,jy, masqueGRISiceshelves(ix,jy)
!!           ENDDO
!!        ENDDO
!!        close(1111)

       do nbd=1,sgd
          write (file_nm, "(A41,I1,A4)") "outputdata/coupler/checkGRISLI_masks_MKI_", nbd, ".txt"
          open(newunit=checkGRISLImasksMKItxt_id,file=file_nm,form='formatted')
          DO ix=1, sgnx(nbd)
             DO jy=1, sgny(nbd)
                write(checkGRISLImasksMKItxt_id,'(2I5,12F15.3)') ix,jy, masqueGRISice(ix,jy,1)
             ENDDO
          ENDDO
          close(checkGRISLImasksMKItxt_id)
       enddo
#endif
        epaisSG(:,:,:) = 0.
        epaisSG(1:sgnx(1),1:sgny(1),1)=epaisglaceGRISnord(:,:)
        topoSG(:,:,:)=topoGRIS(:,:,:)
#endif


        where(epaisSG > 1.0)
           masqueSGice=1.0
        elsewhere
           masqueSGice=0.0
        end where

       return
      end subroutine update_topo_masks
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: aggreg_one_field
!
!>     @brief This subroutine allows aggregating a field from a subgrid variable towards an ECBilt one
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine aggreg_one_field(subgrid_in,ecb_out,moyen_age)

        use input_subgrid2L, only: iEG, jEG

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  moyen_age Do we mean the field after aggregation  = 0 (as in temperature) or sum it = 1 (as in precipitation)
!>    @param[in]  gris_in The grisli variable to be aggregated
!>    @param[out] ecb_out The resulting aggregated field on the ECBilt grid.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer, intent(in) :: moyen_age
        double precision, dimension(:,:,:), intent(in) :: subgrid_in
        ! here ecb_out is expected to be in nlat, nlon, ECBilt std format
        double precision, dimension(:,:), intent(out):: ecb_out

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer :: subgrid_nx, subgrid_ny, subgrid_ng, ecb_nx, ecb_ny, nbd
        double precision, dimension(ubound(ecb_out,dim=1),ubound(ecb_out,dim=2))  :: ecb_outloc

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subgrid_nx = ubound(subgrid_in,1)
        subgrid_ny = ubound(subgrid_in,2)
        subgrid_ng = ubound(subgrid_in,3)
        ecb_nx  = ubound(ecb_out,1) ! nlat
        ecb_ny  = ubound(ecb_out,2) ! nlon

        ecb_outloc(:,:) = 0.

        do nbd=1,subgrid_ng
           call aggregSG2ECB(subgrid_in(:,:,nbd),subgrid_nx,subgrid_ny,ecb_outloc,ecb_nx,ecb_ny,                          &
                                            iEG(:,:,nbd),jEG(:,:,nbd),moyenne=moyen_age)
           ecb_out(:,:) = ecb_out(:,:) + ecb_outloc(:,:)  !the grids are supposed to be disconnected
        enddo

        return
      end subroutine aggreg_one_field

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Cette routine sert a aggreger les variables sous grille de type "topo"
!       sur la grille du modele ECBilt
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Masa Kageyama, Christophe Dumas
!      Date   : ? ? ?
!      Derniere modification : 21 juillet 2008, Didier M. Roche :
!                              Portage de ISM2CLIMBER dans LUDUS
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      SUBROUTINE aggregSG2ECB(champi,nx,ny,champc,ex,ey,iGRID,jGRID,   &
                               masque,moyenne)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  :
!        nx, ny   : taille du champ sous grille
!        ex, ey   : taille du champ ECBilt
!        champi   : champ sous grille
!        NB       : nombre de cases sous grille sur le point de grille ECBilt
!                   considere
!        i(j)GRID : indices de correspondances subgrid <-> ECBilt
!       Variables de sortie :
!        champc   : champ sous grille aggrege sur grille ECBilt
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nx, ny, ex, ey

      REAL, DIMENSION(nx,ny), INTENT(IN)  :: champi, iGRID, jGRID
      REAL, DIMENSION(ex,ey), INTENT(INOUT) :: champc

      REAL, DIMENSION(nx,ny), INTENT(IN), OPTIONAL :: masque
      INTEGER, INTENT(IN), OPTIONAL :: moyenne

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      INTEGER i,n, ii, jj
      INTEGER, DIMENSION(ex,ey) :: NB_lok

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Initialisation dans le cas ou l'on ne fait pas de min/max
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      champc = 0.
      NB_lok = 0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       OPTION 1 : on calcule avec un masque si celui-ci est fourni
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      IF (PRESENT(masque)) THEN

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       OPTION 1a : calcul avec une moyenne ... (type topo)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      IF (PRESENT(moyenne).AND.(moyenne.EQ.1)) THEN
      DO i=1,nx
        DO n=1,ny

            jj = MIN(ex,FLOOR(jGRID(i,n)))
            ii = FLOOR(iGRID(i,n))

            IF (masque(i,n).GT.0.0) THEN
              champc(jj,ii) = champc(jj,ii) + champi(i,n)
              NB_lok(jj,ii) = NB_lok(jj,ii) + 1
            ENDIF

         ENDDO
      ENDDO

      WHERE(NB_lok.NE.0) champc = champc / FLOAT(NB_lok)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       OPTION 1b : calcul sans moyenne ... (type precipitation)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      ELSE
      DO i=1,nx
        DO n=1,ny

            jj = MIN(ex,FLOOR(jGRID(i,n)))
            ii = FLOOR(iGRID(i,n))

            IF (masque(i,n).GT.0.0) THEN
              champc(jj,ii) = champc(jj,ii) + champi(i,n)
            ENDIF

         ENDDO
      ENDDO
      ENDIF

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       OPTION 2 : on calcule brutalement sans masque (aucun fourni)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      ELSE

      DO i=1,nx
        DO n=1,ny

            jj = MIN(ex,FLOOR(jGRID(i,n)))
            ii = FLOOR(iGRID(i,n))

              champc(jj,ii) = champc(jj,ii) + champi(i,n)
              NB_lok(jj,ii) = NB_lok(jj,ii) + 1

         ENDDO
      ENDDO

      IF (PRESENT(moyenne).AND.(moyenne.EQ.1)) THEN
        WHERE(NB_lok.NE.0) champc = champc / FLOAT(NB_lok)
      ENDIF

      ENDIF

      END SUBROUTINE aggregSG2ECB

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!      Auteur : Didier M. Roche
!      Date   : 15 juin 2009
!      Derniere modification : 07 juillet 2010, Didier M. Roche, 14 mars 2016, dmr
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine topoDiffECBsubgrid(champi,nx,ny,champc,ex,ey,iGRID,jGRID, champc2)

        integer, intent(in) :: ex, ey
        integer, dimension(:), intent(in) :: nx,ny

        real, dimension(:,:,:), intent(in)  :: champi            ! field on the subgrid
        real, dimension(:,:,:), intent(in)  :: iGRID, jGRID
        real, dimension(ex,ey), intent(inout) :: champc, champc2 ! fields on ECBilt

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer i,j,n,nd, ii, jj ! indicies
        double precision :: alt_olympus_mons = 21229.0d0

        nd = ubound(nx,dim=1)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        champc  = alt_olympus_mons*(-1.0)
        champc2 = alt_olympus_mons

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       On calcule les min/max de la topographie ECB / la sous grille
!-----|--1--------2---------3---------4---------5---------6---------7-|

        do n=1,nd
           do i=1,nx(n)
              do j=1,ny(n)

                 jj = MIN(ex,FLOOR(jGRID(i,j,n)))
                 ii = FLOOR(iGRID(i,j,n))

! dmr Changement du 07 Juillet 2010 : on sort non plus l'anomalie
! dmr  mais le champs absolu
                 champc(jj,ii)  = max(champi(i,j,n), champc(jj,ii))
                 champc2(jj,ii) = min(champi(i,j,n), champc2(jj,ii))
! Verif
              enddo
           enddo
        enddo

        return
      end subroutine topoDiffECBsubgrid

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
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

      subroutine masks_subgrid_where_flat(grid_minmax_topo, flat_subgrid_mask, threshold)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~        double precision, dimension(:,:),   intent(in) :: subgrid_landmask
       double precision, dimension(:,:,:), intent(in) :: grid_minmax_topo
       double precision,                   intent(in) :: threshold
       double precision, dimension(:,:),   intent(out):: flat_subgrid_mask

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~         integer :: sub_x, sub_y
        integer :: grd_x, grd_y, vert_coord
        double precision, dimension(:,:), allocatable :: local_topo_diff

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   1.0) Get the indicies of the arrays to be dealt with
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~         sub_x = ubound(subgrid_landmask,1)
!~         sub_y = ubound(subgrid_landmask,2)

        grd_x = ubound(flat_subgrid_mask,1)
        grd_y = ubound(flat_subgrid_mask,2)

                ! here is the implicit assumption that grid_minmax_topo and flat_subgrid_mask share the first two axes
        vert_coord = ubound(grid_minmax_topo,3)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   2.0) Create the difference in topography
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         allocate(local_topo_diff(grd_x,grd_y))

         local_topo_diff(:,:) = grid_minmax_topo(:,:,vert_coord) - grid_minmax_topo(:,:,1)

         where(local_topo_diff > threshold)
            flat_subgrid_mask = 1.0
         elsewhere
            flat_subgrid_mask = 0.0
         endwhere

         return
      end subroutine masks_subgrid_where_flat

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( F_PALAEO == 3 )

      subroutine update_forcedISM

        use input_subgrid2L, only: indx_fISM,topoGRIS,masqueGRISice=>masqueGRIS,forcedISM,forcedMSK

        use taillesGrilles, only: sgnx, sgny
        use routageEAU_mod, only: exutiSG, exutjSG, eniSG, enjSG, nbexuts

        integer :: i,j

        topoGRIS(:,:) = max(forcedISM(:,:,indx_fISM),0.0E0)
        masqueGRISice(:,:) = 1.-forcedMSK(:,:,indx_fISM)

! dmr&nb --- MàJ du routage sur nouvelle TOPO et nouveau masque


        !CALL routageEAU(GNX,GNY,topoGRIS,masqueGRISice,nbexuts,exutiSG,exutjSG,eniSG,enjSG)
        !nb Beware the mask should be integers for routageEAU!
        CALL routageEAU(sgnx,sgny,topoGRIS,int(masqueGRISice),nbexuts,exutiSG,exutjSG,eniSG,enjSG,1)

!        do i=1,ubound(eniSG,dim=1)
!        do j=1, ubound(eniSG,dim=2)
!         if (eniSG(i,j) .ne.0) then
!            write(*,*) 'routageEAU exutoires i', i,j, eniSG(i,j)
!         endif
!         if (enjSG(i,j) .ne.0) then
!            write(*,*) 'routageEAU exutoires j', i,j, enjSG(i,j)
!         endif
!        enddo
!        enddo


        return
      end subroutine update_forcedISM

#if ( F_PALAEO_FWF == 1 )

       subroutine update_forced_FWF

          use input_subgrid2L,      only: indx_fISM, forcedTHI, fluxFWF_SG, update_time_fism, fluxFWF_SG_route
          use routageEAU_mod,       only: eniSG, enjSG
          use global_constants_mod, only: days_year360d, solar_day
          use taillesGrilles,       only: sgnx, sgny


          integer :: indx_plusun
          integer :: i,j
          !real :: test_FWF_sum
          real :: test_FWF_sum2

          indx_plusun = indx_fISM+1
          if (indx_plusun.GT.UBOUND(forcedTHI,DIM=3)) then
            indx_plusun = indx_fISM
          endif

          !dmr&nb : calculate the freshwater flux arising from the
          !change in ice-sheet thickness between two timestep. The
          !resulting flux is in m^3.s^-1. The fixed multiplicand
          !40000*40000 arises from the area of the GRISLI grid in m^2
          !fluxFWF_SG(:,:) =  (forcedTHI(:,:,indx_plusun) - forcedTHI(:,:,indx_fISM))*(40000*40000.)                              &
!nb to have a positive fresh water flux to the ocean when the ice sheet
!volume is decreasing
          fluxFWF_SG(:,:) =  (forcedTHI(:,:,indx_fISM) - forcedTHI(:,:,indx_plusun))*(40000*40000.)                              &
!          fluxFWF_SG(:,:) =  (forcedTHI(:,:,indx_fISM) - forcedTHI(:,:,indx_plusun))                                             &
                           / (update_time_fism*days_year360d*solar_day)

          !write(*,*) 'fluxFWF_SG dans transfer_subgrid_ecb_mod ', fluxFWF_SG

          !dmr    : route the water obtained on the GRISLI grid with
          !the river routing network obtained in the routageEAU routine
          !it results in a flux routed to the border of the ice-sheet

          !check: sum total
          !test_FWF_sum=0.0
          test_FWF_sum2=0.0

!          !nb mettre a 0 d abord ?
!          fluxFWF_SG_route(:,:)=0.

!nb at the moment use the routage done at restart
!          do i=1,GNX
!          do j=1,GNY
!           fluxFWF_SG_route(eniSG(i,j), enjSG(i,j)) = fluxFWF_SG_route(eniSG(i,j), enjSG(i,j)) + fluxFWF_SG(i,j) ! in m^3.s^-1
!!           fluxFWF_SG_route(eniSG(i,j), enjSG(i,j)) = fluxFWF_SG_route(eniSG(i,j), enjSG(i,j)) + fluxFWF_SG(i,j)/(40000.*40000.) ! in m.s^-1
!!           if (fluxFWF_SG_route(eniSG(i,j), enjSG(i,j)) .ne. 0) then
!!              write(*,*) 'fluxFWF_SG_route dans transfer_subgrid_ecb_mod ' ,eniSG(i,j),enjSG(i,j),fluxFWF_SG_route(eniSG(i,j), enjSG(i,j))
!!           endif
!           test_FWF_sum=test_FWF_sum+fluxFWF_SG(i,j)
!          enddo
!          enddo
   
          fluxFWF_SG_route(:,:)=fluxFWF_SG(:,:)


          do i=1,sgnx
          do j=1,sgny
          test_FWF_sum2=test_FWF_sum2+fluxFWF_SG_route(i,j)
          enddo
          enddo

          write(*,*) ' '
          !write(*,*) ' TEST FWF_sum in Sv', test_FWF_sum*1e-6
          write(*,*) ' TEST FWF_sum2 in Sv', test_FWF_sum2*1e-6

!          do i=1,ubound(fluxFWF_SG_route,dim=1)
!          do j=1,ubound(fluxFWF_SG_route,dim=2)
!            if (fluxFWF_SG_route(i,j) .ne. 0) then
!              write(*,*) 'fluxFWF_SG_route 1 ' , i,j, fluxFWF_SG_route(i,j)
!           endif
!          enddo
!          enddo

          !write(*,*) 'fluxFWF_SG_route dans transfer_subgrid_ecb_mod ' ,  fluxFWF_SG_route

       end subroutine
#endif


#endif

#if ( CONSEAU == 1 )

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine sert a aggreger les variables ISM de type calving
!       sur la grille du modele CLIO (base sur aggregISM2ECB)
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Masa Kagcyama, Christophe Dumas
!      Date   : ? ? ?
!      Derniere modification : 21 juillet 2008, Didier M. Roche :
!                              Portage de ISM2CLIMBER dans LUDUS
!-----|--1--------2---------3---------4---------5---------6---------7-|

      SUBROUTINE aggregISM2CLIO(champi,nx,ny,champc,cx,cy,iGRID,jGRID,  &
                               masque_oc)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  :
!        nx, ny   : taille du champ ISM
!        cx, cy   : taille du champ CLIO
!        champi   : champ ISM
!        NB       : nombre de cases ISM sur le point de grille ECBilt
!                   considere
!        i(j)GRID : indices de correspondances ISM <-> CLIO
!       Variables de sortie :
!        champc   : champ ISM aggrege sur grille CLIO
!-----|--1--------2---------3---------4---------5---------6---------7-|

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nx, ny, cx, cy

      REAL, DIMENSION(nx,ny), INTENT(IN)  :: champi
      INTEGER, DIMENSION(nx,ny), INTENT(IN)  :: iGRID, jGRID
      REAL, DIMENSION(cx,cy), INTENT(OUT) :: champc

      INTEGER, DIMENSION(cx,cy), INTENT(IN), OPTIONAL :: masque_oc

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

      INTEGER :: i,n,j,k,l,indxi,indxj, radii
      INTEGER, PARAMETER :: radius = 1
      INTEGER iii,jjj
      LOGICAL :: FIXED

      champc(:,:) = 0.

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       OPTION 1 : on calcule avec un masque si celui-ci est fourni
!         ATTENTION : ici le masque est celui de CLIO, pour assurer que
!         la case que l'on considere est bien oceanique
!-----|--1--------2---------3---------4---------5---------6---------7-|

      IF (PRESENT(masque_oc)) THEN

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       OPTION 1b : calcul sans moyenne ... (type precipitation)
!-----|--1--------2---------3---------4---------5---------6---------7-|

      DO i=1,nx
        DO n=1,ny

        IF (masque_oc(iGRID(i,n),jGRID(i,n)).GT.0) THEN
         champc(iGRID(i,n),jGRID(i,n)) =  champc(iGRID(i,n),jGRID(i,n)) &
              + champi(i,n)
        ELSEIF (masque_oc(iGRID(i,n),jGRID(i,n)).LE.0.AND.champi(i,n).GT.0.0) then
!cdmr&mab --- To take into account when there is a mask discrepancy ...
          FIXED=.FALSE.

          radii=1

          DO WHILE ((.NOT.FIXED).AND.(radii.LE.radius))

          DO k=radii*(-1),radii
            DO l=radii*(-1),radii

              indxi =  iGRID(i,n)+k
              indxj =  jGRID(i,n)+l

              IF (indxi.GT.cx) indxi = indxi-cx
              IF (indxj.GT.cy) indxj = indxj-cy

              IF ((masque_oc(indxi,indxj).GT.0).AND.(.NOT.FIXED)) THEN
                champc(indxi,indxj) = champc(indxi,indxj) + champi(i,n)
                FIXED=.TRUE.
              ENDIF

            ENDDO
          ENDDO

          radii = radii+1

          ENDDO ! while not FIXED

            IF (.NOT.FIXED) THEN

             write(*,*) 'ALERT!! ice is lost in the coupling betw. CLIO &
                and GRISLI!! (aggregISM2CLIO)', i,n

            ENDIF ! if not fixed
        ENDIF  ! (no) ocean point but calv flux still
        ENDDO
      ENDDO

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       OPTION 2 : on calcule brutalement sans masque (aucun fourni)
!-----|--1--------2---------3---------4---------5---------6---------7-|
      ELSE

      DO i=1,nx
        DO n=1,ny

         champc(iGRID(i,n),jGRID(i,n)) =  champc(iGRID(i,n),jGRID(i,n)) &
              + champi(i,n)

         ENDDO
      ENDDO

      ENDIF

      END SUBROUTINE aggregISM2CLIO

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#endif

#if ( DOWNSCALING == 2 )
subroutine topoDiffsubgrid(topoin,difftopo)

     use input_subgrid2L,          only: nbpointssg

     double precision, dimension(:,:,:), intent(in)    :: topoin
     double precision, dimension(:,:,:), intent(inout) :: difftopo

     integer :: i,j

     difftopo(:,:,:) = 0.
     do j=1,ubound(topoin,dim=2)
        do i=1,ubound(topoin,dim=1)
           if(nbpointssg(i,j).gt.0) then
              difftopo(i,j,1:nbpointssg(i,j)) = topoin(i,j,1:nbpointssg(i,j))-sum(topoin(i,j,1:nbpointssg(i,j)))/nbpointssg(i,j)
           else
              difftopo(i,j,:) = 0.
           endif
        enddo
     enddo

     return  
end subroutine topoDiffsubgrid
#endif



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      end module transfer_subgrid_ecb_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
