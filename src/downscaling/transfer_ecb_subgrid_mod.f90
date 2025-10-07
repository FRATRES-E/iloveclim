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
!      MODULE: transfer_ecb_grisli_mod
!
!>     @author  Aurelien Quiquet (afq)
!>     @author  Didier M. Roche  (dmr)
!
!>     @brief This module transf_ecb_grisli_mod is handling the transfer of climate variables from the ECBilt grid
!               to the GRISLI grid. It is grid-specific to the particular grids given.
!
!>     @date Creation date: January, 28th, 2016
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       module transfer_ecb_subgrid_mod

        ! -- afq, variables needed by the subgrid, on the ECBilt grid and on vertical levels (i,j,nlevs):
        use sgout_mass_balance_mod, only: temp_profile_smb_ann, totp_profile_smb_ann
#if ( ISM >= 2 || SMB_TYP >= 1 )
        use sgout_mass_balance_mod, only: SMB_iLCM_profile_ann
#endif
        ! -- afq, variables needed by the subgrid, on the sub-grid (i,j,n_point):
#if ( DOWNSCALING == 2 )
        use sgout_mass_balance_mod, only: temp_sg_smb_ann, totp_sg_smb_ann, relhum_sg_smb_ann
#if ( ISM >= 2 || SMB_TYP >= 1 )
        use sgout_mass_balance_mod, only: SMB_iLCM_sg_ann
#if ( REFREEZING == 1 )
        use sgout_mass_balance_mod, only: rain_sg_smb_ann, snow_sg_smb_ann, melt_sg_smb_ann
#endif 
#endif   
#endif
        ! -- afq, variables needed by GRISLI, on the GRISLI grid and 2D:
        use input_subgrid2L, only: tfyear, pfyear, relhumyear
#if ( ISM >= 2 || SMB_TYP >= 1)
        use input_subgrid2L, only: SMB_iLCM
#endif
#if ( ISM >= 2 )
        use input_subgrid2L, only: tfyearnord, pfyearnord, SMB_iLCMnord
#endif
        ! -- afq, function that does the horizontal interpolation:
! afq -- outdated 11/2022        use interpol_ecb_subgrid_mod, only: interpol_one_field
        use interpolate_mod,    only: interpolate
        ! --- afq&dmr, size of the subgrid grids
        use taillesGrilles, only: sgnx, sgny, sgnxm, sgnym, sgd

        use transfer_subgrid_ecb_mod, only: weights_low_2d, index_low_2d

        implicit none


        public  :: transfer_wrapper_subgrid
! afq -- notes on flags: currently (4/7/2016) we need DOWNSTS == 1 for ISM == 3 | ISM == 2
! but on the long go, we should be able to do DOWNSTS == 0 with an ice sheet...
        private  :: write_nc2d_subgrid_init
        public   :: init_subgrid_writing
#if ( DOWNSCALING == 2 )
        private  :: transfer_wrapper_sg
#endif

#if ( 0 )
        ! -- afq, internal routine that horizontally and vertically interpolate from the
        !         ECBilt to the subgrid
        private :: transfer_ecb_subgrid
#endif

        integer, private :: pass_nc
        integer, private :: iyear_nc
        
      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [init_subgrid_writing]
!
!>     @brief 
!
!      DESCRIPTION:
!
!>     *** 
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      function init_subgrid_writing() result(returnValue)
!-----------------------------------------------------------------------
!     *** this data files netCDF writing for subgrid 
!-----------------------------------------------------------------------

       use global_constants_mod, only: str_len

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  ist Blablabla ...
!>    @param[in]  jst Blablabla ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical                :: returnValue

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      character(len=str_len) :: file_nm
      integer nbd
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

# if ( 0 )
       file_nm = "outputdata/downscaling/sgout_profil.nc"
       call write_nc2d_subgrid_init(file_nm)
#endif
#if ( DOWNSCALING == 2 )
       do nbd=1,sgd
          write (file_nm, "(A37,I1,A3)") "outputdata/downscaling/sgout_subgrid_", nbd, ".nc"
          !file_nm = "outputdata/downscaling/sgout_subgrid_"//trim(nbd)//".nc"
          call write_nc2d_subgrid_init(file_nm,nbd)
       enddo
#endif

      returnValue = .TRUE.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      end function init_subgrid_writing

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine transfer_wrapper_subgrid()

        ! -- afq, the topo in the subgrid
        use input_subgrid2L, only: topoSG

#if ( DOWNSCALING == 2 )
        use input_subgrid2L, only: max_nb_points,nbpointsSG,sub_grid_notflat,not_flat_subgrid
#if ( ISM >= 2 || SMB_TYP >= 1 )
        use input_subgrid2L, only: topo_sg,nivo_mer
#endif
        use ecbilt_topography, only: rmount_virt
        use interpolate_mod,   only: interpolate
        use input_subgrid2L,    only: nneigh, index_interpL2G, weights_interpL2G, sumweights_interpL2G
#endif

        integer,parameter :: jEcb  = ubound(temp_profile_smb_ann,dim=1)
        integer,parameter :: iEcb  = ubound(temp_profile_smb_ann,dim=2)
        integer,parameter :: nvert = ubound(temp_profile_smb_ann,dim=3)

        character(len=256) :: file_nm
#if ( DOWNSCALING == 2 )
        double precision, dimension(sgnxm,sgnym,sgd) :: tann_down, prec_down, relhum_down
#if ( ISM >= 2 || SMB_TYP >= 1 )
        double precision, dimension(sgnxm,sgnym,sgd) :: smb_down
        double precision, dimension(iEcb,jEcb) :: temp_sea,prec_sea,smb_sea
        double precision, dimension(sgnxm,sgnym,sgd) :: temp_sea_interp,prec_sea_interp,smb_sea_interp
        integer          :: nsea          ! number of oceanic sub-grid points
#if ( REFREEZING == 1 )
        double precision, dimension(iEcb,jEcb,max_nb_points) :: refreez_ann
#endif
#endif
        logical          :: success
        integer          :: n_point
#endif
        
        double precision, dimension(sgnxm,sgnym,sgd) :: tann_loc, prec_loc, smb_loc, topo_loc
        
        integer :: i,j,nbd
        
        topo_loc(:,:,:) = topoSG(:,:,:)

#if ( 0 )
        write(*,*) "Transfer of temperature"

        call transfer_ecb_subgrid(temp_profile_smb_ann,tann_loc,         &
             &       jEcb,iEcb,nvert,sgnx,sgny,0)

        write(*,*) "....................... temperature done."

        write(*,*) "Transfer of precip"

        call transfer_ecb_subgrid(totp_profile_smb_ann,prec_loc,         &
             &       jEcb,iEcb,nvert,sgnx,sgny,1)

        write(*,*) ".................. precip done."

        write(*,*) "Transfer of SMB"
#if ( ISM >= 2 || SMB_TYP >= 1 )
        call transfer_ecb_subgrid(SMB_iLCM_profile_ann,smb_loc,          &
             &       jEcb,iEcb,nvert,sgnx,sgny,0)

        write(*,*) "............... SMB done."
#endif

        tfyear(:,:)=real(tann_loc(:,:))
        pfyear(:,:)=real(prec_loc(:,:))
#if ( ISM >= 2 || SMB_TYP >= 1 )
        pfyear(:,:)=real(prec_loc(:,:))*1000./910. ! in ice equivalent
        SMB_iLCM(:,:)=real(smb_loc(:,:))
#endif
#endif /* on if(0), we do no longer use the profile variables to force the ISM */

#if ( DOWNSCALING == 2 )

#if ( ISM >= 2 || SMB_TYP >= 1 )

#if ( REFREEZING == 1 )
! afq -- 28/6/17: we modify the SMB to account for the refreezing
        call compute_refreezing(temp_sg_smb_ann,rain_sg_smb_ann,snow_sg_smb_ann,melt_sg_smb_ann,refreez_ann)
        SMB_iLCM_sg_ann(:,:,:) = SMB_iLCM_sg_ann(:,:,:) + refreez_ann(:,:,:) !melt_sg_smb_ann(:,:,:) !SMB_iLCM_sg_ann(:,:,:) !refreez_ann(:,:,:)
#endif

! afq -- this is the native sub-grid field:
! afq -- 15/3/18: new version ---------------------------
! 
! The downsc. subgrid precip are now smooth (almost) so we used them to force GRISLI
! However, we do not have subgrid info where sub_grid_notflat = 0...
! (For ice sheet coupling the sub_grid_notflat has to be 0 only for sea points)
! If we fill the gap with the _profile outputs strong discontinuities appear
! For this reason, we generate a new _profile that is the precip at sea level

        temp_sea(:,:)=0.d0
        prec_sea(:,:)=0.d0
        smb_sea(:,:) =0.d0

        do i= 1, iEcb
           do j= 1, jEcb

                if (sub_grid_notflat(i,j) .gt. 0.d0) then
                
                 nsea=0
                 do n_point=1,nbpointsSG(i,j)
                    if (topo_sg(i,j,n_point).lt.nivo_mer+0.1) then
                       nsea=nsea+1
                       temp_sea(i,j)=temp_sea(i,j)+temp_sg_smb_ann(i,j,n_point)
                       prec_sea(i,j)=prec_sea(i,j)+totp_sg_smb_ann(i,j,n_point)
                       smb_sea(i,j) =smb_sea(i,j) +SMB_iLCM_sg_ann(i,j,n_point)
                    endif
                 enddo
                 if (nsea.eq.0) then
                   temp_sea(i,j)= maxval(temp_sg_smb_ann(i,j,1:nbpointsSG(i,j)))
                   prec_sea(i,j)= minval(totp_sg_smb_ann(i,j,1:nbpointsSG(i,j)))
                   smb_sea(i,j) = minval(SMB_iLCM_sg_ann(i,j,1:nbpointsSG(i,j)))
                   ! in fact in this case we have a continental point
                 else
                   temp_sea(i,j)= temp_sea(i,j)/nsea
                   prec_sea(i,j)= prec_sea(i,j)/nsea
                   smb_sea(i,j) = smb_sea(i,j) /nsea
                 endif
                 
                else
                   temp_sea(i,j)= temp_profile_smb_ann(j,i,1)
                   prec_sea(i,j)= totp_profile_smb_ann(j,i,1)
                   smb_sea(i,j) = SMB_iLCM_profile_ann(j,i,1)
                endif
                
           enddo
        enddo

        do nbd=1,sgd

           success = interpolate(index_interpL2G(:,:,1:sgnx(nbd),1:sgny(nbd),nbd),           &
                weights_interpL2G(:,1:sgnx(nbd),1:sgny(nbd),nbd),                          &
                sumweights_interpL2G(1:sgnx(nbd),1:sgny(nbd),nbd),                       &
                transpose(temp_sea),temp_sea_interp(1:sgnx(nbd),1:sgny(nbd),nbd),            &
                sgnx(nbd),sgny(nbd), 3, nneigh, jEcb, iEcb)

           success = interpolate(index_interpL2G(:,:,1:sgnx(nbd),1:sgny(nbd),nbd),           &
                weights_interpL2G(:,1:sgnx(nbd),1:sgny(nbd),nbd),                          &
                sumweights_interpL2G(1:sgnx(nbd),1:sgny(nbd),nbd),                       &
                transpose(prec_sea),prec_sea_interp(1:sgnx(nbd),1:sgny(nbd),nbd),            &
                sgnx(nbd),sgny(nbd), 3, nneigh, jEcb, iEcb)
           
           success = interpolate(index_interpL2G(:,:,1:sgnx(nbd),1:sgny(nbd),nbd),           &
                weights_interpL2G(:,1:sgnx(nbd),1:sgny(nbd),nbd),                          &
                sumweights_interpL2G(1:sgnx(nbd),1:sgny(nbd),nbd),                       &
                transpose(smb_sea),smb_sea_interp(1:sgnx(nbd),1:sgny(nbd),nbd),              &
                sgnx(nbd),sgny(nbd), 3, nneigh, jEcb, iEcb)
           
        enddo

        call transfer_wrapper_sg(tann_down,prec_down,smb_down,relhum_down)
        
        where (not_flat_subgrid(:,:,:) .le. 0.0d0 + 0.0001d0)
                        tann_down(:,:,:)=temp_sea_interp(:,:,:)
                        prec_down(:,:,:)=prec_sea_interp(:,:,:)
                        smb_down(:,:,:) = smb_sea_interp(:,:,:)
        endwhere

#else
        call transfer_wrapper_sg(tann_down,prec_down,relhum_down)
#endif /* on ISM >= 2 */

        tfyear = real(tann_down)
        pfyear = real(prec_down)
        relhumyear = real(relhum_down)
#if ( ISM >= 2 || SMB_TYP >= 1 ) 
        pfyear = real(prec_down) *1000./910. ! in ice equivalent
        SMB_iLCM = real(smb_down)*1000./910. ! in ice equivalent
#if ( ISM >= 2 )
        tfyearnord(:,:) = tfyear(1:sgnx(1),1:sgny(1),1)
        pfyearnord(:,:) = pfyear(1:sgnx(1),1:sgny(1),1)
        SMB_iLCMnord(:,:) = SMB_iLCM(1:sgnx(1),1:sgny(1),1)
#endif
#endif

#endif /* on DOWNSCALING == 2 */
        
        write(*,*) "AUREL TIME:", iyear_nc, mod(iyear_nc,1)
        
        if (mod(iyear_nc,1).eq.0) then

#if ( 0 )
           file_nm = "outputdata/downscaling/sgout_profil.nc"
#if ( ISM >= 2 || SMB_TYP >= 1 ) 
           call write_nc2d_subgrid(tann_loc,prec_loc,smb_loc,file_nm)
#else
           call write_nc2d_subgrid(tann_loc,prec_loc,prec_loc,file_nm)
#endif
#endif
           
#if ( DOWNSCALING == 2 )

           do nbd=1,sgd

              write (file_nm, "(A37,I1,A3)") "outputdata/downscaling/sgout_subgrid_", nbd, ".nc"
              !file_nm = "outputdata/downscaling/sgout_subgrid_"//trim(nbd)//".nc"
#if ( ISM >= 2 || SMB_TYP >= 1 ) 
              call write_nc2d_subgrid(tann_down(1:sgnx(nbd),1:sgny(nbd),nbd),                  &
                   prec_down(1:sgnx(nbd),1:sgny(nbd),nbd),                                     &
                   smb_down(1:sgnx(nbd),1:sgny(nbd),nbd),                                      &
                   relhum_down(1:sgnx(nbd),1:sgny(nbd),nbd),file_nm)
#else
              call write_nc2d_subgrid(tann_down(1:sgnx(nbd),1:sgny(nbd),nbd),                  &
                   prec_down(1:sgnx(nbd),1:sgny(nbd),nbd),                                     &
                   relhum_down(1:sgnx(nbd),1:sgny(nbd),nbd),file_nm)
#endif
           enddo
#endif
           pass_nc = pass_nc + 1
           
        end if

        iyear_nc = iyear_nc + 1
        

        return
      end subroutine transfer_wrapper_subgrid


#if ( DOWNSTS == 1 )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( DOWNSCALING == 2 )
#if ( ISM >= 2 || SMB_TYP >= 1 )
      subroutine transfer_wrapper_sg(tann_tmp,prec_tmp,smb_tmp,relhum_tmp)
#else
      subroutine transfer_wrapper_sg(tann_tmp,prec_tmp,relhum_tmp)
#endif

        use sgout_mass_balance_mod, only: temp_sg_smb_ann, totp_sg_smb_ann, relhum_sg_smb_ann
#if ( DOWN_T2M == 1 )
        use sgout_mass_balance_mod, only: t2m_sg_smb_ann
#endif
#if ( ISM >= 2 || SMB_TYP >= 1 )
        use sgout_mass_balance_mod, only: SMB_iLCM_sg_ann
#endif
        use input_subgrid2L, only: cast_ecb_subgrid_in_2d

        double precision, dimension(sgnxm,sgnym,sgd), intent(out) :: tann_tmp
        double precision, dimension(sgnxm,sgnym,sgd), intent(out) :: prec_tmp
        double precision, dimension(sgnxm,sgnym,sgd), intent(out) :: relhum_tmp        
#if ( ISM >= 2 || SMB_TYP >= 1 )
        double precision, dimension(sgnxm,sgnym,sgd), intent(out) :: smb_tmp
#endif

#if ( DOWN_T2M == 1 )
        call cast_ecb_subgrid_in_2d(t2m_sg_smb_ann,tann_tmp)
#else
        call cast_ecb_subgrid_in_2d(temp_sg_smb_ann,tann_tmp)
#endif
        call cast_ecb_subgrid_in_2d(totp_sg_smb_ann,prec_tmp)
        call cast_ecb_subgrid_in_2d(relhum_sg_smb_ann,relhum_tmp)        
#if ( ISM >= 2 || SMB_TYP >= 1 )
        call cast_ecb_subgrid_in_2d(SMB_iLCM_sg_ann,smb_tmp)
#endif

        return
      end subroutine transfer_wrapper_sg

#endif

      subroutine write_nc2d_subgrid_init(filename,ngrid)

        use ncio

        character(len=256), intent(in) :: filename
        integer, intent(in)            :: ngrid
        
        double precision, allocatable :: x(:), y(:)

        integer :: ndims, i
        character(len=32), allocatable :: dimnames(:)
        integer, allocatable :: dimlens(:)

        ! Define array sizes and allocate arrays

        allocate(x(sgnx(ngrid)))
        allocate(y(sgny(ngrid)))

        pass_nc = 1
        iyear_nc = 0

! afq -- the extent in X and Y is hard-coded and specific for GRISLI for now
!        this could be changed reading xmin/xmax when reading the subgrid
        x(1) = -5600000.d0
        do i=2,sgnx(ngrid)
          x(i) = x(i-1) + 40000.d0
        enddo

        y(1) = -5200000.d0
        do i=2,sgny(ngrid)
          y(i) = y(i-1) + 40000.d0
        enddo
       
        write(*,*)
        write(*,*) "====== WRITING ======"
        write(*,*)

        ! Create the netcdf file, write global attributes
        call nc_create(filename,overwrite=.TRUE.,netcdf4=.TRUE.)
        call nc_write_attr(filename,"Title","Subgrid downscaling writing ...")
        call nc_write_attr(filename,"Institution", &
                       "Laboratoire de Sciences du Climat et de l'Environnement, ACCLIMATE project")

        call nc_write_dim(filename,"x",x=x,units="m")
        call nc_write_dim(filename,"y",x=y,units="kilometers")
        call nc_write_dim(filename,"time",x=1.0, &
                      units="years",calendar="360_day", unlimited=.TRUE.)


        return
      end subroutine write_nc2d_subgrid_init

#if ( ISM >= 2 || SMB_TYP >= 1 )
      subroutine write_nc2d_subgrid(tann_tmp,prec_tmp,smb_tmp,relhum_tmp,filename)
#else
      subroutine write_nc2d_subgrid(tann_tmp,prec_tmp,relhum_tmp,filename)
#endif    

        use ncio

        character(len=256), intent(in) :: filename

        double precision, dimension(:,:), intent(in) :: tann_tmp
        double precision, dimension(:,:), intent(in) :: prec_tmp
        double precision, dimension(:,:), intent(in) :: relhum_tmp
#if ( ISM >= 2 || SMB_TYP >= 1 )
        double precision, dimension(:,:), intent(in) :: smb_tmp
#endif  

        double precision,   allocatable :: d2D(:,:)

        integer nxloc,nyloc

        nxloc=ubound(tann_tmp,dim=1)
        nyloc=ubound(tann_tmp,dim=2)

        ! Define array sizes and allocate arrays

        allocate(d2D(nxloc,nyloc))

        write(*,*)
        write(*,*) "====== WRITING ======"
        write(*,*)
        write(*,*) "aurel pass_nc", pass_nc
        
        call nc_write(filename,"time",pass_nc,dim1="time",start=[pass_nc],count=[1])

        d2D(:,:) = tann_tmp(:,:)
        call nc_write(filename,"Tann",d2D(:,:),dim1="x",dim2="y",       &
             &       dim3="time",start=[1,1,pass_nc],count=[nxloc,nyloc,1])
        d2D(:,:) = prec_tmp(:,:)
        call nc_write(filename,"Acc",d2D(:,:),dim1="x",dim2="y",       &
             &       dim3="time",start=[1,1,pass_nc],count=[nxloc,nyloc,1])
#if ( ISM >= 2 || SMB_TYP >= 1 )
        d2D(:,:) = smb_tmp(:,:)
        call nc_write(filename,"Bm",d2D(:,:),dim1="x",dim2="y",       &
             &       dim3="time",start=[1,1,pass_nc],count=[nxloc,nyloc,1])
#endif    
        d2D(:,:) = relhum_tmp(:,:)
        call nc_write(filename,"relhum",d2D(:,:),dim1="x",dim2="y",       &
             &       dim3="time",start=[1,1,pass_nc],count=[nxloc,nyloc,1])

        deallocate(d2D)


        return
      end subroutine write_nc2d_subgrid

#endif

#if ( 0 )
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE:
!
!>     @brief This subroutine
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        
      subroutine transfer_ecb_subgrid(varin,varout,nx_ecb,ny_ecb,nz_ecb,nx_ism,ny_ism,positive)

        use input_subgrid2L, only: nneigh, index_interpL2G, weights_interpL2G, sumweights_interpL2G
        use comatm,          only: nlon,nlat

        implicit none

        integer, intent(in) :: nx_ecb,ny_ecb,nz_ecb,nx_ism,ny_ism
        integer, intent(in) :: positive
        double precision, dimension (nx_ecb,ny_ecb,nz_ecb), intent(in) :: varin
        double precision, dimension (nx_ism,ny_ism),        intent(inout) :: varout

        ! -- afq, horizontal interpolation of all vertical levels:
        double precision, dimension (nx_ism,ny_ism,nz_ecb) :: varloc

        ! -- afq, for vertical interpolation:
        integer          :: i,j,nb_lev
        logical          :: success

        do nb_lev = 1, nz_ecb
! afq -- outdated 11/2022           success = interpol_one_field(varin(:,:,nb_lev),varloc(:,:,nb_lev),positive)
           success = interpolate (index_interpL2G,weights_interpL2G, sumweights_interpL2G,             &
                                  varin(:,:,nb_lev),varloc(:,:,nb_lev),                                &
                                  nx_ism, ny_ism, 3, nneigh, nlon, nlat)
        enddo

        do j=1,ny_ism
           do i=1,nx_ism
              varout(i,j) = varloc(i,j,index_low_2d(i,j))*weights_low_2d(i,j)                   &
                          + varloc(i,j,index_low_2d(i,j)+1)*(1.d0 - weights_low_2d(i,j))
           enddo
        enddo

        return
       end subroutine transfer_ecb_subgrid

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#endif /* on if ( 0 ) */
       

#if ( SMB_TYP >= 1 && DOWNSCALING == 2 && REFREEZING == 1 )

      subroutine compute_refreezing(tannin, rainin, snowin, meltin, refreez)

        use global_constants_mod, only: latheat_fus_wat               ! latent heat of fusion 

        double precision,dimension(:,:,:),intent(in)  :: tannin
        double precision,dimension(:,:,:),intent(in)  :: rainin
        double precision,dimension(:,:,:),intent(in)  :: snowin
        double precision,dimension(:,:,:),intent(in)  :: meltin
        double precision,dimension(:,:,:),intent(out) :: refreez
        
        double precision, parameter                   :: tal = 1.     ! thermodynamic active layer
        double precision, parameter                   :: capf = 2.2   ! capillarity factor
        
        ! specific heat capacity of ice:
        double precision, dimension(ubound(tannin,dim=1),ubound(tannin,dim=2),ubound(tannin,dim=3)) :: ishc
        
        ishc(:,:,:) = 152.5+7.122*(tannin(:,:,:)+273.15) ! from Tarasov and Peltier 2002
       !ishc(:,:,:) = 2115.3+7.79293*(tannin(:,:,:)+273.15) ! formule Cat dans GRISLI??

        where(meltin(:,:,:).lt.snowin(:,:,:))
          refreez(:,:,:) = min (                                                                &
               rainin(:,:,:) + meltin(:,:,:),                                                   &
               capf * (snowin(:,:,:)-meltin(:,:,:)) - ishc(:,:,:) * (tal/latheat_fus_wat) * min(tannin(:,:,:), 0d0))
       elsewhere
          refreez(:,:,:) = min (                                                                &
               rainin(:,:,:) + meltin(:,:,:),                                                   &
               - ishc(:,:,:) * (tal/latheat_fus_wat) * min(tannin(:,:,:), 0d0) )
       endwhere
       
      end subroutine compute_refreezing
#endif
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      end module transfer_ecb_subgrid_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
