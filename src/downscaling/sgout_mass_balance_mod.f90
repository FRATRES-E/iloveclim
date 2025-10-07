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
!      MODULE: mass_balance_mod
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module mass_balance_mod is handling the preparation of the arrays and the call of the mass balance submodules
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

      module sgout_mass_balance_mod

       use taillesGrilles,    only: iEcb, jEcb
       use ecbilt_topography, only: nb_levls

#if ( DOWNSCALING == 2 )
       use input_subgrid2L,    only: max_nb_points
#endif
       implicit none
       private

       public :: sgout_smb
       private:: set_smb_variables
       public :: akkuVars
       public :: initakkuVars

! NOTA: top-level public variables should be included here for the means of the computation of the surface mass balance
!       Also be careful that the profile variables below are *not* on the right array for interpolation, they should be
!        transposed at each level.
! afq: could be allocatable. if so, we need a subroutine init_smb (called in emic.f) to allocate what we want

       double precision, dimension(iEcb,jEcb,nb_levls), public ::  temp_profile_smb ! temperature profile ... in K
       double precision, dimension(iEcb,jEcb,nb_levls), public ::  totp_profile_smb ! total water fall (both snow and/or precip) in m.s-1 (dmr -> ?)
#if ( SMB_TYP >= 1 )
       double precision, dimension(iEcb,jEcb,nb_levls), public ::  SMB_iLCM_profile ! smb on each timestep (in m.6th_of_day-1 or m.timestep-1)
#endif


! afq: accumulated yearly:
! dmr : corrected for transposed dimensions, March, 01st, 2016
       double precision, dimension(jEcb,iEcb,nb_levls), public :: temp_profile_smb_ann !for Tann in GRISLI
       double precision, dimension(jEcb,iEcb,nb_levls), public :: totp_profile_smb_ann !for ACC in GRISLI (not used for smb)
#if ( SMB_TYP >= 1 )
 !-- ITM method: (could be between flags)
       double precision, dimension(iEcb,jEcb,nb_levls), private :: shortws_profile_smb ! is computed in set_smb_variables
       double precision, dimension(jEcb,iEcb,nb_levls), public :: SMB_iLCM_profile_ann !for SMB in GRISLI
#endif

! dmr  Compteur de passage dans l accumulateur
       double precision, private :: counter_pass

#if ( DOWNSCALING == 2 )

       double precision, dimension(iEcb,jEcb,max_nb_points), public ::  temp_sg_smb ! temperature profile ... in K
#if ( DOWN_T2M == 1 )
       double precision, dimension(iEcb,jEcb,max_nb_points), public ::  t2m_sg_smb ! temperature profile ... in K
#endif
       double precision, dimension(iEcb,jEcb,max_nb_points), public ::  totp_sg_smb ! total water fall (both snow and/or precip) in m.s-1 (dmr -> ?)
       double precision, dimension(iEcb,jEcb,max_nb_points), public ::  relhum_sg_smb ! relative humidity in ?
#if ( ISM >= 2 || SMB_TYP >= 1 )
       double precision, dimension(iEcb,jEcb,max_nb_points), public ::  SMB_iLCM_sg ! smb on each timestep (in m.6th_of_day-1 or m.timestep-1)
!-- ITM method: (could be between flags)
       double precision, dimension(iEcb,jEcb,max_nb_points), private :: shortws_sg_smb ! is computed in set_smb_variables
#if ( REFREEZING == 1 )
       double precision, dimension(iEcb,jEcb,max_nb_points), public ::  rain_sg_smb
       double precision, dimension(iEcb,jEcb,max_nb_points), public ::  snow_sg_smb
       double precision, dimension(iEcb,jEcb,max_nb_points), public ::  melt_sg_smb
#endif
#endif

! afq: accumulated yearly:
! dmr : corrected for transposed dimensions, March, 01st, 2016
       double precision, dimension(iEcb,jEcb,max_nb_points), public :: temp_sg_smb_ann !for Tann in GRISLI
#if ( DOWN_T2M == 1 )
       double precision, dimension(iEcb,jEcb,max_nb_points), public :: t2m_sg_smb_ann
#endif
       double precision, dimension(iEcb,jEcb,max_nb_points), public :: totp_sg_smb_ann !for ACC in GRISLI (not used for smb)
       double precision, dimension(iEcb,jEcb,max_nb_points), public :: relhum_sg_smb_ann !not for GRISLI, rel. humidity
#if ( SMB_TYP >= 1 )
       double precision, dimension(iEcb,jEcb,max_nb_points), public :: SMB_iLCM_sg_ann !for SMB in GRISLI
#if ( REFREEZING == 1 )
       double precision, dimension(iEcb,jEcb,max_nb_points), public ::  rain_sg_smb_ann
       double precision, dimension(iEcb,jEcb,max_nb_points), public ::  snow_sg_smb_ann
       double precision, dimension(iEcb,jEcb,max_nb_points), public ::  melt_sg_smb_ann
#endif
#endif

#endif


      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: sgout_smb
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

      subroutine sgout_smb()

#if ( SMB_TYP == 0 && ISM >= 2 )
        use ablation_PDD, only: daily_melt_pdd
        use ecbilt_topography, only: rmount_virt
#endif
#if ( SMB_TYP == 1 )
        use ablation_ITM, only: Potential_daily_melt
#elif ( SMB_TYP == 2 )
! afq -- so far the bias correction only works for DOWNSCALING == 2
        use ablation_ITM, only: Potential_daily_melt, Potential_daily_melt_biascorr
        use input_subgrid2L, only: annbias_sg
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! aurel By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! aurel Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( SMB_TYP == 2 )
! afq -- we apply to the melt coef a factor: 1 + biasfactor x bias
!        where bias about 10deg, biasfactor=0.1 gives 100pc change
!        i.e. c_rad = -40 will become -80
       double precision, parameter :: biasfactor = 0.1d0
#endif

       double precision, dimension(iEcb,jEcb,nb_levls)      :: meltrate_profile_smb

#if ( DOWNSCALING == 2 )
       double precision, dimension(iEcb,jEcb,max_nb_points) :: meltrate_sg_smb
#endif

       real    :: returnValue

       integer :: nlev


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the subroutine starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   First needs to transmit the variables to the right array for calculating the SMB
!        Integer given can be used to choose the variables depending on the SMB_type
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       call set_smb_variables(0)

#if ( DOWNSCALING == 2 )
       call set_smb_variables(1)
#endif 

#if ( ISM >= 2 || SMB_TYP >= 1 )

! --- dmr&es Shortwave surface radiation (W.m^-2), Surface air temperature (°C)
!     NOTA BENE: shortw_surf already includes the albedo at surface effect, so it is a NET shortwave surface flux
! --- dmr&es Melt rate (m.day^-1)

#if ( SMB_TYP == 0 && ISM >= 2 )
       call daily_melt_pdd(meltrate_profile_smb,rmount_virt,temp_profile_smb-273.15d0,totp_profile_smb(:,:,:)*(24.0d0*60.*60.))
#elif ( SMB_TYP >= 1 )
! afq -- so far the bias correction only works for DOWNSCALING == 2
       call potential_daily_melt(meltrate_profile_smb,shortws_profile_smb,temp_profile_smb-273.15d0)
#endif
       !smb, on vertical levels:
       !dmr --- Temporary dirty addition of unit change to conform with the SMB unit of m.6th_of_day-1
       SMB_iLCM_profile(:,:,:) = totp_profile_smb(:,:,:)*(24.0d0/6.0d0*60.*60.)-meltrate_profile_smb(:,:,:)/6.0d0

#if ( DOWNSCALING == 2 )
#if ( SMB_TYP == 0 )
! dmr&afq == Following call does not work ... to be updated [2016-04-28]
!~        call daily_melt_pdd(meltrate_sg_smb,rmount_virt,temp_sg_smb-273.15d0,totp_sg_smb(:,:,:)*(24.0d0*60.*60.))
       WRITE(*,*) "Combination of options does not work in SMB routines <STOP>"
       STOP
#elif ( SMB_TYP == 1 )
       call potential_daily_melt(meltrate_sg_smb,shortws_sg_smb,temp_sg_smb-273.15d0)
#elif ( SMB_TYP == 2 )
       call potential_daily_melt_biascorr(meltrate_sg_smb,shortws_sg_smb,temp_sg_smb-273.15d0,annbias_sg,biasfactor)
#endif

       !smb, on vertical levels:
       !dmr --- Temporary dirty addition of unit change to conform with the SMB unit of m.6th_of_day-1
       SMB_iLCM_sg(:,:,:) = totp_sg_smb(:,:,:)*(24.0d0/6.0d0*60.*60.)-meltrate_sg_smb(:,:,:)/6.0d0

#if ( REFREEZING == 1 )
       melt_sg_smb(:,:,:) = meltrate_sg_smb(:,:,:)/(24.*60.*60.) ! afq -- same unit as precip...
#endif

#endif  /* on DOWNSCALING == 2  */
#endif  /* on ISM>=2 */

      end subroutine sgout_smb

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: set_smb_variables
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

      subroutine set_smb_variables(shape_grid)

#define TEST_DOWNSC 1
#define SNOW_AGING 1

       use taillesGrilles, only : sgnx,sgny,sgnxm,sgnym,sgd
        
       use vertDownsc_mod, only: temp_profile, rmount_d !rmount is needed by the ITM
       use comatm,         only: nlat, nlon, nvl, nsh2, iwater  ! needed for the include of comphys, to get the precipitation
       ! -- afq    we need comphys for the precip and comsurf for the albedo
#if ( COMATM == 1 )
       use comphys, only: fswdsfc, torain, tosnow 
       use comsurf_mod, only: albes, nld
#endif

#if ( DOWNSCALING == 2 )
       use input_subgrid2L, only: max_nb_points, torain_sg, tosnow_sg, tsurf_sg, relhum_sg, nbpointsSG  &
     &      , nneigh, index_interpL2G, weights_interpL2G, sumweights_interpL2G 
#if ( DOWN_T2M == 1 )
       use input_subgrid2L, only: tempsg_sg
#endif
#if ( SMB_TYP >= 1 )
       use input_subgrid2L, only: topo_sg, difftopo_sg,epaisglace_sg, snow_age_sg
       use comland_mod, only: albsnow, albkeep
       use interpolate_mod,    only: interpolate
       use input_subgrid2L,     only: cluster_var_subgrid_in_ecbilt_dble
#endif
#endif

#if ( TEST_DOWNSC == 3 )
       use vertDownsc_mod, only: tosnow_d, torain_d
#endif

#if ( COMATM == 0 )
! -- afq    we need comphys for the precip and comsurf for the albedo
#include "comphys.h"  /* torain tosnow */
#include "comsurf.h"  /* albes */
#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  mass_balance_type: enable the choice of mass balance type, needs different variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer, intent(in) :: shape_grid

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer :: n_dim, i, j, nbd
#if ( DOWNSCALING == 2)
        logical                     :: success = .false.
#if ( ISM >= 2 || SMB_TYP >= 1 )
        double precision, parameter :: alb_sg_fact = 0.2d-3 !modif d'albedo par metre
        double precision            :: alb_sg
        double precision,dimension (nlat,nlon)               :: albloc
        double precision,dimension (sgnxm,sgnym,sgd)         :: albloc_interp
        double precision,dimension (nlat,nlon,max_nb_points) :: alb_interp_sg
#if ( SNOW_AGING == 1 )
        double precision, dimension(nlat,nlon,max_nb_points) :: aging
        double precision,parameter            :: agefactor = 1.d-2
        double precision,parameter            :: agereducmax = 0.2d0
#endif
#endif
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the subroutine starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          if ( shape_grid.eq.0 ) then

          temp_profile_smb(:,:,:) = temp_profile(:,:,:)       ! in °K

          do n_dim=1,nb_levls

#if (CPLTYP == 1)
             totp_profile_smb(:,:,n_dim)    = tosnow(:,:,iwater)                ! in m.s-1 I believe ...
#elif (CPLTYP == 2)
             totp_profile_smb(:,:,n_dim)    = tosnow(:,:,iwater) + torain(:,:,iwater)  ! in m.s-1 I believe ...
#else
             !write (*,*) "Default accumulated subgrid: total precipitation"
             totp_profile_smb(:,:,n_dim)    = tosnow(:,:,iwater) + torain(:,:,iwater)  ! in m.s-1 I believe ...
#endif

#if ( TEST_DOWNSC == 3 )
#if (CPLTYP == 1)
           totp_profile_smb(:,:,n_dim)    = tosnow_d(:,:,n_dim)
#elif (CPLTYP == 2)
           totp_profile_smb(:,:,n_dim)    = tosnow_d(:,:,n_dim) + torain_d(:,:,n_dim)
#else
           !write (*,*) "Default accumulated subgrid: total precipitation"
           totp_profile_smb(:,:,n_dim)    = tosnow_d(:,:,n_dim) + torain_d(:,:,n_dim)
#endif
#endif

#if ( SMB_TYP >= 1 )
           shortws_profile_smb(:,:,n_dim) = fswdsfc(:,:)*(1.0d0-albes(:,:)) !could be dependent on elevation

!afq -- we could call top_to_surf_shortwave, but transmissivity should be re-tuned... if wanted, uncomment:
!           call top_to_surf_shortwave(shortws_profile_smb(:,:,n_dim),iEcb,jEcb,fswdtoa0,rmount_d(:,:,n_dim))
!           shortws_profile_smb(:,:,n_dim)=shortws_profile_smb(:,:,n_dim)*(1.0d0-albes(:,:))
#endif

          enddo

        elseif ( shape_grid.eq.1 ) then

#if ( DOWNSCALING == 2 )
          temp_sg_smb(:,:,:) = tsurf_sg (:,:,:) ! in °K
          relhum_sg_smb(:,:,:) = relhum_sg (:,:,:) ! in m/m ?
#if ( DOWN_T2M == 1 )
          t2m_sg_smb(:,:,:) = tempsg_sg (:,:,:) ! in °K
#endif
#if (CPLTYP == 1)
          totp_sg_smb(:,:,:) = tosnow_sg(:,:,:)                    ! in m.s-1 I believe ...
#elif (CPLTYP == 2)
          totp_sg_smb(:,:,:) = tosnow_sg(:,:,:) + torain_sg(:,:,:) ! in m.s-1 I believe ...
#else
          !write (*,*) "Default accumulated subgrid: total precipitation"
          totp_sg_smb(:,:,:) = tosnow_sg(:,:,:) + torain_sg(:,:,:) ! in m.s-1 I believe ...
#endif

#if ( SMB_TYP >= 1 )

#if ( SNOW_AGING == 1 )
          where(totp_sg_smb(:,:,:).gt.1.d-8)
             snow_age_sg(:,:,:) = 0.d0
          elsewhere
             snow_age_sg(:,:,:) = snow_age_sg(:,:,:) + 1.d0
          endwhere

          aging(:,:,:) = min(agefactor * (snow_age_sg(:,:,:) / 6.) , agereducmax)
#endif

          where(albkeep(:,:,nld).gt.0.)
             albloc(:,:) = albkeep(:,:,nld)
          elsewhere
             albloc(:,:) = albes(:,:)
          endwhere

          do nbd=1,sgd

             success = interpolate(index_interpL2G(:,:,1:sgnx(nbd),1:sgny(nbd),nbd),           &
                  weights_interpL2G(:,1:sgnx(nbd),1:sgny(nbd),nbd),                          &
                  sumweights_interpL2G(1:sgnx(nbd),1:sgny(nbd),nbd),                       &
                  transpose(albloc),albloc_interp(1:sgnx(nbd),1:sgny(nbd),nbd),                &
                  sgnx(nbd),sgny(nbd), 3, nneigh, nlon, nlat)

          enddo
          
          call cluster_var_subgrid_in_ecbilt_dble(albloc_interp,alb_interp_sg)
          
          do n_dim=1,max_nb_points
             do j=1,nlon
                do i=1,nlat
                   alb_sg=max(min((1+difftopo_sg(i,j,n_dim)*alb_sg_fact)*alb_interp_sg(i,j,n_dim),albsnow(i)),0.15)
#if ( SNOW_AGING == 1 )
                   if(alb_sg.gt.0.5) then
                      alb_sg=max(alb_sg-aging(i,j,n_dim),0.5)
                   endif
#endif
                   if (epaisglace_sg(i,j,n_dim).gt.10.) then 
                      alb_sg=max(alb_sg,0.5)
                   endif
                   shortws_sg_smb(i,j,n_dim) = fswdsfc(i,j)*(1.0d0-alb_sg)
                enddo
             enddo
          enddo

#if ( REFREEZING == 1 )
          rain_sg_smb (:,:,:) = torain_sg(:,:,:)
          snow_sg_smb (:,:,:) = tosnow_sg(:,:,:)
#endif

#endif 

#else
           ! following is because a call of shape_grid == 1 without DOWNSCALING == 2 is senseless
          write(*,*) "Inconsistent call for set_smb_variables <STOP>"
          STOP
#endif
        else
          write(*,*) "Unrecognized shape_grid <STOP>"
          STOP
        endif

      end subroutine set_smb_variables

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: akkuvars
!
!>     @author  Didier M. Roche  (dmr)
!>     @author  Aurelien Quiquet (afq)
!>     @brief This subroutine accumulates the given variables for end_of_year usage
!
!      DESCRIPTION:
!
!>     Accumulation of the variables for end_of_year use
!!      Creation Date         : 12 Novembre 2008
!!      Last modification     : 15 Septembre 2010, 22 janvier 2016
!!      9 fevrier 2016, afq   : the SMB is the only accumulated var. now
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine akkuVars(k,liatm)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&afq By reference variables ...
!>    @param[in]  k      flag: k = 0 means year ongoing (virtual accumulation year, could be several solar_years ...)
!>    @param[in]  liatm  Number of atmospheric timesteps per day (6 for now, 2016-03-01)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       use global_constants_mod, only: day_length=>solar_day, tzero=>tK_zero_C
       use taillesGrilles,     only: nbmois, nbjours
#if ( ISM >= 2 )
       use input_timerCplGRIS, only: timcplgrisday
!       use input_timerCplGRIS, only: freqdecoup  ! afq
#endif

       implicit none

       integer, intent(in) :: liatm, k

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      double precision, parameter :: cent = 100.0d0
!~       double precision, parameter :: day_length = solar_day

      integer :: n_dim

      if (k.eq.0) then ! annee en cours ...

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Creation du SMB de couplage version 1 :
!       on s'interesse au smb sur chaque niveaux verticaux
!-----|--1--------2---------3---------4---------5---------6---------7-|

! afq  Accumulation sur tous les niveaux verticaux:

         do n_dim=1,nb_levls

            temp_profile_smb_ann(:,:,n_dim) = temp_profile_smb_ann(:,:,n_dim) + transpose(temp_profile_smb(:,:,n_dim))
            totp_profile_smb_ann(:,:,n_dim) = totp_profile_smb_ann(:,:,n_dim) + transpose(totp_profile_smb(:,:,n_dim))
#if ( ISM >= 2 || SMB_TYP >= 1 )
            SMB_iLCM_profile_ann(:,:,n_dim) = SMB_iLCM_profile_ann(:,:,n_dim) + transpose(SMB_iLCM_profile(:,:,n_dim))
#endif

         enddo

#if ( DOWNSCALING == 2 )

            temp_sg_smb_ann(:,:,:) = temp_sg_smb_ann(:,:,:) + temp_sg_smb(:,:,:)
            totp_sg_smb_ann(:,:,:) = totp_sg_smb_ann(:,:,:) + totp_sg_smb(:,:,:)
            relhum_sg_smb_ann(:,:,:) = relhum_sg_smb_ann(:,:,:) + relhum_sg_smb(:,:,:)
#if ( DOWN_T2M == 1 )
            t2m_sg_smb_ann(:,:,:) = t2m_sg_smb_ann(:,:,:) + t2m_sg_smb(:,:,:)
#endif
#if ( ISM >= 2 || SMB_TYP >= 1 )
            SMB_iLCM_sg_ann(:,:,:) = SMB_iLCM_sg_ann(:,:,:) + SMB_iLCM_sg(:,:,:)
#if ( REFREEZING == 1 )
            rain_sg_smb_ann(:,:,:) = rain_sg_smb_ann(:,:,:) + rain_sg_smb(:,:,:)
            snow_sg_smb_ann(:,:,:) = snow_sg_smb_ann(:,:,:) + snow_sg_smb(:,:,:)
            melt_sg_smb_ann(:,:,:) = melt_sg_smb_ann(:,:,:) + melt_sg_smb(:,:,:)
#endif
#endif
#endif
         counter_pass = counter_pass + 1.0d0

      else ! fin de l annee

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Fin de l'annee virtuelle de couplage : on change le format
!-----|--1--------2---------3---------4---------5---------6---------7-|

! dmr le 01er mars 2016
!     Pour simplifier le joyeux bazar précédent en unités, je défini un compteur
!      pour me ficher du nombre de passages en principe ...
!     Variable privée : counter_pass

! afq Temp. annual mean = T accumulated / (nb of atmos. timestep passed through)
        temp_profile_smb_ann(:,:,:) = temp_profile_smb_ann(:,:,:) / counter_pass
! dmr Conversion K -> °C
        temp_profile_smb_ann(:,:,:) = temp_profile_smb_ann(:,:,:) - tzero


! afq torain is in m/s, we want m/an
! dmr Nota: precipitation is cumulative, no need to use a mean
! dmr However, we accumulated over each atmo. timestep, not in days

        totp_profile_smb_ann(:,:,:) = totp_profile_smb_ann(:,:,:) * day_length / liatm

! afq SMB is in m/timestep -> dmr nothing to do!!
        !SMB_iLCM_profile_ann(:,:,:) = SMB_iLCM_profile_ann(:,:,:)

#if ( ISM >= 2 )
! afq we take into account decoupling frequency:
        totp_profile_smb_ann(:,:,:) = totp_profile_smb_ann(:,:,:) / (dble(timCplGRISday)/360.) ! dble(freqdecoup)
        SMB_iLCM_profile_ann(:,:,:) = SMB_iLCM_profile_ann(:,:,:) / (dble(timCplGRISday)/360.) ! dble(freqdecoup)
#endif

#if ( DOWNSCALING == 2 )

! afq Temp. annual mean = T accumulated / (nb of atmos. timestep passed through)
        temp_sg_smb_ann(:,:,:) = temp_sg_smb_ann(:,:,:) / counter_pass
! dmr Conversion K -> °C
        temp_sg_smb_ann(:,:,:) = temp_sg_smb_ann(:,:,:) - tzero

! afq moisture is similar to temperature above:
        relhum_sg_smb_ann(:,:,:) = relhum_sg_smb_ann(:,:,:) / counter_pass
        
#if ( DOWN_T2M == 1 )
! afq Temp. annual mean = T accumulated / (nb of atmos. timestep passed through)
        t2m_sg_smb_ann(:,:,:) = t2m_sg_smb_ann(:,:,:) / counter_pass
! dmr Conversion K -> °C
        t2m_sg_smb_ann(:,:,:) = t2m_sg_smb_ann(:,:,:) - tzero
#endif

! afq torain is in m/s, we want m/an
! dmr Nota: precipitation is cumulative, no need to use a mean
! dmr However, we accumulated over each atmo. timestep, not in days

        totp_sg_smb_ann(:,:,:) = totp_sg_smb_ann(:,:,:) * day_length / liatm

! afq SMB is in m/timestep -> dmr nothing to do!!
        !SMB_iLCM_sg_ann(:,:,:) = SMB_iLCM_sg_ann(:,:,:)

#if ( ISM >= 2 )
! afq we take into account decoupling frequency:
        totp_sg_smb_ann(:,:,:) = totp_sg_smb_ann(:,:,:) / (dble(timCplGRISday)/360.) !dble(freqdecoup)
        SMB_iLCM_sg_ann(:,:,:) = SMB_iLCM_sg_ann(:,:,:) / (dble(timCplGRISday)/360.) !dble(freqdecoup)
#endif

#if ( REFREEZING == 1 && ISM >= 2 )
        rain_sg_smb_ann(:,:,:) = ( rain_sg_smb_ann(:,:,:) * day_length / liatm ) / (dble(timCplGRISday)/360.) !dble(freqdecoup)
        snow_sg_smb_ann(:,:,:) = ( snow_sg_smb_ann(:,:,:) * day_length / liatm ) / (dble(timCplGRISday)/360.) !dble(freqdecoup)
        melt_sg_smb_ann(:,:,:) = ( melt_sg_smb_ann(:,:,:) * day_length / liatm ) / (dble(timCplGRISday)/360.) !dble(freqdecoup)
#elif ( REFREEZING == 1 && SMB_TYP >= 1 )
        rain_sg_smb_ann(:,:,:) = ( rain_sg_smb_ann(:,:,:) * day_length / liatm )
        snow_sg_smb_ann(:,:,:) = ( snow_sg_smb_ann(:,:,:) * day_length / liatm )
        melt_sg_smb_ann(:,:,:) = ( melt_sg_smb_ann(:,:,:) * day_length / liatm )
#endif

#endif

       endif

       return
       end subroutine akkuVars
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: initakkuVars
!
!>     @brief This subroutine is adding a wonderful knew unused functionality
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE initakkuVars

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : entree
!       Variables de sortie : sortie
!-----|--1--------2---------3---------4---------5---------6---------7-|

!$ USE OMP_LIB


       IMPLICIT NONE

!$OMP PARALLEL WORKSHARE
       temp_profile_smb_ann(:,:,:) = 0.d0 !for Tann in GRISLI
       totp_profile_smb_ann(:,:,:) = 0.d0 !for ACC in GRISLI (not used for smb)
#if ( ISM >= 2 || SMB_TYP >= 1 )
       SMB_iLCM_profile_ann(:,:,:) = 0.d0 !for SMB in GRISLI
#endif

#if ( DOWNSCALING == 2 )

       temp_sg_smb_ann(:,:,:) = 0.d0 !for Tann in GRISLI
#if ( DOWN_T2M == 1 )
       t2m_sg_smb_ann(:,:,:) = 0.d0
#endif
       totp_sg_smb_ann(:,:,:) = 0.d0 !for ACC in GRISLI (not used for smb)
       relhum_sg_smb_ann(:,:,:) = 0.d0 !relative humidity, not used for GRISLI
#if ( ISM >= 2 || SMB_TYP >= 1 )
       SMB_iLCM_sg_ann(:,:,:) = 0.d0 !for SMB in GRISLI

#if ( REFREEZING == 1 )
       rain_sg_smb_ann(:,:,:) = 0.d0
       snow_sg_smb_ann(:,:,:) = 0.d0
       melt_sg_smb_ann(:,:,:) = 0.d0
#endif

#endif

#endif
!$OMP END PARALLEL WORKSHARE

       counter_pass                = 0.d0 !counter of accumulation

       END SUBROUTINE initakkuVars

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     end module sgout_mass_balance_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
