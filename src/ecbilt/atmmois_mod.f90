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
!      MODULE: atmmois_mod
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module atmmois_mod is collection moisture related routines initially located in atmphys0.f
!
!>     @date Creation date: Fri Nov 22 08:35:30 CET 2019
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

      module atmmois_mod

!~        use AnotherModule_mod, only: some_variable
!~        use AnotherModule_mod, only: some_otherfunction

       implicit none
       private

!~ #if (ISOATM >= 1 )
!~        public :: ec_moisfields, ec_surfmois, ec_convec_r332, ec_moisture_r332
!~ #else
       public :: ec_moisfields, ec_surfmois, ec_moisture !ec_convec defined elsewhere
!~ #endif

      ! NOTE_avoid_public_variables_if_possible

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ec_moisfields
!
!>     @brief This routine calculates the relative humidity of the moised layer
!
!      DESCRIPTION:
!
!>     *** calculates relative humidity of the moised layer
!!     *** and specific humidity above the surface and at the surface
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ec_moisfields() result(returnValue)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Use of external modules variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comatm,      only: nlat, nlon, iwater
      use comphys,     only: relhum, rmoisg
      use comsurf_mod, only: q10, q10n, nld, ntyps, lwrmois, pgroundn, tempsgn

#if ( DOWNSTS == 1 )
      use comphys,     only: temp2g, temp4g
!~ ###      use comemic_mod, only:

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr added for the temperature downscaling
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
                                       ! Adding the new variables in "_d"
                                       ! These should replace the original
                                       !  "_min" and "_max" variables
      use vertDownsc_mod,   only: q10_d, q10n_d, relhum_d, pground_d,tempsg_d

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr added for the precipitation downscaling
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use ecbilt_topography,only: nb_levls, qmount_virt
      use atmphys_d,        only: ec_ptmoisgp_d
#endif

#if (ISOATM >= 1 )
      use iso_param_mod, ONLY : ieau, ieau18
#endif

#if ( DOWNSTS == 1 )
      use comdyn
#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in] void
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical         :: returnValue

       integer          :: i,j,nn ! indices
       double precision :: ec_qsat,pmount,tmount,qmax,dqmdt


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  New variables for vertical downscaling on nb_ipoints virtual surfaces
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

#if ( DOWNSTS == 1 )

       logical                                 :: success = .false.
       integer                                 :: nb_down ! loop indice over nb_ipoints
       double precision, dimension(nb_levls)   :: tmount_d,qmax_d,dqmdt_d

#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if (ISOATM >= 1 )

! dmr   STD concistency check

      do j=1,nlon
        do i=1,nlat
          if ((rmoisg(i,j,ieau).GT.EPSILON(rmoisg(i,j,ieau))).and.(rmoisg(i,j,ieau18).EQ.0.d0))then
           write(*,*) "PB RMOISG ec_moisfields0 !! ", rmoisg(i,j,ieau), rmoisg(i,j,ieau18),i,j
           write(*,*) "valeur ieau18, iconvn = ", ieau18
           read(*,*)
          endif

        enddo
      enddo
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr  Spatial loop for everybody on all cells at first

      do j=1,nlon
        do i=1,nlat

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  ec_ptmoisgp:
! dmr  Computed on the normal surface
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr    outputs
! dmr              tmount
! dmr              qmax
! dmr
! dmr    inputs
! dmr              dqmdt
! dmr              qmount (provided as include through comphys)
! dmr
! dmr              pmount is unused and is always 0.0d0 here [ deleted ]

          call ec_ptmoisgp(tmount,qmax,i,j,dqmdt)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  New relative humidity computations on the normal surface
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          if (qmax.gt.0d0) then

!~ #if (ISOATM >= 1 )
!~             relhum(i,j)=min(1d0,rmoisg(i,j,ieau)/qmax)
!~ #else
            relhum(i,j)=min(1d0,rmoisg(i,j,iwater)/qmax)
!~ #endif

          else
            relhum(i,j)=0d0
          endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  New "q10n" calculations by surface type on normal surface
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          q10(i,j)= 0.d0
          do nn=1,ntyps
            q10n(i,j,nn)=relhum(i,j) * ec_qsat(pgroundn(i,j,nn),tempsgn(i,j,nn))
          enddo


! *** lwrmois is used in the lwr parameterization

!dqa      lwrmois(i,j)=rmoisg(i,j)**0.3333

!~ #if (ISOATM >= 1 )
!~ !dmr --- In case rmoisg turns negative (may happen in the isotopic
!~ !version I suppressed the upcasting to zero), let's keep it to zero for
!~ !the longwave radiation at least ...

!~           lwrmois(i,j)=max(rmoisg(i,j,ieau),0.0d0)
!~           lwrmois(i,j)=rmoisg(i,j,ieau)
!~ #else
          lwrmois(i,j)=rmoisg(i,j,iwater)
!~ #endif

!
! dmr  [TOBEDONE] The following section could be moved to a standalone
!      ecmoisfields_d subroutine in atmphys_d I think ...
!

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Vertical downscaling FROM HERE BELOW UNTIL END OF SPATIAL LOOP
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( DOWNSTS == 1 )

! dmr [TODO] need to impose to be on land, why is that? is anything depending on it?
         nn = nld
! dmr [NOTA] nn is used in q10n_d below ... move that assertion below?

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Vertical downscaling counter part of ec_ptmoisg
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr [TODO] Update this part to include the definition of the rmount_virt
! dmr         and suppress link to qmound_d that has disappeared!

         success = ec_ptmoisgp_d(temp2g(i,j),temp4g(i,j),                           &
                  geopg(i,j,2),tmount_d,qmax_d,dqmdt_d,qmount_virt(:))

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Below block code addition to add the virtual surfaces computations
!      NOTA: the original calculation is still performed, in order to keep
!       consistency if the sub-grid has only partial coverage of the global one
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  New relative humidity computations on the virtual surfaces
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          do nb_down = 1, nb_levls

            if (qmax_d(nb_down).gt.0.0d0) then

!~ #if (ISOATM >= 1 )
!~               relhum_d(i,j,nb_down)=min(1.0d0,rmoisg(i,j,ieau)/qmax_d(nb_down))
!~ #else
              relhum_d(i,j,nb_down)=min(1.0d0,rmoisg(i,j,iwater)/qmax_d(nb_down))
!~ #endif
            else
              relhum_d(i,j,nb_down)=0d0
            endif ! on qmax_d(nb_down) > 0.0d0

          enddo ! on nb_levls

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  New "q10n" calculations by surface type on virtual surfaces
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          q10_d(i,j,:)= 0.d0

            do nb_down = 1, nb_levls
              q10n_d(i,j,nn,nb_down)=relhum_d(i,j,nb_down) * ec_qsat(pground_d(i,j,nb_down),tempsg_d(i,j,nb_down))
            enddo ! on nb_ipoints

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  CLOSING virtual downscaling SECTION
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#endif /* on DOWNSTS == 1 */

! dmr  Finalizing spatial loop here ...
        enddo
      enddo

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   STD concistency check
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if (ISOATM >= 1 )
      do j=1,nlon
        do i=1,nlat
          if ((rmoisg(i,j,ieau).GT.EPSILON(rmoisg(i,j,ieau))).and.(rmoisg(i,j,ieau18).EQ.0.d0))then
           write(*,*) "PB RMOISG ec_moisfields9 !! ", rmoisg(i,j,ieau),rmoisg(i,j,ieau18),i,j
           write(*,*) "valeur ieau18, iconvn = ", ieau18
           read(*,*)
          endif
        enddo
      enddo
#endif


        returnValue = .true.

        return

      end function ec_moisfields

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ec_surfmois(nn)
!
!>     @brief calculates specific humidity at the surface
!
!      DESCRIPTION:
!
!>
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ec_surfmois(nn) result(returnValue)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Use of external modules variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use global_constants_mod, only: dblp=>dp, ip
      use comatm,               only: nlat, nlon
      use comphys,              only:
      use comrunlabel_mod,      only: irunlabelf
      use comsurf_mod,          only: qsurfn, pgroundn, tsurfn

!~ #if ( DOWNSTS == 1 )
!~       use comemic_mod, only:
!~ #endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Downscaling input variables for the sub-grid computations
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( DOWNSTS == 1 )
      use comsurf_mod,         only: nld
      use vertDownsc_mod,      only: pground_d, tsurfn_d, qsurfn_d
      use ecbilt_topography,   only: nb_levls
#endif

      implicit none


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in] nn land surface type
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer(kind=ip), intent(in) :: nn


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical         :: returnValue

       integer(kind=ip):: i, j, nb_down ! loop indice over nb_ipoints
       real(kind=dblp) :: ec_qsat


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


        do j=1,nlon
          do i=1,nlat

            qsurfn(i,j,nn)=ec_qsat(pgroundn(i,j,nn),tsurfn(i,j,nn))

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Downscaling section
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( DOWNSTS == 1 )

!     --- The other two surface types (sea-ice and ocean) do not admit sub-grid orography.
! dmr --- Moreover, flag_down should only exist if ntyp == nland ...
! dmr --- This section of code is called only if nn == nld anyhow

            if (nn.eq.nld) then
              do nb_down = 1, nb_levls
                 qsurfn_d(i,j,nn,nb_down)= ec_qsat(pground_d(i,j,nb_down),tsurfn_d(i,j,nn,nb_down))
              enddo
            endif ! on nn == nld
#endif

         enddo ! on the spatial loops ...
        enddo

        returnValue = .true.

        return

      end function ec_surfmois

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ec_moisture
!
!>     @brief advection and sources and sinks of atmospheric moisture
!
!      DESCRIPTION:
!
!>
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ec_moisture() result(returnValue)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Use of external modules variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      use global_constants_mod, only: dblp=>dp, ip
      use comatm, only: nlat, nlon, nsh, nsh2, nvl, dtime, grav, nm, ntl, plevel, iwater, nwisos
      use comdyn, only: omegg,pp
      use comphys,only: rmoisg, rowat, tzero, rmoiss, relhum, rlatsub, rlatvap, cormois, vhforg1, thforg1, dysnow, dyrain    &
                      , cpair, dp1, ihavm, ivavm, imsink

       use comsurf_mod, only: tsurf, evap
#if (ISOATM >= 1 )
      use comphys,only: corain, cosnow, temp4g
#endif /* ISOATM */

! apparently nothing depends on comemic
!~ ###       use comemic_mod, only:

#if (ISOATM >= 1 )
       use iso_param_mod,only: neauiso, ieau, ieau18, rsmow
       use iso_alphas,   only: alpha_lv, alpha_sv
       use iso_funcs,    only: rayleigh_dist, mois2ratiosMJ
       use isoatm_mod,   only: datmini

#if ( FRAC_KINETIK == 1 )
       use iso_alphas,   only: alpha_sve
#endif /* FRAC_KINETIC */

#endif /* ISOATM */


! TODO => add properly the following section for downscaling

#if ( DOWNSTS == 1 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr added for the precipitation downscaling
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use ecbilt_topography,only: nb_levls
      use vertDownsc_mod,   only: relhum_d, tsurf_d
#endif

#if ( DOWNSCALING == 2 )
      use input_subgrid2L,  only: sub_grid_notflat,nbpointssg    ! mask : is 0.0 when sub-grid is flat and nothing is to be done
      use vertDownsc_mod,   only: tsurf_d_interp_sg
      use input_subgrid2L,  only: dyrain_sg, dysnow_sg, tsurf_sg, weights_low_sg, index_low_sg, max_nb_points, nneigh        &
                               , index_interpL2G, area_sg, area_sg_onecb, weights_interpL2G, sumweights_interpL2G            &
                               , cluster_nlev_subgrid_in_ecbilt_dble
      use taillesGrilles,   only: sgnx,sgny,sgnxm,sgnym,sgd
      use interpolate_mod,  only: interpolateD
#endif

#if ( DOWNSTS == 1 )
      use vertDownsc_mod,   only: dyrain_d, dysnow_d
#endif

!! USE OMP_LIB

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in] void
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical         :: returnValue

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


#define DEBUGLEV 0
#define CHECKSLEV 0


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer(kind=ip)                       ::  i,j
      real(kind=dblp), dimension(nlat,nlon)  ::  d1, d3, d2
      real(kind=dblp)                        ::  ec_globalmean, qstar
      real(kind=dblp)                        ::  ec_levtempgp, factemv, factems, omegg500, t500, ec_qsat, gm1, gm2

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables with pre-processing dependence
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  0.1  Water isotopes
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( ISOATM >= 1 )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Following section contains variable that are different in iso/non-iso
!        versions.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~       real(kind=dblp), dimension(nlat,nlon,neauiso) :: hdmoisg, hdivmg, vemoisg

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Now variables specific to the isotopic enabled version
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer(kind=ip) :: k
      real(kind=dblp), dimension(nlat,nlon,neauiso) :: rmoiseauav

      real(kind=dblp), dimension(neauiso):: rmoisg_startstep
      real(kind=dblp), dimension(neauiso):: variso, variso2
      real(kind=dblp), dimension(neauiso):: isofluxiso
      real(kind=dblp)                    :: alpha_rayl, fract, rmoisg_semistep, ratio_av, ratio_ap

!~ #else /* ISOATM == 0 */
#endif /* on ISOATM >= 1 */
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  non-iso local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=dblp), dimension(nlat,nlon,nwisos) ::  hdmoisg, hdivmg, vemoisg

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  0.2  Precipitation downscaling
!           Series of variables for looping over the pertaining sub-grid points
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#if ( DOWNSTS == 1 )
      integer(kind=ip)                      :: nl        ! indice on vertical levls
      real(kind=dblp), dimension(nb_levls)  :: qstar_d   ! humidity on vertical levls
      real(kind=dblp), dimension(nb_levls)  :: vemoisg_d ! moisture removal on vertical levls
      real(kind=dblp)                       :: weight_low, dyrain_grisli, dysnow_grisli
      integer(kind=ip)                      :: ind_low
      integer(kind=ip)                      :: n_point

#if ( ISOATM >= 1 )
      real(kind=dblp), dimension(nlat,nlon) :: fract_rain
#endif

#endif
#if ( DOWNSCALING == 2 )
      real(kind=dblp), dimension(max_nb_points) :: vemoisg_grisli, tsurf_grisli,verain_grisli, vesnow_grisli
      real(kind=dblp),dimension(nb_levls,sgnxm,sgnym,sgd)   :: tsurf_d_interp
      real(kind=dblp),dimension(nlat,nlon,max_nb_points) :: tsurf_interp_sg_tmp
      real(kind=dblp), dimension(nb_levls,nlon,nlat) :: tsurf_temp

      integer(kind=ip) :: nbd

      logical                                 :: success = .false.
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Start of main code of the subroutine
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! [TODO] this seems to be always constant => move to common / module?
      factemv=rlatvap*grav*rowat/cpair
      factems=rlatsub*grav*rowat/cpair

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   STD concistency check
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if (ISOATM >= 1 )
      do j=1,nlon
        do i=1,nlat
          if ((rmoisg(i,j,ieau).GT.EPSILON(rmoisg(i,j,ieau))).and.(rmoisg(i,j,ieau18).EQ.0.d0))then
           write(*,*) "PB RMOISG ec_moisture0 !! ", rmoisg(i,j,ieau),dysnow(i,j,ieau18),i,j
           write(*,*) "snow_sum", dysnow(i,j,ieau)+cosnow(i,j,ieau)
           write(*,*) "rain_sum", dyrain(i,j,ieau)+corain(i,j,ieau)
           write(*,*) "valeur ieau18, iconvn = ", ieau18
           read(*,*)
          endif
        enddo
      enddo
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  1.0) Initialization of spectral moisture arrays, horizontal div. of moisture
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if (DOWNSCALING == 2 )

!     il faut aussi interpoler le champ de temperature...
!~ !$OMP PARALLEL DO
      do nl=1, nb_levls
        tsurf_temp(nl,:,:) = transpose(tsurf_d(:,:,nl))
      enddo
!~ !$OMP END PARALLEL DO

      do nbd=1,sgd
         success = interpolateD(index_interpL2G(:,:,1:sgnx(nbd),1:sgny(nbd),nbd),           &
              weights_interpL2G(:,1:sgnx(nbd),1:sgny(nbd),nbd),                           &
              sumweights_interpL2G(1:sgnx(nbd),1:sgny(nbd),nbd),                        &
              tsurf_temp,tsurf_d_interp(:,1:sgnx(nbd),1:sgny(nbd),nbd),                     &
              sgnx(nbd), sgny(nbd), 3, nneigh, nlon, nlat, nb_levls)
      enddo
      
! tsurf_d_interp_sg is dimension(nb_levls,nlat,nlon,max_nb_points)
      call cluster_nlev_subgrid_in_ecbilt_dble(tsurf_d_interp,tsurf_d_interp_sg,nb_levls)

#endif

!~ #if (ISOATM >= 1 )

!~        call ec_rggtosp(rmoisg(:,:,ieau),rmoiss(:,ieau))
!~        call ec_sptogg(rmoiss(:,ieau),rmoisg(:,:,ieau),pp)

!~ !$OMP PARALLEL DO
!~       do j=1,nlon
!~         do i=1,nlat

!~ !       Moisture correction from the grid interpolation process ...
!~ !       Due to the fact that the moisture in any of the components
!~ !       (except 16O so far?) may be zero, then to avoid fractionation,
!~ !       we set everybody to zero for each that is negative

!~           cormois(i,j,ieau) = 0.d0
!~           if (rmoisg(i,j,ieau).lt.0d0) then
!~             cormois(i,j,:)=cormois(i,j,:)-rmoisg(i,j,:)
!~             rmoisg(i,j,:)= 0d0
!~           endif

!~         enddo
!~       enddo
!~ !$OMP END PARALLEL DO

!~ ! *** horizontal divergence of moisture
!~         call ec_trafluxdiv(hdivmg(:,:,ieau),rmoiss(:,ieau), rmoisg(:,:,ieau))

!~ #else /* normal version, no isotopes */

!dmr --- Transfer of the variables from the lat-lon grid to the
!dmr ---   spectral grid
!dmr --- Wild guess: initialization of rmoiss ?
      call ec_rggtosp(rmoisg(:,:,iwater),rmoiss(:,iwater))
!dmr --- Transfer of the variables from the lat-lon grid to the
!dmr ---   spectral grid
!dmr --- Wild guess: initialization of pp
      call ec_sptogg(rmoiss(:,iwater),rmoisg(:,:,iwater),pp)

!~ !$OMP PARALLEL DO
      do j=1,nlon
        do i=1,nlat

          cormois(i,j,iwater)=0d0
          if (rmoisg(i,j,iwater).lt.0d0) then
            cormois(i,j,iwater)=cormois(i,j,iwater)-rmoisg(i,j,iwater)
            rmoisg(i,j,iwater)= 0d0
          endif

        enddo
      enddo
!~ !$OMP END PARALLEL DO

! *** horizontal divergence of moisture

      call ec_trafluxdiv(hdivmg,rmoiss(:,iwater),rmoisg(:,:,iwater))
!~ #endif /* On ISOATM >= 1 */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  1.1) Initialization of spectral moisture arrays, horizontal div. of moisture
!           This section is identical to previous, save the outer loop on isotopes
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if (ISOATM >= 1 )

      do k=ieau+2, neauiso ! skipping useless O16 ...

       call ec_rggtosp(rmoisg(:,:,k),rmoiss(:,k))
       call ec_sptogg(rmoiss(:,k),rmoisg(:,:,k),pp)


!~ !$OMP PARALLEL DO
      do j=1,nlon
        do i=1,nlat

!       Moisture correction from the grid interpolation process ...
!       Due to the fact that the moisture in any of the components
!       (except 16O so far?) may be zero, then to avoid fractionation,
!       we set everybody to zero for each that is negative

          cormois(i,j,k) = 0.d0
          if (rmoisg(i,j,k).lt.0d0) then
            cormois(i,j,ieau+2:neauiso)=cormois(i,j,ieau+2:neauiso)-rmoisg(i,j,ieau+2:neauiso)
            rmoisg(i,j,ieau+2:neauiso)= rmoisg(i,j,1)*(datmini(ieau+2:neauiso)+1.0d0)*rsmow(ieau+2:neauiso)
          endif

        enddo
      enddo
!~ !$OMP END PARALLEL DO

! *** horizontal divergence of moisture
        call ec_trafluxdiv(hdivmg(:,:,k),rmoiss(:,k),rmoisg(:,:,k))

      enddo ! on isotopes, neauiso

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr   1.2) Isotopic ratios for the atmosphere are computed here
!dmr       (out of any loop)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!       CALL mois2ratiosMJ(rmoisg(:,:,:),ration_qatm(:,:,:))

#endif /* On ISOATM >= 1 */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  2.0) Vertical advection of moisture
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      do j=1,nlon
        do i=1,nlat

!~ #if ( ISOATM >= 1 )
!~           vemoisg(i,j,ieau)=0d0
!~           ! [????] init of vemoisg(i,j,ieau+1:neauiso) ?
!~ #else
          vemoisg(i,j,iwater)=0d0
!~ #endif

#if ( DOWNSTS == 1 )
      vemoisg_d(:)     = 0.0d0
#endif

#if ( DOWNSCALING == 2 )
      tsurf_grisli(:)  = 0.0d0
      dyrain_grisli    = 0.0d0
      dysnow_grisli    = 0.0d0
#endif

          omegg500=(omegg(i,j,1)+omegg(i,j,2))/2.d0

          if (omegg500.lt.0.d0) then

            t500=ec_levtempgp(plevel(2),i,j)

            qstar=relhum(i,j)*ec_qsat(plevel(2),t500)


!~ #if ( ISOATM >= 1 )
!~             vemoisg(i,j,ieau)=-omegg500*qstar/(grav*rowat)
!~             vemoisg(i,j,ieau)=min(vemoisg(i,j,ieau),rmoisg(i,j,ieau)/dtime)
!~ #else
            vemoisg(i,j,iwater)=-omegg500*qstar/(grav*rowat)
            vemoisg(i,j,iwater)=min(vemoisg(i,j,iwater),rmoisg(i,j,iwater)/dtime)
!~ #endif /* on ISOATM >= 1 */

#if ( DOWNSTS == 1 )
! [NOTA] The following computation is on the vertical extended grid, not on the sub-grid
!        So far computationally cheap

            qstar_d(:)   = relhum_d(i,j,:)*ec_qsat(plevel(2),t500)
            vemoisg_d(:) = -omegg500*qstar_d(:)/(grav*rowat)

! dmr  Limit the rainout from the convergence to the maximum moisture availability
            do nl = 1, nb_levls
!~ #if ( ISOATM >= 1 )
!~               vemoisg_d(nl)=min(vemoisg_d(nl),rmoisg(i,j,ieau)/dtime)
!~ #else
              vemoisg_d(nl)=min(vemoisg_d(nl),rmoisg(i,j,iwater)/dtime)
!~ #endif

            enddo ! on nb_levls
#endif /* On DOWNSTS == 1 */

#if ( DOWNSCALING == 2 )
         if ( sub_grid_notflat(i,j) > 0.0d0 ) then ! sub-grid is rough ...

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Start loop on sub-grid cells ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          do n_point = 1, nbpointssg(i,j)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Here we need to broadcast vemoisg_d on the sub_grid
!      Also, we need the temperature on the sub_grid to account for the phase
!        change of the water at tzero
!      We do this column-wise to avoid broacasting the full grid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

            ind_low = index_low_sg(i,j,n_point)
            weight_low = weights_low_sg(i,j,n_point)

            vemoisg_grisli(n_point) = vemoisg_d(ind_low)*weight_low + vemoisg_d(ind_low+1)*(1.d0 - weight_low)
            tsurf_grisli(n_point) = tsurf_d_interp_sg(ind_low,i,j,n_point)*weight_low+ tsurf_d_interp_sg(ind_low+1,i,j,n_point)  &
                                    *(1.d0 - weight_low)
            tsurf_sg(i,j,n_point) = tsurf_grisli(n_point)

! afq -- [TODO] To account for elevation desertification we would need to change rmoisg,
!        via relhum_d.... Not implemented yet!

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  (see below)


       ! [NOTA] Implement here the transfer of precipitation on the GRISLI grid
       !        I.e. what needs to be done is the computation of vemoisg-rain and
       !           vemoisg-snow sub-grid scale.

       ! The subroutine on dyrain / dysnow computations uses:
             ! omegg (nlat, nlon, ?) (? = 1,2)
             ! tsurf_grisli   ! surface temperature on GRISLI for the considered n_point
             ! vemoisg_grisli ! moisture content from convergence on GRISLI  for the considered n_point
             ! The output is:
             !   verain_grisli, vesnow_grisli -> on the n_point considered in i,j

             ! The partitioning of available moisture on the sub-grid between rain and snow
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

                if (tsurf_grisli(n_point) .ge. tzero + 0.1d0 ) then
                   dyrain_grisli = dyrain_grisli + vemoisg_grisli(n_point)*area_sg(i,j,n_point)
! dmr&aq Adding the accumulation of snow/rain on the GRISLI grid, using global variable
                   dyrain_sg(i,j,n_point)=dyrain_sg(i,j,n_point) + vemoisg_grisli(n_point)
                else
                   dysnow_grisli = dysnow_grisli + vemoisg_grisli(n_point)*area_sg(i,j,n_point)
! dmr&aq Adding the accumulation of snow/rain on the GRISLI grid, using global variable
                   dysnow_sg(i,j,n_point)=dysnow_sg(i,j,n_point) + vemoisg_grisli(n_point)
                endif

         enddo ! on nbpointsISM

         dyrain_grisli = dyrain_grisli / area_sg_onecb(i,j)

#if ( ISOATM >= 1 )
         dyrain(i,j,ieau) = dyrain(i,j,ieau) + dyrain_grisli
#else
         dyrain(i,j,iwater) = dyrain(i,j,iwater) + dyrain_grisli
#endif
         thforg1(i,j)=thforg1(i,j) + factemv*dyrain_grisli/dp1
         vhforg1(i,j)=vhforg1(i,j) + factemv*dyrain_grisli/dp1

         dysnow_grisli = dysnow_grisli / area_sg_onecb(i,j)

#if ( ISOATM >= 1 )
        dysnow(i,j,ieau) = dysnow(i,j,ieau) + dysnow_grisli
#else
         dysnow(i,j,iwater) = dysnow(i,j,iwater) + dysnow_grisli
#endif
         thforg1(i,j)=thforg1(i,j) + factems*dysnow_grisli/dp1
         vhforg1(i,j)=vhforg1(i,j) + factems*dysnow_grisli/dp1

!~ #if ( ISOATM >= 1 )
!~         vemoisg(i,j,ieau) = dyrain_grisli + dysnow_grisli

!~         if (vemoisg(i,j,ieau).gt.epsilon(vemoisg(i,j,ieau))) then
!~           fract_rain(i,j) = dyrain_grisli / (vemoisg(i,j,ieau))
!~         else
!~           fract_rain(i,j) = 0.0d0
!~         endif
!~ #else
         vemoisg(i,j,iwater) = dyrain_grisli + dysnow_grisli
!~ #endif


         else ! grid is flat ...

#endif /* on DOWNSCALING == 2 */

!       Temperatures are positive, we are dealing with rain ...

            if (tsurf(i,j).ge.tzero + 0.1d0 ) then

!~ #if ( ISOATM >= 1 )

!~ ! --- Currently k = ieau, we are dealing with water ...

!~               dyrain(i,j,ieau) = dyrain(i,j,ieau) + vemoisg(i,j,ieau)
!~               thforg1(i,j)=thforg1(i,j) + factemv*vemoisg(i,j,ieau)/dp1
!~               vhforg1(i,j)=vhforg1(i,j) + factemv*vemoisg(i,j,ieau)/dp1
!~ #else
              dyrain(i,j,iwater) = dyrain(i,j,iwater) + vemoisg(i,j,iwater)
              thforg1(i,j)=thforg1(i,j) + factemv*vemoisg(i,j,iwater)/dp1
              vhforg1(i,j)=vhforg1(i,j) + factemv*vemoisg(i,j,iwater)/dp1
!~ #endif /* ISOATM >= 1 */
            else ! tsurf < tzero

!       Temperatures are negative, we are dealing with snow ...

!~ #if ( ISOATM >= 1 )
!~               dysnow(i,j,ieau) = dysnow(i,j,ieau) + vemoisg(i,j,ieau)
!~               thforg1(i,j)=thforg1(i,j) + factems*vemoisg(i,j,ieau)/dp1
!~               vhforg1(i,j)=vhforg1(i,j) + factems*vemoisg(i,j,ieau)/dp1

!~ #else
              dysnow(i,j,iwater) = dysnow(i,j,iwater) + vemoisg(i,j,iwater)
              thforg1(i,j)=thforg1(i,j) + factems*vemoisg(i,j,iwater)/dp1
              vhforg1(i,j)=vhforg1(i,j) + factems*vemoisg(i,j,iwater)/dp1
!~ #endif
            endif ! on tsurf > tzero

#if ( DOWNSCALING == 2 )
          endif ! on sub_grid_notflat(i,j) > 0.0d0
#endif /* on DOWNSCALING == 2 */

#if ( DOWNSTS == 1 )
! dmr&afq --- Addition of vertically discretized rainout / snowout
            do nl = 1, nb_levls
!       Temperatures are positive, we are dealing with rain ...
            if (tsurf_d(i,j,nl).ge.tzero + 0.1d0 ) then
              dyrain_d(i,j,nl) = dyrain_d(i,j,nl) + vemoisg_d(nl)
            else ! tsurf < tzero
              dysnow_d(i,j,nl) = dysnow_d(i,j,nl) + vemoisg_d(nl)
            endif ! on tsurf > tzero

            enddo
#endif

          endif ! on omegg500 < zero
        enddo ! spatial loop
      enddo ! spatial loop

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  3.0) Horizontal diffusion of moisture
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~ #if ( ISOATM >= 1 )
!~       call ec_hdiff(rmoiss(:,ieau),hdmoisg(:,:,ieau))
!~ #else
      call ec_hdiff(rmoiss(:,iwater),hdmoisg(:,:,iwater))
!~ #endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  3.1) Horizontal diffusion of moisture for isotopes ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( ISOATM >= 1 )

      do k=ieau+2, neauiso ! skipping useless O16 ...
! *** horizontal diffusion of moisture
        call ec_hdiff(rmoiss(:,k),hdmoisg(:,:,k))
      enddo

#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  4.0) Forward time stepping for plain water
!            = Actual modification of the moisture
!              content of the atmosphere
!
!      From here onwards, everything has been calculated for the moisture content
!       update. i.e. rmoisg (initial), hdivmg (divergence flux), hdmoisg
!       (diffusion flux), evap (source from land) and vemoisg (sink from rain/snow)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      do j=1,nlon
        do i=1,nlat

#if ( ISOATM >= 1 )

!        Saving before-advection moisture content & isotopes

        rmoiseauav(i,j,:) = rmoisg(i,j,:)
        
#endif /* on water isotopes */
!~ !       Step in plain water content only

!~         rmoisg(i,j,ieau)=rmoisg(i,j,ieau)+ dtime*(-ihavm*hdivmg(i,j,ieau) - ivavm*vemoisg(i,j,ieau) + hdmoisg(i,j,ieau)     &
!~                                          + imsink*evap(i,j,ieau))

!~ !       Moisture correction for everyone after time-stepping
!~ !       Due to the fact that the moisture in any of the components
!~ !       (except 16O so far?) may be zero, then to avoid fractionation,
!~ !       we set everybody to zero for each that is negative

!~           if (rmoisg(i,j,ieau).lt.0d0) then
!~             cormois(i,j,:)=cormois(i,j,:)-rmoisg(i,j,:)
!~             rmoisg(i,j,:)= 0d0
!~           endif

!~ #else /* no water isotopes */

!       Step in water content

          rmoisg(i,j,iwater) = rmoisg(i,j,iwater) + dtime*(              &
                 -ihavm*hdivmg(i,j,iwater) - ivavm*vemoisg(i,j,iwater) + hdmoisg(i,j,iwater) + imsink*evap(i,j,iwater) )

          if (rmoisg(i,j,iwater).lt.0d0) then
            cormois(i,j,iwater)=cormois(i,j,iwater)-rmoisg(i,j,iwater)
            rmoisg(i,j,iwater)= 0d0
          endif
!~ #endif /* on water isotopes */

        enddo ! end spatial loop
      enddo   ! end spatial loop

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!                 RAYLEIGH TYPE CALCULATIONS FOR ISOTOPES
!
!                         (dmr, 2016.04.05)
!
!       In the previous section, we have defined the amount of water
!        change in the cell due to rain/snow.
!       Here we compute the evolution of the isotopic content using
!        a Rayleigh distillation process
!
!      New version for the water isotopes ...
!       using the Rayleigh distillation formula
!      This is equivalent to the previous version
!       if the fraction of moisture removed is small,
!       and avoid problems in case a large portion
!       is removed.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if (ISOATM == 1 || ISOATM == 3 )
      WRITE(*,*)"This version of the code no longer supports isoatm = 1"
      WRITE(*,*)"And never supported isoatm = 3"
      WRITE(*,*)"<STOP>"
      STOP
#endif

#if (ISOATM >= 1)

      do j=1,nlon
        do i=1,nlat

         if (vemoisg(i,j,ieau).gt.0.0d0) then ! there is snow or rain (or both in Downscaling)

#if ( DOWNSCALING == 2 )
         if ( sub_grid_notflat(i,j) > 0.0d0 ) then ! sub-grid is rough ...

           rmoisg_startstep(:) = rmoiseauav(i,j,:)

         ! two semi-steps: one for rain, one for snow
         if (fract_rain(i,j).gt.0.0d0) then ! there is rain ...

!       Compute fraction of remaining water in the cell after the partial snow/rain step
!        used by all water isotopes

          rmoisg_semistep = rmoisg_startstep(ieau)-ivavm*dtime*vemoisg(i,j,ieau)*fract_rain(i,j)
          fract = rmoisg_semistep / rmoisg_startstep(ieau)

!dmr --- Pathetic case: fract < 0 !
          if (fract.LT.0.0) then
           fract = 0
          endif

          do k=ieau+2, neauiso ! skipping useless O16 ...

!       Store what the ratio of water isotopes in vapor before was
           ratio_av = rmoisg_startstep(k)/rmoisg_startstep(ieau)

!        Determine fractionnation factor type == rain here ...
           alpha_rayl = alpha_lv(temp4g(i,j),k)

!       Compute what the ration should be after step following Rayleigh
!        distillation in the atmospheric vapor
           ratio_ap = rayleigh_dist(fract,ratio_av,alpha_rayl)

!       Compute the new partial step isotopic content of the cell
           rmoisg(i,j,k) = rmoisg_semistep * ratio_ap

!       Compute the dyrain value for the water isotopes
           dyrain(i,j,k) = (rmoisg_startstep(k)-rmoisg(i,j,k)) / (ivavm*dtime)

          enddo ! on isotopes types

          rmoisg_startstep(:)    = rmoisg(i,j,:)
          rmoisg_startstep(ieau) = rmoisg_semistep

         endif! on fract_rain

         if (1.0d0-fract_rain(i,j).gt.0.0d0) then ! there is snow (also or not ...)

!       Compute fraction of remaining water in the cell after the partial snow/rain step
!        used by all water isotopes
          rmoisg_semistep = rmoisg_startstep(ieau)-ivavm*dtime*vemoisg(i,j,ieau)*(1.0d0-fract_rain(i,j))
          fract = rmoisg_semistep / rmoisg_startstep(ieau)

!dmr --- Pathetic case: fract < 0 !
          if (fract.LT.0.0) then
           fract = 0.0d0
          endif

          do k=ieau+2, neauiso ! skipping useless O16 ...

!       Store what the ratio of water isotopes in vapor before was
           ratio_av = rmoisg_startstep(k)/rmoisg_startstep(ieau)

!        Determine fractionnation factor type == rain here ...
#if ( FRAC_KINETIK == 0 )
            alpha_rayl = alpha_sv(temp4g(i,j),k)
#elif ( FRAC_KINETIK == 1 )
            alpha_rayl = alpha_sve(temp4g(i,j),k)
#endif /* FRAC_KINETIC */

!       Compute what the ration should be after step following Rayleigh
!        distillation in the atmospheric vapor
           ratio_ap = rayleigh_dist(fract,ratio_av,alpha_rayl)

!       Compute the new partial step isotopic content of the cell
           rmoisg(i,j,k) = rmoisg_semistep * ratio_ap

!       Compute the dyrain value for the water isotopes
           dysnow(i,j,k) = (rmoisg_startstep(k)-rmoisg(i,j,k)) / (ivavm*dtime)

          enddo ! on isotopes types

         endif ! on fract_rain

         else ! subgrid is flat ...

#endif /* on DOWNSCALING == 2 */

!       Compute fraction of remaining water in the cell after the partial snow/rain step
!        used by all water isotopes
          rmoisg_semistep = rmoiseauav(i,j,ieau)-ivavm*dtime*vemoisg(i,j,ieau)
          fract = rmoisg_semistep / rmoiseauav(i,j,ieau)

!dmr --- Pathetic case: fract < 0 !
          if (fract.LT.0.0) then
           fract = 0
          endif

          do k=ieau+2, neauiso ! skipping useless O16 ...

!       Store what the ratio of water isotopes in vapor before was
           ratio_av = rmoiseauav(i,j,k)/rmoiseauav(i,j,ieau)

!        Determine fractionnation factor type

           if (tsurf(i,j).ge.tzero + 0.1d0 ) then ! rain, use alpha_lv
            alpha_rayl = alpha_lv(temp4g(i,j),k)
           else ! snow, use alpha_sv or alpha_sve
#if ( FRAC_KINETIK == 0 )
            alpha_rayl = alpha_sv(temp4g(i,j),k)
#elif ( FRAC_KINETIK == 1 )
            alpha_rayl = alpha_sve(temp4g(i,j),k)
#endif /* FRAC_KINETIC */
           endif

!       Compute what the ration should be after step following Rayleigh
!        distillation in the atmospheric vapor
           ratio_ap = rayleigh_dist(fract,ratio_av,alpha_rayl)

!       Compute the new partial step isotopic content of the cell
           rmoisg(i,j,k) = rmoisg_semistep * ratio_ap

!       Compute the dyrain and dysnow values for the water isotopes
           if (tsurf(i,j).ge.tzero + 0.1d0 ) then ! rain, use alpha_lv
             dyrain(i,j,k) = (rmoiseauav(i,j,k)-rmoisg(i,j,k)) / (ivavm*dtime)
           else
             dysnow(i,j,k) = (rmoiseauav(i,j,k)-rmoisg(i,j,k)) / (ivavm*dtime)
           endif

          enddo ! on isotopes types

#if ( DOWNSCALING == 2 )
        endif ! on sub_grid_notflat(i,j) > 0.0d0
#endif

         endif   ! there is rain / snow


         do k=ieau+2, neauiso ! skipping useless O16 ...

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!        Finalize the step with the other isotopic fluxes ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! to suppress advection over Antarctica ...
!         if (i.LE.4) hdivmg(i,j,k) = 0.0
! to suppress advection over Antarctica ...

         isofluxiso(k) = dtime*( -ihavm*hdivmg(i,j,k) + hdmoisg(i,j,k) + imsink*evap(i,j,k) )

         rmoisg(i,j,k) = rmoisg(i,j,k) + isofluxiso(k)

!       Moisture correction for everyone after time-stepping
!       Due to the fact that the moisture in any of the components
!       (except 16O so far?) may be zero, then to avoid fractionation,
!       we set everybody to zero for each that is negative


          if (rmoisg(i,j,k).lt.0d0) then
            cormois(i,j,ieau+2:neauiso)=cormois(i,j,ieau+2:neauiso)-rmoisg(i,j,ieau+2:neauiso)
            rmoisg(i,j,ieau+2:neauiso)= rmoisg(i,j,1)*(datmini(ieau+2:neauiso)+1.0d0)*rsmow(ieau+2:neauiso)
          endif

          enddo ! on isotopes types


#if ( DOWNSCALING == 2 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   STD concistency check
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          if ((dysnow(i,j,ieau).GT.EPSILON(dysnow(i,j,ieau))).and.(dysnow(i,j,ieau18).EQ.0.d0))then
           write(*,*) "PB dysnow ec_moisture9 !! ", rmoisg(i,j,ieau),dysnow(i,j,ieau18),i,j
           write(*,*) "snow_sum", dysnow(i,j,ieau)
           write(*,*) "rain_sum", dyrain(i,j,ieau)
           write(*,*) "valeur down = ", sub_grid_notflat(i,j)
           read(*,*)
          endif
          if ((dyrain(i,j,ieau).GT.EPSILON(dyrain(i,j,ieau))).and.(dyrain(i,j,ieau18).EQ.0.d0))then
           write(*,*) "PB dyrain ec_moisture9 !! ", rmoisg(i,j,ieau),dyrain(i,j,ieau18),i,j
           write(*,*) "snow_sum", dysnow(i,j,ieau)
           write(*,*) "rain_sum", dyrain(i,j,ieau)
           write(*,*) "valeur down = ", sub_grid_notflat(i,j)
           read(*,*)
          endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#endif


        enddo ! end of spatial loop
      enddo   ! end of spatial loop

#endif /* on ISOATM >= 1 */

!~ #if ( ISOATM >= 1 )

!~ !      do k=ieau,neauiso
!~         call ec_moisbalance(ieau)
!~ !      enddo

!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!~ ! dmr   STD concistency check
!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!~         do j=1,nlon
!~         do i=1,nlat

!~           if ((rmoisg(i,j,ieau).GT.EPSILON(rmoisg(i,j,ieau))).and.(rmoisg(i,j,ieau18).EQ.0.d0))then
!~            write(*,*) "PB RMOISG ec_moisture-1 !! ", rmoisg(i,j,ieau),dysnow(i,j,ieau18),i,j
!~            write(*,*) "snow_sum", dysnow(i,j,ieau)+cosnow(i,j,ieau)
!~            write(*,*) "rain_sum", dyrain(i,j,ieau)+corain(i,j,ieau)
!~            write(*,*) "valeur ieau18, iconvn = ", ieau18
!~            read(*,*)
!~           endif

!~         enddo
!~         enddo
!~ !-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~ #else
        call ec_moisbalance(iwater)
!~ #endif

        returnValue = .true.

        return

      end function ec_moisture

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


#if (ISOATM >= 1 )

      subroutine longitudo_medius(humidae,humidae_media)
      
      use global_constants_mod, only: dblp=>dp, ip
      use iso_param_mod,only: neauiso, ieau, ieau18, rsmow
      use comphys, only: temp4g
      use iso_funcs, only: rayleigh_dist
      use iso_alphas, only: alpha_lv
      
      real(kind=dblp), dimension(:,:,:), intent(in) :: humidae         ! nlat, nlon, neauiso
      real(kind=dblp),   dimension(:,:), intent(out):: humidae_media
      integer(kind=ip), dimension(:)  , allocatable :: calculare       ! nlon
      logical,          dimension(:,:), allocatable :: personae        ! nlat, nlon
      real(kind=dblp),  dimension(:)  , allocatable :: calor_media     ! nlat
      
      integer :: i, j, k
      
      
!~       allocate(humidae_media(size(humidae,dim=2),size(humidae,dim=3)))
      allocate(  calculare(size(humidae,dim=1)))
      allocate(   personae(size(humidae,dim=1),size(humidae,dim=2)))
      allocate(calor_media(size(humidae,dim=1)))
      
      DO i=LBOUND(humidae,dim=1), UBOUND(humidae,dim=1)
        DO j=LBOUND(humidae,dim=2), UBOUND(humidae,dim=2)
           personae(i,j) = humidae(i,j,ieau).ge.epsilon(0.0d0)        !this is where the water is zero
        ENDDO
      ENDDO
      
!~       write(*,*) personae

      humidae_media(:,:) = 0.0_dblp
      calor_media(:) = 0.0_dblp
      
      calculare(:) = COUNT(personae(:,:),dim=2)
 
!~       DO i=LBOUND(humidae,dim=1), UBOUND(humidae,dim=1)           
!~         write(*,*) i, calculare(i)
!~       ENDDO
      
      DO k=LBOUND(humidae,dim=3), UBOUND(humidae,dim=3)
        DO j=LBOUND(humidae,dim=2), UBOUND(humidae,dim=2)
          DO i=LBOUND(humidae,dim=1), UBOUND(humidae,dim=1)
          
            if (personae(i,j)) then
              humidae_media(i,k) = humidae_media(i,k) + humidae(i,j,k)
            endif

          ENDDO
        ENDDO
      ENDDO
      DO j=LBOUND(humidae,dim=2), UBOUND(humidae,dim=2)
        DO i=LBOUND(humidae,dim=1), UBOUND(humidae,dim=1)
         
!~           if (personae(i,j)) then
            calor_media(i) = calor_media(i) + temp4g(i,j)/size(temp4g,dim=2)
!~           endif

        ENDDO
      ENDDO
      
      DO k=LBOUND(humidae,dim=3), UBOUND(humidae,dim=3)      
        DO i=LBOUND(humidae,dim=1), UBOUND(humidae,dim=1)       
           if (calculare(i).gt.epsilon(0.0_dblp)) then
             humidae_media(i,k) = humidae_media(i,k) / calculare(i)
           endif
        ENDDO
      ENDDO
           

      DO i=LBOUND(humidae,dim=1), UBOUND(humidae,dim=1)-1
        if (humidae_media(i,ieau).gt.epsilon(0.0_dblp)) then
!~           write(*,*) i, (humidae_media(i,ieau18)/humidae_media(i,ieau)/rsmow(ieau18)-1.0)*1000.0
        else
          humidae_media(i,ieau) = 2*epsilon(0.0_dblp) ! I create some water content, but this is *not* used in the water cycle, so no modification of water content ...
          DO k=ieau+1, UBOUND(humidae,dim=3) 
            humidae_media(i,k) = rayleigh_dist(0.1,humidae_media(i+1,ieau18)/humidae_media(i+1,ieau),alpha_lv(calor_media(i+1),k))*humidae_media(i,ieau)
          ENDDO            
!~           write(*,*) "longitudo_media" 
!~           write(*,*) "No water in longitude ...", i, humidae_media(i,ieau), humidae_media(i+1,ieau), humidae_media(i+1,ieau18) &
!~                                                    , calor_media(i), calor_media(i+1)
!~           write(*,*) "RAYLEIGH: alpha", alpha_lv(calor_media(i),ieau18),alpha_lv(calor_media(i+1),ieau18)
!~           write(*,*) "RAYLEIGH: frac ", humidae_media(i,ieau)/humidae_media(i+1,ieau), & 
!~       rayleigh_dist(0.1,humidae_media(i+1,ieau18)/humidae_media(i+1,ieau),alpha_lv(calor_media(i+1),ieau18))
!~           write(*,*) "=========="
!~           READ(*,*)
        endif
      ENDDO

      deallocate(calculare)
      deallocate(personae)
      
!~       READ(*,*)

      return
      end subroutine longitudo_medius
#endif     
 
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module atmmois_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
