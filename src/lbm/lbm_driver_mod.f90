!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2026 Didier M. Roche (a.k.a. dmr) | iLOVECLIM / FRATRES coding group

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!   Style sheet: v1.0.0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [lbm_driver_mod]
!
!>     @brief Grid driver of the LBM land surface: loop over the grid points and grid sums.
!
!      DESCRIPTION:
!
!>     Maps the grid fields of comland_mod onto the point routines of lbm_point_mod, one point at a time, and makes
!!     the sums over the grid that go with them: heat of the snow removed above dsnm (heatsnown, heatsnows), snow sent
!!     to icebergs (ICEBERG = 2). Everything that knows the grid stays here, so that lbm_point_mod does not
!!     (wiso:3a:0003; lbm_landsurface moved verbatim from landmodel_mod, wiso:3a:0000).
!!
!!     Contract (public entry points):
!!       lbm_landsurface(istep) : one land step of the surface (ec_lbm), before the runoff.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module lbm_driver_mod

        use global_constants_mod, only: dblp=>dp, ip

        implicit none
        private

        public :: lbm_landsurface

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE lbm_landsurface(istep)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        !! Land surface temperature and water of every grid point, one point at a time (lbm_landtemp, lbm_landprecip),
        !!   and the grid-wide sums that go with them (heat of the snow removed above dsnm, snow sent to icebergs).
        !! Replaces the whole-field ec_landtemp and ec_landprecip (wiso:3a:0000, dmr&clo).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comatm,        only: nlat, nlon
      use comland_mod,   only: fractl, meltheat, tland, nethfxland, landheat, dsnow, evapl, bmoisg, snowf, rainf, dareas  &
                             , betam, heatsnown, heatsnows
      use input_icemask, only: icemask
      use lbm_point_mod, only: lbm_landtemp, lbm_landprecip, snow_excess_bucket, snow_excess_out, snow_excess_keep

#if ( DOWNSTS == 1 )
      use vertDownsc_mod, only: tland_d, nethfxland_d, meltheat_d
#endif

#if ( ICEBERG == 2 && ISM != 2 )
      use comland_mod,   only: iceb_totl, sntoicebl
      use comemic_mod,   only: nstpyear
      use iceberg_mod,   only: iceberg_info_out_id
#endif

      implicit none

        !> Time step number in the year
      integer, intent(in) :: istep

      integer         :: i, j
      integer         :: snow_excess
      real(kind=dblp) :: snowdiff

      heatsnown = 0.0_dblp
      heatsnows = 0.0_dblp

#if ( ICEBERG == 2 && ISM != 2 )
      if (mod(istep,nstpyear).eq.1) then ! first timestep of the year
         sntoicebl = 0.
         iceb_totl = 0.
      endif
#endif

      do j=1,nlon
        do i=1,nlat

#if ( DOWNSTS == 1 )
          call lbm_landtemp(fractl(i,j),meltheat(i,j),icemask(i,j),tland(i,j),nethfxland(i,j),landheat(i,j),dsnow(i,j,:)     &
                          , tland_d(i,j,:), nethfxland_d(i,j,:), meltheat_d(i,j,:))
#else
          call lbm_landtemp(fractl(i,j),meltheat(i,j),icemask(i,j),tland(i,j),nethfxland(i,j),landheat(i,j),dsnow(i,j,:))
#endif

            ! fate of the snow above dsnm
#if ( ICEBERG < 2 )
          snow_excess = snow_excess_bucket
#elif ( ICEBERG == 2 && ISM != 2 )
          if (icemask(i,j).gt.0.9) then ! only grid cells with ice cover of more than 90%
            snow_excess = snow_excess_out
          else
            snow_excess = snow_excess_bucket
          endif
#else
          snow_excess = snow_excess_keep
#endif

          call lbm_landprecip(landheat(i,j),evapl(i,j,:),bmoisg(i,j,:),dsnow(i,j,:),snowf(i,j,:),rainf(i,j,:),meltheat(i,j)  &
                            , fractl(i,j), snow_excess, snowdiff)

          if (snowdiff.gt.0.0_dblp) then
            if (snow_excess.eq.snow_excess_bucket) then
              if (i.gt.15) then
                heatsnown=heatsnown+snowdiff*dareas(i)*fractl(i,j)/betam
              else
                heatsnows=heatsnows+snowdiff*dareas(i)*fractl(i,j)/betam
              endif
#if ( ICEBERG == 2 && ISM != 2 )
            else
              sntoicebl(i,j)=sntoicebl(i,j) + snowdiff ! accumulate excess snow over ice sheets [m]
#endif
            endif
          endif

        enddo
      enddo

#if ( ICEBERG == 2 && ISM != 2 )
      if (mod(istep,nstpyear).eq.0) then ! last timestep of the year
        do i=1,nlat
          do j=1,nlon
            iceb_totl = iceb_totl + sntoicebl(i,j)*dareas(i)*fractl(i,j)
          end do
        end do
        WRITE(iceberg_info_out_id,*) "yearly sum iceberg flux land-grid [m3]",iceb_totl,istep
      endif
#endif

      return
      end subroutine lbm_landsurface

      end module lbm_driver_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
