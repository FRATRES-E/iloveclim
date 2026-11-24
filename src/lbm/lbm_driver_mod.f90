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
!!     The runoff of each point (overflow of its soil water, lbm_landprecip) is collected by destination (lbm_runoff):
!!     to the outlet ocean cell (eni, enj) with ROUTEAU = 1, or over the ocean basin of its land basin with ROUTEAU = 0
!!     (wiso:3a:0004; the destinations are maps read at set-up, ec_inirunoff; only the receiving side needs areas).
!!
!!     Contract (public entry points):
!!       lbm_landsurface(istep) : one land step of the surface of every point, with its runoff (ec_lbm).
!!       lbm_runoff(istep)      : collection of the runoff to the ocean (runofo), after lbm_landsurface.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module lbm_driver_mod

        use global_constants_mod, only: dblp=>dp, ip

        implicit none
        private

        public :: lbm_landsurface, lbm_runoff

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
                             , betam, heatsnown, heatsnows, bmoism, runofl
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
                            , fractl(i,j), snow_excess, snowdiff, bmoism(i,j), runofl(i,j,:))

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

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE lbm_runoff(istep)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        !! Collection of the runoff of the land points (runofl, from lbm_landprecip) to the ocean (runofo), all species
        !!   ROUTEAU = 1 : each point to its outlet ocean cell (eni, enj), scaled by the land area of the point over the
        !!                 ocean area of the outlet; river-basin sums (flux_per_basin, yearly runoff_per_basin.txt);
        !!                 snow to icebergs with ICEBERG = 2 (yearly)                     (from ec_runoff_2)
        !!   ROUTEAU = 0 : sum over each land basin (ilabas), spread over its ocean basin (iocbas, arocbas)
        !!                                                                                 (from ec_runoff)
        !! Expressions and summation order of the former routines kept (bitwise). wiso:3a:0004, dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comatm,         only: nlat, nlon, nwisos, iwater
      use comland_mod,    only: runofl, runofo, fractl, dareas, epsl
#if ( ROUTEAU == 0 )
      use comland_mod,    only: nbasins, ilabas, iocbas, arocbas
#else
      use comemic_mod,    only: nstpyear
      use routageEAU_mod, only: eni, enj, flux_per_basin, river_basin
#endif
#if ( ICEBERG == 2 && ISM != 2 )
      use comemic_mod,    only: fracto
      use comland_mod,    only: sntoicebl, sntoicebo, iceb_toto
      use iceberg_mod,    only: iceberg_info_out_id
#endif
#if ( WISOATM == 1 )
      use comatm,         only: iwat16
      use comland_mod,    only: lbm_wiso
#endif

      implicit none

#if ( ISM == 1 )
#include "ismecv.com"
#endif

        !> Time step number in the year
      integer, intent(in) :: istep

      integer :: i, j
#if ( ROUTEAU == 0 )
      integer :: ibas
      real(kind=dblp) :: runo(nbasins,nwisos)
#else
      integer :: runoff_per_basin_txt_id
#endif

#if ( ROUTEAU == 0 )
! ***  total runoff for each land basin
      do ibas=1,nbasins
        runo(ibas,:)=0d0
      enddo

      do i=1,nlat
        do j=1,nlon
          if (fractl(i,j).gt.epsl) then
#if  ( ISM == 1 )
            if (flgism) then
              if(rmount_ism(i,j).LT.0.) then
              runo(ilabas(i,j))=runo(ilabas(i,j)) + dareas(i)*runofl(i,j)*fractl(i,j)
              endif
            else
#endif
!     MB land runoff to runoff, meridional convergence, fraction of land per box
!     MB  [m3 step/s] = [m3 step/s] + [m2] [m step/s] [1]
              runo(ilabas(i,j),iwater)=runo(ilabas(i,j),iwater) + dareas(i)*runofl(i,j,iwater)*fractl(i,j)
#if ( WISOATM == 1 )
              if (lbm_wiso) runo(ilabas(i,j),iwat16:nwisos)=runo(ilabas(i,j),iwat16:nwisos)                            &
                                                          + dareas(i)*runofl(i,j,iwat16:nwisos)*fractl(i,j)
#endif
#if  ( ISM == 1 )
            endif
#endif
          endif
        enddo
      enddo

#if ( ISM == 1 )
      if (flgisma) then
       if(runocorra.ne.0.) runo(23)=runocorra+runo(23)
      endif
      if (flgismg) then
       if(runocorrg.ne.0.) runo(5)=runocorrg+runo(5)
      endif
#endif

! ***  distribution of land runoff over the ocean
!     MB runo is in [m3 step/s]... means m3/s per timestep!
      do i=1,nlat
        do j=1,nlon
          runofo(i,j,iwater)=0.
          if (iocbas(i,j).gt.0) then
!     MB   [m step/s] = [m3 step/s]/[m2]
             runofo(i,j,iwater)= runo(iocbas(i,j),iwater)/arocbas(iocbas(i,j))
          endif
#if ( WISOATM == 1 )
          if (lbm_wiso) then
            runofo(i,j,iwat16:nwisos)=0.
            if (iocbas(i,j).gt.0) runofo(i,j,iwat16:nwisos)= runo(iocbas(i,j),iwat16:nwisos)/arocbas(iocbas(i,j))
          endif
#endif
        enddo
      enddo

#else
! ***  ROUTEAU = 1: each land point to its outlet

      runofo = 0.0d0
#if ( ICEBERG == 2 && ISM != 2)
      sntoicebo = 0.
      iceb_toto = 0.
#endif

      do i=1,nlat
        do j=1,nlon
          if ((fractl(i,j).gt.epsl)) then
            if ((eni(i,j).eq.0).or.(enj(i,j).eq.0)) then
              write(*,*) "lbm_runoff: no outlet for the land point ", i, j, eni(i,j), enj(i,j)
              stop 1
            endif
            runofo(eni(i,j),enj(i,j),iwater)=runofo(eni(i,j),enj(i,j),iwater)                                         &
            + runofl(i,j,iwater) * fractl(i,j)/(1.0-fractl(eni(i,j),enj(i,j)))                                         &
            *dareas(i)/dareas(eni(i,j))
            flux_per_basin(river_basin(i,j))=flux_per_basin(river_basin(i,j))                                         &
            + runofl(i,j,iwater) * fractl(i,j) * dareas(i)
#if ( WISOATM == 1 )
            if (lbm_wiso) then
              runofo(eni(i,j),enj(i,j),iwat16:nwisos)                                                                 &
              = runofo(eni(i,j),enj(i,j),iwat16:nwisos)                                                               &
              + runofl(i,j,iwat16:nwisos)                                                                             &
              * fractl(i,j)/(1.0-fractl(eni(i,j),enj(i,j)))                                                           &
              *dareas(i)/dareas(eni(i,j))
            endif
#endif
#if ( ICEBERG == 2 && ISM != 2 )
            if (mod(istep,nstpyear).eq.0) then ! only do conversion at the last timestep of the year
              sntoicebo(eni(i,j),enj(i,j))= sntoicebo(eni(i,j),enj(i,j))                                              &
              +sntoicebl(i,j) * fractl(i,j)/(1.0-fractl(eni(i,j),enj(i,j)))                                           &
              *dareas(i)/dareas(eni(i,j))
            endif
#endif
          endif ! on fractl > epsl
        enddo ! on j, nlon
      enddo ! on i, nlat

#if ( ICEBERG == 2 && ISM != 2 )
      if (mod(istep,nstpyear).eq.0) then ! last timestep of the year
        do i=1,nlat
          do j=1,nlon
            iceb_toto = iceb_toto + sntoicebo(i,j)*dareas(i)*fracto(i,j)
          end do
        end do
        WRITE(iceberg_info_out_id,*) "yearly sum iceberg flux ocean-grid [m3]", iceb_toto, istep
      endif
#endif

      if (mod(istep,nstpyear).eq.0) then ! end of year
        open(newunit=runoff_per_basin_txt_id, file='outputdata/land/runoff_per_basin.txt',                            &
             form='formatted', position='append', status='unknown')
        write(runoff_per_basin_txt_id,'(13E20.10)') flux_per_basin/nstpyear*31104000.0E0/1.0E9
        close(runoff_per_basin_txt_id)
        flux_per_basin = 0.0
      endif
#endif

      return
      end subroutine lbm_runoff

      end module lbm_driver_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
