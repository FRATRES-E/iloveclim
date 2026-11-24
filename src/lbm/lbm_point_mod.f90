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
!      MODULE: [lbm_point_mod]
!
!>     @brief Land surface of the LBM at one grid point: temperature, snow and soil water, and their water isotopes.
!
!      DESCRIPTION:
!
!>     The routines see one point only (no grid, no neighbours): the grid-agnostic part of the LBM (dmr). The loop over
!!     the grid and the grid sums are in lbm_driver_mod. Routines moved verbatim from landmodel_mod (wiso:3a:0003):
!!     lbm_landtemp, lbm_landprecip (dmr's point-wise versions of ec_landtemp, ec_landprecip; wiso:3a:0000) and the
!!     isotope helpers lbm_wiso_ratio, lbm_wiso_mix (wiso:3a:0002, WISOATM = 1).
!!
!!     Contract (public entry points):
!!       lbm_landtemp     : land temperature and melt energy of a point.
!!       lbm_landprecip   : evaporation, snow fall and rain, melt and snow above dsnm of a point (fate of that snow
!!                          set by the caller with SNOW_EXCESS_*), overflow of the soil water (runoff of the point).
!!       lbm_wiso_ratio, lbm_wiso_mix (WISOATM = 1): ratios of a land store, mixing of an inflow.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module lbm_point_mod

        use global_constants_mod, only: dblp=>dp, ip

        implicit none
        private

        public :: lbm_landtemp, lbm_landprecip
        public :: SNOW_EXCESS_BUCKET, SNOW_EXCESS_OUT, SNOW_EXCESS_KEEP
#if ( WISOATM == 1 )
        public :: lbm_wiso_ratio, lbm_wiso_mix
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Fate of the snow above dsnm in lbm_landprecip (wiso:3a:0000)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        integer(ip), parameter :: SNOW_EXCESS_BUCKET = 1   !< into the soil water (bmoisg), the heat counted by the caller
        integer(ip), parameter :: SNOW_EXCESS_OUT    = 2   !< removed from the snow only (icebergs, ICEBERG == 2)
        integer(ip), parameter :: SNOW_EXCESS_KEEP   = 0   !< no cap (ICEBERG == 2 with ISM == 2)

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE lbm_landtemp(fractl,meltheat,icemask,tland,nethfxland,landheat,dsnow,tland_d, nethfxland_d, meltheat_d)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        !! computes surface land temperature given the heat fluxes provided
        !! Revamped version to be location independant

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Downscaling input variables for the sub-grid computations
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( DOWNSTS == 1 )
!      use vertDownsc_mod, only: tland_d, nethfxland_d, meltheat_d
      use ecbilt_topography, only: rmount_virt
#endif

      use comatm, only: iwater, nwisos
      use comland_mod, only: dtland, epsl, lhcap, rdtland, rlatfus, rlhcap, rowat, tzero, betam


      implicit none

        !> Land fraction [1]
      real(kind=dblp), intent(in)    :: fractl
        !> Net heat flux over land [?]
      real(kind=dblp), intent(in)    :: nethfxland
        !> Is this an icesheet (0/1)? [1]
      real(kind=dblp), intent(in)    :: icemask

        !> Landheat flux [?]
      real(kind=dblp), intent(inout) :: landheat

        !> Heat availability to melt snow [?]
      real(kind=dblp)                 , intent(out)   :: meltheat
        !> Computed temperature over the land surface [K]
      real(kind=dblp)                 , intent(inout)   :: tland
        !> Computed temperature over the land surface [K] in downscaling
      real(kind=dblp), dimension(:), optional, intent(inout)   :: tland_d

      real(kind=dblp), dimension(:), optional, intent(inout)   :: nethfxland_d
      real(kind=dblp), dimension(:), optional, intent(inout)   :: meltheat_d


        !> Snow pack [m water equivalent], all species; read only (intent(in): with intent(out) and a strided actual
        !>   argument, gfortran -O3 does not copy the values in, wiso:3a:0000)
      real(kind=dblp), dimension(nwisos), intent(in) :: dsnow

        !> Heat availability to melt ice [?]
      real(kind=dblp) :: iceheat

        meltheat = 0._dblp

        IF (PRESENT(meltheat_d)) then
            meltheat_d(:) = 0.0d0
        ENDIF

          ! If there is a significant land surface
        if (fractl.gt.epsl) then

            !> New temperature over land is netheatflux times landheatcapacity times timestep
          tland = tland + dtland*rlhcap*(nethfxland + landheat)

          IF (PRESENT(tland_d)) then
            tland_d(:) = tland_d(:) + dtland * rlhcap * (nethfxland_d(:)+landheat)
          ENDIF

          if ((dsnow(iwater).gt.0._dblp).and.(tland.gt.tzero)) then
            meltheat = min(dsnow(iwater)/betam,lhcap*(tland-tzero))
            tland    = tland - meltheat/lhcap
            meltheat = meltheat*rdtland

            if (PRESENT(meltheat_d)) then
              meltheat_d(:) = min(dsnow(iwater)/betam,lhcap*(tland_d(:)-tzero))
              tland_d(:)    = tland_d(:) - meltheat_d(:)/lhcap
              meltheat_d(:) = meltheat_d(:)*rdtland
            endif

#if ( IMSK == 1 )
              ! dmr&afq
              ! in the case where there is an ice-sheet below the snow ...
              ! dmr&afq
            if ((icemask.gt.0.9_dblp).and.(tland.gt.tzero)) then
              iceheat    = lhcap*(tland-tzero)
              tland      = tzero
            endif
#endif


          endif ! on snow

        endif ! on fractl > epsl

      return
      end subroutine lbm_landtemp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE lbm_landprecip(landheat,evapl,bmoisg,dsnow,snowf,rainf,meltheat,fractl,snow_excess,snowdiff,bmoism,runofl)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       !! Revamped routine from ec_landprecip, for one grid point
       !! Computes the fate of precipitation on land: evaporation, snow fall and rain, snow melt, snow above dsnm,
       !!   overflow of the soil water above its capacity (runoff of the point; collected by lbm_runoff, wiso:3a:0004)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comatm, only: iwater, nwisos

      use comland_mod, only:  dsnm, betam, dtland, epsl

      use newunit_mod, only: info_id
#if ( WISOATM == 1 )
      use comatm,      only: iwat16
      use comland_mod, only: lbm_wiso
#endif

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        !> Land heat [?]
      real(kind=dblp),                   intent(out) :: landheat
        !> Height of snow removed above dsnm [m]
      real(kind=dblp),                   intent(out) :: snowdiff
        !> Fate of the snow above dsnm: snow_excess_bucket (into bmoisg), snow_excess_out (removed, returned in snowdiff
        !>   only) or snow_excess_keep (no cap)
      integer,                           intent(in)  :: snow_excess
        !> Evaporation over land [m/s]
      real(kind=dblp), dimension(nwisos),intent(in)  :: evapl
        !> Precipitation over land, liquid and solid [m/s]
      real(kind=dblp), dimension(nwisos),intent(in)  :: rainf, snowf
        !> Heat from melting snow/ice [?]
      real(kind=dblp),                   intent(in)  :: meltheat
        !> Fraction of land in the grid cell [1]
      real(kind=dblp),                   intent(in)  :: fractl
        !> Capacity of the soil water (maximum bottom moisture) [m]
      real(kind=dblp),                   intent(in)  :: bmoism
        !> Runoff of the point: overflow of the soil water [m/s], all species (0 off land)
      real(kind=dblp), dimension(nwisos),intent(out) :: runofl

        !> "bottom moisture" ergo water on land [m]
      real(kind=dblp), dimension(nwisos),intent(inout) :: bmoisg
        !> Thickness of snow on land [m water equivalent]
      real(kind=dblp), dimension(nwisos),intent(inout) :: dsnow

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=dblp) ::  dsnomel,dsnovap
#if ( WISOATM == 1 )
        ! ratios of the snow and of the soil water (wiso:3a:0002)
      real(kind=dblp), dimension(nwisos) :: rsn, rso
#endif

      snowdiff = 0.0_dblp
      dsnovap  = 0.0_dblp
      runofl(:) = 0.0_dblp

      landheat = 0d0

          ! if on land then ...
      if (fractl.gt.epsl) then

#if ( WISOATM == 1 )
! dmr&clo   wiso:3a:0002 land isotopes (lbm_wiso): each store keeps a ratio per species, unchanged by what leaves it
!           (evaporation, melt, snow above dsnm: non-fractionating), mixed with what enters it (snow fall, rain, melt
!           water, snow above dsnm); the isotopes are the ratio times the water at the end. The evaporation leaves each
!           store at the ratio it had when the atmosphere computed the flux (ec_la2co, same step, ilan = 1).
        if (lbm_wiso) then
          rsn = lbm_wiso_ratio(dsnow)
          rso = lbm_wiso_ratio(bmoisg)
        endif
#endif

! dmr --- First deal with evaporation from the landsurface

        if (dsnow(iwater).le.0._dblp) then ! there is no snow here
          bmoisg(iwater)=bmoisg(iwater)-dtland*evapl(iwater) ! remove water that goes into evaporation from the land surface
                                                             ![NOTA] no check that there is enough water
        else                               ! there is snow
          dsnovap=dtland*evapl(iwater)           ! how much water needs to be removed in total?
          if (dsnovap.gt.dsnow(iwater)) then     ! if too much evaporation for the content of the snow layer
            bmoisg(iwater) = bmoisg(iwater)-(dsnovap-dsnow(iwater))           ! get the rest from the liquid water
            dsnow(iwater)=0._dblp                                             ! no snow left at all
          else
            dsnow(iwater) = dsnow(iwater)-dsnovap                             ! remove part of the snow layer to evaporation
          endif
        endif

#if ( WISOATM == 1 )
        if (lbm_wiso) then
          rsn = lbm_wiso_mix(rsn, dsnow(iwater), dtland*snowf)
          rso = lbm_wiso_mix(rso, bmoisg(iwater), dtland*rainf)
        endif
#endif
        dsnow(iwater) = dsnow(iwater) + dtland*snowf(iwater)                  ! add snow fall
        bmoisg(iwater)= bmoisg(iwater)+ dtland*rainf(iwater)                  ! add rainfall

! dmr the deal with the melting of snow if necessary

        if (meltheat.gt.0._dblp) then           ! there is energy to melt the snow

          dsnomel = dtland*betam*meltheat       ! potential melting

          if (dsnomel.gt.dsnow(iwater)) then    ! if too much melting for the existing snow
            landheat = landheat + (dsnomel-dsnow(iwater))/(dtland*betam) ! add the energy left to the land temperature
            dsnomel  = dsnow(iwater)                                     ! potential melting the whole layer
          endif

#if ( WISOATM == 1 )
          if (lbm_wiso) rso = lbm_wiso_mix(rso, bmoisg(iwater), rsn*dsnomel)
#endif
          dsnow(iwater) = dsnow(iwater)-dsnomel ! actually melting the snow
          bmoisg(iwater)=bmoisg(iwater)+dsnomel ! acutally adding the snow melt to the water layer

        endif ! on meltheat > 0.0

! *** if snowdepth above a thresshold, remove excessive snow
! *** artificially through the bottom moisture (or to icebergs: snow_excess, set by the caller).
! *** the heat involved is accounted for by the caller (heatsnown, heatsnows)

        if ((snow_excess.ne.snow_excess_keep).and.(dsnow(iwater).gt.dsnm)) then
          snowdiff=dsnow(iwater)-dsnm
          if (snow_excess.eq.snow_excess_bucket) then
#if ( WISOATM == 1 )
            if (lbm_wiso) rso = lbm_wiso_mix(rso, bmoisg(iwater), rsn*snowdiff)
#endif
            bmoisg(iwater)=bmoisg(iwater)+ dsnow(iwater)-dsnm
          endif
          dsnow(iwater) = dsnm
        endif

! [WARNING] This bit is not water conservative. This need to be fixed!!! dmr 2025-07-02
        if (bmoisg(iwater).lt.0d0) then
          if (bmoisg(iwater).lt.-1e-10) then
            write(info_id,*) 'bottom moisture less than zero '
            write(info_id,*) bmoisg(iwater)
            write(info_id,*) dtland*betam*meltheat
            write(info_id,*) dtland*snowf(iwater)
            write(info_id,*) dtland*rainf(iwater)
            write(info_id,*) dtland*evapl(iwater)
            write(info_id,*) dsnow(iwater),dsnovap
          endif
          bmoisg(iwater)=0d0
        endif

#if ( WISOATM == 1 )
        if (lbm_wiso) then
          dsnow(iwat16:nwisos)  = rsn(iwat16:nwisos)*dsnow(iwater)
          bmoisg(iwat16:nwisos) = rso(iwat16:nwisos)*bmoisg(iwater)
        endif
#endif

! dmr&clo   wiso:3a:0004 overflow of the soil water above its capacity (moved from ec_runoff_2 / ec_runoff), at the ratio of
!           the soil water (taken from the store, as before: bitwise identical)
        if (bmoisg(iwater).gt.bmoism) then
#if ( WISOATM == 1 )
          if (lbm_wiso) rso = lbm_wiso_ratio(bmoisg)
#endif
          runofl(iwater) = (bmoisg(iwater)-bmoism)/dtland
          bmoisg(iwater) = bmoism
#if ( WISOATM == 1 )
          if (lbm_wiso) then
            runofl(iwat16:nwisos) = rso(iwat16:nwisos)*runofl(iwater)
            bmoisg(iwat16:nwisos) = rso(iwat16:nwisos)*bmoisg(iwater)
          endif
#endif
        endif

      endif ! on land ...

      return
      end subroutine lbm_landprecip

#if ( WISOATM == 1 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      pure function lbm_wiso_ratio(store) result(r)
        !! Isotope ratios of a land store (amounts / water); 0 for an empty store (its ratio does not matter: no weight)
        !!   wiso:3a:0002
      use comatm, only: iwater, nwisos
      real(kind=dblp), dimension(nwisos), intent(in) :: store
      real(kind=dblp), dimension(nwisos)             :: r
      if (store(iwater).gt.0.0_dblp) then
        r(:) = store(:)/store(iwater)
      else
        r(:) = 0.0_dblp
      endif
      r(iwater) = 1.0_dblp
      end function lbm_wiso_ratio

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      pure function lbm_wiso_mix(r_old, w_old, a_in) result(r)
        !! Ratios of a land store of water w_old (taken as 0 if negative) and ratios r_old after the inflow a_in (amounts
        !!   of all species, a_in(iwater) the water); unchanged without inflow. wiso:3a:0002
      use comatm, only: iwater, nwisos
      real(kind=dblp), dimension(nwisos), intent(in) :: r_old, a_in
      real(kind=dblp),                    intent(in) :: w_old
      real(kind=dblp), dimension(nwisos)             :: r
      real(kind=dblp) :: wp
      wp = max(w_old, 0.0_dblp)
      if (wp + a_in(iwater).gt.0.0_dblp) then
        r(:) = (r_old(:)*wp + a_in(:))/(wp + a_in(iwater))
      else
        r(:) = r_old(:)
      endif
      r(iwater) = 1.0_dblp
      end function lbm_wiso_mix
#endif

      end module lbm_point_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
