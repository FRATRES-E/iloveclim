#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: veget_cell_state_mod
!
!>     @author  Didier M. Roche (dmr), refactoring by <author>
!
!>     @brief   Derived type holding the complete state of a single vegetation grid cell.
!
!>     @details All variables previously accessed as array(lat,lon) inside the computation
!>              routines (veget_sub_mod) are gathered here.  The orchestration loop in
!>              veget.f90 is responsible for:
!>                - extracting a veget_cell_state_t from the global 2-D arrays before
!>                  calling any physics routine,
!>                - writing the updated scalar values back into the 2-D arrays after
!>                  the call returns.
!>
!>              No spatial neighbours are ever accessed inside the physics routines, so
!>              this point-wise design is exact (verified against the full source tree).
!>
!>              The type contains no grid indices (ilat/ilon removed).  The only
!>              grid-dependent quantity needed by the physics is the crop fraction for
!>              the current time step, resolved once per cell by get_crop_fraction() in
!>              the orchestrator and stored in the scalar field crop_fraction.
!>              The full farea(nlat,nlon,mdfor) array remains in veget_mod.
!>
!>     @date    Created : <date>
!>     @date    Last modification : $LastChangedDate$
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module veget_cell_state_mod

        use global_constants_mod, only: dblp=>dp, ip

        implicit none
        private

        public :: veget_cell_state_t
        public :: get_crop_fraction

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        type :: veget_cell_state_t

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Deforestation forcing for the current time step
!  Resolved by get_crop_fraction() in the orchestrator before each physics call;
!  the full farea(nlat,nlon,mdfor) array stays in veget_mod.
!  9.0e+19 is the fill value meaning "no land-use constraint at this cell/time".
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          real(dblp) :: crop_fraction  ! crop fraction for current year           [–]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Local climate forcing (annual means, set from accumulation arrays before each call)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          real(dblp) :: ave_t          ! annual mean temperature                  [°C]
          real(dblp) :: ave_pr         ! annual mean precipitation                [mm yr^-1]
          real(dblp) :: ave_pr05       ! precipitation above daily threshold      [mm yr^-1]
          real(dblp) :: gdd0           ! growing degree days above 0°C            [°C day yr^-1]
          real(dblp) :: fr_ndx         ! freeze index (thaw/freeze ratio)         [–]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Vegetation cover fractions  (st + sg + sd = 1 by construction)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          real(dblp) :: st             ! tree   fraction                          [–]
          real(dblp) :: sg             ! grass  fraction                          [–]
          real(dblp) :: sd             ! desert fraction                          [–]
          real(dblp) :: snlt           ! needleleaf tree fraction                 [–]

!  Potential (equilibrium) fractions — intermediates computed by ccparam
          real(dblp) :: forshare_st    ! potential tree   fraction                [–]
          real(dblp) :: desshare_st    ! potential desert fraction                [–]
          real(dblp) :: nlshare_st     ! potential needleleaf fraction            [–]

!  Reference fractions (used by the deforestation / constant-vegetation scenarios)
          real(dblp) :: stR            ! reference tree   fraction                [–]
          real(dblp) :: sgR            ! reference grass  fraction                [–]
          real(dblp) :: sdR            ! reference desert fraction                [–]
          real(dblp) :: snltR          ! reference needleleaf fraction            [–]
          real(dblp) :: st_const       ! invariant reference tree   fraction      [–]
          real(dblp) :: sd_const       ! invariant reference desert fraction      [–]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  NPP and productivity
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          real(dblp) :: npp            ! total NPP (before CO2 enrichment)        [kgC m^-2 yr^-1]
          real(dblp) :: nppt           ! tree  NPP (after CO2 enrichment)         [kgC m^-2 yr^-1]
          real(dblp) :: nppg           ! grass NPP (after CO2 enrichment)         [kgC m^-2 yr^-1]
          real(dblp) :: pnpp           ! grid-cell net primary productivity       [kgC m^-2 yr^-1]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Leaf Area Index
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          real(dblp) :: lait           ! tree  LAI                                [m^2 m^-2]
          real(dblp) :: laig           ! grass LAI                                [m^2 m^-2]
          real(dblp) :: blai(2)        ! blai(1)=lait, blai(2)=laig (for I/O)    [m^2 m^-2]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Carbon pools — trees (per unit tree area)
!    b1t : leaves   b2t : stems+roots   b3t : fast litter   b4t : slow SOM
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          real(dblp) :: b1t, b2t, b3t, b4t   !                                   [kgC m^-2]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Carbon pools — grasses (per unit grass area)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          real(dblp) :: b1g, b2g, b3g, b4g   !                                   [kgC m^-2]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Carbon pools — grid-cell totals (area-weighted sums, used for diagnostics and I/O)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          real(dblp) :: b1, b2, b3, b4       !                                   [kgC m^-2]
          real(dblp) :: stock                 ! total carbon stock (b1+b2+b3+b4)  [kgC m^-2]
          real(dblp) :: anup                  ! annual carbon uptake              [kgC m^-2 yr^-1]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  FROG_EXP fields — originally guarded by #if ( FROG_EXP > 0 )
!  Included unconditionally; the conditional compilation guard has been removed.
!  These fields support the FROG experiment carbon flux tracking.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          real(dblp) :: Fv              ! integrated vegetation carbon flux       [kgC m^-2]
          real(dblp) :: Fv_t            ! tree  component of Fv                  [kgC m^-2]
          real(dblp) :: Fv_g            ! grass component of Fv                  [kgC m^-2]
          real(dblp) :: r_leaf          ! ratio leaf C / (leaf + wood) C         [–]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Carbon isotope pools — compiled only when CYCC == 2 (full carbon cycle with isotopes)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( CYCC == 2 )
          !  Isotope pools — trees
          real(dblp) :: b1t13, b2t13, b3t13, b4t13   ! 13C pools, trees         [kgC m^-2]
          real(dblp) :: b1t14, b2t14, b3t14, b4t14   ! 14C pools, trees         [kgC m^-2]

          !  Isotope pools — grasses
          real(dblp) :: b1g13, b2g13, b3g13, b4g13   ! 13C pools, grasses       [kgC m^-2]
          real(dblp) :: b1g14, b2g14, b3g14, b4g14   ! 14C pools, grasses       [kgC m^-2]

          !  Isotope totals
          real(dblp) :: bc13               ! grid-cell 13C stock                 [kgC m^-2]
          real(dblp) :: bc14               ! grid-cell 14C stock                 [kgC m^-2]
          real(dblp) :: anup13             ! annual 13C uptake                   [kgC m^-2 yr^-1]

          !  C4 grass
          real(dblp) :: sg4                ! C4 grass fraction within sg         [–]
          real(dblp) :: g4share_st         ! potential C4 grass fraction         [–]

          !  Minimum surface temperature (needed for C4 fraction calculation)
          !  [BUG] Originally never initialised in the source code (hardcoded to 0.0 as
          !  a workaround). Must be properly connected to the atmospheric model.
          real(dblp) :: tatmsmin           ! annual minimum surface temperature   [°C]
#endif

        end type veget_cell_state_t

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      FUNCTION: get_crop_fraction
!
!>     @brief   Resolve the crop fraction for one grid cell at the current model year.
!
!>     @details Encapsulates the index arithmetic previously duplicated inside ccstat
!>              and ccdyn (computation of indxv from irunlabelf, iyear, i0dfor).
!>              Returns the fill value (9.0e+19) when no land-use constraint applies,
!>              matching the convention used in the original farea array.
!>
!>     @param[in]  farea       crop fraction scenario array (nlat,nlon,mdfor)
!>     @param[in]  ilat        latitude  index of the target cell (1..nlat)
!>     @param[in]  ilon        longitude index of the target cell (1..nlon)
!>     @param[in]  irunlabelf  run label year offset (from comrunlabel_mod)
!>     @param[in]  iyear       current model year
!>     @param[in]  i0dfor      year-1 of first record in the scenario file (A.D.)
!>     @param[in]  ndfor       number of records actually available
!>     @param[in]  iscendef    scenario flag (0 = free, +1 = land-use, -1 = constant)
!>     @param[in]  ivegstrt    reference year for the land-use anomaly
!>     @returns    crop fraction [–], or 9.0e+19 if not applicable
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        pure function get_crop_fraction(farea, ilat, ilon,              &
                                        irunlabelf, iyear,              &
                                        i0dfor, ndfor,                  &
                                        iscendef, ivegstrt)             &
                       result(crop_frac)

          integer,    parameter                          :: mdfor_max = 1500
          real(dblp), intent(in) :: farea(:,:,:)         ! (nlat, nlon, mdfor)
          integer(ip),intent(in) :: ilat, ilon
          integer(ip),intent(in) :: irunlabelf, iyear
          integer(ip),intent(in) :: i0dfor, ndfor
          integer(ip),intent(in) :: iscendef, ivegstrt

          real(dblp) :: crop_frac

          integer(ip) :: indxv

          real(dblp), parameter :: FILL = 9.0e+19_dblp   ! no-data sentinel

          crop_frac = FILL

          if (iscendef /= 1) return   ! no land-use scenario active

          indxv = irunlabelf + iyear - i0dfor
          if (indxv < 0) return       ! before the scenario starts

          if (indxv == 0) then
            ! first year of the scenario: reference fractions will be recorded
            ! by the caller; no crop fraction to apply yet
            return
          end if

          indxv = min(indxv, ndfor)   ! clamp to last available record

          crop_frac = real(farea(ilat, ilon, indxv), dblp)

        end function get_crop_fraction



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      end module veget_cell_state_mod

! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
