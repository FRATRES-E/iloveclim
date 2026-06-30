#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: veget_mod
!
!>     @author  ??, Didier M. Roche (dmr), refactoring by <author>
!
!>     @brief   Global state for the VECODE terrestrial vegetation model.
!
!>     @details After refactoring this module holds exactly two categories of data:
!>
!>              (1) Run-control scalars — flags and integer parameters that govern the
!>                  overall behaviour of the model (initialisation mode, scenario type,
!>                  output frequencies, …).  These are set once during initialisation and
!>                  read throughout the run.
!>
!>              (2) Grid-level 2-D arrays (nlat x nlon) — required for restart I/O,
!>                  diagnostic output, and deforestation scenario data.  They are the
!>                  persistent store from which per-cell scalars are extracted at the
!>                  start of each physics call and into which updated values are written
!>                  back afterwards.
!>
!>              Physical parameters have been moved to veget_phys_params_mod.
!>              Per-cell computation state lives in veget_cell_state_mod (veget_cell_state_t).
!>              The former global indices lat/lon and the unused init_flag array
!>              have been removed.
!>
!>     @date    Original : ??, ported to F90 by dmr 06 juin 2018
!>     @date    Last modification : $LastChangedDate$
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module veget_mod

        use global_constants_mod, only: dblp=>dp, silp=>sp, ip
        use comatm,               only: nlat, nlon

        implicit none

        private :: nlat, nlon   ! grid dimensions visible internally only

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  (1) Run-control parameters and flags
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        !  Spectral mesh parameters (unchanged from original)
        integer, parameter :: nm   = 21
        integer, parameter :: nvl  = 3
        integer, parameter :: nsh  = ((nm+1)*(nm+2))/2
        integer, parameter :: nsh2 = 2*nsh

        !  Vegetation I/O unit tag (file-unit offset used in legacy I/O)
        integer, parameter :: iveg = 600

        !  Maximum number of records in the deforestation scenario file
        integer, parameter :: mdfor = 1500

        !  Run-time flags and counters
        integer(ip) :: nstat       ! vegetation model status flag
        integer(ip) :: ieqveg      ! equilibrium vegetation flag
        integer(ip) :: ieqvegwr    ! equilibrium vegetation write flag
        integer(ip) :: iscendef    ! deforestation scenario flag:
                                   !   -1 = constant vegetation,
                                   !    0 = free evolution,
                                   !   +1 = land-use scenario
        integer(ip) :: ivegstrt    ! reference year for the land-use anomaly (A.D.)

        !  Deforestation scenario time axis
        integer(ip) :: i0dfor      ! year-1 of first record in VEGET.dat (A.D.)
        integer(ip) :: ndfor       ! number of records actually present (≤ mdfor)

        !  CO2 vegetation feedback switch
        real(dblp)  :: fco2veg     ! 1.0 → carbon flux from vegetation feeds atm. CO2
                                   ! 0.0 → flux diagnosed but not fed back

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  (2a) Grid-level state arrays — vegetation cover (persistent across time steps)
!       Extracted cell-by-cell before physics calls; written back afterwards.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp), dimension(nlat,nlon) :: st   = 0.0_dblp  ! tree   fraction   [–]
        real(dblp), dimension(nlat,nlon) :: sg   = 0.0_dblp  ! grass  fraction   [–]
        real(dblp), dimension(nlat,nlon) :: sd   = 0.0_dblp  ! desert fraction   [–]
        real(dblp), dimension(nlat,nlon) :: snlt            ! needleleaf fraction [–]

        !  Reference fractions (deforestation / constant-veg scenarios)
        real(dblp), dimension(nlat,nlon) :: stR, sgR, sdR, snltR
        real(dblp), dimension(nlat,nlon) :: st_const        ! tree   reference    [–]
        real(dblp), dimension(nlat,nlon) :: sd_const        ! desert reference    [–]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  (2b) Grid-level carbon pool arrays — persistent state for restart
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp), dimension(nlat,nlon) :: b1t, b2t, b3t, b4t   ! tree  pools   [kgC m^-2]
        real(dblp), dimension(nlat,nlon) :: b1g, b2g, b3g, b4g   ! grass pools   [kgC m^-2]
        real(dblp), dimension(nlat,nlon) :: b1,  b2,  b3,  b4    ! cell totals   [kgC m^-2]

#if ( CYCC == 2 )
        real(dblp), dimension(nlat,nlon) :: b1t13, b2t13, b3t13, b4t13
        real(dblp), dimension(nlat,nlon) :: b1t14, b2t14, b3t14, b4t14
        real(dblp), dimension(nlat,nlon) :: b1g13, b2g13, b3g13, b4g13
        real(dblp), dimension(nlat,nlon) :: b1g14, b2g14, b3g14, b4g14
        real(dblp), dimension(nlat,nlon) :: bc13, bc14
        real(dblp), dimension(nlat,nlon) :: sg4        ! C4 grass fraction        [–]
        real(dblp), dimension(nlat,nlon) :: anup13     ! annual 13C uptake        [kgC m^-2 yr^-1]
#endif

!  FROG_EXP fields — originally guarded by #if ( FROG_EXP > 0 ); now unconditional.
        real(dblp), dimension(nlat,nlon) :: Fv           ! integrated C flux      [kgC m^-2]
        real(dblp), dimension(nlat,nlon) :: Fv_t, Fv_g  ! tree / grass components
        real(dblp), dimension(nlat,nlon) :: r_leaf       ! leaf / (leaf+wood) C   [–]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  (2c) Diagnostic arrays (output / restart)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp), dimension(nlat,nlon)   :: anup        ! annual carbon uptake  [kgC m^-2 yr^-1]
        real(dblp), dimension(nlat,nlon)   :: pnpp        ! net primary prod.     [kgC m^-2 yr^-1]
        real(dblp), dimension(nlat,nlon)   :: stock       ! total carbon stock     [kgC m^-2]
        real(dblp), dimension(nlat,nlon)   :: st_moy      ! time-mean tree frac.  [–]
        real(dblp), dimension(nlat,nlon,2) :: blai        ! LAI (1=tree,2=grass)  [m^2 m^-2]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  (2d) Deforestation scenario data (read once from VEGET.dat)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp), dimension(nlat,nlon,mdfor) :: farea      ! crop fraction R&F  [–]
        real(silp), dimension(mdfor)           :: VegetTime  ! scenario time axis  [yr]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  (2e) NetCDF output metadata (set during initialisation, used by veget_outp)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp), dimension(nlat)      :: phi            ! Gaussian latitudes   [rad]
        integer                          :: numvegvar       ! number of output vars
        real(dblp), dimension(80,20)     :: newvegvar       ! output variable data
        character(len=60)                :: namevegvar(80,6)! output variable names
        real(dblp)                       :: veg_fill_value
        real(dblp)                       :: veg_missing_value

      end module veget_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
