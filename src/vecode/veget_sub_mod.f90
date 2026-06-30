!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2026 iLOVECLIM / LUDUS coding group

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: veget_sub_mod
!
!>     @author  Didier M. Roche (dmr), refactoring by dmr, clo
!
!>     @brief   Point-wise physics routines for the VECODE terrestrial vegetation model.
!
!>     @details All routines operate exclusively on a single veget_cell_state_t value.
!>              No global grid indices (lat, lon) are accessed.  The orchestration loop
!>              in veget.f90 is responsible for extracting cell state before the call
!>              and writing it back afterwards.
!>
!>              Changes from original veget_sub_mod.f90:
!>                - All array(lat,lon) accesses replaced by cell%field scalar accesses.
!>                - Physical parameters moved to veget_phys_params_mod; imported via use.
!>                - fracgr / darea arguments removed (were only used in the deprecated
!>                  CYCC==1 / LOCH carbon cycle path, now deleted).
!>                - PF_CC permafrost blocks removed (abandoned implementation).
!>                - FROG_EXP conditional removed; fields included unconditionally.
!>                - CYCC conditionals retained for the isotope carbon cycle.
!>                - CYCC==1 (LOCH model) block removed (deprecated, commented-out in
!>                  veget.f and marked [DEPRECATED] in original veget_mod).
!>                - ccparam is an internal helper; it is public only so the isotope
!>                  restart-init loop in veget can set residence times without ccstat.
!>                - indxv / farea lookup replaced by cell%crop_fraction resolved by
!>                  get_crop_fraction() in the orchestrator.
!
!>     @date    Original creation : October 17th, 2019 (dmr)
!>     @date    Refactored        : 2026-06-30 (dmr, clo)
!>     @date    Last modification : $LastChangedDate$
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module veget_sub_mod

        use global_constants_mod,  only: dblp=>dp, ip, str_len, stdout
        use veget_cell_state_mod,  only: veget_cell_state_t
        use veget_phys_params_mod

        implicit none
        private

        public :: initcpar, ccstat, ccstatR, ccdyn, ccdynR
#if ( CYCC == 2 )
        public :: ccstat_isotope
        ! ccparam exposed so the isotope-restart init loop can set the per-cell
        ! residence times (t1t..t4g) WITHOUT calling full ccstat, which would
        ! overwrite the b* pools just read from the vegetation restart.
        public :: ccparam
#endif

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      FUNCTION: initcpar
!
!>     @brief  Initialise all physical parameters (veget_phys_params_mod) to their
!>             default values.  Called once at model start-up before any grid cell
!>             is processed.
!
!>     @returns .true. on success
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine initcpar()

#if ( CYCC == 2 )
        ! c13frac, c13frac4 now in veget_phys_params_mod (imported at module level)
        use C_res_mod, only: c13init
#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! --- forest fraction potential ---
        a        = 7000.0_dblp
        bet      = 0.002_dblp
        gamm     = 0.00017_dblp
        gamm2    = 0.00025_dblp
        fmax     = 0.9_dblp
        gdd0_min = 1000.0_dblp
        gdd0_max = 1800.0_dblp
        ades     = 0.0011_dblp
        acr      = 28.0_dblp

        ! --- NPP (Lieth formula) ---
        nppmax = 1.3_dblp
        v1     = 0.000664_dblp
        v2     = 0.119_dblp
        v3     = 3.73_dblp

        ! --- carbon allocation — trees ---
        c1t = 0.046_dblp ;  c2t = 0.58_dblp  ;  c3t = 1.6_dblp
        d1t = 0.22_dblp  ;  d2t = 7.19_dblp  ;  d3t = 5.5_dblp
        e1t = 17.9_dblp  ;  e2t = 167.3_dblp ;  e3t = 15.0_dblp
        k0t = 0.6_dblp
        k2t = 1.0_dblp
        k3t = 0.025_dblp

        ! --- carbon allocation — grasses ---
        c1g = 0.069_dblp ;  c2g = 0.38_dblp  ;  c3g = 1.6_dblp
        d1g = 0.6_dblp   ;  d2g = 0.41_dblp  ;  d3g = 6.0_dblp
        e1g = 0.67_dblp  ;  e2g = 50.5_dblp  ;  e3g = 100.0_dblp
        k0g = 0.2_dblp
        k2g = 0.55_dblp
        k3g = 0.025_dblp
        k4g = 0.025_dblp

        ! --- soil temperature decomposition ---
        ps5   = 0.04_dblp
        soilt = 5.0_dblp

        ! --- needleleaf thresholds ---
        t1tn = 4.0_dblp
        t1td = 1.0_dblp

        ! --- LAI density parameters ---
        deng  = 20.0_dblp
        dentd = 20.0_dblp
        dentn = 6.0_dblp

        ! --- desert propagation densities ---
        deng_prop  = 20.0_dblp
        dentd_prop = 20.0_dblp
        dentn_prop = 6.0_dblp

        ! --- soil water / stomatal resistance ---
        acwd = 100.0_dblp ;  acwt = 100.0_dblp
        acwg = 100.0_dblp ;  acwn = 100.0_dblp
        zrd  = 1.0_dblp   ;  zrt  = 1.0_dblp
        zrg  = 0.6_dblp   ;  zrn  = 0.6_dblp
        rsd  = 0.0_dblp   ;  rst  = 300.0_dblp
        rsg  = 130.0_dblp ;  rsn  = 160.0_dblp

#if ( CYCC == 2 )
        ! --- carbon isotope initialisations ---
        c13frac  = 1.0_dblp - 18.0_dblp / 1000.0_dblp
        c13frac4 = 1.0_dblp -  5.0_dblp / 1000.0_dblp
        c13init  = c13frac
#endif


      end subroutine initcpar

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      FUNCTION: ccparam  (private)
!
!>     @brief  Compute current-year carbon cycle parameters and equilibrium fractions
!>             from the cell's climate forcing.  Writes results into the parameter
!>             module scalars (k1t, t1t, … forshare_st, …) and into cell%nppt/nppg.
!>             Called internally by ccstat, ccstatR, ccdyn, ccdynR.
!
!>     @param[inout] cell   vegetation state for the current grid cell
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine ccparam(cell)

        type(veget_cell_state_t), intent(inout) :: cell

        ! local
        real(dblp) :: npp1, npp2, avefor, differ, pcr
        real(dblp) :: db1, db2, db3

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Potential tree fraction
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! avefor = cell%ave_pr**4
        avefor = cell%ave_pr*cell%ave_pr*cell%ave_pr*cell%ave_pr
        differ = cell%gdd0 - gdd0_min
        db1    = -bet  * differ
        db2    =  gamm * differ
        db3    =  differ * differ

        if (differ < 0.0_dblp) then
          cell%forshare_st = 0.0_dblp
        else
          cell%forshare_st = (1.0_dblp - exp(db1)) * avefor &
                           / (avefor + a * db3 * exp(db2))
        end if
        if (cell%forshare_st > fmax) cell%forshare_st = fmax

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Potential desert fraction
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        cell%desshare_st = 0.0_dblp

        ! northern deserts
        if (cell%gdd0 < 100.0_dblp) then
          cell%desshare_st = 1.0_dblp
        else if (cell%gdd0 < gdd0_min) then
          cell%desshare_st = (gdd0_min - cell%gdd0) / (gdd0_min - 100.0_dblp)
        end if

        ! southern deserts
        if (cell%gdd0 >= gdd0_max) then
          pcr = acr * exp(gamm2 / 2.0_dblp * differ)
          if (cell%ave_pr05 <= pcr) then
            cell%desshare_st = 1.0_dblp
            cell%forshare_st = 0.0_dblp
          else
            db2 = (cell%ave_pr05 - pcr) / exp(gamm2 * differ)
            cell%desshare_st = 1.03_dblp / (1.0_dblp + ades * db2 * db2) - 0.03_dblp
            if (cell%desshare_st < 0.0_dblp) cell%desshare_st = 0.0_dblp
          end if
        end if

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  NPP — Lieth's formula
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        db1  = -v1 * cell%ave_pr
        db2  = -v2 * cell%ave_t
        npp1 = 1.0_dblp - exp(db1)
        npp2 = 1.0_dblp / (1.0_dblp + v3 * exp(db2))
        cell%npp = nppmax * min(npp1, npp2)

        ! CO2 enrichment (betat/betag already include the 1/log(2) factor)
        cell%nppt = cell%npp * (1.0_dblp + betat * log(co2ghg / 280.0_dblp))
        cell%nppg = cell%npp * (1.0_dblp + betag * log(co2ghg / 280.0_dblp))

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Allocation fractions and residence times (functions of NPP)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        k1t = c1t + c2t / (1.0_dblp + c3t * cell%nppt)
        k1g = c1g + c2g / (1.0_dblp + c3g * cell%nppg)

        t1t = d1t + d2t / (1.0_dblp + d3t * cell%nppt)
        t1g = d1g + d2g / (1.0_dblp + d3g * cell%nppg)

        t2t = e1t + e2t / (1.0_dblp + e3t * cell%nppt)
        t2g = e1g + e2g / (1.0_dblp + e3g * cell%nppg)

        ! soil decomposition — temperature-sensitive residence times
        t3t = 16.0_dblp  * exp(-ps5 * (cell%ave_t - soilt))
        t3g = 40.0_dblp  * exp(-ps5 * (cell%ave_t - soilt))
        t4t = 900.0_dblp * exp(-ps5 * (cell%ave_t - soilt))
        t4g = t4t

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Potential needleleaf fraction
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        cell%nlshare_st = (t1t - t1td) / (t1tn - t1td)
        cell%nlshare_st = min(1.0_dblp, max(0.0_dblp, cell%nlshare_st))

#if ( CYCC == 2 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  C4 grass fraction
!  [BUG] cell%tatmsmin is not yet connected to the atmospheric model; it was
!  hardcoded to 0.0 in the original code, which forces g4share_st = 0 always.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (cell%tatmsmin < 12.0_dblp) then
          cell%g4share_st = 0.0_dblp
        else if (cell%tatmsmin < 15.5_dblp) then
          cell%g4share_st = 1.0_dblp - (15.5_dblp - cell%tatmsmin) / 3.5_dblp
        else
          cell%g4share_st = 1.0_dblp
        end if
#endif

      end subroutine ccparam

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: climpar  (private)
!
!>     @brief  Compute LAI, area-weighted carbon stocks, annual uptake and NPP
!>             for the current cell.  Called at the end of ccstat and ccdyn.
!
!>     @param[inout] cell   vegetation state for the current grid cell
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine climpar(cell)

        type(veget_cell_state_t), intent(inout) :: cell

        real(dblp) :: prev_stock   ! total carbon stock at start of step
#if ( CYCC == 2 )
        real(dblp) :: prev_bc13
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  LAI
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        cell%laig  = cell%b1g * deng
        cell%lait  = cell%b1t * (dentn * cell%snlt + dentd * (1.0_dblp - cell%snlt))
        cell%blai(1) = cell%lait
        cell%blai(2) = cell%laig

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Save previous total carbon stock for annual uptake calculation
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        prev_stock = cell%b1 + cell%b2 + cell%b3 + cell%b4
#if ( CYCC == 2 )
        prev_bc13 = cell%bc13
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Area-weighted carbon pool totals
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        cell%b1 = cell%b1t * cell%st + cell%b1g * cell%sg
        cell%b2 = cell%b2t * cell%st + cell%b2g * cell%sg
        cell%b3 = cell%b3t * cell%st + cell%b3g * cell%sg
        cell%b4 = cell%b4t * cell%st + cell%b4g * cell%sg

        ! FROG_EXP: integrated vegetation carbon flux and leaf/wood ratio
        ! (originally guarded by #if FROG_EXP > 0; now unconditional)
        cell%Fv   = cell%Fv + cell%Fv_t * cell%st + cell%Fv_g * cell%sg
        if ((cell%b1 + cell%b2) > 0.0_dblp) then
          cell%r_leaf = (cell%b1 + cell%b2g * cell%sg) / (cell%b1 + cell%b2)
        else
          cell%r_leaf = 0.0_dblp
        end if

#if ( CYCC == 2 )
        cell%bc14 = ( cell%b1t14 + cell%b2t14 + cell%b3t14 + cell%b4t14 ) * cell%st &
                  + ( cell%b1g14 + cell%b2g14 + cell%b3g14 + cell%b4g14 ) * cell%sg

        cell%bc13 = ( cell%b1t13 + cell%b2t13 + cell%b3t13 + cell%b4t13 ) * cell%st &
                  + ( cell%b1g13 + cell%b2g13 + cell%b3g13 + cell%b4g13 ) * cell%sg
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Annual carbon uptake and stock diagnostics
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        cell%stock = cell%b1 + cell%b2 + cell%b3 + cell%b4
        cell%anup  = cell%stock - prev_stock

#if ( CYCC == 2 )
        cell%anup13 = cell%bc13 - prev_bc13
#endif

        cell%pnpp = cell%nppt * cell%st + cell%nppg * cell%sg

      end subroutine climpar

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      FUNCTION: ccstat
!
!>     @brief  Initialise the vegetation state of a cell to its equilibrium
!>             (static run or first call with kveget < 0).
!>             Sets carbon pools, cover fractions, and computes LAI and stocks.
!
!>     @param[inout] cell        vegetation state for the current grid cell
!>     @param[in]    irunlabelf  run label year offset (comrunlabel_mod)
!>     @param[in]    iyear       current model year
!>     @returns      .true. on success
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine ccstat(cell, irunlabelf, iyear)

        use veget_mod, only: iscendef, ivegstrt

        type(veget_cell_state_t), intent(inout) :: cell
        integer(ip),              intent(in)    :: irunlabelf, iyear

        real(dblp) :: tempor1

        call ccparam(cell)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Equilibrium carbon pools (per unit fractional area)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        cell%b1t = k1t * t1t * cell%nppt
        cell%b1g = k1g * t1g * cell%nppg

        cell%b2t = (1.0_dblp - k1t) * t2t * cell%nppt
        cell%b2g = (1.0_dblp - k1g) * t2g * cell%nppg

        cell%b3t = (k0t * cell%b1t / t1t + k2t / t2t * cell%b2t) * t3t
        cell%b3g = (k0g * cell%b1g / t1g + k2g / t2g * cell%b2g) * t3g

        cell%b4t = (k3t / t3t * cell%b3t) * t4t
        cell%b4g = (k4g / t2g * cell%b2g + k3g / t3g * cell%b3g) * t4g

        ! FROG_EXP: initialise flux tracking pools to their equilibrium values
        cell%Fv_t = cell%b4t
        cell%Fv_g = cell%b4g

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Equilibrium cover fractions
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! Land-use scenario: constrain tree fraction against crop fraction
        if ((iscendef == 1) .and. ((irunlabelf + iyear) >= ivegstrt)) then
          if (irunlabelf + iyear == ivegstrt) cell%st_const = cell%st
          tempor1 = 9.0e19_dblp
          if (cell%crop_fraction < 9.0e19_dblp) then
            tempor1 = cell%st_const - cell%crop_fraction
          end if
          cell%st = min(cell%forshare_st, tempor1)
        else
          cell%st = cell%forshare_st
        end if
        if (cell%st < 0.0_dblp) cell%st = 0.0_dblp

        cell%sd   = cell%desshare_st
        cell%snlt = cell%nlshare_st
        if (cell%sd < 0.0_dblp) cell%sd = 0.0_dblp
        cell%sg = 1.0_dblp - cell%st - cell%sd
        if (cell%sg < 0.0_dblp) cell%sg = 0.0_dblp

        call climpar(cell)


      end subroutine ccstat

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      FUNCTION: ccstatR
!
!>     @brief  Compute the equilibrium (reference) cover fractions for a cell
!>             without modifying the actual st/sg/sd state.  Used in parallel with
!>             ccstat when irad == 1 (radiative reference vegetation).
!
!>     @param[inout] cell   vegetation state for the current grid cell
!>     @returns      .true. on success
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine ccstatR(cell)

        type(veget_cell_state_t), intent(inout) :: cell


        call ccparam(cell)

        cell%stR   = cell%forshare_st
        cell%sdR   = cell%desshare_st
        cell%snltR = cell%nlshare_st
        if (cell%sdR < 0.0_dblp) cell%sdR = 0.0_dblp
        cell%sgR = 1.0_dblp - cell%stR - cell%sdR
        if (cell%sgR < 0.0_dblp) cell%sgR = 0.0_dblp

        call climpar(cell)


      end subroutine ccstatR

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      FUNCTION: ccdyn
!
!>     @brief  Advance the vegetation state of a cell by one year (dynamic run).
!>             Updates cover fractions, carbon pools, LAI and diagnostics.
!
!>     @param[inout] cell        vegetation state for the current grid cell
!>     @param[in]    irunlabelf  run label year offset (comrunlabel_mod)
!>     @param[in]    iyear       current model year
!>     @returns      .true. on success
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine ccdyn(cell, irunlabelf, iyear)

#if ( CYCC == 2 )
        ! Isotope pool arrays (B4T14, B3T13, etc.) are now fields of cell% —
        ! no import from veget_iso needed. c13frac/c13frac4 imported at module
        ! level from veget_phys_params_mod.
        use carbone_co2_mod,  only: C14ATM, C14DEC, PA_C
        use mod_sync_time,only: KENDY
        use C_res_mod,    only: c13atm
#endif
        use veget_mod, only: iscendef, ivegstrt

        type(veget_cell_state_t), intent(inout) :: cell
        integer(ip),              intent(in)    :: irunlabelf, iyear

        real(dblp) :: fd, dd, nld, dst, dsd, dstime
        real(dblp) :: temp_st, temp_sg
        real(dblp) :: tempor1, tempor2
        ! saved b4/b3 before correction (conservation law)
        real(dblp) :: b4t_prev, b3t_prev
#if ( CYCC == 2 )
        real(dblp) :: tempor3, tempor4, tempor5, tempor6
        real(dblp) :: G4D
#endif

        call ccparam(cell)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Cover fraction dynamics
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        fd  = cell%forshare_st - cell%st
        dd  = cell%desshare_st - cell%sd
        nld = cell%nlshare_st  - cell%snlt
#if ( CYCC == 2 )
        G4D = cell%g4share_st  - cell%sg4
#endif
        temp_st = cell%st
        temp_sg = cell%sg

        ! forest — exponential filter
        if (abs(fd) < 100.0_dblp * tiny(fd)) then
          dst = 0.0_dblp
        else
          dst = fd * (1.0_dblp - exp(-1.0_dblp / t2t))
        end if
        cell%st = temp_st + dst
        if (cell%st < 0.0_dblp) cell%st = 0.0_dblp

        ! desert — exponential filter with characteristic time adjustment
        dsd     = dd * (1.0_dblp - exp(-1.0_dblp / t2g))
        tempor1 = cell%sd + dsd + cell%st
        if (tempor1 > 0.9_dblp) then
          dstime = (t2g * (1.0_dblp - tempor1) + t2t * (tempor1 - 0.9_dblp)) * 10.0_dblp
          dsd    = dd * (1.0_dblp - exp(-1.0_dblp / dstime))
        end if
        cell%sd = cell%sd + dsd
        if (cell%sd < 0.0_dblp) cell%sd = 0.0_dblp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Constant-vegetation scenario (iscendef == -1)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if ((iscendef == -1) .and. (irunlabelf + iyear - ivegstrt >= 0)) then
          if (irunlabelf + iyear == ivegstrt) then
            cell%st_const = cell%st
            cell%sd_const = cell%sd
          else
            cell%st = cell%st_const
            cell%sd = cell%sd_const
          end if
        end if

        tempor2 = 1.0_dblp - cell%sd
        if (tempor2 < 0.0_dblp) tempor2 = 0.0_dblp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Land-use scenario (iscendef == +1) — scenario version B (Brovkin / EMIC protocol)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if ((iscendef == 1) .and. (irunlabelf + iyear - ivegstrt >= 0)) then
          if (irunlabelf + iyear == ivegstrt) then
            cell%st_const = cell%st
            cell%sd_const = cell%sd
          else
            cell%sd   = cell%sd_const
            tempor2   = 1.0_dblp - cell%sd_const
            cell%st   = cell%st_const
            if (cell%crop_fraction < 9.0e19_dblp) then
              cell%st = cell%st_const - cell%crop_fraction
            end if
            if (cell%st < 0.0_dblp) cell%st = 0.0_dblp
          end if
        end if

        ! update dst for carbon conservation correction below
        dst = cell%st - temp_st

        cell%sg = tempor2 - cell%st
        if (cell%sg < 0.0_dblp) cell%sg = 0.0_dblp

#if ( CYCC == 2 )
        cell%sg4  = cell%g4share_st - G4D * exp(-1.0_dblp / t2g)
#endif
        cell%snlt = cell%nlshare_st - nld * exp(-1.0_dblp / t2t)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Carbon pool conservation — correct b3/b4 for fractional area changes
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! save b4t / b3t before mixing (needed for grass correction below)
        b4t_prev = cell%b4t
        b3t_prev = cell%b3t
#if ( CYCC == 2 )
        tempor3 = cell%b4t14 ;  tempor4 = cell%b3t14
        tempor5 = cell%b4t13 ;  tempor6 = cell%b3t13
#endif

        ! trees advancing into grass area
        if (cell%st > 0.0_dblp) then
          if (dst > 0.0_dblp) then
            if (dst <= temp_sg) then
              cell%b4t = (cell%b4t * temp_st + cell%b4g * dst) / cell%st
              cell%b3t = (cell%b3t * temp_st + cell%b3g * dst) / cell%st
#if ( CYCC == 2 )
              cell%b4t14 = (cell%b4t14 * temp_st + cell%b4g14 * dst) / cell%st
              cell%b3t14 = (cell%b3t14 * temp_st + cell%b3g14 * dst) / cell%st
              cell%b4t13 = (cell%b4t13 * temp_st + cell%b4g13 * dst) / cell%st
              cell%b3t13 = (cell%b3t13 * temp_st + cell%b3g13 * dst) / cell%st
#endif
            else
              cell%b4t = (cell%b4t * temp_st + cell%b4g * temp_sg) / cell%st
              cell%b3t = (cell%b3t * temp_st + cell%b3g * temp_sg) / cell%st
#if ( CYCC == 2 )
              cell%b4t14 = (cell%b4t14 * temp_st + cell%b4g14 * temp_sg) / cell%st
              cell%b3t14 = (cell%b3t14 * temp_st + cell%b3g14 * temp_sg) / cell%st
              cell%b4t13 = (cell%b4t13 * temp_st + cell%b4g13 * temp_sg) / cell%st
              cell%b3t13 = (cell%b3t13 * temp_st + cell%b3g13 * temp_sg) / cell%st
#endif
            end if
          end if
          cell%b2t = cell%b2t * temp_st / cell%st
          cell%b1t = cell%b1t * temp_st / cell%st
#if ( CYCC == 2 )
          cell%b2t14 = cell%b2t14 * temp_st / cell%st
          cell%b1t14 = cell%b1t14 * temp_st / cell%st
          cell%b2t13 = cell%b2t13 * temp_st / cell%st
          cell%b1t13 = cell%b1t13 * temp_st / cell%st
#endif
        end if

        ! grass pools — corrected for tree advance or retreat
        if (cell%sg > 0.0_dblp) then
          if (dst > 0.0_dblp) then
            if (dst <= temp_sg) then
              cell%b4g = cell%b4g * (temp_sg - dst) / cell%sg
              cell%b3g = cell%b3g * (temp_sg - dst) / cell%sg
#if ( CYCC == 2 )
              cell%b4g14 = cell%b4g14 * (temp_sg - dst) / cell%sg
              cell%b3g14 = cell%b3g14 * (temp_sg - dst) / cell%sg
              cell%b4g13 = cell%b4g13 * (temp_sg - dst) / cell%sg
              cell%b3g13 = cell%b3g13 * (temp_sg - dst) / cell%sg
#endif
            else
              cell%b4g = 0.0_dblp ;  cell%b3g = 0.0_dblp
#if ( CYCC == 2 )
              cell%b4g14 = 0.0_dblp ;  cell%b3g14 = 0.0_dblp
              cell%b4g13 = 0.0_dblp ;  cell%b3g13 = 0.0_dblp
#endif
            end if
          else
            ! trees retreating — grass pools gain tree organic matter
            cell%b4g = (cell%b4g * temp_sg - b4t_prev * dst) / cell%sg
            cell%b3g = (cell%b3g * temp_sg - b3t_prev * dst) / cell%sg
#if ( CYCC == 2 )
            cell%b4g14 = (cell%b4g14 * temp_sg - tempor3 * dst) / cell%sg
            cell%b3g14 = (cell%b3g14 * temp_sg - tempor4 * dst) / cell%sg
            cell%b4g13 = (cell%b4g13 * temp_sg - tempor5 * dst) / cell%sg
            cell%b3g13 = (cell%b3g13 * temp_sg - tempor6 * dst) / cell%sg
#endif
          end if
          cell%b2g = cell%b2g * temp_sg / cell%sg
          cell%b1g = cell%b1g * temp_sg / cell%sg
#if ( CYCC == 2 )
          cell%b2g14 = cell%b2g14 * temp_sg / cell%sg
          cell%b1g14 = cell%b1g14 * temp_sg / cell%sg
          cell%b2g13 = cell%b2g13 * temp_sg / cell%sg
          cell%b1g13 = cell%b1g13 * temp_sg / cell%sg
#endif
        end if

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Slow soil organic matter (b4)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        cell%b4t = cell%b4t + k3t / t3t * cell%b3t - cell%b4t / t4t
        cell%b4g = cell%b4g + k4g / t2g * cell%b2g + k3g / t3g * cell%b3g - cell%b4g / t4g

        ! FROG_EXP: fast carbon flux into slow pool (input to soil)
        cell%Fv_t = k3t / t3t * cell%b3t
        cell%Fv_g = k4g / t2g * cell%b2g + k3g / t3g * cell%b3g

#if ( CYCC == 2 )
        cell%b4t14 = cell%b4t14 + k3t / t3t * cell%b3t14 - cell%b4t14 / t4t
        cell%b4g14 = cell%b4g14 + k4g / t2g * cell%b2g14 + k3g / t3g * cell%b3g14 - cell%b4g14 / t4g
        cell%b4t13 = cell%b4t13 + k3t / t3t * cell%b3t13 - cell%b4t13 / t4t
        cell%b4g13 = cell%b4g13 + k4g / t2g * cell%b2g13 + k3g / t3g * cell%b3g13 - cell%b4g13 / t4g
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Fast soil organic matter (b3)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        cell%b3t = cell%b3t + cell%b1t / t1t * k0t + k2t / t2t * cell%b2t - cell%b3t / t3t
        cell%b3g = cell%b3g + cell%b1g / t1g * k0g + k2g / t2g * cell%b2g - cell%b3g / t3g

#if ( CYCC == 2 )
        cell%b3t14 = cell%b3t14 + cell%b1t14 / t1t * k0t + k2t / t2t * cell%b2t14 - cell%b3t14 / t3t
        cell%b3g14 = cell%b3g14 + cell%b1g14 / t1g * k0g + k2g / t2g * cell%b2g14 - cell%b3g14 / t3g
        cell%b3t13 = cell%b3t13 + cell%b1t13 / t1t * k0t + k2t / t2t * cell%b2t13 - cell%b3t13 / t3t
        cell%b3g13 = cell%b3g13 + cell%b1g13 / t1g * k0g + k2g / t2g * cell%b2g13 - cell%b3g13 / t3g
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Leaves biomass (b1)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        cell%b1t = cell%b1t + k1t * cell%nppt - cell%b1t / t1t
        cell%b1g = k1g * cell%nppg * t1g     ! grass leaves: instantaneous equilibrium

#if ( CYCC == 2 )
        cell%b1t14 = cell%b1t14 + k1t * cell%nppt * (C14ATM / PA_C) - cell%b1t14 / t1t
        cell%b1g14 = k1g * cell%nppg * (C14ATM / PA_C) * t1g

        cell%b1t13 = cell%b1t13 + k1t * cell%nppt * c13atm * C13FRAC - cell%b1t13 / t1t
        cell%b1g13 = k1g * cell%nppg * c13atm &
                   * (C13FRAC * (1.0_dblp - cell%sg4) + C13FRAC4 * cell%sg4) * t1g
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Stems and roots (b2)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        cell%b2t = cell%b2t + (1.0_dblp - k1t) * cell%nppt - cell%b2t / t2t
        cell%b2g = cell%b2g + (1.0_dblp - k1g) * cell%nppg - cell%b2g / t2g

#if ( CYCC == 2 )
        cell%b2t14 = cell%b2t14 + (1.0_dblp - k1t) * cell%nppt * (C14ATM / PA_C) - cell%b2t14 / t2t
        cell%b2g14 = cell%b2g14 + (1.0_dblp - k1g) * cell%nppg * (C14ATM / PA_C) - cell%b2g14 / t2g

        cell%b2t13 = cell%b2t13 + (1.0_dblp - k1t) * cell%nppt * c13atm * C13FRAC - cell%b2t13 / t2t
        cell%b2g13 = cell%b2g13 &
                   + (1.0_dblp - k1g) * cell%nppg &
                   * c13atm * (C13FRAC * (1.0_dblp - cell%sg4) + C13FRAC4 * cell%sg4) &
                   - cell%b2g13 / t2g

        ! annual 14C radioactive decay
        if (KENDY == 1) then
          cell%b1t14 = cell%b1t14 * (1.0_dblp - C14DEC)
          cell%b2t14 = cell%b2t14 * (1.0_dblp - C14DEC)
          cell%b3t14 = cell%b3t14 * (1.0_dblp - C14DEC)
          cell%b4t14 = cell%b4t14 * (1.0_dblp - C14DEC)
          cell%b1g14 = cell%b1g14 * (1.0_dblp - C14DEC)
          cell%b2g14 = cell%b2g14 * (1.0_dblp - C14DEC)
          cell%b3g14 = cell%b3g14 * (1.0_dblp - C14DEC)
          cell%b4g14 = cell%b4g14 * (1.0_dblp - C14DEC)
        end if
#endif

        call climpar(cell)


      end subroutine ccdyn

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      FUNCTION: ccdynR
!
!>     @brief  Advance the reference vegetation fractions (stR/sgR/sdR) by one year
!>             without touching the actual carbon pools or cover fractions.
!>             Used in parallel with ccdyn when irad == 1.
!
!>     @param[inout] cell   vegetation state for the current grid cell
!>     @returns      .true. on success
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine ccdynR(cell)

        type(veget_cell_state_t), intent(inout) :: cell

        real(dblp) :: fd, dd, nld, dst, dsd, dstime, tempor1
        real(dblp) :: temp_st, temp_sg

        call ccparam(cell)

        fd  = cell%forshare_st - cell%stR
        dd  = cell%desshare_st - cell%sdR
        nld = cell%nlshare_st  - cell%snltR
        temp_st = cell%stR
        temp_sg = cell%sgR

        ! forest — exponential filter
        if (abs(fd) < 100.0_dblp * tiny(fd)) then
          dst = 0.0_dblp
        else
          dst = fd * (1.0_dblp - exp(-1.0_dblp / t2t))
        end if
        cell%stR = cell%stR + dst
        if (cell%stR < 0.0_dblp) cell%stR = 0.0_dblp

        cell%snltR = cell%nlshare_st - nld * exp(-1.0_dblp / t2t)

        ! desert — exponential filter with characteristic time adjustment
        dsd     = cell%desshare_st - dd * exp(-1.0_dblp / t2g) - cell%sdR
        tempor1 = cell%sdR + dsd + cell%stR
        if (tempor1 > 0.9_dblp) then
          dstime = t2g * (1.0_dblp - tempor1) * 10.0_dblp &
                 + t2t * (tempor1 - 0.9_dblp) * 10.0_dblp
          dsd = cell%desshare_st - dd * exp(-1.0_dblp / dstime) - cell%sdR
        end if
        cell%sdR = cell%sdR + dsd
        if (cell%sdR < 0.0_dblp) cell%sdR = 0.0_dblp

        cell%sgR = 1.0_dblp - cell%stR - cell%sdR
        if (cell%sgR < 0.0_dblp) cell%sgR = 0.0_dblp


      end subroutine ccdynR

#if ( CYCC == 2 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: ccstat_isotope
!
!>     @author  Didier M. Roche (dmr), Veronique Mariotti (vm)
!
!>     @brief   Initialise carbon isotope pools (13C, 14C) for a single grid cell
!>              when no isotope restart is available (new_run_c /= 0).
!>              Called once per cell from the restart initialisation loop in veget.f90,
!>              immediately after ccstat has computed the equilibrium carbon pools and
!>              ccparam has set the residence times t1t..t4g.
!>
!>     @details In the original code (ccstat_isotope.f) this routine operated on the
!>              full 2-D arrays using the last-computed scalar residence times — an
!>              implicit bug since t1t..t4g held the values of the last cell visited.
!>              The point-wise port eliminates that race: ccparam is called by ccstat
!>              before this routine runs, so the cell's own residence times are current.
!>
!>     @param[inout] cell   vegetation state for the current grid cell
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine ccstat_isotope(cell)

        ! c13frac, c13frac4 now in veget_phys_params_mod (imported at module level)
        use carbone_co2_mod, only: c14rstd, c14dec
        use C_res_mod,   only: c13init, c13atm

        type(veget_cell_state_t), intent(inout) :: cell

        real(dblp) :: c14init

        c14init = c14rstd   ! standard modern 14C activity ratio

        ! leaves
        cell%b1t14 = c14init * exp(-c14dec * t1t) * cell%b1t
        cell%b1g14 = c14init * exp(-c14dec * t1g) * cell%b1g
        cell%b1t13 = cell%b1t * c13init * c13atm
        cell%b1g13 = cell%b1g * c13init * c13atm

        ! stems and roots
        cell%b2t14 = cell%b2t * c14init * exp(-c14dec * t2t)
        cell%b2g14 = cell%b2g * c14init * exp(-c14dec * t2g)
        cell%b2t13 = cell%b2t * c13init * c13atm
        cell%b2g13 = cell%b2g * c13init * c13atm

        ! fast litter
        cell%b3t14 = cell%b3t * c14init * exp(-c14dec * t3t)
        cell%b3g14 = cell%b3g * c14init * exp(-c14dec * t3g)
        cell%b3t13 = cell%b3t * c13init * c13atm
        cell%b3g13 = cell%b3g * c13init * c13atm

        ! slow soil organic matter
        cell%b4t14 = cell%b4t * c14init * exp(-c14dec * t4t)
        cell%b4g14 = cell%b4g * c14init * exp(-c14dec * t4g)
        cell%b4t13 = cell%b4t * c13init * c13atm
        cell%b4g13 = cell%b4g * c13init * c13atm

      end subroutine ccstat_isotope

#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of module
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      end module veget_sub_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
