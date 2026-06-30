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
!      SUBROUTINE: veget
!
!>     @author  V. Brovkin, A. Ganopolski, dmr (cleanup 2019), refactoring by dmr, clo
!
!>     @brief   Terrestrial Vegetation annual Model (VECODE) — orchestration layer.
!
!>     @details Called every atmospheric time step by the host model.  On the first
!>              call (iyear == 0) it reads parameter and restart files and initialises
!>              all state.  On subsequent calls it accumulates climate forcing, then
!>              once per vegetation year it:
!>                (1) resolves per-cell climate means,
!>                (2) drives the point-wise physics (ccstat / ccdyn) cell by cell,
!>                (3) writes albedo and land-surface fields back to the atmosphere,
!>                (4) calls ASCII and NetCDF output routines,
!>                (5) writes restart files and resets annual accumulators.
!>
!>              All physics is delegated to veget_sub_mod.  No carbon pool arithmetic
!>              or vegetation dynamics appear here.
!>
!>     @param[in] ist      sub-step index within the atmospheric year (1..nstpyear)
!>     @param[in] jst      unused legacy argument (kept for call-site compatibility)
!>     @param[in] dtime    length of one atmospheric time step                 [s]
!>     @param[in] epss     land-fraction threshold below which cell is skipped [–]
!>     @param[in] patmCO2  atmospheric CO2 concentration                       [ppm]
!>     @param[in] fracgr   land fraction of each grid cell                     [–]
!>     @param[in] darea    area of each latitude band                          [m^2]
!>     @param[in] tempgr   surface temperature field                           [K]
!
!>     @date  Original : V. Brovkin / A. Ganopolski, last mod. 17.11.97
!>     @date  Cleanup  : dmr, 2019-10-17
!>     @date  Refactored to F90, point-wise physics : 2026-06-30 (dmr, clo)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine veget(ist, jst, dtime, epss, patmCO2, fracgr, darea, tempgr)

        use global_constants_mod,  only: dblp=>dp, ip

#if ( CYCC == 2 )
        use carbone_co2_mod,           only: new_run_c
#if ( KC14 == 1 )
        use C_res_mod,             only: cav_la14_b, cav_la_b
#endif
#endif

#if ( ISOATM >= 1 )
        use iso_param_mod
#endif

        use comland_mod,           only: bmoisg, bmoism, albland, forestfr
        use comatm,                only: nlat, nlon, iwater
        use comphys,               only: tosnow, torain
        use comdiag,               only: irad
        use comemic_mod,           only: iyear, nyears, flgveg, iatm,   &
                                         nstpyear, nwrskip
        use comrunlabel_mod,       only: irunlabelf

        use veget_mod
        use veget_phys_params_mod, only: co2ghg, betat, betag, gamm2
        use veget_cell_state_mod,  only: veget_cell_state_t, get_crop_fraction
        use veget_sub_mod,         only: initcpar, ccstat, ccstatR,     &
                                         ccdyn, ccdynR
#if ( CYCC == 2 )
        use veget_sub_mod,         only: ccstat_isotope, ccparam
#endif

#if ( VEG_LUH == 1 )
        use luh_mod,               only: read_luh2   ! external LUH2 reader
#endif

        implicit none

        include 'netcdf.inc'

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Arguments
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer,     intent(in) :: ist, jst
        real(dblp),  intent(in) :: dtime
        real(dblp),  intent(in) :: epss
        real(dblp),  intent(in) :: patmCO2
        real(dblp),  intent(in) :: fracgr(nlat,nlon)
        real(dblp),  intent(in) :: darea(nlat)
        real(dblp),  intent(in) :: tempgr(nlat,nlon)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Local variables — saved across calls (persistent run-control state)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        logical,      save :: flgwrt = .true.
        integer(ip),  save :: kveget, nyvegt, nwrveg, ncumvg, iyr0vg
        real(dblp),   save :: tmxdry, dtrdry, rpfdry, unsdry, prcmin, bmtdry
        integer(ip),  save :: ns0inv(4,0:1)
        character(len=15), save :: titveg

        ! albedo lookup tables (seasonal, read from veget.par)
        real(dblp), save :: albet(4), albeg(4), albed(4), albegc(4), albedc(4)

        ! climate accumulation arrays — annual means built up sub-step by sub-step
        real(dblp), dimension(nlat,nlon), save :: temveg  = 0.0_dblp
        real(dblp), dimension(nlat,nlon), save :: gd0veg  = 0.0_dblp
        real(dblp), dimension(nlat,nlon), save :: gd5veg  = 0.0_dblp
        real(dblp), dimension(nlat,nlon), save :: tpsdry  = 0.0_dblp
        real(dblp), dimension(nlat,nlon), save :: freezeI = 0.0_dblp
        real(dblp), dimension(nlat,nlon), save :: thawI   = 0.0_dblp
        real(dblp), dimension(nlat,nlon), save :: prcday  = 0.0_dblp
        real(dblp), dimension(nlat,nlon,2), save :: prcveg = 0.0_dblp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Local variables — temporaries
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip) :: i, j, k, ns, ii, kk2, nlatd2, ieq, iyrloc

        integer(ip) :: istep
        real(dblp)  :: ttemp, prcrit, ddc, xx, xxalb
        real(dblp)  :: zmoy1, zmoy2, zmoy3, zmoy4
        real(dblp)  :: bmtd00
        real(dblp)  :: ari(nlat/2), rj, rw, dumwei
        integer(ip) :: ilat_g, ird
        character(len=15) :: titv00
        real(dblp) :: fareal(nlat,nlon)

        ! per-cell state (allocated on stack, reused every cell)
        type(veget_cell_state_t) :: cell

        ! file unit integers
        integer :: vegetpar_id, gaussasc_id
        integer :: vegetdat_id, vegetrest_id, vegetinit_id

        ! format for gauss.asc
220     format(f18.10,f17.10)

        !NOTE moved here upwards since co2ghg needs to be initialized at all times
        co2ghg = patmCO2

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Section 0 — First call (iyear == 0): read parameters and initialise state
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (iyear == 0) then

          !--- read veget.par ---
          open(newunit=vegetpar_id, file='veget.par', status='old')
          read(vegetpar_id,*)
          read(vegetpar_id,*)
          read(vegetpar_id,'(A)') titveg
          read(vegetpar_id,*)
          read(vegetpar_id,*) kveget, nyvegt, nwrveg
          read(vegetpar_id,*)
          read(vegetpar_id,*) prcmin, bmtdry, tmxdry, dtrdry, rpfdry
          read(vegetpar_id,*)
          read(vegetpar_id,*) ieqveg, ieqvegwr
          read(vegetpar_id,*)
          read(vegetpar_id,*) iscendef, ivegstrt
          read(vegetpar_id,*)
          read(vegetpar_id,*) (albet(i),  i=1,4)
          read(vegetpar_id,*)
          read(vegetpar_id,*) (albeg(i),  i=1,4)
          read(vegetpar_id,*)
          read(vegetpar_id,*) (albed(i),  i=1,4)
          read(vegetpar_id,*)
          read(vegetpar_id,*) (albegc(i), i=1,4)
          read(vegetpar_id,*)
          read(vegetpar_id,*) (albedc(i), i=1,4)
          read(vegetpar_id,*)
          read(vegetpar_id,*) gamm2
          read(vegetpar_id,*)
          read(vegetpar_id,*) betag, betat
          read(vegetpar_id,*)
          read(vegetpar_id,*) fco2veg
          close(vegetpar_id)

          ! CO2 enrichment factors include 1/log(2) once and for all
          betat = betat / log(2.0_dblp)
          betag = betag / log(2.0_dblp)

          !--- Gaussian latitudes (for NetCDF coordinate output) ---
          open(newunit=gaussasc_id, file='inputdata/gauss.asc', &
               status='old', form='formatted')
          rewind(gaussasc_id)
          ilat_g = nlat / 2
10        continue
            read(gaussasc_id, 220, end=15) rj, rw
            ird = int(rj)
            if (ird == ilat_g) then
              do i = 1, ird
                read(gaussasc_id, 220) ari(i), dumwei
              end do
              goto 20
            else
              goto 10
            end if
15        continue
20        continue
          do i = 1, ilat_g
            phi(i)        = -ari(ilat_g + 1 - i)
            phi(ilat_g+i) =  ari(i)
          end do
          phi(1:nlat) = asin(phi(1:nlat))
          close(gaussasc_id)

          !--- NetCDF output variable metadata (outp_veget.param) ---
          call read_veget_output_param()

          !--- deforestation scenario (VEGET.dat or VEGET.nc) ---
          i0dfor = 0
          ndfor  = 0
          if (iscendef == 1) call read_deforestation_scenario()
          if (iscendef == -1) i0dfor = ivegstrt

          !--- dry-soil parameter ---
          unsdry = 0.0_dblp
          if (dtrdry > 1.0e-6_dblp) unsdry = 1.0_dblp / dtrdry

          !--- early exit if vegetation is disabled ---
          if (nyvegt == 0) then
            nyvegt = nyears + 1
            kveget = abs(kveget)
            flgveg = .false.
            return
          end if
          nwrveg = (nwrveg / nyvegt) * nyvegt
          if (nwrveg == 0) then
            nwrveg  = nyears + 1
            flgwrt  = .false.
          end if
          iyr0vg = 0

          !--- season index permutation (North/South hemisphere swap) ---
          do ns = 1, 4
            ns0inv(ns,0) = 1 + mod(ns+1, 4)
            ns0inv(ns,1) = ns
          end do

          !--- initialise physical parameters ---
          if (kveget /= 0) call initcpar()

          !--- Section 1: initial vegetation state ---

          if (kveget > 0) then
            call read_restart(vegetrest_id, iyr0vg, bmtd00, titv00)
            ! kveget>0 : vegetation read from restart; fall through to section 4
            ! (albedo) like kveget<0. Sections 2, 2b and 3 are skipped via the
            ! iyear/=0 guard and the kveget<0/.or.iyear*kveget/=0 guard.
          else if (kveget < 0) then
            call read_veget_init(vegetinit_id, iyrloc, bmtd00, titv00)
            ! kveget<0 : climate accumulators loaded; fall through to section 3+4.
          end if

        end if  ! iyear == 0
        ! NOTE: when iyear==0 (both kveget>0 and kveget<0) we fall through here,
        ! skipping sections 2 and 2b. Section 3 is additionally skipped for kveget>0
        ! (restart read) via its own guard (kveget<0 .or. iyear*kveget/=0).
        ! Section 4 (albedo) is always executed, initialising forestfr and albland.
        if (iyear /= 0) then

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Section 2 — Accumulate climate forcing (every atmospheric sub-step)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        istep  = (ist - 1) * iatm + jst
        ncumvg = ncumvg + 1

        do j = 1, nlon
          do i = 1, nlat
            ttemp        = tempgr(i,j) - 273.15_dblp
            temveg(i,j)  = temveg(i,j)  + ttemp
            prcday(i,j)  = prcday(i,j)  + dtime * (torain(i,j,iwater) &
                                                   + tosnow(i,j,iwater))
            gd0veg(i,j)  = gd0veg(i,j)  + max(ttemp,      0.0_dblp)
            gd5veg(i,j)  = gd5veg(i,j)  + max(ttemp-5.0_dblp, 0.0_dblp)
            freezeI(i,j) = freezeI(i,j) + abs(min(ttemp, 0.0_dblp))
            thawI(i,j)   = thawI(i,j)   + max(ttemp,      0.0_dblp)
          end do
        end do

        ! end of each model day: accumulate precipitation and dry-soil days
        if (mod(istep, iatm) == 0) then
          prcrit = abs(prcmin)
          do j = 1, nlon
            do i = 1, nlat
              prcveg(i,j,1) = prcveg(i,j,1) + prcday(i,j)
              if (prcday(i,j) >= prcrit) &
                prcveg(i,j,2) = prcveg(i,j,2) + prcday(i,j)
              prcday(i,j) = 0.0_dblp
              if (bmoisg(i,j,iwater) <= bmtdry) &
                tpsdry(i,j) = tpsdry(i,j) + 1.0_dblp
            end do
          end do
        end if

        ! only proceed to physics once per vegetation year
        if (mod(istep, nstpyear*nyvegt) /= 0) return

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Section 2b — Compute annual climate means from accumulators
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        zmoy1 = 0.0_dblp
        if (ncumvg >= 1) zmoy1 = 1.0_dblp / real(ncumvg, dblp)
        zmoy2 = zmoy1 * 360.0_dblp                    ! °C·day → °C·yr^-1
        !dmr The line below is not giving identical results to the one after
        ! zmoy3 = zmoy2 * real(iatm, dblp) * 1000.0_dblp ! m → mm yr^-1
        zmoy3 = zmoy2 * real(iatm*1000, dblp)  ! m → mm yr^-1
        zmoy4 = zmoy2 * real(iatm, dblp)              ! steps → days yr^-1

        kk2 = 2
        if (prcmin < 0.0_dblp) kk2 = 1

        do j = 1, nlon
          do i = 1, nlat
            temveg(i,j)   = temveg(i,j)   * zmoy1
            gd0veg(i,j)   = gd0veg(i,j)   * zmoy2
            gd5veg(i,j)   = gd5veg(i,j)   * zmoy2
            freezeI(i,j)  = freezeI(i,j)  * zmoy2
            thawI(i,j)    = thawI(i,j)    * zmoy1
            prcveg(i,j,1) = prcveg(i,j,1) * zmoy3
            prcveg(i,j,2) = prcveg(i,j,2) * zmoy3
            tpsdry(i,j)   = tpsdry(i,j)   * zmoy4
          end do
        end do

        ! apply dry-soil precipitation reduction
        do j = 1, nlon
          do i = 1, nlat
            prcveg(i,j,2) = prcveg(i,j,kk2) * &
              (1.0_dblp - rpfdry * min(1.0_dblp, &
               max(0.0_dblp, (tpsdry(i,j) - tmxdry) * unsdry)))
          end do
        end do

        end if  ! iyear /= 0  (sections 2 and 2b skipped on first call with kveget<0)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Section 3 — Point-wise vegetation physics
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        co2ghg = patmCO2

        do ieq = 1, ieqveg

          if (kveget < 0 .or. iyear * kveget /= 0) then

            do j = 1, nlon
              do i = 1, nlat

                if (fracgr(i,j) <= epss) cycle

                !--- extract global 2-D arrays into per-cell scalar state ---
                call cell_from_grid(cell, i, j, temveg, gd0veg, &
                                    prcveg, freezeI, thawI)

                !--- resolve crop fraction for this cell and year ---
                cell%crop_fraction = get_crop_fraction(              &
                                       farea, i, j,                  &
                                       irunlabelf, iyear,            &
                                       i0dfor, ndfor,                &
                                       iscendef, ivegstrt)

                !--- physics ---
                if (kveget < 0) then
                  if (irad == 1) call ccstatR(cell)
                  call ccstat(cell, irunlabelf, iyear)
                else
                  do k = 1, kveget
                    if (irad == 1) call ccdynR(cell)
                    call ccdyn(cell, irunlabelf, iyear)
                  end do
                end if

                !--- write updated cell state back to global 2-D arrays ---
                call cell_to_grid(cell, i, j)

              end do
            end do

          end if

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Section 4 — Albedo and land-surface fields
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          if (kveget /= 0) then

#if ( VEG_LUH == 1 )
            ! override computed fractions with prescribed LUH2 land use
            call apply_luh2_override()
#endif
            ! underflow guard
            where (sd < 10.0_dblp * tiny(sd(1,1))) sd = 0.0_dblp
            where (sg < 10.0_dblp * tiny(sg(1,1))) sg = 0.0_dblp
            where (st < 10.0_dblp * tiny(st(1,1))) st = 0.0_dblp
            where (sdR < 10.0_dblp * tiny(sdR(1,1))) sdR = 0.0_dblp
            where (sgR < 10.0_dblp * tiny(sgR(1,1))) sgR = 0.0_dblp
            where (stR < 10.0_dblp * tiny(stR(1,1))) stR = 0.0_dblp

            ! initialise bmoism
            bmoism(:,:) = 0.15_dblp

            nlatd2 = nlat / 2
            do j = 1, nlon
              do i = 6, nlat   ! skip Antarctic ice (rows 1-5)

                if (fracgr(i,j) <= epss) cycle

                ! hemisphere flag (ii=0 southern, ii=1 northern)
                ii = i / nlatd2
                if (ii == 2) ii = 1

                ! transition to tundra/steppe for near-freezing annual T
                xx = (4.0_dblp - temveg(i,j)) * 0.25_dblp
                xx = min(1.0_dblp, max(0.0_dblp, xx))

                ! bright Sahara sand desert flag
                ddc = 0.0_dblp
                if ((i >= 19 .and. i <= 23) .and. &
                    (j <= 12 .or. j >= 62)) ddc = 1.0_dblp
                if (i == 23 .and. j >= 5 .and. j <= 12) ddc = 0.5_dblp

                do ns = 1, 4
                  xxalb = st(i,j) * albet(ns)               &
                        + sd(i,j) * (albed(ns) + ddc*albedc(ns)) &
                        + sg(i,j) * (albeg(ns) + xx*albegc(ns))

                  ! Equatorial rain-forest low-albedo correction
                  if (i == 16 .or. i == 17) xxalb = 0.09_dblp
                  if (i == 15 .or. i == 18) xxalb = 0.5_dblp*xxalb + 0.045_dblp

                  albland(i,j,ns0inv(ns,ii)) = xxalb
                end do

                bmoism(i,j) = (st(i,j) * (1.0_dblp - snlt(i,j))) * 0.25_dblp &
                            +  sd(i,j)                             * 0.10_dblp &
                            +  sg(i,j)                             * 0.15_dblp &
                            + (st(i,j) *               snlt(i,j)) * 0.25_dblp

                forestfr(i,j) = st(i,j)

              end do
            end do

          end if  ! kveget /= 0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Section 5 — ASCII output
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          if (ieqveg /= 1) then
            flgwrt = (mod(ieq, ieqvegwr) == 0)
          end if

          if (flgwrt) call veget_wr(kveget, nyvegt, nwrveg, iyr0vg,    &
                                     iyear, nyears,                      &
                                     temveg, gd0veg, prcveg,            &
                                     prcmin, epss,                       &
                                     fracgr, darea, gd5veg, titveg)

        end do   ! ieq

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Section 6 — Restart write and annual accumulator reset
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if ((iyear /= 0) .and. (mod(istep, nstpyear) == 0)) then
          if (mod(iyear, nwrskip) == 0 .or. iyear == nyears) then
            call write_veget_init(vegetinit_id, iyr0vg, bmtdry, titveg)
            if (kveget /= 0) &
              call write_restart(vegetrest_id, iyr0vg, bmtdry, titveg)
          end if
        end if

        ! reset annual accumulators
        ncumvg = 0
        temveg  = 0.0_dblp
        gd0veg  = 0.0_dblp
        gd5veg  = 0.0_dblp
        prcveg  = 0.0_dblp
        tpsdry  = 0.0_dblp
        freezeI = 0.0_dblp
        thawI   = 0.0_dblp

        ! after first vegetation computation kveget becomes positive
        kveget = abs(kveget)

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  INTERNAL: cell_from_grid
!>  Extract the 2-D global arrays into a veget_cell_state_t scalar for cell (i,j).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine cell_from_grid(c, i, j, temv, gd0v, prcv, frz, thw)

        type(veget_cell_state_t), intent(out) :: c
        integer(ip),              intent(in)  :: i, j
        real(dblp), intent(in) :: temv(nlat,nlon), gd0v(nlat,nlon)
        real(dblp), intent(in) :: prcv(nlat,nlon,2)
        real(dblp), intent(in) :: frz(nlat,nlon), thw(nlat,nlon)

        ! climate forcing
        c%ave_t    = temv(i,j)
        c%gdd0     = gd0v(i,j)
        c%ave_pr   = prcv(i,j,1)
        c%ave_pr05 = prcv(i,j,2)
        if ((frz(i,j) + thw(i,j)) > 0.0_dblp) then
          c%fr_ndx = sqrt(frz(i,j)) / (sqrt(frz(i,j)) + sqrt(thw(i,j)))
        else
          c%fr_ndx = 0.0_dblp
        end if

        ! cover fractions
        c%st   = st(i,j)  ;  c%sg   = sg(i,j)  ;  c%sd   = sd(i,j)
        c%snlt = snlt(i,j)
        c%stR  = stR(i,j) ;  c%sgR  = sgR(i,j) ;  c%sdR  = sdR(i,j)
        c%snltR   = snltR(i,j)
        c%st_const = st_const(i,j)
        c%sd_const = sd_const(i,j)

        ! carbon pools — trees
        c%b1t = b1t(i,j) ;  c%b2t = b2t(i,j)
        c%b3t = b3t(i,j) ;  c%b4t = b4t(i,j)
        ! carbon pools — grasses
        c%b1g = b1g(i,j) ;  c%b2g = b2g(i,j)
        c%b3g = b3g(i,j) ;  c%b4g = b4g(i,j)
        ! cell totals
        c%b1 = b1(i,j) ;  c%b2 = b2(i,j)
        c%b3 = b3(i,j) ;  c%b4 = b4(i,j)
        c%stock = stock(i,j)
        c%anup  = anup(i,j)
        c%pnpp  = pnpp(i,j)
        ! LAI
        c%blai(1) = blai(i,j,1) ;  c%blai(2) = blai(i,j,2)

        ! FROG_EXP fields
        c%Fv    = Fv(i,j)
        c%Fv_t  = Fv_t(i,j)
        c%Fv_g  = Fv_g(i,j)
        c%r_leaf = r_leaf(i,j)

#if ( CYCC == 2 )
        c%b1t13 = b1t13(i,j) ;  c%b2t13 = b2t13(i,j)
        c%b3t13 = b3t13(i,j) ;  c%b4t13 = b4t13(i,j)
        c%b1t14 = b1t14(i,j) ;  c%b2t14 = b2t14(i,j)
        c%b3t14 = b3t14(i,j) ;  c%b4t14 = b4t14(i,j)
        c%b1g13 = b1g13(i,j) ;  c%b2g13 = b2g13(i,j)
        c%b3g13 = b3g13(i,j) ;  c%b4g13 = b4g13(i,j)
        c%b1g14 = b1g14(i,j) ;  c%b2g14 = b2g14(i,j)
        c%b3g14 = b3g14(i,j) ;  c%b4g14 = b4g14(i,j)
        c%bc13  = bc13(i,j)
        c%bc14  = bc14(i,j)
        c%sg4   = sg4(i,j)
        c%anup13 = anup13(i,j)
        c%tatmsmin = 0.0_dblp   ! [BUG] not yet connected to atmospheric model
#endif

      end subroutine cell_from_grid

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  INTERNAL: cell_to_grid
!>  Write updated per-cell scalar state back to the global 2-D arrays.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine cell_to_grid(c, i, j)

        type(veget_cell_state_t), intent(in) :: c
        integer(ip),              intent(in) :: i, j

        ! cover fractions
        st(i,j)   = c%st  ;  sg(i,j)   = c%sg  ;  sd(i,j)   = c%sd
        snlt(i,j) = c%snlt
        stR(i,j)  = c%stR ;  sgR(i,j)  = c%sgR ;  sdR(i,j)  = c%sdR
        snltR(i,j) = c%snltR
        st_const(i,j) = c%st_const
        sd_const(i,j) = c%sd_const

        ! carbon pools
        b1t(i,j) = c%b1t ;  b2t(i,j) = c%b2t
        b3t(i,j) = c%b3t ;  b4t(i,j) = c%b4t
        b1g(i,j) = c%b1g ;  b2g(i,j) = c%b2g
        b3g(i,j) = c%b3g ;  b4g(i,j) = c%b4g
        b1(i,j)  = c%b1  ;  b2(i,j)  = c%b2
        b3(i,j)  = c%b3  ;  b4(i,j)  = c%b4
        stock(i,j) = c%stock
        anup(i,j)  = c%anup
        pnpp(i,j)  = c%pnpp
        blai(i,j,1) = c%blai(1)
        blai(i,j,2) = c%blai(2)

        ! FROG_EXP fields
        Fv(i,j)    = c%Fv
        Fv_t(i,j)  = c%Fv_t
        Fv_g(i,j)  = c%Fv_g
        r_leaf(i,j) = c%r_leaf

#if ( CYCC == 2 )
        b1t13(i,j) = c%b1t13 ;  b2t13(i,j) = c%b2t13
        b3t13(i,j) = c%b3t13 ;  b4t13(i,j) = c%b4t13
        b1t14(i,j) = c%b1t14 ;  b2t14(i,j) = c%b2t14
        b3t14(i,j) = c%b3t14 ;  b4t14(i,j) = c%b4t14
        b1g13(i,j) = c%b1g13 ;  b2g13(i,j) = c%b2g13
        b3g13(i,j) = c%b3g13 ;  b4g13(i,j) = c%b4g13
        b1g14(i,j) = c%b1g14 ;  b2g14(i,j) = c%b2g14
        b3g14(i,j) = c%b3g14 ;  b4g14(i,j) = c%b4g14
        bc13(i,j)  = c%bc13
        bc14(i,j)  = c%bc14
        sg4(i,j)   = c%sg4
        anup13(i,j) = c%anup13
#endif

      end subroutine cell_to_grid

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  INTERNAL: read_restart
!>  Read the full vegetation restart file (veget.rest) into the global 2-D arrays.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine read_restart(fid, iyr0, bmtd, titv)

        integer,          intent(out) :: fid
        integer(ip),      intent(out) :: iyr0
        real(dblp),       intent(out) :: bmtd
        character(len=15),intent(out) :: titv

        integer :: ios, ios2

        open(newunit=fid, file='startdata/veget.rest', status='old', &
             form='unformatted')
        read(fid) iyr0, bmtd, titv
        read(fid) st   
        read(fid) sg   
        read(fid) sd
        read(fid) snlt 
        read(fid) temveg
        read(fid) b1t  
        read(fid) b1g
        read(fid) b2t
        read(fid) b2g
        read(fid) b3t
        read(fid) b3g
        read(fid) b4t
        read(fid) b4g
        read(fid) b1
        read(fid) b2
        read(fid) b3
        read(fid) b4
        read(fid) anup

#if ( CYCC == 2 )
        if (new_run_c == 0) then
          read(fid) b1t13 ;  read(fid) b1g13
          read(fid) b2t13 ;  read(fid) b2g13
          read(fid) b3t13 ;  read(fid) b3g13
          read(fid) b4t13 ;  read(fid) b4g13
          read(fid) b1t14 ;  read(fid) b1g14
          read(fid) b2t14 ;  read(fid) b2g14
          read(fid) b3t14 ;  read(fid) b3g14
          read(fid) b4t14 ;  read(fid) b4g14
#if ( KC14 == 1 )
          read(fid) cav_la14_b
          read(fid) cav_la_b
#endif
        else
          ! No isotope restart available: derive the isotope pools from the carbon
          ! pools just read from the vegetation restart.
          !
          ! We loop cell by cell because ccstat_isotope needs this cell's residence
          ! times t1t..t4g.  Those are set by ccparam.  We must call ccparam ALONE
          ! here, NOT full ccstat: ccstat would recompute the equilibrium b* pools
          ! and overwrite the values just loaded from the restart, zeroing the carbon
          ! diagnostics downstream (ca_la_ini == 0 -> SIGFPE / NaN in eco2).
          ! (This was the original bug: ccstat here clobbered the restart pools.)
          do j = 1, nlon
            do i = 1, nlat
              if (fracgr(i,j) <= epss) cycle
              call cell_from_grid(cell, i, j, temveg, gd0veg, &
                                  prcveg, freezeI, thawI)
              cell%crop_fraction = get_crop_fraction(            &
                                     farea, i, j,                &
                                     irunlabelf, 0_ip,           &
                                     i0dfor, ndfor,              &
                                     iscendef, ivegstrt)
              call ccparam(cell)          ! sets t1t..t4g for this cell (no pool overwrite)
              call ccstat_isotope(cell)   ! uses cell%b* (from restart) + t1t..t4g
              call cell_to_grid(cell, i, j)
            end do
          end do
        end if
#endif

        ! guard against negative slow-pool values from corrupt restarts
        where (b4t < 0.0_dblp) b4t = 0.0_dblp
        where (b4g < 0.0_dblp) b4g = 0.0_dblp
        where (b4  < 0.0_dblp) b4  = 0.0_dblp

        ! reference / scenario fractions
        if (iscendef == 1 .or. iscendef == -1) then
          read(fid, end=99, iostat=ios) st_const
99        if (ios /= 0 .or. irunlabelf == ivegstrt) st_const(:,:) = st(:,:)
          sd_const(:,:) = sd(:,:)
          ios2 = 0
          read(fid, end=98, iostat=ios2) stR
          read(fid, end=98, iostat=ios2) sgR
          read(fid, end=98, iostat=ios2) sdR
        end if
98      if (ios2 /= 0 .or. iscendef /= 1) then
          sdR(:,:) = sd(:,:) ;  sgR(:,:) = sg(:,:) ;  stR(:,:) = st(:,:)
        end if
        close(fid)

        snltR(:,:) = snlt(:,:)

        ! initialise stock from pools
        stock(:,:) = b1(:,:) + b2(:,:) + b3(:,:) + b4(:,:)

      end subroutine read_restart

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  INTERNAL: read_veget_init
!>  Read climate accumulators from veget.init (kveget < 0 path).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine read_veget_init(fid, iyrloc, bmtd, titv)

        integer,           intent(out) :: fid
        integer(ip),       intent(out) :: iyrloc
        real(dblp),        intent(out) :: bmtd
        character(len=15), intent(out) :: titv

        open(newunit=fid, file='startdata/veget.init', status='old', &
             form='unformatted')
        read(fid) temveg
        read(fid) gd0veg
        read(fid) gd5veg
        read(fid) ((prcveg(i,j,1), i=1,nlat), j=1,nlon)
        read(fid) ((prcveg(i,j,2), i=1,nlat), j=1,nlon)
        read(fid) tpsdry
        read(fid) iyrloc, bmtd, titv
        close(fid)

      end subroutine read_veget_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  INTERNAL: write_veget_init
!>  Write climate accumulators to veget.init (restart).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine write_veget_init(fid, iyr0, bmtd, titv)

        integer,           intent(out) :: fid
        integer(ip),       intent(in)  :: iyr0
        real(dblp),        intent(in)  :: bmtd
        character(len=15), intent(in)  :: titv

        open(newunit=fid, file='veget.init', status='unknown', &
             form='unformatted')
        write(fid) temveg
        write(fid) gd0veg
        write(fid) gd5veg
        write(fid) ((prcveg(i,j,1), i=1,nlat), j=1,nlon)
        write(fid) ((prcveg(i,j,2), i=1,nlat), j=1,nlon)
        write(fid) tpsdry
        write(fid) iyear + iyr0, bmtd, titv
        close(fid)

      end subroutine write_veget_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  INTERNAL: write_restart
!>  Write the full vegetation restart file (veget.rest).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine write_restart(fid, iyr0, bmtd, titv)

        integer,           intent(out) :: fid
        integer(ip),       intent(in)  :: iyr0
        real(dblp),        intent(in)  :: bmtd
        character(len=15), intent(in)  :: titv

        open(newunit=fid, file='veget.rest', status='unknown', &
             form='unformatted')
        write(fid) iyear + iyr0, bmtd, titv
        write(fid) st   ;  write(fid) sg   ;  write(fid) sd
        write(fid) snlt ;  write(fid) temveg
        write(fid) b1t  ;  write(fid) b1g
        write(fid) b2t  ;  write(fid) b2g
        write(fid) b3t  ;  write(fid) b3g
        write(fid) b4t  ;  write(fid) b4g
        write(fid) b1   ;  write(fid) b2
        write(fid) b3   ;  write(fid) b4
        write(fid) anup

#if ( CYCC == 2 )
        write(fid) b1t13 ;  write(fid) b1g13
        write(fid) b2t13 ;  write(fid) b2g13
        write(fid) b3t13 ;  write(fid) b3g13
        write(fid) b4t13 ;  write(fid) b4g13
        write(fid) b1t14 ;  write(fid) b1g14
        write(fid) b2t14 ;  write(fid) b2g14
        write(fid) b3t14 ;  write(fid) b3g14
        write(fid) b4t14 ;  write(fid) b4g14
#if ( KC14 == 1 )
        write(fid) cav_la14_b
        write(fid) cav_la_b
#endif
#endif
        if (iscendef == 1) then
          write(fid) st_const
          write(fid) stR ;  write(fid) sgR ;  write(fid) sdR
        else
          write(fid) st
        end if
        close(fid)

      end subroutine write_restart

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  INTERNAL: read_deforestation_scenario
!>  Read farea from VEGET.dat (binary) or VEGET.nc (NetCDF, WRAP_EVOL==2).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine read_deforestation_scenario()

        use newunit_mod, only: info_id
        integer :: k, status, ireg_loc, ivar_id, itime_id

#if ( WRAP_EVOL == 2 )
        status = nf_open('inputdata/VEGET.nc', nf_nowrite, ireg_loc)
        status = nf_inq_dimid(ireg_loc, 'time', itime_id)
        status = nf_inq_dimlen(ireg_loc, itime_id, ndfor)
        if (ndfor > mdfor) stop 'veget: ndfor > mdfor — please adjust mdfor in veget_mod'
        status = nf_inq_varid(ireg_loc, 'time', ivar_id)
        status = nf_get_vara_real(ireg_loc, ivar_id, (/1/), (/ndfor/), VegetTime)
        ivegstrt = int(VegetTime(1))
        status = nf_inq_varid(ireg_loc, 'vegfrac', ivar_id)
        status = nf_get_vara_real(ireg_loc, ivar_id, (/1,1,1/), &
                                  (/nlon,nlat,ndfor/), farea)
        ivegstrt = max(i0dfor, ivegstrt)
        write(info_id,*) 'veget: scenario start=', ivegstrt, ' AD, nrecords=', ndfor
#else
        open(newunit=vegetdat_id, file='inputdata/VEGET.dat', form='unformatted')
        read(vegetdat_id) ndfor, i0dfor
        if (ndfor > mdfor) stop 'veget: ndfor > mdfor — please adjust mdfor in veget_mod'
        if (ndfor <= 0)    stop 'veget: wrong header in VEGET.dat'
        do k = 1, ndfor
          read(vegetdat_id) fareal
          farea(:,:,k) = fareal(:,:)
        end do
        close(vegetdat_id)
        ivegstrt = max(i0dfor, ivegstrt)
#endif

      end subroutine read_deforestation_scenario

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  INTERNAL: read_veget_output_param
!>  Read outp_veget.param for NetCDF output variable metadata.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine read_veget_output_param()

        integer :: fid, i_var, k_col
        character(len=60) :: part1

        open(newunit=fid, file='outp_veget.param')
        read(fid,'(/,/,/,/,/,/,/,/)')
        read(fid,*) numvegvar, veg_fill_value, veg_missing_value
        read(fid,*)
        do i_var = 1, numvegvar
          read(fid,'(A)') part1
          namevegvar(i_var,1) = trim(part1)
          read(fid,*) (namevegvar(i_var,k_col), k_col=2,5)
          read(fid,*) (newvegvar(i_var,k_col),  k_col=1,7)
          do k_col = 1, 6
            newvegvar(i_var,k_col) = newvegvar(i_var,k_col+1)
          end do
        end do
        close(fid)

      end subroutine read_veget_output_param

      end subroutine veget

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
