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
!      MODULE: [wisolbm_restart_mod]
!
!>     @brief netCDF restart and initial state of the water isotopes of the LBM land stores (D65, D68, D69; 3a:0006-0007).
!
!      DESCRIPTION:
!
!>     The water stays in the binary restart of the land (inland*.dat, unchanged). This file holds, in double precision,
!!     for iwat16 .. iwat2h (dimension species), the amounts of
!!       soil      soil water bmoisg
!!       snow      snow on land dsnow
!!       runofo    runoff to the ocean of the last step (as inland*.dat keeps the water one)
!!     and the water fields they go with (water_soil, water_snow, water_runofo), to check the pairing. The isotopes of
!!     the land evaporation of the last step are in the atmospheric restart (coevap + clevap, wisoatm_restart).
!!
!!     Reading (lbm_restart_read = .true. in wisoatm.nml, file lbm_restart_file; land_iso_mode = 'bucket'), in
!!     ec_initlbm after the water restart:
!!       - water of the file bitwise equal to the water just read: amounts taken as they are, the run continues
!!         bit-identically;
!!       - otherwise (land isotopes from another state): the ratios of the file are applied to the water of the run
!!         where the file has water, the cold-start ratios (soil, snow init of wisoatm.nml) elsewhere, with a message.
!!     The land and the atmosphere restarts are read independently (restart_read, lbm_restart_read).
!!
!!     Contract
!!       - wlr_write(fname): writes the current land state; called by emic_write_state with the other restarts
!!         (restartdata/resYYYYYY/wisolbm_restartYYYYYY.nc), land_iso_mode = 'bucket' only.
!!       - wlr_read(fname, r_soil, r_snow): stops if the file is missing or its grid does not match; r_soil, r_snow:
!!         ratios where the file has no water (unpaired case).
!!       - wlr_read_init(fname): cold start from a map ('file', D69, wiso:3a:0007), without a land restart. netCDF on
!!         the model grid (lat, lon coordinates in degrees, either latitude order; checked against the model): d18o_snow
!!         (required), dexc_snow, o17exc_snow, d18o_soil, dexc_soil, o17exc_soil (permil, per meg; optional). Where a
!!         variable is absent or a cell is missing (NaN, |value| > 1e10), the *_init values of wisoatm.nml. Snow at the
!!         snow ratios, soil water and runoff at the soil ratios. Made offline: tools/src/wisolbm_init.
!
!      REFERENCES: ROADMAP D48 (atmospheric restart), D64, D65, D68; ncio (netCDF).
!
!      AUTHORS: dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module wisolbm_restart_mod

#if ( WISOATM == 1 )

        use global_constants_mod, only: ip, dblp=>dp
        use comatm,               only: nlat, nlon, iwater, iwat16, iwat2h, nwisos

        implicit none

        private

        public :: wlr_write, wlr_read, wlr_read_init

        integer(ip), parameter :: NSPEC = iwat2h - iwat16 + 1      !< species in the file (iwat16 .. iwat2h)

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Writing
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wlr_write(fname)

          use ncio,        only: nc_create, nc_write_dim, nc_write, nc_write_attr
          use comland_mod, only: bmoisg, dsnow, runofo
          use comemic_mod, only: iyear

          character(len=*), intent(in) :: fname

          call nc_create(fname, overwrite = .true., netcdf4 = .true.,                                                  &
                         description = 'iLOVECLIM LBM land water isotopes restart (wiso:3a:0006, D65)')
          call nc_write_attr(fname, 'model_year', iyear)
          call nc_write_attr(fname, 'species', 'iwat16 iwat17 iwat18 iwat2h')
          call nc_write_attr(fname, 'ratios', 'isotopologue mass per water mass (rsmow = 1)')
          call nc_write_dim(fname, 'lat', x = 1, dx = 1, nx = int(nlat), long_name = 'ECBilt latitude index')
          call nc_write_dim(fname, 'lon', x = 1, dx = 1, nx = int(nlon), long_name = 'ECBilt longitude index')
          call nc_write_dim(fname, 'species', x = int(iwat16), dx = 1, nx = int(NSPEC),                                &
                            long_name = 'species index: iwat16, iwat17, iwat18, iwat2h')

          call nc_write(fname, 'water_soil',   bmoisg(:,:,iwater), dim1 = 'lat', dim2 = 'lon',                        &
                        long_name = 'water of the soil (bmoisg, m, pairing check)')
          call nc_write(fname, 'water_snow',   dsnow(:,:,iwater), dim1 = 'lat', dim2 = 'lon',                         &
                        long_name = 'water of the snow on land (dsnow, m, pairing check)')
          call nc_write(fname, 'water_runofo', runofo(:,:,iwater), dim1 = 'lat', dim2 = 'lon',                        &
                        long_name = 'water runoff to the ocean (runofo, m/s, pairing check)')
          call nc_write(fname, 'soil',   bmoisg(:,:,iwat16:iwat2h), dim1 = 'lat', dim2 = 'lon', dim3 = 'species',     &
                        long_name = 'isotopologue amount of the soil water')
          call nc_write(fname, 'snow',   dsnow(:,:,iwat16:iwat2h), dim1 = 'lat', dim2 = 'lon', dim3 = 'species',      &
                        long_name = 'isotopologue amount of the snow on land')
          call nc_write(fname, 'runofo', runofo(:,:,iwat16:iwat2h), dim1 = 'lat', dim2 = 'lon', dim3 = 'species',     &
                        long_name = 'isotopologue runoff to the ocean')

        end subroutine wlr_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Reading
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wlr_read(fname, r_soil, r_snow)

          use ncio,        only: nc_read, nc_size
          use comland_mod, only: bmoisg, dsnow, runofo

          character(len=*),                 intent(in) :: fname
          real(dblp), dimension(nwisos),    intent(in) :: r_soil, r_snow    !< ratios where the file has no water

          real(dblp)  :: wb(nlat, nlon), ws(nlat, nlon), wr(nlat, nlon)
          real(dblp)  :: soil(nlat, nlon, NSPEC), snow(nlat, nlon, NSPEC), roff(nlat, nlon, NSPEC)
          logical     :: there, same_b, same_s, same_r
          integer(ip) :: k

          inquire(file = fname, exist = there)
          if (.not. there) then
            write(*,*) "wlr_read: ", trim(fname), " not found (lbm_restart_read = .true. in wisoatm.nml)"
            stop 1
          endif
          if (nc_size(fname, 'lat') /= nlat .or. nc_size(fname, 'lon') /= nlon                                         &
              .or. nc_size(fname, 'species') /= NSPEC) then
            write(*,*) "wlr_read: grid of ", trim(fname), " is not (", nlat, nlon, NSPEC, ")"
            stop 1
          endif

          call nc_read(fname, 'water_soil',   wb)
          call nc_read(fname, 'water_snow',   ws)
          call nc_read(fname, 'water_runofo', wr)
          call nc_read(fname, 'soil',   soil)
          call nc_read(fname, 'snow',   snow)
          call nc_read(fname, 'runofo', roff)

          same_b = all(wb == bmoisg(:,:,iwater))
          same_s = all(ws == dsnow(:,:,iwater))
          same_r = all(wr == runofo(:,:,iwater))

          do k = 1, NSPEC
            if (same_b) then
              bmoisg(:,:,iwat16+k-1) = soil(:,:,k)
            else
              bmoisg(:,:,iwat16+k-1) = ratio_or(soil(:,:,k), wb, r_soil(iwat16+k-1))*bmoisg(:,:,iwater)
            endif
            if (same_s) then
              dsnow(:,:,iwat16+k-1) = snow(:,:,k)
            else
              dsnow(:,:,iwat16+k-1) = ratio_or(snow(:,:,k), ws, r_snow(iwat16+k-1))*dsnow(:,:,iwater)
            endif
            if (same_r) then
              runofo(:,:,iwat16+k-1) = roff(:,:,k)
            else
              runofo(:,:,iwat16+k-1) = ratio_or(roff(:,:,k), wr, r_soil(iwat16+k-1))*runofo(:,:,iwater)
            endif
          enddo

          if (same_b .and. same_s .and. same_r) then
            write(*,*) "wlr_read: ", trim(fname), " read, paired with the land restart (exact restart)"
          else
            write(*,*) "wlr_read: ", trim(fname), " read, WATER DIFFERS from the land restart (soil, snow, runoff",   &
                       " paired: ", same_b, same_s, same_r, "): ratios applied to the water of the run"
          endif

        end subroutine wlr_read

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Cold start from a map (wiso:3a:0007, D69)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wlr_read_init(fname)

          use netcdf
          use comatm,      only: phi
          use comland_mod, only: bmoisg, dsnow, runofo
          use iso_nml_mod, only: iso_ratios, d18o_soil_init, dexc_soil_init, o17exc_soil_init, d18o_snow_init,          &
                                 dexc_snow_init, o17exc_snow_init

          character(len=*), intent(in) :: fname

          real(dblp), parameter :: COORD_TOL = 0.01_dblp     !< degrees
          real(dblp), parameter :: MISS_ABS  = 1.0e10_dblp   !< |value| above: missing
          real(dblp), parameter :: PERMIL = 1.0e-3_dblp, PERMEG = 1.0e-6_dblp

          real(dblp)  :: flat(nlat), flon(nlon), mlat(nlat), mlon(nlon)
          real(dblp)  :: v(nlat, nlon, 6), dflt(6)
          real(dblp)  :: rsn(nwisos), rso(nwisos)
          logical     :: there, north_first, found(6)
          integer(ip) :: ncid, st, i, j, n, nmiss(6)
          character(len=11), parameter :: VNAME(6) = ['d18o_snow  ', 'dexc_snow  ', 'o17exc_snow',                      &
                                                      'd18o_soil  ', 'dexc_soil  ', 'o17exc_soil']

          inquire(file = fname, exist = there)
          if (.not. there) then
            write(*,*) "wlr_read_init: ", trim(fname), " not found (land_init_mode = 'file' in wisoatm.nml)"
            stop 1
          endif
          call check(nf90_open(fname, nf90_nowrite, ncid), 'open')

          ! coordinates against the model grid (ECBilt row 1 is the southernmost)
          call get_coord('lat', flat)
          call get_coord('lon', flon)
          mlat = phi(:)*45.0_dblp/atan(1.0_dblp)
          do j = 1, nlon
            mlon(j) = 360.0_dblp*real(j-1, dblp)/real(nlon, dblp)
          enddo
          north_first = flat(1) > flat(nlat)
          if (north_first) flat = flat(nlat:1:-1)
          if (maxval(abs(flat - mlat)) > COORD_TOL) then
            write(*,*) "wlr_read_init: latitudes of ", trim(fname), " differ from the model grid by ",                &
                       maxval(abs(flat - mlat)), " degrees"
            stop 1
          endif
          if (maxval(abs(modulo(flon - mlon + 180.0_dblp, 360.0_dblp) - 180.0_dblp)) > COORD_TOL) then
            write(*,*) "wlr_read_init: longitudes of ", trim(fname), " differ from the model grid (from 0 E, step ",  &
                       "360/nlon)"
            stop 1
          endif

          ! fields, namelist values where absent or missing
          dflt = [d18o_snow_init, dexc_snow_init, o17exc_snow_init, d18o_soil_init, dexc_soil_init, o17exc_soil_init]
          do n = 1, 6
            call get_field(trim(VNAME(n)), v(:,:,n), found(n))
            if (.not. found(n)) v(:,:,n) = dflt(n)
            nmiss(n) = count(v(:,:,n) /= v(:,:,n) .or. abs(v(:,:,n)) > MISS_ABS)
            where (v(:,:,n) /= v(:,:,n) .or. abs(v(:,:,n)) > MISS_ABS) v(:,:,n) = dflt(n)
          enddo
          st = nf90_close(ncid)
          if (.not. found(1)) then
            write(*,*) "wlr_read_init: d18o_snow not found in ", trim(fname)
            stop 1
          endif

          do j = 1, nlon
            do i = 1, nlat
              rsn = iso_ratios(v(i,j,1)*PERMIL, v(i,j,2)*PERMIL, v(i,j,3)*PERMEG)
              rso = iso_ratios(v(i,j,4)*PERMIL, v(i,j,5)*PERMIL, v(i,j,6)*PERMEG)
              dsnow(i,j,iwat16:iwat2h)  = rsn(iwat16:iwat2h)*dsnow(i,j,iwater)
              bmoisg(i,j,iwat16:iwat2h) = rso(iwat16:iwat2h)*bmoisg(i,j,iwater)
              runofo(i,j,iwat16:iwat2h) = rso(iwat16:iwat2h)*runofo(i,j,iwater)
            enddo
          enddo

          write(*,'(a,a,a)') " wlr_read_init: land isotopes from ", trim(fname), " (cold start, land_init_mode = 'file')"
          do n = 1, 6
            if (found(n)) then
              write(*,'(a,a11,a,i0,a)') "   ", VNAME(n), ": read, ", nmiss(n), " missing cells (namelist value there)"
            else
              write(*,'(a,a11,a)') "   ", VNAME(n), ": absent, namelist value everywhere"
            endif
          enddo

        contains

          subroutine check(stat, what)
            integer, intent(in) :: stat
            character(len=*), intent(in) :: what
            if (stat /= nf90_noerr) then
              write(*,*) "wlr_read_init: ", what, " of ", trim(fname), ": ", trim(nf90_strerror(stat))
              stop 1
            endif
          end subroutine check

          subroutine get_coord(name, val)
            character(len=*), intent(in)  :: name
            real(dblp),       intent(out) :: val(:)
            integer :: varid
            call check(nf90_inq_varid(ncid, name, varid), 'coordinate '//name)
            call check(nf90_get_var(ncid, varid, val), 'coordinate '//name)
          end subroutine get_coord

          subroutine get_field(name, val, ok)       ! file (lon, lat) in Fortran order -> model (lat, lon), south first
            character(len=*), intent(in)  :: name
            real(dblp),       intent(out) :: val(nlat, nlon)
            logical,          intent(out) :: ok
            real(dblp) :: tmp(nlon, nlat)
            integer    :: varid
            ok = nf90_inq_varid(ncid, name, varid) == nf90_noerr
            if (.not. ok) return
            call check(nf90_get_var(ncid, varid, tmp), 'variable '//name)
            val = transpose(tmp)
            if (north_first) val = val(nlat:1:-1, :)
          end subroutine get_field

        end subroutine wlr_read_init

        ! ratio qi/q where q > 0, fallback elsewhere
        pure function ratio_or(qi, q, fallback) result(r)

          real(dblp), intent(in) :: qi(:,:), q(:,:), fallback
          real(dblp)             :: r(size(q, 1), size(q, 2))

          where (q > 0.0_dblp)
            r = qi/q
          elsewhere
            r = fallback
          end where

        end function ratio_or

#endif

      end module wisolbm_restart_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
