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
!      MODULE: [fvt_ncout_mod]
!
!>     @brief Minimal netCDF writer for the FV transport test suite (raw netCDF calls, no io_nc dependency).
!
!      DESCRIPTION:
!>     One file holds several fields on the FV grid at a few times: variables (lon, lat, time), double, with lon/lat
!!     bounds (lat bounds = the FV faces, i.e. from the Gaussian weights), cell area on the unit sphere and minimal
!!     CF-1.8 attributes. Time is in days.
!!
!!     Contract (public entry points):
!!       nco_write_fields(path,g,names,fld,times,title) : fld(nlat,nlon,ntime,nfield), names(nfield).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module fvt_ncout_mod

        use global_constants_mod, only: ip, dblp=>dp, PI=>pi_dp
        use fvt_grid_mod,         only: fvt_grid_t
        use netcdf,               only: nf90_create, nf90_def_dim, nf90_def_var, nf90_put_att, nf90_enddef, nf90_put_var, &
                                        nf90_close, nf90_strerror, NF90_CLOBBER, NF90_NETCDF4, NF90_DOUBLE, NF90_GLOBAL,     &
                                        NF90_NOERR

        implicit none

        private

        public :: nco_write_fields

      contains

        subroutine nco_write_fields(path, g, names, fld, times, title)

          character(len=*), intent(in) :: path
          type(fvt_grid_t), intent(in) :: g
          character(len=*), intent(in) :: names(:)
          real(dblp),       intent(in) :: fld(:,:,:,:)    !< (nlat,nlon,ntime,nfield)
          real(dblp),       intent(in) :: times(:)        !< (ntime) days
          character(len=*), intent(in) :: title

          integer(ip) :: ncid, dlon, dlat, dtim, dbnd, vlon, vlat, vtim, vlonb, vlatb, varea, k, nt, it
          integer(ip) :: vid(size(names))
          real(dblp)  :: rad2deg, bnds(2, max(g%nlat, g%nlon)), buf(g%nlon, g%nlat)

          rad2deg = 180.0_dblp/PI
          nt = size(times)
          call chk(nf90_create(path, ior(NF90_CLOBBER, NF90_NETCDF4), ncid))
          call chk(nf90_def_dim(ncid, 'lon', g%nlon, dlon))
          call chk(nf90_def_dim(ncid, 'lat', g%nlat, dlat))
          call chk(nf90_def_dim(ncid, 'time', nt, dtim))
          call chk(nf90_def_dim(ncid, 'bnds', 2, dbnd))
          call chk(nf90_def_var(ncid, 'lon', NF90_DOUBLE, [dlon], vlon))
          call chk(nf90_put_att(ncid, vlon, 'units', 'degrees_east'))
          call chk(nf90_put_att(ncid, vlon, 'standard_name', 'longitude'))
          call chk(nf90_put_att(ncid, vlon, 'bounds', 'lon_bnds'))
          call chk(nf90_def_var(ncid, 'lat', NF90_DOUBLE, [dlat], vlat))
          call chk(nf90_put_att(ncid, vlat, 'units', 'degrees_north'))
          call chk(nf90_put_att(ncid, vlat, 'standard_name', 'latitude'))
          call chk(nf90_put_att(ncid, vlat, 'bounds', 'lat_bnds'))
          call chk(nf90_def_var(ncid, 'time', NF90_DOUBLE, [dtim], vtim))
          call chk(nf90_put_att(ncid, vtim, 'units', 'days since 2000-01-01 00:00:00'))
          call chk(nf90_put_att(ncid, vtim, 'calendar', '360_day'))
          call chk(nf90_def_var(ncid, 'lon_bnds', NF90_DOUBLE, [dbnd, dlon], vlonb))
          call chk(nf90_def_var(ncid, 'lat_bnds', NF90_DOUBLE, [dbnd, dlat], vlatb))
          call chk(nf90_def_var(ncid, 'cell_area', NF90_DOUBLE, [dlon, dlat], varea))
          call chk(nf90_put_att(ncid, varea, 'long_name', 'cell area on the unit sphere'))
          call chk(nf90_put_att(ncid, varea, 'units', '1'))
          do k = 1, size(names)
            call chk(nf90_def_var(ncid, trim(names(k)), NF90_DOUBLE, [dlon, dlat, dtim], vid(k)))
            call chk(nf90_put_att(ncid, vid(k), 'cell_measures', 'area: cell_area'))
          enddo
          call chk(nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'CF-1.8'))
          call chk(nf90_put_att(ncid, NF90_GLOBAL, 'title', title))
          call chk(nf90_put_att(ncid, NF90_GLOBAL, 'source', 'tools/src/fvtrans_tests/fvt_suite (iLOVECLIM)'))
          call chk(nf90_enddef(ncid))

          call chk(nf90_put_var(ncid, vlon, g%lon*rad2deg))
          call chk(nf90_put_var(ncid, vlat, g%lat*rad2deg))
          call chk(nf90_put_var(ncid, vtim, times))
          bnds(1, 1:g%nlon) = g%lon_face(0:g%nlon-1)*rad2deg
          bnds(2, 1:g%nlon) = g%lon_face(1:g%nlon)*rad2deg
          call chk(nf90_put_var(ncid, vlonb, bnds(:, 1:g%nlon)))
          bnds(1, 1:g%nlat) = g%lat_face(0:g%nlat-1)*rad2deg
          bnds(2, 1:g%nlat) = g%lat_face(1:g%nlat)*rad2deg
          call chk(nf90_put_var(ncid, vlatb, bnds(:, 1:g%nlat)))
          buf = transpose(spread(g%area, 2, g%nlon))
          call chk(nf90_put_var(ncid, varea, buf))
          do k = 1, size(names)
            do it = 1, nt
              buf = transpose(fld(:, :, it, k))
              call chk(nf90_put_var(ncid, vid(k), buf, start=[1, 1, it], count=[g%nlon, g%nlat, 1]))
            enddo
          enddo
          call chk(nf90_close(ncid))

        end subroutine nco_write_fields

        subroutine chk(status)

          integer(ip), intent(in) :: status

          if (status /= NF90_NOERR) then
            write(*, '(a)') 'fvt_ncout_mod: '//trim(nf90_strerror(status))
            error stop 'fvt_ncout_mod: netCDF error'
          endif

        end subroutine chk

      end module fvt_ncout_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
