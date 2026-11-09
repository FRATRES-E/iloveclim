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
!      MODULE: [wiso_restart_mod]
!
!>     @brief netCDF restart of the atmospheric water isotopologues (D34, ROADMAP step 2c, wiso:2c:0002).
!
!      DESCRIPTION:
!
!>     The water itself stays in the binary restart of the atmosphere (inatphy*.dat, unchanged). This file holds, in
!!     double precision, for iwat16 .. iwat2h (dimension species):
!!       vap, rain, snow   isotopologue amounts rmoisg, torain, tosnow
!!       rvap              vapour ratio memory of wiso_vap_mod (state where the water is ~0, D43)
!!     and the water fields they go with (water_vap, water_rain, water_snow), to check the pairing on reading.
!!
!!     Reading (restart_read = .true. in wisoatm.nml, file restart_file), after the atmospheric restart:
!!       - water of the file bitwise equal to the water just read: amounts and memory taken as they are, the run
!!         continues bit-identically;
!!       - otherwise (isotope restart from another state, or a fresh start): the ratios are kept and applied to the
!!         water of the run (vapour: rvap; rain and snow: ratio of the file where its water > 0, rvap elsewhere),
!!         with a message.
!!
!!     Contract
!!       - wiso_restart_write(fname): writes the current state; called by emic_write_state with the other restarts
!!         (restartdata/resYYYYYY/wisoatm_restartYYYYYY.nc).
!!       - wiso_restart_read(fname): stops if the file is missing or its grid does not match.
!
!      REFERENCES: ROADMAP D34, D43; ncio (netCDF).
!
!      AUTHORS: dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module wiso_restart_mod

#if ( WISOATM == 1 )

        use global_constants_mod, only: ip, dblp=>dp
        use comatm,               only: nlat, nlon, iwater, iwat16, iwat2h

        implicit none

        private

        public :: wiso_restart_write, wiso_restart_read

        integer(ip), parameter :: NSPEC = iwat2h - iwat16 + 1      !< species in the file (iwat16 .. iwat2h)

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Writing
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wiso_restart_write(fname)

          use ncio,          only: nc_create, nc_write_dim, nc_write, nc_write_attr
          use comphys,       only: rmoisg, torain, tosnow
          use comemic_mod,   only: iyear
          use wiso_vap_mod,  only: wiso_vap_get_memory

          character(len=*), intent(in) :: fname

          real(dblp) :: rvap(nlat, nlon, NSPEC)

          call wiso_vap_get_memory(rmoisg, rvap)

          call nc_create(fname, overwrite = .true., netcdf4 = .true.,                                                  &
                         description = 'iLOVECLIM ECBilt water isotopologues restart (wiso:2c, D34)')
          call nc_write_attr(fname, 'model_year', iyear)
          call nc_write_attr(fname, 'species', 'iwat16 iwat17 iwat18 iwat2h')
          call nc_write_attr(fname, 'ratios', 'isotopologue mass per water mass (rsmow = 1)')
          call nc_write_dim(fname, 'lat', x = 1, dx = 1, nx = int(nlat), long_name = 'ECBilt latitude index')
          call nc_write_dim(fname, 'lon', x = 1, dx = 1, nx = int(nlon), long_name = 'ECBilt longitude index')
          call nc_write_dim(fname, 'species', x = int(iwat16), dx = 1, nx = int(NSPEC),                                &
                            long_name = 'species index: iwat16, iwat17, iwat18, iwat2h')

          call nc_write(fname, 'water_vap',  rmoisg(:,:,iwater), dim1 = 'lat', dim2 = 'lon',                          &
                        long_name = 'water vapour column amount (pairing check)')
          call nc_write(fname, 'water_rain', torain(:,:,iwater), dim1 = 'lat', dim2 = 'lon',                          &
                        long_name = 'water rain (pairing check)')
          call nc_write(fname, 'water_snow', tosnow(:,:,iwater), dim1 = 'lat', dim2 = 'lon',                          &
                        long_name = 'water snow (pairing check)')
          call nc_write(fname, 'vap',  rmoisg(:,:,iwat16:iwat2h), dim1 = 'lat', dim2 = 'lon', dim3 = 'species',       &
                        long_name = 'isotopologue vapour column amount')
          call nc_write(fname, 'rain', torain(:,:,iwat16:iwat2h), dim1 = 'lat', dim2 = 'lon', dim3 = 'species',       &
                        long_name = 'isotopologue rain')
          call nc_write(fname, 'snow', tosnow(:,:,iwat16:iwat2h), dim1 = 'lat', dim2 = 'lon', dim3 = 'species',       &
                        long_name = 'isotopologue snow')
          call nc_write(fname, 'rvap', rvap, dim1 = 'lat', dim2 = 'lon', dim3 = 'species',                             &
                        long_name = 'vapour ratio memory')

        end subroutine wiso_restart_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Reading
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wiso_restart_read(fname)

          use ncio,          only: nc_read, nc_size
          use comphys,       only: rmoisg, torain, tosnow
          use wiso_vap_mod,  only: wiso_vap_set_memory

          character(len=*), intent(in) :: fname

          real(dblp)  :: wv(nlat, nlon), wr(nlat, nlon), ws(nlat, nlon)
          real(dblp)  :: vap(nlat, nlon, NSPEC), rain(nlat, nlon, NSPEC), snow(nlat, nlon, NSPEC)
          real(dblp)  :: rvap(nlat, nlon, NSPEC)
          logical     :: there, same_vap, same_prc
          integer(ip) :: k

          inquire(file = fname, exist = there)
          if (.not. there) then
            write(*,*) "wiso_restart_read: ", trim(fname), " not found (restart_read = .true. in wisoatm.nml)"
            stop 1
          endif
          if (nc_size(fname, 'lat') /= nlat .or. nc_size(fname, 'lon') /= nlon                                         &
              .or. nc_size(fname, 'species') /= NSPEC) then
            write(*,*) "wiso_restart_read: grid of ", trim(fname), " is not (", nlat, nlon, NSPEC, ")"
            stop 1
          endif

          call nc_read(fname, 'water_vap',  wv)
          call nc_read(fname, 'water_rain', wr)
          call nc_read(fname, 'water_snow', ws)
          call nc_read(fname, 'vap',  vap)
          call nc_read(fname, 'rain', rain)
          call nc_read(fname, 'snow', snow)
          call nc_read(fname, 'rvap', rvap)

          same_vap = all(wv == rmoisg(:,:,iwater))
          same_prc = all(wr == torain(:,:,iwater)) .and. all(ws == tosnow(:,:,iwater))

          if (same_vap) then
            rmoisg(:,:,iwat16:iwat2h) = vap
          else
            do k = 1, NSPEC
              rmoisg(:,:,iwat16+k-1) = rvap(:,:,k)*rmoisg(:,:,iwater)
            enddo
          endif
          if (same_prc) then
            torain(:,:,iwat16:iwat2h) = rain
            tosnow(:,:,iwat16:iwat2h) = snow
          else
            do k = 1, NSPEC
              torain(:,:,iwat16+k-1) = ratio_or(rain(:,:,k), wr, rvap(:,:,k))*torain(:,:,iwater)
              tosnow(:,:,iwat16+k-1) = ratio_or(snow(:,:,k), ws, rvap(:,:,k))*tosnow(:,:,iwater)
            enddo
          endif
          call wiso_vap_set_memory(rvap)

          if (same_vap .and. same_prc) then
            write(*,*) "wiso_restart_read: ", trim(fname), " read, paired with the atmospheric restart (exact restart)"
          else
            write(*,*) "wiso_restart_read: ", trim(fname), " read, WATER DIFFERS from the atmospheric restart",       &
                       " (vapour, precipitation paired: ", same_vap, same_prc, "): ratios applied to the water of the run"
          endif

        end subroutine wiso_restart_read

        ! ratio qi/q where q > 0, fallback elsewhere
        pure function ratio_or(qi, q, fallback) result(r)

          real(dblp), intent(in) :: qi(:,:), q(:,:), fallback(:,:)
          real(dblp)             :: r(size(q, 1), size(q, 2))

          where (q > 0.0_dblp)
            r = qi/q
          elsewhere
            r = fallback
          end where

        end function ratio_or

#endif

      end module wiso_restart_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
