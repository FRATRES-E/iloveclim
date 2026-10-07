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
!      MODULE: [fvt_grid_mod]
!
!>     @brief Finite-volume cells built on a Gaussian grid: faces from the cumulative Gaussian weights.
!
!      DESCRIPTION:
!>     Cell (i,j) is centred on Gaussian latitude i (south to north) and longitude j, as on the ECBilt grid. Its
!!     southern and northern faces are at mu_face(i-1) and mu_face(i), with
!!         mu_face(0) = -1,   mu_face(i) = mu_face(i-1) + w(i),   mu_face(nlat) = +1,
!!     so that the area of every cell on the unit sphere is w(i)*dlon, i.e. exactly the area ECBilt attributes to the
!!     grid point (darea up to the factor radius**2), and a field integrated with the Gaussian weights has the same total
!!     as the FV mass. The poles are faces (zero meridional flux); row 1 and row nlat are rings of nlon cells.
!!     Each Gaussian node lies inside its cell (separation theorem for Gaussian quadrature; checked in the tests).
!!     Longitude j is centred at (j-1)*dlon (j = 1 at 0 degrees, as in the ECBilt FFT), its faces are at
!!     lon_face(j-1) and lon_face(j) = (j-0.5)*dlon.
!!
!!     The construction is symmetric by design: faces are accumulated from the south pole to the equator, the equator
!!     is set to mu = 0 exactly and the northern faces are mirrored, so that mu_face(nlat-i) = -mu_face(i) bitwise.
!!     The FV areas are then w(i)*dlon to round-off (about 1e-16), not bitwise.
!!
!!     Contract (public entry points):
!!       fvt_grid_init(grid,nlat,nlon) : build the grid (nlat even, nlon >= 1).
!!       fvt_grid_free(grid)           : release its arrays.
!!     Everything is on the unit sphere; the caller multiplies areas by radius**2.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module fvt_grid_mod

        use global_constants_mod, only: ip, dblp=>dp, GRID_PI=>pi_dp
        use fvt_legendre_mod,     only: fvt_gauss_nodes

        implicit none

        private

        public :: fvt_grid_t, fvt_grid_init, fvt_grid_free


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Types
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        type :: fvt_grid_t
          integer(ip)             :: nlat = 0          !< number of latitude rows (even)
          integer(ip)             :: nlon = 0          !< number of longitudes
          real(dblp)              :: dlon = 0.0_dblp   !< longitude spacing (radians)
          real(dblp), allocatable :: mu(:)             !< (nlat) sin(latitude) of the cell centres (Gaussian nodes)
          real(dblp), allocatable :: wgt(:)            !< (nlat) Gaussian weights, sum = 2
          real(dblp), allocatable :: lat(:)            !< (nlat) latitude of the cell centres (radians)
          real(dblp), allocatable :: mu_face(:)        !< (0:nlat) sin(latitude) of the faces, -1 and +1 at the poles
          real(dblp), allocatable :: lat_face(:)       !< (0:nlat) latitude of the faces (radians)
          real(dblp), allocatable :: cos_face(:)       !< (0:nlat) cos(latitude) of the faces, 0 at the poles
          real(dblp), allocatable :: dmu(:)            !< (nlat) mu_face(i) - mu_face(i-1), equals wgt to round-off
          real(dblp), allocatable :: area(:)           !< (nlat) cell area on the unit sphere, dmu*dlon
          real(dblp), allocatable :: lon(:)            !< (nlon) longitude of the cell centres (radians)
          real(dblp), allocatable :: lon_face(:)       !< (0:nlon) longitude of the faces (radians)
        end type fvt_grid_t

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Construction and release
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine fvt_grid_init(grid, nlat, nlon)

          type(fvt_grid_t), intent(inout) :: grid
          integer(ip),      intent(in)    :: nlat
          integer(ip),      intent(in)    :: nlon

          integer(ip) :: i, j, nh

          if (nlat < 2 .or. mod(nlat, 2_ip) /= 0) error stop 'fvt_grid_init: nlat must be even and >= 2'
          if (nlon < 1) error stop 'fvt_grid_init: nlon must be >= 1'

          call fvt_grid_free(grid)
          grid%nlat = nlat
          grid%nlon = nlon
          grid%dlon = 2.0_dblp*GRID_PI/real(nlon, dblp)

          allocate(grid%mu(nlat), grid%wgt(nlat), grid%lat(nlat), grid%dmu(nlat), grid%area(nlat))
          allocate(grid%mu_face(0:nlat), grid%lat_face(0:nlat), grid%cos_face(0:nlat))
          allocate(grid%lon(nlon), grid%lon_face(0:nlon))

          call fvt_gauss_nodes(nlat, grid%mu, grid%wgt)
          grid%lat(:) = asin(grid%mu(:))

          ! faces: accumulate the (small) polar weights first, equator exact, mirror to the north
          nh = nlat/2
          grid%mu_face(0) = -1.0_dblp
          do i = 1, nh - 1
            grid%mu_face(i) = grid%mu_face(i-1) + grid%wgt(i)
          enddo
          grid%mu_face(nh) = 0.0_dblp
          do i = 0, nh - 1
            grid%mu_face(nlat-i) = -grid%mu_face(i)
          enddo

          do i = 0, nlat
            grid%lat_face(i) = asin(grid%mu_face(i))
            grid%cos_face(i) = sqrt((1.0_dblp - grid%mu_face(i))*(1.0_dblp + grid%mu_face(i)))
          enddo

          do i = 1, nlat
            grid%dmu(i)  = grid%mu_face(i) - grid%mu_face(i-1)
            grid%area(i) = grid%dmu(i)*grid%dlon
          enddo

          do j = 1, nlon
            grid%lon(j) = real(j - 1, dblp)*grid%dlon
          enddo
          do j = 0, nlon
            grid%lon_face(j) = (real(j, dblp) - 0.5_dblp)*grid%dlon
          enddo

        end subroutine fvt_grid_init

        subroutine fvt_grid_free(grid)

          type(fvt_grid_t), intent(inout) :: grid

          if (allocated(grid%mu))       deallocate(grid%mu)
          if (allocated(grid%wgt))      deallocate(grid%wgt)
          if (allocated(grid%lat))      deallocate(grid%lat)
          if (allocated(grid%mu_face))  deallocate(grid%mu_face)
          if (allocated(grid%lat_face)) deallocate(grid%lat_face)
          if (allocated(grid%cos_face)) deallocate(grid%cos_face)
          if (allocated(grid%dmu))      deallocate(grid%dmu)
          if (allocated(grid%area))     deallocate(grid%area)
          if (allocated(grid%lon))      deallocate(grid%lon)
          if (allocated(grid%lon_face)) deallocate(grid%lon_face)
          grid%nlat = 0
          grid%nlon = 0
          grid%dlon = 0.0_dblp

        end subroutine fvt_grid_free

      end module fvt_grid_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
