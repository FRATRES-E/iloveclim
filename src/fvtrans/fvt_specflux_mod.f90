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
!      MODULE: [fvt_specflux_mod]
!
!>     @brief FV face fluxes from a spectral streamfunction and velocity potential (ECBilt convention).
!
!      DESCRIPTION:
!>     Builds the face volume fluxes expected by fvt_ffsl_mod (units of grid%area, i.e. integral of u.n dl dt /
!!     radius**2) from the spectral coefficients of a streamfunction psi and a velocity potential chi, as exact
!!     integrals of the truncated wind over each face (ROADMAP D22):
!!       - rotational part, u = -(1/a) d(psi)/d(lat), v = (1/(a cos(lat))) d(psi)/d(lon): the flux through a face is
!!         the difference of psi at its two ends, so psi is synthesised at the cell corners and the fluxes are
!!         discretely non-divergent (each cell sums to zero up to round-off);
!!       - divergent part, u = (1/(a cos(lat))) d(chi)/d(lon), v = (1/a) d(chi)/d(lat):
!!           latitude faces: cos(lat)**2 * integral over lon of d(chi)/d(mu), analytic in longitude;
!!           longitude faces: integral over lat of d(chi)/d(lon) / cos(lat), with the integrals of P(n,m)/cos(lat)
!!           over each row precomputed by Gauss-Legendre quadrature in latitude (the integrand is smooth, the
!!           quadrature converges to round-off; checked in test_fvt_specflux).
!!     Spectral convention of ECBilt (as ec_sptogg with the NAG FFT): coefficients c(k,1), c(k,2) for index k = (n,m)
!!     ordered as in fvt_legendre_mod, and
!!         field(lon, mu) = sum over k with m = 0 of c(k,1) P(k)(mu)
!!                        + 2 * sum over k with m > 0 of P(k)(mu) (c(k,1) cos(m lon) - c(k,2) sin(m lon)),
!!     with P the functions of fvt_legendre_mod (the file values of coef.dat, before ECBilt's sqrt(nlon) scaling,
!!     which the NAG normalisation cancels). Checked against ec_sptogg/ec_ddl with the model's coef.dat and nag.f.
!!
!!     Scaling: the caller passes the factors that turn the coefficients into unit-sphere quantities times the
!!     step. For ECBilt's psi(:,3) (non-dimensional, psi_dim = radius**2 * om * psi) and chi(:,3) (dimensional,
!!     m**2/s): scale_psi = umoisr * dt * om, scale_chi = umoisr * dt / radius**2.
!!
!!     Contract (public entry points):
!!       fvt_specflux_init(sf,grid,ntrunc[,nq])  : precompute the Legendre and trigonometric tables.
!!       fvt_specflux_fluxes(sf,grid,psi,chi,scale_psi,scale_chi,fx,fy) : face fluxes (fx(nlat,nlon), fy(0:nlat,nlon)).
!!       fvt_specflux_free(sf)                    : release the tables.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module fvt_specflux_mod

        use global_constants_mod, only: ip, dblp=>dp
        use fvt_legendre_mod,     only: fvt_gauss_nodes, fvt_spec_index, fvt_legendre_eval, fvt_nsh
        use fvt_grid_mod,         only: fvt_grid_t

        implicit none

        private

        public :: fvt_specflux_t, fvt_specflux_init, fvt_specflux_fluxes, fvt_specflux_free

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Module constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip), parameter :: SPF_NQ_DEFAULT = 32   !< quadrature points per row for the longitude-face integrals

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Types
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        type :: fvt_specflux_t
          integer(ip)              :: ntrunc = -1     !< triangular truncation
          integer(ip)              :: nsh    = 0      !< number of (n,m) pairs
          integer(ip), allocatable :: mk(:)           !< (nsh) zonal wavenumber of index k
          real(dblp),  allocatable :: pface(:,:)      !< (0:nlat, nsh) P at the latitude faces (corners)
          real(dblp),  allocatable :: dface(:,:)      !< (1:nlat-1, nsh) dP/dmu at the interior latitude faces
          real(dblp),  allocatable :: prow(:,:)       !< (nlat, nsh) integral over row i of P/cos(lat) dlat
          real(dblp),  allocatable :: cosm(:,:)       !< (0:ntrunc, 0:nlon) cos(m lon_face(j))
          real(dblp),  allocatable :: sinm(:,:)       !< (0:ntrunc, 0:nlon) sin(m lon_face(j))
        end type fvt_specflux_t

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Construction and release
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine fvt_specflux_init(sf, grid, ntrunc, nq)

          type(fvt_specflux_t),  intent(inout) :: sf
          type(fvt_grid_t),      intent(in)    :: grid
          integer(ip),           intent(in)    :: ntrunc
          integer(ip), optional, intent(in)    :: nq          !< quadrature points per row (default SPF_NQ_DEFAULT)

          integer(ip)              :: nlat, nlon, nqq, i, j, k, m, n, q
          integer(ip), allocatable :: nshm(:), ll(:)
          real(dblp),  allocatable :: xq(:), wq(:), muq(:), cq(:), pq(:,:), work(:,:)
          real(dblp)               :: d0, d1, d, half, sgn

          nlat = grid%nlat
          nlon = grid%nlon
          if (ntrunc < 0) error stop 'fvt_specflux_init: ntrunc must be >= 0'
          if (2*ntrunc + 1 > nlon) error stop 'fvt_specflux_init: nlon too small for ntrunc'
          nqq = SPF_NQ_DEFAULT
          if (present(nq)) nqq = nq
          if (nqq < 2 .or. mod(nqq, 2_ip) /= 0) error stop 'fvt_specflux_init: nq must be even and >= 2'

          call fvt_specflux_free(sf)
          sf%ntrunc = ntrunc
          sf%nsh    = fvt_nsh(ntrunc)
          allocate(nshm(0:ntrunc), ll(sf%nsh), sf%mk(sf%nsh))
          call fvt_spec_index(ntrunc, nshm, ll)
          k = 0
          do m = 0, ntrunc
            do n = m, ntrunc
              k = k + 1
              sf%mk(k) = m
            enddo
          enddo

          ! P at all latitude faces (poles included), dP/dmu at the interior ones
          allocate(sf%pface(0:nlat, sf%nsh), sf%dface(1:nlat-1, sf%nsh), work(nlat-1, sf%nsh))
          call fvt_legendre_eval(ntrunc, grid%mu_face(0:nlat), sf%pface)
          call fvt_legendre_eval(ntrunc, grid%mu_face(1:nlat-1), work, sf%dface)

          ! integral of P/cos(lat) over each row, Gauss-Legendre in latitude. The nodes are placed by their angular
          ! distance d from the nearer pole, so that cos(lat) = sin(d) and mu = -+cos(d) stay accurate near the poles
          ! (computing cos(lat) from mu, or mu from lat, loses up to half the digits in the polar rows)
          allocate(xq(nqq), wq(nqq), muq(nqq), cq(nqq), pq(nqq, sf%nsh), sf%prow(nlat, sf%nsh))
          call fvt_gauss_nodes(nqq, xq, wq)
          do i = 1, nlat
            if (i <= nlat/2) then
              ! southern hemisphere: distance from the south pole, d = lat + pi/2
              d0 = pole_distance(grid%mu_face(i-1))
              d1 = pole_distance(grid%mu_face(i))
              sgn = -1.0_dblp
            else
              ! northern hemisphere: distance from the north pole, d = pi/2 - lat
              d0 = pole_distance(-grid%mu_face(i))
              d1 = pole_distance(-grid%mu_face(i-1))
              sgn = 1.0_dblp
            endif
            half = 0.5_dblp*(d1 - d0)
            do q = 1, nqq
              d      = 0.5_dblp*(d0 + d1) + half*xq(q)
              cq(q)  = sin(d)
              muq(q) = sgn*cos(d)
            enddo
            call fvt_legendre_eval(ntrunc, muq, pq, coslat=cq)
            do k = 1, sf%nsh
              sf%prow(i, k) = 0.0_dblp
              do q = 1, nqq
                sf%prow(i, k) = sf%prow(i, k) + half*wq(q)*pq(q, k)/cq(q)
              enddo
            enddo
          enddo

          ! trigonometric tables at the longitude faces
          allocate(sf%cosm(0:ntrunc, 0:nlon), sf%sinm(0:ntrunc, 0:nlon))
          do j = 0, nlon
            do m = 0, ntrunc
              sf%cosm(m, j) = cos(real(m, dblp)*grid%lon_face(j))
              sf%sinm(m, j) = sin(real(m, dblp)*grid%lon_face(j))
            enddo
          enddo

        end subroutine fvt_specflux_init

        subroutine fvt_specflux_free(sf)

          type(fvt_specflux_t), intent(inout) :: sf

          if (allocated(sf%mk))    deallocate(sf%mk)
          if (allocated(sf%pface)) deallocate(sf%pface)
          if (allocated(sf%dface)) deallocate(sf%dface)
          if (allocated(sf%prow))  deallocate(sf%prow)
          if (allocated(sf%cosm))  deallocate(sf%cosm)
          if (allocated(sf%sinm))  deallocate(sf%sinm)
          sf%ntrunc = -1
          sf%nsh    = 0

        end subroutine fvt_specflux_free

        ! Angular distance from the south pole of the latitude whose sine is mu, accurate near mu = -1:
        ! 1 + mu = 2 sin(d/2)**2 (mu_face is built from cumulative Gaussian weights, so 1 + mu is accurate there)
        pure function pole_distance(mu) result(d)

          real(dblp), intent(in) :: mu
          real(dblp)             :: d

          d = 2.0_dblp*asin(sqrt(max(0.0_dblp, 0.5_dblp*(1.0_dblp + mu))))

        end function pole_distance

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Face fluxes
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine fvt_specflux_fluxes(sf, grid, psi, chi, scale_psi, scale_chi, fx, fy)

          type(fvt_specflux_t), intent(in)  :: sf
          type(fvt_grid_t),     intent(in)  :: grid
          real(dblp),           intent(in)  :: psi(:,:)     !< (nsh,2) streamfunction coefficients
          real(dblp),           intent(in)  :: chi(:,:)     !< (nsh,2) velocity potential coefficients
          real(dblp),           intent(in)  :: scale_psi    !< coefficients -> unit-sphere streamfunction * step
          real(dblp),           intent(in)  :: scale_chi    !< coefficients -> unit-sphere velocity potential * step
          real(dblp),           intent(out) :: fx(:,:)      !< (nlat,nlon) east-face volume fluxes
          real(dblp),           intent(out) :: fy(0:,:)     !< (0:nlat,nlon) north-face volume fluxes

          integer(ip) :: nlat, nlon, nt, i, j, m
          real(dblp)  :: psic(0:grid%nlat, 0:grid%nlon)
          real(dblp)  :: ac(0:grid%nlat, 0:sf%ntrunc), bc(0:grid%nlat, 0:sf%ntrunc)
          real(dblp)  :: ad(1:grid%nlat-1, 0:sf%ntrunc), bd(1:grid%nlat-1, 0:sf%ntrunc)
          real(dblp)  :: ar(1:grid%nlat, 0:sf%ntrunc), br(1:grid%nlat, 0:sf%ntrunc)
          real(dblp)  :: s, rm

          nlat = grid%nlat
          nlon = grid%nlon
          nt   = sf%ntrunc
          if (size(psi, 1) /= sf%nsh .or. size(chi, 1) /= sf%nsh .or. size(psi, 2) /= 2 .or. size(chi, 2) /= 2) &
            error stop 'fvt_specflux_fluxes: coefficient arrays have the wrong shape'
          if (size(fx, 1) /= nlat .or. size(fx, 2) /= nlon) error stop 'fvt_specflux_fluxes: fx has the wrong shape'
          if (ubound(fy, 1) /= nlat .or. size(fy, 2) /= nlon) error stop 'fvt_specflux_fluxes: fy has the wrong shape'

          ! Fourier coefficients in longitude of psi at the corners, of d(chi)/d(mu) on the latitude faces and of the
          ! row integrals of chi/cos(lat) (the longitude derivative is applied in the synthesis)
          call spf_fourier(sf, sf%pface, psi, ac, bc)
          call spf_fourier(sf, sf%dface, chi, ad, bd)
          call spf_fourier(sf, sf%prow, chi, ar, br)

          ! rotational part: psi at the corners
          do j = 0, nlon
            do i = 0, nlat
              s = ac(i, 0)
              do m = 1, nt
                s = s + 2.0_dblp*(ac(i, m)*sf%cosm(m, j) - bc(i, m)*sf%sinm(m, j))
              enddo
              psic(i, j) = s
            enddo
          enddo
          do j = 1, nlon
            do i = 1, nlat
              fx(i, j) = -scale_psi*(psic(i, j) - psic(i-1, j))
            enddo
            do i = 1, nlat - 1
              fy(i, j) = scale_psi*(psic(i, j) - psic(i, j-1))
            enddo
          enddo

          ! divergent part, longitude faces: integral over the row of d(chi)/d(lon) / cos(lat) at lon_face(j)
          do j = 1, nlon
            do i = 1, nlat
              s = 0.0_dblp
              do m = 1, nt
                rm = real(m, dblp)
                s  = s + 2.0_dblp*rm*(-ar(i, m)*sf%sinm(m, j) - br(i, m)*sf%cosm(m, j))
              enddo
              fx(i, j) = fx(i, j) + scale_chi*s
            enddo
          enddo

          ! divergent part, latitude faces: cos(lat)**2 * integral over the cell's longitudes of d(chi)/d(mu)
          do j = 1, nlon
            do i = 1, nlat - 1
              s = ad(i, 0)*(grid%lon_face(j) - grid%lon_face(j-1))
              do m = 1, nt
                rm = real(m, dblp)
                s  = s + 2.0_dblp*(ad(i, m)*(sf%sinm(m, j) - sf%sinm(m, j-1))                                    &
                                   + bd(i, m)*(sf%cosm(m, j) - sf%cosm(m, j-1)))/rm
              enddo
              fy(i, j) = fy(i, j) + scale_chi*grid%cos_face(i)**2*s
            enddo
          enddo

          ! pole faces carry no flux
          fy(0, :)    = 0.0_dblp
          fy(nlat, :) = 0.0_dblp

        end subroutine fvt_specflux_fluxes

        ! Fourier coefficients a(i,m), b(i,m) = sum over n of c(k,1|2) * tab(i,k), for the rows i of a table
        subroutine spf_fourier(sf, tab, c, a, b)

          type(fvt_specflux_t), intent(in)  :: sf
          real(dblp),           intent(in)  :: tab(:,:)
          real(dblp),           intent(in)  :: c(:,:)
          real(dblp),           intent(out) :: a(:,0:)
          real(dblp),           intent(out) :: b(:,0:)

          integer(ip) :: k, m

          a = 0.0_dblp
          b = 0.0_dblp
          do k = 1, sf%nsh
            m = sf%mk(k)
            a(:, m) = a(:, m) + c(k, 1)*tab(:, k)
            b(:, m) = b(:, m) + c(k, 2)*tab(:, k)
          enddo

        end subroutine spf_fourier

      end module fvt_specflux_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
