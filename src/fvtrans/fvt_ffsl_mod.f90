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
!      MODULE: [fvt_ffsl_mod]
!
!>     @brief Flux-form semi-Lagrangian transport (Lin and Rood 1996) on the FV Gaussian grid of fvt_grid_mod.
!
!      DESCRIPTION:
!>     Tracer-generic transport over one time step, given the face volume fluxes of that step. Two modes:
!!       - density mode: each field is a density (amount per unit area, e.g. the ECBilt moisture q) transported
!!         independently in flux form, dq/dt = -div(q v). Options B and C of the roadmap.
!!       - ratio mode: a carrier density qc plus ratios r_k = q_k/qc. The carrier is transported as above; each
!!         tracer amount qc*r_k is transported with the carrier's own face fluxes times a reconstructed ratio
!!         (Lin and Rood's consistency of tracer and air mass). A uniform ratio stays uniform to round-off whatever
!!         the flow. Option C' (isotopologue ratios to iwat16), and later 14C/12C.
!!
!!     Scheme (one substep): with X and Y the 1-D flux-form operators and x(.), y(.) the 1-D advective updates
!!         q_y = [q A + Y-fluxes(q)] / [A + Y-volume fluxes]          (advective inner update, keeps constants)
!!         q_x = idem in x
!!         q_new = q + X-divergence(q + (q_y - q)/2) + Y-divergence(q + (q_x - q)/2)
!!     i.e. Lin and Rood (1996), with the inner updates in the advective form of the FV3 dynamical core (to verify)
!!     (division by the updated volume, no subtraction of a divergence term). Exactly conservative (flux form).
!!
!!     Reconstruction (per 1-D sweep): PPM (Colella and Woodward 1984) on the non-uniform mu coordinate in y and
!!     uniform in x, or piecewise linear with MC-limited slopes (van Leer). PPM limiters: none (4th-order edges,
!!     unlimited, for convergence tests), monotone (CW84), positive-definite (own variant in the spirit of Lin 2004:
!!     edges clipped at 0, then the parabola is moved so that its extremum sits on an edge when its interior minimum
!!     is negative).
!!
!!     Grid and fluxes: fields are (nlat, nlon), row 1 south, as in ECBilt. Face volume fluxes are integrals of
!!     u.n dl dt / radius**2 over the face and the step, i.e. in the units of grid%area (unit sphere):
!!       fx(i,j)   through the east face of cell (i,j), at lon_face(j), positive eastward (periodic in j);
!!       fy(i,j)   through the north face of cell (i,j), at mu_face(i), positive northward, i = 0..nlat,
!!                 fy(0,:) = fy(nlat,:) = 0 (pole faces; checked).
!!     Zonal Courant numbers may exceed 1 (integer shift plus fractional flux, as in LR96). The step is split into
!!     nsub equal substeps when the meridional Courant number (swept volume over upwind cell volume) exceeds
!!     FFSL_CY_MAX, or when a 1-D sweep would remove more than FFSL_DIV_MAX of a cell's volume. Across the poles, the
!!     y reconstruction uses the cells of the opposite meridian (j + nlon/2) as ghosts, so nlon must be even.
!!
!!     Contract (public entry points):
!!       fvt_ffsl_density(grid,fx,fy,q,recon,limiter[,nsub])      : q(nlat,nlon,ntr) densities, updated in place.
!!       fvt_ffsl_ratio(grid,fx,fy,qc,r,recon,limiter[,nsub][,limiter_ratio])
!!                                                                : carrier qc(nlat,nlon) and ratios r(nlat,nlon,nr);
!!                                                                  optional separate limiter for the ratios.
!!       fvt_ffsl_nsub(grid,fx,fy) result(nsub)                   : number of substeps the fluxes require.
!!       fvt_recon_1d(n,q,h,recon,limiter,al,ar,a6)               : 1-D reconstruction (exposed for the tests).
!!       FVT_RECON_PPM, FVT_RECON_VL, FVT_LIM_NONE, FVT_LIM_MONO, FVT_LIM_POSDEF : option values.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module fvt_ffsl_mod

        use global_constants_mod, only: ip, dblp=>dp
        use fvt_grid_mod,         only: fvt_grid_t

        implicit none

        private

        public :: fvt_ffsl_density, fvt_ffsl_ratio, fvt_ffsl_nsub, fvt_recon_1d
        public :: FVT_RECON_PPM, FVT_RECON_VL, FVT_LIM_NONE, FVT_LIM_MONO, FVT_LIM_POSDEF

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Module constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip), parameter :: FVT_RECON_PPM  = 1                !< piecewise parabolic
        integer(ip), parameter :: FVT_RECON_VL   = 2                !< piecewise linear, MC-limited slopes
        integer(ip), parameter :: FVT_LIM_NONE   = 0                !< PPM unlimited (4th-order edges)
        integer(ip), parameter :: FVT_LIM_MONO   = 1                !< PPM monotone (Colella and Woodward 1984)
        integer(ip), parameter :: FVT_LIM_POSDEF = 2                !< PPM positive-definite

        real(dblp),  parameter :: FFSL_CY_MAX  = 1.0_dblp           !< max meridional Courant number per substep
        real(dblp),  parameter :: FFSL_DIV_MAX = 0.5_dblp           !< max net volume loss of a cell per 1-D sweep
        integer(ip), parameter :: FFSL_NG      = 2                  !< ghost cells on each side of a 1-D sweep

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Public drivers
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine fvt_ffsl_density(grid, fx, fy, q, recon, limiter, nsub)

          type(fvt_grid_t),      intent(in)    :: grid
          real(dblp),            intent(in)    :: fx(:,:)       !< (nlat,nlon)   east-face volume fluxes
          real(dblp),            intent(in)    :: fy(0:,:)      !< (0:nlat,nlon) north-face volume fluxes
          real(dblp),            intent(inout) :: q(:,:,:)      !< (nlat,nlon,ntr) densities
          integer(ip),           intent(in)    :: recon         !< FVT_RECON_*
          integer(ip),           intent(in)    :: limiter       !< FVT_LIM_* (PPM only)
          integer(ip), optional, intent(out)   :: nsub          !< substeps used

          real(dblp)  :: fxs(grid%nlat, grid%nlon), fys(0:grid%nlat, grid%nlon)
          integer(ip) :: ns, is, k

          call ffsl_check(grid, fx, fy, recon, limiter)
          if (size(q, 1) /= grid%nlat .or. size(q, 2) /= grid%nlon) error stop 'fvt_ffsl_density: q has the wrong shape'

          ns  = fvt_ffsl_nsub(grid, fx, fy)
          fxs = fx/real(ns, dblp)
          fys = fy/real(ns, dblp)
          do is = 1, ns
            do k = 1, size(q, 3)
              call ffsl_substep_density(grid, fxs, fys, q(:, :, k), recon, limiter)
            enddo
          enddo
          if (present(nsub)) nsub = ns

        end subroutine fvt_ffsl_density

        subroutine fvt_ffsl_ratio(grid, fx, fy, qc, r, recon, limiter, nsub, limiter_ratio)

          type(fvt_grid_t),      intent(in)    :: grid
          real(dblp),            intent(in)    :: fx(:,:)       !< (nlat,nlon)   east-face volume fluxes
          real(dblp),            intent(in)    :: fy(0:,:)      !< (0:nlat,nlon) north-face volume fluxes
          real(dblp),            intent(inout) :: qc(:,:)       !< (nlat,nlon) carrier density
          real(dblp),            intent(inout) :: r(:,:,:)      !< (nlat,nlon,nr) ratios to the carrier
          integer(ip),           intent(in)    :: recon         !< FVT_RECON_*
          integer(ip),           intent(in)    :: limiter       !< FVT_LIM_* (PPM only), for the carrier
          integer(ip), optional, intent(out)   :: nsub          !< substeps used
          integer(ip), optional, intent(in)    :: limiter_ratio !< FVT_LIM_* for the ratios (default: limiter)

          real(dblp)  :: fxs(grid%nlat, grid%nlon), fys(0:grid%nlat, grid%nlon)
          integer(ip) :: ns, is, limr

          call ffsl_check(grid, fx, fy, recon, limiter)
          limr = limiter
          if (present(limiter_ratio)) then
            call ffsl_check(grid, fx, fy, recon, limiter_ratio)
            limr = limiter_ratio
          endif
          if (size(qc, 1) /= grid%nlat .or. size(qc, 2) /= grid%nlon) error stop 'fvt_ffsl_ratio: qc has the wrong shape'
          if (size(r, 1) /= grid%nlat .or. size(r, 2) /= grid%nlon) error stop 'fvt_ffsl_ratio: r has the wrong shape'

          ns  = fvt_ffsl_nsub(grid, fx, fy)
          fxs = fx/real(ns, dblp)
          fys = fy/real(ns, dblp)
          do is = 1, ns
            call ffsl_substep_ratio(grid, fxs, fys, qc, r, recon, limiter, limr)
          enddo
          if (present(nsub)) nsub = ns

        end subroutine fvt_ffsl_ratio

        function fvt_ffsl_nsub(grid, fx, fy) result(nsub)

          type(fvt_grid_t), intent(in) :: grid
          real(dblp),       intent(in) :: fx(:,:)
          real(dblp),       intent(in) :: fy(0:,:)
          integer(ip)                  :: nsub

          integer(ip) :: i, j, jw
          real(dblp)  :: cymax, divmax, c

          cymax  = 0.0_dblp
          divmax = 0.0_dblp
          do j = 1, grid%nlon
            jw = modulo(j - 2, grid%nlon) + 1
            do i = 1, grid%nlat
              ! meridional Courant number at the north face, relative to the upwind cell
              if (i < grid%nlat) then
                if (fy(i, j) >= 0.0_dblp) then
                  c = fy(i, j)/grid%area(i)
                else
                  c = -fy(i, j)/grid%area(i+1)
                endif
                cymax = max(cymax, c)
              endif
              ! net volume leaving the cell in each 1-D sweep
              divmax = max(divmax, (fx(i, j) - fx(i, jw))/grid%area(i))
              divmax = max(divmax, (fy(i, j) - fy(i-1, j))/grid%area(i))
            enddo
          enddo
          nsub = max(1_ip, ceiling(cymax/FFSL_CY_MAX, ip), ceiling(divmax/FFSL_DIV_MAX, ip))

        end function fvt_ffsl_nsub

        subroutine ffsl_check(grid, fx, fy, recon, limiter)

          type(fvt_grid_t), intent(in) :: grid
          real(dblp),       intent(in) :: fx(:,:)
          real(dblp),       intent(in) :: fy(0:,:)
          integer(ip),      intent(in) :: recon
          integer(ip),      intent(in) :: limiter

          if (grid%nlat < 4) error stop 'fvt_ffsl: nlat must be >= 4'
          if (mod(grid%nlon, 2_ip) /= 0) error stop 'fvt_ffsl: nlon must be even (pole ghosts)'
          if (size(fx, 1) /= grid%nlat .or. size(fx, 2) /= grid%nlon) error stop 'fvt_ffsl: fx has the wrong shape'
          if (ubound(fy, 1) /= grid%nlat .or. size(fy, 2) /= grid%nlon) error stop 'fvt_ffsl: fy has the wrong shape'
          if (maxval(abs(fy(0, :))) > 0.0_dblp .or. maxval(abs(fy(grid%nlat, :))) > 0.0_dblp) &
            error stop 'fvt_ffsl: fluxes through the pole faces must be zero'
          if (recon /= FVT_RECON_PPM .and. recon /= FVT_RECON_VL) error stop 'fvt_ffsl: unknown reconstruction'
          if (limiter /= FVT_LIM_NONE .and. limiter /= FVT_LIM_MONO .and. limiter /= FVT_LIM_POSDEF) &
            error stop 'fvt_ffsl: unknown limiter'

        end subroutine ffsl_check

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   One substep
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine ffsl_substep_density(g, fx, fy, q, recon, lim)

          type(fvt_grid_t), intent(in)    :: g
          real(dblp),       intent(in)    :: fx(:,:)
          real(dblp),       intent(in)    :: fy(0:,:)
          real(dblp),       intent(inout) :: q(:,:)
          integer(ip),      intent(in)    :: recon
          integer(ip),      intent(in)    :: lim

          real(dblp)  :: px(g%nlat, g%nlon), py(0:g%nlat, g%nlon)
          real(dblp)  :: qx(g%nlat, g%nlon), qy(g%nlat, g%nlon)
          integer(ip) :: i, j, jw

          ! inner advective updates
          call ffsl_yflux(g, fy, q, q, .false., recon, lim, lim, py)
          call ffsl_xflux(g, fx, q, q, .false., recon, lim, lim, px)
          do j = 1, g%nlon
            jw = modulo(j - 2, g%nlon) + 1
            do i = 1, g%nlat
              qy(i, j) = (q(i, j)*g%area(i) + py(i-1, j) - py(i, j))/(g%area(i) + fy(i-1, j) - fy(i, j))
              qx(i, j) = (q(i, j)*g%area(i) + px(i, jw) - px(i, j))/(g%area(i) + fx(i, jw) - fx(i, j))
            enddo
          enddo

          ! outer flux-form update
          qy = 0.5_dblp*(q + qy)
          qx = 0.5_dblp*(q + qx)
          call ffsl_xflux(g, fx, qy, qy, .false., recon, lim, lim, px)
          call ffsl_yflux(g, fy, qx, qx, .false., recon, lim, lim, py)
          do j = 1, g%nlon
            jw = modulo(j - 2, g%nlon) + 1
            do i = 1, g%nlat
              q(i, j) = q(i, j) + (px(i, jw) - px(i, j) + py(i-1, j) - py(i, j))/g%area(i)
            enddo
          enddo

        end subroutine ffsl_substep_density

        subroutine ffsl_substep_ratio(g, fx, fy, qc, r, recon, lim, limr)

          type(fvt_grid_t), intent(in)    :: g
          real(dblp),       intent(in)    :: fx(:,:)
          real(dblp),       intent(in)    :: fy(0:,:)
          real(dblp),       intent(inout) :: qc(:,:)
          real(dblp),       intent(inout) :: r(:,:,:)
          integer(ip),      intent(in)    :: recon
          integer(ip),      intent(in)    :: lim       !< limiter of the carrier
          integer(ip),      intent(in)    :: limr      !< limiter of the ratios

          real(dblp)  :: pcx(g%nlat, g%nlon), pcy(0:g%nlat, g%nlon)       ! carrier fluxes
          real(dblp)  :: prx(g%nlat, g%nlon), pry(0:g%nlat, g%nlon)       ! tracer (carrier * ratio) fluxes
          real(dblp)  :: qcx(g%nlat, g%nlon), qcy(g%nlat, g%nlon)         ! carrier after the inner updates
          real(dblp)  :: qca(g%nlat, g%nlon), qcb(g%nlat, g%nlon)         ! carrier fed to the outer x, y sweeps
          real(dblp)  :: rx(g%nlat, g%nlon), ry(g%nlat, g%nlon)
          real(dblp)  :: qcn(g%nlat, g%nlon), tn
          integer(ip) :: i, j, jw, k

          ! carrier: exactly the density-mode substep, keeping the inner states and the outer fluxes
          call ffsl_yflux(g, fy, qc, qc, .false., recon, lim, lim, pcy)
          call ffsl_xflux(g, fx, qc, qc, .false., recon, lim, lim, pcx)
          do j = 1, g%nlon
            jw = modulo(j - 2, g%nlon) + 1
            do i = 1, g%nlat
              qcy(i, j) = (qc(i, j)*g%area(i) + pcy(i-1, j) - pcy(i, j))/(g%area(i) + fy(i-1, j) - fy(i, j))
              qcx(i, j) = (qc(i, j)*g%area(i) + pcx(i, jw) - pcx(i, j))/(g%area(i) + fx(i, jw) - fx(i, j))
            enddo
          enddo

          ! carrier outer fluxes and new carrier, computed once and shared by all ratios
          qca = 0.5_dblp*(qc + qcy)
          qcb = 0.5_dblp*(qc + qcx)
          call ffsl_xflux(g, fx, qca, qca, .false., recon, lim, lim, pcx)
          call ffsl_yflux(g, fy, qcb, qcb, .false., recon, lim, lim, pcy)
          do j = 1, g%nlon
            jw = modulo(j - 2, g%nlon) + 1
            do i = 1, g%nlat
              qcn(i, j) = qc(i, j) + (pcx(i, jw) - pcx(i, j) + pcy(i-1, j) - pcy(i, j))/g%area(i)
            enddo
          enddo

          do k = 1, size(r, 3)
            ! inner updates of the ratio, weighted by the carrier (a uniform ratio is left unchanged)
            call ffsl_yflux(g, fy, qc, r(:, :, k), .true., recon, lim, limr, pry)
            call ffsl_xflux(g, fx, qc, r(:, :, k), .true., recon, lim, limr, prx)
            do j = 1, g%nlon
              jw = modulo(j - 2, g%nlon) + 1
              do i = 1, g%nlat
                ry(i, j) = ratio_or_keep(qc(i, j)*r(i, j, k)*g%area(i) + pry(i-1, j) - pry(i, j),    &
                                         qcy(i, j)*(g%area(i) + fy(i-1, j) - fy(i, j)), r(i, j, k))
                rx(i, j) = ratio_or_keep(qc(i, j)*r(i, j, k)*g%area(i) + prx(i, jw) - prx(i, j),     &
                                         qcx(i, j)*(g%area(i) + fx(i, jw) - fx(i, j)), r(i, j, k))
              enddo
            enddo
            ! ratios fed to the outer sweeps: carrier-weighted means of the old and inner-updated states,
            ! consistent with the carrier states qca, qcb fed to the same sweeps
            do j = 1, g%nlon
              do i = 1, g%nlat
                ry(i, j) = ratio_or_keep(qc(i, j)*r(i, j, k) + qcy(i, j)*ry(i, j), qc(i, j) + qcy(i, j), r(i, j, k))
                rx(i, j) = ratio_or_keep(qc(i, j)*r(i, j, k) + qcx(i, j)*rx(i, j), qc(i, j) + qcx(i, j), r(i, j, k))
              enddo
            enddo
            call ffsl_xflux(g, fx, qca, ry, .true., recon, lim, limr, prx)
            call ffsl_yflux(g, fy, qcb, rx, .true., recon, lim, limr, pry)
            do j = 1, g%nlon
              jw = modulo(j - 2, g%nlon) + 1
              do i = 1, g%nlat
                tn = qc(i, j)*r(i, j, k) + (prx(i, jw) - prx(i, j) + pry(i-1, j) - pry(i, j))/g%area(i)
                r(i, j, k) = ratio_or_keep(tn, qcn(i, j), r(i, j, k))
              enddo
            enddo
          enddo

          qc = qcn

        end subroutine ffsl_substep_ratio

        ! num/den, or the fallback value where the carrier has vanished
        pure function ratio_or_keep(num, den, fallback) result(rat)

          real(dblp), intent(in) :: num
          real(dblp), intent(in) :: den
          real(dblp), intent(in) :: fallback
          real(dblp)             :: rat

          if (den > 0.0_dblp) then
            rat = num/den
          else
            rat = fallback
          endif

        end function ratio_or_keep

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   1-D sweeps: face fluxes of a density, or of a carrier-weighted ratio
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! Zonal fluxes. If weighted is false, px = flux of the density a (b unused). If weighted is true, px = flux of
        ! a*b where a is the carrier and b the ratio: sum over the swept cells of (carrier content)*(ratio), the
        ! fractional piece using the means of the two reconstructions over that piece. Units of grid%area times a.
        subroutine ffsl_xflux(g, fx, a, b, weighted, recon, lim, limb, px)

          type(fvt_grid_t), intent(in)  :: g
          real(dblp),       intent(in)  :: fx(:,:)
          real(dblp),       intent(in)  :: a(:,:)
          real(dblp),       intent(in)  :: b(:,:)
          logical,          intent(in)  :: weighted
          integer(ip),      intent(in)  :: recon
          integer(ip),      intent(in)  :: lim       !< limiter of a
          integer(ip),      intent(in)  :: limb      !< limiter of b (weighted only)
          real(dblp),       intent(out) :: px(:,:)

          integer(ip) :: n, i, j, k, m, jj, iu
          real(dblp)  :: h(1-FFSL_NG:g%nlon+FFSL_NG), qe(1-FFSL_NG:g%nlon+FFSL_NG)
          real(dblp)  :: al_a(g%nlon), ar_a(g%nlon), a6_a(g%nlon)
          real(dblp)  :: al_b(g%nlon), ar_b(g%nlon), a6_b(g%nlon)
          real(dblp)  :: c, f, s, ma, mb

          n = g%nlon
          h = 1.0_dblp
          do i = 1, g%nlat
            call periodic_extend(n, a, i, qe)
            call fvt_recon_1d(n, qe, h, recon, lim, al_a, ar_a, a6_a)
            if (weighted) then
              call periodic_extend(n, b, i, qe)
              call fvt_recon_1d(n, qe, h, recon, limb, al_b, ar_b, a6_b)
            endif
            do j = 1, n
              c = fx(i, j)/g%area(i)
              if (c >= 0.0_dblp) then
                ! upwind is to the west: whole cells j, j-1, ..., then a fraction of cell j-k (its eastern part)
                k = int(c, ip)
                f = c - real(k, dblp)
                s = 0.0_dblp
                do m = 0, k - 1
                  jj = modulo(j - 1 - m, n) + 1
                  s  = s + cell_content(a(i, jj), b(i, jj), weighted)
                enddo
                iu = modulo(j - 1 - k, n) + 1
                if (f > 0.0_dblp) then
                  ma = mean_right(al_a(iu), ar_a(iu), a6_a(iu), f)
                  mb = 1.0_dblp
                  if (weighted) mb = mean_right(al_b(iu), ar_b(iu), a6_b(iu), f)
                  s = s + f*ma*mb
                endif
                px(i, j) = s*g%area(i)
              else
                ! upwind is to the east: whole cells j+1, ..., j+k, then a fraction of cell j+1+k (its western part)
                k = int(-c, ip)
                f = -c - real(k, dblp)
                s = 0.0_dblp
                do m = 1, k
                  jj = modulo(j - 1 + m, n) + 1
                  s  = s + cell_content(a(i, jj), b(i, jj), weighted)
                enddo
                iu = modulo(j + k, n) + 1
                if (f > 0.0_dblp) then
                  ma = mean_left(al_a(iu), ar_a(iu), a6_a(iu), f)
                  mb = 1.0_dblp
                  if (weighted) mb = mean_left(al_b(iu), ar_b(iu), a6_b(iu), f)
                  s = s + f*ma*mb
                endif
                px(i, j) = -s*g%area(i)
              endif
            enddo
          enddo

        end subroutine ffsl_xflux

        ! Meridional fluxes, same conventions as ffsl_xflux; |Courant| <= 1 is guaranteed by the substepping.
        subroutine ffsl_yflux(g, fy, a, b, weighted, recon, lim, limb, py)

          type(fvt_grid_t), intent(in)  :: g
          real(dblp),       intent(in)  :: fy(0:,:)
          real(dblp),       intent(in)  :: a(:,:)
          real(dblp),       intent(in)  :: b(:,:)
          logical,          intent(in)  :: weighted
          integer(ip),      intent(in)  :: recon
          integer(ip),      intent(in)  :: lim       !< limiter of a
          integer(ip),      intent(in)  :: limb      !< limiter of b (weighted only)
          real(dblp),       intent(out) :: py(0:,:)

          integer(ip) :: n, i, j
          real(dblp)  :: h(1-FFSL_NG:g%nlat+FFSL_NG), qe(1-FFSL_NG:g%nlat+FFSL_NG)
          real(dblp)  :: al_a(g%nlat), ar_a(g%nlat), a6_a(g%nlat)
          real(dblp)  :: al_b(g%nlat), ar_b(g%nlat), a6_b(g%nlat)
          real(dblp)  :: f, ma, mb

          n = g%nlat
          call pole_extend_h(g, h)
          py(0, :) = 0.0_dblp
          py(n, :) = 0.0_dblp
          do j = 1, g%nlon
            call pole_extend(g, a, j, qe)
            call fvt_recon_1d(n, qe, h, recon, lim, al_a, ar_a, a6_a)
            if (weighted) then
              call pole_extend(g, b, j, qe)
              call fvt_recon_1d(n, qe, h, recon, limb, al_b, ar_b, a6_b)
            endif
            do i = 1, n - 1
              if (fy(i, j) >= 0.0_dblp) then
                ! northward: northern part of cell i
                f  = fy(i, j)/g%area(i)
                ma = mean_right(al_a(i), ar_a(i), a6_a(i), f)
                mb = 1.0_dblp
                if (weighted) mb = mean_right(al_b(i), ar_b(i), a6_b(i), f)
              else
                ! southward: southern part of cell i+1
                f  = -fy(i, j)/g%area(i+1)
                ma = mean_left(al_a(i+1), ar_a(i+1), a6_a(i+1), f)
                mb = 1.0_dblp
                if (weighted) mb = mean_left(al_b(i+1), ar_b(i+1), a6_b(i+1), f)
              endif
              py(i, j) = fy(i, j)*ma*mb
            enddo
          enddo

        end subroutine ffsl_yflux

        pure function cell_content(a, b, weighted) result(s)

          real(dblp), intent(in) :: a
          real(dblp), intent(in) :: b
          logical,    intent(in) :: weighted
          real(dblp)             :: s

          if (weighted) then
            s = a*b
          else
            s = a
          endif

        end function cell_content

        ! Row i with FFSL_NG periodic ghost cells on each side
        subroutine periodic_extend(n, q, i, qe)

          integer(ip), intent(in)  :: n
          real(dblp),  intent(in)  :: q(:,:)
          integer(ip), intent(in)  :: i
          real(dblp),  intent(out) :: qe(1-FFSL_NG:n+FFSL_NG)

          integer(ip) :: j

          do j = 1 - FFSL_NG, n + FFSL_NG
            qe(j) = q(i, modulo(j - 1, n) + 1)
          enddo

        end subroutine periodic_extend

        ! Column j with two ghost cells beyond each pole, taken from the opposite meridian (rows 1, 2 and nlat, nlat-1)
        subroutine pole_extend(g, q, j, qe)

          type(fvt_grid_t), intent(in)  :: g
          real(dblp),       intent(in)  :: q(:,:)
          integer(ip),      intent(in)  :: j
          real(dblp),       intent(out) :: qe(1-FFSL_NG:g%nlat+FFSL_NG)

          integer(ip) :: jo, n, m

          n  = g%nlat
          jo = modulo(j - 1 + g%nlon/2, g%nlon) + 1
          qe(1:n) = q(:, j)
          do m = 1, FFSL_NG
            qe(1-m) = q(m, jo)
            qe(n+m) = q(n+1-m, jo)
          enddo

        end subroutine pole_extend

        subroutine pole_extend_h(g, h)

          type(fvt_grid_t), intent(in)  :: g
          real(dblp),       intent(out) :: h(1-FFSL_NG:g%nlat+FFSL_NG)

          integer(ip) :: n, m

          n = g%nlat
          h(1:n) = g%dmu(:)
          do m = 1, FFSL_NG
            h(1-m) = g%dmu(m)
            h(n+m) = g%dmu(n+1-m)
          enddo

        end subroutine pole_extend_h

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Reconstruction
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! Sub-cell profile of cells 1..n from cell means q (with FFSL_NG ghosts each side) and widths h.
        ! Profile over cell j, xi in [0,1] from its left edge: p(xi) = al + xi*(ar - al + a6*(1 - xi)), mean q(j).
        subroutine fvt_recon_1d(n, q, h, recon, lim, al, ar, a6)

          integer(ip), intent(in)  :: n
          real(dblp),  intent(in)  :: q(1-FFSL_NG:n+FFSL_NG)     !< cell means, with ghosts
          real(dblp),  intent(in)  :: h(1-FFSL_NG:n+FFSL_NG)     !< cell widths, with ghosts
          integer(ip), intent(in)  :: recon
          integer(ip), intent(in)  :: lim
          real(dblp),  intent(out) :: al(n)
          real(dblp),  intent(out) :: ar(n)
          real(dblp),  intent(out) :: a6(n)

          real(dblp)  :: d(0:n+1), e(0:n)
          real(dblp)  :: dl, dr, z1, z2, hs
          integer(ip) :: j
          logical     :: limit_slopes

          ! slopes (CW84 eq. 1.7), MC-limited for van Leer and monotone PPM (CW84 eq. 1.8)
          limit_slopes = (recon == FVT_RECON_VL) .or. (lim == FVT_LIM_MONO)
          do j = 0, n + 1
            dl = q(j) - q(j-1)
            dr = q(j+1) - q(j)
            d(j) = h(j)/(h(j-1) + h(j) + h(j+1))*((2.0_dblp*h(j-1) + h(j))/(h(j+1) + h(j))*dr       &
                                                + (h(j) + 2.0_dblp*h(j+1))/(h(j-1) + h(j))*dl)
            if (limit_slopes) then
              if (dl*dr > 0.0_dblp) then
                d(j) = sign(min(abs(d(j)), 2.0_dblp*abs(dl), 2.0_dblp*abs(dr)), d(j))
              else
                d(j) = 0.0_dblp
              endif
            endif
          enddo

          if (recon == FVT_RECON_VL) then
            do j = 1, n
              al(j) = q(j) - 0.5_dblp*d(j)
              ar(j) = q(j) + 0.5_dblp*d(j)
              a6(j) = 0.0_dblp
            enddo
            return
          endif

          ! edge values between cells j and j+1 (CW84 eq. 1.6, exact for cubic profiles)
          do j = 0, n
            hs = h(j-1) + h(j) + h(j+1) + h(j+2)
            z1 = (h(j-1) + h(j))/(2.0_dblp*h(j) + h(j+1))
            z2 = (h(j+2) + h(j+1))/(2.0_dblp*h(j+1) + h(j))
            e(j) = q(j) + h(j)/(h(j) + h(j+1))*(q(j+1) - q(j))                                       &
                 + (2.0_dblp*h(j+1)*h(j)/(h(j) + h(j+1))*(z1 - z2)*(q(j+1) - q(j))                   &
                    - h(j)*z1*d(j+1) + h(j+1)*z2*d(j))/hs
          enddo

          do j = 1, n
            al(j) = e(j-1)
            ar(j) = e(j)
            select case (lim)
              case (FVT_LIM_MONO)
                call limit_mono(q(j), al(j), ar(j))
              case (FVT_LIM_POSDEF)
                call limit_posdef(q(j), al(j), ar(j))
              case default
                continue
            end select
            a6(j) = 6.0_dblp*q(j) - 3.0_dblp*(al(j) + ar(j))
          enddo

        end subroutine fvt_recon_1d

        ! Colella and Woodward (1984) eq. 1.10
        pure subroutine limit_mono(q, al, ar)

          real(dblp), intent(in)    :: q
          real(dblp), intent(inout) :: al
          real(dblp), intent(inout) :: ar

          real(dblp) :: da, a6

          if ((ar - q)*(q - al) <= 0.0_dblp) then
            al = q
            ar = q
          else
            da = ar - al
            a6 = 6.0_dblp*q - 3.0_dblp*(al + ar)
            if (da*a6 > da*da) then
              al = 3.0_dblp*q - 2.0_dblp*ar
            else if (-da*da > da*a6) then
              ar = 3.0_dblp*q - 2.0_dblp*al
            endif
          endif

        end subroutine limit_mono

        ! Positive definite: edges clipped at zero; if the parabola still has a negative interior minimum, its
        ! extremum is moved to the lower edge (zero slope there), or the cell is flattened.
        pure subroutine limit_posdef(q, al, ar)

          real(dblp), intent(in)    :: q
          real(dblp), intent(inout) :: al
          real(dblp), intent(inout) :: ar

          real(dblp) :: da, a6, pmin

          if (q <= 0.0_dblp) then
            al = q
            ar = q
            return
          endif
          al = max(al, 0.0_dblp)
          ar = max(ar, 0.0_dblp)
          da = ar - al
          a6 = 6.0_dblp*q - 3.0_dblp*(al + ar)
          ! interior minimum only if the parabola is convex (a6 < 0) with its vertex inside the cell
          if (a6 < 0.0_dblp .and. abs(da) < -a6) then
            pmin = q + 0.25_dblp*da*da/a6 + a6/12.0_dblp
            if (pmin < 0.0_dblp) then
              if (q < al .and. q < ar) then
                al = q
                ar = q
              else if (ar > al) then
                ar = 3.0_dblp*q - 2.0_dblp*al
              else
                al = 3.0_dblp*q - 2.0_dblp*ar
              endif
            endif
          endif

        end subroutine limit_posdef

        ! Mean of the profile over the rightmost fraction f of the cell
        pure function mean_right(al, ar, a6, f) result(m)

          real(dblp), intent(in) :: al
          real(dblp), intent(in) :: ar
          real(dblp), intent(in) :: a6
          real(dblp), intent(in) :: f
          real(dblp)             :: m

          m = ar - 0.5_dblp*f*((ar - al) - (1.0_dblp - 2.0_dblp*f/3.0_dblp)*a6)

        end function mean_right

        ! Mean of the profile over the leftmost fraction f of the cell
        pure function mean_left(al, ar, a6, f) result(m)

          real(dblp), intent(in) :: al
          real(dblp), intent(in) :: ar
          real(dblp), intent(in) :: a6
          real(dblp), intent(in) :: f
          real(dblp)             :: m

          m = al + 0.5_dblp*f*((ar - al) + (1.0_dblp - 2.0_dblp*f/3.0_dblp)*a6)

        end function mean_left

      end module fvt_ffsl_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
