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

!   Style sheet: v1.0.0 (applied by analogy: standalone programs are not yet covered)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      PROGRAM: [test_fvt_ffsl]
!
!>     @brief Checks of fvt_ffsl_mod (ROADMAP step 1.1, patch wiso:1a:0002). The full test suite is wiso:1a:0003.
!
!      DESCRIPTION:
!>     Unit sphere (radius 1), time in days. Fluxes of the rotational flows come from the streamfunction at the cell
!!     corners (discretely non-divergent); the divergent flow has analytic face integrals.
!!     Pass/fail checks (invariants):
!!       1. PPM edges exact for cubic profiles on the non-uniform mu grid and on a uniform grid;
!!       2. a uniform density stays uniform under non-divergent flow (solid-body rotation, alpha = 0, pi/4, pi/2);
!!       3. mass conservation to round-off, all schemes;
!!       4. no negative values with the monotone, positive-definite and van Leer options;
!!       5. ratio mode under divergent flow: a uniform ratio stays uniform, carrier and tracer masses conserved.
!!     Reported only (accuracy, for review): Williamson et al. (1992) cosine bell after one revolution (12 days),
!!     normalised l1, l2, linf errors, extrema, substeps used; and the ratio bounds of options C and C' under the
!!     divergent flow.
!!     Stops with exit code 1 if a pass/fail check fails.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      program test_fvt_ffsl

        use global_constants_mod, only: ip, dblp=>dp, PI=>pi_dp
        use fvt_grid_mod,         only: fvt_grid_t, fvt_grid_init
        use fvt_ffsl_mod,         only: fvt_ffsl_density, fvt_ffsl_ratio, fvt_ffsl_nsub, fvt_recon_1d,                &
                                        FVT_RECON_PPM, FVT_RECON_VL, FVT_LIM_NONE, FVT_LIM_MONO, FVT_LIM_POSDEF

        implicit none

        integer(ip), parameter :: NLAT = 32, NLON = 64           !< T21 grid
        real(dblp),  parameter :: U0   = 2.0_dblp*PI/12.0_dblp   !< one revolution in 12 days (radius 1)
        real(dblp),  parameter :: DT   = 1.0_dblp/6.0_dblp       !< 4 h, the ECBilt time step
        integer(ip), parameter :: NSTEP = 72                     !< 12 days
        integer(ip), parameter :: NSCHEME = 4
        integer(ip), parameter :: RECONS(NSCHEME) = [FVT_RECON_PPM, FVT_RECON_PPM, FVT_RECON_PPM, FVT_RECON_VL]
        integer(ip), parameter :: LIMS(NSCHEME)   = [FVT_LIM_NONE, FVT_LIM_MONO, FVT_LIM_POSDEF, FVT_LIM_MONO]
        character(len=8), parameter :: NAMES(NSCHEME) = ['PPM-none', 'PPM-mono', 'PPM-pdef', 'vanLeer ']
        real(dblp),  parameter :: TOL_ROUND = 1.0e-12_dblp       !< round-off tolerance (relative)
        real(dblp),  parameter :: TOL_STEP  = 1.0e-14_dblp       !< round-off growth per substep (uniform ratio)
        real(dblp),  parameter :: TOL_NEG   = 1.0e-14_dblp       !< tolerated negative values (relative to the maximum)

        type(fvt_grid_t) :: g
        integer(ip)      :: nfail

        nfail = 0
        call fvt_grid_init(g, NLAT, NLON)

        call check_edges()
        call check_rotation()
        call check_stability()
        call check_divergent()

        write(*, '(a)') repeat('-', 100)
        if (nfail > 0) then
          write(*, '(i0,a)') nfail, ' check(s) FAILED'
          error stop 1
        endif
        write(*, '(a)') 'all pass/fail checks passed'

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Reporting
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine report(label, value, tol)

          character(len=*), intent(in) :: label
          real(dblp),       intent(in) :: value
          real(dblp),       intent(in) :: tol

          character(len=4) :: verdict

          if (value <= tol) then
            verdict = 'ok'
          else
            verdict = 'FAIL'
            nfail = nfail + 1
          endif
          write(*, '(a4,1x,a,t72,es10.3,a,es8.1)') verdict, label, value, '  tol ', tol

        end subroutine report

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Flows and initial fields (unit sphere)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! Solid-body rotation, Williamson et al. (1992) case 1: psi = -u0 (sin(lat) cos(alpha) - cos(lon) cos(lat) sin(alpha))
        pure function psi_sbr(lon, mu, coslat, alpha) result(psi)

          real(dblp), intent(in) :: lon
          real(dblp), intent(in) :: mu
          real(dblp), intent(in) :: coslat
          real(dblp), intent(in) :: alpha
          real(dblp)             :: psi

          psi = -U0*(mu*cos(alpha) - cos(lon)*coslat*sin(alpha))

        end function psi_sbr

        ! Face volume fluxes over dt from the streamfunction at the corners (exactly non-divergent per cell)
        subroutine fluxes_sbr(alpha, fx, fy)

          real(dblp), intent(in)  :: alpha
          real(dblp), intent(out) :: fx(NLAT, NLON)
          real(dblp), intent(out) :: fy(0:NLAT, NLON)

          real(dblp)  :: psic(0:NLAT, 0:NLON)
          integer(ip) :: i, j

          do j = 0, NLON
            do i = 0, NLAT
              psic(i, j) = psi_sbr(g%lon_face(j), g%mu_face(i), g%cos_face(i), alpha)
            enddo
          enddo
          do j = 1, NLON
            do i = 1, NLAT
              fx(i, j) = -DT*(psic(i, j) - psic(i-1, j))
            enddo
            do i = 0, NLAT
              fy(i, j) = DT*(psic(i, j) - psic(i, j-1))
            enddo
          enddo
          fy(0, :)    = 0.0_dblp
          fy(NLAT, :) = 0.0_dblp

        end subroutine fluxes_sbr

        ! Divergent flow: solid-body rotation (alpha) plus velocity potential chi = c1 sin(lat) + c2 cos(lat) cos(lon)
        subroutine fluxes_div(alpha, c1, c2, fx, fy)

          real(dblp), intent(in)  :: alpha
          real(dblp), intent(in)  :: c1
          real(dblp), intent(in)  :: c2
          real(dblp), intent(out) :: fx(NLAT, NLON)
          real(dblp), intent(out) :: fy(0:NLAT, NLON)

          integer(ip) :: i, j
          real(dblp)  :: lw, le

          call fluxes_sbr(alpha, fx, fy)
          do j = 1, NLON
            lw = g%lon_face(j-1)
            le = g%lon_face(j)
            do i = 1, NLAT
              ! u = -c2 sin(lon): integral of u dlat along the east face
              fx(i, j) = fx(i, j) - DT*c2*sin(le)*(g%lat_face(i) - g%lat_face(i-1))
            enddo
            do i = 1, NLAT - 1
              ! v = c1 cos(lat) - c2 sin(lat) cos(lon): integral of v cos(lat) dlon along the north face
              fy(i, j) = fy(i, j) + DT*g%cos_face(i)*(c1*g%cos_face(i)*(le - lw) - c2*g%mu_face(i)*(sin(le) - sin(lw)))
            enddo
          enddo

        end subroutine fluxes_div

        ! Cosine bell (Williamson et al. 1992), radius 1/3, centre (3 pi/2, 0), height 1; with step = .true. the
        ! indicator of the same disc (sharp field). Cell means by 4 x 4 Gauss-Legendre quadrature in (mu, lon).
        subroutine cell_means_bell(step, h)

          logical,    intent(in)  :: step
          real(dblp), intent(out) :: h(NLAT, NLON)

          real(dblp), parameter :: XG(4) = [-0.861136311594052575_dblp, -0.339981043584856265_dblp, &
                                             0.339981043584856265_dblp,  0.861136311594052575_dblp]
          real(dblp), parameter :: WG(4) = [ 0.347854845137453857_dblp,  0.652145154862546143_dblp, &
                                             0.652145154862546143_dblp,  0.347854845137453857_dblp]
          real(dblp)  :: mu, lon, rr, s, lonc, latc, rb
          integer(ip) :: i, j, a, b

          lonc = 1.5_dblp*PI
          latc = 0.0_dblp
          rb   = 1.0_dblp/3.0_dblp
          do j = 1, NLON
            do i = 1, NLAT
              s = 0.0_dblp
              do a = 1, 4
                mu = 0.5_dblp*(g%mu_face(i-1) + g%mu_face(i)) + 0.5_dblp*g%dmu(i)*XG(a)
                do b = 1, 4
                  lon = g%lon(j) + 0.5_dblp*g%dlon*XG(b)
                  rr  = acos(max(-1.0_dblp, min(1.0_dblp, sin(latc)*mu + cos(latc)*sqrt(1.0_dblp - mu*mu)*cos(lon - lonc))))
                  if (rr < rb) then
                    if (step) then
                      s = s + 0.25_dblp*WG(a)*WG(b)
                    else
                      s = s + 0.25_dblp*WG(a)*WG(b)*0.5_dblp*(1.0_dblp + cos(PI*rr/rb))
                    endif
                  endif
                enddo
              enddo
              h(i, j) = s
            enddo
          enddo

        end subroutine cell_means_bell

        pure function total(h) result(t)

          real(dblp), intent(in) :: h(NLAT, NLON)
          real(dblp)             :: t

          integer(ip) :: i

          t = 0.0_dblp
          do i = 1, NLAT
            t = t + g%area(i)*sum(h(i, :))
          enddo

        end function total

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   1. Reconstruction
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine check_edges()

          integer(ip), parameter :: NG = 2
          real(dblp)  :: xf(-NG:NLAT+NG), h(1-NG:NLAT+NG), q(1-NG:NLAT+NG)
          real(dblp)  :: al(NLAT), ar(NLAT), a6(NLAT), err
          integer(ip) :: i

          write(*, '(a)') repeat('-', 100)
          write(*, '(a)') '1. PPM edge interpolation (unlimited)'
          ! non-uniform: the mu faces, extended beyond the poles with mirrored widths
          xf(0:NLAT) = g%mu_face(:)
          do i = 1, NG
            xf(-i)      = xf(1-i) - g%dmu(i)
            xf(NLAT+i)  = xf(NLAT+i-1) + g%dmu(NLAT+1-i)
          enddo
          call cubic_means(xf, h, q)
          call fvt_recon_1d(NLAT, q, h, FVT_RECON_PPM, FVT_LIM_NONE, al, ar, a6)
          err = max(maxval(abs(al - xf(0:NLAT-1)**3)), maxval(abs(ar - xf(1:NLAT)**3)))
          call report('non-uniform grid (mu faces): max |edge - x**3|, cubic profile', err, 1.0e-14_dblp)
          ! uniform
          do i = -NG, NLAT + NG
            xf(i) = -1.0_dblp + 2.0_dblp*real(i, dblp)/real(NLAT, dblp)
          enddo
          call cubic_means(xf, h, q)
          call fvt_recon_1d(NLAT, q, h, FVT_RECON_PPM, FVT_LIM_NONE, al, ar, a6)
          err = max(maxval(abs(al - xf(0:NLAT-1)**3)), maxval(abs(ar - xf(1:NLAT)**3)))
          call report('uniform grid: max |edge - x**3|, cubic profile', err, 1.0e-14_dblp)

        end subroutine check_edges

        subroutine cubic_means(xf, h, q)

          real(dblp), intent(in)  :: xf(-2:NLAT+2)
          real(dblp), intent(out) :: h(-1:NLAT+2)
          real(dblp), intent(out) :: q(-1:NLAT+2)

          integer(ip) :: i

          do i = -1, NLAT + 2
            h(i) = xf(i) - xf(i-1)
            q(i) = (xf(i)**4 - xf(i-1)**4)/(4.0_dblp*h(i))
          enddo

        end subroutine cubic_means

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   2-4. Solid-body rotation
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine check_rotation()

          real(dblp), parameter :: ALPHAS(3) = [0.0_dblp, 0.25_dblp*PI, 0.5_dblp*PI]
          character(len=6), parameter :: ANAMES(3) = ['0     ', 'pi/4  ', 'pi/2  ']
          real(dblp)  :: fx(NLAT, NLON), fy(0:NLAT, NLON)
          real(dblp)  :: q(NLAT, NLON, 2), h0(NLAT, NLON), dq(NLAT, NLON)
          real(dblp)  :: m0, l1, l2, linf, econst, emass, w(NLAT, NLON)
          integer(ip) :: ia, is, n, ns, i
          character(len=100) :: lab

          write(*, '(a)') repeat('-', 100)
          write(*, '(a)') '2-4. Solid-body rotation, T21 grid, dt = 4 h, 12 days (cosine bell, radius 1/3)'
          call cell_means_bell(.false., h0)
          m0 = total(h0)
          do i = 1, NLAT
            w(i, :) = g%area(i)
          enddo

          do ia = 1, 3
            call fluxes_sbr(ALPHAS(ia), fx, fy)
            write(*, '(a,a,a,i0,a,f6.2)') 'alpha = ', trim(ANAMES(ia)), ': substeps ', fvt_ffsl_nsub(g, fx, fy), &
                                          ', max zonal Courant ', maxval(abs(fx)/spread(g%area, 2, NLON))
            do is = 1, NSCHEME
              q(:, :, 1) = 1.0_dblp
              q(:, :, 2) = h0
              do n = 1, NSTEP
                call fvt_ffsl_density(g, fx, fy, q, RECONS(is), LIMS(is), ns)
              enddo
              econst = maxval(abs(q(:, :, 1) - 1.0_dblp))
              emass  = abs(total(q(:, :, 2)) - m0)/m0
              dq     = q(:, :, 2) - h0
              l1     = sum(abs(dq)*w)/sum(abs(h0)*w)
              l2     = sqrt(sum(dq*dq*w)/sum(h0*h0*w))
              linf   = maxval(abs(dq))/maxval(abs(h0))
              write(lab, '(a,a,a)') '  ', NAMES(is), ': uniform field max |q - 1|'
              call report(trim(lab), econst, TOL_ROUND)
              write(lab, '(a,a,a)') '  ', NAMES(is), ': relative mass change'
              call report(trim(lab), emass, TOL_ROUND)
              if (LIMS(is) /= FVT_LIM_NONE .or. RECONS(is) == FVT_RECON_VL) then
                ! round-off negatives (cancellation in the flux sums) are tolerated, relative to the maximum
                write(lab, '(a,a,a)') '  ', NAMES(is), ': -min(q)/max(q0) (no negative values)'
                call report(trim(lab), max(0.0_dblp, -minval(q(:, :, 2)))/maxval(h0), TOL_NEG)
              endif
              write(*, '(6x,a,3(a,f7.4),2(a,es10.2))') NAMES(is), '  l1 ', l1, '  l2 ', l2, '  linf ', linf, &
                                                      '  min ', minval(q(:, :, 2)), '  max-1 ', maxval(q(:, :, 2)) - maxval(h0)
            enddo
          enddo

        end subroutine check_rotation

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   4b. Stability at small Courant numbers (many substeps)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! Regression test: PPM reconstructed on the non-uniform mu coordinate in y was unstable over the pole when the
        ! Courant numbers per substep are small. Forcing about 60 substeps per step, the unlimited and positive-definite
        ! schemes must not grow.
        subroutine check_stability()

          real(dblp), parameter :: OUT_SMALL = 0.05_dblp
          real(dblp)  :: fx(NLAT, NLON), fy(0:NLAT, NLON), q(NLAT, NLON, 1), h0(NLAT, NLON)
          integer(ip) :: is, n, ns
          character(len=100) :: lab

          write(*, '(a)') repeat('-', 100)
          call fluxes_sbr(0.5_dblp*PI, fx, fy)
          write(*, '(a,i0,a)') '4b. Stability: solid-body rotation alpha = pi/2, 12 days, ', &
                               fvt_ffsl_nsub(g, fx, fy, OUT_SMALL), ' substeps per step'
          call cell_means_bell(.false., h0)
          do is = 1, 3, 2
            q(:, :, 1) = h0
            do n = 1, NSTEP
              call fvt_ffsl_density(g, fx, fy, q, RECONS(is), LIMS(is), ns, out_max=OUT_SMALL)
            enddo
            write(lab, '(a,a,a)') '  ', NAMES(is), ': max|q| / max|q0| (no growth)'
            call report(trim(lab), maxval(abs(q(:, :, 1)))/maxval(h0), 1.0_dblp)
          enddo

        end subroutine check_stability

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   5. Divergent flow: ratio mode (C') against independent densities (C)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine check_divergent()

          real(dblp), parameter :: C1 = 0.05_dblp, C2 = 0.03_dblp     !< velocity potential amplitudes (1/day)
          real(dblp), parameter :: R_UNIF = 0.7_dblp
          real(dblp), parameter :: R_LOW = 0.2_dblp, R_HIGH = 1.0_dblp
          real(dblp)  :: fx(NLAT, NLON), fy(0:NLAT, NLON)
          real(dblp)  :: qc(NLAT, NLON), r(NLAT, NLON, 3), d(NLAT, NLON, 4), bell(NLAT, NLON), step(NLAT, NLON)
          real(dblp)  :: mc0, mt0, rr(NLAT, NLON)
          integer(ip) :: is, n, ns, rec, lim, limr
          character(len=8)   :: name
          character(len=100) :: lab

          write(*, '(a)') repeat('-', 100)
          write(*, '(a)') "5. Divergent flow (rotation pi/4 + chi = 0.05 sin(lat) + 0.03 cos(lat) cos(lon)), 12 days"
          call fluxes_div(0.25_dblp*PI, C1, C2, fx, fy)
          call cell_means_bell(.false., bell)
          call cell_means_bell(.true., step)
          write(*, '(a,i0)') 'substeps ', fvt_ffsl_nsub(g, fx, fy)

          do is = 1, NSCHEME + 1
            ! configuration NSCHEME+1: positive-definite carrier with monotone ratios (C' only; C as PPM-pdef)
            if (is <= NSCHEME) then
              rec  = RECONS(min(is, NSCHEME))
              lim  = LIMS(min(is, NSCHEME))
              limr = lim
              name = NAMES(min(is, NSCHEME))
            else
              rec  = FVT_RECON_PPM
              lim  = FVT_LIM_POSDEF
              limr = FVT_LIM_MONO
              name = 'pdef+mon'
            endif
            ! C': carrier 1, ratios uniform / bell / sharp step between R_LOW and R_HIGH
            qc = 1.0_dblp
            r(:, :, 1) = R_UNIF
            r(:, :, 2) = bell
            r(:, :, 3) = R_LOW + (R_HIGH - R_LOW)*step
            mc0 = total(qc)
            mt0 = total(qc*r(:, :, 2))
            ! C: the same tracers as independent densities
            d(:, :, 1) = qc
            d(:, :, 2) = qc*r(:, :, 1)
            d(:, :, 3) = qc*r(:, :, 2)
            d(:, :, 4) = qc*r(:, :, 3)
            do n = 1, NSTEP
              call fvt_ffsl_ratio(g, fx, fy, qc, r, rec, lim, ns, limiter_ratio=limr)
              call fvt_ffsl_density(g, fx, fy, d, rec, lim, ns)
            enddo
            ! the uniform ratio is exact in structure; round-off accumulates with the number of substeps
            write(lab, '(a,a,a)') '  ', name, " C': uniform ratio max |r - 0.7|"
            call report(trim(lab), maxval(abs(r(:, :, 1) - R_UNIF)), max(TOL_ROUND, TOL_STEP*real(NSTEP*ns, dblp)))
            write(lab, '(a,a,a)') '  ', name, " C': relative carrier mass change"
            call report(trim(lab), abs(total(qc) - mc0)/mc0, TOL_ROUND)
            write(lab, '(a,a,a)') '  ', name, " C': relative tracer mass change (qc*bell)"
            call report(trim(lab), abs(total(qc*r(:, :, 2)) - mt0)/mt0, TOL_ROUND)
            write(lab, '(a,a,a)') '  ', name, " C : relative tracer mass change (bell)"
            call report(trim(lab), abs(total(d(:, :, 3)) - mt0)/mt0, TOL_ROUND)
            rr = d(:, :, 2)/d(:, :, 1)
            write(*, '(6x,a,a,2(f8.4,1x),a,es10.2)') name, "  carrier min max ", minval(qc), maxval(qc), &
                                                   "  C uniform ratio max |q/qc - 0.7| ", maxval(abs(rr - R_UNIF))
            rr = d(:, :, 4)/d(:, :, 1)
            write(*, '(6x,a,a,2(es10.2,1x),a,2(es10.2,1x))') name, "  step ratio, undershoot / overshoot of [0.2,1]:  C' ", &
                     R_LOW - minval(r(:, :, 3)), maxval(r(:, :, 3)) - R_HIGH, "  C ", R_LOW - minval(rr), maxval(rr) - R_HIGH
          enddo

        end subroutine check_divergent

      end program test_fvt_ffsl

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
