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
!      PROGRAM: [test_fvt_specflux]
!
!>     @brief Checks of the face fluxes of fvt_specflux_mod (ROADMAP step 1.1, patches wiso:1d:0001, wiso:1f:0002).
!
!      DESCRIPTION:
!>     Usage: test_fvt_specflux (no argument, no input file)
!!     Random spectral psi and chi with a red spectrum (fixed seed). Pass/fail checks, at T21 on a 32 x 64 grid and
!!     at T42 on a 64 x 128 grid:
!!       1. rotational fluxes are discretely non-divergent (each cell sums to zero up to round-off);
!!       2. divergent fluxes: the outflow of each cell equals the integral of the divergence over the cell
!!          (Gauss theorem; divergence integrated independently by 16 x 16 quadrature);
!!       3. the row integrals have converged in the quadrature (default against 2x the points);
!!       4. pole faces carry no flux;
!!       5. a uniform field transported by fvt_ffsl with the rotational fluxes stays uniform.
!!     The spectral convention itself (against ECBilt's transforms) is checked by test_fvt_specflux_ecbilt.
!!     Stops with exit code 1 if a check fails.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      program test_fvt_specflux

        use fvt_kinds_mod,        only: ip=>fvt_ip, dblp=>fvt_dp
        use fvt_legendre_mod,     only: fvt_gauss_nodes, fvt_legendre_eval, fvt_nsh
        use fvt_grid_mod,         only: fvt_grid_t, fvt_grid_init
        use fvt_specflux_mod,     only: fvt_specflux_t, fvt_specflux_init, fvt_specflux_fluxes
        use fvt_ffsl_mod,         only: fvt_ffsl_density, FVT_RECON_PPM, FVT_LIM_POSDEF

        implicit none

        real(dblp), parameter :: TOL_ROUND = 1.0e-12_dblp

        integer(ip)            :: nfail

        nfail = 0
        call check_fluxes(21_ip, 32_ip, 64_ip)
        call check_fluxes(42_ip, 64_ip, 128_ip)

        write(*, '(a)') repeat('-', 100)
        if (nfail > 0) then
          write(*, '(i0,a)') nfail, ' check(s) FAILED'
          error stop 1
        endif
        write(*, '(a)') 'all checks passed'

      contains

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
          write(*, '(a4,1x,a,t76,es10.3,a,es8.1)') verdict, label, value, '  tol ', tol

        end subroutine report

        ! Random coefficients with amplitude ~ 1/(n(n+1)) (red spectrum), zero imaginary part for m = 0
        subroutine random_coefficients(nt, seed, c)

          integer(ip), intent(in)  :: nt
          integer(ip), intent(in)  :: seed
          real(dblp),  intent(out) :: c(:,:)

          integer(ip) :: k, m, n
          real(dblp)  :: r(2)
          integer     :: nseed
          integer, allocatable :: s(:)

          call random_seed(size=nseed)
          allocate(s(nseed))
          s = int(seed) + 37*[(int(k), k = 1, nseed)]
          call random_seed(put=s)
          k = 0
          do m = 0, nt
            do n = m, nt
              k = k + 1
              call random_number(r)
              c(k, :) = (2.0_dblp*r - 1.0_dblp)/real(max(1_ip, n*(n + 1)), dblp)
              if (m == 0) c(k, 2) = 0.0_dblp
            enddo
          enddo
          c(1, :) = 0.0_dblp

        end subroutine random_coefficients

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   1-5. Face fluxes
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine check_fluxes(nt, nl, nlo)

          integer(ip), intent(in) :: nt
          integer(ip), intent(in) :: nl
          integer(ip), intent(in) :: nlo

          type(fvt_grid_t)     :: g
          type(fvt_specflux_t) :: sf, sf2
          real(dblp), allocatable :: psi(:,:), chi(:,:), zero(:,:), fx(:,:), fy(:,:), fx2(:,:), fy2(:,:), q(:,:,:)
          real(dblp), allocatable :: divint(:,:)
          real(dblp)  :: div, fmax, e, cfl
          integer(ip) :: i, j, jw, n, ns
          character(len=40) :: tag

          write(tag, '(a,i0,a,i0,a,i0)') 'T', nt, ', ', nl, ' x ', nlo
          write(*, '(a)') repeat('-', 100)
          write(*, '(a)') '1-5. Face fluxes, '//trim(tag)
          call fvt_grid_init(g, nl, nlo)
          call fvt_specflux_init(sf, g, nt)
          call fvt_specflux_init(sf2, g, nt, nq=64_ip)
          allocate(psi(fvt_nsh(nt), 2), chi(fvt_nsh(nt), 2), zero(fvt_nsh(nt), 2))
          allocate(fx(nl, nlo), fy(0:nl, nlo), fx2(nl, nlo), fy2(0:nl, nlo), divint(nl, nlo), q(nl, nlo, 1))
          call random_coefficients(nt, 11_ip, psi)
          call random_coefficients(nt, 23_ip, chi)
          zero = 0.0_dblp

          ! 2. rotational part alone
          call fvt_specflux_fluxes(sf, g, psi, zero, 1.0_dblp, 0.0_dblp, fx, fy)
          fmax = max(maxval(abs(fx)), maxval(abs(fy)))
          e = 0.0_dblp
          do j = 1, nlo
            jw = modulo(j - 2, nlo) + 1
            do i = 1, nl
              div = fx(i, j) - fx(i, jw) + fy(i, j) - fy(i-1, j)
              e = max(e, abs(div))
            enddo
          enddo
          call report('rotational: max |cell net outflow| / max |face flux|', e/fmax, 1.0e-14_dblp)

          ! 3. divergent part alone: net outflow = integral of the divergence over the cell
          call fvt_specflux_fluxes(sf, g, zero, chi, 0.0_dblp, 1.0_dblp, fx, fy)
          call cell_divergence(g, nt, chi, divint)
          fmax = max(maxval(abs(fx)), maxval(abs(fy)))
          e = 0.0_dblp
          do j = 1, nlo
            jw = modulo(j - 2, nlo) + 1
            do i = 1, nl
              div = fx(i, j) - fx(i, jw) + fy(i, j) - fy(i-1, j)
              e = max(e, abs(div - divint(i, j)))
            enddo
          enddo
          call report('divergent: max |net outflow - integral of div| / max |face flux|', e/fmax, TOL_ROUND)

          ! 4. quadrature convergence of the row integrals (default 32 against 64 points)
          call fvt_specflux_fluxes(sf2, g, zero, chi, 0.0_dblp, 1.0_dblp, fx2, fy2)
          call report('row integrals: max |fx(nq=32) - fx(nq=64)| / max |fx|', maxval(abs(fx - fx2))/maxval(abs(fx)), &
                      1.0e-14_dblp)

          ! 5. pole faces
          call fvt_specflux_fluxes(sf, g, psi, chi, 1.0_dblp, 1.0_dblp, fx, fy)
          call report('pole faces: max |fy(0,:)|, |fy(nlat,:)|', max(maxval(abs(fy(0, :))), maxval(abs(fy(nl, :)))), &
                      0.0_dblp)

          ! 6. uniform field under the rotational fluxes, scaled to a zonal Courant number of about 1 at the equator
          call fvt_specflux_fluxes(sf, g, psi, zero, 1.0_dblp, 0.0_dblp, fx, fy)
          cfl = maxval(abs(fx(nl/2, :)))/g%area(nl/2)
          call fvt_specflux_fluxes(sf, g, psi, zero, 1.0_dblp/cfl, 0.0_dblp, fx, fy)
          q = 1.0_dblp
          do n = 1, 50
            call fvt_ffsl_density(g, fx, fy, q, FVT_RECON_PPM, FVT_LIM_POSDEF, ns)
          enddo
          call report('uniform field, 50 fvt_ffsl steps, rotational fluxes: max |q - 1|', &
                      maxval(abs(q - 1.0_dblp)), TOL_ROUND)
          write(*, '(6x,a,i0)') 'substeps per step: ', ns

        end subroutine check_fluxes

        ! Integral over each cell of the unit-sphere Laplacian of chi, by 16 x 16 Gauss-Legendre quadrature in (lat, lon)
        subroutine cell_divergence(g, nt, chi, divint)

          type(fvt_grid_t), intent(in)  :: g
          integer(ip),      intent(in)  :: nt
          real(dblp),       intent(in)  :: chi(:,:)
          real(dblp),       intent(out) :: divint(:,:)

          integer(ip), parameter :: NQ = 16
          real(dblp)  :: xq(NQ), wq(NQ), mu(NQ), p(NQ, fvt_nsh(nt)), lat, half, hl, lon, s, f, lam
          integer(ip) :: i, j, a, b, k, m, n

          call fvt_gauss_nodes(NQ, xq, wq)
          do i = 1, g%nlat
            half = 0.5_dblp*(g%lat_face(i) - g%lat_face(i-1))
            do a = 1, NQ
              mu(a) = sin(0.5_dblp*(g%lat_face(i) + g%lat_face(i-1)) + half*xq(a))
            enddo
            call fvt_legendre_eval(nt, mu, p)
            do j = 1, g%nlon
              hl = 0.5_dblp*g%dlon
              s = 0.0_dblp
              do a = 1, NQ
                lat = 0.5_dblp*(g%lat_face(i) + g%lat_face(i-1)) + half*xq(a)
                do b = 1, NQ
                  lon = g%lon(j) + hl*xq(b)
                  f = 0.0_dblp
                  k = 0
                  do m = 0, nt
                    do n = m, nt
                      k = k + 1
                      lam = -real(n*(n + 1), dblp)
                      if (m == 0) then
                        f = f + lam*chi(k, 1)*p(a, k)
                      else
                        f = f + 2.0_dblp*lam*p(a, k)*(chi(k, 1)*cos(m*lon) - chi(k, 2)*sin(m*lon))
                      endif
                    enddo
                  enddo
                  ! area element on the unit sphere: cos(lat) dlat dlon
                  s = s + wq(a)*wq(b)*half*hl*cos(lat)*f
                enddo
              enddo
              divint(i, j) = s
            enddo
          enddo

        end subroutine cell_divergence

      end program test_fvt_specflux

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
