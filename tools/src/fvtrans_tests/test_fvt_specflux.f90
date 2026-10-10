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
!      PROGRAM: [test_fvt_specflux]
!
!>     @brief Checks of fvt_specflux_mod (ROADMAP step 1.1, patch wiso:1d:0001).
!
!      DESCRIPTION:
!>     Usage: test_fvt_specflux <path to coef.dat>
!!     Random spectral psi and chi with a red spectrum (fixed seed). Pass/fail checks:
!!       1. the spectral convention of fvt_specflux_mod reproduces ECBilt's transforms: the model winds u800, v800
!!          (ec_geowin) and udivg, vdivg (ec_divwin), computed here exactly as in atmdyn_mod with coef.dat and the
!!          NAG FFT of nag.f, against the same winds synthesised with fvt_legendre_mod at the Gaussian points;
!!       2. rotational fluxes are discretely non-divergent (each cell sums to zero up to round-off);
!!       3. divergent fluxes: the outflow of each cell equals the integral of the divergence over the cell
!!          (Gauss theorem; divergence integrated independently by 16 x 16 quadrature);
!!       4. the row integrals have converged in the quadrature (default against 2x the points);
!!       5. pole faces carry no flux;
!!       6. a uniform field transported by fvt_ffsl with the rotational fluxes stays uniform.
!!     The same checks 2-6 at T42 on a 64 x 128 grid (no coef.dat needed there).
!!     Stops with exit code 1 if a check fails.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      program test_fvt_specflux

        use global_constants_mod, only: ip, dblp=>dp, str_len
        use fvt_legendre_mod,     only: fvt_gauss_nodes, fvt_spec_index, fvt_legendre_eval, fvt_nsh
        use fvt_grid_mod,         only: fvt_grid_t, fvt_grid_init
        use fvt_specflux_mod,     only: fvt_specflux_t, fvt_specflux_init, fvt_specflux_fluxes
        use fvt_ffsl_mod,         only: fvt_ffsl_density, FVT_RECON_PPM, FVT_LIM_POSDEF

        implicit none

        real(dblp), parameter :: TOL_ROUND = 1.0e-12_dblp

        character(len=str_len) :: coef_file
        integer(ip)            :: nfail

        nfail = 0
        if (command_argument_count() < 1) then
          write(*, '(a)') 'usage: test_fvt_specflux <path to coef.dat>'
          error stop 2
        endif
        call get_command_argument(1, coef_file)

        call check_convention(trim(coef_file))
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
! dmr&clo   1. Spectral convention against the ECBilt transforms (T21)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine check_convention(path)

          character(len=*), intent(in) :: path

          integer(ip), parameter :: NT = 21, NL = 32, NLO = 64
          integer(ip), parameter :: NS = (NT + 1)*(NT + 2)/2
          real(dblp),  parameter :: RADIUS = 6.37e6_dblp
          real(dblp),  parameter :: OM = 4.0_dblp*3.14159265358979323846_dblp/(24.0_dblp*3600.0_dblp)
          integer(ip) :: nshm_f(0:NT), ll_f(NS), u, ios, i, j, k, m, n
          real(dblp)  :: pp(NL, NS), pd(NL, NS), pw(NL, NS)
          real(dblp)  :: psi(NS, 2), chi(NS, 2), vv(NS, 2)
          real(dblp)  :: dpsdl(NL, NLO), dpsdm(NL, NLO), dxdl(NL, NLO), dxdm(NL, NLO)
          real(dblp)  :: um(NL, NLO), vm(NL, NLO), ud(NL, NLO), vd(NL, NLO)
          real(dblp)  :: u2(NL, NLO), v2(NL, NLO), ud2(NL, NLO), vd2(NL, NLO)
          real(dblp)  :: trigi(2*NLO), wgg(NL*NLO), dummy(NL, NLO)
          real(dblp)  :: p(NL, NS), dp(NL, NS), mu(NL), w(NL), cosfi(NL), lon, sq, fa, fb, fda, fdb, gl
          integer     :: ifail

          write(*, '(a)') repeat('-', 100)
          write(*, '(a)') '1. Spectral convention against ec_geowin / ec_divwin (T21, coef.dat, NAG FFT)'
          open(newunit=u, file=path, form='unformatted', access='sequential', status='old', action='read', iostat=ios)
          if (ios /= 0) error stop 'test_fvt_specflux: cannot open coef.dat'
          read(u) nshm_f, ll_f
          read(u) pp
          read(u) pd
          read(u) pw
          close(u)
          ! ECBilt's NAG normalisation (ec_initdyn)
          sq = sqrt(real(NLO, dblp))
          pp = pp*sq
          pd = pd*sq
          dummy = 1.0_dblp
          ifail = 0
          call c06fqf(int(NL), int(NLO), dummy, 'i', trigi, wgg, ifail)

          call random_coefficients(NT, 11_ip, psi)
          call random_coefficients(NT, 23_ip, chi)

          ! model side, as ec_geowin and ec_divwin (chi is used directly, as x in ec_divwin)
          call model_ddl(NT, psi, vv)
          call model_sptogg(vv, dpsdl, pp)
          call model_sptogg(psi, dpsdm, pd)
          call model_ddl(NT, chi, vv)
          call model_sptogg(vv, dxdl, pp)
          call model_sptogg(chi, dxdm, pd)
          call fvt_gauss_nodes(NL, mu, w)
          cosfi = sqrt((1.0_dblp - mu)*(1.0_dblp + mu))
          do j = 1, NLO
            do i = 1, NL
              um(i, j) = -RADIUS*OM*dpsdm(i, j)*cosfi(i)
              vm(i, j) =  RADIUS*OM*dpsdl(i, j)/cosfi(i)
              ud(i, j) = dxdl(i, j)/(RADIUS*cosfi(i))
              vd(i, j) = dxdm(i, j)*cosfi(i)/RADIUS
            enddo
          enddo

          ! fvt side: the convention documented in fvt_specflux_mod, synthesised point by point
          call fvt_legendre_eval(NT, mu, p, dp)
          do j = 1, NLO
            lon = real(j - 1, dblp)*2.0_dblp*3.14159265358979323846_dblp/real(NLO, dblp)
            do i = 1, NL
              fa = 0.0_dblp
              fb = 0.0_dblp
              fda = 0.0_dblp
              fdb = 0.0_dblp
              k = 0
              do m = 0, NT
                do n = m, NT
                  k = k + 1
                  if (m == 0) then
                    fb = fb + psi(k, 1)*dp(i, k)
                    fdb = fdb + chi(k, 1)*dp(i, k)
                  else
                    ! d/dlon of (c1 cos - c2 sin) and d/dmu of the same
                    fa  = fa  + 2.0_dblp*real(m, dblp)*p(i, k)*(-psi(k, 1)*sin(m*lon) - psi(k, 2)*cos(m*lon))
                    fb  = fb  + 2.0_dblp*dp(i, k)*(psi(k, 1)*cos(m*lon) - psi(k, 2)*sin(m*lon))
                    fda = fda + 2.0_dblp*real(m, dblp)*p(i, k)*(-chi(k, 1)*sin(m*lon) - chi(k, 2)*cos(m*lon))
                    fdb = fdb + 2.0_dblp*dp(i, k)*(chi(k, 1)*cos(m*lon) - chi(k, 2)*sin(m*lon))
                  endif
                enddo
              enddo
              u2(i, j)  = -RADIUS*OM*fb*cosfi(i)
              v2(i, j)  =  RADIUS*OM*fa/cosfi(i)
              ud2(i, j) = fda/(RADIUS*cosfi(i))
              vd2(i, j) = fdb*cosfi(i)/RADIUS
            enddo
          enddo
          gl = maxval(abs(um))
          call report('u800  (rotational, from psi): max |fvt - model| / max |model|', maxval(abs(u2 - um))/gl, TOL_ROUND)
          gl = maxval(abs(vm))
          call report('v800  (rotational, from psi): max |fvt - model| / max |model|', maxval(abs(v2 - vm))/gl, TOL_ROUND)
          gl = maxval(abs(ud))
          call report('udivg (divergent, from chi):  max |fvt - model| / max |model|', maxval(abs(ud2 - ud))/gl, TOL_ROUND)
          gl = maxval(abs(vd))
          call report('vdivg (divergent, from chi):  max |fvt - model| / max |model|', maxval(abs(vd2 - vd))/gl, TOL_ROUND)

        end subroutine check_convention

        ! ec_ddl: zonal derivative in spectral space
        subroutine model_ddl(nt, as, dadl)

          integer(ip), intent(in)  :: nt
          real(dblp),  intent(in)  :: as(:,:)
          real(dblp),  intent(out) :: dadl(:,:)

          integer(ip) :: k, m, n

          k = 0
          do m = 0, nt
            do n = m, nt
              k = k + 1
              dadl(k, 1) = -real(m, dblp)*as(k, 2)
              dadl(k, 2) =  real(m, dblp)*as(k, 1)
            enddo
          enddo

        end subroutine model_ddl

        ! ec_sptogg, verbatim algorithm (inverse Legendre transform, then NAG c06fqf)
        subroutine model_sptogg(as, agg, pploc)

          real(dblp), intent(in)  :: as(:,:)
          real(dblp), intent(out) :: agg(:,:)
          real(dblp), intent(in)  :: pploc(:,:)

          integer(ip), parameter :: NT = 21
          integer(ip) :: nlat, nlon, k, k1, k2, m, mr, mi
          real(dblp)  :: trigi(2*size(agg, 2)), wgg(size(agg))
          real(dblp)  :: dummy(size(agg, 1), size(agg, 2))
          integer     :: ifail

          nlat = size(agg, 1)
          nlon = size(agg, 2)
          agg = 0.0_dblp
          k2 = NT + 1
          do k = 1, k2
            agg(:, 1) = agg(:, 1) + as(k, 1)*pploc(:, k)
          enddo
          do m = 1, NT
            mr = m + 1
            mi = nlon + 1 - m
            k1 = k2 + 1
            k2 = k2 + NT + 1 - m
            do k = k1, k2
              agg(:, mr) = agg(:, mr) + as(k, 1)*pploc(:, k)
              agg(:, mi) = agg(:, mi) - as(k, 2)*pploc(:, k)
            enddo
          enddo
          dummy = 1.0_dblp
          ifail = 0
          call c06fqf(int(nlat), int(nlon), dummy, 'i', trigi, wgg, ifail)
          call c06fqf(int(nlat), int(nlon), agg, 'r', trigi, wgg, ifail)

        end subroutine model_sptogg

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   2-6. Face fluxes
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
          write(*, '(a)') '2-6. Face fluxes, '//trim(tag)
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
