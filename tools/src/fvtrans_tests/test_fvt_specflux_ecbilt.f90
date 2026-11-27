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
!      PROGRAM: [test_fvt_specflux_ecbilt]
!
!>     @brief Spectral convention of fvt_specflux_mod against ECBilt's transforms (ROADMAP step 1.1, wiso:1d:0001).
!
!      DESCRIPTION:
!>     Usage: test_fvt_specflux_ecbilt <path to coef.dat>
!!     Random spectral psi and chi with a red spectrum (fixed seed). Pass/fail check: the model winds u800, v800
!!     (ec_geowin) and udivg, vdivg (ec_divwin), computed here exactly as in atmdyn_mod with coef.dat and the model's
!!     FFT (NAG, nag.f), against the same winds synthesised with fvt_legendre_mod at the Gaussian points (T21).
!!     ECBilt-side test (split from test_fvt_specflux in wiso:1f:0002): it follows the model's transforms, so it
!!     changes with them (e.g. when the NAG FFT is replaced). The face-flux checks are in test_fvt_specflux.
!!     Stops with exit code 1 if the check fails.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      program test_fvt_specflux_ecbilt

        use global_constants_mod, only: ip, dblp=>dp, str_len
        use fvt_legendre_mod,     only: fvt_gauss_nodes, fvt_legendre_eval

        implicit none

        real(dblp), parameter :: TOL_ROUND = 1.0e-12_dblp

        character(len=str_len) :: coef_file
        integer(ip)            :: nfail

        nfail = 0
        if (command_argument_count() < 1) then
          write(*, '(a)') 'usage: test_fvt_specflux_ecbilt <path to coef.dat>'
          error stop 2
        endif
        call get_command_argument(1, coef_file)

        call check_convention(trim(coef_file))

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

      end program test_fvt_specflux_ecbilt

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
