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
!      PROGRAM: [test_fvt_grid]
!
!>     @brief Checks of fvt_legendre_mod and fvt_grid_mod (ROADMAP step 1.1, patch wiso:1:0001).
!
!      DESCRIPTION:
!>     1. T21 / 32 latitudes against inputdata/coef.dat (path given as first argument): Gaussian nodes and weights,
!!        spectral index tables, pp, pd, pw.
!!     2. For several truncations: orthonormality of P(n,m) under Gaussian quadrature, dP/dmu against a 4th-order
!!        finite difference, values at the poles.
!!     3. For the same grids: FV faces (exact poles and equator, bitwise symmetry, monotonicity, interlacing with the
!!        Gaussian nodes), cell areas against the Gaussian weights, total area 4 pi.
!!     Prints one line per check and stops with exit code 1 if any check fails.
!!     coef.dat is big-endian with 4-byte record markers: compile with -fconvert=big-endian (as the model).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      program test_fvt_grid

        use global_constants_mod, only: ip, dblp=>dp, str_len, PI=>pi_dp
        use fvt_legendre_mod,     only: fvt_gauss_nodes, fvt_spec_index, fvt_legendre_eval, fvt_legendre_tables, fvt_nsh
        use fvt_grid_mod,         only: fvt_grid_t, fvt_grid_init, fvt_grid_free

        implicit none

        integer(ip), parameter :: NTEST = 5
        integer(ip), parameter :: TRUNCS(NTEST) = [21, 42, 63, 106, 170]  !< triangular truncations tested
        integer(ip), parameter :: NLATS(NTEST)  = [32, 64, 96, 160, 256]  !< matching Gaussian grids
        integer(ip), parameter :: NLONS(NTEST)  = [64, 128, 192, 320, 512]

        character(len=str_len) :: coef_file
        integer(ip)            :: nfail, itest

        nfail = 0
        if (command_argument_count() < 1) then
          write(*, '(a)') 'usage: test_fvt_grid <path to coef.dat>'
          error stop 2
        endif
        call get_command_argument(1, coef_file)

        call check_coef_dat(trim(coef_file))
        do itest = 1, NTEST
          call check_legendre(TRUNCS(itest), NLATS(itest))
          call check_grid(NLATS(itest), NLONS(itest))
        enddo

        write(*, '(a)') repeat('-', 72)
        if (nfail > 0) then
          write(*, '(i0,a)') nfail, ' check(s) FAILED'
          error stop 1
        endif
        write(*, '(a)') 'all checks passed'

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Reporting
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine report(label, value, tol)

          character(len=*), intent(in) :: label
          real(dblp),       intent(in) :: value   !< measured error
          real(dblp),       intent(in) :: tol     !< pass if value <= tol

          character(len=4) :: verdict

          if (value <= tol) then
            verdict = 'ok'
          else
            verdict = 'FAIL'
            nfail = nfail + 1
          endif
          write(*, '(a4,1x,a,t62,es10.3,a,es8.1)') verdict, label, value, '  tol ', tol

        end subroutine report

        subroutine report_flag(label, ok)

          character(len=*), intent(in) :: label
          logical,          intent(in) :: ok

          if (ok) then
            write(*, '(a4,1x,a)') 'ok', label
          else
            write(*, '(a4,1x,a)') 'FAIL', label
            nfail = nfail + 1
          endif

        end subroutine report_flag

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   1. Against coef.dat (T21)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine check_coef_dat(path)

          character(len=*), intent(in) :: path

          integer(ip), parameter :: NT = 21, NL = 32
          integer(ip) :: nshm_f(0:NT), ll_f(fvt_nsh(NT)), nshm(0:NT), ll(fvt_nsh(NT))
          real(dblp)  :: pp_f(NL, fvt_nsh(NT)), pd_f(NL, fvt_nsh(NT)), pw_f(NL, fvt_nsh(NT))
          real(dblp)  :: pp(NL, fvt_nsh(NT)), pd(NL, fvt_nsh(NT)), pw(NL, fvt_nsh(NT))
          real(dblp)  :: mu(NL), w(NL)
          integer(ip) :: u, ios, k
          real(dblp)  :: epp, epd, epw

          write(*, '(a)') repeat('-', 72)
          write(*, '(a)') 'T21 / 32 latitudes against '//path
          open(newunit=u, file=path, form='unformatted', access='sequential', status='old', action='read', iostat=ios)
          if (ios /= 0) then
            call report_flag('open coef.dat', .false.)
            return
          endif
          read(u) nshm_f, ll_f
          read(u) pp_f
          read(u) pd_f
          read(u) pw_f
          close(u)

          call fvt_spec_index(NT, nshm, ll)
          call report_flag('spectral index tables nshm, ll identical', all(nshm == nshm_f) .and. all(ll == ll_f))

          ! nodes and weights: mu = P(1,0)/sqrt(3), w = 2 pw(:,1) since P(0,0) = 1
          call fvt_gauss_nodes(NL, mu, w)
          call report('Gaussian nodes   max |mu - pp(:,2)/sqrt(3)|', maxval(abs(mu - pp_f(:, 2)/sqrt(3.0_dblp))), 1.0e-14_dblp)
          call report('Gaussian weights max |w - 2 pw(:,1)|', maxval(abs(w - 2.0_dblp*pw_f(:, 1))), 1.0e-14_dblp)

          call fvt_legendre_tables(NT, NL, pp, pd, pw)
          ! errors relative to the largest value of each function (pd reaches about 1e3 near the poles)
          epp = 0.0_dblp
          epw = 0.0_dblp
          do k = 1, fvt_nsh(NT)
            epp = max(epp, maxval(abs(pp(:, k) - pp_f(:, k)))/maxval(abs(pp_f(:, k))))
            epw = max(epw, maxval(abs(pw(:, k) - pw_f(:, k)))/maxval(abs(pw_f(:, k))))
          enddo
          call report('pp  max |P - pp| / max|pp|  (per function)', epp, 1.0e-12_dblp)
          ! pd(:,1) is identically zero (P(0,0) = 1): absolute error for k = 1, relative for k > 1
          epd = epd_nonzero(pd, pd_f)
          call report('pd  max |dP/dmu - pd| / max|pd|  (per function)', epd, 1.0e-12_dblp)
          call report('pw  max |(w/2)P - pw| / max|pw|  (per function)', epw, 1.0e-12_dblp)

        end subroutine check_coef_dat

        function epd_nonzero(a, b) result(err)

          real(dblp), intent(in) :: a(:,:)
          real(dblp), intent(in) :: b(:,:)
          real(dblp)             :: err

          integer(ip) :: k

          err = maxval(abs(a(:, 1) - b(:, 1)))
          do k = 2, size(a, 2)
            err = max(err, maxval(abs(a(:, k) - b(:, k)))/maxval(abs(b(:, k))))
          enddo

        end function epd_nonzero

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   2. Legendre functions at any truncation
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine check_legendre(nt, nl)

          integer(ip), intent(in) :: nt
          integer(ip), intent(in) :: nl

          real(dblp), allocatable :: pp(:,:), pd(:,:), pw(:,:)
          real(dblp), allocatable :: xm(:), pm2(:,:), pm1(:,:), pp1(:,:), pp2(:,:), fd(:,:), pp_pole(:,:)
          integer(ip), allocatable :: nshm(:), ll(:)
          integer(ip) :: nsh, k1, k2, m, ka, kb, i, k
          real(dblp)  :: eorth, s, ederiv, epole, h, scale
          character(len=40) :: tag

          nsh = fvt_nsh(nt)
          write(tag, '(a,i0,a,i0,a)') 'T', nt, ' / ', nl, ' lat'
          write(*, '(a)') repeat('-', 72)
          write(*, '(a)') trim(tag)//': Legendre functions'

          allocate(pp(nl, nsh), pd(nl, nsh), pw(nl, nsh), nshm(0:nt), ll(nsh))
          call fvt_spec_index(nt, nshm, ll)
          call fvt_legendre_tables(nt, nl, pp, pd, pw)

          ! orthonormality: sum_i (w_i/2) P_a P_b = delta_ab for equal m (exact quadrature, degree <= 2 nt < 2 nl)
          eorth = 0.0_dblp
          k2 = 0
          do m = 0, nt
            k1 = k2 + 1
            k2 = k2 + nshm(m)
            do ka = k1, k2
              do kb = ka, k2
                s = sum(pw(:, ka)*pp(:, kb))
                if (ka == kb) s = s - 1.0_dblp
                eorth = max(eorth, abs(s))
              enddo
            enddo
          enddo
          call report('orthonormality max |<P_a,P_b> - delta_ab|', eorth, 1.0e-12_dblp)

          ! dP/dmu against a 4th-order centred difference, at the Gaussian nodes (|mu| < 1)
          h = 1.0e-4_dblp/real(nt, dblp)
          allocate(xm(nl), pm2(nl, nsh), pm1(nl, nsh), pp1(nl, nsh), pp2(nl, nsh), fd(nl, nsh))
          xm = pp(:, 2)/sqrt(3.0_dblp)
          call fvt_legendre_eval(nt, xm - 2.0_dblp*h, pm2)
          call fvt_legendre_eval(nt, xm - h, pm1)
          call fvt_legendre_eval(nt, xm + h, pp1)
          call fvt_legendre_eval(nt, xm + 2.0_dblp*h, pp2)
          fd = (8.0_dblp*(pp1 - pm1) - (pp2 - pm2))/(12.0_dblp*h)
          ederiv = 0.0_dblp
          do k = 2, nsh
            scale  = maxval(abs(pd(:, k)))
            ederiv = max(ederiv, maxval(abs(pd(:, k) - fd(:, k)))/scale)
          enddo
          call report('dP/dmu vs 4th-order finite difference (rel.)', ederiv, 1.0e-7_dblp)

          ! poles: P(n,0)(+-1) = sqrt(2n+1) (+-1)**n, P(n,m>0)(+-1) = 0 exactly
          allocate(pp_pole(2, nsh))
          call fvt_legendre_eval(nt, [-1.0_dblp, 1.0_dblp], pp_pole)
          epole = 0.0_dblp
          do k = 1, nsh
            if (k <= nshm(0)) then
              epole = max(epole, abs(pp_pole(2, k) - sqrt(real(2*ll(k) + 1, dblp))))
              epole = max(epole, abs(pp_pole(1, k) - sqrt(real(2*ll(k) + 1, dblp))*real((-1)**ll(k), dblp)))
            else
              do i = 1, 2
                epole = max(epole, abs(pp_pole(i, k)))
              enddo
            endif
          enddo
          call report('pole values max error', epole, 1.0e-12_dblp)

        end subroutine check_legendre

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   3. FV grid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine check_grid(nl, nlon)

          integer(ip), intent(in) :: nl
          integer(ip), intent(in) :: nlon

          type(fvt_grid_t) :: g
          integer(ip)      :: i
          logical          :: ok
          character(len=40) :: tag

          write(tag, '(i0,a,i0)') nl, ' x ', nlon
          write(*, '(a)') 'grid '//trim(tag)//': FV cells'
          call fvt_grid_init(g, nl, nlon)

          ok = g%mu_face(0) == -1.0_dblp .and. g%mu_face(nl) == 1.0_dblp .and. g%mu_face(nl/2) == 0.0_dblp
          call report_flag('faces exact at the poles and equator', ok)
          ok = .true.
          do i = 0, nl
            ok = ok .and. g%mu_face(nl-i) == -g%mu_face(i)
          enddo
          call report_flag('faces symmetric (bitwise)', ok)
          ok = g%cos_face(0) == 0.0_dblp .and. g%cos_face(nl) == 0.0_dblp
          call report_flag('cos(latitude) zero at the pole faces', ok)
          ok = .true.
          do i = 1, nl
            ok = ok .and. g%mu_face(i-1) < g%mu(i) .and. g%mu(i) < g%mu_face(i)
          enddo
          call report_flag('each Gaussian node strictly inside its cell', ok)
          call report('max |dmu - w| (FV area vs Gaussian weight)', maxval(abs(g%dmu - g%wgt)), 1.0e-14_dblp)
          call report('|sum(area)*nlon - 4 pi|', abs(sum(g%area)*real(nlon, dblp) - 4.0_dblp*PI), 1.0e-13_dblp)
          call report('|lon_face(nlon) - lon_face(0) - 2 pi|', abs(g%lon_face(nlon) - g%lon_face(0) - 2.0_dblp*PI), &
                      1.0e-13_dblp)

          call fvt_grid_free(g)

        end subroutine check_grid

      end program test_fvt_grid

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
