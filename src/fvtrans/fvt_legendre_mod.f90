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
!      MODULE: [fvt_legendre_mod]
!
!>     @brief Gaussian quadrature and associated Legendre functions in the ECBilt spectral convention, at any truncation.
!
!      DESCRIPTION:
!>     Independent generator for what ECBilt reads from inputdata/coef.dat (T21 only), so that spectral fields (psi, chi)
!!     can be synthesised at latitudes other than the Gaussian ones (FV cell faces and corners) and at other truncations.
!!
!!     Convention (recovered from coef.dat and reproduced to round-off, see test_fvt_grid):
!!       - triangular truncation ntrunc, spectral index k ordered by m outer (0..ntrunc), n inner (m..ntrunc);
!!         nshm(m) = ntrunc+1-m, ll(k) = n;
!!       - P(n,m)(mu) = sqrt((2n+1)(n-m)!/(n+m)!) * P_n^m(mu), *without* the Condon-Shortley phase and without a
!!         sqrt(2) for m > 0: the mean of P(n,m)**2 over mu in [-1,1] is 1;
!!       - pd = dP/dmu (plain mu derivative, not (1-mu**2) dP/dmu);
!!       - pw = (w/2) * P, w the Gaussian weights (sum of w = 2);
!!       - latitude index 1 is the southernmost row (mu < 0).
!!     ECBilt multiplies pp, pd by sqrt(nlon) and pw by 1/sqrt(nlon) after reading (NAG FFT normalisation); the tables
!!     returned here are the file values, before that scaling.
!!
!!     Contract (public entry points):
!!       fvt_gauss_nodes(nlat,mu,w)              : Gaussian nodes (south to north) and weights (sum 2).
!!       fvt_spec_index(ntrunc,nshm,ll)          : ECBilt spectral index tables for triangular truncation ntrunc.
!!       fvt_legendre_eval(ntrunc,mu,p[,dpdmu])  : P(n,m) at any mu in [-1,1]; dP/dmu for |mu| < 1 only.
!!       fvt_legendre_tables(ntrunc,nlat,pp,pd,pw): coef.dat equivalent at the nlat Gaussian latitudes.
!!       fvt_nsh(ntrunc)                         : number of (n,m) pairs, (ntrunc+1)(ntrunc+2)/2.
!!
!>     Recurrences: diagonal P(m,m) = sqrt((2m+1)/(2m)) s P(m-1,m-1), s = sqrt((1-mu)(1+mu)); then the standard
!!     three-term recurrence in n for fully normalised functions (the normalisation constant differs from the usual
!!     4-pi one by a factor independent of n, so the coefficients are the same). Suitable up to truncations of a few
!!     hundred (s**m underflow near the poles is harmless: the true values are below the double precision range).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module fvt_legendre_mod

        use global_constants_mod, only: ip, dblp=>dp, LEG_PI=>pi_dp

        implicit none

        private

        public :: fvt_gauss_nodes, fvt_spec_index, fvt_legendre_eval, fvt_legendre_tables, fvt_nsh

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Module constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip), parameter :: LEG_NEWTON_MAX = 100                                           !< max Newton iterations per node
        real(dblp),  parameter :: LEG_NEWTON_TOL = 1.0e-15_dblp                                  !< tolerance on a Newton update

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Gaussian quadrature
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine fvt_gauss_nodes(nlat, mu, w)

          integer(ip), intent(in)  :: nlat       !< number of nodes (even)
          real(dblp),  intent(out) :: mu(nlat)   !< nodes, sin(latitude), ascending (south to north)
          real(dblp),  intent(out) :: w(nlat)    !< weights, sum = 2

          integer(ip) :: i, iter, k
          real(dblp)  :: x, dx, p0, p1, p2, dp

          if (nlat < 2 .or. mod(nlat, 2_ip) /= 0) error stop 'fvt_gauss_nodes: nlat must be even and >= 2'

          ! Newton on P_nlat for the northern half (x > 0), mirrored to the south: exact symmetry by construction
          do i = 1, nlat/2
            x = cos(LEG_PI*(real(i, dblp) - 0.25_dblp)/(real(nlat, dblp) + 0.5_dblp))
            do iter = 1, LEG_NEWTON_MAX
              call leg_pn(nlat, x, p0, p1)
              ! p1 = P_nlat(x), p0 = P_nlat-1(x)
              dp = real(nlat, dblp)*(p0 - x*p1)/((1.0_dblp - x)*(1.0_dblp + x))
              dx = p1/dp
              x  = x - dx
              if (abs(dx) <= LEG_NEWTON_TOL) exit
            enddo
            if (iter > LEG_NEWTON_MAX) error stop 'fvt_gauss_nodes: Newton iteration did not converge'
            ! weight from the converged node
            call leg_pn(nlat, x, p0, p1)
            dp = real(nlat, dblp)*(p0 - x*p1)/((1.0_dblp - x)*(1.0_dblp + x))
            p2 = 2.0_dblp/((1.0_dblp - x)*(1.0_dblp + x)*dp*dp)
            ! i = 1 is the node closest to the north pole
            k = nlat + 1 - i
            mu(k) = x
            w(k)  = p2
            mu(i) = -x
            w(i)  = p2
          enddo

        end subroutine fvt_gauss_nodes

        ! Unnormalised Legendre polynomials P_{n-1}(x) and P_n(x) by the Bonnet recurrence
        subroutine leg_pn(n, x, pnm1, pn)

          integer(ip), intent(in)  :: n
          real(dblp),  intent(in)  :: x
          real(dblp),  intent(out) :: pnm1, pn

          integer(ip) :: j
          real(dblp)  :: pm2

          pnm1 = 1.0_dblp
          pn   = x
          do j = 2, n
            pm2  = pnm1
            pnm1 = pn
            pn   = (real(2*j - 1, dblp)*x*pnm1 - real(j - 1, dblp)*pm2)/real(j, dblp)
          enddo

        end subroutine leg_pn

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Spectral indexing (ECBilt triangular convention)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        pure function fvt_nsh(ntrunc) result(nsh)

          integer(ip), intent(in) :: ntrunc
          integer(ip)             :: nsh

          nsh = ((ntrunc + 1)*(ntrunc + 2))/2

        end function fvt_nsh

        subroutine fvt_spec_index(ntrunc, nshm, ll)

          integer(ip), intent(in)  :: ntrunc
          integer(ip), intent(out) :: nshm(0:ntrunc)   !< number of n values for each m
          integer(ip), intent(out) :: ll(:)            !< total wavenumber n of index k (size fvt_nsh(ntrunc))

          integer(ip) :: k, m, n

          if (size(ll) /= fvt_nsh(ntrunc)) error stop 'fvt_spec_index: ll has the wrong size'

          k = 0
          do m = 0, ntrunc
            nshm(m) = ntrunc + 1 - m
            do n = m, ntrunc
              k = k + 1
              ll(k) = n
            enddo
          enddo

        end subroutine fvt_spec_index

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Associated Legendre functions
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine fvt_legendre_eval(ntrunc, mu, p, dpdmu)

          integer(ip),          intent(in)  :: ntrunc
          real(dblp),           intent(in)  :: mu(:)         !< points, sin(latitude), in [-1,1]
          real(dblp),           intent(out) :: p(:,:)        !< (size(mu), nsh) P(n,m)(mu)
          real(dblp), optional, intent(out) :: dpdmu(:,:)    !< (size(mu), nsh) dP(n,m)/dmu, requires |mu| < 1

          integer(ip) :: i, k, kdiag, m, n, npts
          real(dblp)  :: x, s, pmm, a, b, c, rn, rm

          npts = size(mu)
          if (size(p, 1) /= npts .or. size(p, 2) /= fvt_nsh(ntrunc)) error stop 'fvt_legendre_eval: p has the wrong shape'
          if (any(abs(mu) > 1.0_dblp)) error stop 'fvt_legendre_eval: |mu| > 1'
          if (present(dpdmu)) then
            if (size(dpdmu, 1) /= npts .or. size(dpdmu, 2) /= fvt_nsh(ntrunc)) &
              error stop 'fvt_legendre_eval: dpdmu has the wrong shape'
            ! dP/dmu is singular at the poles for m = 1: callers never need it there (zero-flux pole faces)
            if (any(abs(mu) >= 1.0_dblp)) error stop 'fvt_legendre_eval: dP/dmu requested at a pole'
          endif

          do i = 1, npts
            x   = mu(i)
            s   = sqrt((1.0_dblp - x)*(1.0_dblp + x))
            pmm = 1.0_dblp
            k   = 0
            do m = 0, ntrunc
              rm = real(m, dblp)
              if (m > 0) pmm = sqrt((2.0_dblp*rm + 1.0_dblp)/(2.0_dblp*rm))*s*pmm
              kdiag = k + 1
              ! n = m
              k = k + 1
              p(i, k) = pmm
              ! n = m+1
              if (m < ntrunc) then
                k = k + 1
                p(i, k) = sqrt(2.0_dblp*rm + 3.0_dblp)*x*pmm
              endif
              ! n >= m+2
              do n = m + 2, ntrunc
                rn = real(n, dblp)
                a  = sqrt((4.0_dblp*rn*rn - 1.0_dblp)/(rn*rn - rm*rm))
                b  = sqrt(((rn - 1.0_dblp)**2 - rm*rm)/(4.0_dblp*(rn - 1.0_dblp)**2 - 1.0_dblp))
                k  = k + 1
                p(i, k) = a*(x*p(i, k-1) - b*p(i, k-2))
              enddo
              ! derivative: (1-mu**2) dP(n,m)/dmu = -n mu P(n,m) + c(n,m) P(n-1,m)
              if (present(dpdmu)) then
                do n = m, ntrunc
                  rn = real(n, dblp)
                  if (n == m) then
                    dpdmu(i, kdiag) = -rn*x*p(i, kdiag)/(s*s)
                  else
                    c = sqrt((2.0_dblp*rn + 1.0_dblp)*(rn*rn - rm*rm)/(2.0_dblp*rn - 1.0_dblp))
                    dpdmu(i, kdiag+n-m) = (-rn*x*p(i, kdiag+n-m) + c*p(i, kdiag+n-m-1))/(s*s)
                  endif
                enddo
              endif
            enddo
          enddo

        end subroutine fvt_legendre_eval

        subroutine fvt_legendre_tables(ntrunc, nlat, pp, pd, pw)

          integer(ip), intent(in)  :: ntrunc
          integer(ip), intent(in)  :: nlat
          real(dblp),  intent(out) :: pp(:,:)   !< (nlat, nsh) P at the Gaussian latitudes
          real(dblp),  intent(out) :: pd(:,:)   !< (nlat, nsh) dP/dmu at the Gaussian latitudes
          real(dblp),  intent(out) :: pw(:,:)   !< (nlat, nsh) (w/2) P, Legendre analysis weights

          real(dblp)  :: mu(nlat), w(nlat)
          integer(ip) :: k

          call fvt_gauss_nodes(nlat, mu, w)
          call fvt_legendre_eval(ntrunc, mu, pp, pd)
          if (size(pw, 1) /= nlat .or. size(pw, 2) /= fvt_nsh(ntrunc)) error stop 'fvt_legendre_tables: pw has the wrong shape'
          do k = 1, fvt_nsh(ntrunc)
            pw(:, k) = 0.5_dblp*w(:)*pp(:, k)
          enddo

        end subroutine fvt_legendre_tables

      end module fvt_legendre_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
