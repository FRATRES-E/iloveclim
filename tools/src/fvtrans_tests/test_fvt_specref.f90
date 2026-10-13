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
!      PROGRAM: [test_fvt_specref]
!
!>     @brief Checks of the spectral baseline fvt_specref_mod (ROADMAP step 1.1, patch wiso:1d:0002).
!
!      DESCRIPTION:
!>     On a field that is a single spherical harmonic Y(n=7,m=3) (exactly representable at T21), plus a constant:
!!       1. grid -> spectral -> grid round trip is exact (ec_rggtosp / ec_sptogg);
!!       2. ec_trafluxdiv under solid-body rotation about the polar axis (u = u0 cos(lat), v = 0, div = 0) gives
!!          u0 dq/dlon, i.e. -u0 m (c2 cos(m lon) + c1 sin(m lon)) P (the analytic derivative of the field);
!!       3. ec_hdiff gives -n(n+1)/tdifq times the harmonic (the constant is not diffused);
!!       4. one full spr_step conserves the Gaussian-weighted total (no clipping occurs here).
!!     Also at T42 on 64 x 128 for 1-3. Stops with exit code 1 if a check fails.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      program test_fvt_specref

        use global_constants_mod, only: ip, dblp=>dp
        use fvt_legendre_mod,     only: fvt_legendre_eval, fvt_nsh
        use fvt_grid_mod,         only: fvt_grid_t, fvt_grid_init
        use fvt_specref_mod,      only: spr_t, spr_init, spr_step, spr_trafluxdiv, spr_hdiff, spr_rggtosp, spr_sptogg

        implicit none

        integer(ip) :: nfail

        nfail = 0
        call check(32_ip, 64_ip)
        call check(64_ip, 128_ip)
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

        subroutine check(nl, nlo)

          integer(ip), intent(in) :: nl
          integer(ip), intent(in) :: nlo

          integer(ip), parameter :: NN = 7, MM = 3
          real(dblp),  parameter :: C1 = 0.3_dblp, C2 = -0.2_dblp, U0 = 0.5_dblp, TDIFQ = 10.0_dblp, BG = 2.0_dblp
          type(fvt_grid_t) :: g
          type(spr_t)      :: sr
          real(dblp), allocatable :: p(:,:), y(:,:), dy(:,:), q(:,:), q2(:,:), u(:,:), v(:,:), dv(:,:), tend(:,:)
          real(dblp), allocatable :: qs(:,:)
          integer(ip) :: i, j, k, kk, m, n
          real(dblp)  :: lon, m0, m1
          character(len=40) :: tag

          call fvt_grid_init(g, nl, nlo)
          call spr_init(sr, g)
          write(tag, '(a,i0,a,i0,a,i0)') 'T', sr%nt, ', ', nl, ' x ', nlo
          write(*, '(a)') repeat('-', 100)
          write(*, '(a)') 'Spectral baseline, '//trim(tag)
          allocate(p(nl, fvt_nsh(sr%nt)), y(nl, nlo), dy(nl, nlo), q(nl, nlo), q2(nl, nlo), u(nl, nlo), v(nl, nlo))
          allocate(dv(nl, nlo), tend(nl, nlo), qs(fvt_nsh(sr%nt), 2))
          call fvt_legendre_eval(sr%nt, g%mu, p)
          ! index of (n, m)
          kk = 0
          k = 0
          do m = 0, sr%nt
            do n = m, sr%nt
              k = k + 1
              if (n == NN .and. m == MM) kk = k
            enddo
          enddo
          ! y = 2 P (c1 cos(m lon) - c2 sin(m lon)) (ECBilt convention) and its longitude derivative
          do j = 1, nlo
            lon = g%lon(j)
            do i = 1, nl
              y(i, j)  = 2.0_dblp*p(i, kk)*(C1*cos(MM*lon) - C2*sin(MM*lon))
              dy(i, j) = 2.0_dblp*p(i, kk)*real(MM, dblp)*(-C1*sin(MM*lon) - C2*cos(MM*lon))
            enddo
          enddo
          q = BG + y

          ! 1. round trip
          call spr_rggtosp(sr, q, qs)
          call spr_sptogg(sr, qs, q2, sr%pp)
          call report('round trip: max |q - sptogg(rggtosp(q))| / max |q|', maxval(abs(q2 - q))/maxval(abs(q)), 1.0e-13_dblp)

          ! 2. advection tendency under rotation about the polar axis
          do i = 1, nl
            u(i, :) = U0*cos(g%lat(i))
          enddo
          v  = 0.0_dblp
          dv = 0.0_dblp
          call spr_trafluxdiv(sr, qs, q, u, v, dv, tend)
          call report('ec_trafluxdiv: max |tendency - u0 dq/dlon| / max |u0 dq/dlon|', &
                      maxval(abs(tend - U0*dy))/maxval(abs(U0*dy)), 1.0e-12_dblp)

          ! 3. diffusion tendency (round-off leaking into the high modes is amplified by n(n+1), up to ~1800 at T42)
          call spr_hdiff(sr, qs, TDIFQ, tend)
          call report('ec_hdiff: max |tendency + n(n+1)/tdifq y| / max |n(n+1)/tdifq y|', &
                      maxval(abs(tend + real(NN*(NN + 1), dblp)/TDIFQ*y))/maxval(abs(real(NN*(NN + 1), dblp)/TDIFQ*y)), &
                      1.0e-11_dblp)

          ! 4. a full step conserves the Gaussian-weighted total
          m0 = 0.0_dblp
          do i = 1, nl
            m0 = m0 + g%wgt(i)*sum(q(i, :))
          enddo
          q2 = q
          call spr_step(sr, q2, u, v, dv, 0.05_dblp, TDIFQ)
          m1 = 0.0_dblp
          do i = 1, nl
            m1 = m1 + g%wgt(i)*sum(q2(i, :))
          enddo
          call report('spr_step: relative change of the Gaussian-weighted total', abs(m1 - m0)/m0, 1.0e-13_dblp)

        end subroutine check

      end program test_fvt_specref

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
