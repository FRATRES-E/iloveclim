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
!      MODULE: [fvt_testcases_mod]
!
!>     @brief Standard 2-D transport test cases on the sphere for the FV transport (flows, initial fields, diagnostics).
!
!      DESCRIPTION:
!>     Unit sphere (R = 1), time in days, period TC_PERIOD = 12 days for every flow.
!!     Flows:
!!       TC_FLOW_SBR       solid-body rotation, Williamson et al. (1992) case 1, angle alpha;
!!       TC_FLOW_DEFORM    non-divergent deformational flow, Lauritzen et al. (2012) eqs. 18-20;
!!       TC_FLOW_DIVERGENT divergent deformational flow, Lauritzen et al. (2012) eqs. 21-22.
!!     The non-divergent flows give face fluxes from the streamfunction at the cell corners (discretely
!!     non-divergent); the divergent flow from Gauss-Legendre quadrature of the winds along each face.
!!     Initial fields (cell means by 6 x 6 Gauss-Legendre quadrature in mu and longitude):
!!       TC_IC_WG_BELL      Williamson cosine bell, radius 1/3, centre (3 pi/2, 0), height 1;
!!       TC_IC_HILLS        Lauritzen Gaussian hills (h_max 0.95, b 5), centres (5 pi/6, 0), (7 pi/6, 0);
!!       TC_IC_BELLS        Lauritzen cosine bells, radius 1/2, background 0.1, amplitude 0.9;
!!       TC_IC_CYLINDERS    Lauritzen slotted cylinders, radius 1/2, values 0.1 and 1;
!!       TC_IC_CORRELATED   correlated cosine bells, -0.8 bells**2 + 0.9 (Lauritzen and Thuburn 2012).
!!     Diagnostics: error norms (l1 as Williamson, l2, linf, phi_min, phi_max as Lauritzen), filament diagnostic l_f,
!!     mixing diagnostics l_r, l_u, l_o (Lauritzen and Thuburn 2012, as summarised in Lauritzen et al. 2012, app. C;
!!     the closest point on the curve is found numerically on the initial range instead of the closed-form root).
!!
!!     Contract (public entry points):
!!       tc_fluxes(g,flow,alpha,t,dt,fx,fy)      : face volume fluxes over dt with the winds at time t.
!!       tc_cell_means(g,ic,phi)                 : cell means of an initial field.
!!       tc_norms(g,phi,phit,dphi0) result(nrm)  : error norms of phi against the exact phit.
!!       tc_filament(g,phi,phi0,tau,lf)          : l_f(tau) in percent.
!!       tc_mixing(g,chi,xi,lr,lu,lo)            : mixing diagnostics of the pair (cosine bells, correlated bells).
!!       tc_total(g,phi)                         : area integral of phi over the unit sphere.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module fvt_testcases_mod

        use global_constants_mod, only: ip, dblp=>dp, PI=>pi_dp
        use fvt_legendre_mod,     only: fvt_gauss_nodes
        use fvt_grid_mod,         only: fvt_grid_t

        implicit none

        private

        public :: TC_PERIOD, TC_FLOW_SBR, TC_FLOW_DEFORM, TC_FLOW_DIVERGENT
        public :: TC_IC_WG_BELL, TC_IC_HILLS, TC_IC_BELLS, TC_IC_CYLINDERS, TC_IC_CORRELATED
        public :: tc_norms_t, tc_fluxes, tc_cell_means, tc_norms, tc_filament, tc_mixing, tc_total

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Module constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp),  parameter :: TC_PERIOD = 12.0_dblp           !< period of every flow (days)

        integer(ip), parameter :: TC_FLOW_SBR       = 1
        integer(ip), parameter :: TC_FLOW_DEFORM    = 2
        integer(ip), parameter :: TC_FLOW_DIVERGENT = 3

        integer(ip), parameter :: TC_IC_WG_BELL    = 1
        integer(ip), parameter :: TC_IC_HILLS      = 2
        integer(ip), parameter :: TC_IC_BELLS      = 3
        integer(ip), parameter :: TC_IC_CYLINDERS  = 4
        integer(ip), parameter :: TC_IC_CORRELATED = 5

        real(dblp),  parameter :: TC_TAU_TOL = 1.0e-10_dblp       !< tolerance of the filament thresholds
        integer(ip), parameter :: TC_NQ_CELL = 6                  !< quadrature points per direction, cell means
        integer(ip), parameter :: TC_NQ_FACE = 8                  !< quadrature points along a face (divergent flow)

        ! Lauritzen et al. (2012) initial conditions
        real(dblp),  parameter :: LON1 = 5.0_dblp*PI/6.0_dblp, LON2 = 7.0_dblp*PI/6.0_dblp
        real(dblp),  parameter :: HILL_HMAX = 0.95_dblp, HILL_B = 5.0_dblp
        real(dblp),  parameter :: BELL_R = 0.5_dblp, BELL_BG = 0.1_dblp, BELL_AMP = 0.9_dblp
        real(dblp),  parameter :: CYL_BG = 0.1_dblp, CYL_C = 1.0_dblp
        real(dblp),  parameter :: CORR_A = -0.8_dblp, CORR_B = 0.9_dblp

        ! mixing diagnostics: initial ranges of chi (bells) and xi (correlated bells)
        real(dblp),  parameter :: MIX_CHI_MIN = 0.1_dblp, MIX_CHI_MAX = 1.0_dblp
        real(dblp),  parameter :: MIX_XI_MIN = CORR_A*MIX_CHI_MAX**2 + CORR_B
        real(dblp),  parameter :: MIX_XI_MAX = CORR_A*MIX_CHI_MIN**2 + CORR_B

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Types
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        type :: tc_norms_t
          real(dblp) :: l1     = 0.0_dblp   !< sum |phi - phit| A / sum |phit| A
          real(dblp) :: l2     = 0.0_dblp   !< sqrt(sum (phi - phit)**2 A / sum phit**2 A)
          real(dblp) :: linf   = 0.0_dblp   !< max |phi - phit| / max |phit|
          real(dblp) :: phimin = 0.0_dblp   !< (min phi - min phit) / (max phi0 - min phi0)
          real(dblp) :: phimax = 0.0_dblp   !< (max phi - max phit) / (max phi0 - min phi0)
        end type tc_norms_t

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Flows
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine tc_fluxes(g, flow, alpha, t, dt, fx, fy)

          type(fvt_grid_t), intent(in)  :: g
          integer(ip),      intent(in)  :: flow
          real(dblp),       intent(in)  :: alpha       !< rotation angle (TC_FLOW_SBR only)
          real(dblp),       intent(in)  :: t           !< time at which the winds are evaluated (days)
          real(dblp),       intent(in)  :: dt          !< step (days)
          real(dblp),       intent(out) :: fx(:,:)     !< (nlat,nlon)
          real(dblp),       intent(out) :: fy(0:,:)    !< (0:nlat,nlon)

          real(dblp)  :: psic(0:g%nlat, 0:g%nlon), xq(TC_NQ_FACE), wq(TC_NQ_FACE), s, x
          integer(ip) :: i, j, k

          select case (flow)
            case (TC_FLOW_SBR, TC_FLOW_DEFORM)
              do j = 0, g%nlon
                do i = 0, g%nlat
                  psic(i, j) = psi_point(flow, alpha, g%lon_face(j), g%mu_face(i), g%cos_face(i), t)
                enddo
              enddo
              do j = 1, g%nlon
                do i = 1, g%nlat
                  fx(i, j) = -dt*(psic(i, j) - psic(i-1, j))
                enddo
                do i = 0, g%nlat
                  fy(i, j) = dt*(psic(i, j) - psic(i, j-1))
                enddo
              enddo
            case (TC_FLOW_DIVERGENT)
              call fvt_gauss_nodes(TC_NQ_FACE, xq, wq)
              do j = 1, g%nlon
                do i = 1, g%nlat
                  ! integral of u dlat along the east face
                  s = 0.0_dblp
                  do k = 1, TC_NQ_FACE
                    x = 0.5_dblp*(g%lat_face(i-1) + g%lat_face(i)) + 0.5_dblp*(g%lat_face(i) - g%lat_face(i-1))*xq(k)
                    s = s + wq(k)*u_div(g%lon_face(j), x, t)
                  enddo
                  fx(i, j) = dt*0.5_dblp*(g%lat_face(i) - g%lat_face(i-1))*s
                enddo
                do i = 0, g%nlat
                  ! integral of v cos(lat) dlon along the north face
                  s = 0.0_dblp
                  do k = 1, TC_NQ_FACE
                    x = 0.5_dblp*(g%lon_face(j-1) + g%lon_face(j)) + 0.5_dblp*g%dlon*xq(k)
                    s = s + wq(k)*v_div(x, g%cos_face(i), t)
                  enddo
                  fy(i, j) = dt*g%cos_face(i)*0.5_dblp*g%dlon*s
                enddo
              enddo
            case default
              error stop 'tc_fluxes: unknown flow'
          end select
          fy(0, :)       = 0.0_dblp
          fy(g%nlat, :)  = 0.0_dblp

        end subroutine tc_fluxes

        ! Streamfunction (unit sphere) with u = -d(psi)/d(lat), v = d(psi)/d(lon) / cos(lat)
        pure function psi_point(flow, alpha, lon, mu, coslat, t) result(psi)

          integer(ip), intent(in) :: flow
          real(dblp),  intent(in) :: alpha
          real(dblp),  intent(in) :: lon
          real(dblp),  intent(in) :: mu
          real(dblp),  intent(in) :: coslat
          real(dblp),  intent(in) :: t
          real(dblp)              :: psi

          real(dblp) :: lonp

          if (flow == TC_FLOW_SBR) then
            psi = -(2.0_dblp*PI/TC_PERIOD)*(mu*cos(alpha) - cos(lon)*coslat*sin(alpha))
          else
            lonp = lon - 2.0_dblp*PI*t/TC_PERIOD
            psi  = (10.0_dblp/TC_PERIOD)*sin(lonp)**2*coslat**2*cos(PI*t/TC_PERIOD) - (2.0_dblp*PI/TC_PERIOD)*mu
          endif

        end function psi_point

        pure function u_div(lon, lat, t) result(u)

          real(dblp), intent(in) :: lon
          real(dblp), intent(in) :: lat
          real(dblp), intent(in) :: t
          real(dblp)             :: u

          real(dblp) :: lonp

          lonp = lon - 2.0_dblp*PI*t/TC_PERIOD
          u = -(5.0_dblp/TC_PERIOD)*sin(0.5_dblp*lonp)**2*sin(2.0_dblp*lat)*cos(lat)**2*cos(PI*t/TC_PERIOD) &
              + (2.0_dblp*PI/TC_PERIOD)*cos(lat)

        end function u_div

        pure function v_div(lon, coslat, t) result(v)

          real(dblp), intent(in) :: lon
          real(dblp), intent(in) :: coslat
          real(dblp), intent(in) :: t
          real(dblp)             :: v

          real(dblp) :: lonp

          lonp = lon - 2.0_dblp*PI*t/TC_PERIOD
          v = (2.5_dblp/TC_PERIOD)*sin(lonp)*coslat**3*cos(PI*t/TC_PERIOD)

        end function v_div

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Initial fields
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine tc_cell_means(g, ic, phi)

          type(fvt_grid_t), intent(in)  :: g
          integer(ip),      intent(in)  :: ic
          real(dblp),       intent(out) :: phi(:,:)    !< (nlat,nlon)

          real(dblp)  :: xq(TC_NQ_CELL), wq(TC_NQ_CELL), mu, lon, s
          integer(ip) :: i, j, a, b

          call fvt_gauss_nodes(TC_NQ_CELL, xq, wq)
          do j = 1, g%nlon
            do i = 1, g%nlat
              s = 0.0_dblp
              do a = 1, TC_NQ_CELL
                mu = 0.5_dblp*(g%mu_face(i-1) + g%mu_face(i)) + 0.5_dblp*g%dmu(i)*xq(a)
                do b = 1, TC_NQ_CELL
                  lon = g%lon(j) + 0.5_dblp*g%dlon*xq(b)
                  s   = s + wq(a)*wq(b)*ic_point(ic, lon, asin(mu))
                enddo
              enddo
              phi(i, j) = 0.25_dblp*s
            enddo
          enddo

        end subroutine tc_cell_means

        function ic_point(ic, lon, lat) result(phi)

          integer(ip), intent(in) :: ic
          real(dblp),  intent(in) :: lon
          real(dblp),  intent(in) :: lat
          real(dblp)              :: phi

          real(dblp) :: r1, r2, d1, d2, x, y, z

          select case (ic)
            case (TC_IC_WG_BELL)
              r1 = gc_dist(lon, lat, 1.5_dblp*PI, 0.0_dblp)
              phi = 0.0_dblp
              if (r1 < 1.0_dblp/3.0_dblp) phi = 0.5_dblp*(1.0_dblp + cos(3.0_dblp*PI*r1))
            case (TC_IC_HILLS)
              x = cos(lat)*cos(lon)
              y = cos(lat)*sin(lon)
              z = sin(lat)
              phi = HILL_HMAX*exp(-HILL_B*((x - cos(LON1))**2 + (y - sin(LON1))**2 + z**2))                &
                  + HILL_HMAX*exp(-HILL_B*((x - cos(LON2))**2 + (y - sin(LON2))**2 + z**2))
            case (TC_IC_BELLS, TC_IC_CORRELATED)
              r1 = gc_dist(lon, lat, LON1, 0.0_dblp)
              r2 = gc_dist(lon, lat, LON2, 0.0_dblp)
              phi = BELL_BG
              if (r1 < BELL_R) phi = BELL_BG + BELL_AMP*0.5_dblp*(1.0_dblp + cos(PI*r1/BELL_R))
              if (r2 < BELL_R) phi = BELL_BG + BELL_AMP*0.5_dblp*(1.0_dblp + cos(PI*r2/BELL_R))
              if (ic == TC_IC_CORRELATED) phi = CORR_A*phi*phi + CORR_B
            case (TC_IC_CYLINDERS)
              r1 = gc_dist(lon, lat, LON1, 0.0_dblp)
              r2 = gc_dist(lon, lat, LON2, 0.0_dblp)
              d1 = abs(modulo(lon - LON1 + PI, 2.0_dblp*PI) - PI)
              d2 = abs(modulo(lon - LON2 + PI, 2.0_dblp*PI) - PI)
              phi = CYL_BG
              if (r1 <= BELL_R .and. d1 >= BELL_R/6.0_dblp) phi = CYL_C
              if (r2 <= BELL_R .and. d2 >= BELL_R/6.0_dblp) phi = CYL_C
              if (r1 <= BELL_R .and. d1 < BELL_R/6.0_dblp .and. lat < -5.0_dblp/12.0_dblp*BELL_R) phi = CYL_C
              if (r2 <= BELL_R .and. d2 < BELL_R/6.0_dblp .and. lat >  5.0_dblp/12.0_dblp*BELL_R) phi = CYL_C
            case default
              error stop 'tc_cell_means: unknown initial field'
          end select

        end function ic_point

        pure function gc_dist(lon, lat, lonc, latc) result(r)

          real(dblp), intent(in) :: lon
          real(dblp), intent(in) :: lat
          real(dblp), intent(in) :: lonc
          real(dblp), intent(in) :: latc
          real(dblp)             :: r

          r = acos(max(-1.0_dblp, min(1.0_dblp, sin(latc)*sin(lat) + cos(latc)*cos(lat)*cos(lon - lonc))))

        end function gc_dist

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Diagnostics
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        pure function tc_total(g, phi) result(tot)

          type(fvt_grid_t), intent(in) :: g
          real(dblp),       intent(in) :: phi(:,:)
          real(dblp)                   :: tot

          integer(ip) :: i

          tot = 0.0_dblp
          do i = 1, g%nlat
            tot = tot + g%area(i)*sum(phi(i, :))
          enddo

        end function tc_total

        function tc_norms(g, phi, phit, dphi0) result(nrm)

          type(fvt_grid_t), intent(in) :: g
          real(dblp),       intent(in) :: phi(:,:)
          real(dblp),       intent(in) :: phit(:,:)
          real(dblp),       intent(in) :: dphi0      !< max - min of the initial field
          type(tc_norms_t)             :: nrm

          nrm%l1     = tc_total(g, abs(phi - phit))/tc_total(g, abs(phit))
          nrm%l2     = sqrt(tc_total(g, (phi - phit)**2)/tc_total(g, phit**2))
          nrm%linf   = maxval(abs(phi - phit))/maxval(abs(phit))
          nrm%phimin = (minval(phi) - minval(phit))/dphi0
          nrm%phimax = (maxval(phi) - maxval(phit))/dphi0

        end function tc_norms

        subroutine tc_filament(g, phi, phi0, tau, lf)

          type(fvt_grid_t), intent(in)  :: g
          real(dblp),       intent(in)  :: phi(:,:)
          real(dblp),       intent(in)  :: phi0(:,:)
          real(dblp),       intent(in)  :: tau(:)
          real(dblp),       intent(out) :: lf(:)

          real(dblp)  :: a0, a1, w(g%nlat, g%nlon)
          integer(ip) :: k

          ! the tolerance keeps background cells (exactly tau = 0.1 up to round-off in the quadrature) on one side
          w = spread(g%area, 2, g%nlon)
          do k = 1, size(tau)
            a0 = sum(w, mask=phi0 >= tau(k) - TC_TAU_TOL)
            a1 = sum(w, mask=phi >= tau(k) - TC_TAU_TOL)
            lf(k) = 0.0_dblp
            if (a0 > 0.0_dblp) lf(k) = 100.0_dblp*a1/a0
          enddo

        end subroutine tc_filament

        subroutine tc_mixing(g, chi, xi, lr, lu, lo)

          type(fvt_grid_t), intent(in)  :: g
          real(dblp),       intent(in)  :: chi(:,:)   !< cosine bells
          real(dblp),       intent(in)  :: xi(:,:)    !< correlated cosine bells
          real(dblp),       intent(out) :: lr         !< real mixing
          real(dblp),       intent(out) :: lu         !< range-preserving unmixing
          real(dblp),       intent(out) :: lo         !< overshooting

          real(dblp)  :: d, atot, c, x, fline
          integer(ip) :: i, j

          lr = 0.0_dblp
          lu = 0.0_dblp
          lo = 0.0_dblp
          atot = 0.0_dblp
          do j = 1, g%nlon
            do i = 1, g%nlat
              c = chi(i, j)
              x = xi(i, j)
              d = mix_distance(c, x)*g%area(i)
              atot = atot + g%area(i)
              fline = MIX_XI_MAX + (MIX_XI_MIN - MIX_XI_MAX)*(c - MIX_CHI_MIN)/(MIX_CHI_MAX - MIX_CHI_MIN)
              if (c >= MIX_CHI_MIN .and. c <= MIX_CHI_MAX .and. x >= fline .and. x <= CORR_A*c*c + CORR_B) then
                lr = lr + d
              else if (c >= MIX_CHI_MIN .and. c <= MIX_CHI_MAX .and. x >= MIX_XI_MIN .and. x <= MIX_XI_MAX) then
                lu = lu + d
              else
                lo = lo + d
              endif
            enddo
          enddo
          lr = lr/atot
          lu = lu/atot
          lo = lo/atot

        end subroutine tc_mixing

        ! Normalised distance from (c, x) to the closest point of the curve x = a c**2 + b on the initial range of c:
        ! coarse sampling, then golden-section refinement around the best sample
        pure function mix_distance(c, x) result(d)

          real(dblp), intent(in) :: c
          real(dblp), intent(in) :: x
          real(dblp)             :: d

          integer(ip), parameter :: NSAMP = 200, NGOLD = 80
          real(dblp),  parameter :: GR = 0.618033988749894848_dblp
          real(dblp)  :: h, cb, lb, ca, cc, c1, c2, f1, f2, cs
          integer(ip) :: k

          h  = (MIX_CHI_MAX - MIX_CHI_MIN)/real(NSAMP, dblp)
          cb = MIX_CHI_MIN
          lb = huge(1.0_dblp)
          do k = 0, NSAMP
            cs = MIX_CHI_MIN + real(k, dblp)*h
            if (mix_l2(cs, c, x) < lb) then
              lb = mix_l2(cs, c, x)
              cb = cs
            endif
          enddo
          ca = max(MIX_CHI_MIN, cb - h)
          cc = min(MIX_CHI_MAX, cb + h)
          c1 = cc - GR*(cc - ca)
          c2 = ca + GR*(cc - ca)
          f1 = mix_l2(c1, c, x)
          f2 = mix_l2(c2, c, x)
          do k = 1, NGOLD
            if (f1 < f2) then
              cc = c2
              c2 = c1
              f2 = f1
              c1 = cc - GR*(cc - ca)
              f1 = mix_l2(c1, c, x)
            else
              ca = c1
              c1 = c2
              f1 = f2
              c2 = ca + GR*(cc - ca)
              f2 = mix_l2(c2, c, x)
            endif
          enddo
          d = sqrt(min(lb, f1, f2))

        end function mix_distance

        pure function mix_l2(cs, c, x) result(l2)

          real(dblp), intent(in) :: cs
          real(dblp), intent(in) :: c
          real(dblp), intent(in) :: x
          real(dblp)             :: l2

          l2 = ((c - cs)/(MIX_CHI_MAX - MIX_CHI_MIN))**2 + ((x - (CORR_A*cs*cs + CORR_B))/(MIX_XI_MAX - MIX_XI_MIN))**2

        end function mix_l2

      end module fvt_testcases_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
