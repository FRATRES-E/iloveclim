!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright for initial code: ECBilt authors (atmdyn_mod, atmphys_mod, atmmois_mod), algorithms reproduced here
!      (ec_rggtosp, ec_sptogg, ec_ddl, ec_lap, ec_trafluxdiv, ec_hdiff, ec_moisbalance, transport part of ec_moisture)

!   Further development:
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
!      MODULE: [fvt_specref_mod]
!
!>     @brief Standalone copy of ECBilt's current moisture transport (spectral baseline for the FV transport suite).
!
!      DESCRIPTION:
!>     Reproduces, for one tracer and transport only (no evaporation, precipitation or vertical exchange), the chain
!!     of ec_moisture with ihavm = 1 and umoisr = 1 (ROADMAP §1.5):
!!       1. grid -> spectral -> grid round trip (ec_rggtosp, ec_sptogg): truncation to T(ntrunc);
!!       2. negative values set to zero, the added amount accumulated in cormois;
!!       3. ec_trafluxdiv: u dq/dlon / cos(lat) + v cos(lat) dq/dmu (derivatives of the spectral field of step 1,
!!          before clipping) + q div (grid field after clipping), then truncated and global mean removed (vv(1) = 0);
!!       4. ec_hdiff: spectral Laplacian -n(n+1) q / tdifq (ec_lap has no 1/radius**2: mode n is damped in
!!          tdifq / (n(n+1)) days), global mean removed;
!!       5. forward Euler: q = q + dt (-advection + diffusion), negative values set to zero into cormois;
!!       6. ec_moisbalance: the cormois additions removed again latitude band by latitude band, from each pole to the
!!          equator, passing any deficit to the next band (weights pw(:,1)).
!!     Units: unit sphere, time in days; the winds u, v (rad/day) and the divergence (1/day) are point values at the
!!     Gaussian grid points, as u800 + udivg, v800 + vdivg, divg in the model. The Legendre tables are the coef.dat
!!     equivalent of fvt_legendre_tables (any truncation), scaled as in ec_initdyn; the FFT is ECBilt's NAG copy
!!     (src/ecbilt/nag.f).
!!
!!     Contract (public entry points):
!!       spr_init(sr,g)                     : tables for the grid g, truncation (nlon - 1)/3 (T21 on 64 longitudes).
!!       spr_step(sr,q,u,v,div,dt,tdifq[,clipped]) : one ec_moisture transport step; tdifq <= 0 switches the
!!                                            diffusion off; optional amount clipped into cormois.
!!       spr_trafluxdiv, spr_hdiff, spr_rggtosp, spr_sptogg : the individual operators (exposed for the tests).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module fvt_specref_mod

        use global_constants_mod, only: ip, dblp=>dp
        use fvt_legendre_mod,     only: fvt_legendre_tables, fvt_spec_index, fvt_nsh
        use fvt_grid_mod,         only: fvt_grid_t

        implicit none

        private

        public :: spr_t, spr_init, spr_step, spr_trafluxdiv, spr_hdiff, spr_rggtosp, spr_sptogg

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Types
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        type :: spr_t
          integer(ip)             :: nlat = 0, nlon = 0, nt = -1, nsh = 0
          integer(ip), allocatable :: nshm(:)       !< (0:nt)
          real(dblp),  allocatable :: pp(:,:)       !< (nlat,nsh) P * sqrt(nlon)
          real(dblp),  allocatable :: pd(:,:)       !< (nlat,nsh) dP/dmu * sqrt(nlon)
          real(dblp),  allocatable :: pw(:,:)       !< (nlat,nsh) (w/2) P / sqrt(nlon)
          real(dblp),  allocatable :: rm(:)         !< (nsh) zonal wavenumber m
          real(dblp),  allocatable :: lap(:)        !< (nsh) -n(n+1) (ec_lap, rinhel(:,0))
          real(dblp),  allocatable :: cosfi(:)      !< (nlat) cos(latitude) of the Gaussian points
          real(dblp),  allocatable :: trigd(:), trigi(:)
        end type spr_t

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Initialisation (ec_initdyn equivalent)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine spr_init(sr, g)

          type(spr_t),      intent(inout) :: sr
          type(fvt_grid_t), intent(in)    :: g

          integer(ip), allocatable :: ll(:)
          integer(ip) :: k, m, n
          real(dblp)  :: sq
          real(dblp), allocatable :: dummy(:,:), wgg(:)
          integer     :: ifail

          sr%nlat = g%nlat
          sr%nlon = g%nlon
          sr%nt   = (g%nlon - 1)/3
          sr%nsh  = fvt_nsh(sr%nt)
          allocate(sr%nshm(0:sr%nt), ll(sr%nsh), sr%rm(sr%nsh), sr%lap(sr%nsh))
          call fvt_spec_index(sr%nt, sr%nshm, ll)
          k = 0
          do m = 0, sr%nt
            do n = m, sr%nt
              k = k + 1
              sr%rm(k)  = real(m, dblp)
              sr%lap(k) = -real(n*(n + 1), dblp)
            enddo
          enddo
          allocate(sr%pp(sr%nlat, sr%nsh), sr%pd(sr%nlat, sr%nsh), sr%pw(sr%nlat, sr%nsh))
          call fvt_legendre_tables(sr%nt, sr%nlat, sr%pp, sr%pd, sr%pw)
          sq = sqrt(real(sr%nlon, dblp))
          sr%pp = sr%pp*sq
          sr%pd = sr%pd*sq
          sr%pw = sr%pw/sq
          allocate(sr%cosfi(sr%nlat))
          sr%cosfi = cos(g%lat)

          allocate(sr%trigd(2*sr%nlon), sr%trigi(2*sr%nlon), dummy(sr%nlat, sr%nlon), wgg(sr%nlat*sr%nlon))
          dummy = 1.0_dblp
          ifail = 0
          call c06fpf(int(sr%nlat), int(sr%nlon), dummy, 'i', sr%trigd, wgg, ifail)
          dummy = 1.0_dblp
          ifail = 0
          call c06fqf(int(sr%nlat), int(sr%nlon), dummy, 'i', sr%trigi, wgg, ifail)

        end subroutine spr_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   One transport step of ec_moisture
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine spr_step(sr, q, u, v, div, dt, tdifq, clipped)

          type(spr_t), intent(in)    :: sr
          real(dblp),  intent(inout) :: q(:,:)      !< (nlat,nlon) grid field
          real(dblp),  intent(in)    :: u(:,:)      !< zonal wind at the grid points (rad/day)
          real(dblp),  intent(in)    :: v(:,:)      !< meridional wind at the grid points (rad/day)
          real(dblp),  intent(in)    :: div(:,:)    !< divergence at the grid points (1/day)
          real(dblp),  intent(in)    :: dt          !< step (days)
          real(dblp),  intent(in)    :: tdifq       !< diffusion time scale (days); <= 0: no diffusion
          real(dblp), optional, intent(out) :: clipped  !< Gaussian-weighted total of cormois (amount added by clipping
                                                       !! in steps 2 and 5, removed again by ec_moisbalance)

          real(dblp)  :: qs(sr%nsh, 2), cormois(sr%nlat, sr%nlon), hdivmg(sr%nlat, sr%nlon), hdmoisg(sr%nlat, sr%nlon)
          integer(ip) :: i, j

          ! 1. round trip (ec_rggtosp, ec_sptogg)
          call spr_rggtosp(sr, q, qs)
          call spr_sptogg(sr, qs, q, sr%pp)

          ! 2. clipping into cormois
          cormois = 0.0_dblp
          do j = 1, sr%nlon
            do i = 1, sr%nlat
              if (q(i, j) < 0.0_dblp) then
                cormois(i, j) = cormois(i, j) - q(i, j)
                q(i, j) = 0.0_dblp
              endif
            enddo
          enddo

          ! 3. advection (ec_trafluxdiv), 4. diffusion (ec_hdiff)
          call spr_trafluxdiv(sr, qs, q, u, v, div, hdivmg)
          call spr_hdiff(sr, qs, tdifq, hdmoisg)

          ! 5. forward Euler, clipping into cormois
          do j = 1, sr%nlon
            do i = 1, sr%nlat
              q(i, j) = q(i, j) + dt*(-hdivmg(i, j) + hdmoisg(i, j))
              if (q(i, j) < 0.0_dblp) then
                cormois(i, j) = cormois(i, j) - q(i, j)
                q(i, j) = 0.0_dblp
              endif
            enddo
          enddo

          if (present(clipped)) then
            clipped = 0.0_dblp
            do i = 1, sr%nlat
              clipped = clipped + 2.0_dblp*sqrt(real(sr%nlon, dblp))*sr%pw(i, 1)*sum(cormois(i, :))
            enddo
          endif

          ! 6. hole filling (ec_moisbalance)
          call spr_moisbalance(sr, q, cormois)

        end subroutine spr_step

        subroutine spr_trafluxdiv(sr, ctrasp, ctra, u, v, div, tfdiv)

          type(spr_t), intent(in)  :: sr
          real(dblp),  intent(in)  :: ctrasp(:,:)
          real(dblp),  intent(in)  :: ctra(:,:)
          real(dblp),  intent(in)  :: u(:,:)
          real(dblp),  intent(in)  :: v(:,:)
          real(dblp),  intent(in)  :: div(:,:)
          real(dblp),  intent(out) :: tfdiv(:,:)

          real(dblp)  :: vv(sr%nsh, 2), dcdl(sr%nlat, sr%nlon), dcdm(sr%nlat, sr%nlon)
          integer(ip) :: i, j

          call spr_ddl(sr, ctrasp, vv)
          call spr_sptogg(sr, vv, dcdl, sr%pp)
          call spr_sptogg(sr, ctrasp, dcdm, sr%pd)
          do j = 1, sr%nlon
            do i = 1, sr%nlat
              tfdiv(i, j) = dcdl(i, j)*u(i, j)/sr%cosfi(i) + dcdm(i, j)*v(i, j)*sr%cosfi(i) + ctra(i, j)*div(i, j)
            enddo
          enddo
          call spr_rggtosp(sr, tfdiv, vv)
          vv(1, 1) = 0.0_dblp
          call spr_sptogg(sr, vv, tfdiv, sr%pp)

        end subroutine spr_trafluxdiv

        subroutine spr_hdiff(sr, qs, tdifq, hdmg)

          type(spr_t), intent(in)  :: sr
          real(dblp),  intent(in)  :: qs(:,:)
          real(dblp),  intent(in)  :: tdifq
          real(dblp),  intent(out) :: hdmg(:,:)

          real(dblp) :: hd(sr%nsh, 2), difq

          difq = 0.0_dblp
          if (tdifq > 0.0_dblp) difq = 1.0_dblp/tdifq
          hd(:, 1) = difq*sr%lap*qs(:, 1)
          hd(:, 2) = difq*sr%lap*qs(:, 2)
          hd(1, 1) = 0.0_dblp
          call spr_sptogg(sr, hd, hdmg, sr%pp)

        end subroutine spr_hdiff

        ! ec_moisbalance, same order of operations
        subroutine spr_moisbalance(sr, q, cormois)

          type(spr_t), intent(in)    :: sr
          real(dblp),  intent(inout) :: q(:,:)
          real(dblp),  intent(inout) :: cormois(:,:)

          integer(ip) :: i, j, nlat, nlon
          real(dblp)  :: gmc, gmm, gfac

          nlat = sr%nlat
          nlon = sr%nlon
          do i = 1, nlat/2
            gmc = sum(cormois(i, :))
            gmm = sum(q(i, :))
            if (gmm > 0.0_dblp .and. gmc > 0.0_dblp) then
              if (gmm > gmc) then
                gfac = gmc/gmm
                do j = 1, nlon
                  q(i, j) = q(i, j) - gfac*q(i, j)
                enddo
              else
                gfac = gmm/gmc
                gmc  = (gmc - gmm)*sr%pw(i, 1)/(sr%pw(i+1, 1)*real(nlon, dblp))
                do j = 1, nlon
                  q(i, j) = 0.0_dblp
                  cormois(i, j) = gfac*cormois(i, j)
                  cormois(i+1, j) = cormois(i+1, j) + gmc
                enddo
              endif
            endif
          enddo
          do i = nlat, 1 + nlat/2, -1
            gmc = sum(cormois(i, :))
            gmm = sum(q(i, :))
            if (gmm > 0.0_dblp .and. gmc > 0.0_dblp) then
              if (gmm > gmc) then
                gfac = 1.0_dblp - gmc/gmm
                do j = 1, nlon
                  q(i, j) = q(i, j)*gfac
                enddo
              else
                gfac = gmm/gmc
                gmc  = (gmc - gmm)*sr%pw(i, 1)/(sr%pw(i-1, 1)*real(nlon, dblp))
                do j = 1, nlon
                  q(i, j) = 0.0_dblp
                  cormois(i, j) = gfac*cormois(i, j)
                  cormois(i-1, j) = cormois(i-1, j) + gmc
                enddo
              endif
            endif
          enddo

        end subroutine spr_moisbalance

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Transforms (ec_ddl, ec_sptogg, ec_rggtosp)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine spr_ddl(sr, as, dadl)

          type(spr_t), intent(in)  :: sr
          real(dblp),  intent(in)  :: as(:,:)
          real(dblp),  intent(out) :: dadl(:,:)

          dadl(:, 1) = -sr%rm*as(:, 2)
          dadl(:, 2) =  sr%rm*as(:, 1)

        end subroutine spr_ddl

        subroutine spr_sptogg(sr, as, agg, pploc)

          type(spr_t), intent(in)  :: sr
          real(dblp),  intent(in)  :: as(:,:)
          real(dblp),  intent(out) :: agg(:,:)
          real(dblp),  intent(in)  :: pploc(:,:)

          real(dblp)  :: trig(2*sr%nlon), wgg(sr%nlat*sr%nlon)
          integer(ip) :: k, k1, k2, m, mr, mi
          integer     :: ifail

          agg = 0.0_dblp
          k2 = sr%nshm(0)
          do k = 1, k2
            agg(:, 1) = agg(:, 1) + as(k, 1)*pploc(:, k)
          enddo
          do m = 1, sr%nt
            mr = m + 1
            mi = sr%nlon + 1 - m
            k1 = k2 + 1
            k2 = k2 + sr%nshm(m)
            do k = k1, k2
              agg(:, mr) = agg(:, mr) + as(k, 1)*pploc(:, k)
              agg(:, mi) = agg(:, mi) - as(k, 2)*pploc(:, k)
            enddo
          enddo
          trig = sr%trigi
          ifail = 0
          call c06fqf(int(sr%nlat), int(sr%nlon), agg, 'r', trig, wgg, ifail)

        end subroutine spr_sptogg

        subroutine spr_rggtosp(sr, agg, as)

          type(spr_t), intent(in)  :: sr
          real(dblp),  intent(in)  :: agg(:,:)
          real(dblp),  intent(out) :: as(:,:)

          real(dblp)  :: store(sr%nlat, sr%nlon), trig(2*sr%nlon), wgg(sr%nlat*sr%nlon)
          integer(ip) :: i, k, k1, k2, m, mr, mi
          integer     :: ifail

          store = agg
          trig  = sr%trigd
          ifail = 0
          call c06fpf(int(sr%nlat), int(sr%nlon), store, 'r', trig, wgg, ifail)
          as = 0.0_dblp
          k2 = sr%nshm(0)
          do k = 1, k2
            do i = 1, sr%nlat
              as(k, 1) = as(k, 1) + store(i, 1)*sr%pw(i, k)
            enddo
          enddo
          do m = 1, sr%nt
            mr = m + 1
            mi = sr%nlon + 1 - m
            k1 = k2 + 1
            k2 = k2 + sr%nshm(m)
            do k = k1, k2
              do i = 1, sr%nlat
                as(k, 1) = as(k, 1) + store(i, mr)*sr%pw(i, k)
                as(k, 2) = as(k, 2) + store(i, mi)*sr%pw(i, k)
              enddo
            enddo
          enddo

        end subroutine spr_rggtosp

      end module fvt_specref_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
