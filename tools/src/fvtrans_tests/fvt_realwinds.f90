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
!      PROGRAM: [fvt_realwinds]
!
!>     @brief FV transport with the real ECBilt 800 hPa winds (ROADMAP step 1.1, patch wiso:1d:0003).
!
!      DESCRIPTION:
!>     Usage: fvt_realwinds <dump file> <outdir>
!!     Reads a dump of the model state taken just before ec_moisture (temporary diagnostic, not in the model):
!!       header : nlat, nlon, nsh2 (int32), dtime, om, radius, umoisr (float64);
!!       per step: istep (int32), psi(:,3), chi(:,3) (nsh2 each), u800 + udivg(:,:,3), v800 + vdivg(:,:,3),
!!                 divg(:,:,3), rmoisg(:,:,iwater) (nlat x nlon each), all float64, big-endian, stream access.
!!     For every step:
!!       1. winds consistency: u, v and the divergence synthesised from psi and chi with the fvt_specflux_mod
!!          convention against the model's own grid fields (relative maximum differences);
!!       2. face fluxes from fvt_specflux_mod (scale_psi = umoisr dt om, scale_chi = umoisr dt / radius**2): zonal and
!!          meridional Courant numbers and outflow fractions per cell, number of FV substeps for out_max = 1 and 0.75;
!!       3. one transport step of the dumped moisture: FV (positive-definite PPM, density mode) and the current
!!          spectral chain (fvt_specref_mod, winds and divergence times umoisr, tdifq = 1000 d and without
!!          diffusion): amount clipped into cormois, increments, minimum after FV.
!!     Outputs: <outdir>/realwinds_steps.csv (one line per step) and <outdir>/realwinds_rows.csv (statistics per
!!     month (30 days) and latitude row).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      program fvt_realwinds

        use, intrinsic :: iso_fortran_env, only: int32, real64
        use global_constants_mod, only: ip, dblp=>dp, str_len
        use fvt_legendre_mod,     only: fvt_legendre_eval, fvt_nsh
        use fvt_grid_mod,         only: fvt_grid_t, fvt_grid_init
        use fvt_specflux_mod,     only: fvt_specflux_t, fvt_specflux_init, fvt_specflux_fluxes
        use fvt_ffsl_mod,         only: fvt_ffsl_density, fvt_ffsl_nsub, FVT_RECON_PPM, FVT_LIM_POSDEF
        use fvt_specref_mod,      only: spr_t, spr_init, spr_step

        implicit none

        integer(ip), parameter :: STEPS_PER_MONTH = 180     !< 30 days of 4 h steps
        real(dblp),  parameter :: TDIFQ_S = 1000.0_dblp*86400.0_dblp

        character(len=str_len) :: dumpfile, outdir
        integer(int32)         :: h_nlat, h_nlon, h_nsh2, istep
        real(real64)           :: dtime, om, radius, umoisr
        integer(ip)            :: nlat, nlon, nt, nsh, u_in, u_st, u_rw, ios, n, i, j, jw, k, m, mon, nmon
        type(fvt_grid_t)       :: g
        type(fvt_specflux_t)   :: sf
        type(spr_t)            :: sr
        real(dblp), allocatable :: psi(:), chi(:), ug(:,:), vg(:,:), dg(:,:), qg(:,:)
        real(dblp), allocatable :: p(:,:), dp(:,:), cosm(:,:), sinm(:,:), lapn(:)
        real(dblp), allocatable :: us(:,:), vs(:,:), ds(:,:), fx(:,:), fy(:,:)
        real(dblp), allocatable :: qf(:,:,:), qs(:,:), qs0(:,:), w2(:,:)
        real(dblp), allocatable :: rmax_cx(:,:), rsum_cx(:,:), rmax_cy(:,:), rsum_cy(:,:), rmax_out(:,:)
        integer(ip), allocatable :: rcount(:)
        real(dblp)  :: cx, cy, out, cxmax, cymax, outmax, eu, ev, ed, msum, clipped, clipped0, cosl
        real(dblp)  :: rmsf, rmss, rmss0, rmsd, rmsd0, wsum, dmass_f, dmass_s
        integer(ip) :: ns1, ns075, nsf

        if (command_argument_count() /= 2) then
          write(*, '(a)') 'usage: fvt_realwinds <dump file> <outdir>'
          error stop 2
        endif
        call get_command_argument(1, dumpfile)
        call get_command_argument(2, outdir)

        open(newunit=u_in, file=trim(dumpfile), access='stream', form='unformatted', status='old', action='read', &
             iostat=ios)
        if (ios /= 0) error stop 'fvt_realwinds: cannot open the dump file'
        read(u_in) h_nlat, h_nlon, h_nsh2, dtime, om, radius, umoisr
        nlat = int(h_nlat, ip)
        nlon = int(h_nlon, ip)
        nt   = (nlon - 1)/3
        nsh  = fvt_nsh(nt)
        if (int(h_nsh2, ip) /= 2*nsh) error stop 'fvt_realwinds: nsh2 does not match the grid'
        write(*, '(a,i0,a,i0,a,i0,a,f8.1,a,f5.2)') 'fvt_realwinds: grid ', nlat, ' x ', nlon, ', T', nt, ', dt ', dtime, &
                                                   ' s, umoisr ', umoisr

        call fvt_grid_init(g, nlat, nlon)
        call fvt_specflux_init(sf, g, nt)
        call spr_init(sr, g)
        allocate(psi(2*nsh), chi(2*nsh), ug(nlat, nlon), vg(nlat, nlon), dg(nlat, nlon), qg(nlat, nlon))
        allocate(p(nlat, nsh), dp(nlat, nsh), cosm(0:nt, nlon), sinm(0:nt, nlon), lapn(nsh))
        allocate(us(nlat, nlon), vs(nlat, nlon), ds(nlat, nlon), fx(nlat, nlon), fy(0:nlat, nlon))
        allocate(qf(nlat, nlon, 1), qs(nlat, nlon), qs0(nlat, nlon), w2(nlat, nlon))
        call fvt_legendre_eval(nt, g%mu, p, dp)
        do j = 1, nlon
          do m = 0, nt
            cosm(m, j) = cos(real(m, dblp)*g%lon(j))
            sinm(m, j) = sin(real(m, dblp)*g%lon(j))
          enddo
        enddo
        k = 0
        do m = 0, nt
          do n = m, nt
            k = k + 1
            lapn(k) = -real(n*(n + 1), dblp)
          enddo
        enddo
        do i = 1, nlat
          w2(i, :) = g%wgt(i)
        enddo
        nmon = 12
        allocate(rmax_cx(nmon, nlat), rsum_cx(nmon, nlat), rmax_cy(nmon, nlat), rsum_cy(nmon, nlat), &
                 rmax_out(nmon, nlat), rcount(nmon))
        rmax_cx = 0.0_dblp
        rsum_cx = 0.0_dblp
        rmax_cy = 0.0_dblp
        rsum_cy = 0.0_dblp
        rmax_out = 0.0_dblp
        rcount = 0

        open(newunit=u_st, file=trim(outdir)//'/realwinds_steps.csv', status='replace', action='write')
        write(u_st, '(a)') 'istep,month,err_u,err_v,err_div,cx_max,cy_max,out_max_nonshift,nsub_1,nsub_075,'// &
                           'clipped_frac,clipped_frac_nodif,rms_incr_fv,rms_incr_spec,rms_incr_spec_nodif,'// &
                           'rms_diff,rms_diff_nodif,qmin_fv,mass_change_fv,mass_change_spec,nsub_fv'

        do
          read(u_in, iostat=ios) istep, psi, chi, ug, vg, dg, qg
          if (ios /= 0) exit
          mon = min(nmon, (int(istep, ip) - 1)/STEPS_PER_MONTH + 1)

          ! 1. winds from psi and chi against the model fields
          call synth_winds(psi, chi, us, vs, ds)
          eu = maxval(abs(us - ug))/maxval(abs(ug))
          ev = maxval(abs(vs - vg))/maxval(abs(vg))
          ed = maxval(abs(ds - dg))/maxval(abs(dg))

          ! 2. face fluxes, Courant numbers, substeps
          call fvt_specflux_fluxes(sf, g, reshape(psi, [nsh, 2]), reshape(chi, [nsh, 2]), umoisr*dtime*om,           &
                                   umoisr*dtime/radius**2, fx, fy)
          cxmax = 0.0_dblp
          cymax = 0.0_dblp
          outmax = 0.0_dblp
          do j = 1, nlon
            jw = modulo(j - 2, nlon) + 1
            do i = 1, nlat
              cx  = max(abs(fx(i, j)), abs(fx(i, jw)))/g%area(i)
              cy  = max(abs(fy(i, j)), abs(fy(i-1, j)))/g%area(i)
              out = (max(fx(i, j), 0.0_dblp) + max(-fx(i, jw), 0.0_dblp) + max(fy(i, j), 0.0_dblp)                   &
                     + max(-fy(i-1, j), 0.0_dblp))/g%area(i)
              cxmax = max(cxmax, cx)
              cymax = max(cymax, cy)
              if (cx < 1.0_dblp) outmax = max(outmax, out)
              rmax_cx(mon, i) = max(rmax_cx(mon, i), cx)
              rsum_cx(mon, i) = rsum_cx(mon, i) + cx
              rmax_cy(mon, i) = max(rmax_cy(mon, i), cy)
              rsum_cy(mon, i) = rsum_cy(mon, i) + cy
              if (cx < 1.0_dblp) rmax_out(mon, i) = max(rmax_out(mon, i), out)
            enddo
          enddo
          rcount(mon) = rcount(mon) + 1
          ns1   = fvt_ffsl_nsub(g, fx, fy)
          ns075 = fvt_ffsl_nsub(g, fx, fy, 0.75_dblp)

          ! 3. one transport step of the dumped moisture
          msum = sum(w2*qg)
          qf(:, :, 1) = qg
          call fvt_ffsl_density(g, fx, fy, qf, FVT_RECON_PPM, FVT_LIM_POSDEF, nsf)
          qs  = qg
          call spr_step(sr, qs, umoisr*ug/radius, umoisr*vg/radius, umoisr*dg, dtime, TDIFQ_S, clipped)
          qs0 = qg
          call spr_step(sr, qs0, umoisr*ug/radius, umoisr*vg/radius, umoisr*dg, dtime, 0.0_dblp, clipped0)
          wsum  = sum(w2)
          rmsf  = sqrt(sum(w2*(qf(:, :, 1) - qg)**2)/wsum)
          rmss  = sqrt(sum(w2*(qs - qg)**2)/wsum)
          rmss0 = sqrt(sum(w2*(qs0 - qg)**2)/wsum)
          rmsd  = sqrt(sum(w2*(qf(:, :, 1) - qs)**2)/wsum)
          rmsd0 = sqrt(sum(w2*(qf(:, :, 1) - qs0)**2)/wsum)
          dmass_f = (sum(spread(g%area, 2, nlon)*qf(:, :, 1)) - sum(spread(g%area, 2, nlon)*qg))                    &
                    /sum(spread(g%area, 2, nlon)*qg)
          dmass_s = (sum(w2*qs) - msum)/msum

          write(u_st, '(i0,a,i0,18(a,es14.6),a,i0)') istep, ',', mon, ',', eu, ',', ev, ',', ed, ',', cxmax, ',', cymax, &
                ',', outmax, ',', real(ns1, dblp), ',', real(ns075, dblp), ',', clipped/msum, ',', clipped0/msum,         &
                ',', rmsf, ',', rmss, ',', rmss0, ',', rmsd, ',', rmsd0, ',', minval(qf), ',', dmass_f, ',', dmass_s,       &
                ',', nsf
          if (mod(istep, 360_int32) == 0) write(*, '(a,i0)') '  step ', istep
        enddo
        close(u_in)
        close(u_st)

        open(newunit=u_rw, file=trim(outdir)//'/realwinds_rows.csv', status='replace', action='write')
        write(u_rw, '(a)') 'month,row,lat,cx_max,cx_mean,cy_max,cy_mean,out_max_nonshift'
        do mon = 1, nmon
          if (rcount(mon) == 0) cycle
          do i = 1, nlat
            cosl = real(rcount(mon)*nlon, dblp)
            write(u_rw, '(i0,a,i0,6(a,es14.6))') mon, ',', i, ',', g%lat(i)*180.0_dblp/acos(-1.0_dblp), ',',          &
                  rmax_cx(mon, i), ',', rsum_cx(mon, i)/cosl, ',', rmax_cy(mon, i), ',', rsum_cy(mon, i)/cosl, ',',   &
                  rmax_out(mon, i)
          enddo
        enddo
        close(u_rw)
        write(*, '(a)') 'fvt_realwinds: done'

      contains

        ! Winds and divergence at the Gaussian points from the spectral psi (non-dimensional) and chi (m**2/s), with the
        ! fvt_specflux_mod convention: field = sum P (c1 cos(m lon) - c2 sin(m lon)) (factor 2 for m > 0)
        subroutine synth_winds(ps, cs, u, v, d)

          real(dblp), intent(in)  :: ps(:)
          real(dblp), intent(in)  :: cs(:)
          real(dblp), intent(out) :: u(:,:)
          real(dblp), intent(out) :: v(:,:)
          real(dblp), intent(out) :: d(:,:)

          integer(ip) :: ii, jj, kk, mm, nn
          real(dblp)  :: f, c1, c2, e1, e2, cl, psl, psm, chl, chm, lap

          do jj = 1, nlon
            do ii = 1, nlat
              psl = 0.0_dblp
              psm = 0.0_dblp
              chl = 0.0_dblp
              chm = 0.0_dblp
              lap = 0.0_dblp
              kk = 0
              do mm = 0, nt
                do nn = mm, nt
                  kk = kk + 1
                  c1 = ps(kk)
                  c2 = ps(kk + nsh)
                  e1 = cs(kk)
                  e2 = cs(kk + nsh)
                  if (mm == 0) then
                    psm = psm + c1*dp(ii, kk)
                    chm = chm + e1*dp(ii, kk)
                    lap = lap + lapn(kk)*e1*p(ii, kk)
                  else
                    f   = 2.0_dblp*real(mm, dblp)*p(ii, kk)
                    psl = psl + f*(-c1*sinm(mm, jj) - c2*cosm(mm, jj))
                    chl = chl + f*(-e1*sinm(mm, jj) - e2*cosm(mm, jj))
                    psm = psm + 2.0_dblp*dp(ii, kk)*(c1*cosm(mm, jj) - c2*sinm(mm, jj))
                    chm = chm + 2.0_dblp*dp(ii, kk)*(e1*cosm(mm, jj) - e2*sinm(mm, jj))
                    lap = lap + 2.0_dblp*lapn(kk)*p(ii, kk)*(e1*cosm(mm, jj) - e2*sinm(mm, jj))
                  endif
                enddo
              enddo
              cl = cos(g%lat(ii))
              u(ii, jj) = -radius*om*psm*cl + chl/(radius*cl)
              v(ii, jj) =  radius*om*psl/cl + chm*cl/radius
              d(ii, jj) = lap/radius**2
            enddo
          enddo

        end subroutine synth_winds

      end program fvt_realwinds

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
