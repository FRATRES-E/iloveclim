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

!   Style sheet: v1.0.0 (standalone program, by analogy, D25)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      PROGRAM: [fvt_ratiocheck]
!
!>     @brief Offline checks of the isotope transport (ROADMAP 2c:0001) on the real ECBilt winds of a wind dump.
!
!      DESCRIPTION:
!
!>     Usage: fvt_ratiocheck <wind_dump.bin>
!!
!!     The model moves water in density mode (D29). With isotopes (D30, D32) it moves water as the carrier of
!!     fvt_ffsl_ratio, the ratios iwat17, iwat18, iwat2h on the carrier's fluxes (monotone limiter), and iwat16 as a
!!     separate density. Built with the model's compile options, this program checks on one year of real winds:
!!       C1  water from fvt_ffsl_ratio is bitwise identical to water from fvt_ffsl_density (D4 at the root), both
!!           from the dumped state at every step, and on a free transport-only run of the whole year;
!!       C2  a uniform ratio stays uniform to round-off (T2 offline);
!!       C3  a structured ratio field (smooth and step) stays within its initial range at every step (T4 offline,
!!           including the polar rows, O23), and every new ratio stays within the range of the old ratios over the
!!           cells that can reach it in the step (local bound: rows i +- nsub, columns j +- (nsub + the largest zonal
!!           shift of the row)), both from the dumped state at every step and on the free run;
!!       C4  "isotope mass" sum(area q R) is conserved: to round-off for a uniform ratio; to 1e-6 over the year for
!!           the structured ratios, whose transport is not exactly conservative (clip of fvt_ffsl_mod, ~1e-7 a year);
!!       C5  iwat16 in density mode is bitwise identical to water (same call, same data).
!!     One line per check, exit code 1 on failure.
!
!      AUTHORS: dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      program fvt_ratiocheck

        use, intrinsic :: iso_fortran_env, only: int32, real64
        use global_constants_mod, only: ip, dblp=>dp, str_len
        use fvt_legendre_mod,     only: fvt_nsh
        use fvt_grid_mod,         only: fvt_grid_t, fvt_grid_init
        use fvt_specflux_mod,     only: fvt_specflux_t, fvt_specflux_init, fvt_specflux_fluxes
        use fvt_ffsl_mod,         only: fvt_ffsl_density, fvt_ffsl_ratio, FVT_RECON_PPM, FVT_LIM_POSDEF, FVT_LIM_MONO

        implicit none

        integer(ip), parameter :: NR = 3                                !< ratios: uniform, smooth, step

        type(fvt_grid_t)       :: g
        type(fvt_specflux_t)   :: sf
        character(len=str_len) :: dumpfile
        integer(int32)         :: h_nlat, h_nlon, h_nsh2, istep
        real(real64)           :: dtime, om, radius, umoisr
        integer(ip)            :: nlat, nlon, nt, nsh, u_in, ios, nstep, i, j, k, nsub_d, nsub_r, nfail
        real(dblp), allocatable :: psi(:), chi(:), ug(:,:), vg(:,:), dg(:,:), qg(:,:)
        real(dblp), allocatable :: fx(:,:), fy(:,:), area(:,:)
        real(dblp), allocatable :: qd(:,:,:), qr(:,:), rr(:,:,:)        ! per-step check from the dumped state
        real(dblp), allocatable :: fd(:,:,:), fr(:,:), frr(:,:,:), f16(:,:,:), r0(:,:,:)   ! free run
        real(dblp)             :: rmin0(NR), rmax0(NR), exc(NR), mass_r0(NR), dmass(NR), uni_dev
        real(dblp)             :: lexc_step(NR), lexc_free(NR)
        real(dblp), allocatable :: rold(:,:,:)
        integer(ip)            :: nsub_f
        integer(ip)            :: nmis_step, nmis_free, nmis_16, nsub_diff

        nfail = 0
        call get_command_argument(1, dumpfile)
        if (len_trim(dumpfile) == 0) error stop 'usage: fvt_ratiocheck <wind_dump.bin>'
        open(newunit=u_in, file=trim(dumpfile), access='stream', form='unformatted', status='old', action='read', &
             iostat=ios)
        if (ios /= 0) error stop 'fvt_ratiocheck: cannot open the dump'
        read(u_in) h_nlat, h_nlon, h_nsh2, dtime, om, radius, umoisr
        nlat = int(h_nlat, ip)
        nlon = int(h_nlon, ip)
        nt   = (nlon - 1)/3
        nsh  = fvt_nsh(nt)
        if (int(h_nsh2, ip) /= 2*nsh) error stop 'fvt_ratiocheck: nsh2 does not match the grid'

        call fvt_grid_init(g, nlat, nlon)
        call fvt_specflux_init(sf, g, nt)
        allocate(psi(2*nsh), chi(2*nsh), ug(nlat, nlon), vg(nlat, nlon), dg(nlat, nlon), qg(nlat, nlon))
        allocate(fx(nlat, nlon), fy(0:nlat, nlon), area(nlat, nlon))
        allocate(qd(nlat, nlon, 1), qr(nlat, nlon), rr(nlat, nlon, NR))
        allocate(fd(nlat, nlon, 1), fr(nlat, nlon), frr(nlat, nlon, NR), f16(nlat, nlon, 1), r0(nlat, nlon, NR))
        allocate(rold(nlat, nlon, NR))
        do i = 1, nlat
          area(i, :) = g%area(i)
        enddo

        ! dmr&clo   initial ratios: uniform; smooth (latitude and longitude); step (a sharp band, with polar caps)
        do j = 1, nlon
          do i = 1, nlat
            r0(i, j, 1) = 0.988_dblp
            r0(i, j, 2) = 1.0_dblp + 0.02_dblp*g%mu(i) + 0.01_dblp*cos(g%lon(j))
            r0(i, j, 3) = merge(0.95_dblp, 1.05_dblp, abs(g%mu(i)) > 0.8_dblp .or. (j > nlon/4 .and. j <= nlon/2))
          enddo
        enddo
        do k = 1, NR
          rmin0(k) = minval(r0(:, :, k))
          rmax0(k) = maxval(r0(:, :, k))
        enddo

        nstep = 0
        nmis_step = 0
        nmis_free = 0
        nmis_16 = 0
        nsub_diff = 0
        exc = 0.0_dblp
        lexc_step = 0.0_dblp
        lexc_free = 0.0_dblp
        dmass = 0.0_dblp
        uni_dev = 0.0_dblp

        do
          read(u_in, iostat=ios) istep, psi, chi, ug, vg, dg, qg
          if (ios /= 0) exit
          nstep = nstep + 1

          ! dmr&clo   fluxes exactly as atm_fvmois_step
          call fvt_specflux_fluxes(sf, g, reshape(psi, [nsh, 2]), reshape(chi, [nsh, 2]), umoisr*dtime*om,       &
                                   umoisr*dtime/radius**2, fx, fy)

          ! C1, from the dumped state
          qd(:, :, 1) = qg
          qr = qg
          rr = r0
          call fvt_ffsl_density(g, fx, fy, qd, FVT_RECON_PPM, FVT_LIM_POSDEF, nsub_d)
          call fvt_ffsl_ratio(g, fx, fy, qr, rr, FVT_RECON_PPM, FVT_LIM_POSDEF, nsub_r, limiter_ratio=FVT_LIM_MONO)
          nmis_step = nmis_step + count(qd(:, :, 1) /= qr)
          if (nsub_d /= nsub_r) nsub_diff = nsub_diff + 1
          do k = 2, NR
            lexc_step(k) = max(lexc_step(k), local_excess(r0(:, :, k), rr(:, :, k), nsub_r)/(rmax0(k) - rmin0(k)))
          enddo

          ! free transport-only run over the year
          if (nstep == 1) then
            fd(:, :, 1) = qg
            fr = qg
            f16(:, :, 1) = qg
            frr = r0
            do k = 1, NR
              mass_r0(k) = sum(area*fr*frr(:, :, k))
            enddo
          endif
          call fvt_ffsl_density(g, fx, fy, fd, FVT_RECON_PPM, FVT_LIM_POSDEF)
          rold = frr
          call fvt_ffsl_ratio(g, fx, fy, fr, frr, FVT_RECON_PPM, FVT_LIM_POSDEF, nsub_f, limiter_ratio=FVT_LIM_MONO)
          do k = 2, NR
            lexc_free(k) = max(lexc_free(k), local_excess(rold(:, :, k), frr(:, :, k), nsub_f)/(rmax0(k) - rmin0(k)))
          enddo
          call fvt_ffsl_density(g, fx, fy, f16, FVT_RECON_PPM, FVT_LIM_POSDEF)
          nmis_free = nmis_free + count(fd(:, :, 1) /= fr)
          nmis_16   = nmis_16 + count(f16(:, :, 1) /= fd(:, :, 1))
          uni_dev = max(uni_dev, maxval(abs(frr(:, :, 1) - r0(1, 1, 1))))
          do k = 1, NR
            exc(k) = max(exc(k), (rmin0(k) - minval(frr(:, :, k)))/(rmax0(k) - rmin0(k) + tiny(1.0_dblp)),             &
                                 (maxval(frr(:, :, k)) - rmax0(k))/(rmax0(k) - rmin0(k) + tiny(1.0_dblp)))
            dmass(k) = max(dmass(k), abs(sum(area*fr*frr(:, :, k)) - mass_r0(k))/mass_r0(k))
          enddo
        enddo
        close(u_in)

        write(*,'(a,i0,a,i0,a,i0)') 'fvt_ratiocheck: ', nstep, ' steps, grid ', nlat, ' x ', nlon
        call check('C1 water ratio mode == density mode, per step (cells differing)', nmis_step == 0 .and. nsub_diff == 0, &
                   real(nmis_step + nsub_diff, dblp))
        call check('C1 water ratio mode == density mode, free run (cells differing)', nmis_free == 0,                     &
                   real(nmis_free, dblp))
        call check('C2 uniform ratio stays uniform (max |R - R0|)', uni_dev <= 1.0e-13_dblp, uni_dev)
        call check('C3 smooth ratio within initial range (excursion / range)', exc(2) <= 1.0e-12_dblp, exc(2))
        call check('C3 step ratio within initial range (excursion / range)',   exc(3) <= 1.0e-12_dblp, exc(3))
        call check('C3 smooth ratio, local bound per step (excursion / range)', lexc_step(2) <= 1.0e-12_dblp,        &
                   lexc_step(2))
        call check('C3 step ratio, local bound per step (excursion / range)',   lexc_step(3) <= 1.0e-12_dblp,        &
                   lexc_step(3))
        call check('C3 smooth ratio, local bound free run (excursion / range)', lexc_free(2) <= 1.0e-12_dblp,        &
                   lexc_free(2))
        call check('C3 step ratio, local bound free run (excursion / range)',   lexc_free(3) <= 1.0e-12_dblp,        &
                   lexc_free(3))
        call check('C4 isotope mass, uniform ratio (max rel. change)', dmass(1) <= 1.0e-12_dblp, dmass(1))
        call check('C4 isotope mass, smooth ratio (max rel. change)',  dmass(2) <= 1.0e-6_dblp, dmass(2))
        call check('C4 isotope mass, step ratio (max rel. change)',    dmass(3) <= 1.0e-6_dblp, dmass(3))
        call check('C5 iwat16 density == water (cells differing)', nmis_16 == 0, real(nmis_16, dblp))

        if (nfail > 0) then
          write(*,'(a,i0,a)') 'fvt_ratiocheck: ', nfail, ' FAILED'
          error stop 1
        endif
        write(*,'(a)') 'fvt_ratiocheck: all passed'

      contains

        subroutine check(name, ok, val)

          character(len=*), intent(in) :: name
          logical,          intent(in) :: ok
          real(dblp),       intent(in) :: val

          if (ok) then
            write(*,'(a,a,es12.4)') 'PASS  ', name, val
          else
            write(*,'(a,a,es12.4)') 'FAIL  ', name, val
            nfail = nfail + 1
          endif

        end subroutine check

        ! Largest distance of rnew outside the range of rold over the cells that can reach each cell in one step
        function local_excess(rold, rnew, ns) result(e)

          real(dblp),  intent(in) :: rold(:,:)
          real(dblp),  intent(in) :: rnew(:,:)
          integer(ip), intent(in) :: ns
          real(dblp)              :: e

          integer(ip) :: i, j, hj, ii, m
          real(dblp)  :: lo, hi

          e = 0.0_dblp
          do i = 1, nlat
            hj = ns + ceiling(maxval(abs(fx(i, :)))/g%area(i), ip)
            do j = 1, nlon
              lo = rold(i, j)
              hi = rold(i, j)
              do ii = max(1_ip, i - ns), min(nlat, i + ns)
                do m = -hj, hj
                  lo = min(lo, rold(ii, modulo(j - 1 + m, nlon) + 1))
                  hi = max(hi, rold(ii, modulo(j - 1 + m, nlon) + 1))
                enddo
              enddo
              e = max(e, lo - rnew(i, j), rnew(i, j) - hi)
            enddo
          enddo

        end function local_excess

      end program fvt_ratiocheck

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
