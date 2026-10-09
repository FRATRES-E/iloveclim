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
!      PROGRAM: [fvt_suite]
!
!>     @brief Evaluation suite of the FV transport (ROADMAP step 1.1, patch wiso:1c:0001).
!
!      DESCRIPTION:
!>     Usage: fvt_suite nlat nlon nstep timing outdir
!!       nlat, nlon : Gaussian grid (nlat even, nlon even)
!!       nstep      : time steps per period (12 days; even); the ECBilt step (4 h) is nstep = 72
!!       timing     : 'mid' (winds at mid-step, second order) or 'start' (winds at the start of the step, as the
!!                    model provides them)
!!       outdir     : output directory (must exist)
!!     Cases (exact solution at t = T is the initial state in all of them):
!!       A. Williamson et al. (1992) solid-body rotation of a cosine bell, alpha = 0, pi/4, pi/2 (density mode);
!!       B. Lauritzen et al. (2012) non-divergent deformational flow: Gaussian hills, cosine bells, slotted cylinders,
!!          correlated cosine bells (density mode = mixing ratio, the density stays 1);
!!       C. Lauritzen et al. (2012) divergent flow, density 1 at t = 0: the same fields as mixing ratios, transported
!!          either as independent densities rho and rho*phi (option C, phi = rho*phi / rho) or as ratios to the
!!          carrier rho (option C').
!!     Schemes: PPM unlimited, monotone, positive-definite, van Leer; in C' also positive-definite carrier with
!!     monotone ratios.
!!     Outputs: outdir/summary.csv (one line per metric, appended), and one netCDF file per case, mode and scheme
!!     with the fields at t = 0, T/2, T. Stops with exit code 1 if mass is not conserved to round-off.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      program fvt_suite

        use global_constants_mod, only: ip, dblp=>dp, PI=>pi_dp, str_len
        use fvt_grid_mod,         only: fvt_grid_t, fvt_grid_init
        use fvt_ffsl_mod,         only: fvt_ffsl_density, fvt_ffsl_ratio, FVT_RECON_PPM, FVT_RECON_VL,                  &
                                        FVT_LIM_NONE, FVT_LIM_MONO, FVT_LIM_POSDEF
        use fvt_testcases_mod,    only: TC_PERIOD, TC_FLOW_SBR, TC_FLOW_DEFORM, TC_FLOW_DIVERGENT, TC_IC_WG_BELL,       &
                                        TC_IC_HILLS, TC_IC_BELLS, TC_IC_CYLINDERS, TC_IC_CORRELATED, tc_norms_t,         &
                                        tc_fluxes, tc_cell_means, tc_norms, tc_filament, tc_mixing, tc_total
        use fvt_ncout_mod,        only: nco_write_fields

        implicit none

        integer(ip), parameter :: NSCHEME = 5
        integer(ip), parameter :: S_RECON(NSCHEME) = [FVT_RECON_PPM, FVT_RECON_PPM, FVT_RECON_PPM, FVT_RECON_VL, FVT_RECON_PPM]
        integer(ip), parameter :: S_LIM(NSCHEME)   = [FVT_LIM_NONE, FVT_LIM_MONO, FVT_LIM_POSDEF, FVT_LIM_MONO, FVT_LIM_POSDEF]
        integer(ip), parameter :: S_LIMR(NSCHEME)  = [FVT_LIM_NONE, FVT_LIM_MONO, FVT_LIM_POSDEF, FVT_LIM_MONO, FVT_LIM_MONO]
        character(len=9), parameter :: S_NAME(NSCHEME) = ['ppm_none ', 'ppm_mono ', 'ppm_pdef ', 'vanleer  ', 'pdef_mono']
        integer(ip), parameter :: NTAU = 19
        real(dblp),  parameter :: TOL_MASS = 1.0e-12_dblp       !< relative mass conservation tolerance

        type(fvt_grid_t)       :: g
        integer(ip)            :: nlat, nlon, nstep, ucsv, nfail, it
        real(dblp)             :: dt, tau(NTAU)
        logical                :: midstep
        character(len=str_len) :: arg, outdir, timing, prefix

        call read_args()
        call fvt_grid_init(g, nlat, nlon)
        dt = TC_PERIOD/real(nstep, dblp)
        do it = 1, NTAU
          tau(it) = 0.10_dblp + 0.05_dblp*real(it - 1, dblp)
        enddo
        nfail = 0
        call open_csv()
        write(prefix, '(i0,a,i0,a,i0,a,a)') nlat, ',', nlon, ',', nstep, ',', trim(timing)

        write(*, '(a,i0,a,i0,a,i0,a,a)') 'fvt_suite: grid ', nlat, ' x ', nlon, ', ', nstep, ' steps per period, winds at ', &
                                         trim(timing)
        call case_sbr()
        call case_deform()
        call case_divergent()
        close(ucsv)

        if (nfail > 0) then
          write(*, '(i0,a)') nfail, ' mass conservation failure(s)'
          error stop 1
        endif
        write(*, '(a)') 'fvt_suite: done'

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Arguments and output
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine read_args()

          integer(ip) :: ios

          if (command_argument_count() /= 5) then
            write(*, '(a)') 'usage: fvt_suite nlat nlon nstep mid|start outdir'
            error stop 2
          endif
          call get_command_argument(1, arg)
          read(arg, *, iostat=ios) nlat
          if (ios /= 0) error stop 'fvt_suite: bad nlat'
          call get_command_argument(2, arg)
          read(arg, *, iostat=ios) nlon
          if (ios /= 0) error stop 'fvt_suite: bad nlon'
          call get_command_argument(3, arg)
          read(arg, *, iostat=ios) nstep
          if (ios /= 0 .or. nstep < 2 .or. mod(nstep, 2_ip) /= 0) error stop 'fvt_suite: nstep must be even and >= 2'
          call get_command_argument(4, timing)
          if (trim(timing) /= 'mid' .and. trim(timing) /= 'start') error stop 'fvt_suite: timing must be mid or start'
          midstep = trim(timing) == 'mid'
          call get_command_argument(5, outdir)

        end subroutine read_args

        subroutine open_csv()

          logical :: exists

          inquire(file=trim(outdir)//'/summary.csv', exist=exists)
          open(newunit=ucsv, file=trim(outdir)//'/summary.csv', position='append', action='write')
          if (.not. exists) write(ucsv, '(a)') 'nlat,nlon,nstep,timing,case,mode,scheme,tracer,metric,value'

        end subroutine open_csv

        subroutine csv(cas, mode, scheme, tracer, metric, value)

          character(len=*), intent(in) :: cas
          character(len=*), intent(in) :: mode
          character(len=*), intent(in) :: scheme
          character(len=*), intent(in) :: tracer
          character(len=*), intent(in) :: metric
          real(dblp),       intent(in) :: value

          write(ucsv, '(12a,es22.14e3)') trim(prefix), ',', trim(cas), ',', trim(mode), ',', trim(scheme), ',', &
                                          trim(tracer), ',', trim(metric), ',', value

        end subroutine csv

        subroutine csv_norms(cas, mode, scheme, tracer, nrm)

          character(len=*), intent(in) :: cas
          character(len=*), intent(in) :: mode
          character(len=*), intent(in) :: scheme
          character(len=*), intent(in) :: tracer
          type(tc_norms_t), intent(in) :: nrm

          call csv(cas, mode, scheme, tracer, 'l1', nrm%l1)
          call csv(cas, mode, scheme, tracer, 'l2', nrm%l2)
          call csv(cas, mode, scheme, tracer, 'linf', nrm%linf)
          call csv(cas, mode, scheme, tracer, 'phimin', nrm%phimin)
          call csv(cas, mode, scheme, tracer, 'phimax', nrm%phimax)

        end subroutine csv_norms

        subroutine check_mass(label, m, m0)

          character(len=*), intent(in) :: label
          real(dblp),       intent(in) :: m
          real(dblp),       intent(in) :: m0

          if (abs(m - m0)/abs(m0) > TOL_MASS) then
            write(*, '(a,a,es10.3)') 'FAIL mass conservation ', label, abs(m - m0)/abs(m0)
            nfail = nfail + 1
          endif

        end subroutine check_mass

        function ncname(cas, mode, scheme) result(path)

          character(len=*), intent(in) :: cas
          character(len=*), intent(in) :: mode
          character(len=*), intent(in) :: scheme
          character(len=str_len)       :: path

          write(path, '(a,a,a,a,a,a,a,a,i0,a,a,a)') trim(outdir), '/', trim(cas), '_', trim(mode), '_', trim(scheme), &
                                                    '_n', nlat, '_', trim(timing), '.nc'

        end function ncname

        ! Time at which the winds of step n (0-based) are evaluated
        pure function t_wind(n) result(t)

          integer(ip), intent(in) :: n
          real(dblp)              :: t

          t = real(n, dblp)*dt
          if (midstep) t = t + 0.5_dblp*dt

        end function t_wind

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   A. Solid-body rotation
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine case_sbr()

          real(dblp), parameter :: ALPHAS(3) = [0.0_dblp, 0.25_dblp*PI, 0.5_dblp*PI]
          character(len=8), parameter :: CNAME(3) = ['sbr_a00 ', 'sbr_a45 ', 'sbr_a90 ']
          real(dblp)  :: fx(nlat, nlon), fy(0:nlat, nlon), q(nlat, nlon, 1), q0(nlat, nlon), fld(nlat, nlon, 3, 1)
          integer(ip) :: ia, is, n, ns, nsmax

          call tc_cell_means(g, TC_IC_WG_BELL, q0)
          do ia = 1, 3
            do is = 1, NSCHEME - 1
              q(:, :, 1) = q0
              fld(:, :, 1, 1) = q0
              nsmax = 0
              do n = 0, nstep - 1
                call tc_fluxes(g, TC_FLOW_SBR, ALPHAS(ia), t_wind(n), dt, fx, fy)
                call fvt_ffsl_density(g, fx, fy, q, S_RECON(is), S_LIM(is), ns)
                nsmax = max(nsmax, ns)
                if (n + 1 == nstep/2) fld(:, :, 2, 1) = q(:, :, 1)
              enddo
              fld(:, :, 3, 1) = q(:, :, 1)
              call check_mass(trim(CNAME(ia))//' '//S_NAME(is), tc_total(g, q(:, :, 1)), tc_total(g, q0))
              call csv_norms(CNAME(ia), 'density', S_NAME(is), 'wg_bell', tc_norms(g, q(:, :, 1), q0, 1.0_dblp))
              call csv(CNAME(ia), 'density', S_NAME(is), 'wg_bell', 'nsub_max', real(nsmax, dblp))
              call nco_write_fields(ncname(CNAME(ia), 'density', S_NAME(is)), g, ['wg_bell'], fld,                       &
                                    [0.0_dblp, 0.5_dblp*TC_PERIOD, TC_PERIOD], 'Williamson case 1, '//trim(CNAME(ia)))
            enddo
            write(*, '(a,a)') '  done ', trim(CNAME(ia))
          enddo

        end subroutine case_sbr

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   B. Non-divergent deformational flow
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine case_deform()

          integer(ip), parameter :: NT = 4
          integer(ip), parameter :: ICS(NT) = [TC_IC_HILLS, TC_IC_BELLS, TC_IC_CYLINDERS, TC_IC_CORRELATED]
          character(len=10), parameter :: TNAME(NT) = ['hills     ', 'bells     ', 'cylinders ', 'correlated']
          real(dblp)  :: fx(nlat, nlon), fy(0:nlat, nlon), q(nlat, nlon, NT), q0(nlat, nlon, NT), fld(nlat, nlon, 3, NT)
          real(dblp)  :: lf(NTAU), lr, lu, lo, tu(NT), to(NT)
          integer(ip) :: k, is, n, ns, nsmax
          character(len=16) :: met

          do k = 1, NT
            call tc_cell_means(g, ICS(k), q0(:, :, k))
          enddo
          do is = 1, NSCHEME - 1
            q = q0
            fld(:, :, 1, :) = q0
            nsmax = 0
            tu = 0.0_dblp
            to = 0.0_dblp
            do n = 0, nstep - 1
              call tc_fluxes(g, TC_FLOW_DEFORM, 0.0_dblp, t_wind(n), dt, fx, fy)
              call fvt_ffsl_density(g, fx, fy, q, S_RECON(is), S_LIM(is), ns)
              nsmax = max(nsmax, ns)
              if (n + 1 == nstep/2) fld(:, :, 2, :) = q
              call track_bounds(q, q0, tu, to)
            enddo
            fld(:, :, 3, :) = q
            do k = 1, NT
              call check_mass('deform '//S_NAME(is)//' '//TNAME(k), tc_total(g, q(:, :, k)), tc_total(g, q0(:, :, k)))
              call csv_norms('deform', 'density', S_NAME(is), TNAME(k),                                              &
                             tc_norms(g, q(:, :, k), q0(:, :, k), maxval(q0(:, :, k)) - minval(q0(:, :, k))))
              call csv('deform', 'density', S_NAME(is), TNAME(k), 'trans_under', tu(k))
              call csv('deform', 'density', S_NAME(is), TNAME(k), 'trans_over', to(k))
            enddo
            call csv('deform', 'density', S_NAME(is), 'all', 'nsub_max', real(nsmax, dblp))
            ! filaments (cosine bells) and mixing at T/2
            call tc_filament(g, fld(:, :, 2, 2), q0(:, :, 2), tau, lf)
            do k = 1, NTAU
              write(met, '(a,f4.2)') 'lf_', tau(k)
              call csv('deform', 'density', S_NAME(is), 'bells', trim(met), lf(k))
            enddo
            call tc_mixing(g, fld(:, :, 2, 2), fld(:, :, 2, 4), lr, lu, lo)
            call csv('deform', 'density', S_NAME(is), 'pair', 'lr', lr)
            call csv('deform', 'density', S_NAME(is), 'pair', 'lu', lu)
            call csv('deform', 'density', S_NAME(is), 'pair', 'lo', lo)
            call nco_write_fields(ncname('deform', 'density', S_NAME(is)), g, TNAME, fld,                             &
                                  [0.0_dblp, 0.5_dblp*TC_PERIOD, TC_PERIOD], 'Lauritzen et al. 2012, non-divergent flow')
            write(*, '(a,a)') '  done deform ', trim(S_NAME(is))
          enddo

        end subroutine case_deform

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   C. Divergent flow: option C (independent densities) and option C' (carrier and ratios)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine case_divergent()

          integer(ip), parameter :: NT = 3
          integer(ip), parameter :: ICS(NT) = [TC_IC_BELLS, TC_IC_CYLINDERS, TC_IC_CORRELATED]
          character(len=10), parameter :: TNAME(NT+1) = ['rho       ', 'bells     ', 'cylinders ', 'correlated']
          real(dblp)  :: fx(nlat, nlon), fy(0:nlat, nlon)
          real(dblp)  :: phi0(nlat, nlon, NT), d(nlat, nlon, NT+1), qc(nlat, nlon), r(nlat, nlon, NT)
          real(dblp)  :: fld(nlat, nlon, 3, NT+1), m0(NT+1), rho0(nlat, nlon), lr, lu, lo
          real(dblp)  :: phit(nlat, nlon, NT+1), tu(NT+1), to(NT+1)
          integer(ip) :: k, is, n, ns, nsmax, imode
          character(len=7) :: mode

          do k = 1, NT
            call tc_cell_means(g, ICS(k), phi0(:, :, k))
          enddo
          rho0 = 1.0_dblp
          m0(1) = tc_total(g, rho0)
          do k = 1, NT
            m0(k+1) = tc_total(g, phi0(:, :, k))
          enddo

          do imode = 1, 2
            do is = 1, NSCHEME
              ! the pdef carrier with monotone ratios exists in ratio mode only
              if (imode == 1 .and. is == NSCHEME) cycle
              if (imode == 1) then
                mode = 'density'
                d(:, :, 1)  = rho0
                d(:, :, 2:) = phi0
              else
                mode = 'ratio'
                qc = rho0
                r  = phi0
              endif
              fld(:, :, 1, 1)  = rho0
              fld(:, :, 1, 2:) = phi0
              nsmax = 0
              tu = 0.0_dblp
              to = 0.0_dblp
              do n = 0, nstep - 1
                call tc_fluxes(g, TC_FLOW_DIVERGENT, 0.0_dblp, t_wind(n), dt, fx, fy)
                if (imode == 1) then
                  call fvt_ffsl_density(g, fx, fy, d, S_RECON(is), S_LIM(is), ns)
                else
                  call fvt_ffsl_ratio(g, fx, fy, qc, r, S_RECON(is), S_LIM(is), ns, limiter_ratio=S_LIMR(is))
                endif
                nsmax = max(nsmax, ns)
                if (n + 1 == nstep/2) call store_div(imode, d, qc, r, fld(:, :, 2, :))
                call store_div(imode, d, qc, r, phit)
                call track_bounds(phit(:, :, 2:), fld(:, :, 1, 2:), tu(2:), to(2:))
              enddo
              call store_div(imode, d, qc, r, fld(:, :, 3, :))
              ! conservation of rho and of rho*phi
              do k = 1, NT + 1
                if (k == 1) then
                  call check_mass('divergent '//trim(mode)//' '//S_NAME(is)//' rho', tc_total(g, fld(:, :, 3, 1)), m0(1))
                else
                  call check_mass('divergent '//trim(mode)//' '//S_NAME(is)//' '//TNAME(k),                          &
                                  tc_total(g, fld(:, :, 3, 1)*fld(:, :, 3, k)), m0(k))
                endif
              enddo
              do k = 1, NT + 1
                call csv_norms('divergent', mode, S_NAME(is), TNAME(k),                                              &
                               tc_norms(g, fld(:, :, 3, k), fld(:, :, 1, k),                                         &
                                        max(maxval(fld(:, :, 1, k)) - minval(fld(:, :, 1, k)), 1.0_dblp)))
              enddo
              do k = 2, NT + 1
                call csv('divergent', mode, S_NAME(is), TNAME(k), 'trans_under', tu(k))
                call csv('divergent', mode, S_NAME(is), TNAME(k), 'trans_over', to(k))
              enddo
              call csv('divergent', mode, S_NAME(is), 'all', 'nsub_max', real(nsmax, dblp))
              call tc_mixing(g, fld(:, :, 2, 2), fld(:, :, 2, 4), lr, lu, lo)
              call csv('divergent', mode, S_NAME(is), 'pair', 'lr', lr)
              call csv('divergent', mode, S_NAME(is), 'pair', 'lu', lu)
              call csv('divergent', mode, S_NAME(is), 'pair', 'lo', lo)
              call nco_write_fields(ncname('divergent', mode, S_NAME(is)), g, TNAME, fld,                             &
                                    [0.0_dblp, 0.5_dblp*TC_PERIOD, TC_PERIOD],                                       &
                                    'Lauritzen et al. 2012, divergent flow, mixing ratios ('//trim(mode)//' mode)')
              write(*, '(a,a,a,a)') '  done divergent ', trim(mode), ' ', trim(S_NAME(is))
            enddo
          enddo

        end subroutine case_divergent

        ! Largest excursions so far outside the initial range of each field, relative to that range
        subroutine track_bounds(q, q0, tu, to)

          real(dblp), intent(in)    :: q(:,:,:)
          real(dblp), intent(in)    :: q0(:,:,:)
          real(dblp), intent(inout) :: tu(:)      !< (min q0 - min q) / range, maximum over time (0 if none)
          real(dblp), intent(inout) :: to(:)      !< (max q - max q0) / range, maximum over time (0 if none)

          integer(ip) :: k
          real(dblp)  :: rng

          do k = 1, size(q, 3)
            rng   = maxval(q0(:, :, k)) - minval(q0(:, :, k))
            tu(k) = max(tu(k), (minval(q0(:, :, k)) - minval(q(:, :, k)))/rng)
            to(k) = max(to(k), (maxval(q(:, :, k)) - maxval(q0(:, :, k)))/rng)
          enddo

        end subroutine track_bounds

        ! Density and mixing ratios of the divergent case, from either mode
        subroutine store_div(imode, d, qc, r, out)

          integer(ip), intent(in)  :: imode
          real(dblp),  intent(in)  :: d(:,:,:)      !< density mode: rho, rho*phi_k
          real(dblp),  intent(in)  :: qc(:,:)       !< ratio mode: carrier rho
          real(dblp),  intent(in)  :: r(:,:,:)      !< ratio mode: phi_k
          real(dblp),  intent(out) :: out(:,:,:)    !< rho, phi_k

          integer(ip) :: k

          if (imode == 1) then
            out(:, :, 1) = d(:, :, 1)
            do k = 2, size(d, 3)
              out(:, :, k) = d(:, :, k)/d(:, :, 1)
            enddo
          else
            out(:, :, 1)  = qc
            out(:, :, 2:) = r
          endif

        end subroutine store_div

      end program fvt_suite

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
