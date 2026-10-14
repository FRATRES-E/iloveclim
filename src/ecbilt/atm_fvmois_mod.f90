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
!      MODULE: [atm_fvmois_mod]
!
!>     @brief ECBilt side of the finite-volume moisture transport (ROADMAP step 1.2, option B): mode switch and set-up.
!
!      DESCRIPTION:
!>     Holds the transport mode of the water vapour, read in the ECBilt namelist group moipar:
!!         adv_mode = 'spectral'   current chain (round trip, clipping, ec_trafluxdiv, ec_hdiff, ec_moisbalance)
!!         adv_mode = 'fv'         flux-form semi-Lagrangian transport of src/fvtrans (D20-D24), fluxes from the
!!                                 spectral psi(:,3) and chi(:,3) times umoisr (D22)
!!     and, in FV mode, the FV grid and the spectral-to-flux tables built once from the model's own dimensions and
!!     truncation. At set-up the FV cells are checked against the model grid: Gaussian weights against pw(:,1) of
!!     coef.dat (to its precision) and the Gaussian latitudes against those of gauss.asc (10 decimals, O25).
!!
!!     In FV mode, ec_moisture calls atm_fvmois_step on the column water instead of the spectral chain: one
!!     positive-definite PPM step (D23, out_max = 1, D26) with the winds of the current step (O18), no explicit
!!     diffusion. Each step is checked for conservation and positivity; a summary (substeps, worst mass error,
!!     minimum) is written to standard output once per model year of steps. The default 'spectral' leaves the
!!     model unchanged.
!!
!!     Contract (public entry points):
!!       adv_mode                 : namelist variable (moipar), default 'spectral'.
!!       atm_fvmois_check_mode()  : validate adv_mode, after the namelist is read; stops on an unknown value.
!!       atm_fvmois_init()        : FV mode only: build and check the FV grid and flux tables (needs coef.dat read).
!!       atm_fvmois_is_on()       : .true. in FV mode.
!!       atm_fvmois_step(q)       : FV mode only: transport q(nlat,nlon) over one dtime (column amount, density mode).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module atm_fvmois_mod

        use global_constants_mod, only: ip, dblp=>dp
        use fvt_grid_mod,         only: fvt_grid_t, fvt_grid_init
        use fvt_specflux_mod,     only: fvt_specflux_t, fvt_specflux_init, fvt_specflux_fluxes
        use fvt_ffsl_mod,         only: fvt_ffsl_density, FVT_RECON_PPM, FVT_LIM_POSDEF

        implicit none

        private

        public :: ADV_MODE_LEN, adv_mode, atm_fvmois_check_mode, atm_fvmois_init, atm_fvmois_is_on, atm_fvmois_step


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Constants and module state
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip),          parameter :: ADV_MODE_LEN = 16               !< length of adv_mode
        real(dblp),           parameter :: FVM_TOL_WGT  = 1.0e-12_dblp     !< weights vs coef.dat (file precision 2e-13)
        real(dblp),           parameter :: FVM_TOL_MU   = 1.0e-9_dblp      !< sin(latitude) vs gauss.asc (10 decimals)
        real(dblp),           parameter :: FVM_TOL_MASS = 1.0e-12_dblp     !< relative global mass change per step, and
                                                                           !< negatives relative to the global mean
        integer(ip),          parameter :: FVM_MAX_WARN = 10               !< warnings written before going silent

        character(len=ADV_MODE_LEN)     :: adv_mode     = 'spectral'       !< moisture transport: 'spectral' | 'fv'
        logical,              save      :: fvm_on       = .false.          !< .true. in FV mode (after check_mode)
        type(fvt_grid_t),     save      :: fvm_grid                        !< FV cells on the model Gaussian grid
        type(fvt_specflux_t), save      :: fvm_sf                          !< spectral psi, chi -> face flux tables
        real(dblp), allocatable, save   :: fvm_fx(:,:)                     !< (nlat,nlon) east-face fluxes
        real(dblp), allocatable, save   :: fvm_fy(:,:)                     !< (0:nlat,nlon) north-face fluxes
        real(dblp), allocatable, save   :: fvm_q(:,:,:)                    !< (nlat,nlon,1) work copy of the tracer
        real(dblp), allocatable, save   :: fvm_area(:,:)                   !< (nlat,nlon) cell areas (unit sphere)

! dmr&clo   Diagnostics, reset after each summary
        integer(ip),          save      :: fvm_nstep    = 0                !< steps since the last summary
        integer(ip),          save      :: fvm_nsub_sum = 0                !< sum of the substeps
        integer(ip),          save      :: fvm_nsub_max = 0                !< largest number of substeps
        real(dblp),           save      :: fvm_dm_max   = 0.0_dblp         !< largest relative mass change
        real(dblp),           save      :: fvm_q_min    = huge(1.0_dblp)   !< smallest value after transport
        integer(ip),          save      :: fvm_nwarn    = 0                !< warnings written so far
        integer(ip),          save      :: fvm_report   = 0                !< steps per summary (one model year)


      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Mode switch
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine atm_fvmois_check_mode()

          select case (trim(adv_mode))
          case ('spectral')
            fvm_on = .false.
          case ('fv')
            fvm_on = .true.
          case default
            write(*,*) "atm_fvmois: unknown adv_mode = ", trim(adv_mode), " in namelist group moipar (spectral | fv)"
            stop 1
          end select
          write(*,*) "atm_fvmois: moisture transport adv_mode = ", trim(adv_mode)

        end subroutine atm_fvmois_check_mode

        function atm_fvmois_is_on() result(on)

          logical :: on

          on = fvm_on

        end function atm_fvmois_is_on

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Set-up of the FV grid and flux tables, checked against the model grid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine atm_fvmois_init()

          use comatm, only: nlat, nlon, nm, phi, dtime
          use comdyn, only: pw

          real(dblp)  :: err_wgt, err_mu
          integer(ip) :: j

          if (.not. fvm_on) return

          call fvt_grid_init(fvm_grid, nlat, nlon)
          call fvt_specflux_init(fvm_sf, fvm_grid, nm)

          ! pw(:,1) is (w/2) * P(0,0) / sqrt(nlon) after the scaling in ec_iatmdyn, with P(0,0) = 1
          err_wgt = maxval(abs(fvm_grid%wgt(:) - 2.0_dblp*sqrt(real(nlon, dblp))*pw(:,1)))
          err_mu  = maxval(abs(fvm_grid%mu(:) - sin(phi(:))))
          write(*,'(a,es10.2,a,es10.2)') " atm_fvmois: FV grid vs model, max |dw| = ", err_wgt,                     &
                                         ", max |d sin(lat)| = ", err_mu
          if (err_wgt > FVM_TOL_WGT .or. err_mu > FVM_TOL_MU) then
            write(*,*) "atm_fvmois: FV grid does not match the model Gaussian grid (coef.dat, gauss.asc)"
            stop 1
          endif

          allocate(fvm_fx(nlat, nlon), fvm_fy(0:nlat, nlon), fvm_q(nlat, nlon, 1), fvm_area(nlat, nlon))
          do j = 1, nlon
            fvm_area(:, j) = fvm_grid%area(:)
          enddo
          fvm_report = nint(360.0_dblp*86400.0_dblp/dtime)

        end subroutine atm_fvmois_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   One transport step of the column water, with the winds of the current step
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine atm_fvmois_step(q)

          use comatm,  only: nsh, dtime, om, radius
          use comdyn,  only: psi, chi
          use comphys, only: umoisr

          real(dblp), intent(inout) :: q(:,:)   !< (nlat,nlon) moisture (column amount), positive on entry

          real(dblp)  :: mass0, mass1, dm, qmin, qneg
          integer(ip) :: nsub

          if (.not. fvm_on) error stop 'atm_fvmois_step: called with adv_mode /= fv'

          ! face fluxes of the 800 hPa wind reduced by umoisr (D22), ec_trafluxdiv's wind
          call fvt_specflux_fluxes(fvm_sf, fvm_grid, reshape(psi(:,3), [nsh, 2]), reshape(chi(:,3), [nsh, 2]),      &
                                   umoisr*dtime*om, umoisr*dtime/radius**2, fvm_fx, fvm_fy)

          fvm_q(:,:,1) = q(:,:)
          mass0 = sum(fvm_area*q)
          call fvt_ffsl_density(fvm_grid, fvm_fx, fvm_fy, fvm_q, FVT_RECON_PPM, FVT_LIM_POSDEF, nsub)
          q(:,:) = fvm_q(:,:,1)
          mass1 = sum(fvm_area*q)

          dm = 0.0_dblp
          if (mass0 > 0.0_dblp) dm = abs(mass1 - mass0)/mass0
          qmin = minval(q)
          qneg = -FVM_TOL_MASS*mass0/sum(fvm_area)   ! round-off negatives (about -1e-22) are not reported

          fvm_nstep    = fvm_nstep + 1
          fvm_nsub_sum = fvm_nsub_sum + nsub
          fvm_nsub_max = max(fvm_nsub_max, nsub)
          fvm_dm_max   = max(fvm_dm_max, dm)
          fvm_q_min    = min(fvm_q_min, qmin)

          if ((dm > FVM_TOL_MASS .or. qmin < qneg) .and. fvm_nwarn < FVM_MAX_WARN) then
            fvm_nwarn = fvm_nwarn + 1
            write(*,'(a,es10.2,a,es10.2)') " atm_fvmois: WARNING relative mass change ", dm, ", minimum ", qmin
            if (fvm_nwarn == FVM_MAX_WARN) write(*,*) "atm_fvmois: further warnings suppressed"
          endif

          if (fvm_nstep >= fvm_report) then
            write(*,'(a,i0,a,f6.3,a,i0,a,es10.2,a,es10.2)') " atm_fvmois: ", fvm_nstep, " steps, substeps mean ",     &
                  real(fvm_nsub_sum, dblp)/real(fvm_nstep, dblp), " max ", fvm_nsub_max, ", max rel. mass change ",     &
                  fvm_dm_max, ", min ", fvm_q_min
            fvm_nstep    = 0
            fvm_nsub_sum = 0
            fvm_nsub_max = 0
            fvm_dm_max   = 0.0_dblp
            fvm_q_min    = huge(1.0_dblp)
          endif

        end subroutine atm_fvmois_step

      end module atm_fvmois_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
