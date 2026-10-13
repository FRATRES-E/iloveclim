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
!!     wiso:1e:0001 provides the switch and the set-up only: adv_mode = 'fv' builds and checks the tables, then stops
!!     (the transport step comes with wiso:1e:0002). The default 'spectral' leaves the model unchanged.
!!
!!     Contract (public entry points):
!!       adv_mode                 : namelist variable (moipar), default 'spectral'.
!!       atm_fvmois_check_mode()  : validate adv_mode, after the namelist is read; stops on an unknown value.
!!       atm_fvmois_init()        : FV mode only: build and check the FV grid and flux tables (needs coef.dat read).
!!       atm_fvmois_is_on()       : .true. in FV mode.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module atm_fvmois_mod

        use global_constants_mod, only: ip, dblp=>dp
        use fvt_grid_mod,         only: fvt_grid_t, fvt_grid_init
        use fvt_specflux_mod,     only: fvt_specflux_t, fvt_specflux_init

        implicit none

        private

        public :: ADV_MODE_LEN, adv_mode, atm_fvmois_check_mode, atm_fvmois_init, atm_fvmois_is_on


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Constants and module state
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip),          parameter :: ADV_MODE_LEN = 16               !< length of adv_mode
        real(dblp),           parameter :: FVM_TOL_WGT  = 1.0e-12_dblp     !< weights vs coef.dat (file precision 2e-13)
        real(dblp),           parameter :: FVM_TOL_MU   = 1.0e-9_dblp      !< sin(latitude) vs gauss.asc (10 decimals)

        character(len=ADV_MODE_LEN)     :: adv_mode     = 'spectral'       !< moisture transport: 'spectral' | 'fv'
        logical,              save      :: fvm_on       = .false.          !< .true. in FV mode (after check_mode)
        type(fvt_grid_t),     save      :: fvm_grid                        !< FV cells on the model Gaussian grid
        type(fvt_specflux_t), save      :: fvm_sf                          !< spectral psi, chi -> face flux tables


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

          use comatm, only: nlat, nlon, nm, phi
          use comdyn, only: pw

          real(dblp) :: err_wgt, err_mu

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

          write(*,*) "atm_fvmois: adv_mode = fv is not available yet (transport step in wiso:1e:0002)"
          stop 1

        end subroutine atm_fvmois_init

      end module atm_fvmois_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
