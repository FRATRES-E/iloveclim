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
!      MODULE: [atm_wbudget_mod]
!
!>     @brief Diagnostics of the atmospheric water budget of ECBilt (ROADMAP step 1.3): global budget and clipping maps.
!
!      DESCRIPTION:
!>     Pure diagnostics: nothing here changes the model state (results are bit-identical with and without them).
!!     Every atmospheric step, between the start of ec_moisture and the end of ec_convec_two, the global integrals of
!!         W    column water (rmoisg) at the start and at the end,
!!         E    evaporation applied (imsink * evap),
!!         P    precipitation applied (ivavm * (dyrain + dysnow + corain + cosnow)), large-scale and convective,
!!         T    horizontal transport applied as a tendency (-ihavm * hdivmg; zero with adv_mode = 'fv'),
!!         D    horizontal diffusion (hdmoisg; zero with adv_mode = 'fv'),
!!         C1   water added by the clipping in ec_moisture (cormois before ec_moisbalance),
!!         C2   water removed by ec_moisbalance,
!!         C3   water added by the clipping in ec_convec_two,
!!     are accumulated, and the residual  R = W_end - W_start - (E - P + T + D + C1 - C2 + C3)  measures everything
!!     not accounted for (the FV transport, which conserves mass to round-off, and any non-conservation of the
!!     spectral round trip or of ec_moisbalance end up there). All terms are global means in mm of water.
!!     One line per model month is written to outputdata/globals/wbudget<fini>.txt (terms summed over the month,
!!     W at the end of the month), at the last step of the month (30 days of iatm steps); an incomplete month (start
!!     or restart in the middle of a month) is written when the calendar month changes.
!!
!!     The clipping maps go to the user field user_t1 (dumt1, otherwise unused): level 1 = C1, level 2 = C2,
!!     level 3 = C3, as rates in cm/yr like the precipitation outputs, averaged over the last iatm steps (one day)
!!     so that the daily sampling of the outputs sees a daily mean.
!!
!!     Contract (public entry points), all called from ECBilt:
!!       atm_wb_step_begin()                   : before ec_moisture (month bookkeeping, W at the start).
!!       atm_wb_transport(hdivmg, hdmoisg)     : in ec_moisture after the forward step, the applied T and D tendencies.
!!       atm_wb_moisbalance_pre()              : just before ec_moisbalance (C1, state before).
!!       atm_wb_moisbalance_post()             : just after ec_moisbalance (C2).
!!       wb_clip_conv(nlat,nlon)               : public accumulator, incremented by ec_convec_two where it clips (C3).
!!       atm_wb_step_end()                     : after ec_convec_two (E, P, W at the end, residual, maps, output).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module atm_wbudget_mod

        use global_constants_mod, only: ip, dblp=>dp, days_month360d_i
        use comatm,               only: nlat, nlon

        implicit none

        private

        public :: wb_clip_conv
        public :: atm_wb_step_begin, atm_wb_transport, atm_wb_moisbalance_pre, atm_wb_moisbalance_post, atm_wb_step_end


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Constants and module state
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip), parameter :: WB_NTERM = 11                !< number of accumulated terms
        integer(ip), parameter :: WB_E   = 1                   !< evaporation
        integer(ip), parameter :: WB_PLS = 2                   !< large-scale precipitation (dyrain + dysnow)
        integer(ip), parameter :: WB_PCV = 3                   !< convective precipitation (corain + cosnow)
        integer(ip), parameter :: WB_T   = 4                   !< transport tendency
        integer(ip), parameter :: WB_D   = 5                   !< diffusion tendency
        integer(ip), parameter :: WB_C1  = 6                   !< clipping in ec_moisture
        integer(ip), parameter :: WB_C2  = 7                   !< removal by ec_moisbalance
        integer(ip), parameter :: WB_C3  = 8                   !< clipping in ec_convec_two
        integer(ip), parameter :: WB_R   = 9                   !< residual
        integer(ip), parameter :: WB_RAB = 10                  !< sum of |residual| per step
        integer(ip), parameter :: WB_DW  = 11                  !< W_end - W_start
        real(dblp),  parameter :: WB_MM  = 1000.0_dblp         !< m -> mm

        real(dblp), dimension(nlat,nlon), save :: wb_clip_conv = 0.0_dblp   !< C3 of the current step (m), set by convection

        logical,                            save :: wb_ready   = .false.    !< weights and buffers set
        real(dblp), dimension(nlat),        save :: wb_wgt                  !< latitude weights for global means
        real(dblp),                         save :: wb_w0      = 0.0_dblp   !< W at the start of the step (m)
        real(dblp), dimension(WB_NTERM),    save :: wb_step                 !< terms of the current step (m)
        real(dblp), dimension(WB_NTERM),    save :: wb_month   = 0.0_dblp   !< terms summed over the month (m)
        integer(ip),                        save :: wb_nstep   = 0          !< steps in the current month
        integer(ip),                        save :: wb_cur_month = -1       !< month being accumulated
        integer(ip),                        save :: wb_cur_year  = -1       !< year being accumulated
        real(dblp),                         save :: wb_w_last  = 0.0_dblp   !< W at the end of the last step (m)
        integer(ip),                        save :: wb_unit    = -1         !< output unit (newunit, negative once open)
        logical,                            save :: wb_open    = .false.    !< output file opened
        real(dblp), dimension(nlat,nlon),   save :: wb_map_c1, wb_map_c2    !< C1, C2 maps of the current step (m)
        real(dblp), dimension(nlat,nlon),   save :: wb_qpre                 !< rmoisg before ec_moisbalance
        real(dblp), allocatable,            save :: wb_ring(:,:,:,:)        !< (nlat,nlon,3,iatm) last iatm steps, cm/yr
        integer(ip),                        save :: wb_iring   = 0          !< last filled slot of the ring
        integer(ip),                        save :: wb_nring   = 0          !< filled slots (<= iatm)


      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Step bracket
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine atm_wb_step_begin()

          use comdyn,      only: pw
          use comphys,     only: rmoisg
          use comatm,      only: iwater
          use comemic_mod, only: imonth, iyear, iatm

          if (.not. wb_ready) then
            wb_wgt(:) = 2.0_dblp*sqrt(real(nlon, dblp))*pw(:,1)/(2.0_dblp*real(nlon, dblp))   ! sum over cells = 1
            allocate(wb_ring(nlat, nlon, 3, max(iatm, 1)))
            wb_ring = 0.0_dblp
            wb_ready = .true.
          endif

          ! the calendar month changed with steps still pending (a month that was not complete, e.g. after a restart
          ! in the middle of a month): write them; complete months are written at their last step (atm_wb_step_end)
          if (imonth /= wb_cur_month .or. iyear /= wb_cur_year) then
            if (wb_nstep > 0) call wb_write_month()
            wb_cur_month = imonth
            wb_cur_year  = iyear
            wb_month     = 0.0_dblp
            wb_nstep     = 0
          endif

          wb_step      = 0.0_dblp
          wb_map_c1    = 0.0_dblp
          wb_map_c2    = 0.0_dblp
          wb_clip_conv = 0.0_dblp
          wb_w0        = wb_gmean(rmoisg(:,:,iwater))

        end subroutine atm_wb_step_begin

        subroutine atm_wb_step_end()

          use comphys,     only: rmoisg, dyrain, dysnow, corain, cosnow, ivavm, imsink
          use comsurf_mod, only: evap
          use comatm,      only: iwater, dtime
          use comemic_mod, only: iatm

          real(dblp) :: w1, rfac

          w1 = wb_gmean(rmoisg(:,:,iwater))

          wb_step(WB_E)   = real(imsink, dblp)*dtime*wb_gmean(evap(:,:,iwater))
          wb_step(WB_PLS) = real(ivavm, dblp)*dtime*wb_gmean(dyrain(:,:,iwater) + dysnow(:,:,iwater))
          wb_step(WB_PCV) = real(ivavm, dblp)*dtime*wb_gmean(corain(:,:,iwater) + cosnow(:,:,iwater))
          wb_step(WB_C3)  = wb_gmean(wb_clip_conv)
          wb_step(WB_DW)  = w1 - wb_w0
          wb_step(WB_R)   = wb_step(WB_DW) - (wb_step(WB_E) - wb_step(WB_PLS) - wb_step(WB_PCV) + wb_step(WB_T)         &
                            + wb_step(WB_D) + wb_step(WB_C1) - wb_step(WB_C2) + wb_step(WB_C3))
          wb_step(WB_RAB) = abs(wb_step(WB_R))

          wb_month  = wb_month + wb_step
          wb_nstep  = wb_nstep + 1
          wb_w_last = w1

          ! last step of a complete month: write it now, so that the last month of a run is not lost
          if (wb_nstep == days_month360d_i*iatm) then
            call wb_write_month()
            wb_month = 0.0_dblp
            wb_nstep = 0
          endif

          ! clipping maps (cm/yr), mean over the last iatm steps, into user_t1
          rfac = 100.0_dblp*86400.0_dblp*360.0_dblp/dtime
          wb_iring = modulo(wb_iring, size(wb_ring, 4)) + 1
          wb_nring = min(wb_nring + 1, size(wb_ring, 4))
          wb_ring(:,:,1,wb_iring) = wb_map_c1*rfac
          wb_ring(:,:,2,wb_iring) = wb_map_c2*rfac
          wb_ring(:,:,3,wb_iring) = wb_clip_conv*rfac
          call wb_set_user_field()

        end subroutine atm_wb_step_end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Hooks inside ec_moisture
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine atm_wb_transport(hdivmg, hdmoisg)

          use comphys, only: ihavm
          use comatm,  only: dtime

          real(dblp), intent(in) :: hdivmg(:,:)    !< (nlat,nlon) horizontal flux divergence applied (m/s)
          real(dblp), intent(in) :: hdmoisg(:,:)   !< (nlat,nlon) horizontal diffusion applied (m/s)

          wb_step(WB_T) = -real(ihavm, dblp)*dtime*wb_gmean(hdivmg)
          wb_step(WB_D) = dtime*wb_gmean(hdmoisg)

        end subroutine atm_wb_transport

        subroutine atm_wb_moisbalance_pre()

          use comphys, only: rmoisg, cormois
          use comatm,  only: iwater

          wb_map_c1(:,:) = cormois(:,:,iwater)
          wb_step(WB_C1) = wb_gmean(wb_map_c1)
          wb_qpre(:,:)  = rmoisg(:,:,iwater)

        end subroutine atm_wb_moisbalance_pre

        subroutine atm_wb_moisbalance_post()

          use comphys, only: rmoisg
          use comatm,  only: iwater

          wb_map_c2(:,:) = wb_qpre(:,:) - rmoisg(:,:,iwater)
          wb_step(WB_C2) = wb_gmean(wb_map_c2)

        end subroutine atm_wb_moisbalance_post

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Internal helpers
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        function wb_gmean(f) result(gm)

          real(dblp), intent(in) :: f(:,:)   !< (nlat,nlon) field
          real(dblp)             :: gm

          integer(ip) :: i

          gm = 0.0_dblp
          do i = 1, nlat
            gm = gm + wb_wgt(i)*sum(f(i,:))
          enddo

        end function wb_gmean

        subroutine wb_set_user_field()

          use comphys, only: dumt1

          integer(ip) :: k

          do k = 1, 3
            dumt1(:,:,k) = sum(wb_ring(:,:,k,1:wb_nring), dim=3)/real(wb_nring, dblp)
          enddo

        end subroutine wb_set_user_field

        subroutine wb_write_month()

          use comemic_mod, only: fini, irunlabel

          real(dblp) :: m(WB_NTERM)

          if (.not. wb_open) then
            open(newunit=wb_unit, file='outputdata/globals/wbudget'//fini//'.txt', form='formatted',               &
                 status='replace', action='write')
            write(wb_unit, '(a,i0)') '# ECBilt water budget (atm_wbudget_mod), global means in mm, irunlabel = ',    &
                                    irunlabel
            write(wb_unit, '(a)') '# R = dW - (E - Pls - Pcv + T + D + C1 - C2 + C3); |R| = sum over steps of |R_step|'
            write(wb_unit, '(a)') '# year month nstep    W_end        dW           E            Pls          Pcv'//   &
                                  '          T            D            C1           C2           C3           R'//   &
                                  '            |R|'
            wb_open = .true.
          endif

          m = wb_month*WB_MM
          write(wb_unit, '(i6,i4,i6,13es13.4)') wb_cur_year, wb_cur_month, wb_nstep, wb_w_last*WB_MM, m(WB_DW),       &
                m(WB_E), m(WB_PLS), m(WB_PCV), m(WB_T), m(WB_D), m(WB_C1), m(WB_C2), m(WB_C3), m(WB_R), m(WB_RAB)
          flush(wb_unit)

        end subroutine wb_write_month

      end module atm_wbudget_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
