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
!      MODULE: [fvt_kinds_mod]
!
!>     @brief Kinds and constants of the FV transport library (src/fvtrans), so that it depends on no other module.
!
!      DESCRIPTION:
!>     The library is self-contained: it uses only this module and the intrinsic ones, needs no preprocessor and no
!!     include file. The kinds are those of the host model (iLOVECLIM, global_constants_mod: ip = 4, dp = 8 with
!!     gfortran), so arrays of the host's double precision are passed without conversion; the names carry the fvt_
!!     prefix so that a caller can use both modules without a clash.
!!
!!     Contract (public entities):
!!       fvt_ip : integer kind (default 32-bit integer).
!!       fvt_dp : real kind, IEEE double precision.
!!       fvt_pi : pi in fvt_dp.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module fvt_kinds_mod

        implicit none

        private

        public :: fvt_ip, fvt_dp, fvt_pi

        integer,          parameter :: fvt_ip = selected_int_kind(9)            !< integer kind (4 bytes)
        integer,          parameter :: fvt_dp = selected_real_kind(15, 307)     !< real kind (8 bytes, IEEE double)
        real(kind=fvt_dp), parameter :: fvt_pi = 4.0_fvt_dp*atan(1.0_fvt_dp)  !< pi

      end module fvt_kinds_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
