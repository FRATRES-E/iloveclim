!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2026 iLOVECLIM / LUDUS coding group

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: restart_cc_mod
!
!>     @author  Nathaelle Bouttes (nb), refactoring by dmr, clo
!
!>     @brief   Read / write the global carbon-reservoir restart file (rest_cc.dat).
!
!>     @details Holds a single public routine, restart_cc(choix), that serialises the
!>              run-level carbon scalars (atmospheric CO2, 13C, 14C, ocean/land totals,
!>              alkalinity) to a direct-access binary file.  choix == 0 writes, choix == 1
!>              reads; on read it also restores the derived initial values (PA0_C, c13atm,
!>              the "before" 14C reservoirs) and rescales ocean alkalinity by the volume
!>              ratio when BATHY >= 1.
!>
!>              Refactoring notes (relative to restart_cc.f90):
!>                - wrapped the free subroutine in a module (restart_cc_mod).
!>                - file unit obtained with open(newunit=…) instead of the hand-rolled
!>                  INQUIRE search loop.
!>                - kinds via global_constants_mod (dblp, ip); the record length is sized
!>                  from kind(PA0_C) exactly as before.
!>                - `use carbone_co2` -> `use carbone_co2_mod`.
!>                - character lengths written as character(len=*) parameters.
!>                - I/O logic and record layout are byte-for-byte identical to the original
!>                  so existing rest_cc.dat files remain readable.
!>
!>     @date    Original   : 26 February 2019 (nb)
!>     @date    Refactored : 2026-06-29 (dmr, clo)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module restart_cc_mod

        use global_constants_mod, only: dblp=>dp, ip, stdout

        implicit none
        private

        public :: restart_cc

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: restart_cc
!
!>     @brief   Write (choix==0) or read (choix==1) the carbon-reservoir restart.
!>     @param[in] choix   0 = write rest_cc.dat, 1 = read startdata/rest_cc.dat
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine restart_cc(choix)

        use C_res_mod, only: ca_oc_rest, ca_la_rest, cav_oc, cav_la,     &
                             ca13_oc_rest, ca13_la_rest, cav_oc13,       &
                             cav_la13, c13atm, c13atm_rest, alk_oc_rest
#if ( KC14 == 1 )
        use C_res_mod, only: cav_oc14_rest, cav_la14_rest, cav_oc14,     &
                             cav_la14, cav_oc14_b, cav_la14_b,           &
                             cav_oc14_b_rest, cav_la14_b_rest,           &
                             cav_oc_b_rest, cav_la_b_rest,               &
                             cav_oc_b, cav_la_b
#endif

        use carbone_co2_mod, only: PA0_C, PA_C, C14ATM, C14ATM_rest

        use marine_bio_mod,  only: OALK_ini
#if ( BATHY >= 1 )
        use loveclim_transfer_mod, only: OVOL, OVOL_prev
#endif

        integer(ip), intent(in) :: choix

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Local variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip) :: fich_num, nrecl

        character(len=*), parameter :: fich_res_name     = 'rest_cc.dat'
        character(len=*), parameter :: fich_res_name_old = 'startdata/rest_cc.dat'

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Record length: number of scalars written x storage size of one (kind of PA0_C).
!  7 variables in the base case, 14 when 14C is active.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( KC14 == 1 )
        nrecl = 14 * kind(PA0_C)
#else
        nrecl = 7  * kind(PA0_C)
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  choix == 0 : write the restart
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (choix == 0) then

          ca_oc_rest   = cav_oc
          ca_la_rest   = cav_la
          ca13_oc_rest = cav_oc13
          ca13_la_rest = cav_la13
          c13atm_rest  = c13atm
          alk_oc_rest  = OALK_ini
#if ( KC14 == 1 )
          cav_oc14_rest   = cav_oc14
          cav_la14_rest   = cav_la14
          C14ATM_rest     = C14ATM
          cav_oc14_b_rest = cav_oc14_b
          cav_la14_b_rest = cav_la14_b
          cav_oc_b_rest   = cav_oc_b
          cav_la_b_rest   = cav_la_b
#endif

          open(newunit=fich_num, file=fich_res_name, status='unknown',  &
               access='direct', recl=nrecl, action='write')
#if ( KC14 == 1 )
          write(unit=fich_num, rec=1)                                   &
                ca_oc_rest, ca_la_rest, PA_C,                           &
                ca13_oc_rest, ca13_la_rest, c13atm_rest,                &
                alk_oc_rest,                                            &
                cav_oc14_rest, cav_la14_rest, C14ATM_rest,              &
                cav_oc14_b_rest, cav_la14_b_rest,                       &
                cav_oc_b_rest, cav_la_b_rest
#else
          write(unit=fich_num, rec=1)                                   &
                ca_oc_rest, ca_la_rest, PA_C,                           &
                ca13_oc_rest, ca13_la_rest, c13atm_rest,                &
                alk_oc_rest
#endif
          close(unit=fich_num)

          write(stdout,*) 'write carbon values in rest_cc.dat'
          write(stdout,*) ca_oc_rest, ca_la_rest, ca13_oc_rest, ca13_la_rest
          write(stdout,*) PA_C, c13atm_rest, alk_oc_rest

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  choix == 1 : read the restart
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        else if (choix == 1) then

          write(stdout,*) 'Reading carbon restart from file:'
          write(stdout,*) fich_res_name_old

          open(newunit=fich_num, file=fich_res_name_old, status='old',  &
               access='direct', recl=nrecl, action='read')
#if ( KC14 == 1 )
          read(unit=fich_num, rec=1)                                    &
                ca_oc_rest, ca_la_rest, PA_C,                           &
                ca13_oc_rest, ca13_la_rest, c13atm_rest,                &
                alk_oc_rest,                                            &
                cav_oc14_rest, cav_la14_rest, C14ATM_rest,              &
                cav_oc14_b_rest, cav_la14_b_rest,                       &
                cav_oc_b_rest, cav_la_b_rest
#else
          read(unit=fich_num, rec=1)                                    &
                ca_oc_rest, ca_la_rest, PA_C,                           &
                ca13_oc_rest, ca13_la_rest, c13atm_rest,                &
                alk_oc_rest
#endif
          close(unit=fich_num)

          PA0_C  = PA_C
          c13atm = c13atm_rest
#if ( KC14 == 1 )
          C14ATM     = C14ATM_rest
          cav_oc14   = cav_oc14_rest
          cav_la14   = cav_la14_rest
          cav_oc14_b = cav_oc14_b_rest
          cav_la14_b = cav_la14_b_rest
          cav_oc_b   = cav_oc_b_rest
          cav_la_b   = cav_la_b_rest
#endif
          write(stdout,*) 'initialisation of atmospheric carbon'
          write(stdout,*) PA0_C, c13atm, C14ATM

          !  Restore ocean alkalinity; rescale by volume ratio under variable bathymetry.
          OALK_ini = alk_oc_rest
#if ( BATHY >= 1 )
          OALK_ini = alk_oc_rest * OVOL_prev / OVOL
#endif

        end if

      end subroutine restart_cc

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      end module restart_cc_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
