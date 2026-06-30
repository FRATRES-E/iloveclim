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
!      MODULE: c14_prod_mod
!
!>     @author  Veronique Mariotti (vm), refactoring by dmr, clo
!
!>     @brief   Cosmogenic 14C production and atmospheric 14C update.
!
!>     @details Two public routines, active only when KC14 == 1:
!>                c14_prod  : compute the cosmogenic 14C production rate PRODC14 and the
!>                            associated mass / mole steady-state limits.  With KC14P == 1
!>                            the production is interpolated from the prod_C14.dat record
!>                            file (Be10-derived); with KC14P == 0 it is a fixed rate.
!>                c14atm_dp : advance the atmospheric 14C activity C14ATM by one year
!>                            (radioactive decay then cosmogenic production), once per
!>                            year (KENDY == 1).
!>
!>              Refactoring notes (relative to ec14.f):
!>                - both subroutines wrapped in a module (c14_prod_mod); the whole module
!>                  body remains under #if ( KC14 == 1 ) as in the original file.
!>                - kinds via global_constants_mod (dblp, ip); literals carry _dblp.
!>                - input unit via open(newunit=…); diagnostics go to stdout.
!>                - `use carbone_co2` -> `use carbone_co2_mod`.
!>                - production-rate constants kept verbatim (numerically identical).
!>
!>     @date    Original   : 09 October 2012 (vm), last mod. 10 January 2013
!>     @date    Refactored : 2026-06-29 (dmr, clo)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( KC14 == 1 )

      module c14_prod_mod

        use global_constants_mod, only: dblp=>dp, ip, stdout

        implicit none
        private

        public :: c14_prod, c14atm_dp

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: c14_prod
!
!>     @brief   Compute cosmogenic 14C production rate PRODC14 (ppm 14C yr^-1) and the
!>              derived steady-state mass (LIMC14MASSE) and mole (N14LIM) limits.
!
!>     @details With KC14P == 1, PRODC14 is interpolated in time from prod_C14.dat, whose
!>              PC14M column is normalised so preindustrial production equals 1 and whose
!>              TPSC14 column lists years (negative = BP).  With KC14P == 0, a fixed rate
!>              tuned to give DC14atm = 0 is used.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine c14_prod()

        use carbone_co2_mod, only: PRODC14, MPRODC14, LIMC14MASSE,       &
                                   N14LIM, C14DEC
#if ( KC14P == 1 )
        use carbone_co2_mod, only: NC14max, PC14M, KTIME, n, NC14, NYR,  &
                                   NYR0, NYR01, NYRSR, PC14VAR, TPSC14
#endif

#if ( KC14P == 1 )
        integer(ip) :: bel10dat_id
#endif

#if ( KC14P == 1 )
        !  prod_C14.dat: PC14M normalised so preindustrial production = 1;
        !  TPSC14 in years, negative for BP (e.g. 21 ka BP -> -21000).
        open(newunit=bel10dat_id, file='inputdata/prod_C14.dat', status='unknown')
        write(stdout,*) 'reading prod_C14.dat, NC14max = ', NC14max
        do n = 1, NC14max
          read(bel10dat_id,*) TPSC14(n), PC14M(n)
          write(stdout,*) n, TPSC14(n), PC14M(n)
        end do
        NC14 = n - 1
        close(bel10dat_id)

        !  Resolve current calendar year (real time vs fixed).
        if (KTIME == 1) then
          NYR0 = NYRSR + NYR
        else
          NYR0 = NYRSR
        end if
        NYR01 = -NYR0

        !  Linear interpolation of the read production series.
        if (NYR0 <= TPSC14(1)) then
          PC14VAR = PC14M(1)
        else if (NYR0 >= TPSC14(NC14)) then
          PC14VAR = PC14M(NC14)
        else
          do n = 1, NC14 - 1
            if ((NYR0 >= TPSC14(n)) .and. (NYR0 < TPSC14(n+1)))          &
              PC14VAR = PC14M(n) + (PC14M(n+1) - PC14M(n))               &
                      * (NYR0 - TPSC14(n)) / (TPSC14(n+1) - TPSC14(n))
          end do
        end if

        PRODC14 = 2.625e-12_dblp * 40918.0_dblp / 48296.75_dblp * PC14VAR

#elif ( KC14P == 0 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Fixed cosmogenic production (Masarik & Beer 1999 baseline, iLOVECLIM-tuned).
!  The 40918/48296.75 ratio rescales the CLIMBER PI global carbon budget to
!  iLOVECLIM's; the 1.010 factor matches atmospheric 14C in the PI setting.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        PRODC14 = 3.216e-12_dblp * 1.010_dblp   ! ppm 14C yr^-1
#endif

        mPRODC14    = PRODC14 / 0.40_dblp * 1.0e15_dblp   ! gC14 yr^-1
        limC14masse = mPRODC14 / C14DEC                   ! gC14   (steady state)
        n14lim      = (mPRODC14 / 14.0_dblp) / C14DEC     ! mol C14 (steady state)

      end subroutine c14_prod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: c14atm_dp
!
!>     @brief   Advance atmospheric 14C activity by one year: radioactive decay followed
!>              by cosmogenic production.  Acts only on the last day of the year (KENDY==1).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine c14atm_dp()

        use carbone_co2_mod, only: C14DEC, PRODC14, C14ATM
        use mod_sync_time,   only: KENDY

        call c14_prod()

        if (KENDY == 1) then
          write(stdout,*) 'C14ATM before prod', C14ATM, C14DEC
          C14ATM = (1.0_dblp - C14DEC) * C14ATM
          write(stdout,*) 'C14ATM after decay', C14ATM
          C14ATM = C14ATM + PRODC14
          write(stdout,*) 'C14ATM after prod', C14ATM
        end if

      end subroutine c14atm_dp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      end module c14_prod_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#endif /* KC14 == 1 */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
