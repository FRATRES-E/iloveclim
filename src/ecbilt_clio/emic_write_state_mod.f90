!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of iLOVECLIM/COUPLER
!!      iLOVECLIM/COUPLER is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
!!      as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!!
!!      Foobar is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!!      of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!!
!!      You should have received a copy of the GNU General Public License along with Foobar.
!!      If not, see <http://www.gnu.org/licenses/>.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"
!dmr -- Adding the choice of components through the pre-processing options

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: emic_write_state_mod
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module emic_write_state_mod is handling the coherent restart state for all subcomponents of ILOVECLIM
!
!>     @date Creation date: October, 07th, 2019
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     On the long run, each piece/component should be moved to a different subroutine, not in one big subroutine as now
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module emic_write_state_mod


       implicit none
       private

       public :: ec_writestate

      ! NOTE_avoid_public_variables_if_possible

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ec_writestate
!
!>     @brief This subroutine / function is adding a wonderful knew unused functionality
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ec_writestate(ist,ndaytot) result(returnValue)

#define IFORT_USAGE 0
#define PGIFORTRAN 1

      USE comatm,      only:
      use comemic_mod, only: iatm, nstpyear, iyear, nwrskip, nyears, irunlabel
      use comunit,     only: iuo

#if ( OCYCC == 1 && CLIO_OUT_NEWGEN == 1 )
use WRTE_RESTART_IO_NC, only: WRTE_RESTART_OCYCC
#endif

#if ( MEDUSA == 1 )
      USE flux_from_sediments_mod, only : restart_flx_from_sediments
#endif

!dmr --- [DEPRECATED]#if ( VAMPER_KEY >= 1 )
!dmr --- [DEPRECATED]      USE Params, ONLY: V_RESTART_FILE, V2_RESULTS_FILE, V_RESULTS_FILE
!dmr --- [DEPRECATED]      USE Error
!dmr --- [DEPRECATED]#endif

#if ( IFORT_USAGE == 1 )

!     dmr --- Ajout pour la creation de redemarrages en repertoires ...
      USE IFPORT, ONLY : SYSTEMQQ, MAKEDIRQQ, IERRNO,
     &     RENAMEFILEQQ, ENOENT
!     dmr ---

#endif

#if ( PATH >= 1 )
      use path_mod, only: restart_fnm, path_write_rest
#endif

#if ( CORAL == 1 )
      use coral_mod, only: restart_coral
#endif

      use newunit_mod, only: newunit_id, wisocpl_restart_id
      use landmodel_mod, only: ec_wrendland

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  ist
!>    @param[in]  ndaytot
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical                :: returnValue

       CHARACTER(len=255) :: cwd
       integer istep,ist, ndaytot
       character*6 chf

!     dmr ---
       CHARACTER*300 ::  REPERTOIRE, COMMANDE, file_move, file_dest
       INTEGER (KIND=4) :: ETAT
       LOGICAL :: RESULT, existant

#if ( IFORT_USAGE == 0 )
       integer :: resultat
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



!dmr --- set timing ...

      istep=ist*iatm

      if (mod(istep,nstpyear).eq.0) then
         if (mod(iyear,nwrskip).eq.0.or.iyear.eq.nyears) then

            write(chf,1) iyear+irunlabel
 1          format(i6.6)


!dmr --- Location where we are putting this particular restart state
            REPERTOIRE=''//'restartdata/'//'res'//chf

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Let's create that directory if it does not exist
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr --- intel fortran version
#if ( IFORT_USAGE == 1 )

            INQUIRE(DIRECTORY=TRIM(REPERTOIRE), EXIST=existant)

            IF (.NOT. existant) THEN
               RESULT = MAKEDIRQQ(TRIM(REPERTOIRE))
               WRITE(*,*) 'Created ', TRIM(REPERTOIRE), RESULT
            ENDIF

#else

!dmr --- pgi fortran version
#if ( PGIFORTRAN == 1 )

                call system('test -d '//TRIM(REPERTOIRE))
                call system('mkdir '//TRIM(REPERTOIRE))

!dmr --- gnu fortran version
#else
            call EXECUTE_COMMAND_LINE('test -d '//TRIM(REPERTOIRE), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
                call EXECUTE_COMMAND_LINE('mkdir '//TRIM(REPERTOIRE),wait=.TRUE., exitstat=resultat)
                if (resultat .ne. 0) then
                  WRITE(*,*) 'Problem in directory creation ', TRIM(REPERTOIRE), resultat
                else
                  WRITE(*,*) 'Created ', TRIM(REPERTOIRE), resultat
                endif
            endif
#endif /* on PGI FORTRAN */
#endif /* on IFORT_USAGE */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Write section: ATMOSPHERE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr --- Dynamics

            open(newunit=newunit_id,file=''//TRIM(REPERTOIRE)//'/inatdyn'//chf//'.dat',form='unformatted')
            call ec_wrenddyn
            close(newunit_id)

!dmr --- [DELETE] I do not think we need to call this twice?
!BdB-02-2019 - I think you need this for restarting, it is inatphy
             open(newunit=newunit_id,file=''//TRIM(REPERTOIRE)//'/inatphy'//chf//'.dat',form='unformatted')
             call ec_wrendphy
             close(newunit_id)


!dmr --- Water isotopes, just need to copy a file

#if ( ISOATM >= 1 )
          COMMANDE='cp '//'wisoatm_restart.dat '//TRIM(REPERTOIRE)//'/.'

!dmr --- intel fortran version
#if ( IFORT_USAGE == 1 )

            WRITE(*,*) "test ", TRIM(COMMANDE)
            ETAT = SYSTEMQQ(TRIM(COMMANDE))
            WRITE(*,*) 'Resultat du cp : ', ETAT
#else

!dmr --- pgi fortran version
#if ( PGIFORTRAN == 1 )
            call system(TRIM(COMMANDE))

!dmr --- gnu fortran version
#else
            call EXECUTE_COMMAND_LINE(TRIM(COMMANDE), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du cp : ', resultat
            endif
#endif /* on PGI FORTRAN */
#endif /* on IFORT_USAGE */

#endif /* on ISOATM */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Write section: LAND MODEL LBM
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr --- Land surface variables

            open(newunit=newunit_id,file=''//TRIM(REPERTOIRE)//'/inland'//chf//'.dat',form='unformatted')
            call ec_wrendland
            close(newunit_id)

!dmr --- water isotopes on land, just need to copy a file

#if ( ISOLBM >= 1 )
          COMMANDE='cp '//'wisolbm_restart.dat '//TRIM(REPERTOIRE)//'/.'

!dmr --- intel fortran version
#if ( IFORT_USAGE == 1 )
            WRITE(*,*) "test ", TRIM(COMMANDE)
            ETAT = SYSTEMQQ(TRIM(COMMANDE))
            WRITE(*,*) 'Resultat du cp : ', ETAT
#else

!dmr --- pgi fortran version
#if ( PGIFORTRAN == 1 )
            call system(TRIM(COMMANDE))
#else

!dmr --- gnu fortran version
            call EXECUTE_COMMAND_LINE(TRIM(COMMANDE), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du cp : ', resultat
            endif
#endif /* on PGI FORTRAN */
#endif /* on IFORT_USAGE */

#endif /* ISOLBM */


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Write section: COUPLER
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr--- For coupler variables
            open(newunit=newunit_id ,file=''//TRIM(REPERTOIRE)//'/incoup'//chf//'.dat',form='unformatted')

!dmr--- For water isotopes variables
#if ( ISOATM >= 2 && ISOOCN >= 1 )
            open(newunit=wisocpl_restart_id,file=''//TRIM(REPERTOIRE)//'/wisocpl_restart.dat', form='unformatted')
#endif

            call ec_wrendcoup

            close(newunit_id)
#if ( ISOATM >= 2 && ISOOCN >= 1 )
            close(wisocpl_restart_id)
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Write section: ICEBERGS
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!PB#if (CALVFLUX == 1 )
!PB!mab: coupling GRISLI - ICEBERGS
!PB         open(iuo+95,file=''//TRIM(REPERTOIRE)//'/mass_icb'//chf//'.dat',form='unformatted')
!PB          call ec_wrendcalv
!PB          close(iuo+95)
!PB#endif

!dmr --- intel fortran version

#if ( IFORT_USAGE == 1 )
!dmr [DEPRECATED]          if (flgicb) then
#if ( ICEBERG_MODEL > 0 )

            file_move='resicb.om'
            file_dest=''//TRIM(REPERTOIRE)//'/'//trim(file_move)

            INQUIRE(FILE=TRIM(file_move), EXIST=existant)

            IF (existant) THEN
               RESULT = RENAMEFILEQQ(TRIM(file_move),TRIM(file_dest))
               WRITE(*,*) 'Moved ', TRIM(file_move), RESULT
            ENDIF

!PB            file_move='resvol.om'
!PB            file_dest=''//TRIM(REPERTOIRE)//'/'//trim(file_move)

!PB            INQUIRE(FILE=TRIM(file_move), EXIST=existant)

!PB            IF (existant) THEN
!PB               RESULT = RENAMEFILEQQ(TRIM(file_move),TRIM(file_dest))
!PB               WRITE(*,*) 'Moved ', TRIM(file_move), RESULT
!PB            ENDIF

!dmr [DEPRECATED]         endif

#endif
!dmr --- pgi & gnu fortran version

#else


!dmr --- pgi fortran version

#if ( PGIFORTRAN == 1 )
!dmr [DEPRECATED]          if (flgicb) then
#if ( ICEBERG_MODEL > 0 )
            file_move='resicb.om'
            file_dest=''//TRIM(REPERTOIRE)//'/'//trim(file_move)

            call system("mv "//TRIM(file_move)//" "//TRIM(file_dest))

            file_move='resvol.om'
            file_dest=''//TRIM(REPERTOIRE)//'/'//trim(file_move)

            call system("mv "//TRIM(file_move)//" "//TRIM(file_dest))

!dmr [DEPRECATED]          endif
#endif
#else

#if ( ICEBERG_MODEL > 0 )
            write(*,*) "Attention partie de code non implentée !!"
            read(*,*)
#endif

#endif
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Write section: ICE SHEET CALVING [DEPRECATED?]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!PB#if (CALVFLUX >= 1 )
!PB          open(iuo+95,file=''//TRIM(REPERTOIRE)//'/mass_calvgris'//chf//'.dat',form='unformatted')
!PB          call ec_wrcalvgris
!PB          close(iuo+95)

!PB#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Write section: ICE-SHEET / ICEBERG COUPLING [DEPRECATED?]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!afq, deprecated -- #if (WATGRISCONS == 1 )
!afq, deprecated -- !PB          open(newunit=newunit_id
!afq, deprecated --          open(iuo+95,file=''//TRIM(REPERTOIRE)//'/mass_grisrunECB'//chf//'.dat',form='unformatted')
!afq, deprecated --           call ec_wrgrisrunECB
!afq, deprecated --           close(iuo+95)
!afq, deprecated -- #endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Write section: CLIO OCEAN MODEL
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr physical ocean variables: just need to move a file ...

            file_move='res'//chf//'.om'
            file_dest=''//TRIM(REPERTOIRE)//'/'//trim(file_move)

!dmr --- intel fortran version

#if ( IFORT_USAGE == 1 )

            INQUIRE(FILE=TRIM(file_move), EXIST=existant)

            IF (existant) THEN
               RESULT = RENAMEFILEQQ(TRIM(file_move),TRIM(file_dest))
               WRITE(*,*) 'Moved ', TRIM(file_move), RESULT
            ELSE
               WRITE(*,*) "I CANNOT MOVE FILE:", TRIM(file_move)
               WRITE(*,*) "WHAT THE HELL IS GOING ON?"
               WRITE(*,*) "I STOP DOING THIS STUPID BUSINESS ...."
               STOP
            ENDIF

#else

!dmr --- pgi fortran version

#if ( PGIFORTRAN == 1 )

          call system("mv "//TRIM(file_move)//" "//TRIM(file_dest))

!dmr --- gnu fortran version

#else

          call EXECUTE_COMMAND_LINE("mv "//TRIM(file_move)//" "//TRIM(file_dest), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du mv : ', resultat
            endif

#endif /* on PGI FORTRAN */
#endif /* on IFORT_USAGE */

!#if ( BATHY >= 1 )
!nb by default
!dmr Save the masks of current run for evolving bathy: just need to move a file ...
            file_move='res_masks'//chf//'.om'
            file_dest=''//TRIM(REPERTOIRE)//'/'//trim(file_move)

!dmr --- intel fortran version

#if ( IFORT_USAGE == 1 )

            INQUIRE(FILE=TRIM(file_move), EXIST=existant)

            IF (existant) THEN
               RESULT = RENAMEFILEQQ(TRIM(file_move),TRIM(file_dest))
               WRITE(*,*) 'Moved ', TRIM(file_move), RESULT
            ELSE
               WRITE(*,*) "I CANNOT MOVE FILE:", TRIM(file_move)
               WRITE(*,*) "WHAT THE HELL IS GOING ON?"
               WRITE(*,*) "I STOP DOING THIS STUPID BUSINESS ...."
               STOP
            ENDIF

#else

!dmr --- pgi fortran version

#if ( PGIFORTRAN == 1 )

          call system("mv "//TRIM(file_move)//" "//TRIM(file_dest))

!dmr --- gnu fortran version

#else

          call EXECUTE_COMMAND_LINE("mv "//TRIM(file_move)//" "//TRIM(file_dest), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du mv : ', resultat
            endif

#endif /* on PGI FORTRAN */
#endif /* on IFORT_USAGE */
!#endif /* on BATHY */
!dmr --- file where the water conservation inbalance for CLIO is reported -- just a copy

            COMMANDE='cp '//'correcw.dat '//TRIM(REPERTOIRE)//'/.'

!dmr --- intel fortran version

#if ( IFORT_USAGE == 1 )
            WRITE(*,*) "test ", TRIM(COMMANDE)
            ETAT = SYSTEMQQ(TRIM(COMMANDE))
            WRITE(*,*) 'Resultat du cp : ', ETAT
#else

!dmr --- pgi fortran version

#if ( PGIFORTRAN == 1 )

            call system(TRIM(COMMANDE))

#else

!dmr --- gnu fortran version

            call EXECUTE_COMMAND_LINE(TRIM(COMMANDE), wait=.TRUE.
     &          , exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du cp : ', resultat
            endif

#endif /* on PGI FORTRAN */
#endif /* on IFORT_USAGE */

!dmr water isotopes in the CLIO ocean if any .... just need to copy the restart file

#if ( ISOOCN >= 1 )
          COMMANDE='cp '//'wisoocn_restart.dat '//TRIM(REPERTOIRE)//'/.'

!dmr --- intel fortran version

#if ( IFORT_USAGE == 1 )
            WRITE(*,*) "test ", TRIM(COMMANDE)
            ETAT = SYSTEMQQ(TRIM(COMMANDE))
            WRITE(*,*) 'Resultat du cp : ', ETAT
#else

!dmr --- pgi fortran version

#if ( PGIFORTRAN == 1 )

            call system(TRIM(COMMANDE))

#else

!dmr --- gnu fortran version

            call EXECUTE_COMMAND_LINE(TRIM(COMMANDE), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du cp : ', resultat
            endif

#endif /* on PGI FORTRAN */
#endif /* on IFORT_USAGE */

#endif /* on ISOOCN */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Write section: VECODE VEGETATION MODEL
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!     MB  --- veget.init/.rest is created in veget.f

            file_move='veget.init'
            file_dest=''//TRIM(REPERTOIRE)//'/'//TRIM(file_move)

#if ( IFORT_USAGE == 1 )
            INQUIRE(FILE=TRIM(file_move), EXIST=existant)

            IF (existant) THEN
               RESULT = RENAMEFILEQQ(TRIM(file_move),TRIM(file_dest))
               WRITE(*,*) 'Moved ', TRIM(file_move), RESULT
            ENDIF
#else
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!

#if ( PGIFORTRAN == 1 )

          call system("mv "//TRIM(file_move)//" "//TRIM(file_dest))

#else

          call EXECUTE_COMMAND_LINE("mv "//TRIM(file_move)//" "//TRIM(file_dest), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du mv : ', resultat
            endif

#endif
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#endif

            file_move='veget.rest'
            file_dest=''//TRIM(REPERTOIRE)//'/'//TRIM(file_move)

#if ( IFORT_USAGE == 1 )

            INQUIRE(FILE=TRIM(file_move), EXIST=existant)

            IF (existant) THEN
               RESULT = RENAMEFILEQQ(TRIM(file_move),TRIM(file_dest))
               WRITE(*,*) 'Moved ', TRIM(file_move), RESULT
            ENDIF
#else
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
#if ( PGIFORTRAN == 1 )

          call system("mv "//TRIM(file_move)//" "//TRIM(file_dest))

#else

          call EXECUTE_COMMAND_LINE("mv "//TRIM(file_move)//" "//TRIM(file_dest), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du mv : ', resultat
            endif

#endif
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#endif

!cnb write restart for total carbon in each reservoir
#if ( CYCC == 2 )
            CALL  restart_cc(0)

            file_move='rest_cc.dat'
            file_dest=''//TRIM(REPERTOIRE)//'/'//TRIM(file_move)

#if ( IFORT_USAGE == 1 )

            INQUIRE(FILE=TRIM(file_move), EXIST=existant)

            IF (existant) THEN
               RESULT = RENAMEFILEQQ(TRIM(file_move),TRIM(file_dest))
               WRITE(*,*) 'Moved ', TRIM(file_move), RESULT
            ENDIF
#else

#if ( PGIFORTRAN == 1 )

          call system("mv "//TRIM(file_move)//" "//TRIM(file_dest))

#else
          call EXECUTE_COMMAND_LINE("mv "//TRIM(file_move)//" "//TRIM(file_dest), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du mv : ', resultat
            endif
#endif

#endif
#endif
!cnb end restart CC

!cnb write restart for corals
#if ( CORAL == 1 )
            CALL  restart_coral(0)

            file_move='rest_coral.dat'
            file_dest=''//TRIM(REPERTOIRE)//'/'//TRIM(file_move)

#if ( IFORT_USAGE == 1 )

            INQUIRE(FILE=TRIM(file_move), EXIST=existant)

            IF (existant) THEN
               RESULT = RENAMEFILEQQ(TRIM(file_move),TRIM(file_dest))
               WRITE(*,*) 'Moved ', TRIM(file_move), RESULT
            ENDIF
#else

#if ( PGIFORTRAN == 1 )

          call system("mv "//TRIM(file_move)//" "//TRIM(file_dest))

#else
          call EXECUTE_COMMAND_LINE("mv "//TRIM(file_move)//"
"//TRIM(file_dest), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du mv : ', resultat
            endif
#endif

#endif
#endif
!cnb end restart corals


#if ( OCYCC == 1 )
            CALL  restart_mb(0)

            file_move='rest_mb.dat'
            file_dest=''//TRIM(REPERTOIRE)//'/'//TRIM(file_move)

#if ( CLIO_OUT_NEWGEN == 1 )
            CALL WRTE_RESTART_OCYCC()
#endif

!tbd cnb write restart for total carbon in each reservoir
!tbd             CALL  restart_cc(0)


#if ( IFORT_USAGE == 1 )

    INQUIRE(FILE=TRIM(file_move), EXIST=existant)

    IF (existant) THEN
       RESULT = RENAMEFILEQQ(TRIM(file_move),TRIM(file_dest))
       WRITE(*,*) 'Moved ', TRIM(file_move), RESULT
    ENDIF
#else
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( PGIFORTRAN == 1 )

  call system("mv "//TRIM(file_move)//" "//TRIM(file_dest))

#else
  call EXECUTE_COMMAND_LINE("mv "//TRIM(file_move)//" "//TRIM(file_dest), wait=.TRUE., exitstat=resultat)
    if ( resultat .ne. 0 ) then
      WRITE(*,*) 'Resultat du mv : ', resultat
    endif
#endif
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#endif
#endif

!cnb write restart sediments
#if ( MEDUSA == 1 )
    CALL restart_flx_from_sediments(0)

    file_move='rest_sed.dat'
    file_dest=''//TRIM(REPERTOIRE)//'/'//TRIM(file_move)

#if ( IFORT_USAGE == 1 )

            INQUIRE(FILE=TRIM(file_move), EXIST=existant)

            IF (existant) THEN
               RESULT = RENAMEFILEQQ(TRIM(file_move),TRIM(file_dest))
               WRITE(*,*) 'Moved ', TRIM(file_move), RESULT
            ENDIF
#else
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( PGIFORTRAN == 1 )

          call system("mv "//TRIM(file_move)//" "//TRIM(file_dest))

#else
          call EXECUTE_COMMAND_LINE("mv "//TRIM(file_move)//" "//TRIM(file_dest), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du mv : ', resultat
            endif
#endif
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#endif
#endif
!cnb end sediments

#if ( ISM == 2 )
            COMMANDE='mv -f '//'*-grestart* '//TRIM(REPERTOIRE)//'/.'
#if ( IFORT_USAGE == 1 )
            ETAT = SYSTEMQQ(TRIM(COMMANDE))
            WRITE(*,*) 'Moved ISM: ', ETAT
#else
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!

#if ( PGIFORTRAN == 1 )

            call system(TRIM(COMMANDE))

#else

            call EXECUTE_COMMAND_LINE(TRIM(COMMANDE), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du mv : ', resultat
            endif

#endif
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#endif
#endif

#if ( PATH >= 1 )
          call path_write_rest()

          COMMANDE='cp '//trim(restart_fnm)//' '//TRIM(REPERTOIRE)//'/.'

#if ( IFORT_USAGE == 1 )
            WRITE(*,*) "test ", TRIM(COMMANDE)
            ETAT = SYSTEMQQ(TRIM(COMMANDE))
            WRITE(*,*) 'Resultat du cp : ', ETAT
#else
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
#if ( PGIFORTRAN == 1 )

            call system(TRIM(COMMANDE))

#else
            call EXECUTE_COMMAND_LINE(TRIM(COMMANDE), wait=.TRUE., exitstat=resultat)
            if ( resultat .ne. 0 ) then
              WRITE(*,*) 'Resultat du cp : ', resultat
            endif

#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#endif
#endif

         endif
      endif


      returnValue = .TRUE.
      return

      end function ec_writestate

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module emic_write_state_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
