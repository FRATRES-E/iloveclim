!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2020-2021 Didier M. Roche (a.k.a. dmr)

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

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#include "choixcomposantes.h"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+---


       MODULE MOD_SYNC_TIME

        USE global_constants_mod, ONLY: dblp=>dp, ip, days_year360d, solar_day

        INTEGER(kind=ip) :: KENDY ! MaJ dans sync_timer
        INTEGER(kind=ip), POINTER :: NYR ! Pointeur dans fait_pointer_CC
        INTEGER(kind=ip), POINTER :: NMONTH ! Pointeur dans fait_pointer_CC
        INTEGER(kind=ip), POINTER :: NDAY ! Pointeur dans fait_pointer_CC nb of days
        INTEGER(kind=ip), POINTER :: TSTOC_jour ! Pointeur dans fait_pointer_CC = nb appel par an

        REAL(kind=dblp), PARAMETER :: TDAY=solar_day, TYER=days_year360d*TDAY
        REAL(KIND=dblp)            :: TSTOC ! TSTOC en secondes ... (compatibilite OCYC)

        INTEGER(kind=ip) :: KMON ! MaJ dans sync_timer ! Ajout temporaire pas dans CLIMBER !!
        INTEGER(kind=ip) :: KWEEK !if last day of week KWEEK=1, else 0
        INTEGER(kind=ip), parameter :: joursmois=30, joursan=360
        INTEGER(kind=ip), parameter :: joursemaine=5 ! number of days per week



       CONTAINS
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+---
!      Cette routine sync_timer sert a synchroniser les variables
!       de temps dans la simulation entre ECbilt et OCYCC
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche, Nathaelle Bouttes
!      Date   : 15 fevrier 2010
!      Derniere modification : 08 octobre 2024
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+---

       SUBROUTINE sync_timer(jour,mois)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  :
!       Variables de sortie :
!-----|--1--------2---------3---------4---------5---------6---------7-|

       implicit none

       integer(kind=ip), intent(in) :: jour,mois

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|


       IF (((mois-1)*joursmois + jour).EQ.joursan) THEN
        KENDY=1
       ELSE
        KENDY=0
       ENDIF

       IF (jour.EQ.joursmois) THEN
        KMON=1
       ELSE
        KMON=0
       ENDIF

! We have 5 days per week, and 6 weeks per month and 12 month per year
! (360 days/year)
       IF (MODULO(((mois-1)*joursmois + jour),joursemaine).EQ.0) THEN
        KWEEK=1
       ELSE
        KWEEK=0
       ENDIF

       return
       END SUBROUTINE sync_timer

       subroutine fait_pointer_CC_TIME

       use comemic_mod, only: iyear, nocstpyear, imonth, iday

       IMPLICIT NONE

       NYR => iyear
       NMONTH => imonth
       NDAY => iday
       TSTOC_jour => nocstpyear
       TSTOC = REAL(TSTOC_jour*(TDAY/TYER)*TDAY)

       end subroutine fait_pointer_CC_TIME

       END MODULE MOD_SYNC_TIME

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

