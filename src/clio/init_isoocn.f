!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine <squelette> sert a ???
!       et caetera, et caetera
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche 
!      Date   : 13 octobre 2011
!      Derniere modification : 25 septembre 2018
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE init_isoocn

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : ~
!       Variables de sortie : scal ocean tracers, of nsmax size
!-----|--1--------2---------3---------4---------5---------6---------7-|

!! START_OF_USE_SECTION

#if ( ISOOCN >= 1 )

      use iso_param_mod, only : rsmow 
      use para0_mod,     only : isoocn_restart, ocnw17, ocnw18, ocnw2h, 
     &                          owisostrt, owisostop, oc2atwisoindx
      use bloc0_mod,     only : scal ! ocean tracers, of nsmax size

!! END_OF_USE_SECTION

!! START_OF_INCLUDE_SECTION

! [TRNOTYPE] #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"

!! END_OF_INCLUDE_SECTION


!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       integer :: iz

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      If we do not perform a restart ...
!-----|--1--------2---------3---------4---------5---------6---------7-|
       if ( isoocn_restart.EQ.0) then


         
#if ( LGMSWITCH == 1 )

! dmr  --- INIT SECTION WHEN NO ISO restart AT LGM
! dmr  ---   choosing meand18O = +1.0E-3, no excess ...

         call oc2atwisoindx(ocnw18,iz)
         scal(:,:,:,ocnw18) = (1.0+1.0E-3)*rsmow(iz)

! dmr  ---
!       Then we assume a d-excess of 0 and a d17excess of zero
!     
!       d-excess is defined as : d-excess = dD - 8*d18O
!         choosing d-excess = 0.0 implies dD = 8 * d18O
!         so : (Rd / Rdsmow -1.0)*1000. = d18O * 8.0
!         meaning : Rd = ( d18O/1000.* 8.0 + 1.0 ) * Rdsmow
! dmr   ---

         call oc2atwisoindx(ocnw2h,iz)
         scal(:,:,:,ocnw2h) = ( 1.0E-3 * 8.0 + 1.0 ) * rsmow(iz)


! dmr   --- 
!
!       17Oexcess is defined as : 
!           17Oexcess = ln(d17O/1000.+1.0) - 0.528*ln(d18O/1000.+1.0)
!         choosing 17Oexcess = 0.0 thus implies : 
!            R17O = ((d18O/1000.+1.0) **0.528) * R17smow
! dmr   ---

         call oc2atwisoindx(ocnw17,iz)
         scal(:,:,:,ocnw17) = ( (1.0E-3 + 1.0) **0.528 ) * rsmow(iz)


#else

!dmr --- INIT SECTION WHEN NO ISO restart AT PD
!dmr ---  i.e. all water isotopes at rsmow ...

         call oc2atwisoindx(ocnw17,iz)
         scal(:,:,:,ocnw17) = rsmow(iz)

         call oc2atwisoindx(ocnw18,iz)
         scal(:,:,:,ocnw18) = rsmow(iz)

         call oc2atwisoindx(ocnw2h,iz)
         scal(:,:,:,ocnw2h) = rsmow(iz)

#endif

       else ! isoocn_restart != 0
       
!-----|--1--------2---------3---------4---------5---------6---------7-|
!      If we do perform a restart ...
!-----|--1--------2---------3---------4---------5---------6---------7-|
         OPEN(unit=20,FILE='startdata/wisoocn_restart.dat',STATUS='old'
     &,FORM='unformatted')
           READ(20) scal(:,:,:,owisostrt:owisostop)
         CLOSE(20)
       endif

#endif
       END SUBROUTINE init_isoocn
