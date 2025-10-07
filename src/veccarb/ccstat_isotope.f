!     dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:27 CET 2009
#include "choixcomposantes.h"
!     dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:27 CET 2009

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      INIT of isotopes in vecarb when no restart in place for carb
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Dideir M. Roche, Veronique Mariotti (vm)
!      Date   : 09 octobre 2012
!      Derniere modification : 07 novembre 2012, dmr, vm
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE CCSTAT_ISOTOPE()

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : 
!       Variables de sortie : 
!-----|--1--------2---------3---------4---------5---------6---------7-|

       USE veget_iso, ONLY: b1t14, b1g14, b1t13, b1g13, b2t14, b2g14
     >              , b2t13, b2g13, b3t14, b3g14, b3t13, b3g13, b4t14
     >              , b4g14, b4t13, b4g13

       USE carbone_co2, ONLY: c14rstd, c14dec
       USE C_res_mod, ONLY: c13init, c13atm

#if ( COMATM == 1 )
       use veget_mod
#endif

       IMPLICIT NONE

#if ( COMATM == 0 )
#include "veget.h"
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

      REAL :: c14init = c14rstd

c initialization of isotopes ratios

      b1t14(:,:)=c14init*exp(-c14dec*t1t)*b1t(:,:)
      b1g14(:,:)=c14init*exp(-c14dec*t1g)*b1g(:,:)

      b1t13(:,:)=b1t(:,:)*c13init*c13atm
      b1g13(:,:)=b1g(:,:)*c13init*c13atm

c   stems and roots biomass

      b2t14(:,:)=b2t(:,:)*c14init*exp(-c14dec*t2t)
      b2g14(:,:)=b2g(:,:)*c14init*exp(-c14dec*t2g)

      b2t13(:,:)=b2t(:,:)*c13init*c13atm
      b2g13(:,:)=b2g(:,:)*c13init*c13atm

c   litter

      b3t14(:,:)=b3t(:,:)*c14init*exp(-c14dec*t3t)
      b3g14(:,:)=b3g(:,:)*c14init*exp(-c14dec*t3g)

      b3t13(:,:)=b3t(:,:)*c13init*c13atm
      b3g13(:,:)=b3g(:,:)*c13init*c13atm

c mortmass and soil organic matter

      b4t14(:,:)=b4t(:,:)*c14init*exp(-c14dec*t4t)
      b4g14(:,:)=b4g(:,:)*c14init*exp(-c14dec*t4g)

      b4t13(:,:)=b4t(:,:)*c13init*c13atm
      b4g13(:,:)=b4g(:,:)*c13init*c13atm

      return

      END SUBROUTINE CCSTAT_ISOTOPE
