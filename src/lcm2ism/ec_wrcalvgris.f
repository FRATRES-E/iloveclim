!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to save mass that has not been used to produce icebergs so that
! it can be used after restart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if ( CALVFLUX > 0 )

      SUBROUTINE ec_wrcalvgris

!! START_OF_USE_SECTION

      USE icb_coupl, ONLY: calvCLIO

!! END_OF_USE_SECTION

!! START_OF_INCLUDE_SECTION

#include "para.com"
#include "comunit.h"

!! END_OF_INCLUDE_SECTION

      write(iuo+95) calvCLIO

      return

      END SUBROUTINE ec_wrcalvgris

#endif
