!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to save mass that has not been used to produce icebergs so that 
! it can be used after restart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if ( WATGRISCONS == 1 )
      SUBROUTINE ec_wrgrisrunECB

!afq, deprecated --       USE runoff_coupl, ONLY: runoflECB

!afq, deprecated -- #include "comunit.h"

!afq, deprecated --      write(iuo+95) runoflECB
      
      return 
         
      END SUBROUTINE ec_wrgrisrunECB
#endif
