!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:54 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:54 CET 2009

      FUNCTION sol_cfc(pt,ps,kn)
!_ ---------------------------------------------------------------------
!_ RCS lines preceded by "c_ "
!_ ---------------------------------------------------------------------
!_
!_ $Source: /usr/people/severijn/.cvs-repository/EMIC/clio/sources/sol_cfc.f,v $
!_ $Revision: 1.1.1.1 $
!_ $Date: 2001/05/21 14:42:20 $   ;  $State: Exp $
!_ $Author: severijn $ ;  $Locker:  $
!_
!_ ---------------------------------------------------------------------
!_ $Log: sol_cfc.f,v $
!_ Revision 1.1.1.1  2001/05/21 14:42:20  severijn
!_ Import of ECBILT/CLIO
!_
!_ Revision 1.2  1998/07/17 07:37:02  jomce
!_ Fixed slight bug in units conversion: converted 1.0*e-12 to 1.0e-12
!_ following warning from Matthew Hecht at NCAR.
!_
!_ Revision 1.1  1998/07/07 15:22:00  orr
!_ Initial revision
!_
!_ ---------------------------------------------------------------------
!_
!jmc  REAL FUNCTION sol_cfc(pt,ps,kn)
!-------------------------------------------------------------------
!
!     CFC 11 and 12 Solubilities in seawater
!     ref: Warner & Weiss (1985) , Deep Sea Research, vol32
!
!     pt:       temperature (degre Celcius)
!     ps:       salinity    (o/oo)
!     kn:       11 = CFC-11, 12 = CFC-12
!     sol_cfc:  in mol/m3/pptv
!               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm

!
!     J-C Dutay - LSCE
!-------------------------------------------------------------------

      REAL    pt, ps,ta,d,sol_cfc
      REAL a1 ( 11: 12), a2 ( 11: 12), a3 ( 11: 12), a4 ( 11: 12)
      REAL b1 ( 11: 12), b2 ( 11: 12), b3 ( 11: 12)

      INTEGER kn

!c
!c coefficient for solubility in  mol/l/atm
!c ----------------------------------------
!
!     for CFC 11
!     ----------
      a1 ( 11) = -229.9261
      a2 ( 11) =  319.6552
      a3 ( 11) =  119.4471
      a4 ( 11) =   -1.39165
      b1 ( 11) =   -0.142382
      b2 ( 11) =    0.091459
      b3 ( 11) =   -0.0157274
!    
!     for CFC/12
!     ----------
      a1 ( 12) = -218.0971
      a2 ( 12) =  298.9702
      a3 ( 12) =  113.8049
      a4 ( 12) =   -1.39165
      b1 ( 12) =   -0.143566
      b2 ( 12) =    0.091015
      b3 ( 12) =   -0.0153924
!

      ta       = ( pt + 273.16)* 0.01
      d    = ( b3 ( kn)* ta + b2 ( kn))* ta + b1 ( kn)
!
!
      sol_cfc 
     $    = exp ( a1 ( kn)
     $    +       a2 ( kn)/ ta 
     $    +       a3 ( kn)* alog ( ta )
     $    +       a4 ( kn)* ta * ta  + ps* d )
!
!     conversion from mol/(l * atm) to mol/(m^3 * atm) 
!     ------------------------------------------------
      sol_cfc = 1000. * sol_cfc
!
!     conversion from mol/(m^3 * atm) to mol/(m3 * pptv) 
!     --------------------------------------------------
      sol_cfc = 1.0e-12 * sol_cfc

      END 
