!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009

      FUNCTION sc_cfc(t,kn)
!_ ---------------------------------------------------------------------
!_ RCS lines preceded by "c_ "
!_ ---------------------------------------------------------------------
!_
!_ $Source: /usr/people/severijn/.cvs-repository/EMIC/clio/sources/sc_cfc.f,v $
!_ $Revision: 1.1.1.1 $
!_ $Date: 2001/05/21 14:42:20 $   ;  $State: Exp $
!_ $Author: severijn $ ;  $Locker:  $
!_
!_ ---------------------------------------------------------------------
!_ $Log: sc_cfc.f,v $
!_ Revision 1.1.1.1  2001/05/21 14:42:20  severijn
!_ Import of ECBILT/CLIO
!_
!_ Revision 1.1  1998/07/07 15:22:00  orr
!_ Initial revision
!_
!_ ---------------------------------------------------------------------
!_
!jmc  REAL FUNCTION sc_cfc(t,kn)
!---------------------------------------------------
!     CFC 11 and 12 schmidt number
!     as a fonction of temperature.
!
!     ref: Zheng et al (1998), JGR, vol 103,No C1
!
!     t: temperature (degree Celcius)
!     kn: = 11 for CFC-11,  12 for CFC-12
!
!     J-C Dutay - LSCE
!---------------------------------------------------
      IMPLICIT NONE
      INTEGER kn
      REAL  a1 ( 11: 12), a2 ( 11: 12), a3 ( 11: 12), a4 ( 11: 12)
      REAL t, sc_cfc
!
!   coefficients with t in degre Celcius
!   ------------------------------------
      a1(11) = 3501.8
      a2(11) = -210.31
      a3(11) =    6.1851
      a4(11) =   -0.07513
!
      a1(12) = 3845.4
      a2(12) = -228.95
      a3(12) =    6.1908
      a4(12) =   -0.067430
!

      sc_cfc = a1(kn) + a2(kn) * t + a3(kn) *t*t  
     &         + a4(kn) *t*t*t

      RETURN
      END
