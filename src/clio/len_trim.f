!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:49 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:49 CET 2009

      FUNCTION len_trim(STRING)
!_ ---------------------------------------------------------------------
!_ RCS lines preceded by "c_ "
!_ ---------------------------------------------------------------------
!_
!_ $Source: /usr/people/severijn/.cvs-repository/EMIC/clio/sources/len_trim.f,v $
!_ $Revision: 1.1.1.1 $
!_ $Date: 2001/05/21 14:42:20 $   ;  $State: Exp $
!_ $Author: severijn $ ;  $Locker:  $
!_
!_ ---------------------------------------------------------------------
!_ $Log: len_trim.f,v $
!_ Revision 1.1.1.1  2001/05/21 14:42:20  severijn
!_ Import of ECBILT/CLIO
!_
!_ Revision 1.1  1998/07/21 16:55:54  jomce
!_ Initial revision
!_
!_ ---------------------------------------------------------------------
!_
!jmc        FUNCTION LEN_TRIM(STRING)

        CHARACTER*(*) STRING
        INTEGER*4 I 
        INTEGER*4 LEN_TRIM

        DO I=LEN(STRING),1,-1
                LEN_TRIM=I
                IF (STRING(I:I).NE.' ') RETURN
        enddo

        END
