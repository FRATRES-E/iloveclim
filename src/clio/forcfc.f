!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:46 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:46 CET 2009

!_ ---------------------------------------------------------------------
!_ RCS lines preceded by "c_ "
!_ ---------------------------------------------------------------------
!_
!_ $Source: /usr/people/severijn/.cvs-repository/EMIC/clio/sources/forcfc.f,v $
!_ $Revision: 1.1.1.1 $
!_ $Date: 2001/05/21 14:42:20 $   ;  $State: Exp $
!_ $Author: severijn $ ;  $Locker:  $
!_
!_ ---------------------------------------------------------------------
!_ $Log: forcfc.f,v $
!_ Revision 1.1.1.1  2001/05/21 14:42:20  severijn
!_ Import of ECBILT/CLIO
!_
!_ Revision 1.1  1998/07/07 15:22:00  orr
!_ Initial revision
!_
!_ ---------------------------------------------------------------------
!_
      SUBROUTINE cfc_interp(pnorth, psouth, rlat, imt, jmt, patm) 
!     SUBROUTINE cfc_interp(pnorth, psouth, xintp, rlat, imt, jmt, patm) 

!     ==================================================================
!     ROUTINE to interpolate CFC-11 and CFC-12 DATA for OCMIP-2 model runs

!     Two stations (41S and 45N), each representative of its own hemisphere,
!     except between 10S and 10N WHERE values are interpolated linearly
!     thus there are 3 regions:
!     (1) 90s to 10s, WHERE values take on the value at the station at 41s;
!     (2) 10n to 90n, WHERE values take on the value at the station at 45n; and
!     (3) 10s to 10n, WHERE values are interpolated

!     ---------------------------------------------------------
!     James Orr, LSCE/CEA-CNRS, Saclay, France, 10 June 1998
!     Core algortihm of this ROUTINE by Jean-Claude Dutay, LSCE
!     ---------------------------------------------------------

!     INPUT: 
!     ------
!     Pnorth           REAL   2-member array of CFC-11 and CFC-12 
!                             atmospheric concentrations (mixing ratio) at
!                             the northern station (45N) at time t.
!     Psouth           REAL   2-member array of CFC-11 and CFC-12 
!                             atmospheric concentrations (mixing ratio) at
!                             the southern station (41S) at time t.
!     rlat             REAL   2-D array of latitudes, needed to compute xintp
!     imt           INTEGER   max. dimension in direction i   
!     jmt           INTEGER   max. DIMENSION in direction j

!     OUTPUT:
!     -------
!     Patm             REAL   spatial array for computed atmospheric 
!                             CFC-11 and CFC-12 
!
!     LOCAL:
!     ------
!     xintp            REAL   2-D array of interpolation factors computed on 
!                             1st pass through this routine
!     ientry        INTEGER   Scalar counter to indicate number of calls to
!                             this routine by main PROGRAM
!     ys, yn           REAL   Scaler Latitude limits, within which one does 
!                             a linear interpolation (i.e., 10S to 10N)

!     ==================================================================
      IMPLICIT NONE

      INTEGER nt, ijmax
      PARAMETER (nt=2, ijmax=200*200)

      INTEGER imt, jmt 
      INTEGER ientry
      INTEGER i, j, ij, n     

      REAL ys, yn 
      REAL rlat(imt,jmt), patm(imt,jmt,nt)

      REAL xintp(ijmax) 
      REAL Pnorth(nt), Psouth(nt)

      SAVE ys, yn, ientry, xintp

      DATA ys, yn / -10., 10./
      DATA ientry /0/

      ientry = ientry + 1

!     Test to see if ijmax is large enough
      IF (ijmax .LT. imt*jmt) THEN
          PRINT *," cfc_interp: ERROR -> ijmax must be at least "
     $             ,        imt*jmt
          STOP
      ENDIF

!     IF Block to be executed ONLY on first call to this routine
      IF (ientry .EQ. 1) THEN
          DO j = 1,jmt
            DO i = 1,imt
              ij = (j-1)*imt + i
              IF (rlat(i,j) .GE. yn) THEN
                  xintp(ij) = 1.
              ELSE IF (rlat(i,j) .LE. ys) THEN
                  xintp(ij) = 0.
              ELSE 
                  xintp(ij) =  (rlat(i,j) - ys)/(yn-ys)
              ENDIF 
            enddo
          enddo 
      ENDIF 

!     Block to be executed every pass
      DO n=1,nt
        DO j=1,jmt
          DO i=1,imt
            ij = (j-1)*imt + i
            Patm(i,j,n) =         xintp(ij)  * Pnorth(n) 
     $                   + (1.0 - xintp(ij)) * Psouth(n)
          enddo
        enddo
      enddo

      RETURN
      END
