!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:21 CET 2009

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       module O2SAT_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       use global_constants_mod, only: dblp=>dp, ip
       use declars_mod, only: lt, noc_cbr
       use marine_bio_mod, only: jprod

       implicit none 
       
       real(kind=dblp), dimension(lt,jprod,noc_cbr) :: O2_sat_thistime

       contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       subroutine O2sat_in_photic_zone()
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       
         use loveclim_transfer_mod, ONLY: TM, SM, MGT
          
          integer(kind=ip) :: i,j,n
          
          do n=1,NOC_CBR
            do j=1,JPROD          
              do i=1,LT
                if (MGT(i,j,n).eq.1) then          
                   O2_sat_thistime(i,j,n) = O2sat(TM(i,j,n),SM(i,j,n))
                endif
              enddo
            enddo           
          enddo              
          
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|       
       end subroutine O2sat_in_photic_zone
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        real(kind=dblp) function O2sat(T,S)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!
!  By: J. Bendtsen
!  Last modification: 02.02.99
!  -- f90 by dmr, 2020-06-15
!
!  Purpose: calculate the O2 saturation concentration
!
!       T     : TEMPERATURE (CELCIUS DEGREES)
!       S     : SALINITY (0/00)
!       O2sat : O2 SATURATION (mumol/kg)
!
        
        real(kind=dblp), intent(in):: T, S
        real(kind=dblp)            :: A, xK100!, O2sat                 
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        xK100   = (T+273.15)/100.

        A      = -177.7888+255.5907/xK100+146.4813*ALOG(xK100)-22.204*xK100
        A      = A + S*(-0.037362+xK100*(0.016504-0.0020564*xK100))

        O2sat = exp(A) *1.0e3/22.392 !  units=mumol/kg 

        end function O2sat

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Version taken from the MOCSY PACKAGE, subroutine gasx.f90
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!>    Compute O2 saturation concentration of surface seawater (mol/m3) at 1 atm (Garcia & Gordon, 1992)
      subroutine o2sato(T, S, N, o2sat_molm3)
  !    Purpose:
  !    Compute O2 saturation concentration of surface seawater (mol/m3) at 1 atm (Garcia & Gordon, 1992)
  !
  !    ********************************************************************
  !    Computes the oxygen saturation concentration at 1 atm total pressure
  !    in mol/m^3 given sea surface temperature T (deg C) and salinity S (permil) 
  !
  !    From: Garcia & Gordon (1992) Oxygen solubility in seawater: better fitting equations,
  !          Limnol. Oceanogr., 37(6), 1307-1312.
  !          This routine uses:
  !          - equation (8) on page 1310
  !          - coefficients from Table 1, column 2 (Benson & Krause, [cm3/dm3], i.e, same as [ml/L])
  !
  !    *** NOTE: The "A3*ts^2" term in the equation (8) in the paper is a TYPO.    ***
  !    *** It shouldn't be there. It is not used in this routine.                  ***
  !
  !    'o2sat' is fit between T(freezing) <= T <= 40(deg C)  and  0 <= S <= 42 permil
  !
  !    CHECK VALUE:  T = 10.0 deg C, S = 35.0 permil, 
  !    o2sat_molm3 = 0.282015 mol/m^3
  !    ********************************************************************

      

  !> number of records
      INTEGER(kind=ip), intent(in) :: N

! INPUT variables
  !> surface temperature [C]
      REAL(kind=dblp), INTENT(in), DIMENSION(N) :: T
  !> surface salinity [psu]
      REAL(kind=dblp), INTENT(in), DIMENSION(N) :: S
!!!f2py optional , depend(temp) :: n=len(temp)
!f2py depend(temp) :: n

! OUTPUT variables:
  !> O2 saturation concentration of seawater [mol/m3] 
      REAL(kind=dblp), INTENT(out), DIMENSION(N) :: o2sat_molm3

! LOCAL variables:
      REAL(kind=dblp)  :: A0, A1, A2, A3, A4, A5
      REAL(kind=dblp)  :: B0, B1, B2, B3
      REAL(kind=dblp)  :: C0
      REAL(kind=dblp)  :: tmp
      REAL(kind=dblp)  :: o2sat_mlL 
      REAL(kind=dblp)  :: tt, tk, ts, ts2, ts3, ts4, ts5
      INTEGER(kind=ip) :: i
  
      DATA A0/ 2.00907_dblp   /, A1/ 3.22014_dblp   /, A2/ 4.05010_dblp /,  &
           A3/ 4.94457_dblp   /, A4/-2.56847E-1_dblp/, A5/ 3.88767_dblp /
      DATA B0/-6.24523E-3_dblp/, B1/-7.37614E-3_dblp/, B2/-1.03410E-2_dblp/, B3/-8.17083E-3_dblp/
      DATA C0/-4.88682E-7_dblp/

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Main Code starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
       DO i = 1, N
         tt  = 298.15_dblp - T(i)
         tk  = 273.15_dblp + T(i)
         ts  = LOG(tt/tk)

         ts2 = ts**2
         ts3 = ts**3
         ts4 = ts**4
         ts5 = ts**5

!     O2 saturation concentration (ml/L) 
         tmp  = A0 + A1*ts + A2*ts2 + A3*ts3 + A4*ts4 + A5*ts5 + S(i)*(B0 + B1*ts + B2*ts2 + B3*ts3) + C0*(S(i)*S(i))
         o2sat_mlL = EXP(tmp)

!     Convert from ml/L to mol/m^3
         o2sat_molm3(i) = o2sat_mlL / 22391.6_dblp*1000.0_dblp
       END DO
   
       RETURN
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|       
       END SUBROUTINE o2sato
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|       
       end module O2SAT_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
