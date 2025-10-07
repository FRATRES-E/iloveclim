! module de declaration de variables
! ocean carbonate
! cnb

module carbonate_mod

USE declars_mod

implicit none


!***********************************************************************
!      COMMON /OCEANS/

!* VARIABLES
!*************
REAL, dimension(JT) :: CO3sat_ar,CO3sat_ca
REAL, dimension(LT,JT,NOC_CBR) :: sat_ar, sat_ca

end module carbonate_mod
