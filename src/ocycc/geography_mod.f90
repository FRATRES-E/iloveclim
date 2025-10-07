      MODULE geography_mod

      IMPLICIT NONE

      !dmr --- Pour test de compilation
      INTEGER, PARAMETER :: IT=32
      INTEGER, PARAMETER :: IX=IT+1
      INTEGER, PARAMETER :: NS=64
      !dmr --- Pour test de compilation

      !************* GEOGRAPHY **************************************
      REAL, dimension(IT) :: FIT, SINT, COST, DX
      REAL, dimension(IX) :: FIX, SINU, COSU
      REAL :: DY, DPHI, SQA
      REAL :: SQO
      REAL, TARGET :: SQL
      REAL, dimension(IT,NS), TARGET :: SQRA
      REAL, dimension(IT,NS) :: SQRO


      END MODULE geography_mod
