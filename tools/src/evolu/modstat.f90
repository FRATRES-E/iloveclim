MODULE MODSTAT
  IMPLICIT NONE
  REAL, PARAMETER :: PI=3.14159265

  ! MEAN(x)			 Mittelwert
  ! STDDEV(x)		 Standard Abweichung vom Mittelwert
  ! MEAN_FD(ff,dd)	 Mittelwert fuer Wind Vektoren
  ! STDDEV_FD(ff,dd) Standard Abweichung fuer Wind Vektoren
  ! COVAR(x,y)		 Kovarianz 
  ! VAR(x)			 Varianz
  ! PotReg(x,y)		 Calculate Regression
  
  CONTAINS

    REAL FUNCTION MEAN(M)
      REAL, DIMENSION(:), INTENT(IN) :: M
      INTEGER :: I,N

      N=0
      MEAN = 0.
      DO I=1,SIZE(M)
         IF (M(I).NE.-1) THEN
            MEAN = MEAN + M(I)
            N=N+1
         ENDIF
      ENDDO
      MEAN = MEAN/N
    END FUNCTION MEAN
    
    REAL FUNCTION STDDEV(M)
      REAL, DIMENSION(:), INTENT(IN) :: M
      REAL :: MITTEL
      INTEGER :: N,I

      N=0
      MITTEL = MEAN(M)
      STDDEV = 0.
      DO I=1,SIZE(M)
         IF (M(I).NE.-1) THEN
            STDDEV = STDDEV + (M(I)-MITTEL)**2
            N=N+1
         ENDIF
      ENDDO
      STDDEV = SQRT(STDDEV/(N-1))
    END FUNCTION STDDEV

	! MEAN for Wind (speed, direction)
    FUNCTION MEAN_FD(FF,DD)
      REAL, DIMENSION(:), INTENT(IN) :: FF,DD
      REAL :: MEAN_FD(2)
      REAL :: u,v
      INTEGER :: I,N
      
      N=0
      u = 0.
      v = 0.
      MEAN_FD = 0.
      DO I=1,SIZE(FF)
         IF (FF(I).NE.-1.AND.DD(I).NE.-1) THEN
            u = u + FF(I)*SIN(DD(I)*PI/180.)
            v = v + FF(I)*COS(DD(I)*PI/180.)
            N=N+1
         ENDIF
      ENDDO
      MEAN_FD(1) = SQRT((u/N)**2+(v/N)**2)
      MEAN_FD(2) = MOD((360. + ATAN2(v/N,u/N)*180./PI),360.)
    END FUNCTION MEAN_FD

	! STDDEV for Wind (speed, direction)
    FUNCTION STDDEV_FD(FF,DD)
      REAL, DIMENSION(:), INTENT(IN) :: FF,DD
      REAL :: STDDEV_FD(2),MITTEL(2)
      REAL :: u,v
      INTEGER :: I,N
      
      N=0
      u = 0.
      v = 0.
      STDDEV_FD = 0.
      MITTEL = 0.
      DO I=1,SIZE(FF)
         IF (FF(I).NE.-1.AND.DD(I).NE.-1) THEN
            MITTEL(1) = MITTEL(1) + FF(I)*SIN(DD(I)*PI/180.)
            MITTEL(2) = MITTEL(2) + FF(I)*COS(DD(I)*PI/180.)
            N=N+1
         ENDIF
      ENDDO
      MITTEL = MITTEL /N
      
      DO I=1,SIZE(FF)
         IF (FF(I).NE.-1.AND.DD(I).NE.-1) THEN
            u = u + (FF(I)*SIN(DD(I)*PI/180.)-MITTEL(1))**2
            v = v + (FF(I)*COS(DD(I)*PI/180.)-MITTEL(2))**2
            N=N+1
         ENDIF
      ENDDO
      STDDEV_FD(1) = SQRT((u/(N-1))**2+(v/(N-1))**2)
      STDDEV_FD(2) = MOD((360. + ATAN2(v/N,u/N)*180./PI),360.)
    END FUNCTION STDDEV_FD

    REAL FUNCTION COVAR(X,Y)
      REAL, DIMENSION(:), INTENT(IN) :: X,Y
      REAL :: XM,YM
      INTEGER :: N,I

      IF (SIZE(X).EQ.SIZE(Y)) THEN
         N=0
         XM = MEAN(X)
         YM = MEAN(Y)
         COVAR = 0 
         DO I=1,SIZE(X)
            IF (X(I).NE.-1.AND.Y(I).NE.-1) THEN
               COVAR = COVAR + (X(I)-XM)*(Y(I)-YM)
               N=N+1
            ENDIF
         ENDDO
         COVAR = COVAR/(N-1)
      ELSE
         COVAR=-999.
      ENDIF
    END FUNCTION COVAR
    
    REAL FUNCTION VAR(X)
      REAL, DIMENSION(:), INTENT(IN) :: X
      REAL :: XM
      INTEGER :: I,N

      VAR=0
      N=0
      XM = MEAN(X)
      DO I=1,SIZE(X)
         IF (X(I).NE.-1) THEN
            VAR = VAR + (X(I)-XM)**2
            N=N+1
         ENDIF
      ENDDO
      VAR = VAR/(N-1)
    END FUNCTION VAR
    
    FUNCTION PotReg(X,Y)
      REAL, DIMENSION(:), INTENT(IN) :: X,Y
      REAL, ALLOCATABLE :: PotReg(:),XL(:),YL(:)
      REAL :: AL,A,B
      INTEGER :: I,N
      
      !? >0?
      N = SIZE(X)
      allocate(PotReg(N),XL(N),YL(N))
      IF (all(X.NE.-1).AND.all(Y.NE.-1)) THEN
         !no prob
         XL = LOG10(X)
         YL = LOG10(Y)
      ELSE
         DO I=1,N
            IF (X(I).NE.-1.AND.Y(I).NE.-1) THEN
               XL(I) = LOG10(X(I))
               YL(I) = LOG10(Y(I))
            ENDIF
         ENDDO
      ENDIF
      
      B = COVAR(XL,YL)/VAR(XL)
      AL = MEAN(YL)-B*MEAN(XL)
      A = 10**AL
      WRITE(*,*) 'PotReg: ',A,'*VV^(',B,')'
      PotReg = A*X**(B)
      
    END FUNCTION PotReg
END MODULE MODSTAT
