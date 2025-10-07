!     dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:27 CET 2009
#include "choixcomposantes.h"
!     dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:27 CET 2009

#if ( BIOM_GEN >= 1 )
      
      MODULE DIAGNOSE_BIOME_MOD

      USE global_constants_mod, ONLY: silp=>sp, big_dp => alt_olympus_mons

      IMPLICIT NONE

      REAL(KIND=silp), parameter :: undef = big_dp*-1._silp

      CONTAINS


! dmr --- Input variables:
!                         min_T(size_Y,size_X,nbrecordmax)
!                         max_T(size_Y,size_X,nbrecordmax)
!                         vegmapout(size_Y,size_X,nvgmax,nbrecordmax)
!                               from which:
!                                           nvg ==  1 => tree fraction (and undef for not land)
!                                           nvg ==  4 => needle leaf tree fraction
!                                           nvg == 20 => GDD5
!                                           nvg ==  2 => grass fraction
!                                           nvg ==  3 => desert fraction
!                                           nvg ==  9 => GDD0
!                                           nvg == 19 => icemask % [0-1] theoretically


!~       SUBROUTINE determine_biome(min_T,max_T,vegmapout,size_Y, size_X, nvgmax, nbrecordmax,nbrecordeff,biomeout)
      SUBROUTINE determine_biome(min_T,max_T,treefrac, needletreefrac, grassfrac, desertfrac, GDD0_val, GDD5_val, icemask_frac &
                                ,size_Y, size_X, nbrecordmax,nbrecordeff,biomeout)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~       INTEGER, INTENT(IN) :: size_Y, size_X, nvgmax, nbrecordmax, nbrecordeff
      INTEGER, INTENT(IN) :: size_Y, size_X, nbrecordmax, nbrecordeff
      REAL(KIND=4), DIMENSION(size_Y, size_X, nbrecordmax), INTENT(IN) :: min_T, max_T, treefrac, needletreefrac, grassfrac   &
                             , desertfrac, GDD0_val, GDD5_val, icemask_frac
      REAL(KIND=4), DIMENSION(size_Y, size_X, nbrecordmax), INTENT(OUT):: biomeout
!~       REAL(KIND=4), DIMENSION(size_Y, size_X, nvgmax, nbrecordmax), INTENT(IN) :: vegmapout

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      INTEGER i,j,k
      REAL(KIND=4) :: minit, maxit
!~       REAL(kind=8) :: cpu_start, cpu_finish

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      biomeout(:,:,:) = undef

!~       call cpu_time(cpu_start)


      do k=1,nbrecordeff
        do j=1,size_X
          do i=1,size_Y      
            
          
       minit = min_T(i,j,k)
       maxit = max_T(i,j,k)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       NO LAND POINT
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      IF (icemask_frac(i,j,k).GT.undef) THEN ! GOTO 900

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       TREE DOMINANCE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      IF (treefrac(i,j,k).GT.50.) THEN                                                                                   ! Majorite d'arbres

        IF ((treefrac(i,j,k)*(1.-needletreefrac(i,j,k)/100.)).GE.(treefrac(i,j,k)*(needletreefrac(i,j,k)/100.))) THEN       ! Arbres a feuilles dominant

        IF (minit.GT.18.0) then  ! minT > 18C
          biomeout(i,j,k) = 1.0                                         ! --- TROPICAL FOREST

          ELSEIF ((minit.GT.-2.0).AND.(minit.LT.18.0) ) THEN  ! -2C < minT < 18C
            biomeout(i,j,k) = 3.0                                       ! --- TEMPERATE FOREST

          ELSE ! en theorie, pas possible !
            biomeout(i,j,k) = 6.0                                       ! --- BOREAL FOREST
          ENDIF

        ELSEIF ((treefrac(i,j,k)*(1.-needletreefrac(i,j,k)/100.)).LT.(treefrac(i,j,k)*(needletreefrac(i,j,k)/100.))) THEN   ! Arbres a aiguilles dominant

           IF  ((minit.GT.-19.0).AND.(minit.LT.-15.0)) THEN  ! -19C < minT < -15C
             biomeout(i,j,k) = 5.0                                      ! --- COLD CONIFEROUS FOREST

           ELSEIF (((minit.GT.-15.0).AND.(minit.LT.-2.0)).AND.(GDD5_val(i,j,k).LT.1200.)) THEN ! -15C < minT < -2C && gdd5 < 1200
             biomeout(i,j,k) = 5.0                                      ! --- COLD CONIFEROUS FOREST

           ELSE
             IF ((minit.LT.-15.0)) THEN ! -35C < minT < -15C
                  biomeout(i,j,k) = 6.0                                 ! --- BOREAL FOREST

             ELSEIF ((minit.GT.-2.0).AND.(minit.LT.7.0)) THEN ! -2C < minT < 7C
                  biomeout(i,j,k) = 4.0                                 ! --- COLD MIXED FOREST

             ELSE
                  biomeout(i,j,k) = 2.0                                 ! --- WARM MIXED FOREST

             ENDIF

           ENDIF

        ELSEIF ((minit.GT.5.0).AND.(minit.LT.18)) THEN   ! 5C < minT < 18C
           biomeout(i,j,k) = 2.0                                        ! --- WARM MIXED FOREST

        ELSEIF ((minit.GT.-2.0).AND.(minit.LT.5.0)) THEN ! -2C < minT < 5C
           biomeout(i,j,k) = 4.0                                        ! --- COLD MIXED FOREST

        ELSEIF (minit.LT.-2.0) THEN  ! minT < -2C
           biomeout(i,j,k) = 6.0                                        ! --- BOREAL FOREST

        ENDIF

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       GRASS DOMINANCE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      ELSEIF (grassfrac(i,j,k).GT.50) THEN                            ! Majorite de gazon
        IF (minit.GT.18) THEN    ! minT > 18
          biomeout(i,j,k) = 7.0                                         ! --- TROPICAL SAVANNAH

        ELSEIF ( maxit.GT.32.0 ) THEN ! maxT > 32C
          biomeout(i,j,k) = 8.0                                         ! --- WARM GRASSLAND

        ELSEIF ((maxit.LT.32.0).AND.(GDD5_val(i,j,k).GT.500)) THEN  ! maxT < 32C && gdd5 > 500
          biomeout(i,j,k) = 9.0                                         ! --- COLD GRASSLAND

        ELSEIF ((GDD5_val(i,j,k).LT.500)) THEN                      ! gdd5 < 500
          biomeout(i,j,k) = 10.0                                        ! --- TUNDRA

        ENDIF

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       DESERT DOMINANCE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      ELSEIF ((grassfrac(i,j,k)+treefrac(i,j,k)).GT.30.0) THEN     ! arbres + gazon > 30%
          biomeout(i,j,k) = 11.0                                        ! --- SEMI-DESERTIC

      ELSEIF ((desertfrac(i,j,k).GT.50.0).AND.((minit.GT.10.0).OR.(maxit.GT.20.0)) ) THEN ! desert > 50% && minT > 10C ou maxT > 20C
          biomeout(i,j,k) = 13.0                                        ! --- WARM DESERT

      ELSEIF ((icemask_frac(i,j,k).GT.0.5).OR.(GDD0_val(i,j,k).LE.500).OR.(GDD5_val(i,j,k).LE.300)) THEN ! calottes glace > 50% && gdd0 <= 100 || gdd0 <500 || gdd5 < 300
          biomeout(i,j,k) = 12.                                         ! --- COLD DESERT

      ELSEIF ((desertfrac(i,j,k).GT.80.0) ) THEN ! desert > 80% 
          biomeout(i,j,k) = 14.0                                        ! --- OTHER DESERT
          
      ENDIF

      IF (((icemask_frac(i,j,k).GT.undef).AND.(biomeout(i,j,k).LE.undef))) THEN
        PRINT*, "Biome non assigne en case : ", i,j,k
        PRINT*, "arbres : ", treefrac(i,j,k)
        PRINT*, "gazon : ", grassfrac(i,j,k)
        PRINT*, "desert : ", desertfrac(i,j,k)
        PRINT*, "minT :", minit
        PRINT*, "maxT :", maxit
        PRINT*, "Broadleaf :", (treefrac(i,j,k)*(1.-needletreefrac(i,j,k)/100.))
        PRINT*, "Needleleaf :", (treefrac(i,j,k)*(needletreefrac(i,j,k)/100.))
        PRINT*, "GDD5 :", GDD5_val(i,j,k)
        PRINT*, "GDD0 :", GDD0_val(i,j,k)
        PRINT*, "icemask :", icemask_frac(i,j,k)
      ENDIF


      ENDIF

      ENDDO
      ENDDO
      ENDDO
!~       call cpu_time(cpu_finish)
!~       print '("Time = ",f6.3," seconds.")',cpu_finish-cpu_start

      END SUBROUTINE determine_biome

      END MODULE DIAGNOSE_BIOME_MOD

#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
