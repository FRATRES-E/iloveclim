!     dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:27 CET 2009
#include "choixcomposantes.h"
!     dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:27 CET 2009

#if ( CLM_INDICES >= 1 )

      MODULE CLIMATE_INDICES_MOD

!     [NOTA]: all temperature variables are assumed to be in K


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      USE global_constants_mod, ONLY: dblp=>dp, silp=>sp, sip,  freezeT => tK_zero_C, big_dp => alt_olympus_mons
      USE comatm, ONLY: size_Y=>nlat, size_X=>nlon

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: DAILYSTEP_FOR_CLIM_INDICES, GLOBAL_RE_INIT, GLOBAL_FINALIZE

#if ( CLM_INDICES >= 2 )      
      PUBLIC :: SET_VEG_VARS
#endif      


      REAL(KIND=dblp), PARAMETER :: bigN_dp = -1.0*big_dp, d_zero = 0.0_dblp, d_five = 5.0_dblp


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       CLIMATE INDICES DEFINITIONS
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       INTEGER(KIND=sip),DIMENSION(size_Y,size_X) :: frost_days_0   ! 01
       REAL(KIND=dblp),  DIMENSION(size_Y,size_X) :: TmaX           ! 02
       REAL(KIND=dblp),  DIMENSION(size_Y,size_X) :: TmiN           ! 03
       REAL(KIND=dblp),  DIMENSION(size_Y,size_X) :: TmaX_y         ! 04
       REAL(KIND=dblp),  DIMENSION(size_Y,size_X) :: TmiN_y         ! 05
       REAL(KIND=dblp),  DIMENSION(size_Y,size_X) :: GDD_zero       ! 06
       REAL(KIND=dblp),  DIMENSION(size_Y,size_X) :: GDD_five       ! 07
       
#if ( CLM_INDICES >= 2 )
       REAL(KIND=dblp),  DIMENSION(size_Y,size_X) :: treefrac       ! 08       
       REAL(KIND=dblp),  DIMENSION(size_Y,size_X) :: needletreefrac ! 09
       REAL(KIND=dblp),  DIMENSION(size_Y,size_X) :: grassfrac      ! 10
       REAL(KIND=dblp),  DIMENSION(size_Y,size_X) :: desertfrac     ! 11
#endif

#if ( CLM_INDICES == 1 )              
       INTEGER(KIND=sip), PARAMETER :: number_indices = 7
#elif ( CLM_INDICES == 2 )
       INTEGER(KIND=sip), PARAMETER :: number_indices = 11
#endif      

       INTEGER(KIND=sip), PARAMETER :: frost_days_0_key=1
       INTEGER(KIND=sip), PARAMETER :: TmaX_key=2
       INTEGER(KIND=sip), PARAMETER :: TmiN_key=3
       INTEGER(KIND=sip), PARAMETER :: TmaX_y_key=4
       INTEGER(KIND=sip), PARAMETER :: TmiN_y_key=5
       INTEGER(KIND=sip), PARAMETER :: GDD_zero_key=6
       INTEGER(KIND=sip), PARAMETER :: GDD_five_key=7

#if ( CLM_INDICES >= 2 )       
       INTEGER(KIND=sip), PARAMETER :: treefrac_key=8
       INTEGER(KIND=sip), PARAMETER :: needletreefrac_key=9
       INTEGER(KIND=sip), PARAMETER :: grassfrac_key = 10
       INTEGER(KIND=sip), PARAMETER :: desertfrac_key = 11
#endif       

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       FILE NAMES DEFINITIONS AND HANDLING
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       CHARACTER(LEN=8),  PARAMETER :: indices_suff="_indices"
       CHARACTER(LEN=4),  PARAMETER :: binaary_suff=".grd"
       CHARACTER(LEN=4),  PARAMETER :: control_suff=".ctl"
       CHARACTER(LEN=27), PARAMETER :: directoire_stuff="outputdata/climate_indices/"
       
#if ( CLM_INDICES == 1 )          
       CHARACTER(LEN=4), DIMENSION(number_indices), PARAMETER ::                               &
                                             indices_filenames= (/"FRD0","TMMX","TMMN","TYMX","TYMN","GDD0","GDD5"/)
#elif ( CLM_INDICES == 2 )                   
       CHARACTER(LEN=4), DIMENSION(number_indices), PARAMETER ::                               &
                                             indices_filenames= (/"FRD0","TMMX","TMMN","TYMX","TYMN","GDD0","GDD5",   & 
                                                                  "TRFR","NTFR","GRFR","DSFR"/)
#endif

       CHARACTER(LEN=43),DIMENSION(number_indices)            :: binaary_filenames, control_filenames

      CONTAINS

! ---

      SUBROUTINE GLOBAL_RE_INIT()

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       THINGS TO DO ONCE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


        binaary_filenames(:) = directoire_stuff//indices_filenames(:)//indices_suff//binaary_suff
        control_filenames(:) = directoire_stuff//indices_filenames(:)//indices_suff//control_suff

        CALL MONTHLY_RE_INIT()
        CALL YEARLY_RE_INIT()


      END SUBROUTINE GLOBAL_RE_INIT

! ---

! ---

      SUBROUTINE GLOBAL_FINALIZE()

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! write the CTL files

        CALL WRITE_CTL_HEADER(control_filenames(frost_days_0_key),binaary_filenames(frost_days_0_key)         &
                            , indices_filenames(frost_days_0_key),'y','i')
        CALL WRITE_CTL_HEADER(control_filenames(TmaX_key)        ,binaary_filenames(TmaX_key)                 &
                            , indices_filenames(TmaX_key)        ,'m','f')
        CALL WRITE_CTL_HEADER(control_filenames(TmiN_key)        ,binaary_filenames(TmiN_key)                 &
                            , indices_filenames(TmiN_key)        ,'m','f')
        CALL WRITE_CTL_HEADER(control_filenames(TmaX_y_key)      ,binaary_filenames(TmaX_y_key)               &
                            , indices_filenames(TmaX_y_key)      ,'y','f')
        CALL WRITE_CTL_HEADER(control_filenames(TmiN_y_key)      ,binaary_filenames(TmiN_y_key)               &
                            , indices_filenames(TmiN_y_key)      ,'y','f')
        CALL WRITE_CTL_HEADER(control_filenames(GDD_zero_key)    ,binaary_filenames(GDD_zero_key)             &
                            , indices_filenames(GDD_zero_key)    ,'y','f')
        CALL WRITE_CTL_HEADER(control_filenames(GDD_five_key)    ,binaary_filenames(GDD_five_key)             &
                            ,indices_filenames(GDD_five_key)    ,'y','f')
                            
#if ( CLM_INDICES >= 2 )   

        CALL WRITE_CTL_HEADER(control_filenames(treefrac_key)      ,binaary_filenames(treefrac_key)           &
                            , indices_filenames(treefrac_key)      ,'y','f')
        CALL WRITE_CTL_HEADER(control_filenames(needletreefrac_key),binaary_filenames(needletreefrac_key)     &
                            ,indices_filenames(needletreefrac_key) ,'y','f')
        CALL WRITE_CTL_HEADER(control_filenames(grassfrac_key)     ,binaary_filenames(grassfrac_key)          &
                            ,indices_filenames(grassfrac_key)      ,'y','f')                            
        CALL WRITE_CTL_HEADER(control_filenames(desertfrac_key)    ,binaary_filenames(desertfrac_key)         &
                            ,indices_filenames(desertfrac_key)     ,'y','f')                                    

#if ( BIOM_GEN == 1 )
       CALL COMPUTE_AND_WRITE_BIOMES()
#endif

#endif                            

      END SUBROUTINE GLOBAL_FINALIZE

! ---

      SUBROUTINE YEARLY_RE_INIT()

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         frost_days_0(:,:) = 0.0_dblp
         GDD_zero(:,:)     = 0.0_dblp
         GDD_five(:,:)     = 0.0_dblp
         TmaX_y(:,:)       = bigN_dp
         TmiN_y(:,:)       = big_dp

      END SUBROUTINE YEARLY_RE_INIT

      SUBROUTINE MONTHLY_RE_INIT()

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         TmaX(:,:) = bigN_dp
         TmiN(:,:) = big_dp

      END SUBROUTINE MONTHLY_RE_INIT

! ---


      SUBROUTINE DAILYSTEP_FOR_CLIM_INDICES()

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     NOTA: currently called every day, not at each atmospheric timestep, since there is no diurnal variability in ECBilt
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        use comoutlocal_mod, only: tsurf1
        use comemic_mod, only: iday, imonth, iyear
        use global_constants_mod, only: months_year_i, days_month360d_i, ip

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER(kind=sip) :: i,j,c_month, c_year

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         ! WRITE(*,*) 'i_day, i_month == ', iday, imonth, iyear

         DO j=1, size_X
         DO i=1, size_Y
            IF (tsurf1(i,j) < freezeT) THEN
              frost_days_0(i,j) = frost_days_0(i,j) + 1
            ENDIF
            TmiN(i,j) = min(TmiN(i,j),tsurf1(i,j))
            TmaX(i,j) = max(TmaX(i,j),tsurf1(i,j))
            GDD_zero(i,j) = GDD_zero(i,j)+ max( (tsurf1(i,j) - freezeT ), d_zero )
            GDD_five(i,j) = GDD_five(i,j)+ max( (tsurf1(i,j) - freezeT - d_five), d_zero )
         ENDDO
         ENDDO

         if (iday == days_month360d_i) then
           c_month = (iyear-1)*months_year_i+imonth
           CALL DO_END_MONTH(c_month)
         endif

         if (imonth == months_year_i .and. iday == days_month360d_i) then
           c_year = iyear
           CALL DO_END_YEAR(c_year)
         endif

      END SUBROUTINE DAILYSTEP_FOR_CLIM_INDICES


! ---


      SUBROUTINE DO_END_MONTH(nb_month)

         INTEGER(KIND=sip), INTENT(IN) :: nb_month

         INTEGER(KIND=sip) :: i,j

         CALL WRITE_ONE_VAR_D(binaary_filenames(TmiN_key),nb_month,TmiN)
         CALL WRITE_ONE_VAR_D(binaary_filenames(TmaX_key),nb_month,TmaX)

         DO j=1, size_X
         DO i=1, size_Y
            TmiN_y(i,j) = min(TmiN_y(i,j),TmiN(i,j))
            TmaX_y(i,j) = max(TmaX_y(i,j),TmaX(i,j))
         ENDDO
         ENDDO

         CALL MONTHLY_RE_INIT()

      END SUBROUTINE DO_END_MONTH


! ---


      SUBROUTINE DO_END_YEAR(nb_year)

         INTEGER(KIND=sip), INTENT(IN) :: nb_year

         CALL WRITE_ONE_VAR_I(binaary_filenames(frost_days_0_key),nb_year,frost_days_0)
         CALL WRITE_ONE_VAR_D(binaary_filenames(GDD_zero_key)    ,nb_year,GDD_zero)
         CALL WRITE_ONE_VAR_D(binaary_filenames(GDD_five_key)    ,nb_year,GDD_five)
         CALL WRITE_ONE_VAR_D(binaary_filenames(TmaX_y_key)      ,nb_year,TmaX_y)
         CALL WRITE_ONE_VAR_D(binaary_filenames(TmiN_y_key)      ,nb_year,TmiN_y)

         CALL YEARLY_RE_INIT()

      END SUBROUTINE DO_END_YEAR

#if ( CLM_INDICES >= 2 )
! ---

      SUBROUTINE SET_VEG_VARS(tree,needle,grass,desert,nb_year)

         INTEGER(KIND=sip), INTENT(IN) :: nb_year
         REAL(KIND=dblp),  DIMENSION(:,:), INTENT(in) :: tree  
         REAL(KIND=dblp),  DIMENSION(:,:), INTENT(in) :: needle
         REAL(KIND=dblp),  DIMENSION(:,:), INTENT(in) :: grass
         REAL(KIND=dblp),  DIMENSION(:,:), INTENT(in) :: desert
         

         treefrac(:,:) = tree(:,:)
         needletreefrac(:,:) = needle(:,:)
         grassfrac(:,:) = grass(:,:)
         desertfrac(:,:) = desert(:,:)

         CALL DO_END_YEAR_VEG(nb_year)

      END SUBROUTINE SET_VEG_VARS


! ---

      SUBROUTINE DO_END_YEAR_VEG(nb_year)

         INTEGER(KIND=sip), INTENT(IN) :: nb_year

         CALL WRITE_ONE_VAR_D(binaary_filenames(treefrac_key)    ,nb_year,treefrac)
         CALL WRITE_ONE_VAR_D(binaary_filenames(needletreefrac_key)    ,nb_year,needletreefrac)
         CALL WRITE_ONE_VAR_D(binaary_filenames(grassfrac_key)      ,nb_year,grassfrac)
         CALL WRITE_ONE_VAR_D(binaary_filenames(desertfrac_key)      ,nb_year,desertfrac)

         CALL YEARLY_RE_INIT()

      END SUBROUTINE DO_END_YEAR_VEG


! ---

#endif

      SUBROUTINE WRITE_ONE_VAR_D(filename,curr_record,var_to_write)

        USE global_constants_mod, only: ip
        USE file_libs, only: get_fID, release_fID

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER(KIND=sip), INTENT(IN) :: curr_record
        REAL(KIND=dblp), DIMENSION(:,:), INTENT(IN) :: var_to_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER(kind=sip) :: i,j
        INTEGER(kind=ip)  :: file_ID

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        file_ID = get_fID()

        OPEN(file_ID,file=filename, access='direct', form='unformatted', status='unknown',recl=size_X*size_Y*1*silp)

        WRITE(file_ID,REC=curr_record) ((REAL(var_to_write(i,j),KIND=silp),j=1,size_X),i=1,size_Y)

        CLOSE(file_ID)

        call release_fID(file_ID)

      END SUBROUTINE WRITE_ONE_VAR_D


! ---

      SUBROUTINE WRITE_ONE_VAR_S_ALLREC(filename,var_to_write)

        USE global_constants_mod, only: ip
        USE file_libs, only: get_fID, release_fID

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        CHARACTER(LEN=*), INTENT(IN) :: filename
        REAL(KIND=silp), DIMENSION(:,:,:), INTENT(IN) :: var_to_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER(kind=sip) :: i,j,k
        INTEGER(kind=ip)  :: file_ID

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        file_ID = get_fID()

        OPEN(file_ID,file=filename, access='direct', form='unformatted', status='unknown',recl=size_X*size_Y*1*silp)
  
        DO k=LBOUND(var_to_write,DIM=3),UBOUND(var_to_write,DIM=3)
           WRITE(file_ID,REC=k) ((var_to_write(i,j,k),j=1,size_X),i=1,size_Y)
        ENDDO

        CLOSE(file_ID)

        call release_fID(file_ID)

      END SUBROUTINE WRITE_ONE_VAR_S_ALLREC

! ---


      SUBROUTINE WRITE_ONE_VAR_I(filename,curr_record,var_to_write)

        USE global_constants_mod, only: ip
        USE file_libs, only: get_fID, release_fID

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER(KIND=sip), INTENT(IN) :: curr_record
        INTEGER(KIND=sip), DIMENSION(:,:), INTENT(IN) :: var_to_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER(kind=sip) :: i,j
        INTEGER(kind=ip)  :: file_ID

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        file_ID = get_fID()

        OPEN(unit=file_ID,file=filename, access='direct', form='unformatted', status='unknown',recl=size_X*size_Y*1*sip)

        WRITE(unit=file_ID,REC=curr_record) ((var_to_write(i,j),j=1,size_X),i=1,size_Y)

        CLOSE(unit=file_ID)

        call release_fID(file_ID)

      END SUBROUTINE WRITE_ONE_VAR_I


! ---

      SUBROUTINE READ_ONE_VAR_S_ALLREC(filename,var_to_read,nb_records)

        USE global_constants_mod, only: ip
        USE file_libs, only: get_fID, release_fID

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER(KIND=ip), INTENT(out) :: nb_records
        REAL(KIND=silp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: var_to_read

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER(kind=sip) :: i,j,k
        INTEGER(kind=ip)  :: file_ID, file_size
        INTEGER(kind=ip), PARAMETER  :: size_rec = size_X*size_Y*1*silp
!~         REAL(KIND=silp), DIMENSION(size_Y,size_X) :: local_var_read

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        file_ID = get_fID()
        
        INQUIRE(FILE=filename, SIZE=file_size)

        nb_records = file_size/size_rec

        allocate(var_to_read(size_Y,size_X,nb_records))
        
        OPEN(file_ID,file=filename, access='direct', form='unformatted', status='old',recl=size_rec)

        DO k = 1, nb_records
           READ(file_ID,REC=k) ((var_to_read(i,j,k),j=1,size_X),i=1,size_Y)
!~            var_to_read(:,:,k) = REAL(local_var_read(:,:),kind=dblp)
        ENDDO

        CLOSE(file_ID)

        call release_fID(file_ID)

      END SUBROUTINE READ_ONE_VAR_S_ALLREC


! ---


      SUBROUTINE WRITE_CTL_HEADER(ctrlfilename, binfilename, varname, freq, dtype, undef)

        USE global_constants_mod, ONLY: ip, str_len
        USE file_libs,            ONLY: get_fID, release_fID
        USE comemic_mod,          ONLY: iyear
        USE global_constants_mod, ONLY: months_year_i

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        CHARACTER(LEN=*),INTENT(IN) :: ctrlfilename, binfilename, varname
        CHARACTER,       INTENT(IN) :: freq,dtype
        REAL(KIND=dblp), INTENT(IN), OPTIONAL :: undef
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER(kind=sip)     :: i,j
        INTEGER(kind=ip)      :: file_ID
        CHARACTER(LEN=str_len):: timing_string


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        file_ID = get_fID()

        OPEN(unit=file_ID,file=ctrlfilename, form='formatted', action='write', status='unknown')

        WRITE(unit=file_ID,fmt="(A)") "DSET ^"//binfilename
        WRITE(unit=file_ID,fmt="(A)") "TITLE "//"Climate indices derived from run ... of iLOVECLIM"
        if ( dtype == "i" ) then
          WRITE(unit=file_ID,fmt="(A)") "UNDEF 0"
        else
          IF (PRESENT(undef)) THEN
             WRITE(unit=file_ID,fmt="(A,F10.1)") "UNDEF ",undef
          ELSE
             WRITE(unit=file_ID,fmt="(A)") "UNDEF 9999.0"
          ENDIF
        endif
        WRITE(unit=file_ID,fmt="(A)") "OPTIONS big_endian"
        WRITE(unit=file_ID,fmt="(A)") "XDEF      64  LINEAR     0.000     5.625"
        WRITE(unit=file_ID,fmt="(A)") "YDEF      32  LEVELS   -85.761   -80.269   -74.745   -69.213   -63.679   -58.143"
        WRITE(unit=file_ID,fmt="(A)") "-52.607   -47.070   -41.532   -35.995   -30.458   -24.920   -19.382   -13.844"
        WRITE(unit=file_ID,fmt="(A)") "-8.307    -2.769     2.769     8.307    13.844    19.382    24.920    30.458"
        WRITE(unit=file_ID,fmt="(A)") "35.995    41.532    47.070    52.607    58.143    63.679    69.213    74.745"
        WRITE(unit=file_ID,fmt="(A)") "80.269    85.761"
        WRITE(unit=file_ID,fmt="(A)") "ZDEF       1  LEVELS 1 1"
        if ( freq == "m" ) then
           WRITE(timing_string,fmt="(A,I6,A)") "TDEF   ", iyear*months_year_i ,"  LINEAR  0Z15jan0001  1mo"
        else
           WRITE(timing_string,fmt="(A,I6,A)") "TDEF   ", iyear ,"  LINEAR  0Z1jul0001  1yr"
        endif
        WRITE(unit=file_ID,fmt="(A)") TRIM(timing_string)
        WRITE(unit=file_ID,fmt="(A)") "VARS 1"
        if ( dtype == "i" ) then
          WRITE(unit=file_ID,fmt="(A)") varname//"  1 -1,40,2,-1  "//" blaaaa"
        else
          WRITE(unit=file_ID,fmt="(A)") varname//"  1 99          "//" blaaaa"
        endif
        WRITE(unit=file_ID,fmt="(A)") "ENDVARS"
        CLOSE(file_ID)

        call release_fID(file_ID)

      END SUBROUTINE WRITE_CTL_HEADER


! ---

#if ( CLM_INDICES >= 2  && BIOM_GEN == 1)

      SUBROUTINE COMPUTE_AND_WRITE_BIOMES()
      
      USE global_constants_mod, ONLY: ip, eps_dp
      USE DIAGNOSE_BIOME_MOD,   ONLY: determine_biome
      USE input_icemask,        ONLY: icemask
      USE comsurf_mod,          ONLY: fractn, nld
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         REAL(KIND=silp), DIMENSION(:,:,:), ALLOCATABLE :: tree_fraction, needle_fraction, grass_fraction, desert_fraction   &
                                                         , gdd0_values, gdd5_values, tmin_values, tmax_values, biomes_values &
                                                         , imask_fraction
                                                         
         INTEGER(kind=ip) :: nbrecs
         INTEGER(kind=sip):: i,j,k
         
         CHARACTER(LEN=12), PARAMETER :: biome_filenm="BIOM_indices"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---    READ THE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         CALL READ_ONE_VAR_S_ALLREC(binaary_filenames(TmiN_y_key)        ,tmin_values     ,nbrecs)
         CALL READ_ONE_VAR_S_ALLREC(binaary_filenames(TmaX_y_key)        ,tmax_values     ,nbrecs)
         CALL READ_ONE_VAR_S_ALLREC(binaary_filenames(treefrac_key)      ,tree_fraction   ,nbrecs)
         CALL READ_ONE_VAR_S_ALLREC(binaary_filenames(needletreefrac_key),needle_fraction ,nbrecs)
         CALL READ_ONE_VAR_S_ALLREC(binaary_filenames(grassfrac_key)     ,grass_fraction  ,nbrecs)
         CALL READ_ONE_VAR_S_ALLREC(binaary_filenames(desertfrac_key)    ,desert_fraction ,nbrecs)
         CALL READ_ONE_VAR_S_ALLREC(binaary_filenames(GDD_zero_key)      ,gdd0_values     ,nbrecs)
         CALL READ_ONE_VAR_S_ALLREC(binaary_filenames(GDD_five_key)      ,gdd5_values     ,nbrecs) 


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---    CREATE MISSING VARIABLES AND CHANGE UNITS
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


         ALLOCATE(imask_fraction(size_Y, size_X, nbrecs)) ! for now undefined ...
         
         DO k=LBOUND(imask_fraction,DIM=3),UBOUND(imask_fraction,DIM=3)
           DO j=LBOUND(imask_fraction,DIM=2),UBOUND(imask_fraction,DIM=2)
             DO i=LBOUND(imask_fraction,DIM=1),UBOUND(imask_fraction,DIM=1)
                IF (fractn(i,j,nld).gt.eps_dp) then
                  imask_fraction(i,j,k) = icemask(i,j)*100._silp
                ELSE
                  imask_fraction(i,j,k) = big_dp * (-1._silp)
                ENDIF
             ENDDO
           ENDDO
         ENDDO
         
         
         tmin_values(:,:,:) = tmin_values(:,:,:) - freezeT
         tmax_values(:,:,:) = tmax_values(:,:,:) - freezeT
         tree_fraction(:,:,:) = tree_fraction(:,:,:) * 100._silp
         grass_fraction(:,:,:) = grass_fraction(:,:,:) * 100._silp
         needle_fraction(:,:,:) = needle_fraction(:,:,:) * 100._silp
         desert_fraction(:,:,:) = desert_fraction(:,:,:) * 100._silp
         
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---    CREATE BIOME ARRAY AND COMPUTE BIOMES VALUES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         
         ALLOCATE(biomes_values(size_Y, size_X, nbrecs))

         CALL determine_biome(tmin_values, tmax_values, tree_fraction, needle_fraction, grass_fraction, desert_fraction       &
                            , gdd0_values, gdd5_values, imask_fraction, size_Y, size_X, nbrecs, nbrecs, biomes_values)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---    WRITE OUT BIOME RESULTS
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         CALL WRITE_ONE_VAR_S_ALLREC(""//directoire_stuff//biome_filenm//binaary_suff,biomes_values)
         CALL WRITE_CTL_HEADER(""//directoire_stuff//biome_filenm//control_suff, biome_filenm//binaary_suff                   &
                               , "biom",'y','f', big_dp*(-1._dblp))

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---    GARBAGE COLLECTION
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
   
         DEALLOCATE(tree_fraction)
         DEALLOCATE(needle_fraction)
         DEALLOCATE(grass_fraction)
         DEALLOCATE(desert_fraction)
         DEALLOCATE(gdd0_values)
         DEALLOCATE(gdd5_values)
         DEALLOCATE(tmin_values)
         DEALLOCATE(tmax_values)
         DEALLOCATE(imask_fraction)
         DEALLOCATE(biomes_values)


      END SUBROUTINE COMPUTE_AND_WRITE_BIOMES

#endif


      END MODULE CLIMATE_INDICES_MOD
#endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
