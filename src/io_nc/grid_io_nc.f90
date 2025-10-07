!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2020-2021 Didier M. Roche (a.k.a. dmr)
!   Copyright 2021      Pepijn Bakker   (a.k.a. PB)

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#include "choixcomposantes.h"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( CLIO_OUT_NEWGEN == 1 )

      MODULE GRID_IO_NC

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       USE global_constants_mod, ONLY: dblp=>dp, silp=>sp, sip,  freezeT => tK_zero_C, big_dp => alt_olympus_mons, ip
      
       USE IO_NC_MOD, ONLY: IO_NC_FILE, IO_GRID_VAR, IO_NC_AXIS
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      IMPLICIT NONE
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   History
! dmr           Change from 0.0.0: Created a first version that can test write a CLIO file
! PB            Change from 0.1.0: Added a first support for the iceberg variables
! dmr           Change from 0.2.0: Added a first support 3D ocean files (with time)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
      CHARACTER(LEN=5), PARAMETER :: version_mod ="0.3.0"
      

      PRIVATE
      PUBLIC :: DAILYSTEP_IO_NC, GLOBAL_RE_INIT, GLOBAL_FINALIZE


        REAL(KIND=dblp), PARAMETER :: bigN_dp = -1.0*big_dp, d_zero = 0.0_dblp

!--- UGLY HACK FOR CLIO GENERATION TO BE REPLACED !!!

        INTEGER(KIND=ip), PARAMETER :: size_x = 120, size_y = 65, size_class = 10, size_tdepth = 20      

        REAL(KIND=silp), DIMENSION(size_y) :: y_axis_values =                                                                     &
                 (/  -79.5, -76.5, -73.5, -70.5, -67.5, -64.5, -61.5, -58.5, -55.5, -52.5, -49.5, -46.5, -43.5, -40.5, -37.5,     &
                     -34.5, -31.5, -28.5, -25.5, -22.5, -19.5, -16.5, -13.5, -10.5, -7.5, -4.5, -1.5, 1.5, 4.5, 7.5, 10.5, 13.5,  &
                      16.5, 19.5, 22.5, 25.5, 28.5, 31.5, 34.5, 37.5, 40.5, 43.5, 46.5, 49.5, 52.5, 55.5, 58.5, 61.5, 64.5, 67.5, &
                      70.5, 73.5, 76.5, 79.5, 82.5, 85.5, 88.5, 91.5, 94.5, 97.5, 100.5, 103.5, 106.5, 109.5, 112.5 /)

        REAL(KIND=silp), DIMENSION(size_x) :: x_axis_values =                                                                     &
                   (/ 28.5, 31.5, 34.5, 37.5, 40.5, 43.5, 46.5, 49.5, 52.5, 55.5, 58.5, 61.5, 64.5, 67.5, 70.5, 73.5, 76.5, 79.5, &
                      82.5, 85.5, 88.5, 91.5, 94.5, 97.5, 100.5, 103.5, 106.5, 109.5, 112.5, 115.5, 118.5, 121.5, 124.5, 127.5,   &
                      130.5, 133.5, 136.5, 139.5, 142.5, 145.5, 148.5, 151.5, 154.5, 157.5, 160.5, 163.5, 166.5, 169.5, 172.5,    &
                      175.5, 178.5, 181.5, 184.5, 187.5, 190.5, 193.5, 196.5, 199.5, 202.5, 205.5, 208.5, 211.5, 214.5, 217.5,    &
                      220.5, 223.5, 226.5, 229.5, 232.5, 235.5, 238.5, 241.5, 244.5, 247.5, 250.5, 253.5, 256.5, 259.5, 262.5,    &
                      265.5, 268.5, 271.5, 274.5, 277.5, 280.5, 283.5, 286.5, 289.5, 292.5, 295.5, 298.5, 301.5, 304.5, 307.5,    &
                      310.5, 313.5, 316.5, 319.5, 322.5, 325.5, 328.5, 331.5, 334.5, 337.5, 340.5, 343.5, 346.5, 349.5, 352.5,    &
                      355.5, 358.5, 361.5, 364.5, 367.5, 370.5, 373.5, 376.5, 379.5, 382.5, 385.5 /)

        REAL(KIND=silp), DIMENSION(size_tdepth) :: tdepth_axis_values =                                                           &
                   (/ -5126.18, -4385.22, -3661.11, -2963.25, -2307.36, -1717.9, -1225.11, -850.19, -588.88, -415.07, -299.29,    &
                      -219.86, -163.28, -121.52, -89.75, -64.96, -45.2, -29.17, -15.98, -5.0 /)

! axis for the iceberg 'classes'
        INTEGER(KIND=ip), DIMENSION(size_class) :: class_axis_values =                                                            &
                   (/  33, 100, 166, 233, 300, 366, 450, 550, 700, 900 /)



!--- UGLY HACK FOR CLIO GENERATION TO BE REPLACED !!!



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       DEFINITION OF OUTPUT VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        TYPE(IO_NC_FILE), TARGET                :: test_CLIO_file
        TYPE(IO_NC_AXIS), TARGET                :: CLIO_axis_ptlat, CLIO_axis_ptlon, CLIO_axis_time, CLIO_axis_class,&
                                                   CLIO_axis_tdepth
 
!---    CLIO variables
        TYPE(IO_GRID_VAR)                         :: mu_SST
        REAL(dblp), DIMENSION(:,:), ALLOCATABLE   :: my_SST

        TYPE(IO_GRID_VAR)                         :: mu_fwruno
        REAL(dblp), DIMENSION(:,:), ALLOCATABLE   :: my_fwruno

        TYPE(IO_GRID_VAR)                         :: mu_3D_temp
        REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE :: my_3D_temp

        TYPE(IO_GRID_VAR)                         :: mu_fsolcn
        REAL(dblp), DIMENSION(:,:), ALLOCATABLE   :: my_fsolcn


#if ( REMIN == 1 )
        TYPE(IO_GRID_VAR)                         :: mu_3D_remin 
        REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE :: my_3D_remin
#endif

#if ( NEOD > 0 )
        TYPE(IO_GRID_VAR)                       :: mu_neodymium143
        REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE :: my_neodymium143

        TYPE(IO_GRID_VAR)                       :: mu_neodymium144
        REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE :: my_neodymium144

        TYPE(IO_GRID_VAR)                       :: mu_epsNd
        REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE :: my_epsNd
#endif 

#if ( OCYCC == 1 )
        TYPE(IO_GRID_VAR)                           :: mu_fPOC_clio
        REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE   :: my_fPOC_clio

        TYPE(IO_GRID_VAR)                           :: mu_fCAL_clio
        REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE   :: my_fCAL_clio
#endif

! eclermont - added
#if ( OOISO == 1 )
        TYPE(IO_GRID_VAR)                         :: mu_phyto_clio 
        REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE :: my_phyto_clio
        
        TYPE(IO_GRID_VAR)                         :: mu_respO2
        REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE :: my_respO2
        
        TYPE(IO_GRID_VAR)                         :: mu_prodO2
        REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE :: my_prodO2
        
        TYPE(IO_GRID_VAR)                       :: mu_reminO2
        REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE :: my_reminO2
        
        TYPE(IO_GRID_VAR)                       :: mu_NCP
        REAL(dblp), DIMENSION(:,:,:), ALLOCATABLE :: my_NCP

        TYPE(IO_GRID_VAR)                       :: mu_FtotalO
        REAL(dblp), DIMENSION(:,:), ALLOCATABLE :: my_FtotalO
        
        TYPE(IO_GRID_VAR)                       :: mu_F17O
        REAL(dblp), DIMENSION(:,:), ALLOCATABLE :: my_F17O
        
        TYPE(IO_GRID_VAR)                       :: mu_F18O
        REAL(dblp), DIMENSION(:,:), ALLOCATABLE :: my_F18O
#endif

#if ( ICEBERG > 0 )
!---    ICEBERG variables
        TYPE(IO_GRID_VAR)                       :: mu_vol_icb
        REAL(dblp), DIMENSION(:,:), ALLOCATABLE :: my_vol_icb
        TYPE(IO_GRID_VAR)                       :: mu_heat_icb
        REAL(dblp), DIMENSION(:,:), ALLOCATABLE :: my_heat_icb
        TYPE(IO_GRID_VAR)                       :: mu_dvol_icb
        REAL(dblp), DIMENSION(:,:), ALLOCATABLE :: my_dvol_icb        
        TYPE(IO_GRID_VAR)                       :: mu_hiceb_class_icb
        REAL(dblp), DIMENSION(:,:,:),ALLOCATABLE:: my_hiceb_class_icb  
        TYPE(IO_GRID_VAR)                       :: mu_wiceb_class_icb
        REAL(dblp), DIMENSION(:,:,:),ALLOCATABLE:: my_wiceb_class_icb 
        TYPE(IO_GRID_VAR)                       :: mu_pond_class_icb
        REAL(dblp), DIMENSION(:,:,:),ALLOCATABLE:: my_pond_class_icb         
        TYPE(IO_GRID_VAR)                       :: mu_uiceb_class_icb
        REAL(dblp), DIMENSION(:,:,:),ALLOCATABLE:: my_uiceb_class_icb  
        TYPE(IO_GRID_VAR)                       :: mu_viceb_class_icb
        REAL(dblp), DIMENSION(:,:,:),ALLOCATABLE:: my_viceb_class_icb                                      
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       FILE NAMES DEFINITIONS AND HANDLING
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      CONTAINS

! ---

      SUBROUTINE GLOBAL_RE_INIT()


        use para0_mod, only: imax, jmax, kmax
        use iceberg_mod, only: numclass        
        
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

        ! initialize netCDF files, including axes attributes and initial coordinates variables
        ! for now, do a relatively ugly direct import. Will do better later on.

        call test_CLIO_file%init(FileName="CLIO3_NewGen.nc")
        call CLIO_axis_ptlat%init("ptlat", y_axis_values, OPT_AxisUnit="degrees_north")
        call CLIO_axis_ptlon%init("ptlon", x_axis_values, OPT_AxisUnit="degrees_east")
        call CLIO_axis_class%init("class", class_axis_values, OPT_AxisUnit="iceberg_size_class")
        call CLIO_axis_tdepth%init("tdepth", tdepth_axis_values, OPT_AxisUnit="meters")
        call CLIO_axis_time%init("time",OPT_isTime=.true., OPT_AxisUnit="years since 0001-01-01", OPT_calendar="360_day")
        
        call CLIO_axis_ptlat%wrte("CLIO3_NewGen.nc")
        call CLIO_axis_ptlon%wrte("CLIO3_NewGen.nc")
        call CLIO_axis_class%wrte("CLIO3_NewGen.nc")
        call CLIO_axis_tdepth%wrte("CLIO3_NewGen.nc")
        call CLIO_axis_time%wrte("CLIO3_NewGen.nc")     
        
        allocate(my_SST(imax-2,jmax))
        allocate(my_fwruno(imax-2,jmax))
        allocate(my_3D_temp(imax-2,jmax,kmax))

#if ( REMIN == 1 )
        allocate(my_3D_remin(imax-2,jmax,kmax))
#endif

#if ( NEOD > 0 ) 
        allocate(my_neodymium143(imax-2,jmax,kmax))
        allocate(my_neodymium144(imax-2,jmax,kmax))
        allocate(my_epsNd(imax-2,jmax,kmax))
#endif

#if ( OCYCC == 1 )
         allocate(my_fPOC_clio(imax-2,jmax,kmax))
         allocate(my_fCAL_clio(imax-2,jmax,kmax))
#endif

! ajout par eclermont
#if ( OOISO == 1 )
        allocate(my_phyto_clio(imax-2,jmax,kmax))
        allocate(my_respO2(imax-2,jmax,kmax))
        allocate(my_prodO2(imax-2,jmax,kmax))
        allocate(my_reminO2(imax-2,jmax,kmax))
        allocate(my_NCP(imax-2,jmax,kmax))
        allocate(my_FtotalO(imax-2,jmax))
        allocate(my_F17O(imax-2,jmax))
        allocate(my_F18O(imax-2,jmax))
#endif

#if ( ICEBERG > 0 )        
        allocate(my_vol_icb(imax-2,jmax))
        allocate(my_heat_icb(imax-2,jmax))
        allocate(my_dvol_icb(imax-2,jmax))        
        allocate(my_hiceb_class_icb(imax-2,jmax,numclass))        
        allocate(my_wiceb_class_icb(imax-2,jmax,numclass))        
        allocate(my_pond_class_icb(imax-2,jmax,numclass))        
        allocate(my_uiceb_class_icb(imax-2,jmax,numclass))        
        allocate(my_viceb_class_icb(imax-2,jmax,numclass))              
        
#endif        
        
        call mu_SST%init("SST",test_CLIO_file,"ptlon ptlat time")
        call mu_SST%show()

        call mu_fwruno%init("fwruno",test_CLIO_file,"ptlon ptlat time")
        call mu_fwruno%show()

        call mu_3D_temp%init("temp",test_CLIO_file,"ptlon ptlat tdepth time")
        call mu_3D_temp%show()

#if ( REMIN == 1 )
        call mu_3D_remin%init("remin",test_CLIO_file,"ptlon ptlat tdepth time")
        call mu_3D_remin%show()
#endif

#if ( NEOD > 0 ) 
        call mu_neodymium143%init("neodymium143",test_CLIO_file,"ptlon ptlat tdepth time")
        call mu_neodymium143%show()
        call mu_neodymium144%init("neodymium144",test_CLIO_file,"ptlon ptlat tdepth time")
        call mu_neodymium144%show()     
        call mu_epsNd%init("epsNd",test_CLIO_file,"ptlon ptlat tdepth time")
        call mu_epsNd%show()        
#endif   

#if ( OCYCC == 1 )
        call mu_fPOC_clio%init("fPOC",test_CLIO_file,"ptlon ptlat tdepth time")
        call mu_fPOC_clio%show()

        call mu_fCAL_clio%init("fCAL",test_CLIO_file,"ptlon ptlat tdepth time")
        call mu_fCAL_clio%show()
#endif

! Ajout par eclermont 
#if ( OOISO == 1 )
        call mu_phyto_clio%init("phyto_growth",test_CLIO_file,"ptlon ptlat tdepth time")
        call mu_phyto_clio%show()

        call mu_respO2%init("resp_O2",test_CLIO_file,"ptlon ptlat tdepth time")
        call mu_respO2%show()
                
        call mu_prodO2%init("prod_O2",test_CLIO_file,"ptlon ptlat tdepth time")
        call mu_prodO2%show()
        
        call mu_reminO2%init("Remin",test_CLIO_file,"ptlon ptlat tdepth time")
        call mu_reminO2%show()

        call mu_NCP%init("NCP",test_CLIO_file,"ptlon ptlat tdepth time")
        call mu_NCP%show()

        call mu_FtotalO%init("Ototal_flux",test_CLIO_file,"ptlon ptlat time")
        call mu_FtotalO%show()

        call mu_F17O%init("O17_flux",test_CLIO_file,"ptlon ptlat time")
        call mu_F17O%show()

        call mu_F18O%init("O18_flux",test_CLIO_file,"ptlon ptlat time")
        call mu_F18O%show()
#endif

#if ( ICEBERG > 0 )           
        call mu_vol_icb%init("vol_icb",test_CLIO_file,"ptlon ptlat time")
        call mu_vol_icb%show()
        call mu_heat_icb%init("heat_icb",test_CLIO_file,"ptlon ptlat time")
        call mu_heat_icb%show()
        call mu_dvol_icb%init("dvol_icb",test_CLIO_file,"ptlon ptlat time")
        call mu_dvol_icb%show()                
        call mu_hiceb_class_icb%init("hiceb_class",test_CLIO_file,"ptlon ptlat class time")
        call mu_hiceb_class_icb%show() 
        call mu_wiceb_class_icb%init("wiceb_class",test_CLIO_file,"ptlon ptlat class time")
        call mu_wiceb_class_icb%show() 
        call mu_pond_class_icb%init("pond_class",test_CLIO_file,"ptlon ptlat class time")
        call mu_pond_class_icb%show() 
        call mu_uiceb_class_icb%init("uiceb_class",test_CLIO_file,"ptlon ptlat class time")
        call mu_uiceb_class_icb%show() 
        call mu_viceb_class_icb%init("viceb_class",test_CLIO_file,"ptlon ptlat class time")
        call mu_viceb_class_icb%show() 
#endif        
        
        ! re-initialize monthly accumulation variables
        CALL MONTHLY_RE_INIT()

        ! re-initialize annual accumulation variables
        CALL YEARLY_RE_INIT()




      END SUBROUTINE GLOBAL_RE_INIT

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

        ! write summary netCDF IO file?

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

!--- Temporary 

        my_SST(:,:) = d_zero
        my_fwruno(:,:) = d_zero        
        my_3D_temp(:,:,:) = d_zero

#if ( REMIN == 1 )
        my_3D_remin(:,:,:) = d_zero
#endif

#if ( NEOD > 0 )
        my_neodymium143(:,:,:) = d_zero
        my_neodymium144(:,:,:) = d_zero
        my_epsNd(:,:,:) = d_zero
#endif  

#if ( OCYCC == 1 )
        my_fPOC_clio(:,:,:) = d_zero
        my_fCAL_clio(:,:,:) = d_zero
#endif        
        
! ajout par eclermont
#if ( OOISO == 1 )
        my_phyto_clio(:,:,:) = d_zero
        my_respO2(:,:,:) = d_zero
        my_prodO2(:,:,:) = d_zero
        my_reminO2(:,:,:) = d_zero
        my_NCP(:,:,:) = d_zero
        my_FtotalO(:,:) = d_zero
        my_F17O(:,:) = d_zero
        my_F18O(:,:) = d_zero
#endif 

#if ( ICEBERG > 0 )         
        my_vol_icb(:,:) = d_zero
        my_heat_icb(:,:) = d_zero
        my_dvol_icb(:,:) = d_zero
        my_hiceb_class_icb(:,:,:) = d_zero
        my_wiceb_class_icb(:,:,:) = d_zero
        my_pond_class_icb(:,:,:) = d_zero
        my_uiceb_class_icb(:,:,:) = d_zero
        my_viceb_class_icb(:,:,:) = d_zero
#endif        
        

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


      END SUBROUTINE MONTHLY_RE_INIT


      SUBROUTINE DAILYSTEP_IO_NC()

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     NOTA: currently called every day as name indicates
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        use comemic_mod, only: iday, imonth, iyear
        use global_constants_mod, only: months_year_i, days_month360d_i, ip, days_year360d


!--- Temporary 
        use bloc0_mod, only: ks2, scal
        use ice_mod, only: fwruno
        
#if ( REMIN == 1 )
        use bloc0_mod, only: kremin_clio
#endif
        
#if ( ICEBERG > 0 )                 
        use iceberg_mod, only: vol_icb, dVol_icb, heat_icb, hiceb_class, wiceb_class, pond_icb_class, uiceb_class, viceb_class
#endif        

#if ( OCYCC == 1 )
        use bloc0_mod, only: fPOC_flx_clio, fCAL_flx_clio
#endif

! Ajout par eclermont
#if ( OOISO == 1 )
        use bloc0_mod, only: phyto_clio, respO2_clio, prodO2_clio, reminO2_clio, NCP_clio
        use marine_bio_mod, only: FOO2
#endif
                   
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER(kind=ip) :: i,j, c_month, c_year

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!--- Temporary, mean of the variable along timesteps

!---     CLIO SST
         my_SST(:,:) = my_SST(:,:) + scal(2:UBOUND(scal,dim=1)-1,:,ks2,1)/days_year360d ! tryout with SST
         my_fwruno(:,:) = my_fwruno(:,:) + fwruno(2:UBOUND(fwruno,dim=1)-1,:)/days_year360d ! tryout with SST
         my_3D_temp(:,:,:) = my_3D_temp(:,:,:) + scal(2:UBOUND(scal,dim=1)-1,:,:,1)/days_year360d ! tryout with SST

#if ( REMIN == 1 )
         my_3D_remin(:,:,:) = my_3D_remin(:,:,:) + kremin_clio(2:UBOUND(kremin_clio,dim=1)-1,:,:)/days_year360d
#endif

#if ( NEOD > 0 )
!---     CLIO neodymium
         my_neodymium143(:,:,:) = my_neodymium143(:,:,:) + scal(2:UBOUND(scal,dim=1)-1,:,:,18)/days_year360d
         my_neodymium144(:,:,:) = my_neodymium144(:,:,:) + scal(2:UBOUND(scal,dim=1)-1,:,:,19)/days_year360d
#endif 

#if ( ICEBERG > 0 )   
!---     ICEBERGS total iceberg volume
         my_vol_icb(:,:) = my_vol_icb(:,:) + vol_icb/days_year360d
!---     ICEBERGS total iceberg heat to clio
         my_heat_icb(:,:) = my_heat_icb(:,:) + sum(heat_icb,DIM=3)/days_year360d
!---     ICEBERGS melting to clio
         my_dvol_icb(:,:) = my_dvol_icb(:,:) + sum(dVol_icb,DIM=3)/days_year360d
!---     ICEBERGS mean iceberg height per class to clio
         my_hiceb_class_icb(:,:,:) = my_hiceb_class_icb(:,:,:) + hiceb_class/days_year360d
!---     ICEBERGS mean iceberg width per class to clio
         my_wiceb_class_icb(:,:,:) = my_wiceb_class_icb(:,:,:) + wiceb_class/days_year360d
!---     ICEBERGS mean number of iceberg per class to clio
         my_pond_class_icb(:,:,:) = my_pond_class_icb(:,:,:) + pond_icb_class/days_year360d
!---     ICEBERGS mean iceberg zonal velocity per class to clio
         my_uiceb_class_icb(:,:,:) = my_uiceb_class_icb(:,:,:) + uiceb_class/days_year360d
!---     ICEBERGS mean iceberg meridional velocity per class to clio
         my_viceb_class_icb(:,:,:) = my_viceb_class_icb(:,:,:) + viceb_class/days_year360d

#endif         


#if ( OCYCC == 1 )
         my_fPOC_clio(:,:,:) = my_fPOC_clio(:,:,:) + fPOC_flx_clio(:,:,:)/days_year360d
         my_fCAL_clio(:,:,:) = my_fCAL_clio(:,:,:) + fCAL_flx_clio(:,:,:)/days_year360d
#endif        

! Ajout par eclermont
#if ( OOISO == 1)
         my_phyto_clio(:,:,:) = my_phyto_clio(:,:,:) + phyto_clio(2:UBOUND(phyto_clio,dim=1)-1,:,:)/days_year360d
         my_respO2(:,:,:)=my_respO2(:,:,:) + respO2_clio(2:UBOUND(respO2_clio,dim=1)-1,:,:)/days_year360d
         my_prodO2(:,:,:) = my_prodO2(:,:,:) + prodO2_clio(2:UBOUND(prodO2_clio,dim=1)-1,:,:)/days_year360d
         my_reminO2(:,:,:) = my_reminO2(:,:,:) + reminO2_clio(2:UBOUND(reminO2_clio,dim=1)-1,:,:)/days_year360d
         my_NCP(:,:,:) = my_NCP(:,:,:) + NCP_clio(2:UBOUND(NCP_clio,dim=1)-1,:,:)/days_year360d

         my_FtotalO(:,:) = my_FtotalO(:,:) + TRANSPOSE(FOO2(2:UBOUND(FOO2,dim=1)-1,:,1))
         my_F18O(:,:) = my_F18O(:,:) + TRANSPOSE(FOO2(2:UBOUND(FOO2,dim=1)-1,:,4))
         my_F17O(:,:) = my_F17O(:,:) + TRANSPOSE(FOO2(2:UBOUND(FOO2,dim=1)-1,:,3))
#endif   

         if (iday == days_month360d_i) then
           c_month = (iyear-1)*months_year_i+imonth
           CALL DO_END_MONTH(c_month)
         endif

         if (imonth == months_year_i .and. iday == days_month360d_i) then
           c_year = iyear
           CALL DO_END_YEAR(c_year)
         endif

      END SUBROUTINE DAILYSTEP_IO_NC


! ---

      SUBROUTINE DO_END_MONTH(nb_month)

         INTEGER(KIND=ip), INTENT(IN) :: nb_month

         INTEGER(KIND=ip) :: i,j

         ! write out the monthly fields

         CALL MONTHLY_RE_INIT()

      END SUBROUTINE DO_END_MONTH

! ---

      SUBROUTINE DO_END_YEAR(nb_year)

!--- Temporary 
         use bloc0_mod, only: ks2, tms, kmax
         
         use IO_NC_MOD, only: undef_dblp

#if ( NEOD > 0 )
         use neodymium_mod, only: epsilonNd ! to compute the epsilonNd from Nd144 and Nd143
#endif

         INTEGER(KIND=ip), INTENT(IN) :: nb_year
         INTEGER(KIND=ip) classes, k

         ! write out the yearly fields
         
!--- Temporary        
         WHERE(tms(2:UBOUND(tms,dim=1)-1,:,ks2).LT.EPSILON(tms(1,1,ks2)))

            my_SST(:,:) = undef_dblp
            my_fwruno(:,:) = undef_dblp
            
#if ( OOISO == 1 )
            my_FtotalO(:,:)=undef_dblp
            my_F17O(:,:)=undef_dblp
            my_F18O(:,:)=undef_dblp
#endif 

#if ( ICEBERG > 0 )              
            my_vol_icb(:,:) = undef_dblp
            my_heat_icb(:,:) = undef_dblp
            my_dvol_icb(:,:) = undef_dblp
#endif           
            
         ENDWHERE

         DO k=1,kmax
           WHERE(tms(2:UBOUND(tms,dim=1)-1,:,k).LT.EPSILON(tms(1,1,k)))
            my_3D_temp(:,:,k) = undef_dblp
#if ( REMIN == 1 )
            my_3D_remin(:,:,k) = undef_dblp
#endif

! Ajouter par eclermont : 
#if ( OOISO == 1 )
            my_phyto_clio(:,:,k) = undef_dblp
            my_respO2(:,:,k) = undef_dblp
            my_prodO2(:,:,k) =undef_dblp
            my_reminO2(:,:,k) =undef_dblp
            my_NCP(:,:,k) =undef_dblp
#endif

           ENDWHERE
         ENDDO


#if ( ICEBERG > 0 )              
         do classes=LBOUND(my_hiceb_class_icb,dim=3), UBOUND(my_hiceb_class_icb,dim=3)
            WHERE(tms(2:UBOUND(tms,dim=1)-1,:,ks2).LT.EPSILON(tms(1,1,ks2)))
                my_hiceb_class_icb(:,:,classes) = undef_dblp
                my_wiceb_class_icb(:,:,classes) = undef_dblp
                my_pond_class_icb(:,:,classes) = undef_dblp
                my_uiceb_class_icb(:,:,classes) = undef_dblp
                my_viceb_class_icb(:,:,classes) = undef_dblp
            ENDWHERE
         enddo            
#endif            

#if ( NEOD > 0 ) 
         my_epsNd(:,:,:) = epsilonNd(my_neodymium143(:,:,:),my_neodymium144(:,:,:))             
         DO k=1,kmax
            WHERE(tms(2:UBOUND(tms,dim=1)-1,:,k).LT.EPSILON(tms(1,1,k)))
                my_neodymium143(:,:,k) = undef_dblp
                my_neodymium144(:,:,k) = undef_dblp
                my_epsNd(:,:,k) = undef_dblp
            ENDWHERE
         enddo            
#endif            

         call mu_SST%wrte(my_SST,nb_year)
         call mu_fwruno%wrte(my_fwruno,nb_year)
         call mu_3D_temp%wrte(my_3D_temp,nb_year)

#if ( REMIN == 1 )
         call mu_3D_remin%wrte(my_3D_remin,nb_year)
#endif

#if ( NEOD > 0 ) 
         call mu_neodymium143%wrte(my_neodymium143,nb_year)
         call mu_neodymium144%wrte(my_neodymium144,nb_year)
         call mu_epsNd%wrte(my_epsNd,nb_year)         
#endif 

#if ( ICEBERG > 0 )                       
         call mu_vol_icb%wrte(my_vol_icb,nb_year)
         call mu_heat_icb%wrte(my_heat_icb,nb_year)
         call mu_dvol_icb%wrte(my_dvol_icb,nb_year)
         call mu_hiceb_class_icb%wrte(my_hiceb_class_icb,nb_year)
         call mu_wiceb_class_icb%wrte(my_wiceb_class_icb,nb_year)
         call mu_pond_class_icb%wrte(my_pond_class_icb,nb_year)
         call mu_uiceb_class_icb%wrte(my_uiceb_class_icb,nb_year)
         call mu_viceb_class_icb%wrte(my_viceb_class_icb,nb_year)
         
#endif         

#if ( OCYCC == 1 ) 
         call mu_fCAL_clio%wrte(my_fCAL_clio,nb_year)
         call mu_fPOC_clio%wrte(my_fPOC_clio,nb_year)
#endif

! Ajout par eclermont
#if ( OOISO == 1)
         call mu_phyto_clio%wrte(my_phyto_clio,nb_year)
         call mu_respO2%wrte(my_respO2,nb_year)
         call mu_prodO2%wrte(my_prodO2,nb_year)
         call mu_FtotalO%wrte(my_FtotalO, nb_year)
         call mu_F17O%wrte(my_F17O, nb_year)
         call mu_F18O%wrte(my_F18O, nb_year)
         call mu_reminO2%wrte(my_reminO2, nb_year)
         call mu_NCP%wrte(my_NCP, nb_year)
#endif

         CALL YEARLY_RE_INIT()

      END SUBROUTINE DO_END_YEAR

! ---

      END MODULE GRID_IO_NC
#endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
