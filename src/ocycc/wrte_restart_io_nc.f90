!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Copyright 2020 Didier M. Roche (a.k.a. dmr)

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

#if ( CLIO_OUT_NEWGEN == 1 && OCYCC == 1 )

      MODULE WRTE_RESTART_IO_NC

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       USE global_constants_mod, ONLY: dblp=>dp, silp=>sp, sip, big_dp => alt_olympus_mons, ip

       USE IO_NC_MOD, ONLY: IO_NC_FILE, IO_GRID_VAR, IO_NC_AXIS


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: WRTE_RESTART_OCYCC


        REAL(KIND=dblp), PARAMETER :: bigN_dp = -1.0*big_dp, d_zero = 0.0_dblp

!--- UGLY HACK FOR CLIO GENERATION TO BE REPLACED !!!

        INTEGER(KIND=ip), PARAMETER :: size_x = 120, size_y = 65, size_z = 20

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


        REAL(KIND=silp), DIMENSION(size_z) :: z_axis_values =                                                                     &
                 (/ -5.0, -15.98, -29.17, -45.2, -64.96, -89.75, -121.52, -163.28, -219.86, -299.29, -415.07, -588.88, -850.19,   &
                    -1225.11, -1717.9, -2307.36, -2963.25, -3661.11, -4385.22, -5126.18 /)


!--- UGLY HACK FOR CLIO GENERATION TO BE REPLACED !!!

      CONTAINS

! ---
      SUBROUTINE WRTE_RESTART_OCYCC

         use marine_bio_mod, only: OPO4, ONO3,  OSI, OALK,  ODIC, ODOC, OPOC, OC13, OC14, ODOCS, ODOC13, ODOCS13                  &
                                , FOPO4, FONO3, FOSI, FOALK, FODIC, FODOC, FODOCS, FOC13, FODOC13
         use C_res_mod,      only: C13ATM
         use carbone_co2,    only: PA_C, PA0_C
         use mbiota_mod,     only: PHYTO_M, ZOO_M


         integer(kind=ip), parameter :: len_char = 7, nb_vars = 26

         character(len=len_char), dimension(:), allocatable :: variables_to_write
         integer(kind=ip)            :: var

         allocate(variables_to_write(nb_vars))
         variables_to_write = [Character(len=len_char) :: "OPO4", "ONO3",  "OSI", "OALK",  "ODIC", "ODOC", "OPOC", "OC13", "OC14" &
                                              , "PHYTO_M", "ZOO_M", "ODOCS", "ODOC13", "ODOCS13", "PA0_C", "PA_C", "C13ATM"       &
                                              , "FOPO4", "FONO3", "FOSI", "FOALK", "FODIC", "FODOC", "FODOCS", "FOC13", "FODOC13" &
                              ]
!~ #if ( OOISO == 1 )
!~              OO2(:,:,:,1),                                              &
!~ #else
!~              OO2,                                                       &
!~ #endif
!~ #if ( OOISO == 1 )
!~              FOO2(:,:,1),                                               &
!~ #else
!~              FOO2,                                                      &
!~ #endif

!~ #if ( OXNITREUX == 1 )
!~                      ON2O ,FON2O,                                       &
!~ #endif
!~ #if ( KC14 == 1 )
!~                      C14ATM, cav_oc14_b, cav_oc_b,                      &
!~ #endif
!~ #if ( OOISO ==1 )
!~                      FODOCS13, OO2(:,:,:,2:NISOO2), FOO2(:,:,2:NISOO2)
!~ #else
!~                      FODOCS13
!~ #endif

         call WRTE_VAR_IN_FILE(OPO4,trim(variables_to_write(1))//"_rest.nc",trim(variables_to_write(1)))
         call WRTE_VAR_IN_FILE(ONO3,trim(variables_to_write(2))//"_rest.nc",trim(variables_to_write(2)))

      END SUBROUTINE WRTE_RESTART_OCYCC


      SUBROUTINE WRTE_VAR_IN_FILE(array_to_wrte,filename,varname)

        use loveclim_transfer_mod, only: MGT
        use IO_NC_MOD, only: undef_dblp
        use ncio, only: nc_read

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


        real(kind=dblp), dimension(:,:,:),         intent(in) :: array_to_wrte
        CHARACTER(len=*),                          intent(in) :: filename, varname

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        TYPE(IO_NC_FILE), TARGET                       :: restart_mb_file
        TYPE(IO_NC_AXIS), TARGET                       :: CLIO_axis_ptlat, CLIO_axis_ptlon, CLIO_axis_depth

!---    CLIO variables
        TYPE(IO_GRID_VAR)                              :: var_towrte
        REAL(kind=dblp), DIMENSION(:,:,:), ALLOCATABLE :: tbl_towrte
        REAL(kind=dblp), DIMENSION(:,:,:), ALLOCATABLE :: tbl_tochck

        INTEGER(kind=silp)                             :: i,j,k
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       THINGS TO DO ONCE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! initialize netCDF files, including axes attributes and initial coordinates variables
        ! for now, do a relatively ugly direct import. Will do better later on.

        call restart_mb_file%init(FileName=filename)
        call CLIO_axis_ptlat%init("ptlat", y_axis_values, OPT_AxisUnit="degrees_north")
        call CLIO_axis_ptlon%init("ptlon", x_axis_values, OPT_AxisUnit="degrees_east")
        call CLIO_axis_depth%init("depth", z_axis_values, OPT_AxisUnit="meters_depth")


        call CLIO_axis_ptlat%wrte(filename)
        call CLIO_axis_ptlon%wrte(filename)
        call CLIO_axis_depth%wrte(filename)

        ! [NOTA] The netCDF file expects lon, lat, depth, but our internal variable is lat, depth, lon
        !         -> needs re-writing, swapping axes, and reducing the longitude because of CLIO overlap (122/120)
        !         -> en termes d'axes : 1,2,3 -> 3*,1,2
        ! [NOTA] The internal format is the exact reverse of netCDF. If you ask for lon, lat, depth, you obtain
        !            depth, lat, lon in the netCDF. Hence the choice made here

        allocate(tbl_towrte(UBOUND(array_to_wrte,DIM=3)-2,UBOUND(array_to_wrte,DIM=1),UBOUND(array_to_wrte,DIM=2)))
        allocate(tbl_tochck(UBOUND(array_to_wrte,DIM=3)-2,UBOUND(array_to_wrte,DIM=1),UBOUND(array_to_wrte,DIM=2)))

        call var_towrte%init(varname,restart_mb_file, "ptlon ptlat depth")
        call var_towrte%show()

        do k=2,UBOUND(array_to_wrte,DIM=3)-1   ! longitude
          do j=1,UBOUND(array_to_wrte,DIM=2)   ! depth
            do i=1,UBOUND(array_to_wrte,DIM=1) ! latitude
              tbl_towrte(k-1,i,j) = array_to_wrte(i,j,k)
              if (MGT(i,j,k).EQ.0) then
                 tbl_towrte(k-1,i,j) = undef_dblp
              endif
            enddo
          enddo
        enddo

         call var_towrte%wrte(tbl_towrte(:,:,:))
! ---
#if ( 0 )
! --- Sanity check .... (could be removed)
         call nc_read(filename,varname,tbl_tochck(:,:,:))
         WRITE(*,*) "SANITY CHECK ::"//TRIM(varname), MINVAL(tbl_towrte(:,:,:)-tbl_tochck(:,:,:)), MAXVAL(tbl_towrte(:,:,:)-tbl_tochck(:,:,:))
#endif
         deallocate(tbl_towrte)
         deallocate(tbl_tochck)

      END SUBROUTINE WRTE_VAR_IN_FILE


      END MODULE WRTE_RESTART_IO_NC
#endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
