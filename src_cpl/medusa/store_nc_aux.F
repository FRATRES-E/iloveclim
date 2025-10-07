************************************************************************
      SUBROUTINE STORE_NC_AUX(filename, title_string)
************************************************************************

! Add surface area information
      USE mod_defines_medusa
      USE mod_execontrol_medusa, ONLY: ABORT_MEDUSA
      USE mod_gridparam
      USE mod_indexparam
      USE mod_materialcharas

      USE mod_seafloor_central, ONLY: N_COLUMNS_USED, IJ_COLUMNS_USED,
     &                              COLUMN_IJ2N, COLUMN_N2IJ,
     &                              COLUMN_N2XY, COLUMN_AREA4N

      USE mod_netcdfparam
      use mod_netcdfinc, only: HANDLE_ERRORS

      IMPLICIT NONE


#include "netcdf.inc" 

      CHARACTER(LEN=*), INTENT(IN)           :: filename
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: title_string

      INTEGER :: n_grid_seafloor
      INTEGER :: iflag

      DOUBLE PRECISION, DIMENSION(idnw:idnb) :: xzdn


      INTEGER, DIMENSION(3) :: dim, start, count

      INTEGER :: istatus, k
      INTEGER :: ncid, iccf
      INTEGER :: i, j, n

! Dimensions and dimension variables
!  'lev' (levels)
      INTEGER :: dim_lev      ! dimension ID
      INTEGER ::  id_lev      ! ID of dimension variable
! 'col' (columns)
      INTEGER :: dim_col      ! dimension ID
      INTEGER ::  id_col      ! ID of dimension variable
! 'ilon' ('longitudinal' index of the iLOVECLIM grid)
      INTEGER :: dim_ilon
      INTEGER ::  id_ilon
! 'jlat' ('latitudinal' index of the iLOVECLIM grid)
      INTEGER :: dim_jlat
      INTEGER ::  id_jlat

#ifdef STORENCAUX_POLYGONS
! 'ngeped' (grid element polygon edge index index)
      INTEGER, PARAMETER :: ngeped = 5
      INTEGER :: dim_ngeped
      INTEGER ::  id_ngeped
#endif
      DOUBLE PRECISION, PARAMETER :: badlola = 1D36


! Data variables
      INTEGER ::  id_xzdn      ! Depth coordinates
      INTEGER ::  id_col4ij    ! column index n for grid element (i,j)
      INTEGER ::  id_lon4ij    ! longitude of grid element (i,j)
      INTEGER ::  id_lat4ij    ! latitude of grid element (i,j)
      INTEGER ::  id_i4col     ! index i(n)
      INTEGER ::  id_j4col     ! index j(n)
#ifdef STORENCAUX_POLYGONS
      INTEGER ::  id_lon4ijgep ! longitudes of grid elt polygon edges
      INTEGER ::  id_lat4ijgep ! latitudes of grid elt polygon edges
      INTEGER ::  id_lon4ngep  ! longitudes of grid elt polygon edges
      INTEGER ::  id_lat4ngep  ! latitudes of grid elt polygon edges
#endif
      INTEGER ::  id_area      ! surface areas of columns (n)

      DOUBLE PRECISION :: xlon_ij
      DOUBLE PRECISION :: xlat_ij
#ifdef STORENCAUX_POLYGONS
      DOUBLE PRECISION, DIMENSION(ngeped)      :: xlongep, xlatgep
#endif

      INTEGER :: nlons = -1, nlats = -1
      DOUBLE PRECISION                         :: sfc_area


! STORE_NC_AUX is outside the module, and so has to use the interface
      CALL N_COLUMNS_USED(n_grid_seafloor)
      CALL IJ_COLUMNS_USED(nix_columns = nlons, njy_columns = nlats)


! If n_grid_seafloor == -1, the basic setup has not yet been
! done, and it does not make any sense to proceed
      IF (n_grid_seafloor == -1) THEN
        WRITE(jp_stderr,*) '[STORE_NC_AUX] error: ' //
     &    'Sea-floor setup not yet done -- aborting!'
        CALL ABORT_MEDUSA()
      ENDIF


!-----------------------
! Create the data file
!-----------------------
      istatus = NF_CREATE(TRIM(filename), NF_SHARE, ncid)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


!-----------------------
! Define dimensions and dimension variables
!-----------------------
!
! Columns on the string ('col')              dim_col, n_grid_seafloor
! Levels ('lev' )                            dim_lev, ndn
! 'Longitudinal' index of iLOVECLIM ('ilon') dim_ilon, nlons
! 'Latitudinal' index of iLOVECLIM ('jlat')  dim_jlat, nlats
!
! (Id's of dimension variables for dim_xxx are id_xxx)

! Columns on the string
      istatus = NF_DEF_DIM(ncid,
     &                    'col', n_grid_seafloor,
     &                    dim_col)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_DEF_VAR(ncid, 'col', NF_INT,
     &                    1, dim_col,
     &                    id_col)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

! Levels
      istatus = NF_DEF_DIM(ncid,
     &                    'lev', ndn,
     &                    dim_lev)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_DEF_VAR(ncid, 'lev', NF_INT,
     &                    1, dim_lev,
     &                    id_lev)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lev, 'positive',
     &                       4, 'down')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


! 'Longitudinal' index of the iLOVECLIM grid
      istatus = NF_DEF_DIM(ncid, 'ilon',
     &                    nlons, dim_ilon)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

      istatus = NF_DEF_VAR(ncid, 'ilon', NF_INT,
     &                    1, dim_ilon, id_ilon)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

      istatus = NF_PUT_ATT_TEXT(ncid, id_ilon, 'long_name',
     &                       20, 'LVCL_Longitude_Index')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


! 'Latitudinal' index of the iLOVECLIM grid
      istatus = NF_DEF_DIM(ncid, 'jlat',
     &                    nlats, dim_jlat)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

      istatus = NF_DEF_VAR(ncid, 'jlat', NF_INT,
     &                    1, dim_jlat, id_jlat)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

      istatus = NF_PUT_ATT_TEXT(ncid, id_jlat, 'long_name',
     &                         19, 'LVCL_Latitude_Index')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


#ifdef STORENCAUX_POLYGONS
! Polygon edge index of grid elements
      istatus = NF_DEF_DIM(ncid, 'ngeped',
     &                    ngeped, dim_ngeped)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

      istatus = NF_DEF_VAR(ncid, 'ngeped', NF_INT,
     &                    1, dim_ngeped, id_ngeped)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

      istatus = NF_PUT_ATT_TEXT(ncid, id_ngeped, 'long_name',
     &                         24, 'Grid_El_Polyg_Edge_index')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
#endif

         
!-----------------------
! Define data variables
!-----------------------

! xzdn(lev)
      istatus = NF_DEF_VAR(ncid, vsn_xzdn, NF_DOUBLE,
     &                    1, dim_lev,
     &                    id_xzdn)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_xzdn,
     &                         'long_name', vll_xzdn, vln_xzdn)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_xzdn,
     &                         'units', ul_m, un_m)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


! col4ij(ilon,jlat)
      dim(1) = dim_ilon
      dim(2) = dim_jlat
      
      istatus = NF_DEF_VAR(ncid, 'col4ij', NF_INT,
     &                    2, dim(1:2),
     &                    id_col4ij)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_col4ij, 'long_name',
     &                         17, 'Column_index_of_ij')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_INT(ncid, id_col4ij,
     &                        '_FillValue', NF_INT,
     &                        1, -1)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_INT(ncid, id_col4ij,
     &                        'valid_min', NF_INT,
     &                        1, 1)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_INT(ncid, id_col4ij,
     &                        'valid_max', NF_INT,
     &                        1, n_grid_seafloor)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


! lon4ij(ilon,jlat)
      dim(1) = dim_ilon
      dim(2) = dim_jlat
      
      istatus = NF_DEF_VAR(ncid, 'lon4ij', NF_DOUBLE,
     &                    2, dim(1:2),
     &                    id_lon4ij)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lon4ij, 'long_name',
     &                         9, 'Longitude')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lon4ij, 'units',
     &                         ul_degE, un_degE)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


! lat4ij(ilon,jlat)
      dim(1) = dim_ilon
      dim(2) = dim_jlat

      istatus = NF_DEF_VAR(ncid, 'lat4ij', NF_DOUBLE,
     &                    2, dim(1:2),
     &                    id_lat4ij)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lat4ij, 'long_name',
     &                         8, 'Latitude')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lat4ij, 'units',
     &                         ul_degN, un_degN)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


! i4col(col)
      istatus = NF_DEF_VAR(ncid, 'i4col', NF_INT,
     &                    1, dim_col,
     &                    id_i4col)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_i4col, 'long_name',
     &                   33, 'Longitude_Index_from_Column_Index')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_INT(ncid, id_i4col,
     &                        '_FillValue', NF_INT,
     &                        1, -1)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_INT(ncid, id_i4col,
     &                        'valid_min', NF_INT,
     &                        1, 1)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_INT(ncid, id_i4col,
     &                        'valid_max', NF_INT,
     &                        1, nlons)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


! j4col(col)
      istatus = NF_DEF_VAR(ncid, 'j4col', NF_INT,
     &                    1, dim_col,
     &                    id_j4col)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_j4col, 'long_name',
     &                   32, 'Latitude_Index_from_Column_Index')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_INT(ncid, id_j4col,
     &                        '_FillValue', NF_INT,
     &                        1, -1)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_INT(ncid, id_j4col,
     &                        'valid_min', NF_INT,
     &                        1, 1)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_INT(ncid, id_j4col,
     &                        'valid_max', NF_INT,
     &                        1, nlats)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

#ifdef STORENCAUX_POLYGONS
! lon4ijgep(ilon,jlat,ngeped)
      dim(1) = dim_ilon
      dim(2) = dim_jlat
      dim(3) = dim_ngeped

      istatus = NF_DEF_VAR(ncid, 'lon4ijgep', NF_DOUBLE,
     &                    3, dim(1:3),
     &                    id_lon4ijgep)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lon4ijgep, 'long_name',
     &                         9, 'Longitude')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lon4ijgep, 'units',
     &                         ul_degE, un_degE)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_DOUBLE(ncid, id_lon4ijgep,
     &                        '_FillValue', NF_DOUBLE,
     &                        1, badlola)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


! lat4ijgep(ilon,jlat,ngeped)
      dim(1) = dim_ilon
      dim(2) = dim_jlat
      dim(3) = dim_ngeped

      istatus = NF_DEF_VAR(ncid, 'lat4ijgep', NF_DOUBLE,
     &                    3, dim(1:3),
     &                    id_lat4ijgep)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lat4ijgep, 'long_name',
     &                         8, 'Latitude')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lat4ijgep, 'units',
     &                         ul_degN, un_degN)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_DOUBLE(ncid, id_lat4ijgep,
     &                        '_FillValue', NF_DOUBLE,
     &                        1, badlola)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

! lon4ngep(col,ngeped)
      dim(1) = dim_col
      dim(2) = dim_ngeped

      istatus = NF_DEF_VAR(ncid, 'lon4ngep', NF_DOUBLE,
     &                    2, dim(1:2),
     &                    id_lon4ngep)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lon4ngep, 'long_name',
     &                         9, 'Longitude')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lon4ngep, 'units',
     &                         ul_degE, un_degE)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_DOUBLE(ncid, id_lon4ngep,
     &                        '_FillValue', NF_DOUBLE,
     &                        1, badlola)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


! lat4ngep(col,ngeped)
      dim(1) = dim_col
      dim(2) = dim_ngeped
      
      istatus = NF_DEF_VAR(ncid, 'lat4ngep', NF_DOUBLE,
     &                    2, dim(1:2),
     &                    id_lat4ngep)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lat4ngep, 'long_name',
     &                         8, 'Latitude')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_lat4ngep, 'units',
     &                         ul_degN, un_degN)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_DOUBLE(ncid, id_lat4ngep,
     &                        '_FillValue', NF_DOUBLE,
     &                        1, badlola)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
#endif

! area(col)
      istatus = NF_DEF_VAR(ncid, 'area', NF_DOUBLE,
     &                    1, dim_col,
     &                    id_area)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_area, 'long_name',
     &                   12, 'Surface_Area')
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      istatus = NF_PUT_ATT_TEXT(ncid, id_area, 'units',
     &                         ul_m2, un_m2)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


!----------------------
! Put global attributes
!----------------------

      IF(PRESENT(title_string)) THEN
        istatus = NF_PUT_ATT_TEXT(ncid, NF_GLOBAL, 'title',
     &                          LEN_TRIM(title_string),title_string)
        IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      ENDIF


!--------------------
! End define mode
!--------------------
      
      istatus = NF_ENDDEF(ncid)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


      DO k = 1, n_grid_seafloor
        i = k
        istatus = NF_PUT_VAR1_INT(ncid, id_col, i, k)
        IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      ENDDO

      DO k = idnt, idnb
        i = k-idnt+1
        istatus = NF_PUT_VAR1_INT(ncid, id_lev, i, k)
        IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      ENDDO

      DO k = 1, nlons
        i = k
        istatus = NF_PUT_VAR1_INT(ncid, id_ilon, i, k)
        IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      ENDDO

      DO k = 1, nlats
        i = k
        istatus = NF_PUT_VAR1_INT(ncid, id_jlat, i, k)
        IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      ENDDO

#ifdef STORENCAUX_POLYGONS
      DO k = 1, ngeped
        i = k
        istatus = NF_PUT_VAR1_INT(ncid, id_ngeped, i, k)
        IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      ENDDO
#endif

! xzdn(lev)
      CALL GRID_DEF(xzdn)
      istatus = NF_PUT_VAR_DOUBLE(ncid, id_xzdn, xzdn)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)


! col4ij
      DO i = 1, nlons
        start(1) = i
        DO j = 1, nlats
          start(2) = j
          CALL COLUMN_IJ2N(i, j, iflag, k)
          IF (iflag == 0) THEN
            istatus = NF_PUT_VAR1_INT(ncid, id_col4ij,
     &                               start(1:2), k)
            IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
          ENDIF
        ENDDO
      ENDDO


! lon4ij, lat4ij

      DO i = 1, nlons
        start(1) = i
        DO j = 1, nlats
          start(2) = j
          CALL COLUMN_IJ2N(i, j, iflag, k)
          CALL COLUMN_N2XY(k, iflag, xlon_ij, xlat_ij)
          istatus = NF_PUT_VAR1_DOUBLE(ncid, id_lon4ij,
     &                                start(1:2), xlon_ij)
          IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
          istatus = NF_PUT_VAR1_DOUBLE(ncid, id_lat4ij,
     &                                start(1:2), xlat_ij)
          IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
        ENDDO
      ENDDO

! icol, jcol
      DO k = 1, n_grid_seafloor
        CALL COLUMN_N2IJ(k, iflag, i, j)
        IF (iflag == 0) THEN
          istatus = NF_PUT_VAR1_INT(ncid, id_i4col, k, i)
          IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
          istatus = NF_PUT_VAR1_INT(ncid, id_j4col, k, j)
          IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
        ELSE
          WRITE(jp_stderr,*) 'Inconsistency regarding COLUMN_N2IJ call:'
          WRITE(jp_stderr,*) ' called with k = ',k
          WRITE(jp_stderr,*) ' returned iflag = ', iflag
          CALL ABORT_MEDUSA()
        ENDIF
      ENDDO


#ifdef STORENCAUX_POLYGONS
!Longitudes and latitudes of GridElementPolygons
      DO j = 1, nlats
        IF(j == 1) THEN
          xlatgep(1) = xlat_ij(1)
          xlatgep(2) = 90D0
          xlatgep(3) = 90D0
          xlatgep(4) = xlat_ij(1)
          xlatgep(5) = xlat_ij(2)
        ELSEIF(j == nlats) THEN
          xlatgep(1) = xlat_ij(nlats)
          xlatgep(2) = xlat_ij(nlats-1)
          xlatgep(3) = xlat_ij(nlats)
          xlatgep(4) = -90D0
          xlatgep(5) = -90D0
        ELSE
          xlatgep(1) = xlat_ij(j)
          xlatgep(2) = xlat_ij(j-1)
          xlatgep(3) = xlat_ij(j)
          xlatgep(4) = xlat_ij(j+1)
          xlatgep(5) = badlola
        ENDIF

        DO i = 1, nlons
          IF((i == 1).AND.(j==nlats)) THEN
            xlongep(1) = xlon_ij(1, nlats-1)
            xlongep(2) = xlon_ij(1, nlats  )
            xlongep(3) = xlon_ij(1, nlats-1) - 360D0/(DBLE(nlons))
            xlongep(4) = xlon_ij(1, nlats-1) - 360D0/(DBLE(nlons))
            xlongep(5) = xlon_ij(1, nlats-1)
          ELSEIF((i == nlons).AND.(j == 1)) THEN
            xlongep(1) = xlon_ij(nlons, 2) + 360D0/(DBLE(nlons))
            xlongep(2) = xlon_ij(nlons, 2) + 360D0/(DBLE(nlons))
            xlongep(3) = xlon_ij(nlons, 2)
            xlongep(4) = xlon_ij(nlons, 2)
            xlongep(5) = xlon_ij(nlons, 1)
          ELSEIF(j == 1) THEN
            xlongep(1) = xlon_ij(i+1, 2)
            xlongep(2) = xlon_ij(i+1, 2)
            xlongep(3) = xlon_ij(i  , 2)
            xlongep(4) = xlon_ij(i  , 2)
            xlongep(5) = xlon_ij(i  , 1)
          ELSEIF(j == nlats) THEN
            xlongep(1) = xlon_ij(i  , nlats-1)
            xlongep(2) = xlon_ij(i  , nlats  )
            xlongep(3) = xlon_ij(i-1, nlats-1)
            xlongep(4) = xlon_ij(i-1, nlats-1)
            xlongep(5) = xlon_ij(i  , nlats-1)
          ELSE
            xlongep(1) = xlon_ij(i, j-1)
            xlongep(2) = xlon_ij(i, j  )
            xlongep(3) = xlon_ij(i, j+1)
            xlongep(4) = xlon_ij(i, j  )
            xlongep(5) = badlola
          ENDIF

! First save the corner coordinates for the ij-array-elements
          start(1) = i
          start(2) = j
          start(3) = 1

          count(1) = 1
          count(2) = 1
          count(3) = ngeped
          
          istatus = NF_PUT_VARA_DOUBLE(ncid, id_lon4ijgep,
     &                                start(1:3), count(1:3), xlongep)
          IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

          istatus = NF_PUT_VARA_DOUBLE(ncid, id_lat4ijgep,
     &                                start(1:3), count(1:3), xlatgep)
          IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

! Then check whether there is a column associated with the
! current ij-array-element
          CALL COLUMN_IJ2N(i, j, iflag, k)
! if so, also save the corner coordinates to the associated
! k-array-element
          IF (iflag == 0) THEN

            start(1) = k
            start(2) = 1

            count(1) = 1
            count(2) = ngeped

            istatus = NF_PUT_VARA_DOUBLE(ncid, id_lon4ngep,
     &                                  start(1:2), count(1:2), xlongep)
            IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

            istatus = NF_PUT_VARA_DOUBLE(ncid, id_lat4ngep,
     &                                  start(1:2), count(1:2), xlatgep)
            IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

          ENDIF

        ENDDO
      ENDDO
#endif

! Surface areas of grid_elements
      DO k = 1, n_grid_seafloor
        CALL COLUMN_AREA4N(k, iflag, sfc_area)
        istatus = NF_PUT_VAR1_DOUBLE(ncid, id_area, k, sfc_area)
        IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)
      ENDDO


      istatus = NF_CLOSE(ncid)
      IF (istatus /= NF_NOERR) CALL HANDLE_ERRORS(istatus)

      RETURN
      END SUBROUTINE STORE_NC_AUX
