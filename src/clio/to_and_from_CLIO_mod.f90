!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2022 Didier M. Roche (a.k.a. dmr)

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


      module to_and_from_CLIO

      implicit none
      
      
      private
      
      public :: get_indexes_C, get_lonlat_C, get_indexes_C_wmask, transfer_to_CLIO
      
!L30      parameter ( dxo = 1.5, dyo = 1.5 )      
!L30  latini = -79.5

      real, parameter :: demitour = 180.0, pi = 4.0 * atan(1.0), latini=-81.0, lonini=24.0, dxo = 3.0, dyo = 3.0
      real, parameter :: lonAAini = 90., xwpoln = -111. + 360., latini2=-79.5, lonini2=25.5, radian = pi / demitour

!~ CL30      lonini2 =  24.75
!~ CL30      latini2 = -78.75

     
      real            :: lon3, lat3, xx, yy, latm1, latm2, latt, lont, zzz, clat, clon, slat, slon, lat, lon, latAA, lonAA

      integer, parameter :: iocn = 122, jocn = 65, ijocn = iocn*jocn, jsepar = 50, isepar = 82
      integer            :: jsep(iocn), jeq, inumo, jnumo
      integer            :: lat_dat_id, lon_dat_id


      real, dimension(iocn,jocn) :: zlatt, zlont
      
      logical :: is_initd = .false.


      contains 

      subroutine mod_init 
      
! --- reading helper files

       integer :: i, j


      open (newunit=lat_dat_id,file='outputdata/ocean/lat.dat')
      do j=1,jocn
         read(lat_dat_id,'(122( F10.5))' ) (zlatt(i,j),i=1,iocn)
!L30         write(lat_dat_id,'(242( F10.5))' ) (zlatt(i,j),i=1,imax)
      enddo
      close (lat_dat_id)

      open (newunit=lon_dat_id,file='outputdata/ocean/long.dat')
      do j=1,jocn
         read(lon_dat_id,'(122( F10.5))' ) (zlont(i,j),i=1,iocn)
!L30         write(lon_dat_id,'(242( F10.5))' ) (zlont(i,j),i=1,imax)

      enddo
      write(lon_dat_id,*)
      close (lon_dat_id)

      
!L30      parameter ( iocn = 242, jocn = 128, ijocn = iocn*jocn )
!L30      parameter ( isepar = 180, jsepar=98 )
      
      
!-- calcul de la limite selon j pour la separation des 2 grilles: jsep(i)
      jeq = nint(0.5 - latini / dyo) + 1
      if (jsepar.eq.jocn) then
        do i=1,iocn
          jsep(i) = jsepar
        enddo
      else
!-- limite (incluse ds AA) entre les 2 grilles : x = 296.E - y.N, jeq<j<jsepar
         do i=1,iocn
            xx = lonini + dxo * float(i-1)
            yy = 296. - xx
            j = nint( (yy - latini) / dyo  ) + 1
            jsep(i) = min(max(j,jeq),jsepar+1)
         enddo
      endif
         
      is_initd = .true.
      
      end subroutine mod_init 
      

      function get_lonlat_C(io,jo) result(lonlat)
      
          integer, intent(in) ::  io, jo
          real, dimension(2)  :: lonlat

!~         determiner sur la grille du modele oceanique pour la maille ocn ijo, 
    
    
          if (.not. is_initd) then
             call mod_init
          endif
    
    
!     lonini est la longitude minimale pour les points vecteurs (position3)
!     latini est la latitude minimale pour les points vecteurs (position3)
         lon3 = lonini + (io-1)*dxo
         lat3 = latini + (jo-1)*dyo

         if ( jo.ge.jsep(io) ) then
            latm1 = (-111.0+360+90-lon3)*pi/demitour
            latm2 = (-111.0+360+90-lon3-dyo)*pi/demitour
         else
            latm1 = (latini + (jo-1)*dyo)*pi/demitour
            latm2 = (latini + (jo-1)*dyo + dyo)*pi/demitour
         endif

      
         if ( jo.ge.jsep(io) ) then
  	        latt = lat3 ! + ran(iseed )*dyo
            zzz=0.0 ! for now, could something else ... was ran(iseed )

            lont = lon3 - (asin((1.-zzz)*sin(latm1)+zzz*sin(latm2))-latm1)/radian
!~          write(*,*) lat3,latt,lon3,lont
         else
	        lont = lon3 ! + ran(iseed )*dxo
            zzz=0.0  ! ran(iseed )
            latt = lat3 + (asin((1.-zzz)*sin(latm1)+zzz*sin(latm2))-latm1)/radian
!~          write(*,*) lat3,latt,lon3,lont
         endif


! from mod_geo ...
!     io -> inumo
!     jo -> jnumo

        inumo = io
        jnumo = jo

!     lont -> lon
!     latt -> lat

        lon = lont
        lat = latt

!--   clon est le cosinus de la longitide, slon le sinus
!-- lonini -> lonini2 /      latini -> latini2

!-- Passage de la grille fictive (batit a partir de la carte du modele 
!-- en ajoutant 3o en longitude et latitude quand on croit en indice)
!-- a la grille tournee "AA".
!-- cela signifie pour les points j> jsep (limite entre WW et AA). 
!-- Pour cela, il faut 
!-- modifier les lat et lon (en lat' et lon') du fait que la grille soit tournee
!-- et utiliser les formules de passage de la grille AA en coordonnees 
!-- geographique, en pensant a repasser a la reference de Greenwich
!-- pour les longitudes.
 
         if ( jnumo.ge.jsep(inumo) ) then
!--      2/ modifier les lon et lat en lonAA et latAA dans la grille AA
            lonAA= lonAAini + lat 
            latAA= xwpoln + 90. - lon 
!
!--      passage de la grille AA en WW non tournee
            clon= cos(lonAA*radian)
            slon= sin(lonAA*radian)
            clat= cos(latAA*radian)
            slat= sin(latAA*radian)
!
            lon= ( atan2(clat*slon,slat) )/radian
            lat= ( -asin(clat*clon) )/radian
!
!--      re-normalisation de la longitude (0o pour le meridien de Greenwich)
            lon= mod(xwpoln + lon + 360., 360.) 
         endif
!
!     ecriture de la longitude entre -180 et 180
!     if (lon.gt.180.) then
!        lon= lon-360.
!     endif
!     ecriture de la longitude entre 0 et 360
         if (lon.gt.360.) then
            lon= lon-360.
         endif

         lonlat(1) = lon
         lonlat(2) = lat

      end function get_lonlat_C

        function get_indexes_C(lon,lat) result(index_close)
          real, intent(in) :: lon, lat
          integer          :: index_close(2)

          if (.not. is_initd) then
             call mod_init
          endif

         index_close = find_closest_2D(lon*radian, lat*radian,zlont, zlatt)
     
       end function get_indexes_C

        function get_indexes_C_wmask(lon,lat,mask_array,value_mask) result(index_close)
          real, intent(in) :: lon, lat
          integer, dimension(:,:), intent(in) :: mask_array
          integer, intent(in) :: value_mask
          
          integer          :: index_close(2)

          if (.not. is_initd) then
             call mod_init
          endif

         index_close = find_closest_2D_masked(lon*radian, lat*radian,zlont, zlatt,mask_array,value_mask)
     
       end function get_indexes_C_wmask

       elemental pure function haversine_distance(lon1, lat1, lon2, lat2) result (h_dist)
       
         real, intent(in) :: lat1, lat2, lon1, lon2
         real             :: dlat, dlon, q, h_dist
       
         dlat = lat2 - lat1
         dlon = lon2 - lon1
         q = sin(dlat/2)**2 + (cos(lat1) * cos(lat2) * (sin(dlon/2)**2))
         h_dist = 2 * atan2(sqrt(q), sqrt(1-q))

       end function haversine_distance

       function find_closest_2D(tagt_lon,tagt_lat,lons_array,lats_array) result(indx_loc)

          real, intent(in)                  :: tagt_lon, tagt_lat
          real, dimension(:,:), intent(in)  :: lons_array, lats_array

          real, dimension(:,:), allocatable :: work_array
          
          integer :: indx_loc(2), size_array_1, size_array_2

          size_array_1 = UBOUND(lons_array, DIM=1)
          size_array_2 = UBOUND(lons_array, DIM=2)
          allocate(work_array(size_array_1, size_array_2))

          work_array = haversine_distance(tagt_lon, tagt_lat, lons_array, lats_array)

          indx_loc = minloc(work_array)

       end function find_closest_2D

       function find_closest_2D_masked(tagt_lon,tagt_lat,lons_array,lats_array, mask_array, value_mask) result(indx_loc)

          real, intent(in)                    :: tagt_lon, tagt_lat
          real, dimension(:,:), intent(in)    :: lons_array, lats_array
          integer, dimension(:,:), intent(in) :: mask_array
          integer, intent(in)                 :: value_mask

          real, dimension(:,:), allocatable   :: work_array, masked_wkarray
                    
          integer :: indx_loc(2), size_array_1, size_array_2

          size_array_1 = UBOUND(lons_array, DIM=1)
          size_array_2 = UBOUND(lons_array, DIM=2)
          
          allocate(work_array(size_array_1, size_array_2))
          allocate(masked_wkarray(size_array_1, size_array_2))

          work_array = haversine_distance(tagt_lon, tagt_lat, lons_array, lats_array)
          
          !dmr&afq --- Penalty when not belonging to the right area
          masked_wkarray = merge(work_array,work_array+1000000.,mask_array==value_mask)

          indx_loc = minloc(masked_wkarray)

       end function find_closest_2D_masked

       function transfer_to_CLIO(input_field, input_lat, input_lon, clio_mask) result(clio_field)
       
          integer, dimension(:,:),       intent(in)           :: input_field
          real   , dimension(:,:),       intent(in)           :: input_lat, input_lon
          integer, dimension(iocn,jocn), intent(in), optional :: clio_mask
       
          integer, dimension(iocn,jocn)                       :: clio_field
          integer                                             :: i,j
          integer, dimension(2)                               :: index_clio
          
          if ( PRESENT(clio_mask)) then
            clio_field = merge(0, 99 , clio_mask .eq. 0)

          
          do j=1, UBOUND(input_field,dim=2)
          do i=1, UBOUND(input_field,dim=1)
             index_clio = get_indexes_C(input_lon(i,j), input_lat(i,j))            
            if ( clio_mask(index_clio(1), index_clio(2)) .gt. 0 ) then 
                clio_field(index_clio(1), index_clio(2)) = input_field(i,j)
            endif
          enddo
          enddo
          
          
          else
          
          do j=1, UBOUND(input_field,dim=2)
          do i=1, UBOUND(input_field,dim=1)
            index_clio = get_indexes_C(input_lon(i,j), input_lat(i,j))            
            clio_field(index_clio(1), index_clio(2)) = input_field(i,j)
          enddo
          enddo
                       
          endif ! on the optional argument       
       end function transfer_to_CLIO

      
      end module to_and_from_CLIO

!  for flattening / unflattening
!         j = int((ijo-1)/iocn)+1
!         i = ijo-(j-1)*iocn

!         ijo = (jo-1)*iocn + io
