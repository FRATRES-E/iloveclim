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


      module delta_latlon

      use global_constants_mod, only: dblp=>dp, ip

      implicit none

      private

      public :: get_Nlatlon

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      elemental pure subroutine get_Nlatlon(uu,vv,deltaT,oLat,oLon,nLat,nLon)

        use global_constants_mod, only: deg_to_rad, rad_to_deg

        real(kind=dblp), intent(in) :: uu, vv, oLat, oLon
        real(kind=dblp), intent(out):: nLat, nLon
        integer(kind=ip), intent(in) :: deltaT
        real(kind=dblp) :: bearing_rad, d_togo, oLat_rad, oLon_rad, nLat_rad, nLon_rad


        bearing_rad = get_bearing(uu,vv)
        d_togo = distance_to_go(uu,vv,deltaT)

        oLat_rad = oLat * deg_to_rad
        oLon_rad = oLon * deg_to_rad

        call get_newlatlon(oLat_rad,oLon_rad,d_togo,bearing_rad,nLat_rad,nLon_rad)

        nLat = nLat_rad * rad_to_deg
        nLon = nLon_rad * rad_to_deg

      end subroutine get_Nlatlon


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      elemental pure function get_bearing(uu,vv) result(bearing_rad)

        real(kind=dblp), intent(in) :: uu, vv
        real(kind=dblp)             :: bearing_rad

        bearing_rad = atan2(vv,uu)

      end function get_bearing


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      elemental pure function distance_to_go(uu,vv,deltaT) result(ToGo)

        real(kind=dblp), intent(in) :: uu,vv
        integer(kind=ip), intent(in) :: deltaT

        real(kind=dblp)             :: normUV, ToGo

        normUV = sqrt(uu**2+vv**2)
        ToGo   = normUV * deltaT

      end function distance_to_go


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      elemental pure subroutine get_newlatlon(oLat,oLon,distance,bearing,nLat,nLon)

      use global_constants_mod, only: radius_earth

        ! nota all angles are in radians here
        ! nota need R = radius of Earth

        real(kind=dblp), intent(in) :: oLat, oLon, distance, bearing
        real(kind=dblp), intent(out):: nLat, nLon

        nLat = asin( sin(oLat)*cos(distance/radius_earth) + cos(oLat)*sin(distance/radius_earth)*cos(bearing))
        nLon = oLon + atan2(sin(bearing)*sin(distance/radius_earth)*cos(oLat),cos(distance/radius_earth)-sin(oLat)*sin(nLat))

      end subroutine get_newlatlon


      end module delta_latlon
