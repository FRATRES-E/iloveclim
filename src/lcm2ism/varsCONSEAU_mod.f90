!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of the iLOVECLIM coupled climate model under the LUDUS framework.
!!      global_constants_mod is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
!!      License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
!!      version.
!!
!!      global_constants_mod is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
!!      warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!!
!!      You should have received a copy of the GNU General Public License along with Foobar.
!!      If not, see <http://www.gnu.org/licenses/>.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module varsCONSEAU_mod

       use taillesGrilles, only: iecb, jecb, CNX, CNY, sgnxm, sgnym, sgd

       implicit none
       public

! afq -- the GRISLI variables related to water conservation, yearly accumulated in GRISLI

       ! Generic multigrid volume variation:
       double precision                          :: trendgris
       !surface mass balance(either positive or negative):
       double precision, dimension(sgnxm,sgnym,sgd)        :: smbgris
       !basal melting rates (either grounded or shelves):
       double precision, dimension(sgnxm,sgnym,sgd)        :: bmeltgris
       !calving rate:
       double precision, dimension(sgnxm,sgnym,sgd)        :: calgris


       ! the volume variation of the ice sheets in the North (South can be added when ready):
       double precision                , pointer :: trendgrisnord
       !surface mass balance(either positive or negative):
       double precision, dimension(:,:), pointer :: smbgrisnord
       !basal melting rates (either grounded or shelves):
       double precision, dimension(:,:), pointer :: bmeltgrisnord
       !calving rate:
       double precision, dimension(:,:), pointer :: calgrisnord

! afq&dmr --
       double precision, dimension(iEcb, jEcb) :: flux_gris2ecb
       double precision, dimension(CNX,CNY)    :: calgrisCLIO, calgrisCLIOms


!---dmr&afq
      real(kind=8), dimension(CNX,CNY)         :: fwatconseau


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     end module varsCONSEAU_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
