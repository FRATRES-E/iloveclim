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

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: global_constants_mod
!
!>     @author  Didier M. Roche (dmr)
!
!
!>     @brief This module is putting together in Fortran90 form the gather&scater routines
!
!>     @date Creation date: June, 24th, 2020
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module gather_and_scater_mod

      implicit none

      private
            
      public:: gather, scater
      
      contains
      

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      subroutine gather(n,a,b,index_t)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use global_constants_mod, only: dblp=>dp, ip      

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       
      real(kind=dblp), dimension(n), intent(out)::  a
      integer(kind=ip), dimension(n), intent(in)::  index_t
      real(kind=dblp), dimension(*), intent(in) ::  b


      integer(kind=ip) :: n, ji
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|  

      do ji=1,n
         a(ji)=b(index_t(ji))
      enddo

      return
      end subroutine gather


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|      
      subroutine scater(n,a,index_t,b)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use global_constants_mod, only: dblp=>dp, ip           

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=dblp), dimension(n), intent(in) :: b
      integer(kind=ip), dimension(n), intent(in):: index_t
      real(kind=dblp), dimension(*), intent(out):: a
      
      
      integer(kind=ip) :: n, ji
 
      do ji=1,n
         a(index_t(ji))=b(ji)
      enddo

      return
      end subroutine scater

      end module gather_and_scater_mod
