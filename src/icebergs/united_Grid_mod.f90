!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Copyright 2021 Didier M. Roche (a.k.a. dmr)

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


       MODULE united_Grid_mod 


       USE global_constants_mod, only: pi_dp, deg_to_rad, dblp=>dp, ip, rad_to_deg

       IMPLICIT NONE

       PRIVATE

       PUBLIC :: test, getUnifiedGrid_indx, grid_elmnt, init_module

       LOGICAL :: is_initialized = .FALSE.

       REAL(KIND=dblp),  PARAMETER  :: lambda_poleN =  69.0_dblp * deg_to_rad
       REAL(KIND=dblp),  PARAMETER  :: lambda_minAA =  296._dblp ! longitude of index point 1 in AA
       REAL(KIND=dblp),  PARAMETER  :: lambda_minWW =  25.5_dblp ! longitude of index point 1 in WW
       REAL(KIND=dblp),  PARAMETER  :: lambda_refAA =  91.5_dblp ! longitude of index point 28 in AA
       REAL(KIND=dblp),  PARAMETER  :: theta_minWW  = -79.5_dblp ! latitude  of index point 1 in WW
       REAL(KIND=dblp),  PARAMETER  :: dlambda      =   3.0_dblp ! model resolution in longitude
       REAL(KIND=dblp),  PARAMETER  :: dtheta       =   3.0_dblp ! model resolution in longitude
       
       REAL(KIND=dblp),  PARAMETER  :: dblePI       =   2.0_dblp * pi_dp
       
       INTEGER(KIND=ip), PARAMETER :: imax  = 122
       INTEGER(KIND=ip), PARAMETER :: jmn_AA = 28   ! index in WW of the first AA cell
       INTEGER(KIND=ip), PARAMETER :: jmx_WW = 50   ! northernmost possible index for WW
       INTEGER(KIND=ip), PARAMETER :: jeq_AA = 105  ! index of the "equator" of AA
       INTEGER(KIND=ip), PARAMETER :: ist_AA = 28   ! index of the "equator" of AA

       REAL(KIND=dblp), DIMENSION(imax) :: jsep

       TYPE :: grid_elmnt
          INTEGER(KIND=ip) :: X_idx
          INTEGER(KIND=ip) :: Y_idx
          LOGICAL          :: is_inWW
       END TYPE


!   with lambda, lambda_p the longitudes
!   with theta, theta_p the latitudes
!   relationships to be included are:
!         theta  = arcsin(cos(lambda_p)*cos(theta_p))
!         lambda =  atan2(sin(lambda_p)*cos(theta_p), -sin(theta_p))

!         theta_p  = -arcsin(cos(lambda-lambda_PN)*cos(theta)
!         lambda_p = atan2(sin(lambda-lambda_PN)*cos(theta),sin(theta)) + pi/2

       CONTAINS

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   PUBLIC FUNCTIONS AND SUBROUTINES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE init_module()

         if (.NOT. is_initialized) then
            call init_jsepar()
            is_initialized = .TRUE.
         endif

       END SUBROUTINE init_module


       SUBROUTINE test(lambda_in,theta_in,lambda_out,theta_out)
       
         REAL(KIND=dblp), intent(in)  :: lambda_in, theta_in
         REAL(KIND=dblp), intent(out) :: lambda_out, theta_out
       
         REAL(KIND=dblp) :: lambda_prime, theta_prime
       
         if (.NOT. is_initialized) then
            call init_jsepar()
            is_initialized = .TRUE.
         endif
       
         lambda_prime = lambda_p(lambda_in,theta_in)
         theta_prime  = theta_p(lambda_in,theta_in)
!~          WRITE(*,*) "prime coordinates lambda, theta == ", lambda_prime*rad_to_deg, theta_prime*rad_to_deg
         lambda_out = lambda(lambda_prime,theta_prime)
         theta_out  = theta(lambda_prime,theta_prime)
       
       END SUBROUTINE test

       FUNCTION getUnifiedGrid_indx(lambda,theta,is_WW,is_degree) result(this_elmnt)
         
         REAL(KIND=dblp), INTENT(in)    :: lambda, theta
         LOGICAL, OPTIONAL, INTENT(out) :: is_WW
         LOGICAL, OPTIONAL, INTENT(in)  :: is_degree
         
         REAL(KIND=dblp)                :: lambdap, thetap
         TYPE(grid_elmnt)               :: this_elmnt         
         INTEGER(KIND=ip), DIMENSION(2) :: indexs_AA
         
         if (PRESENT(is_degree)) then
           if (is_degree) then 
              this_elmnt = which_grid(lambda*deg_to_rad,theta*deg_to_rad)
           else
              this_elmnt = which_grid(lambda,theta)
           endif
         else ! assume radians
              this_elmnt = which_grid(lambda,theta)
         endif
         
!         write(*,*) "back in Unified grid", this_elmnt%is_inWW
         if (.not. this_elmnt%is_inWW) then
            lambdap = lambda_p(lambda,theta)
            thetap  = theta_p(lambda,theta)
            indexs_AA = getAA_indexes(lambdap,thetap)
!            WRITE(*,*) "element in AA indexes == ", indexs_AA(1), indexs_AA(2)
            this_elmnt%X_idx = indexs_AA(1)
            this_elmnt%Y_idx = indexs_AA(2)
         endif

         if ( PRESENT(is_WW) ) then   
            is_WW = this_elmnt%is_inWW
         endif

       
       END FUNCTION getUnifiedGrid_indx


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   PRIVATE FUNCTIONS AND SUBROUTINES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE init_jsepar()


       INTEGER(KIND=ip) :: i, j
       REAL(KIND=dblp) :: xx, yy

!- limit (included in AA) between the 2 grids : x = 296.E - y.N, jeq<j<jsepar
          DO i=1,imax
            xx = lambda_minWW + dlambda * REAL(i-1,kind=dblp)
            yy = lambda_minAA - xx
            j = nint( (yy - theta_minWW) / dtheta ) + 1
            jsep(i) = min(max(j,jmn_AA),jmx_WW+1)
          ENDDO

       END SUBROUTINE init_jsepar
       
       FUNCTION which_grid(lambda,theta) result(this_elmnt)

         REAL(KIND=dblp), INTENT(in)    :: lambda, theta
         INTEGER(KIND=ip), DIMENSION(2) :: indxs_WW
         LOGICAL                        :: is_WW
         TYPE(grid_elmnt)               :: this_elmnt
       
         indxs_WW = getWW_indexes(lambda, theta)
!         WRITE(*,*) "Indexes obtained == ", indxs_WW(1), indxs_WW(2), jsep(indxs_WW(2))
         
         if (indxs_WW(2).GE.jsep(indxs_WW(1))) then
            is_WW = .FALSE.
         else
            is_WW = .TRUE.
         endif

!         WRITE(*,*) "Which grid == ", is_WW
         this_elmnt = grid_elmnt(X_idx = indxs_WW(1), Y_idx = indxs_WW(2), is_inWW = is_WW)
         
       END FUNCTION which_grid


       FUNCTION getWW_indexes(lambda,theta) RESULT(indexes_WW)

         REAL(KIND=dblp), INTENT(in)    :: lambda, theta
         INTEGER(KIND=ip), DIMENSION(2) :: indexes_WW
       
         indexes_WW(1) = mod(nint((mod(lambda+dblePI,dblePI)*rad_to_deg-lambda_minWW)/dlambda) + imax,imax)
         indexes_WW(2) = nint((theta*rad_to_deg-theta_minWW)/dtheta) + 1
          

       END FUNCTION getWW_indexes


       FUNCTION getAA_indexes(lambda_p,theta_p) RESULT(indexes_AA)

         REAL(KIND=dblp), INTENT(in)    :: lambda_p, theta_p
         INTEGER(KIND=ip), DIMENSION(2) :: indexes_AA
       
!         write(*,*) "Received coordsAA ==", lambda_p*rad_to_deg, theta_p*rad_to_deg
       
         indexes_AA(1) = jeq_AA - nint(theta_p*rad_to_deg/dlambda)
         indexes_AA(2) = ist_AA + nint((lambda_p*rad_to_deg-lambda_refAA)/dtheta)
          
       END FUNCTION getAA_indexes


!dmr ---  | --------------------------------------------------------- |  ---

       REAL(KIND=dblp) FUNCTION lambda_p(lambda, theta)
       
         REAL(KIND=dblp), intent(in)  :: lambda, theta
       
         REAL(KIND=dblp)              :: lambda_2
       
         lambda_2 = lambda - lambda_poleN

         lambda_p = atan2(sin(lambda_2)*cos(theta),sin(theta)) + pi_dp

       END FUNCTION lambda_p

!dmr ---  | --------------------------------------------------------- |  ---

       REAL(KIND=dblp) FUNCTION theta_p(lambda, theta)
       
         REAL(KIND=dblp), intent(in)  :: lambda, theta
         REAL(KIND=dblp)              :: lambda_2
       
         lambda_2 = lambda - lambda_poleN
  
         theta_p = -asin(cos(theta)*cos(lambda_2))
       
       END FUNCTION theta_p

!dmr ---  | --------------------------------------------------------- |  ---

       REAL(KIND=dblp) FUNCTION theta(lambda_p, theta_p)
       
         REAL(KIND=dblp), intent(in)  :: lambda_p, theta_p

         theta = asin(cos(lambda_p-pi_dp)*cos(theta_p))


       END FUNCTION theta

!dmr ---  | --------------------------------------------------------- |  ---

       REAL(KIND=dblp) FUNCTION lambda(lambda_p, theta_p)
       
         REAL(KIND=dblp), intent(in)  :: lambda_p, theta_p
              
         lambda = atan2(sin(lambda_p-pi_dp)*cos(theta_p),-sin(theta_p))+lambda_poleN
       
       END FUNCTION lambda

!dmr ---  | --------------------------------------------------------- |  ---

       END MODULE united_Grid_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
