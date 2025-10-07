!
#define OFFLINE 0
! 0 run coupled
! 1 run uncoupled but use grid (nlat=32,nlon=64)
! 2 run uncoupled but single site run (nlat=1,nlon=1)
!
#define RESTART 0
! WARNING: this switch can only be implemented in offline run
!          when coupled, must be set to "0"
!          should be corrected in future
! 0 use data from initialization ("fresh" run)
! 1 use data from previously saved run
!
#define POROSITY 1
! 0 use constant porosity throughout depth
! 1 use depth porosity function
!
#define SNOW 1
! 0 no snow
! 1 snow allowed
