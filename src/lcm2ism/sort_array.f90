!--------------------------------------------------------------------!
! afq -- module that contains a sorting routine                      !
!     -- input X is 2D, sorted on the first dimension X(1,:)         !
!     -- finding the minimum has a cost, might be improved if needed !
!--------------------------------------------------------------------!
  
module sort_array

  implicit none

  private :: FindMinimum_dbl
  private :: Swap_dbl
  public :: sort_dbl

contains

  subroutine FindMinimum_dbl(x, istart, iend, xmini)

    double precision,intent(in) :: x(:) ! the array to be sorted
    integer, intent(in) :: istart,iend  ! the range of x considered
    integer, intent(out) :: xmini       ! position of the minimum

    double precision :: minival         ! value of the minimum
    integer :: xsize                    ! size of the input array
    integer :: i                        ! loop integer
    
    xsize = ubound (x , dim = 1)

    if ((xsize.lt.istart).or.(xsize.lt.iend).or.(iend.lt.istart)) then
       write(*,*) "In sort_array: invalid dimensions"
       write(*,*) "nx,start,end:", xsize, istart, iend
       STOP
    end if

    minival  = x(istart)                ! assume the first is the min
    xmini = istart                      ! record its position
    do i = istart+1, iend               ! start with next elements
       if (x(i) < minival) then         ! if x(i) less than the min?
          minival  = x(i)               ! Yes, a new minimum found
          xmini = i                     ! record its position
       end if
    end do

    return
  end subroutine FindMinimum_dbl

  
  subroutine Swap_dbl(a,b)

    double precision, intent(inout) :: a(:), b(:)
    double precision,dimension(ubound(a,dim=1)) :: Temp

    Temp = a
    a    = b
    b    = Temp

    return
  end subroutine Swap_dbl

  
  subroutine sort_dbl(x,nxlim)

    ! afq -- we sort on the first dimension, but dim of x can be gt 2 ...
    
    double precision, intent(inout) :: x(:,:) ! the array to be sorted
    integer, intent(in) :: nxlim              ! sort on the first nxlim elmt

    integer :: xsize                    ! size of the input array
    integer :: i                        ! loop integer
    integer :: ipos                     ! position in the array
    
!    xsize = ubound (x , 1)
    
!    do i = 1, xsize-1                         ! except for the last
    do i = 1, nxlim-1                          ! except for the last
       !call FindMinimum_dbl(x(:,1), i, xsize, ipos) ! find min from this to last
       call FindMinimum_dbl(x(1,:), i, nxlim, ipos) ! find min from this to last
       call  Swap_dbl(x(:,i), x(:,ipos))            ! swap this and the minimum
    end do

    return
  end subroutine sort_dbl

  

end module sort_array


  
