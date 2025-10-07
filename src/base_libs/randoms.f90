      module randoms

      implicit none

      private

      public :: get_random_array


      logical, private :: is_seeded = .false.

      contains

      subroutine get_random_array(rand_array,naim,maim,sort)

      integer, intent(out), dimension(:) :: rand_array
      integer, intent(in), optional      :: naim, maim
      logical, intent(in), optional      :: sort

      real(kind=8)  :: u
      integer       :: i
      integer       :: n=1, m=360, ns
      integer       :: clock

      integer, dimension(:), allocatable :: seed


      if (present(naim)) n = naim
      if (present(maim)) m = maim

      if (.not. is_seeded) then

         call random_seed(size = ns)
         allocate(seed(ns))

         call system_clock(count=clock)

         seed = clock + (/ (i - 1, i = 1, ns) /)
         call random_seed(put=seed)

         deallocate(seed)

         is_seeded = .true.
      endif

      do i = lbound(rand_array,dim=1), ubound(rand_array,dim=1)

        call random_number(u)

        rand_array(i) = n + floor((m+1-n)*u)  ! We want to choose one from m-n+1 integers

      enddo

      if (present(sort) .and. sort ) then
         call quicksort(rand_array,lbound(rand_array,dim=1), ubound(rand_array,dim=1))
      endif

      end subroutine get_random_array


      ! quicksort.f -*-f90-*-
      ! Author: t-nissie
      ! License: GPLv3
      ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
      !!
      recursive subroutine quicksort(a, first, last)
        implicit none
        integer  a(*), x, t
        integer first, last
        integer i, j

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
           do while (a(i) < x)
              i=i+1
           end do
           do while (x < a(j))
              j=j-1
           end do
           if (i >= j) exit
           t = a(i);  a(i) = a(j);  a(j) = t
           i=i+1
           j=j-1
        end do
        if (first < i-1) call quicksort(a, first, i-1)
        if (j+1 < last)  call quicksort(a, j+1, last)
      end subroutine quicksort

      end module randoms
