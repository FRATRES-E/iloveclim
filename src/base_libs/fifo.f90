! --- dmr Code initially taken from: http://computer-programming-forum.com/49-fortran/a9765f89de6aa523.htm
! --- dmr  modified to make it run since it was only a showcase but not an actual code

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      This module contains all for a FIFO to be used in and for
!       the iLOVECLIM model and more.
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche
!      Date   : 10 novembre 2013
!      Derniere modification : 20 novembre 2013, Didier M. Roche
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      MODULE FIFO_mod       !declare a module full of the FIFO queuing
                            !routines
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      IMPLICIT NONE

      Type FIFO
         integer :: head=0, tail=0, size=0   ! initialize as well as declare
         logical :: full=.true., empty=.true.
         real, allocatable, dimension(:) :: data      ! array is dynamic - can be resized
      end type FIFO

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      CONTAINS

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine q_put (queue, item)
!     assignment subroutines - all the user says is 'queue = item'

         type(FIFO) :: queue
         real item
         integer, parameter :: default_increase = 50  !initialize to 50
         real, allocatable, dimension(:) :: temp_data

         ! WRITE(*,*) "[DEBUG] Adding item to queue ... "
         if (queue%full) then
            ! WRITE(*,*) "[DEBUG] Increasing queue size ... "
            queue%size=queue%size+default_increase

            if (allocated(queue%data)) then

! Need to resize the array queue%data while keeping its contained values
! dmr That is not possible directly in FORTRAN see the MOVE_ALLOC in FORTRAN 2003 in
! IDRIS FORTRAN lessons
!
              allocate(temp_data(queue%size))
              temp_data(1:SIZE(queue%data)) = queue%data(:)
              call MOVE_ALLOC(TO=queue%data,FROM=temp_data)
            else ! not allocated so far ...
              allocate(queue%data(queue%size))
            endif
            queue%full=.false.

         endif

         queue%head = queue%head + 1
         if (queue%head.ge.queue%size) then
          queue%full=.true.
          ! WRITE(*,*) "[DEBUG] Queue is full ... "
         endif
!         if (queue%head.eq.queue%tail) queue%full=.true.

         queue%data(queue%head) = item
         queue%empty=.false.
         return

      end subroutine q_put

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function q_get (queue) ! returns next item in queue

         type(FIFO) :: queue

         real :: q_get

         if (queue%empty) then
!           CODE FOR EMPTY QUEUE REQUEST
         endif

         queue%tail = queue%tail + 1
         if (queue%tail.gt.queue%size) queue%tail = 1 ! Really? Should never happen ... dmr
         if (queue%head.eq.queue%tail) then
           queue%empty=.true.
         endif

         q_get = queue%data(queue%tail)

         ! Need code here when the queue head is after size but space
         ! exist on the left ...
         IF (queue%full) then
           call pack_data(queue)
         endif

         return

      end function q_get

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine pack_data(queue)
        type(FIFO) :: queue


        integer :: effective_size
        real, allocatable, dimension(:) :: temp_store

        effective_size = queue%head - queue%tail + 1
        if (effective_size.eq.0) then
          queue%tail = 1
          queue%head = 1
          queue%empty=.true.
          queue%full = .false.
        else
          if (effective_size.LT.queue%size) then
           ALLOCATE(temp_store(effective_size))
           temp_store(:) = queue%data(queue%tail:queue%head)
           queue%tail = 1
           queue%head = effective_size
           queue%data(queue%tail:queue%head) = temp_store(:)
           DEALLOCATE(temp_store)
           queue%full = .false.
          else
           WRITE(*,*) "[DEBUG]: NOT POSSIBLE THAT effective_size > size"
           READ(*,*)
          endif
        endif
      end subroutine pack_data

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      END MODULE FIFO_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
