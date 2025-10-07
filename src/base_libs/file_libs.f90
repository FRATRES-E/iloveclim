!!!
!
!This program and/or (associated) package(s) are developped within the AC²ME project, 
!  NWO project n°864.09.013
!
!Copyright 2012, D.M. Roche, didier.roche<AT>vu.nl
!Release under GPLv3 license.
!
!This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!!!

!!!
! This is version 0.2.1 of the software (see info at http://semver.org)
!
! Contributors list is:
! Didier M. Roche / dmr / Main development, initial release
!!!

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      This module contains file related libs functions and sub-
!       routines for the iLOVECLIM model and more.
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche
!      Date   : 19 novembre 2013
!      Derniere modification : 20 novembre 2013, Didier M. Roche
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       MODULE file_libs


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       List of module content:
!
! dmr   Function returning an unsued file ID
!      INTEGER FUNCTION get_fID()
!
! dmr   Function releasing a used file ID
!      SUBROUTINE release_fID()
!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       USE FIFO_mod, ONLY: FIFO, q_put, q_get

       IMPLICIT NONE

       Type fileDescriptor

         integer :: id = 0, lstr = 0, f_status = 99   ! initialize as well as declare
         integer :: lstrostat = 0 ! length of string in open_status variable
         logical :: isOpened = .false., isFormatted = .true.
         character(len=:), allocatable :: fileName ! will be allocated as
                                                   ! allocate(character(len=lstr) :: fileName)
         character(len=:), allocatable :: open_status
         integer :: rec_len = 0
       end type fileDescriptor

       INTEGER, PRIVATE, SAVE :: glo_id_last_search
       TYPE(FIFO), PRIVATE, SAVE :: pile_fID

       CONTAINS

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE open_f(f_desc,f_name,o_stat)

         character(len=*), intent(in)           :: f_name
         TYPE(fileDescriptor), intent(inout)    :: f_desc
         character(len=*), optional, intent(in) :: o_stat

         character(len=11) :: formattd


        formattd = "UNFORMATTED"

         if (f_desc%id.eq.0) then

           if (present(o_stat)) then
             call set_new_file(f_name,f_desc,oo_stat = o_stat)
           else
             call set_new_file(f_name,f_desc)
           endif
         endif

         if (.NOT.f_desc%isOpened) then

           if ( f_desc%isFormatted ) then
             formattd = "FORMATTED"
           endif
           if (f_desc%rec_len.GT.0.0) then

             OPEN(unit=f_desc%id,file=f_desc%fileName, FORM=formattd    &
                 ,iostat=f_desc%f_status,recl=f_desc%rec_len            &
                 ,status='unknown', access='direct')

           else

             OPEN(unit=f_desc%id,file=f_desc%fileName, FORM=formattd    &
                 ,iostat=f_desc%f_status,status='unknown')

           endif

           if ( f_desc%f_status.eq.0 ) then ! file is correctly opened
             f_desc%isOpened = .true.
           else ! there was an error, file is not opened ...
             ! write an error generation ??
           endif
         else
           ! do nothing, file is already opened
         endif

       END SUBROUTINE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE close_f(f_desc)

         TYPE(fileDescriptor), intent(inout) :: f_desc

         if ( f_desc%isOpened ) then
           CLOSE(unit=f_desc%id)
           call release_fID(f_desc%id)
           f_desc%id = 0
           f_desc%f_status = 99
           f_desc%isOpened = .false.
           ! I keep the file name in f_desc allocated if no formal destruction of the variable
         endif 

       END SUBROUTINE 


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE set_new_file(name_file,f_desc,oo_stat)

         character(len=*), intent(in)           :: name_file
         TYPE(fileDescriptor), intent(inout)    :: f_desc
         character(len=*), optional, intent(in) :: oo_stat

         integer :: llen
         character(len=:), allocatable :: lnamefile, lo_stat

         if (f_desc%id.ne.0) then 
           ! strange, I am called, but the variable has been assigned already
           ! design an error message? 
           WRITE(*,*) "I am called, but id already assigned ... why?"
         else
           f_desc%id = get_fID()
           llen = len(name_file)
! dmr gfortran does not handle well the allocation of character strings.
! dmr workaround suggested is the following (from: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=51055#c0)
!           allocate(character(len=llen) :: lnamefile)
           lnamefile = repeat(' ', llen)
! dmr end wrokaround
           lnamefile = adjustl(name_file) ! remove head blanks
           llen = len_trim(lnamefile)
           f_desc%lstr = llen
! dmr gfortran does not handle well the allocation of character strings.
! dmr workaround suggested is the following (from: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=51055#c0)
!           allocate(character(len=llen) :: f_desc%fileName)
           f_desc%fileName = repeat(' ', llen)
! dmr end wrokaround
           f_desc%fileName(1:f_desc%lstr) = lnamefile(1:f_desc%lstr)
           deallocate(lnamefile)

           if ( present(oo_stat) ) then

             llen = len(oo_stat)

! dmr gfortran does not handle well the allocation of character strings.
! dmr workaround suggested is the following (from: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=51055#c0)
             lo_stat = repeat(' ', llen)
! dmr end wrokaround

             lo_stat = adjustl(oo_stat) ! remove head blanks
             llen = len_trim(lo_stat)
             f_desc%lstrostat = llen

! dmr gfortran does not handle well the allocation of character strings.
! dmr workaround suggested is the following (from: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=51055#c0)
             f_desc%open_status = repeat(' ', llen)
! dmr end wrokaround
             f_desc%open_status(1:f_desc%lstrostat) = oo_stat(1:f_desc%lstrostat)

           endif

         endif

       END SUBROUTINE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Function returning an unsued file ID
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       INTEGER FUNCTION get_fID()

       IMPLICIT NONE

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr    Local variables
!-----|--1--------2---------3---------4---------5---------6---------7-|

       INTEGER, PARAMETER:: start_search = 100
       INTEGER :: l
       LOGICAL :: I_OPENED

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr   Beginning of main code
!-----|--1--------2---------3---------4---------5---------6---------7-|

       get_fID = 0

       if (.NOT. pile_fID%empty) then
         get_fID = NINT(q_get(pile_fID))

       else ! no previously released numbers ...

         l = MAX(start_search,glo_id_last_search)

         DO WHILE (get_fID.EQ.0)
           INQUIRE(l,OPENED=I_OPENED)
           IF (I_OPENED) THEN
            l = l + 1
           ELSE
            get_fID = l
           ENDIF
         ENDDO

         glo_id_last_search = get_fID

       endif

       RETURN

       END FUNCTION get_fID


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Function formally releasing a used file ID
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       SUBROUTINE release_fID(id_num)

         INTEGER id_num

         CALL q_put(pile_fID,REAL(id_num))

       END SUBROUTINE

       END MODULE file_libs
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
