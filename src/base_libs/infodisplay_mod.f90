!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!
!   Copyright 2026 FRATRES-E (https://github.com/FRATRES-E)
!     FRamework for fAst TRansient Earth-system Studies and Education

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

      module infodisplay_mod

       !! version: v1.0
       !! display: public private protected
       !! proc_internals: true


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [infodisplay_mod]

!!     @author  Didier M. Roche  (dmr)
!!     @date Creation date: April, 03rd, 2026

!!     @brief This module [infodisplay_mod] is designed to standardized communication of information to the screen user

!>
!>     DESCRIPTION : Here add the long_description of the module ...
!>        - Subroutine [NAME] : This subroutine is for blabla
!>        - Formula: $$ \frac{d\lambda}{dt} , \frac{d\phi}{dt},  \frac{dz}{dt} $$
!>     @reference References: papers or other documents to be cited... [Site the website if possible](https://iloveclim.eu)
!>     @date Last modification: SLastChangedDate$
!>     @author Last modified by : tfa, tsa!

!!     REVISION HISTORY:
!!        2026-04-03 - Initial Version
!!        TODO_YYYY-MM-DD - TODO_describe_appropriate_changes to be done or discussed - TODO_name
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       use global_constants_mod, only: str_len, stderr, ip
       use face, only: colorize

       implicit none
       private

       ! people -- PUBLIC variables, functions, subroutines.
       public :: write_em, write_im

       ! people -- Variables only required for this module
       integer(kind=ip), PARAMETER :: disp_len = 103         !! Maximum display text length
       integer(kind=ip), PARAMETER :: ori_len  = 11
       integer(kind=ip), PARAMETER :: msg_len  = 52

       ! ADDING DATE and TIME
       character(len=8)     :: date
       character(len=10)    :: time
       character(len=5)     :: zone
       integer,dimension(8) :: values

       character(len=5), parameter :: premesg="!    "
       character(len=7), parameter :: infomsg="[INFO ]", errrmsg="[ERROR]", verbmsg="[VERB+]"
       character(len=3), parameter :: sepseq=" | "
       character(len=9), parameter :: preori="CPNT <-- "
       character(len=3), parameter :: threeb="   "
       character(len=2), parameter :: postmg=" |"

       !dmr total length: 5+7+3+msg_len+3+9+ori_len+8+2 = 37 + ori_len + msg_len

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!     | [INFO ] | Mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm | CPNT <-- SSSSSSSSSSSS   18:42:06 |
!     | [ERROR] |
!     | [VERB+] |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! Source origin    : 8 characters
! Message max size : 30 characters

      contains


! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SUBROUTINE PART >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! **********************************************************************************************************************************
      SUBROUTINE write_em(message,origin)
! **********************************************************************************************************************************

!!      AUTHOR : dmr
!!      DESCRIPTION: Subroutine should be used to display error messages to the screen

       character(len=*), intent(in) :: message !! the message to be broadcasted on screen
       character(len=*), intent(in) :: origin  !! some information on the origin of the message

       ! Local variables

       character(len=disp_len) :: msgtodisplay

       ! Begin of the subroutine

       msgtodisplay = assemble_mesg(message,origin,0)

       WRITE(stderr,'(A)') colorize(msgtodisplay, color_fg='red', style='underline_on')

      END SUBROUTINE write_em
! **********************************************************************************************************************************

! **********************************************************************************************************************************
      SUBROUTINE write_im(message,origin)
! **********************************************************************************************************************************

!!      AUTHOR : dmr
!!      DESCRIPTION: Subroutine should be used to display error messages to the screen

       character(len=*), intent(in) :: message !! the message to be broadcasted on screen
       character(len=*), intent(in) :: origin  !! some information on the origin of the message

       ! Local variables

       character(len=disp_len) :: msgtodisplay

       ! Begin of the subroutine

       msgtodisplay = assemble_mesg(message,origin,1)

       WRITE(stderr,'(A)') colorize(msgtodisplay, color_fg='yellow')

      END SUBROUTINE write_im
! **********************************************************************************************************************************


! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  FUNCTIONS PART >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!===================================================================================================================================
      function assemble_mesg(message, origin, typemsg) result(theassembledmessage)
!===================================================================================================================================

!>      AUTHOR : dmr
!>      DESCRIPTION: Function to assemble the message to be displayed in a consistently formatted way

!>      Input variable :
!>         - message = the message to be displayed - will be cut to 77 characters if necessary
!>         - origin  = an info on the origin of the message
!>         - typemsg = the type of the message (error, info etc.)
!>      Output variable : the message assembled into a proper format


       character(len=*), intent(in)  :: message, origin
       integer,          intent(in)  :: typemsg
       ! Local variables

       character(len=disp_len)       :: theassembledmessage
       character(len=7)              :: disptypemsg
       character(len=msg_len)        :: messagetolength
       character(len=ori_len)        :: origintolength
       character, parameter          :: local_fillchar=" "


       integer(kind=ip)              :: msg_efflen, ori_efflen, eff_width
       ! Main body of function

       SELECT CASE (typemsg)
         CASE (0)
           disptypemsg = errrmsg
         CASE (1)
           disptypemsg = infomsg
         CASE (2)
           disptypemsg = verbmsg
       END SELECT

       call date_and_time(date,time,zone,values) ! retrieve the current date and time

       !dmr this section got some inspiration from https://github.com/eengl/fortran-strings.git

       msg_efflen = LEN_TRIM(message)
       eff_width = msg_len - msg_efflen
       if (msg_efflen.ge.msg_len) then
         messagetolength(1:msg_len) = message(1:msg_len)
       else if (eff_width.LE.3) then
         messagetolength(1:ori_len) = " "//TRIM(message)//repeat(local_fillchar,eff_width-1)
       else
         if(mod(eff_width,2).eq.0)then
            messagetolength=repeat(local_fillchar,eff_width/2)//TRIM(message)//repeat(local_fillchar,eff_width/2)
         else
            messagetolength=repeat(local_fillchar,(eff_width/2)-1)//TRIM(message)//repeat(local_fillchar,(eff_width/2)-2)
         endif
       endif

       ori_efflen = LEN_TRIM(origin)
       eff_width = ori_len - ori_efflen
       if (ori_efflen.ge.ori_len) then
         origintolength(1:ori_len) = origin(1:ori_len)
       else if (eff_width.LE.3) then
         origintolength(1:ori_len) = " "//TRIM(origin)//repeat(local_fillchar,eff_width-1)
       else
         if(mod(eff_width,2).eq.0)then
            origintolength=repeat(local_fillchar,eff_width/2)//TRIM(origin)//repeat(local_fillchar,eff_width/2)
         else
            origintolength=repeat(local_fillchar,(eff_width/2)-1)//TRIM(origin)//repeat(local_fillchar,(eff_width/2)-2)
         endif
       endif

       theassembledmessage = "" // premesg // disptypemsg // sepseq // messagetolength                               &
                               // sepseq // preori // origintolength // threeb                                       &
                              // time(1:2)//":"//time(3:4)//":"//time(5:6) // postmg

      end function assemble_mesg
!===================================================================================================================================

end module infodisplay_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
