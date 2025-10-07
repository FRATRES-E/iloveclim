!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of iLOVECLIM/LUDUS
!!      unix_like_libs is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
!!      as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!!
!!      unix_like_libs is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!!      of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!!
!!      You should have received a copy of the GNU General Public License along with unix_like_libs.
!!      If not, see <http://www.gnu.org/licenses/>.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: asa183_mod
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module is a wrapper for the asa183 program
!
!>     @date Creation date: April, 17th, 2020
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module asa183_mod

       implicit none
       private

       public:: rand

      ! NOTE_avoid_public_variables_if_possible_if_not_parameter
      logical :: is_seeded = .false.
      integer, parameter            :: number_ints=3
      integer(kind=selected_int_kind(8)), dimension(number_ints) :: seed

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      FUNCTION: rand
!
!>     @brief Function generating a random number, handles the setting of the seed etc.
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function rand() result(rand_num)


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       real                                :: rand_num
       integer(kind=selected_int_kind(16)) :: count, i
       character(len=1024)           :: string
       character(len=:), allocatable ::full_string
       integer                       :: len_string, step

       !dmr added for date_and_time management
       character(8)  :: date
       character(10) :: timing, reversetiming
       character(5)  :: zone
       integer :: date_i
       real :: timing_r
       integer :: zone_i
       integer,dimension(8) :: values
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       if ( .not. is_seeded ) then

         ! number here has no importance really. 
         ! modified to use a portable timing function
         ! dmr
         
         ! retrieve time with date_and_time (portable)         

         call date_and_time(DATE=date,TIME=timing,ZONE=zone,VALUES=values)
         count = 0
         read(date,'(1i8)') date_i
         read(zone,'(1i5)') zone_i
         write(reversetiming,'(*(A))') (timing(i:i),i=len(timing),1,-1)
         read(reversetiming,*) timing_r
         count=date_i*10+NINT(timing_r*10000)+zone_i

         ! call time(count)
         ! count = time()
       
         write(string,'(I0)') count
       
         len_string=LEN(trim(string))
         allocate(character(len=len_string) :: full_string)
         full_string=trim(string)
         step = len_string/(number_ints)
       
         do i=1,len_string-step,step
           read(full_string(i:i+step-1),'(I5)') seed(i/step+1)
         enddo
         
         is_seeded=.true.
      endif
   
      rand_num = r4_random(seed(1),seed(2),seed(3))
       
      end function rand

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


function r4_random ( s1, s2, s3 )

!*****************************************************************************80
!
!! R4_RANDOM returns a pseudorandom number between 0 and 1.
!
!  Discussion:
!
!    This function returns a pseudo-random number rectangularly distributed
!    between 0 and 1.   The cycle length is 6.95E+12.  (See page 123
!    of Applied Statistics (1984) volume 33), not as claimed in the
!    original article.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2008
!
!  Author:
!
!    Original FORTRAN77 original version by Brian Wichman, David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Brian Wichman, David Hill,
!    Algorithm AS 183: An Efficient and Portable Pseudo-Random
!    Number Generator,
!    Applied Statistics,
!    Volume 31, Number 2, 1982, pages 188-190.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) S1, S2, S3, three values used as the
!    seed for the sequence.  These values should be positive
!    integers between 1 and 30,000.
!
!    Output, real ( kind = 4 ) R4_RANDOM, the next value in the sequence.
!
  implicit none

  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) s3
  real ( kind = 4 ) r4_random

  s1 = mod ( 171 * s1, 30269 )
  s2 = mod ( 172 * s2, 30307 )
  s3 = mod ( 170 * s3, 30323 )
 
  r4_random = mod ( real ( s1, kind = 4 ) / 30269.0E+00 &
                  + real ( s2, kind = 4 ) / 30307.0E+00 &
                  + real ( s3, kind = 4 ) / 30323.0E+00, 1.0E+00 )

  return
end function r4_random
function r4_uni ( s1, s2 )

!*****************************************************************************80
!
!! R4_UNI returns a pseudorandom number between 0 and 1.
!
!  Discussion:
!
!    This function generates uniformly distributed pseudorandom numbers
!    between 0 and 1, using the 32-bit generator from figure 3 of
!    the article by L'Ecuyer.
!
!    The cycle length is claimed to be 2.30584E+18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2008
!
!  Author:
!
!    Original Pascal version by Pierre L'Ecuyer.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Pierre LEcuyer,
!    Efficient and Portable Combined Random Number Generators,
!    Communications of the ACM,
!    Volume 31, Number 6, June 1988, pages 742-751.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) S1, S2, two values used as the
!    seed for the sequence.  On first call, the user should initialize
!    S1 to a value between 1 and 2147483562;  S2 should be initialized
!    to a value between 1 and 2147483398.
!
!    Output, real ( kind = 4 ) R4_UNI, the next value in the sequence.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 4 ) r4_uni
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) z

  k = s1 / 53668
  s1 = 40014 * ( s1 - k * 53668 ) - k * 12211
  if ( s1 < 0 ) then
    s1 = s1 + 2147483563
  end if

  k = s2 / 52774
  s2 = 40692 * ( s2 - k * 52774 ) - k * 3791
  if ( s2 < 0 ) then
    s2 = s2 + 2147483399
  end if

  z = s1 - s2
  if ( z < 1 ) then
    z = z + 2147483562
  end if

  r4_uni = real ( z, kind = 4 ) / 2147483563.0E+00

  return
end function r4_uni
function r8_random ( s1, s2, s3 )

!*****************************************************************************80
!
!! R8_RANDOM returns a pseudorandom number between 0 and 1.
!
!  Discussion:
!
!    This function returns a pseudo-random number rectangularly distributed
!    between 0 and 1.   The cycle length is 6.95E+12.  (See page 123
!    of Applied Statistics (1984) volume 33), not as claimed in the
!    original article.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2008
!
!  Author:
!
!    FORTRAN77 original version by Brian Wichman, David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Brian Wichman, David Hill,
!    Algorithm AS 183: An Efficient and Portable Pseudo-Random
!    Number Generator,
!    Applied Statistics,
!    Volume 31, Number 2, 1982, pages 188-190.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) S1, S2, S3, three values used as the
!    seed for the sequence.  These values should be positive
!    integers between 1 and 30,000.
!
!    Output, real ( kind = 8 ) R8_RANDOM, the next value in the sequence.
!
  implicit none

  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) s3
  real ( kind = 8 ) r8_random

  s1 = mod ( 171 * s1, 30269 )
  s2 = mod ( 172 * s2, 30307 )
  s3 = mod ( 170 * s3, 30323 )
 
  r8_random = mod ( real ( s1, kind = 8 ) / 30269.0D+00 &
                  + real ( s2, kind = 8 ) / 30307.0D+00 &
                  + real ( s3, kind = 8 ) / 30323.0D+00, 1.0D+00 )

  return
end function r8_random
function r8_uni ( s1, s2 )

!*****************************************************************************80
!
!! R8_UNI returns a pseudorandom number between 0 and 1.
!
!  Discussion:
!
!    This function generates uniformly distributed pseudorandom numbers
!    between 0 and 1, using the 32-bit generator from figure 3 of
!    the article by L'Ecuyer.
!
!    The cycle length is claimed to be 2.30584E+18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2008
!
!  Author:
!
!    Original Pascal original version by Pierre L'Ecuyer
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Pierre LEcuyer,
!    Efficient and Portable Combined Random Number Generators,
!    Communications of the ACM,
!    Volume 31, Number 6, June 1988, pages 742-751.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) S1, S2, two values used as the
!    seed for the sequence.  On first call, the user should initialize
!    S1 to a value between 1 and 2147483562;  S2 should be initialized
!    to a value between 1 and 2147483398.
!
!    Output, real ( kind = 8 ) R8_UNI, the next value in the sequence.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uni
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) z

  k = s1 / 53668
  s1 = 40014 * ( s1 - k * 53668 ) - k * 12211
  if ( s1 < 0 ) then
    s1 = s1 + 2147483563
  end if

  k = s2 / 52774
  s2 = 40692 * ( s2 - k * 52774 ) - k * 3791
  if ( s2 < 0 ) then
    s2 = s2 + 2147483399
  end if

  z = s1 - s2
  if ( z < 1 ) then
    z = z + 2147483562
  end if

  r8_uni = real ( z, kind = 8 ) / 2147483563.0D+00

  return
end function r8_uni
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module asa183_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
