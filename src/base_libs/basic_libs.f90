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
! This is version 0.5.0 of the software (see info at http://semver.org)
! 
! Contributors list is:
! Didier M. Roche / dmr / Main development, initial release
!!!

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      This module contains basic libs functions and subroutines for
!       the iLOVECLIM model and more.
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche (dmr)
!      Date   : 05 juillet 2012
!      Derniere modification : 26 novembre 2015, dmr
!-----|--1--------2---------3---------4---------5---------6---------7-|


       MODULE basic_libs


!-----|--1--------2---------3---------4---------5---------6---------7-|
!       List of module content:
!
! dmr   A method to emulate the linux/unix 'basename'
!      CHARACTER(LEN=256) FUNCTION BASENAME(pathname)
!
! dmr  Function to get the argument of some command line option
!      INTEGER FUNCTION get_argchar(longopt,shortopt,descopt,actual_arg)

! dmr A Subroutine providing a primitive progress bar
!      SUBROUTINE PROGRESS(string,ndone,ntotal)

!
!-----|--1--------2---------3---------4---------5---------6---------7-|

       IMPLICIT NONE

       CHARACTER(LEN=9), PUBLIC, PARAMETER :: emic_root = "EMIC_ROOT"

! dmr for the progress bar

       CHARACTER(LEN=255), PRIVATE, SAVE :: oldprog
       REAL(KIND=4), PRIVATE, SAVE :: oldtime

       CONTAINS


!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr   A method to emulate the linux/unix 'basename' 
!-----|--1--------2---------3---------4---------5---------6---------7-|
       CHARACTER(LEN=256) FUNCTION BASENAME(pathname)

       CHARACTER(LEN=*), INTENT(IN)  :: pathname

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr    Local variables
!-----|--1--------2---------3---------4---------5---------6---------7-|
       CHARACTER(LEN=1), PARAMETER  :: slash = "/"
       INTEGER(KIND=4) :: location

!-----|--1--------2---------3---------4---------5---------6---------7-|
! dmr   Beginning of main code
!-----|--1--------2---------3---------4---------5---------6---------7-|
       location = INDEX(pathname, slash, BACK=.TRUE., KIND=4)

       IF (location.EQ.0) THEN
         BASENAME = TRIM(pathname)
       ELSE
         BASENAME = pathname(location+1:LEN_TRIM(pathname))
       ENDIF

       END FUNCTION BASENAME

!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|

       INTEGER FUNCTION get_argchar(longopt,shortopt,descopt,actual_arg)

        IMPLICIT NONE

        CHARACTER(LEN=*), INTENT(IN) :: longopt
        CHARACTER(LEN=*), INTENT(IN) :: shortopt
        CHARACTER(LEN=*), INTENT(IN) :: descopt
        CHARACTER(LEN=256), INTENT(OUT) :: actual_arg

!dmr --- Local Variables

        INTEGER ::  narg, cptArg=1 !#of arg & counter of arg
        LOGICAL :: stop_handle = .FALSE.
        CHARACTER(LEN=256) :: l_arg

        !Check if any arguments are found
        narg=COMMAND_ARGUMENT_COUNT()

        IF (narg.GT.0) THEN

         DO WHILE((cptArg.LE.narg).AND.(.NOT.stop_handle))

         CALL GET_COMMAND_ARGUMENT(cptArg,l_arg)

         IF ((TRIM(adjustl(l_arg)).EQ."--help").OR.(TRIM(adjustl(l_arg))&
             .EQ."-h")) THEN
           WRITE(*,*) "Try using the c.l. option: "                     &
                     ,TRIM(adjustl(shortopt))
           get_argchar = -1
           RETURN
         ELSEIF(                                                        &
                (TRIM(adjustl(l_arg)).EQ.TRIM(adjustl(longopt)))        &
            .OR.(TRIM(adjustl(l_arg)).EQ.TRIM(adjustl(shortopt)))       &
               ) THEN
           cptArg = cptArg + 1 
           IF (cptArg.LE.narg) THEN
             CALL GET_COMMAND_ARGUMENT(cptArg,actual_arg)
             get_argchar = 0
           ELSE
             WRITE(*,*) "Switch: ", TRIM(adjustl(l_arg))                &
       , " should have a character string argument behind, as "         &
       , TRIM(adjustl(descopt))
             get_argchar = -1
             RETURN
           ENDIF
         ELSE
           get_argchar = -1
           WRITE(*,*) "Unkown argument: ", TRIM(adjustl(l_arg))
         ENDIF
           cptArg = cptArg + 1
         ENDDO
        ELSE
          get_argchar = -1
          WRITE(*,*) "Can't handle less than one argument"
        ENDIF

       END FUNCTION get_argchar

       SUBROUTINE PROGRESS(string,ndone,ntotal)

       IMPLICIT NONE

       character*(*) string
       character*255 prog ! ,oldprog
       double precision tl
       integer ndone,ntotal,i
!      save oldprog,oldtime

!      if (ndone.eq.0) oldtime=hires_time()
       if (ndone.eq.0) oldtime=secnds(0.0_4)
!      tl=hires_time()-oldtime
       tl=secnds(oldtime)
       if (tl.lt.0) tl=0
!     if (ndone.gt.0) tl=(1.0*ntotal/ndone)*tl-tl

!     When finished, print the total time taken rather
!     than just 00:00!
!      if (ndone.eq.ntotal) tl=secnds(oldtime)
       write(prog,'(a25,1x,''['')') string
       do i=1,40
          prog(27+i:27+i)=' '
       enddo
       write(prog(43:51),'(f7.1,''%'')') 100.0*ndone/ntotal
       do i=1,40
          if ((1.0*ndone/ntotal).gt.(1.0*i/40)) then
             if (prog(27+i:27+i).eq.' ') prog(27+i:27+i)='#'
          endif
       enddo
       prog(67:67)=']'
       write(prog(70:72),'(i2.2,'':'')')int(tl/3600)
       write(prog(73:75),'(i2.2,'':'')')int((tl-int(tl/3600)*3600)/60)
       write(prog(76:77),'(i2.2)')int((tl-int(tl/60)*60))
       if (prog.ne.oldprog) write(0,'(a,a,$)') prog(1:77),char(13)
       oldprog=prog
       if (ndone.eq.ntotal) write(0,*)
       return
       END SUBROUTINE PROGRESS

       END MODULE basic_libs
