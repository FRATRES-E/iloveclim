      SUBROUTINE out_routageEAU(iEcb,jEcb,eni,enj,direction)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!
!-----|--1--------2---------3---------4---------5---------6---------7-|
      IMPLICIT NONE

!-----|--1--------2---------3---------4---------5---------6---------7-|

      INTEGER, INTENT(IN) :: iEcb, jEcb
      INTEGER, DIMENSION(iEcb,jEcb), INTENT(IN) :: eni, enj

!-----|--1--------2---------3---------4---------5---------6---------7-|
!-----|--1--------2---------3---------4---------5---------6---------7-|
      INTEGER, DIMENSION(iEcb,jEcb), INTENT(OUT) :: direction

      REAL*4, DIMENSION(jEcb,iEcb) :: directionF

      INTEGER i,j, indx_bas
      INTEGER, PARAMETER :: fileid=777


      WRITE(*,*) "Routage fini !"

!-----|--1--------2---------3---------4---------5---------6---------7-|
!-----|--1--------2---------3---------4---------5---------6---------7-|
      direction(:,:) = 0

      DO i=1, iEcb
       DO j=1, jEcb
        if ((eni(i,j).ne.0).and.(enj(i,j).ne.0)) then

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Vers ou ? Par bassin ... 
!-----|--1--------2---------3---------4---------5---------6---------7-|

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Eurasian Arctic
!-----|--1--------2---------3---------4---------5---------6---------7-|
        if((eni(i,j).GE.29).AND.((enj(i,j).GE.3).AND.(enj(i,j).LE.34))) &
             direction(i,j) = 1
         
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       N. American Arctic
!-----|--1--------2---------3---------4---------5---------6---------7-|
        if((eni(i,j).GE.29).AND.((enj(i,j).GE.35).AND.(enj(i,j).LE.53)))&
             direction(i,j) = 2
         
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Greenland Arctic
!-----|--1--------2---------3---------4---------5---------6---------7-|
        if((eni(i,j).GE.31).AND.((enj(i,j).GE.54).AND.(enj(i,j).LE.62)))&
             direction(i,j) = 3
         
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Hudson Bay
!-----|--1--------2---------3---------4---------5---------6---------7-|

        if (                                                            &
             ((eni(i,j).GE.27).AND.(eni(i,j).LE.28))                    &
        .AND.((enj(i,j).GE.49).AND.(enj(i,j).LE.51))                    &
           )                                                            &
             direction(i,j) = 4

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Mediterranean Sea
!-----|--1--------2---------3---------4---------5---------6---------7-|
        if (                                                            &
             ((eni(i,j).GE.22).AND.(eni(i,j).LE.24))                    &
        .AND.((enj(i,j).GE.1).AND.(enj(i,j).LE.7))                      &
           )                                                            &
             direction(i,j) = 5

        if ((eni(i,j).eq.25).AND.((enj(i,j).GE.3).AND.(enj(i,j).LE.8))) &
             direction(i,j) = 5

        if ((eni(i,j).eq.24).AND.(enj(i,j).EQ.8))                       &
             direction(i,j) = 5

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Caspian Sea
!-----|--1--------2---------3---------4---------5---------6---------7-|

       if ((enj(i,j).EQ.10).AND.((eni(i,j).GE.23).AND.(eni(i,j).LE.25)))&
             direction(i,j) = 6

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       North Indian Ocean
!-----|--1--------2---------3---------4---------5---------6---------7-|

       if (((eni(i,j).GE.17).AND.(eni(i,j).LE.21)).AND.                 &
          ((enj(i,j).GE.6).AND.(enj(i,j).LE.18)))                       &
               direction(i,j) = 7

       if ((eni(i,j).EQ.22).AND.((enj(i,j).GE.8).AND.(enj(i,j).LE.18))) &
              direction(i,j) = 7

       if ((eni(i,j).EQ.17).AND.(enj(i,j).EQ.18))                       &
              direction(i,j) = 7

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       North Pacific Ocean 
!-----|--1--------2---------3---------4---------5---------6---------7-|

       if (                                                             &
           (direction(i,j).EQ.0).AND.                                   &
           ( ((eni(i,j).GE.17).AND.(eni(i,j).LE.28))                    &
        .AND.((enj(i,j).GE.18).AND.(enj(i,j).LE.47)) )                  &
          )  direction(i,j) = 8

       if (((eni(i,j).GE.17).AND.(eni(i,j).LE.19)).AND.                 &
          ((enj(i,j).GE.48).AND.(enj(i,j).LE.49)))                      &
             direction(i,j) = 8

       if (((eni(i,j).GE.17).AND.(eni(i,j).LE.18)).AND.                 &
            (enj(i,j).EQ.50) )                                          &
             direction(i,j) = 8

       if ((eni(i,j).EQ.17).AND.(enj(i,j).EQ.50))                       &
              direction(i,j) = 8

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Gulf of Mexico, Caribbean Sea
!-----|--1--------2---------3---------4---------5---------6---------7-|

       if (                                                             &
           (direction(i,j).EQ.0).AND.                                   &
           ( ((eni(i,j).GE.17).AND.(eni(i,j).LE.22))                    &
        .AND.((enj(i,j).GE.48).AND.(enj(i,j).LE.50)) )                  &
          )  direction(i,j) = 9

       if (                                                             &
           (direction(i,j).EQ.0).AND.                                   &
           ( ((eni(i,j).GE.17).AND.(eni(i,j).LE.20))                    &
        .AND.((enj(i,j).GE.51).AND.(enj(i,j).LE.54)) )                  &
          )  direction(i,j) = 9

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       North Atlantic Ocean
!-----|--1--------2---------3---------4---------5---------6---------7-|

       if (                                                             &
           (direction(i,j).EQ.0).AND.                                   &
           ( ((eni(i,j).GE.17).AND.(eni(i,j).LE.32))                    &
        .AND.((enj(i,j).GE.47).OR.(enj(i,j).LE.10)) )                   &
          )  direction(i,j) = 10

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       South Atlantic Ocean
!-----|--1--------2---------3---------4---------5---------6---------7-|

       if (((eni(i,j).GE.10).AND.(eni(i,j).LE.16)).AND.                 &
          ((enj(i,j).GE.54).OR.(enj(i,j).LE.4)))                        &
             direction(i,j) = 11

       if ((eni(i,j).LE.9).AND.                                         &
          ((enj(i,j).GE.55).OR.(enj(i,j).LE.4)))                        &
             direction(i,j) = 11

       if (((eni(i,j).GE.6).AND.(eni(i,j).LE.9)).AND.                   &
          ((enj(i,j).GE.53).AND.(enj(i,j).LE.54)))                      &
             direction(i,j) = 11

       if ((enj(i,j).EQ.54).AND.                                        &
          (eni(i,j).LE.4)       )                                       &
             direction(i,j) = 11

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       South Indian Ocean
!-----|--1--------2---------3---------4---------5---------6---------7-|

       if ((eni(i,j).LE.16).AND.                                        &
          ((enj(i,j).GE.5).AND.(enj(i,j).LE.19)))                       &
             direction(i,j) = 12

       if ((eni(i,j).LE.15).AND.                                        &
          ((enj(i,j).GE.20).AND.(enj(i,j).LE.22)))                      &
             direction(i,j) = 12

       if (((eni(i,j).GE.13).AND.(eni(i,j).LE.14)).AND.                 &
          ((enj(i,j).GE.23).AND.(enj(i,j).LE.24)))                      &
             direction(i,j) = 12

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       South Pacific: the unassigned rest ...
!-----|--1--------2---------3---------4---------5---------6---------7-|

       if (   (direction(i,j).EQ.0).AND.                                &
          (((eni(i,j).GE.1).AND.(eni(i,j).LE.16)).AND.                  &
          ((enj(i,j).GE.19).AND.(enj(i,j).LE.53)))                      &
          )                                                             &
             direction(i,j) = 13


       if (direction(i,j).EQ.0) WRITE(*,*) "not assigned", eni(i,j),     &
         &   enj(i,j)

        endif

       ENDDO
      ENDDO

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       For the output, convert it to REAL*4 variable type
!-----|--1--------2---------3---------4---------5---------6---------7-|

      do i=1,iEcb
        do j=1,jEcb
          directionF(j,i) = direction(i,j)*1.0E0
        enddo
      enddo

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Output the river basins to a runoff file...
!-----|--1--------2---------3---------4---------5---------6---------7-|

      open(fileid, CONVERT='LITTLE_ENDIAN',                             &
             file='outputdata/land/runoff_basin.dat',                   &
             form='unformatted', access='direct', status='unknown',     &
             recl=Size(directionF)*Kind(directionF(1,1)))

      write(fileid,REC=1) directionF
      close(fileid)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Create the .ctl that goes with it...
!-----|--1--------2---------3---------4---------5---------6---------7-|

      open(fileid,file='outputdata/land/runoff_basin.ctl')
      write(fileid,fmt="('dset   ^runoff_basin.dat')")
      write(fileid,fmt="('undef ',1p,e12.4)") 0.0d0
      write(fileid,fmt="('title Automatic output ECBilt var1')")
      write(fileid,fmt="('xdef ',i3,' linear ',2f7.3)") 64,0.00,5.625
      write(fileid,fmt="('ydef ',i3,' levels')") 32
      write(fileid,fmt="(' -85.7606 -80.2688 -74.7445 -69.2130 -63.6786'&
      )")
      write(fileid,fmt="(' -58.1430 -52.6065 -47.0696 -41.5325 -35.9951'&
      )")
      write(fileid,fmt="(' -30.4576 -24.9199 -19.3822 -13.8445 -8.30670'&
      )")
      write(fileid,fmt="(' -2.76890 2.76890 8.30670 13.8445 19.3822')")
      write(fileid,fmt="('  24.9199 30.4576 35.9951 41.5325 47.0696')")
      write(fileid,fmt="('  52.6065 58.1430 63.6786 69.2130 74.7445')")
      write(fileid,fmt="('  80.2688 85.7606')")
      write(fileid,fmt="('zdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
      write(fileid,fmt="('tdef ',i4,' linear 1jan0001  1YR')") 1
      write(fileid,fmt="('vars 1')")
      write(fileid,fmt="('riverbasin    1    99 River basins on ECBilt g&
      rid')")
      write(fileid,fmt="('endvars')")
      close(fileid)
      


      END SUBROUTINE out_routageEAU
