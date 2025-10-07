! +--------------------+------------------------+----------------------------+
! subroutine to calculate the daily mass of calving - regarding the ice volume 
! flux that is given to this routine by GRISLI
!
! Author: Marianne Bugelmayer
! Date: 01. March 2012
! Last Modification: idem
! +--------------------+------------------------+----------------------------+
! subroutine that is called by thersf.f once a day to calculate the daily extra
! runoff
! +--------------------+------------------------+----------------------------+

      SUBROUTINE calvThersf(cday,sumday,irunlabelclio)

      USE icb_coupl, ONLY: MonIcb,calvCLIO,mass_icb, ocn_mask,RO
     &                     ,mass_day

      use const_mod

#include "type.com"
#include "para.com"
#include "bloc.com"
#include "reper.com"
#include "iceberg.com"
!#include "ice.com"

!-------------------------------------------------------------------------------+
! Subroutine arguments variables
!-------+-----------------+------------------------+----------------------------+
      INTEGER :: cday,sumday
!--------------------------------------------------------------------------------

! variables used only in this subroutine:   
      INTEGER, PARAMETER :: mon=12,kll=10
      INTEGER :: ii,jj,monend,actmon,cday1
     &           ,mass_ex1,pond_c(kll),month, testvar,lij
     &           ,num_berg(sumday),blij,calcmon,yearsec
     &           ,irunlabelclio

      REAL    :: mass_class(kll),mass_ex(kll)
     &          ,tmc,hiceb_c(kll),wiceb_c(kll)
     &          ,mass_icb1(imax,jmax)
     &          ,mass_add(imax,jmax,sumday),prct_class(kll)
     &          ,prct(mon),calvCLIO1(imax,jmax)
     &          ,mass_tot,mass_tot_icb,mass_calvgris(imax,jmax)
      CHARACTER*26 :: fname
      LOGICAL :: lexist,lexist2
      CHARACTER*6 ::cvinit 

!-------------------------------------------------------------------------------+
! definitions:
!-------+-----------------+------------------------+----------------------------+
       
       if(mod(cday,360).eq.0) testvar=int(cday/360)
       if(mod(cday,360).ne.0) testvar=int(cday/360)+1
       if(mod(cday,30).eq.0) then
         month = (int(cday/30)-(testvar-1)*12)
       else 
         month = (int(cday/30)-(testvar-1)*12)+1
       endif
! if cday=1,year =1,cday1=1,cday = 361,year=2,cday1=1,...       
       cday1 = cday-360*(testvar-1)
       calcmon =  30*(month-1)+1+360*(testvar-1)

!      mass_tot = 0d0


!--- percentage of mass that is available per month (after Martin&Adcroft 2010) 
!!      prct(1) = 0.09
!!      prct(2) = 0.12
!!      prct(3) = 0.13
!!      prct(4) = 0.18
!!      prct(5) = 0.19
!!      prct(6) = 0.06
!!      prct(7) = 0.03
!!      prct(8) = 0.01
!!      prct(9) = 0.02
!!      prct(10)= 0.04
!!      prct(11)= 0.08
!!      prct(12)= 0.05
!--- percentage of mass that is available per month
!--- (after NRC publications Archive,Report SR-2002-25, 2005) 
      prct(1) = 0.02
      prct(2) = 0.07
      prct(3) = 0.14
      prct(4) = 0.2
      prct(5) = 0.23
      prct(6) = 0.16
      prct(7) = 0.11
      prct(8) = 0.04
      prct(9) = 0.01
      prct(10)= 0.01
      prct(11)= 0.0
      prct(12)= 0.01

! yearsec=1 because volume flux from GRISLI is given in m3/yr
      yearsec = 1

!+-----------------+---------------------+-----------------------+--------------
! needed parameters:                                                            !
!+-----------------+---------------------+-----------------------+--------------+

! mab:imax,jmax in param0.com (clio/sources) defined as imax=122,jmax=65

! after the first year the initialization file won't be called anymore, so
! the calculation happens from the beginning of the 2nd year on
      if ( testvar .eq. 1 ) then
        
        if (cday .eq. 1) then
          write(cvinit,'(I6.6)') irunlabelclio-1
        
          inquire(file='startdata/mass_calvgris'//cvinit//'.dat',
     &    exist=lexist2)
          if(lexist2) then
            write(*,*) 'calvThersf, lexist2 :',lexist2,cvinit
            open(1050,file='startdata/mass_calvgris'//cvinit//'.dat',
     &      form='unformatted')
            read(1050) mass_calvgris
          else 
            write(*,*) 'calvThersf, lexist2 :',lexist2,cvinit
            mass_calvgris(:,:)=0d0
          endif
          close(1050)  
        
          calvCLIO (:,:) = calvCLIO(:,:) + mass_calvgris(:,:) 
          print*, 'calvCLIO,anfang:', sum(calvCLIO)
          print*, sum(calvCLIO)/1E6/(86400.*360.)
        endif

! to calculate the monthly available amount of ice
        if(cday .eq. calcmon) then
!      write(*,*) 'calculate monthly mass,month:',month,prct(month)
          do jj=1,jmax
            do ii=1,imax
              MonIcb(ii,jj,month)=calvCLIO(ii,jj)*RO*prct(month)
            enddo
           enddo
        endif

      else

       if(cday .eq. calcmon) then
!      write(*,*) 'calculate monthly mass,month:',month,prct(month)
         do jj=1,jmax
           do ii=1,imax
              MonIcb(ii,jj,month)=calvCLIO(ii,jj)*RO*prct(month)
           enddo
         enddo
       endif

      endif
!+-----------------+---------------------+-----------------------+--------------+
! to calculate the total daily mass of ice at each grid cell 
!+-----------------+---------------------+-----------------------+--------------+
      mass_day(:,:) = 0d0
    
      do jj=1,jmax
        do ii=1,imax
          if(MonIcb(ii,jj,month).ge.0.0)then
            mass_day(ii,jj)=(MonIcb(ii,jj,month)/30.0)
          else
            mass_day(ii,jj)=0.0
          endif !MonIcb smaller 0d0
        enddo
      enddo


      calvCLIO1 = calvCLIO

      WHERE(ocn_mask.EQ.0)
        calvCLIO1 = -99999.E0
      ENDWHERE

      if(testvar .eq. 1) then
      open(35,file='outputdata/ocean/calvmass.ctl')
      write(35,fmt="('dset   ^calvdir_bin.dat')")    
      write(35,fmt="('options big_endian')")
      write(35,fmt="('undef ',1p,e12.4)") -99999.
      write(35,fmt="('title unused mass - calving')")
      write(35,fmt="('xdef ',i3,' linear ',2f7.2)") 122,25.5,3.00
      write(35,fmt="('ydef ',i3,' linear ',2f7.2)") 65,-81.0,3.0
      write(35,fmt="('zdef ',i3,' linear ', 2i3)") 1,1,1
      write(35,fmt="('tdef ',i5,' linear 1jan0001  1YR')") testvar 
      write(35,fmt="('vars 1')") 
      write(35,fmt="('mass       0  99 unused mass   kg')")
      write(35,fmt="('endvars')")       
      close(35)
      endif

      open(15,file='outputdata/icebergs/calvdir_grads.dat'
     &,form='unformatted',status='replace')
      write(15) calvCLIO1
      close(15)

      
      if ( cday .eq. sumday) call ec_wrcalvdir_grads

!+-----------------+---------------------+-----------------------+--------------+
! end of calculations                                                           !
!+-----------------+---------------------+-----------------------+--------------+
      return

      END SUBROUTINE calvThersf

!+-----------------+---------------------+-----------------------+--------------+
!+-----------------+---------------------+-----------------------+--------------+


