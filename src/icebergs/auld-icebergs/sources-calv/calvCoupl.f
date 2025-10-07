! +--------------------+------------------------+----------------------------+
! subroutine to calculate the iceberg distribution (size classes, amount of 
! bergs of each class) regarding the ice volume flux that is given to this 
! routine by GRISLI
!
! Author: Marianne Bugelmayer
! Date: 09.September 2011
! Last Modification: 05.12.2011
! +--------------------+------------------------+----------------------------+

! +--------------------+------------------------+----------------------------+
! Initializations....

      SUBROUTINE init_calv(irunlabelclio)

      USE icb_coupl, ONLY: mass_rest, calvCLIO
     &                          ,mass_icb,MonIcb,RO


      use const_mod

#include "type.com"
#include "para.com"
#include "bloc.com"
#include "reper.com"
#include "iceberg.com"

      INTEGER, PARAMETER :: mon=12
      INTEGER :: irunlabelclio,wout
      REAL :: prct(mon),lon(imax,jmax),lat(imax,jmax)
     &  ,lon1(imax,jmax),lat1(imax,jmax),mass_calvgris(imax,jmax)
      LOGICAL :: lexist,lexist2
      CHARACTER*6 ::cvinit

      print*,'in init_calv'

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


! file that includes mass that hasn't been used in the previous run
! -> mass conservation

      write(cvinit,'(I6.6)') irunlabelclio

      inquire(file='startdata/mass_icb'//cvinit//'.dat',exist=lexist)
      inquire(file='startdata/mass_calvgris'//cvinit//'.dat',
     &  exist=lexist2)
!      inquire(file='startdata/mass_icb.dat',exist=lexist)
      if(lexist) then
	    write(*,*) 'lexist = TRUE'
      else
        write(*,*) 'lexist = FALSE'
      endif   
      if(lexist2) then
	    write(*,*) 'lexist2 = TRUE'
      else
        write(*,*) 'lexist2 = FALSE'
      endif

      if(lexist) then
        open(1050,file='startdata/mass_icb'//cvinit//'.dat',
!        open(1050,file='startdata/mass_icb.dat',
     &  form='unformatted')
        read(1050) mass_rest
        WHERE(mass_rest .eq. -99999.E0)
          mass_rest = 0d0
        ENDWHERE
!        do jj=1,jmax
!           write(16,*) (mass_rest(ii,jj),ii=i,imax)
!        enddo
      else 
        mass_rest(:,:)=0d0
      endif
      close(1050) 

      if(lexist2) then
        open(1050,file='startdata/mass_calvgris'//cvinit//'.dat',
     &  form='unformatted')
        read(1050) mass_calvgris
      else 
        mass_calvgris(:,:)=0d0
      endif
      close(1050)  
 
 

!+-----------------+---------------------+-----------------------+--------------
! needed parameters:                                                            !
!+-----------------+---------------------+-----------------------+--------------+
! calvCLIO in m3/yr
! mass_rest in m3

        do jj=1,jmax
           do ii=1,imax
            mass_icb(ii,jj)=(calvCLIO(ii,jj) + mass_calvgris(ii,jj))
     &         + mass_rest(ii,jj)
           enddo
        enddo
       
       open(30,file='outputdata/icebergs/init_calv.dat',status='replace')
        do jj=1,jmax
          do ii=1,imax
	        write(30,'(3F24.8)') mass_icb(ii,jj),calvCLIO(ii,jj)
     &       ,mass_calvgris(ii,jj),mass_rest(ii,jj)
          enddo
        enddo
	    close(30)

        mass_rest(:,:)=0d0


      END SUBROUTINE init_calv


! +--------------------+------------------------+----------------------------+
! subroutine that is called by iceberg.f once a day to generate the bergs 
! according to the amount of mass available (size/mass as in Bigg et al,1997)
! +--------------------+------------------------+----------------------------+

      SUBROUTINE calvCoupl(cday,sumday)

      USE icb_coupl, ONLY: MonIcb, mass_rest,startbox
     &                     ,calvCLIO,mass_icb, ocn_mask,RO

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
      INTEGER :: ii,jj,monend,actmon,cday1,sice_c(kll)
     &           ,mass_ex1,pond_c(kll),month, testvar,lij
     &           ,num_berg(2),blij,calcmon,yearsec

      REAL    :: mass_class(kll),mass_ex(kll)
     &          ,tmc,hiceb_c(kll),wiceb_c(kll)
     &          ,mass_icb1(imax,jmax)
     &          ,mass_add(imax,jmax),prct_class(kll)
     &          ,mass_day(imax,jmax),prct(mon),calvCLIO1(imax,jmax)
     &          ,mass_tot,mass_tot_icb
      CHARACTER*26 :: fname
 

!-------------------------------------------------------------------------------+
! definitions:
!-------+-----------------+------------------------+----------------------------+
c         do jj=1,jmax
c           write(16,*) (mass_icb(ii,jj),ii=i,imax)
c        enddo 

       if(mod(cday,360).eq.0) testvar=int(cday/360)
       if(mod(cday,360).ne.0) testvar=int(cday/360)+1
       if(mod(cday,30).eq.0) then
         month = (int(cday/30)-(testvar-1)*12)
       else 
         month = (int(cday/30)-(testvar-1)*12)+1
       endif
       cday1 = cday - 360*(testvar-1)
       calcmon =  30*(month-1)+1+360*(testvar-1)
!       write(*,*) 'calvCoupl:day,month,year: ',cday,month,testvar

      mass_tot = 0d0
      mass_tot_icb = 0d0
      

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

! percentage of mass depending on the size class of the berg
      prct_class(1) = 0.15
      prct_class(2) = 0.15
      prct_class(3) = 0.2
      prct_class(4) = 0.15
      prct_class(5) = 0.08
      prct_class(6) = 0.07
      prct_class(7) = 0.05
      prct_class(8) = 0.05
      prct_class(9) = 0.05
      prct_class(10) =0.05

! mass of the different classes
!     mass_class(1) = 0.473E+9
!      mass_class(2) = 3.70E+9
!      mass_class(3) = 12.6E+9
!      mass_class(4) = 29.8E+9
!      mass_class(5) = 52.3E+9
!      mass_class(6) = 75.4E+9
!      mass_class(7) = 123.0E+9
!      mass_class(8) = 177.0E+9
!      mass_class(9) = 315.0E+9
!      mass_class(10)= 492.0E+9

! height of bergs depending on classes
      hiceb_c(1) = 67
      hiceb_c(2) = 133
      hiceb_c(3) = 200
      hiceb_c(4) = 267
      hiceb_c(5) = 300
      hiceb_c(6) = 300
      hiceb_c(7) = 300
      hiceb_c(8) = 300
      hiceb_c(9) = 300
      hiceb_c(10)= 300

! width of bergs depending on classes
      wiceb_c(1) = 67
      wiceb_c(2) = 133
      wiceb_c(3) = 200
      wiceb_c(4) = 267
      wiceb_c(5) = 333
      wiceb_c(6) = 400
      wiceb_c(7) = 500
      wiceb_c(8) = 600
      wiceb_c(9) = 800
      wiceb_c(10)= 1000

! yearsec=1 because volume flux from GRISLI is given in m3/yr
      yearsec = 1

!+-----------------+---------------------+-----------------------+--------------
! needed parameters:                                                            !
!+-----------------+---------------------+-----------------------+--------------+
#if (1)
      IF(cday1 .eq. 1) THEN
       open(100,file='outputdata/ism/sumcalvcoupl.txt',
     & form='formatted',access='append')
       write(100,'(I15,5F25.10)') cday,sum(calvCLIO(:,:))
     &  ,sum(calvCLIO(:,:)*0.91)/1E6/(86400.*360.)
     &  ,sum(mass_icb(:,:)*0.91)/1E6/(86400.*360.)
       close(100)

       ENDIF
#endif     

      IF(cday1 .eq. 1) disticeb(:,:)=0d0       

! mab:imax,jmax in param0.com (clio/sources) defined as imax=122,jmax=65

! after the first year the initialization file won't be called anymore, so
! the calculation happens from the beginning of the 2nd year on

      if(testvar .gt. 1 .and. cday1 .eq. 1 .and. month .eq. 1) then
       
        do ii=1,imax
          do jj=1,jmax
             mass_tot=mass_tot+calvCLIO(ii,jj)
             mass_tot_icb=mass_tot_icb+mass_icb(ii,jj)
          enddo
        enddo
 
        open(105,file='mass_tot_calv.dat')
        write(105,*)'mass_tot calv:', mass_tot, mass_tot_icb
        close(105)
      endif

      
! to calculate the monthly available amount of ice
      if(cday .eq. calcmon) then
        do jj=1,jmax
          do ii=1,imax
            if(testvar .eq. 1) then
              MonIcb(ii,jj,month)=mass_icb(ii,jj)*prct(month)
            else
             mass_icb(ii,jj) = 0d0
             MonIcb(ii,jj,month)=(calvCLIO(ii,jj)*yearsec)*prct(month)
            endif 
          enddo
         enddo
      endif

!+-----------------+---------------------+-----------------------+--------------+
! to calculate the total daily mass of ice at each grid cell 
!+-----------------+---------------------+-----------------------+--------------+
      mass_day(:,:) = 0d0
      if (cday.eq.1) mass_add(:,:) = 0d0
      if (cday.eq.1) num_berg = 0d0
    
      do jj=1,jmax
        do ii=1,imax
          if(cday .eq. 1) then
            mass_day(ii,jj)=(MonIcb(ii,jj,month)/30.0)
          else
            mass_day(ii,jj)=(MonIcb(ii,jj,month)/30.0)+mass_add(ii,jj)
            if(mass_add(ii,jj).lt. 0.0) then           
              write(*,*) 'mass_add smaller 0!!!'
              if(MonIcb(ii,jj,month).ge.0.0)then
                mass_day(ii,jj)=(MonIcb(ii,jj,month)/30.0)
              else
               if(ii.eq.1.and.jj.eq.1) write(*,*)'MonIcb smaller 0.0!!!'
               mass_day(ii,jj)=0.0
              endif !MonIcb smaller 0d0
            endif   !mass_add gt 0d0 
          endif     !cday eq 1
        enddo
      enddo

! mab&dmr : the following call modifies lmx !!
      CALL rearrange_lmx_arrays(lmx)

! to separate the total amount of ice into the 10 size classes (Bigg,1997)
      mass_ex=0d0
! to be able to make the loops for the actual bergs as lij is the number of ALL
!the bergs produced       

! mab: HERE LMX IS MODIFIED
      if(cday.eq.1) then
        if (lmx .eq. 0) then
          id_max_used = 0
          id_abs_table(:) = 0
          lij=0
          blij=1
        else
!          write(*,*) 'in calvCoupl, lmx: ',lmx
          lij=lmx
          blij=lmx
          id_max_used = lmx
          DO i=1,id_max_used
           id_abs_table(i) = i
          ENDDO
        endif
      else      
        if(lmx.ne.0) then
          lij=lmx
          blij=lij 
          id_max_used = id_abs_table(lij)
!! mab & dmr pay attention, lmx can go to zero after some years without 
!! bergs but id_abs_table SHOULD NOT be set to zero
!! The correct fix would be to rely the first if (cday.eq.1) on
!! the first day of first year after GRISLI coupling. That is easy with
!! parameters that are set in input_timer_CPL_GRISLI.f90
         else
          blij=1
          lij=0    
          id_max_used = 0
          id_abs_table(:) = 0
         endif
      endif

      if (mod(cday,360) .eq. 0) print*,'in calvCoupl,id_max '
     &    ,id_max_used,lmx

      do jj=1,jmax
         do ii=1,imax
            if(mass_day(ii,jj) .gt. 0.0) then
              do kl=kll,1,-1
                tmc=0d0
                mass_class(kl) = 
     &         (hiceb_c(kl)*(1+remsim)*wiceb_c(kl)*1.5*wiceb_c(kl))
                if(kl.eq.kll) then
                  tmc=mass_day(ii,jj)*prct_class(kl)
                else
                  tmc=mass_day(ii,jj)*prct_class(kl)+mass_ex(kl+1)
                  mass_ex(kl+1)=0d0
                  mass_ex(kl)=0d0
                endif
                if(tmc .ge. mass_class(kl)) then
                  pond_c(kl)=aint(tmc/mass_class(kl))
                  sice_c(kl)=kl
                  if(pond_c(kl).lt.0) then
                    write(*,*)'help!pond_c,tmc,mass_class:',pond_c(kl)
     &              ,tmc,mass_class(kl), tmc/mass_class(kl)
                    write(*,*)'mass_day,mass_ex(kl+1):',mass_day(ii,jj)
     &              ,mass_ex(kl+1)
                  endif
                  mass_ex(kl)=tmc-mass_class(kl)*pond_c(kl)
                  if(mass_ex(kl).lt.0.0) write(*,*)'mass_ex smaller 0!'
                  if(mass_ex(kl).gt.mass_class(kl)) then
                    write(*,*) 'ALERT!! ',tmc,mass_ex(kl),mass_class(kl)
                  endif
                else
                  mass_ex(kl)=tmc
                  pond_c(kl)=0
                endif
                
                if(kl.eq.1.and.mass_ex(kl).gt.0.0) then
                  mass_add(ii,jj)=mass_ex(kl)
                  mass_rest(ii,jj)=mass_ex(kl)
                  mass_ex(:)=0d0 
                endif    
!                where (mass_rest .EQ. -99999.0) 
!                  mass_rest = 0d0
!                endwhere 
!                if(kl.eq.1) then
!                  mass_rest(ii,jj)=mass_add(ii,jj)
!                endif
             enddo  !kl=10,1,-1
           else
             pond_c=0
           endif
! to pass the defined size and calculated amount of icebergs of one size class
           
           do kli=1,kll 
            if(pond_c(kli).gt.0) then
! mab : lij corresponds to lmx + adding new bergs
               lij=lij+1
               siceb(lij)=sice_c(kli)
               pond_icb(lij)=pond_c(kli)
               hiceb(lij)=hiceb_c(kli)
               wiceb(lij)=wiceb_c(kli)
               xn(lij)=25.5+3.0*(ii-1)
               yn(lij)=-80.5+3.0*(jj-1)
               vol_orig(lij)= wiceb(lij)*wiceb(lij)*1.5
     &                        *hiceb(lij)*(1+remsim)*pond_icb(lij)
! mab&dmr We have a new berg, we need to generate its id
               id_max_used = id_max_used + 1 
               id_abs_table(lij) = id_max_used
               startxn(lij) = 25.5+3.0*(ii-1)
               startyn(lij) = -80.5+3.0*(jj-1)
               starthiceb(lij) = hiceb_c(kli)
               startwiceb(lij) = wiceb_c(kli)
               flagstartwritten(lij) = 0
!               startbox(ii,jj,testvar) = startbox(ii,jj,testvar) 
               startbox(ii,jj) = startbox(ii,jj) 
     &                        + 1.0
!     &                        + 1.0*real(pond_c(kli))

             endif
          enddo ! kli le 10
        enddo   ! ii=1,imax
      enddo     ! jj=1,jmax


!mab: for writing out of startposition
!           if (cday .ge. (sumday-720))
!     &       open (2055,file='outputdata/icebergs/startpos.out'
!            open (2055,file='outputdata/icebergs/startpos.out'
!     &      ,form='formatted',access='append')


! same steps as in deficeb to check the position of the icebergs
! each grid point (i,j) can have 10 types of bergs -> i*j*nr of berg

      if(lij .gt. 0) then
        do ll=blij,lij
!mab: to calculate the i,j where the iceberg is situated at the moment
          i=1+nint((xn(ll)-xwi1)/dxwi)
          j=1+nint((yn(ll)-ywj1)/dywj)
          call iceb_dist(ll,i,j,xn(ll),yn(ll))
        enddo

! taken from deficeb.f: 
! mab: icebergs that are outside of the grid are set to kiceb=-1 that means that 
! they are not taken into account in icebdyn or icebtraj and their volumen is, 
! in this case, LOST !!!
! mab: pay ATTENTION as this is bad concerning the MASS CONSERVATION -> need a SOLUTION !!!!!

        do ll=blij,lij
          if(ll.gt.0) then
            i=1+nint((xn(ll)-xwi1)/dxwi)
            j=1+nint((yn(ll)-ywj1)/dywj)
            if (i.lt.1.or.i.gt.imax.or.j.lt.1.or.j.gt.jmax) then
              kiceb(ll)=-1
              write(99,*)'ALERT!!kiceb in calvCoupl = -1-volume
     &          is lost!!i,j'
            write(*,*)'ALERT!!kiceb in calvCoupl = -1-volume is lost!!'
            write(*,*) i,nint((xn(ll)-xwi1)/dxwi),j,nint((yn(ll)-ywj1)/dywj),
     &         xn(ll),yn(ll),ll,wiceb(ll),hiceb(ll)
            else if (tms(i,j,ks2).eq.0) then
!             write(16,*)'tms,i,j,ks2,: ',tms(i,j,ks2),i,j,ks2
              kiceb(ll)=-1
              write(99,*)'ALERT!!kiceb in calvCoupl=-1-volume is
     &        lost!!-tms'
            else if( pond_icb(ll).eq. 0) then
              kiceb(ll)=-1
            endif
          endif
   
! taken from deficeb.f: 
! Definition of the niveau in the sea in which the iceberg is situated

! mab: ks2 = 20; zw is negativ and somehow related to dz, both variables (ks2, zw)
! mab: are defined in clio/sources-> defgrid.f

           kiceb(ll)=ks2
           do k=ks2,2,-1
              if (-zw(kiceb(ll)).le.hiceb(ll)) then
                kiceb(ll)=k-1
              endif
           enddo

!mab: writing out startpositions of last 2 years
!           if (cday .ge. (sumday-720)) then
!             write(2055,'(1I10,2F8.2)') id_abs_table(ll)
!     &                    ,real(startxn(ll)),real(startyn(ll))
!           endif
        enddo
      endif ! if lij gt 0

!mab: for writing out startpositions of last 2 years
!      if (cday .ge. (sumday-720)) close(2055)    
!      close(2055)    

      num_berg(2) = lij
      lmx=num_berg(2)
      num_berg(1)=lmx 
       if(num_berg(2) .ne. 0) then 
          if( mod(cday,30) .eq. 0) then
             write(*,*) 'num_berg,id_max_used,lmx: ',num_berg(2)
     &         ,id_max_used,lmx
          endif
       endif


       WHERE(ocn_mask.EQ.0)
         mass_rest = -99999.E0
         startbox = -99999.E0
       ENDWHERE
       IF (cday .eq. sumday) THEN
         ijm = Size(startbox)*Kind(startbox(1,1))
!         ijm = Size(startbox)*Kind(startbox(1,1,1))
         open(1051,file='outputdata/icebergs/startdistrib-bin.out',
     &      recl=ijm,access='direct')
         irec = 1
         write(1051,rec=irec) startbox
         close(1051)
       ENDIF

!+-----------------+---------------------+-----------------------+--------------+
! end of calculations                                                           !
!+-----------------+---------------------+-----------------------+--------------+

! create a .ctl file for grads to plat the rest mass

      if(testvar .eq. 1) then
      open(35,file='outputdata/ism/mass.ctl')
      write(35,fmt="('dset   ^massbin_grads.dat')")    
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

      open(155,file='outputdata/icebergs/mass_grads.dat'
     &,form='unformatted',status='replace')
      write(155) mass_rest
      close(155)


      return

      END SUBROUTINE calvCoupl

!+-----------------+---------------------+-----------------------+--------------+
!+-----------------+---------------------+-----------------------+--------------+

! subroutine to write out results

      SUBROUTINE ec_wrendcalv

!      implicit none

      USE icb_coupl, ONLY : mass_rest

      implicit none

#include "para.com"
#include "comunit.h" 

! mass that has to be read in at the next run of the iceberg modul
     
      write(iuo+95) mass_rest
     
      call ec_wrcalv_grads
     
      return

      END SUBROUTINE ec_wrendcalv


