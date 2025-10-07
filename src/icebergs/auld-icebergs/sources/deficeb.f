#include "choixcomposantes.h"

      subroutine deficeb(irunlabelclio)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  modif : 17/02/00
!  modif : 25/09/2006 JONO, APW restructuring for iceberg armada releases
!  modif : 27/07/2011 mab: added possibility of coupling GRISLI to iceberg
!                          added flag_calv_coupl
!- modif 24/10/19: cleaned use declaration sections              

      use comemic_mod, only: flag_snow_to_flux
      use global_constants_mod, only: dblp=>dp, ip             
      use const_mod, only:
      use para0_mod, only:
      use para_mod, only:
      use bloc0_mod, only: ks2, tpstot, zw, nitrun, dts, tms
      use bloc_mod, only: refexp, numit, nstart
      use ice_mod, only: imax, jmax, uwind, vwind
      use iceberg_mod, only: ndbug, nfricb, d_cotes,
     >    necriflag, nbtot_write_max, nfrqicb, icbwrout,
     >    repuls, dticeb, icbexp, linit, flag_read_positb,
     >    wind10u, wind10v, seaicevel, remsim, necriture,
     >    nbricb_plus, nbricb_plus_n, nbricb_plus_s, nbricb,
     >    icbmax, urep, vrep, flag_calv_coupl,flag_output_icb,
     >    it_icb, numtracki, filenumber,fonte_icb, vol_icb,
     >    uiceb, viceb,numiticb,nbricb_n,nbricb_s,nsource,lmx,
     >    nbtot_write_iceb, id_max_used, id_abs_table, srcx,
     >    srcy,srch,srcw,debk,berg_mk, t_lach, pond, xn, yn,
     >    hiceb, wiceb, kiceb, pond_icb, siceb, vol_orig, xn_old,
     >    yn_old, xn0, yn0,hiceb0,wiceb0,pond_icb0, nitrunicb,
     >    flag_write_iceb, fmticb, disticeb

      use reper_mod, only: xwi1, dxwi, ywj1, dywj


      implicit none

!--variables locales :
      character(len=8) :: fmtinf
      character(len=30) :: fmtflx, fmtexp
      character(len=50) :: ccfile

      real(kind=dblp) :: fonte1(500,500),fonte2(500,500),fonte0(500,500)
      logical :: flag_restart,flag_read_source
      logical :: rexist=.FALSE.
      logical :: plfound
      integer(kind=ip) :: nstart_icb,nbwires
! mab: debug
      integer(kind=ip) :: xpl1,xpl2,ypl1,ypl2, initstart, initend,
     >                 irunlabelclio, nfriqicbb = 0

!PB variables added after imposing implicit none
      integer(kind=ip) :: icb_source,icb_restart,nmax,j,i,l,iyear,ii,jj,
     >                 icdist_lost, kounter,k,nnfin,nndeb,nbrbloc,jm,
     >                 nncf,nfrc
      real(kind=dblp) :: spvicb, xvi1, yvj1

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|     
!--instructions "data" :

      write(6,*) 'deficeb : Initialisation of Icebergs'

 1000 format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
 1111 format(3(F7.1,1X,F7.3,1X),I3,A)
 1222 format (4(F7.2,1x),1x,I3,F7.2,2(1P1E12.2,1x))

! JONO_out
        open(1059,file='outputdata/icebergs/JONO.out',status='unknown')

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!      Read positb.init run parameters
!      and define setup of the run.                                          |
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-----------------------------------------------------------------------

      open (30, file="positb.init", status='old', form='formatted')
      write(6,*) 'read positb'
      read(30,*)
      read(30,*)
      read(30,*) nstart_icb
      read (30,*)
      read (30,*)
      read (30,*) ndbug, nmax, nfricb, d_cotes, necriflag
      read (30,*) nbtot_write_max,nfrqicb,icbwrout
      read (30,*)
      read (30,*) repuls, dticeb
      read (30,*)
      read (30,'(2A)') icbexp
      read (30,*)
      read (30,*) linit
      read (30,*)


!       APW comment write to mouchard
!      write(99,*) 'JA nstart_icb setup script:',nstart_icb
!      write(99,*)'icbexp',icbexp,
!     &   'reading ',linit,' bergs from positb'
!      write(99,*) 'dticeb(seconds/icb_timestep(day)):',dticeb
!      write(99,*) 'ndbug (iter betw debug writing etc):',ndbug
!      write(99,*) 'nmax (a redundant nr of iterations?):',nmax
!      write(99,*) 'nfricb (days between write?)',nfricb
!      write(99,*) 'd_cotes (minimum distance to coast):',d_cotes
!      write(99,*) 'necriflag (writing stuff (redundant)):',
!     &            necriflag
!      write(99,*) 'nbtot_write_max (maximum number of ',
!     &            'icebergs tracked):',nbtot_write_max
!      write(99,*) 'nfrqicb (old nfricb (redundant)):',nfrqicb
!      write(99,*) 'repuls:',repuls
!       APW end of commenting

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! nstart_icb = 1 : read production sources from source_iceb.param (when lacking a restart file)
! nstart_icb = 2 : read sources and icebergs from restart file resicb.om
! nstart_icb = 3 : read armada from positb.init, no iceberg sources but parameterised

      if (nstart_icb.eq.1) then
	flag_snow_to_flux = .FALSE.
	flag_read_source = .TRUE.
	flag_restart = .FALSE.
	flag_read_positb = .FALSE.
	flag_calv_coupl = .FALSE.
	flag_output_icb = .FALSE.
      endif

      if (nstart_icb.eq.2) then
	flag_snow_to_flux = .FALSE.
	flag_read_source = .FALSE.
	flag_restart = .TRUE.
	flag_read_positb = .FALSE.
	flag_calv_coupl = .FALSE.
	flag_output_icb = .FALSE.
      endif

      if (nstart_icb.eq.3) then
	flag_snow_to_flux = .TRUE.
	flag_read_source = .FALSE.
	flag_restart = .FALSE.
	flag_read_positb = .TRUE.
	flag_calv_coupl = .FALSE.
	flag_output_icb = .FALSE.
      endif

      if (nstart_icb.eq.4) then
	flag_snow_to_flux = .TRUE.
	flag_read_source = .FALSE.
	flag_restart = .TRUE.
	flag_read_positb = .FALSE.
	flag_calv_coupl = .TRUE.
	flag_output_icb = .TRUE.
      endif

      print*,'flag_calv_coupl:' ,flag_calv_coupl

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Some Initialisations of the run.                                          |
!-----------------------------------------------------------------------

!- initialisation :
         it_icb=0
         numtracki = 1
         filenumber = 1
         wind10u = 0
         wind10v = 0
         seaicevel = 0
      print*,'deficeb: ',filenumber, numtracki,wind10u,wind10v,seaicevel
!--variable en plus :
      remsim=0.1494
      necriture=0

!- origine, pt vitesse:
! JONO: what happens here?
! looks like trafo from middle of a grid to the (topleft?) corner
! scalar to vector i guess
      xvi1=xwi1-0.5*dxwi
      yvj1=ywj1-0.5*dywj
!      write(99,*)'JA write xwi1,dxwi,ywj1,dywj',xwi1,dxwi,ywj1,dywj

!-initialisation des tableaux de font et de volume de glace
!-a faire au debut de chaque troncon
! JONO oh, initialising fonte_icb is done here!
! thus, after each restart first day of fonte is missing...?
! well, it should be already put in ocean,
! and this init happens before the first day right?
      do 204 j=1,jmax
        do 204 i=1,imax
          fonte_icb(i,j)=0.
          vol_icb(i,j)=0.
 204  continue

      nbricb_plus=0
      nbricb_plus_n=0
      nbricb_plus_s=0
! JONO_nbricb lets initialise nbricb here, just to make sure
! (only relevant when not restarting, nor reading bergs from positb, thus when
! starting from source_iceb.param... (check there)

! REDUNDANT      if (.not.flag_restart) then
      nbricb=0
      do 165 l=1,icbmax
         uiceb(l)=0.
         viceb(l)=0.
165    continue
! REDUNDANT      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! JONO_says looks like calculation of vector to kick (repuls) stranded bergs
!- calcul vitesses de repulsion:

! mab - initialization of repulsion velocity
      urep=0
      vrep=0
      do 170 j=2,jmax
       do 170 i=2,imax
! mab: tms is the land mask defined in ism/inputdata->tms.dat
! mab: ks2 = 20, defined in clio/sources/defgrid that means surface layer?!
        urep(i,j)=repuls*(tms(i,j,ks2)
     &         -tms(i-1,j-1,ks2)+tms(i,j-1,ks2)-tms(i-1,j,ks2))
        vrep(i,j)=repuls*(tms(i,j,ks2)
     &         -tms(i-1,j-1,ks2)+tms(i-1,j,ks2)-tms(i,j-1,ks2))
 170  continue



!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-flag_restart: read restart file: resicb.om
!----------------------------------
      if (flag_restart) then

!       write(99,*)'flagresTRUE open restart file: resicb.om'
        inquire(file='startdata/resicb.om',exist=rexist)
        if (rexist) then
          write(*,*) 'RESICB.OM', rexist
          if (.not. flag_calv_coupl) then
            open(59,file='startdata/resicb.om',status='old',
     &      form='unformatted')
            write(99,*)'APW read restart file: resicb.om'
            read(59) numiticb,tpstot, refexp, icbexp
! mab: HERE LMX IS READ IN
            read(59)nbricb,nbricb_n,nbricb_s,nsource,lmx,nbtot_write_iceb
            id_max_used = lmx
            DO i=1,id_max_used
             id_abs_table(i) = i
            ENDDO
!           write(99,*)'nbricb_n,s,nbricb,nsource,lmx',nbricb_n,nbricb_s,
!     &            nbricb,nsource,lmx
            do 115 j=1,nsource
              read(59) srcx(j),srcy(j),srch(j),srcw(j),debk(j),berg_mk(j),
     &               t_lach(j),pond(j)
!             write(99,*) 'APW:read sources',srcx(j),srcy(j),
!     &               srch(j),srcw(j),debk(j),berg_mk(j),t_lach(j),
!     &               pond(j)
 115        continue
            do 800 l=1,lmx
              read(59) xn(l),yn(l),uiceb(l),viceb(l),hiceb(l),wiceb(l),
     &             kiceb(l),pond_icb(l)
              write(99,*) 'JA 10prcnt;read resicb.om l:',
     &              l,
     &        ' kiceb(l):',kiceb(l)
!     &               xn(l),yn(l),uiceb(l),
!     &               viceb(l),hiceb(l),wiceb(l),kiceb(l),pond_icb(l)
 800       continue
           write(99,*)'APW resicb.om has been read'

! mab: for calving!
         elseif (flag_calv_coupl) then
           open(59,file='startdata/resicb.om',status='old',
     &      form='unformatted')
           read(59) refexp,nbwires,iyear,numiticb
           write(*,*) nbwires
           do l=1,nbwires
            read(59) xn(l),yn(l),siceb(l),kiceb(l),
     &            uiceb(l),viceb(l),
     &            hiceb(l),wiceb(l),pond_icb(l)
            vol_orig(l)=wiceb(l)*wiceb(l)*1.5
     &                 *hiceb(l)*(1+remsim)*pond_icb(l)
           enddo
           lmx=nbwires
           xn_old=0d0
           yn_old=0d0
           do l=1,lmx
!mab - counting the bergs passing by the grid cells from the last run
          if(xn(l).ne.xn_old(l).or.yn(l).ne.yn_old(l)) then
            plfound = .false.
            xpl1=25.5+3.0*92
            xpl2=25.5+3.0*120
            ypl1=-81.0+3.0*29
            ypl2=-81.0+3.0*65
            if(xn(l).ge.xpl1 .and. xn(l) .le. xpl2 .and.
     &        yn(l) .ge. ypl1 .and. yn(l) .le. ypl2 ) then
              do while (.not. plfound)
              do ii=93,120
                xpl1=25.5+3.0*(ii-1)
                xpl2=25.5+3.0*(ii)
                do jj= 30,65
                  ypl1=-81.0+3.0*(jj-1)
                  ypl2=-81.0+3.0*(jj)
                  if( xn(l) .ge. xpl1 .and. xn(l) .lt. xpl2 ) then
                    if( yn(l) .ge. ypl1 .and. yn(l) .lt. ypl2 ) then
                      disticeb(ii-1,jj-1)=disticeb(ii-1,jj-1)+1.0
                      plfound = .true.
!                      write(*,*) 'disticeb: ',int(disticeb(ii-1,jj-1))
                    endif
                  endif
                enddo
              enddo
              enddo
              plfound = .false.
              xn_old(l) = xn(l)
              yn_old(l) = yn(l)
            else
             icdist_lost=icdist_lost+1
!               write(*,*) xn(l),yn(l)
            endif
          endif
         enddo   !l=1,lmx
         endif   ! flag_calv_coupl
         close(59)
       else
!         write(*,*) 'no restart file for icebergs'
         write(*,*) 'RESICB.OM', rexist
! when not reading restart file:
         nbricb_n=0
         nbricb_s=0
! JONO dont know why exactly numiticb should differ when restarting, but preserving this
       numiticb=numit
       endif
      else

! when not reading restart file:
       nbricb_n=0
       nbricb_s=0
! JONO dont know why exactly numiticb should differ when restarting, but preserving this
       numiticb=numit

      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Read rest mass from last run
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

       print*,'test calv_coupl :', flag_calv_coupl
#if (CALVFLUX == 1)
      if (flag_calv_coupl) then
         if(.not. rexist) lmx=0
         call init_calv(irunlabelclio)
      endif
#endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--- flag_read_positb
! Read iceberg armada from the rest of positb.init:
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      if (flag_read_positb) then

!- check :
        if (linit.gt.icbmax) then
          write(6,'(A,2I6)') 'Pb. ss dim: linit,icbmax=', linit, icbmax
          stop
        endif

! JONO APW instead of reading iceberg variables xn(l),..,... from positb,
! we put these vars into basic production vars xn0(1-linit),..0,..
! READ 30=positb.init icebergs
! (pond_icb is a factor to simulate a number of bergs as one particle)
        do 150 l=1,linit
          read(30,*) xn0(l),yn0(l),hiceb0(l),wiceb0(l),pond_icb0(l)
 150    continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-taille du tableau des icebergs
! lmx: number of lines in iceberg array  is set to match
! linit: number of bergs in positb.init, plus
! nbricb: number of bergs read from restart file (if it was read)

        if(.not.flag_restart) then
          lmx=0
          nbricb=0
        endif

!        write(99,*)'APW write lmx,linit,initstart,initend',
!     &     lmx,linit,initstart,initend

! JONO APW reading the vars in positb as basic production vars
        kounter=1
        do 160 l=initstart,initend
! MaB - positb.init is already read in why is that happening here
! MaB - and not just read in like that in line 227?
!          read(30,*) xn(l),yn(l),hiceb(l),wiceb(l),pond_icb(l)
      xn(l)=xn0(kounter)
	  yn(l)=yn0(kounter)
	  hiceb(l)=hiceb0(kounter)
	  wiceb(l)=wiceb0(kounter)
	  pond_icb(l)=pond_icb0(kounter)
	  kounter=kounter+1

!          write(99,*) 'APW: read icebergs',l,
!     &                xn(l),yn(l),hiceb(l),wiceb(l),pond_icb(l)
!          pond_icb(l)=1

! JONO instead of:         if(yn(l).ge.0.)nbricb_plus_n=nbricb_plus_n+1
! adding directly to the commons to avoid confusing initialising later
! not even sure if nbrplus_cum is a very usefull var
!  (it is described as the number of new bergs in a month in iceb_moy

! JONO APW doin this in iceberg.f           if(yn(l).lt.0.) then
! JONO APW doin this in iceberg.f               nbricb_s=nbricb_s+1
! JONO APW doin this in iceberg.f               nbrplus_cum_s=nbrplus_cum_s+1
! JONO APW doin this in iceberg.f           endif
! JONO APW doin this in iceberg.f           if(yn(l).ge.0.) then
! JONO APW doin this in iceberg.f               nbricb_n=nbricb_n+1
! JONO APW doin this in iceberg.f               nbrplus_cum_n=nbrplus_cum_n+1
! JONO APW doin this in iceberg.f           endif

 160    continue
        write(6,*)'end of positb.init'

      endif
      close (30)


! JONO 2006 Usefull loop to remove melted bergs? probably redundant
!     if(flaginit2.and.(.not.flaginit)) then
!      do 161 k=1,linit
!       do 202 l=1,lmx
!         if ((xn(l).eq.999).AND.(kiceb(l).gt.(-1))) kiceb(l)=(-1)
!         if (kiceb(l).eq.(-2)) then
!               kflag=l
!               read (30,*) xn(l),yn(l),hiceb(l),wiceb(l)
!               pond_icb(l)=1
!               goto 201
!          endif
!         write(6,*)'qwert'
!202    continue
!201    continue
!
!       if (kflag.eq.0) then
!         kflag=l
!         lmx=lmx+1
!       endif
!       nbricb_plus= nbricb_plus+1
!       if(yn(l).lt.0.)nbricb_plus_s=nbricb_plus_s+1
!       if(yn(l).ge.0.)nbricb_plus_n=nbricb_plus_n+1
!161   continue
!     endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!---verification distance cotes
!-------------------------------------------------------
! JONO_says d_cotes is read from positb and used by iceb_dist
! xwi, ywj are centres of gridcells, dxwi,j are grid-widths (?)

      if (flag_read_positb) then
      do 162 l=initstart,initend
          i=1+nint((xn(l)-xwi1)/dxwi)
          j=1+nint((yn(l)-ywj1)/dywj)
          call iceb_dist(l,i,j,xn(l),yn(l))
 162    continue
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! JONO_says looks like check of start position and depth of positb bergs

!- 2.verification des positions initiales
       if (flag_read_positb) then
         do 200 l=initstart,initend
          i=1+nint((xn(l)-xwi1)/dxwi)
          j=1+nint((yn(l)-ywj1)/dywj)
          write(6,*)'JA flag_read_positb;',
     & 'i,j,xn(l),yn(l),xwi1,ywj1,dxwi,dywj',
     & i,j,xn(l),yn(l),xwi1,ywj1,dxwi,dywj

           if (i.lt.1.or.i.gt.imax.or.j.lt.1.or.j.gt.jmax) then
             kiceb(l)=-1
             write(6,'(A,3I4,2F11.4)') 'Iceb. out : l,i,j,x,y=',
     &                                  l,i,j, xn(l),yn(l)
             goto 200
           else if (tms(i,j,ks2).eq.0) then
             kiceb(l)=-1
             write(6,'(A,3I4,2F11.4)') 'Iceb. Contin.: l,i,j,x,y=',
     &                                  l,i,j, xn(l),yn(l)
             goto 200
           endif

           write(99,*)'JA 10prcnt;read_positb;b4: l,lmx,initstart,initend:',
     &        l,lmx,initstart,initend,
     &        ' kiceb(l):',kiceb(l)
!-- Def. du niveau dans lequel se trouve l iceberg
            kiceb(l)=ks2
              do 120 k=ks2,2,-1
	    write(99,*)'JA l,ks2,k,kiceb(l),-zw(kiceb(l),hiceb(l)',
     & l,ks2,k,kiceb(l),-zw(kiceb(l)),hiceb(l)
                if (-zw(kiceb(l)).le.hiceb(l)) then
                 kiceb(l)=k-1
		 write(99,*)'kiceb(l)',kiceb(l)
                else
                 goto 200
                endif
 120          continue
 200     continue
            do l=1,initend
	      write(99,*) 'JA 10prcnt;read resicb.om l:',
     &              l,
     &        ' kiceb(l): ',kiceb(l)
            enddo
        endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!---flag_read_source
!---reading source_iceb.param (located in workdir?) if so flagged
!---unnecessary when using a restart file with source info!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-initiation de la generation des icebergs
!-nsricb=nombre de source d iceberg
!-srcx=position source en x,srcw=taille iceb

! JONO CHECK the j- to nsource before using this
      if (flag_read_source) then
!       write (99,*) 'reading iceberg source file source_iceb.param:'
       open (32,file="source_iceb.param", status='old')
       read (32,*) nsource
!       write (99,*) 'number of sources:',nsource
       do 163 j=1,nsource
!mab         read(32,*)srcx(j),srcy(j),srcw(j),srch(j),t_lach(j),
!mab    &              pond(j),Debk(j),berg_mk(j)
         read(32,1222) srcx(j),srcy(j),srcw(j),srch(j),t_lach(j),
     &                  pond(j),Debk(j),berg_mk(j)
! JONO check this line?:
         t_lach(j)=numiticb+t_lach(j)

!         write(99,*) 'source:',j
         write(99,*) srcx(j),srcy(j),srcw(j),srch(j),
     &                   Debk(j),berg_mk(j)
163    continue
       close(32)
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! flag_snow_to_flux:
! reset sources when using parameterised icebergs
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! nsrmax is set to 20k in iceberg  com
! setting debk s to zero just to make sure (? debk is never used??)
! (repulsion speeds for positb bergs?)

      if (flag_snow_to_flux) then

! JONO_says redundant?
        do 105 j=1,jmax
        do 105 i=1,imax
         urep(i,j) = 0.
         vrep(i,j) = 0.
 105   continue

! JONO minimising the sources to make sure CHECK division thru zero (wiceb,h) in icb_dyn and icb_traj
! INSTEAD of THIS, we will flag away gener_iceb in iceberg.f
!       do 110 j=1,nsrmax
!        srcx(j)=0.
!        srcy(j)=0.
!        srch(j)=0.00000001
!        srcw(j)=0.00000001
!        debk(j)=0.
! 110   continue

      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-- 3.2. Outputs. (JONO_says havent looked into this)
!-----------------

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! looks like allowing iceberg writing as function of iteration number
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! nfricb=30 means writing once per month?
!
! JONO_says nitrun is a common, read(10 from run.param
! clio.f: nlast=numit+nitrun, and nstart=numit+1 so nitrun is number of total iterations of the run
      nitrunicb=nitrun*(int(dts(ks2)/dticeb))
      if (nfricb.eq.0) then
        flag_write_iceb = .FALSE.
        nfricb =nitrun+1
      endif
      if (nfricb.le.nitrun) flag_write_iceb = .TRUE.


!-ouverture de icebtraj.out et nbreiceb.out
!------------------------------------------
!-entete
        open(1055,file='outputdata/icebergs.nbreiceb.out',status='unknown')
! JONO_dbg debugging writing to fort.54:
! naming, opening and closing this file (was54) to 1054=icbtraj.out in gener_iceb.f
! commenting similar lines here in deficeb.f (accesed just once)
! renaming and closing 1054 here from icbtraj.out (unformatted) to icbtraj.header (formatted)
        open(1054,file='outputdata/icebergs/icbtraj.header',status='unknown')
        write(1054,*) refexp,icbexp
! JONO_out
!        open(1059,file='outputdata/icebergs/JONO.out',status='unknown')
        write(1059,*) 'refexp,icbexp', refexp,icbexp
        write(1059,*) 'nstart_icb',nstart_icb
        write(1059,*) 'lecture ndbug,nmax,nfricb', ndbug,nmax,nfricb
! mab: causes error while debugging!!!
        write(1059,*) 'd_cotes,necriflag,nbtot_write_max,nfriqicbb',
     &            d_cotes,necriflag,nbtot_write_max,nfriqicbb

        write(6,*)'ouverture des fichiers (icbtraj en nbreiceb)'
!-calcul du nombre de blocs a ecrire
        nnfin=(nstart-1+nitrunicb)/nfricb

! JONO_says the saga continues
      nfricb=max(1,nfricb)
! JONO_end

        if (nstart.eq.1)then
          nndeb=-1
        else
          nndeb=(nstart-1)/nfricb
        endif

        nbrbloc=nnfin-nndeb
        write(1054,*) nbrbloc
! JONO_dbg see comments above
        close(1054)
        nbtot_write_iceb=0.

! JONO 06 cutting this from an obsolete nstart_icb.eq.2 (read positb bergs) loop
!-ecriture des productions sur icbtraj.out
! JONO_dbg see comments above, i think this subroutine is accessed just once
! JONO_dbg         if (flag_write_iceb) then
! JONO_dbg             write(1054) l,numiticb,real(xn(l)),real(yn(l)),
! JONO_dbg     &                  real(hiceb(l)),real(wiceb(l))
! JONO_dbg             nbtot_write_iceb = nbtot_write_iceb+1
! JONO_dbg           endif



!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! looks like allowing iceberg writing as function of iteration number END
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|



        jm=1+nitrun/nfricb
        fmticb = '(2(2F8.3,2F9.2,2F9.4))'

!----------------------------------------------------
!definition des numeros de fichiers temporaire et autres

!       nftemp1=51
!       nftemp2=52
!       nftemp3=53
!       nfterm=50
!       nfmoy=63
!       fmttrm = '(6(2F9.3))'

!----------------------------------------------------
!-temporaire pour la sortie des xys des icebergs
!-on sort un fichier contenant les xyhw sans distinction
!-et un autre contenant la taille du tableau a chaque pas.
!-un programme reconverti le tout pour utiliser courb inchange.

!       ccfile ='tempxyhw.out'
!       write(6,'(2A)') 'Ecriture sur fichier : ', ccfile
!       open(nftemp1, file=ccfile, status='unknown',form='unformatted')

!       write (nftemp1) linit
!       write (nftemp1) (xn(l),l=1,linit)
!       write (nftemp1) (yn(l),l=1,linit)
!       write (nftemp1) (hiceb(l),l=1,linit)
!       write (nftemp1) (wiceb(l),l=1,linit)
!       write (nftemp1) (uiceb(l),l=1,linit)
!       write (nftemp1) (viceb(l),l=1,linit)

!        write(nftemp1,fmticb) (xn(l),yn(l),hiceb(l)
!     &                        ,wiceb(l),uiceb(l),viceb(l),l=1,linit)

!        if ((necriflag.eq.2).or.(necriflag.eq.3)) then
!          ccfile ='temptaille.out'
!          write(6,'(2A)') 'Ecriture sur fichier : ', ccfile
!          open(nftemp2, file=ccfile, status='unknown',form=
!     &       'unformatted')

!        write(nftemp2,'(I4)')  linit

!          write (nftemp2) linit

          ccfile ='tempterme.out'

!         write(6,'(2A)') 'Ecriture sur fichier : ', ccfile
!         open(nftemp3, file=ccfile, status='unknown',form=
!    &        'unformatted')
!        endif
!---


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Ecriture de verification :
        spvicb = -99.

!     if (flaginit2) then
!--Ecriture sur fichier ASCII :
        ccfile = 'uvrep_1.fmt'
        nncf = 11
        fmtexp = '(1P125E10.2)'
        write(6,'(2A)') 'Ecriture sur fichier : ', ccfile(:nncf)
!-
        open(37, file=ccfile(:nncf), status='unknown')
        write(37,1000) fmtexp, spvicb, imax, -jmax, -2, nfrc
        write(37,1111) xvi1, dxwi, yvj1, dywj,
     &                  1., 1., 0
        write(37,'(4A)') 'Vitesse de repulsion (m/s)'
        write(37,*) 'Composante u'
        do 390 j=jmax,1,-1
          write(37,fmtexp) (urep(i,j),i=1,imax)
 390    continue
        write(37,*)
!--
        write(37,'(4A)') 'Vitesse de repulsion (m/s)'
        write(37,*) 'Composante v'
        do 395 j=jmax,1,-1
          write(37,fmtexp) (vrep(i,j),i=1,imax)
 395    continue
        write(37,*)
!--
        close(37)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-- rest of outputting was commented
!     endif
!--Fin de la partie "ecriture sur fichier".
!     open(38,file='cote.out',status='unknown')

!     do 399 j=jmax,1,-1
!         write(38,'(250I5)') (tms(i,j,ks2),i=1,imax)
!399    continue
!      close(38)

!---pour convertir tableau fonte 1.5 en 3 (driess 08/2001)

!     open (45,file='tab1.out',status='unknown')
!     open (46,file='tab2.out',status='unknown')
!     open (47,file='tab5.out',status='unknown')
!     read(45,*)
!     read(45,*)
!     read(45,*)
!     read(45,*)
!     fonte_av=0.
!     fonte_ap=0.
!     do 398  j=jmax,1,-1

!        read(45,'(1P250E8.1)') (fonte0(i,j),i=1,imax)
!398  continue
!
!      do 397 j=1,jmax
!        do 397 i=1,imax
!           fonte1(i,j)=fonte0(i,j)*
!    &            (cmx(i,j,0)*cmy(i,j,0)*dx*dy)
!           fonte_av=fonte_av+fonte1(i,j)
!397     continue

!     do 198 j=jmax,1,-1
!        write(47,'(1P250E8.1)') (fonte1(i,j),i=1,imax)
!198  continue
!
!     do 396 i=1,(imax-1),2
!        i2=int((i+1)/2)
!        j2=1
!            fonte2(i2,j2)=(fonte1(i,1)+fonte1(i+1,1))
!            fonte2(121,j2)=fonte2(1,j2)
!            fonte2(122,j2)=fonte2(2,j2)
!        i2=int((i+1)/2)
!        j2=65
!            fonte2(i2,j2)=(fonte1(i,jmax)+fonte1(i+1,jmax))
!            fonte2(121,j2)=fonte2(1,j2)
!            fonte2(122,j2)=fonte2(2,j2)
!396  continue

!     do 391 j=2,(jmax-2),2
!           j2=int((j+2)/2)
!        do 394 i=1,(imax-1),2
!           i2=int((i+1)/2)
!              fonte2(i2,j2)=(fonte1(i,j)+fonte1(i,j+1)+
!    &          fonte1(i+1,j)+fonte1(i+1,j+1))
!394     continue
!     fonte2(121,j2)=fonte2(1,j2)
!     fonte2(122,j2)=fonte2(2,j2)
!391  continue

!     do 393 j2=65,1,-1
!        write(46,'(1P125E8.1)') (fonte2(i2,j2),i2=1,122)
!        do 392 i2=1,122
!          fonte_ap=fonte_ap+fonte2(i2,j2)
!392     continue
!393  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--- Section to read fixed wind field,
! redundant since ecbilt wind is now coupled to icebergs
! JONO_icb_fix_wind
! fixed wind for iceberg forcing is read and outputed to clio3.out here
! leaving this untouched for checking purposes,
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!      write(99,*)'reading fixed wind field (just for checking purposes)'
!--- Lecture du fichier de vent en surface, units = m/s

! but [windu, windv] in icebdyn.f will no longer be a function of [uwind,vwind]
! instead now a function of ecbiltvar  [utot(i,j,3),vtot3] -> [utot10(i,j),vtot10]
! -coupling0-> [sumu10(i,j),sumv10] -ec_co2oc-> [wind10_u(ii,jj),wind10_v] (iceberg var)
!dmr @-@ iceb0
!dmr ---      open(29, file='wind_uv.om', status='old', form='unformatted')
!dmr ---      read(29)
!dmr ---      read(29)
!dmr ---      read(29)
!dmr ---      read(29)
!dmr ---      read(29) uwind
!dmr ---      read(29)
!dmr ---      read(29)
!dmr ---      read(29) vwind
!dmr ---      close(29)
!dmr @-@ iceb0
!      write(99,*) '(lat,lon),fixed Surface wind read on "wind_uv.om"'
! JONO_out
!      do 164 i=1,nlat
!        write(6,*) ((i,j,uwind(i,j),vwind(i,j)), j=1,nlon)
!      enddo
! 164  continue
      write(6,*)  uwind(10,10), vwind(10,10)

!--Lecture du fichier de vent (Annual Mean), units = m/s
!      open(28, file='vvent.om', status='old', form='unformatted')
!      read(28)
!      read(28)
!      read(28)
!      read(28)
!      read(28) vabq
!      close(28)
!      write(6,*) 'Annual Mean Wind Speed read on "vvent.om"'
!     write(6,*) 'vabq(10,10)=', vabq(10,10)

      return

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine deficeb -
      end subroutine deficeb

