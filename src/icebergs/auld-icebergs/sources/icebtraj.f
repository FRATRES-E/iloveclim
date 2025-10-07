      subroutine icebtraj(nit)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  modif : 06/12/99
c that ---^ was assigning 999 to xn and n (iteration) to yn of ELIMINATEd bergs?
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c 2006,2007: modifications by JONO (and WIEA)
c (volume melt separated per layer & clean up)
c English comments are interpretations of JONO
c- modif 24/10/19: cleaned use declaration sections                            
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c in subroutine icebtraj, all icebergs with kiceb(l) > -1 MELT,
c       and their position, size, layer depth and status kiceb(l) is updated:
c	*) position is determined,
c	*) iceb_dyn is called to calculate movement of the bergs,
c       *) landed or off-the-grid icebergs are pushed to sea or ELIMINATEd,
c	*) basal and lateral MELTING is calculated and,
c	*) melted volume {m3/day} put in fonte_icb(i,j) and dVol_icb(i,j,k),
c	*) icebergs smaller then 5m. are ELIMINATEd ( <--> kiceb(l) = -1 ),
c	*) iceberg size, layer depth and status kiceb(l) is updated.
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c January 2007: clearing up kiceb(l):
c   kiceb(l) 1 to 20 (=ks2=surface layer) is iceberg depth in layer number
c   kiceb(l)=21=ks2+1 is a berg that melted completely and
c              is ready for elimination
c   kiceb(l)=0 means grounded; too shallow for keel.
c              feb 2007: bergs can no longer ground (repulsed from shallows)
c              Only new icebergs can be grounded, treat in iceberg or deficeb.
c   kiceb(l)=-1 means iceberg has been eliminated:
c		speed set to zero;
c		not allowed into icebtraj/icebdyn;
c		but still written to track.out in iceberg.f;
c		after which it is set to -2;
c		nbricb number of icebergs updated;
c		If it came from a valid ocean gridcell,
c		its volume was put into that surface layer, (CHECK new bergs)
c		else WARNING and write volume
c   kiceb(l)=-2 means the array spot is free:
c		not written anymore;
c		xn(l)=999. (--> =yn=hiceb=wiceb=uiceb=viceb in iceb_out)
c--
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- removed this, and doing that when kiceb is set to -2 in iceberg.f:
Ccmofif 24/11/99
C          yn(l) = n
Cc yn no longer records the iter.nr. of elimination (consequences not CHECKed)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c c7: march07 getting rid off xpcorr. remsim part of melt piped to surface layer
c thus make sure you dont produce grounded bergs
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip             
      use const_mod, only: untour,zero
      use para0_mod, only:
      use para_mod, only:
      use bloc0_mod, only: tms, kfs, u, v, unsdx, unsdy, zw,dz,ks2,scal
      use bloc_mod, only: smx, smy
      use ice_mod, only: imax, jmax, kmax
      use iceberg_mod, only: it_icbl, dticeb, nbricb_moins,
     >    nbricb_moins_s, nbricb_moins_n,flag_calv_coupl,dVol_icb,
     >    lmx,kiceb,wiceb,hiceb,remsim, pond_icb, uiceb,viceb,
     >    xn, yn, fonte_icb, repuls, wind10_u, wind10_v, nbricb,
     >    nbricb_plus, nbricb_n, nbricb_plus_n, nbricb_plus_s,
     >    nbricb_s,icbmax
      use reper_mod, only: xwi1, dxwi, ywj1, dywj

!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

!PB#include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "reper.com"
! [SCRPTCOM] #include "iceberg.com"
! [SCRPTCOM] #include "ice.com"

!! END_OF_INCLUDE_SECTION

#define dmrmabversion 1

c--variables locales equivalente :
c- a way to read sea temperature scal(i,j,k,1) into t(i,j,k)(?)
cdbRundo      real*8 t(imax,jmax,kmax)

!dmr [NOEQUI]        dimension t(imax,jmax,kmax)
!dmr [NOEQUI]        equivalence (scal(1,1,1,1), t(1,1,1))

c--variables locales :
      real(kind=dblp), dimension(kmax) :: un, vn

cdbR declare all as reals
c7 back to double precisions; reals made no difference for ratio fonte/dvol
      real(kind=dblp) :: deltat
      real(kind=dblp) :: ti,tkc,fbcf,flcf1,flcf2,flcf3,cfst1,cfst2
      real(kind=dblp) :: V_old(icbmax),V_new(icbmax)
      real(kind=dblp) :: volmeltday(icbmax,kmax)
      real(kind=dblp) :: vol_in(icbmax),vol_diff(icbmax)
      real(kind=dblp) :: vol_rest(imax,jmax)
c- Repulsed bergs are put at fron, d_edge degrees from the grid edge:
      real(kind=dblp) :: d_edge

c- melting variables
      real(kind=dblp) :: temp, vit, wwicb, hnew, tocean
      real(kind=dblp) :: flat, flvag, flvag_dt, fl_sum, fl_dt
      real(kind=dblp) :: dwiceb, dziceb, fb, fb_dt, fb_dt_k

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- 23-2-07 db30.CG
c- local variables for keeping track of volumes and melts:
c- manipulate writing frequency in WRITING VOLUMES loop below
c------
c- icebergs will no longer get grounded, but get repulsed from shallows
c- therefore grounded bergs can only exist at production site
c- deal with this in deficeb (and iceberg.f), but monitor just in case
c------
c- remove 'cwa' to write more warnings and alerts
c--------- INTEGERS
c- number of OOB=out off bounds bergs eliminated by icebtraj:
      integer(kind=ip) ::  i_OOB
c- number of grounded icebergs entering icebtraj (from grounded production!):
      integer(kind=ip) ::  i_grounded
c- number of rolling bergs today
      integer(kind=ip) ::  i_rolled
c--------- VOLUMES
c- total iceberg volume entering icebtraj (sum{V_old}, so not counting OOB) :
      real(kind=dblp) ::  vol_all
c- total grounded volume entering icebtraj today:
      real(kind=dblp) :: vol_grounded
      real(kind=dblp) :: vol_grounded_south
c- grounded volume after melt:
      real(kind=dblp) :: vol_grounded_after
      real(kind=dblp) :: vol_grounded_after_south
c- volume lost OOB=out of bounds is on land or off the grid
      real(kind=dblp) :: vol_lost_OOB
      real(kind=dblp) :: vol_lost_OOB_south
c- volumes of rolling bergs
      real(kind=dblp) :: old_volume, new_volume
c- volume lost by rolling bergs (small bug, repare for ane and didier)
      real(kind=dblp) :: vol_lost_roll
c--------- MELTS
c- total melted volume per gridcell, integrated over layers today:
      real(kind=dblp) :: dvol_grid(imax,jmax)
c- all melt today (integrated over grid):
      real(kind=dblp) :: dvol
c- global melt today per layer
      real(kind=dblp) :: dvol_k(ks2)
c--------- TOTAL time-integral in writing loop
c- variables will be summed (see "WRITING VOLUMES")
      real(kind=dblp) :: yr_vmt,yr_vmg,yr_vmgs
      real(kind=dblp) :: yr_vlOOB,yr_vlOOBs,yr_vlr
      real(kind=dblp) :: yr_vmtk_20,yr_vmtk_19,yr_vmtk_18,yr_vmtk_17,
     >                 yr_vmtk_16,yr_vmtk_15,yr_vmtk_14,yr_vmtk_13,
     >                 yr_vmtk_12,yr_vmtk_11,yr_vmtk_10
      real(kind=dblp) :: fonte_icbtr_dp
      integer(kind=ip) ::  yr_ioob,yr_ir,yr_ig
      logical :: vexist

!PB variables added after imposing implicit none
      integer(kind=ip) ::n,l,i,j,ksubm,k,inew,jnew,kl1,kl2,nit,
     >                   flag_write
      real(kind=dblp) ::  yr_va,xi,yj,xx,yy,xx1,yy1,htot,xs, xnew,
     >                    ynew,xfron,yfron,wwiceb,tip_of_the_iceberg,
     >                    ratio,cfstab,sum_jono3,fonte_icbtr
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      flag_write=1
      d_edge=0.000001

c      write(*,*) ' in icebtraj: ',it_icbl,it_icb
      n=it_icbl
      tkc=273.15d0
      if (kmax.eq.15) tkc=0.
c mab: dticeb = 86400
      deltat = dticeb
C      deltat = dts(ks2)

c-      ti is iceberg temperature
      ti=-4.+tkc
c-      fbcf and flcfi are basal and lateral constants for melt parameterisation
c mab: basal turbulent melting rate coefficient
      fbcf=0.58/86400
c mab:buoyant convection melting rate coefficient
      flcf1=7.62*0.001/86400
c mab:buoyant convection melting rate coefficient
      flcf2=1.29*0.001/86400
c mab: wave erosion melting rate
      flcf3=0.5/86400
c mab: coefficients to define if berg can roll over
      cfst1=0.92
      cfst2=58.32

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      nbricb_moins=0
      nbricb_moins_s=0
      nbricb_moins_n=0

c- 23-2-07 db30.CG
c--
      vol_all=0.
c-
      i_OOB=0
      vol_lost_OOB=0.
      vol_lost_OOB_south=0.
c-
      i_rolled=0
      vol_lost_roll=0.
c-
      i_grounded=0
      vol_grounded=0.
      vol_grounded_after=0.
      vol_grounded_south=0.
      vol_grounded_after_south=0.

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        if(n.eq.1) then
        write(99,*)'(first iteration) n=1',n
        yr_vlOOB=0.
        yr_ioob=0
        yr_vlOOBs=0.
        yr_vlr=0.
        yr_ir=0
        yr_va=0.
        yr_vmg=0.
        yr_vmgs=0.
        yr_vmt=0.
        fonte_icbtr_dp=0.

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c mab: to take into account the ice volume of foregoing run when using restart
         if (flag_calv_coupl) then
         inquire(file='startdata/resvol.om',exist=vexist)
         if (vexist) then
           write(*,*) 'RESVOL :', vexist
       open(59,file='startdata/resvol.om',status='old'
     &    ,form='unformatted')
           read(59) dVol_icb
           close(59)
         else
           write(*,*) 'RESVOL :', vexist
         endif
       endif
      endif !     (if n.eq.1)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-
c- main loop over all icebergs
c-
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c mab: lmx=linit (=nr of icebergs defined in positb.init)
!mab      open(250,file='volmeltday.dat',position='append')

      volmeltday(:,:) = 0d0
      V_new(:) = 0d0
      vol_in(:) = 0d0
      vol_diff(:)=0d0

      do 800 l=1,lmx
c- Negative kiceb means there is no berg to be moved or melted
c- thus TERMINATEd icebergs do not enter this routine again.
       if (kiceb(l).le.-1) goto 800



!mab       if(mod(nit,360).eq.0) then
!mab       write(250,'(1I8,3F40.15)')l,wiceb(l),hiceb(l)
!mab            endif

      vol_in(l) = wiceb(l)*wiceb(l)*1.5*hiceb(l)*(remsim+1)*pond_icb(l)

      if (vol_in(l) .le. 0d0) then
c--- ELIMINATEd:
          kiceb(l)=-1
          uiceb(l) = 0.
          viceb(l) = 0.
        goto 800
      endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- 1) position iceberg au pas de temps n+1
c     initial position (and depth) of iceberg "l"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-
        i=1+nint((xn(l)-xwi1)/dxwi)
        j=1+nint((yn(l)-ywj1)/dywj)
c- Here the (int) grid numbers of the icebergs position (xn,yn) are determined
c- feedomL.f: 1er point scalaire (1,1) lat,long en degre: xwi1=25.5, ywj1=-79.5

c---------------------------------
c- eliminate OOB bergs
        if ( (i.lt.1.or.i.gt.imax.or.j.lt.1.or.j.gt.jmax).or.
     &	     (tms(i,j,ks2).eq.0) ) then
c- db30.CG
          i_OOB=i_OOB+1
	  vol_lost_OOB=vol_lost_OOB+
     &      pond_icb(l)*(1+remsim)*hiceb(l)*wiceb(l)*1.5*wiceb(l)
          if(yn(l).lt.0.) then
            vol_lost_OOB_south=vol_lost_OOB_south+
     &        pond_icb(l)*(1+remsim)*hiceb(l)*wiceb(l)*1.5*wiceb(l)
          endif
c-
cmab	  write(99,*)'ALERT icebtraj:',
cmab     &      'icb(',l,') entered main loop with (i,j)=(',i,',',j,
cmab     &      ') OUT OF BOUNDS (land or off grid).',
cmab     &      ' ... Eliminating ... Volume lost=',
cmab     &       pond_icb(l)*(1+remsim)*hiceb(l)*wiceb(l)*1.5*wiceb(l),
cmab     &      'Total eliminated volume bergs entering icebtraj OOB=',
cmab     &       vol_lost_OOB

c--- ELIMINATEd:
          kiceb(l)=-1
	  if(i.lt.1.or.i.gt.imax.or.j.lt.1.or.j.gt.jmax) kiceb(l)=-2
          uiceb(l) = 0.
          viceb(l) = 0.
	  nbricb_moins=nbricb_moins+(1*pond_icb(l))
          if(j.le.int(jmax/2))nbricb_moins_s=nbricb_moins_s
     &                +(1*pond_icb(l))
          if(j.gt.(int(jmax/2)+1))nbricb_moins_n=nbricb_moins_n
     &                    +(1*pond_icb(l))
c--- EXIT loop, go to next berg
          goto 800
        endif
c----------------------------------

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-
        xi=xwi1+dxwi*(i-1.5)
        yj=ywj1+dywj*(j-1.5)
        xx=(xn(l)-xi)/dxwi
        yy=(yn(l)-yj)/dywj
c- (xi,yj) are the left bottom corners of the grids (i,j),
c- (presuming xwi1,ywi1 are in the center.)
c- (xx,yy) are distance (0 to 1) to the berg; inverse 'interpolation weights'
c---

        ksubm=max(kiceb(l),kfs(i,j))
c- ksubm is the (layer-)depth of iceberg bottom. kfs is ocean bottom layer.
c- grounded bergs get ksubm=kfs=ocean depth (redundant?)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c
c  2 ) Movement of iceberg "l"
c
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- 2.a) Ocean speed per layer interpolated at icb position
c-      ps.  u(.,.,.) (look like corner values: vectors!)
c- 	Calcul des vit. oc. a la posit del iceb aux diff niveaux:
c-----------------------------------------------------------------------

        do 140 k=ksubm,ks2
          yy1=1.+tms(i,j-1,k)*(tms(i,j+1,k)*yy-1.)
          xx1=1.+tms(i-1,j,k)*(tms(i+1,j,k)*xx-1.)

          un(k) = (1.-yy1)*( xx*u(i+1,j,k)+(1.-xx)*u(i,j,k) )
     &          +  yy1*( xx*u(i+1,j+1,k)+(1.-xx)*u(i,j+1,k) )
          vn(k) = (1.-yy)*( xx1*v(i+1,j,k)+(1.-xx1)*v(i,j,k) )
     &          +  yy*( xx1*v(i+1,j+1,k)+(1.-xx1)*v(i,j+1,k) )
 140    continue


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Icebergs qui touchent le fond
c- Grounded Icebergs <--> kiceb(l)=0 and speed is set to zero and SKIP to 200,
c-   where they do melt, kiceb is reevaluated at 'mise a jour', below.
c-   note: side melt affects above water part ->
c-   accentuated for grounded bergs and not counted in dVol_icb.
c
c- clean.2: changing this to avoid thick bergs superglued into shallow waters.
c instead of skipping to 200, calculate dynamics,
c then in 'check validity of new position', below
c dont allow j-movement into shallow grid; and dont allow i-movement,

        if (kiceb(l).lt.kfs(i,j)) then
CdbG          write(99,*)'WARNING icebtraj:',
CdbG     &     ' grounded iceberg detected, setting kiceb',
CdbG     &      kiceb(l),'to 0 Note: due to shallow coldness?,',
CdbG     &     ' chronically grounded bergs experience slow? melt.',
CdbG     &     ' Height of grounded berg(',l,'): (1+remsim)*hiceb=',
CdbG     &      (1+remsim)*hiceb(l), 'depth(',kfs(i,j),')=',zw(kfs(i,j))
          kiceb(l)=0
          uiceb(l)=0.
          viceb(l)=0.
c- 23-2-07 db30.CG
c total number and volume of grounded bergs entering icebtraj today
          i_grounded=i_grounded+1
          vol_grounded=vol_grounded+
     &      pond_icb(l)*1.5*wiceb(l)*wiceb(l)*hiceb(l)*(1+remsim)
          if (yn(l).lt.0) then
	    vol_grounded_south=vol_grounded_south+
     &      pond_icb(l)*1.5*wiceb(l)*wiceb(l)*hiceb(l)*(1+remsim)
	  endif
	  goto 200
        endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- 2.b) calculate acceleration of (non-stranded) iceberg l by calling icebdyn
c- 	calcul de la vitesse de l'iceberg "l" :
c-----------------------------------------------------------------------

        call icebdyn(i,j,l,un,vn,xx,yy)

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2.c) New position of iceberg "l":
c        calcul de la nouvelle position de l'iceberb "l" :
c-----------------------------------------------------------------------

        xnew=xn(l)+uiceb(l)*deltat*dxwi*unsdx*smx(i,j,0)
        ynew=yn(l)+viceb(l)*deltat*dywj*unsdy*smy(i,j,0)
        inew=1+nint((xnew-xwi1)/dxwi)
        jnew=1+nint((ynew-ywj1)/dywj)

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2.d) Check validity of new position
c-----------------------------------------------------------------------


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- DRIESS said: if bergs jump more than one gridcell, they are ELIMINATEd
c-----------------------------------------------------------------------
        if(inew.lt.1.or.inew.gt.imax.or.jnew.lt.1.or.jnew.gt.jmax)then
!mab: commented otherwise too much mouchard output
!          write(99,*)'WARNING9 icebtraj.f: '
!     &      ,((wiceb(l)*wiceb(l)*1.5*hiceb(l)*(remsim+1)*pond_icb(l))),
!     &      'is volume of FLYING? off-grid ELIMINATEd Icebergs.(',l,
!     &      ') at (inew,jnew): (',inew,',',jnew,
!     &      ') { (imax,jmax): (',imax,',',jmax,') }',
!     &      ' Volume added to surface layer of old grid-cell: ',
!     &      '(i,j): (',i,',',j,'); [x,y]=[',xn(l),',',yn(l),'] ',
!     &      'ps. iter.n=',n,' tms(i,j,ks2)=',tms(i,j,ks2),
!     &      'tms(inew,jnew,ks2)=',tms(inew,jnew,ks2),
!     &      ',H', hiceb(l),'k',kiceb(l)
!mab: end of commented part
c- put volume to fonte_ resp. dVol(surface) in old grid cell and EXIT to 800
c-mab: remsim = 0.1494
          fonte_icb(i,j)=fonte_icb(i,j)+(wiceb(l)*1.5*wiceb(l)
     &                    *1.5*hiceb(l)*(remsim+1)*pond_icb(l))
          dVol_icb(i,j,ks2)=dVol_icb(i,j,ks2)+(wiceb(l)*wiceb(l)
     &                    *1.5*hiceb(l)*(remsim+1)*pond_icb(l))
          volmeltday(l,ks2)=volmeltday(l,ks2)+(wiceb(l)*wiceb(l)
     &                    *1.5*hiceb(l)*(remsim+1)*pond_icb(l))
c----------------

c--- ELIMINATEd:
          kiceb(l)=-1
          uiceb(l) = 0.
          viceb(l) = 0.
	  nbricb_moins=nbricb_moins+(1*pond_icb(l))
          if(j.le.int(jmax/2))nbricb_moins_s=nbricb_moins_s
     &                +(1*pond_icb(l))
          if(j.gt.(int(jmax/2)+1))nbricb_moins_n=nbricb_moins_n
     &                    +(1*pond_icb(l))
c--- EXIT loop, go to next berg
          goto 800
        endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c test if bergs moved onto LAND-cells, then put at (fron= is) water edge,
c at intersection of move-direction and 'fron' and
c REPULSEd orthogonnally from fron.
c jan 2007: same for icebergs moving into waters that are too shallow
c-----------------------------------------------------------------------

        if ((tms(inew,jnew,ks2).eq.0).or.
     &      (kiceb(l).lt.kfs(inew,jnew)) ) then
c-----------------------------------------------------------------------
c---
          if (inew.lt.i) then
            xfron=xwi1+dxwi*(i-1)-0.5*dxwi
            ynew=yn(l)+
     &        ((ynew-yn(l))/(xnew-xn(l)-d_edge))*(xfron-xn(l))
            xnew=xfron+d_edge
            inew=i
            jnew=1+nint((ynew-ywj1)/dywj)
cdbug	    write(99,*)'icebtraj REPULSE grounding or stranding berg',
cdbug     &        'Arr.p.Est:old tms,new tms,old kfs,new kfs:',
cdbug     &        tms(i,j,ks2),tms(inew,jnew,ks2),kfs(i,j),kfs(inew,jnew),
cdbug     &        ',k,l,H,xnew,xfr,ynew,yn,repuls:',
cdbug     &        kiceb(l),l,hiceb(l),xnew,xfron,ynew,yn(l),repuls
            uiceb(l) = repuls
          else
            if (inew.gt.i) then
              xfron=xwi1+dxwi*(i-1)+0.5*dxwi
              ynew=yn(l)+
     &          ((ynew-yn(l))/(xnew-xn(l)+d_edge))*(xfron-xn(l))
              xnew=xfron-d_edge
              inew=i
              jnew=1+nint((ynew-ywj1)/dywj)
cdbug	    write(99,*)'icebtraj REPULSE grounding or stranding berg',
cdbug     &        'Arr.p.West:old tms,new tms,old kfs,new kfs:',
cdbug     &        tms(i,j,ks2),tms(inew,jnew,ks2),kfs(i,j),kfs(inew,jnew),
cdbug     &        ',k,l,H,xnew,xfr,ynew,yn,repuls:',
cdbug     &        kiceb(l),l,hiceb(l),xnew,xfron,ynew,yn(l),repuls
              uiceb(l) = -repuls
            endif
          endif
c---
c- xnew and ynew could be changed, so have to check again
c---
          if ((tms(inew,jnew,ks2).eq.0).or.
     &        (kiceb(l).lt.kfs(inew,jnew)) ) then
            if (jnew.lt.j) then
              yfron=ywj1+dywj*(j-1)-0.5*dywj
cdbug              write(99,*)'dbug4',(yfron-yn(l)),kiceb(l),kfs(inew,jnew),
cdbug     &          (ynew-yn(l)-d_edge),d_edge,yfron,repuls
c dbug4_jan07 -2.00000000000000         11          14  -2.52844755413137
c  epsil:1.000000000000000E-010 (setting d_edge=0.0001)
c 51.0000000000000       3.000000000000000E-003
              xnew=xn(l)+
     &          ((xnew-xn(l))/(ynew-yn(l)-d_edge))*(yfron-yn(l))
              ynew=yfron+d_edge
              jnew=j
              inew=1+nint((xnew-xwi1)/dxwi)
cdbug       	     write(99,*)'icebtraj REPULSE grounding or stranding berg',
cdbug     &        ' Arr.p.Nord:old tms,new tms,old kfs,new kfs:'
cdbug             write(99,*) tms(i,j,ks2),tms(inew,jnew,ks2)
cdbug             write(99,*) kfs(i,j),kfs(inew,jnew)
cdbug             write(99,*) 'icebergNR:',l,',k:',kiceb(l),'h:',hiceb(l)
cdbug             write(99,*) ',xnew,xn,ynew,yfron,repuls:'
cdbug             write(99,*) xnew,xn(l),ynew,yfron,repuls
               viceb(l) = repuls
            else
              if (jnew.gt.j) then
                yfron=ywj1+dywj*(j-1)+0.5*dywj
                xnew=xn(l)+
     &            ((xnew-xn(l))/(ynew-yn(l)+d_edge))*(yfron-yn(l))
                ynew=yfron-d_edge
                jnew=j
                inew=1+nint((xnew-xwi1)/dxwi)
cdbug  	      write(99,*)'icebtraj REPULSE grounding or stranding berg',
cdbug     &         ' Arr.p.Sud:old tms,new tms,old kfs,new kfs:',
cdbug     &          tms(i,j,ks2),tms(inew,jnew,ks2),kfs(i,j),kfs(inew,jnew),
cdbug     &         ',k,l,H,xnew,xn,ynew,yfron,repuls:',
cdbug     &          kiceb(l),l,hiceb(l),xnew,xn(l),ynew,yfron,repuls
                viceb(l) = -repuls
              endif
            endif
          endif
c---

c-----------------------------------------------------------------------
c- if the icebergs are still on land after repulsion,
c- (or if their old grid cell is also land;) they are ELIMINATEd,
c- 25-1-07: adding their volume to the surface layer of their old gridcell,
c-      unless that also is land; writing lost volume

          if (tms(inew,jnew,ks2).eq.0.) then
c- putting volume in old grid-cell if it was not land
            if ( (tms(i,j,ks2).eq.0.).or.
     &           (xn(l).gt.998).or.(yn(l).gt.998)) then
!mab: commented otherwise too much mouchard output
!              write(99,*)'ALERT icebtraj.f:'
!	      write(99,*)
!     &       (wiceb(l)*wiceb(l)*1.5*hiceb(l)*(remsim+1)*pond_icb(l)),
!     &       'is lost volume of ELIMINATEd iceberg(', l , ') ',
!     &       'STUCK on land from land-grid-cell (',
!     &        i, ',' ,j, '); [x,y]old =',xn(l),',',yn(l),
!     &       'ps. tms(i,j,ks2)=',tms(i,j,ks2),
!     &       'tms(inew,jnew,ks2)=',tms(inew,jnew,ks2)
!mab: commented part ends here
            else
              fonte_icb(i,j)=fonte_icb(i,j)+
     &         (wiceb(l)*wiceb(l)*1.5*hiceb(l)*(remsim+1)*pond_icb(l))
              dVol_icb(i,j,(ks2))=dVol_icb(i,j,(ks2))+
     &         (wiceb(l)*wiceb(l)*1.5*hiceb(l)*(remsim+1)*pond_icb(l))
              volmeltday(l,ks2)=volmeltday(l,ks2)+
     &         (wiceb(l)*wiceb(l)*1.5*hiceb(l)*(remsim+1)*pond_icb(l))

!mab: commented otherwise too much mouchard output
!              write(99,*)'WARNING icebtraj.f: ',
!     &       'iceberg(',l,')stuck on land;eliminated;volume was put in',
!     &       ' surface layer (ks2) of old grid-cell: (',
!     &        i, ',' ,j, '); [x,y]=',xn(l),',',yn(l),
!     &       'ps. tms(i,j,ks2)=',tms(i,j,ks2),
!     &       'tms(inew,jnew,ks2)=',tms(inew,jnew,ks2)
!mab: commented part ends here
            endif
c--- ELIMINATEd:
            kiceb(l)=-1
            uiceb(l) = 0.
            viceb(l) = 0.
	    nbricb_moins=nbricb_moins+(1*pond_icb(l))
            if(j.le.int(jmax/2))nbricb_moins_s=nbricb_moins_s
     &               +(1*pond_icb(l))
            if(j.gt.(int(jmax/2)+1))nbricb_moins_n=nbricb_moins_n
     &               +(1*pond_icb(l))

c            goto 800
c---

c---
c- if not tms(inew,jnew,ks2).eq.0.,
c- stranding/grounding berg stopped at front and repulsed {read from positb}
          else
	      xn(l)=xnew
            yn(l)=ynew
          endif

c---
c- warning for chronicly grounded bergs
          if (kiceb(l).lt.kfs(inew,jnew)) then
!mab: commented otherwise too much mouchard output
!	    write(99,*) 'WARNING icebtraj:',
!     &        'grounding iceberg still grounded at grid front',
!     &        'note: possible slow melt',
!     &        ': k:',kiceb(l),' i,inew;j,jnew:',i,inew,j,jnew,
!     &        tms(i,j,ks2),tms(inew,jnew,ks2),kfs(i,j),kfs(inew,jnew),
!     &        'l,H,xnew,ynew at front:',l,hiceb(l),xnew,ynew
!mab: commented part ends here
            kiceb(l)=0
	    uiceb=0.
	    viceb=0.
          endif
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c------ else bergs were not stranding on land-cells or grounding in shallows

        else
          xn(l)=xnew
          yn(l)=ynew
        endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- 2.f)	untour=360d0, so this must be the circular boundary condition.
c----------------------------------------------------------------------

        xn(l)=xwi1+dxwi+mod(xn(l)-(xwi1+dxwi)+untour,untour)
C        xn(l)=xwi1+0.5*dxwi+mod(xn(l)-(xwi1+0.5*dxwi)+untour,untour)


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

 200  continue

c from here on chronicly grounded bergs (kiceb(l)=0) are also treated
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3) 	Desintegration et fonte de l iceberg "l".
c 	iceberg MELT and DiSINTEGRATION
c-----------------------------------------------------------------------
c JONO quoting Bigg: "neglecting the significant train
c                     (envelope) of cold meltwater will accentuate melting"
c => The lateral melt only affecting the 'front' sides of the berg
c    can be seen as a (crude) cold melt water train!
c    Waves can also be thought to be weakened at the 'back' of the berg...
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c check for zero hiceb, because of divisions:
        if(hiceb(l).le.0.) then
	  write(99,*)'ALERT icebtraj; hiceb<=0'
	  goto 790
	endif
c---

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- 3.1)	calcul du volume de l iceberg "l"
c- initial volume of iceberg l (side is 1.5*front...)
c----------------------------------------------------------------------
c to calculate fonte_icb, new iceberg volume is subtracted from this initial vol
        V_old(l)=wiceb(l)*wiceb(l)*1.5*hiceb(l)*(remsim+1)

c db30.CG total volume of live (melting) icebergs
        vol_all=vol_all+V_old(l)


c- 3.2a) BASAL MELT
c----------------------------------------------------------------------
c basal melt based on ice-shelf parameterisation (3.17: Weeks and Campbell 1973)
c Rb = 0.58 |Vw-Vi|^0.8 (Tw-Ti) / L^0.2 (m/day)
c
cmab: why **0.4 if above the equation is with 0.8???
!dmr [NOEQUI]        temp=(t(i,j,ksubm)-ti)
        temp=(scal(i,j,ksubm,1)-ti)
        vit=((un(ksubm)-uiceb(l))*(un(ksubm)-uiceb(l))+
     &       (vn(ksubm)-viceb(l))*(vn(ksubm)-viceb(l)))**0.4
c- un(k) is layer k waterspeed, interpolated to icb position (see above)
        wwiceb=wiceb(l)*1.5

!        write(250,*)'fb :'
!        write(250,'(5F20.10,I5)')vit,un(ksubm),uiceb(l),vn(ksubm),viceb(l),ksubm
cdmr&mab --- I do not understand why the following equation pertains to wwiceb only (dmr)
        if(wwiceb .gt. 0d0) then
          fb=fbcf*vit*temp/(wwiceb**0.2)
        else
          fb=0d0
        endif
cdbR !!!
c        fb=0.
cdbR-

c-- no negative melt:
cdbr comment: once NEGATIVE in mouchard(anes non-grounding test positb)
        if (fb.lt.0.) then
!mab: commented otherwise too much mouchard output
!         write(99,*) 'WARNING icebtraj.f; ',
!     &      ' daily basal melt fb=',fb*deltat,'NEGATIVE icb l=',l,
!     &      ' Ocean Temperature',t(i,j,ksubm),
!     &      'at iceberg-bottom depth layer ksubm:', ksubm,
!     &      'in grid-cell: (' ,i, ',' ,j, '). Reset fb to ZERO '
!mab: commented part ends here
          fb=0.
        endif
c--
        fb_dt=fb*deltat

c correct for melt exceeding iceberg-size
        if ((fb_dt).gt.((1+remsim)*hiceb(l))) then
!mab: commented otherwise too much mouchard output
!          write(99,*) 'WARNING icebtraj.f; Corrected for ',
!     &    'basal melt exceeding berg(',l,')s ',
!     &    'total height:',hiceb(l)*(1+remsim),
!     &    ' at (i,j)=(',i,',',j,'iteration:',n,
!     &    ' Daily basal melt=',fb*deltat,
!     &    ' icb-velocity(u,v):',uiceb(l),viceb(l),' SST:',t(i,j,ksubm)
!mab: commented part ends here
          fb_dt=(1+remsim)*hiceb(l)
        endif

cdbr
cdbr comment: not true before it 61...
c7        if ((fb_dt).gt.hiceb(l)) then
c7          write(99,*) 'cdbr WARNING icebtraj.f; fb_dt gt hiceb ',
c7     &    'basal melt exceeding berg(',l,')s ',
c7     &    'submarine height:',hiceb(l)*(1+remsim),
c7     &    ' at (i,j)=(',i,',',j,'iteration:',n,
c7     &    ' Daily basal melt=',fb*deltat,
c7     &    ' icb-velocity(u,v):',uiceb(l),viceb(l),' SST:',t(i,j,ksubm)
c7        endif
cdbr-

c- hnew is used to find the new, buoyancy adjusted, height. see NOTE below
c CHECKed        hnew=hiceb(l)-(fb*deltat)
         hnew=(1+remsim)*hiceb(l)-fb_dt
cdbR
c7        write(99,*)'fb_dt',fb_dt,
c7     &    '(1+r)hiceb(',l,'):',(1+remsim)*hiceb(l),'hnew:',hnew
c7        fb_dvol=1.5*wiceb(l)*wiceb(l)*fb_dt
c7        fb_dvol_sum=fb_dvol_sum+fb_dvol
c7        write(99,*)'fb_dvol=surface_area*fb_dt=',fb_dvol,
c7     &    'cumulative:',fb_dvol_sum
cdbR-

c--- 3.2b) Lateral wave erosion in surface layer is a function of sea-state
c-         It is extrapolated over total iceberg height (1+remsim)hiceb
c----------------------------------------------------------------------
c       Rvague = 0.5 Ss    (m/day!!) , with
c       Ss = -5 + (32+2|Va|)^0.5 (=1 if Va=0)
c ("Both above and below waterline")
c Bigg et al '97, quoting El-Tahan et al '87
c-
        xs = -5 + sqrt( 32+2*
     &  sqrt(wind10_u(i,j)*wind10_u(i,j)+wind10_v(i,j)*wind10_v(i,j)))
        flvag=flcf3*xs
cdbR
c        flvag=0.
cdbR

c- we2top wave erosion to top layer
c- separate flvag_dt (local var in this k-loop; dont touch flvag),
c- but keep fl_sum->fl_dt for compensation purposes
          flvag_dt=flvag*deltat
c- we2top end

cdbr
cdbr comment:
        if ((flvag_dt).lt.0.) then
          write(99,*) 'cdbr WARNING icebtraj.f; flvag_dt negative'
!mab: commented otherwise too much mouchard output
!     &    ,'berg(',l,')s h,w:',hiceb(l),wiceb(l),
!     &    ' at (i,j)=(',i,',',j,'iteration:',n,
!     &    ' Daily wave melt= flvag_dt',flvag_dt,
!     &    'times total heigt:',hiceb(l)*(1+remsim),'=',
!     &     flvag_dt*hiceb(l)*(1+remsim),
!     &    ' icb-velocity(u,v):',uiceb(l),viceb(l),' SST:',t(i,j,ksubm)
!mab: commented part ends here
cdmr&mab --- I do not understand why flvag_dt < 0.0 != 0.0 (dmr)
cdmr&mab --- Let's do it!
           flvag_dt = 0.0d0
        endif

cdbr-
cdbr
cdbr comment:
        if ((flvag_dt).gt.wiceb(l)) then
          write(99,*) 'cdbr WARNING icebtraj.f; flvag_dt>wiceb'
!mab: commented otherwise too much mouchard output
!     &    ,'berg(',l,')s h,w:',hiceb(l),wiceb(l),
!     &    ' at (i,j)=(',i,',',j,'iteration:',n,
!     &    ' Daily wave melt= flvag_dt',flvag_dt,
!     &    'times total heigt:',hiceb(l)*(1+remsim),'=',
!     &     flvag_dt*hiceb(l)*(1+remsim),
!     &    ' icb-velocity(u,v):',uiceb(l),viceb(l),' SST:',t(i,j,ksubm)
!mab: commented part ends here
cdmr&mab --- I do not understand why flvag_dt > wiceb(l)  != wiceb(l) (dmr)
cdmr&mab --- Let's do it!
           flvag_dt = wiceb(l)
        endif


cdbr-

c optionally, one could limit wave erosion to the part
c of the berg above (and including) the surface layer
c-option      dwiceb_fv_ave=deltat*flvag*
c-option     &    ( (dz(ks2)+hiceb(l)*remsim) / (hiceb(l)*(1+remsim) )
c one would have to separate lateral melt from wave erosion
c in the scheme below...

c----------------------------------------------------------------------
c--- 3.3) Calculate Lateral melt and total melted volume per layer.
c----------------------------------------------------------------------
        dwiceb=0.

c     do not leave the skip, it freezes grounded bergs...
cdbr check grounding bergs as cause
cdbr comment: doesnt help RATIO, but could still be flawed
cmab        if (kiceb(l).lt.1) then
cmab          write(99,*)'cdbr WARNING icebtraj: kiceb<1 entering melt',
cmab     &      'at:',n,'in (',i,j,') hiceb=',hiceb(l)
c7        ,'SKIPPING MELT'
c7          goto 710
cmab        endif
cdbr-

c----------------------------------------------------------------------
c- loop over layers (bottom up)
        do 700 k=ksubm,ks2
c----------------------------------------------------------------------

          dziceb=dz(k)
c- iceberg height per layer
       if (k.eq.kiceb(l)) dziceb=hiceb(l)+zw(k+1) ! dmr&mab: works only if zw is negative!! checked!
c- tip of the iceberg in bottom layer

cdbr
       if(dziceb.lt.0) then  !should set dziceb to zero??? --- dmr&mab
        write(99,*)'ALERT! dziceb(',l,'):',dziceb,'(k:',k,') at',n,
     & 'in',i,j,'kiceb:',kiceb(l)
!mab: commented otherwise too much mouchard output
!       write(99,*)'dz(k):',dz(k),'hiceb:',hiceb(l),'zw(k+1)',zw(k+1)
!mab: commented part ends here
       endif
cdbr

!### : dmr&mab moved out of the loop !!!
!###c- we2top wave erosion to top layer
!###c- separate flvag_dt (local var in this k-loop; dont touch flvag),
!###c- but keep fl_sum->fl_dt for compensation purposes
!###          flvag_dt=flvag*deltat
!###c- we2top end

c--- 3.2c) Lateral melt ("Buoyant convection")
c----------------------------------------------------------------------
c Lateral melt with a first and second order term according to El-Tahan 87:
c Rlatt=7.62e-3 Tw +1.29e-3 Tw^2 (m/day!!)
c-
!dmr [NOEQUI]           tocean=t(i,j,k)-tkc
          tocean=scal(i,j,k,1)-tkc

!### dmr&mab --- below shall NEVER happen anyway t and scal(...,1) are IDENTICAL in PHYSICAL MEMORY!
!### cdbR          tocean=scal(i,j,k,1)-tkc
!###          if(t(i,j,k)-scal(i,j,k,1).ne.0.) then
!###!mab: commented otherwise too much mouchard output
!###!            write(99,*)'RED ALERT icebtraj; ',
!###!     &      't(',i,j,k,'):',t(i,j,k),'.ne. scal(ijk,1):',scal(i,j,k,1)
!###!mab: commented part ends here
!###          endif

          if(tocean.gt.0d0) then
            flat=flcf1*tocean+flcf2*tocean*tocean
          else
            flat=0.
          endif
cmab
cmab
c- no freeze-on allowed:
cmab	  if (flat.lt.0.) then
cmab            print*,'flat .lt. 0!',tocean
cmab!mab: commented otherwise too much mouchard output
cmab!            write(99,*)'WARNING icebtraj:',
cmab!     &      'flat negative; set to zero; Tocean:',tocean,'k:',k
cmab!mab: commented part ends here
cmab!            flat=0.
cmab!          endif

cdbR          flat=max(flat,zero)

cdbR !!!
c          flat=0.
cdbR-

c--- 3.2d) Total lateral melt per layer fl_dt into dVol_icb(ijk)
c----------------------------------------------------------------------
c--- total lateral melt is fl_dt {m/day} (<wiceb):
          fl_sum=flat+flvag

cdbr
cdbr comment:
!mab: commented otherwise too much mouchard output
!        if ((fl_sum).lt.0.) then
!          write(99,*) 'cdbr WARNING icebtraj.f; fl_sum negative',
!     &    'berg(',l,')s h,w:',hiceb(l),wiceb(l),
!     &    ' at (i,j)=(',i,',',j,'iteration:',n,
!     &    ' Daily wave melt= flvag_dt',flvag_dt,
!     &    'times total heigt:',hiceb(l)*(1+remsim),'=',
!     &     flvag_dt*hiceb(l)*(1+remsim),
!     &    ' Daily lat melt= flat*deltat',flat*deltat,
!     &    ' icb-velocity(u,v):',uiceb(l),viceb(l),' SST:',t(i,j,ksubm)
!        endif
!mab: commented part ends here
cdbr-

          fl_dt=fl_sum*deltat

          if(fl_dt.gt.wiceb(l)) then
!mab: commented otherwise too much mouchard output
!	    write(99,*)'cdbr WARNING icebtraj: Corrected for ',
!     &      'lateral melt exceeds wiceb(',l,')=',wiceb,'in layer:',k,
!     &      'flat=',flat,'dT=',tocean,'flvag=',flvag,'xs=',xs
!mab: commented part ends here
cdbr comment: never mentioned
            fl_dt=wiceb(l)
!###c- flvag_dt is local var of k-loop; shouldnt exceed wiceb
!###            if (flvag_dt.gt.wiceb(l) ) then
!###              flvag_dt=wiceb(l)/deltat
!###            endif
c-
cdbr:dangerous 2do in k-loop	    flvag=wiceb(l)/deltat-flat
c changing flvag, only used in dbug write, below. Beware if used later!
          endif

c- weighted sum over layers for averaging lateral melt {m/day}
c- to determine new width after k-loop:

cdmr&mab --- [fl_dt] = m ; [dziceb] = m
cdmr&mab --- WTF next line?
          dwiceb=dwiceb+dziceb*fl_dt ! ?????????????????

cdmr&mab --- ALTERNATIVE:
cdmr&mab ---  dwiceb=fl_dt


cdbRR          wiceb(l)=wiceb(l)-fl_dt*(dziceb/hiceb(l))

c--- chop basal melt into bottom layers ----------------------
          fb_dt_k=fb_dt
cdbR comment:
c7            write(99,*)
c7            write(99,*) 'c7 icebtraj.f; fb_dt chopped',
c7     &    ' berg(',l,')s h,w:',kiceb(l),hiceb(l),wiceb(l),
c7     &    ' at (i,j)=(',i,',',j,'iteration:',n,
c7     &    ' k,ksubm,kiceb:',k,ksubm,kiceb(l),
c7     &    ' Daily basal melt b4 chop= fb_dt=',fb_dt,
c7     &    'dziceb=',dziceb,
c7     &    'fb_dt-dziceb=',fb_dt-dziceb,'fb_dt_k=',fb_dt_k
cdbR-


c- fb_dt_k is daily basal melt for layer k
          if (fb_dt.gt.dziceb.and.k.lt.ks2) then
c- limit it if melting into next layer (surface layer excluded)

cdbr comment: dziceb = NEGATIVE!? kiceb klopte niet
c7            write(99,*) 'cdbr WARNING icebtraj.f; fb_dt chopped',
c7     &    ' berg(',l,')s h,w:',hiceb(l),wiceb(l),
c7     &    ' at (i,j)=(',i,',',j,'iteration:',n,
c7     &    ' Daily basal melt b4 chop= fb_dt=',fb_dt,
c7     &    'this layers part=fb_dt_k=dziceb=',dziceb,
c7     &    'fb_dt-dziceb=',fb_dt-dziceb
cdbr-
            fb_dt_k=dziceb
	  endif
c- calculate leftovers for next layer:

      fb_dt=fb_dt-fb_dt_k

cdbr comment:
c7          if (fb_dt .gt. 0.) then
c7           write(99,*) 'cdbr icebtraj.f; fb_dt leftover>0',
c7     &         fb_dt
c7          endif
           if (fb_dt .lt. 0.) then
             write(99,*) 'cdbr WARNING icebtraj.f; fb_dt leftover<0!?',
     &         fb_dt
            fb_dt = 0d0
           endif
cdbr


c--- Write relative contributions
c7 outdated, check before use
c7          if( (fb_dt_k/(fb_dt_k+(flat+flvag)*deltat).gt.0.96).or.
c7     &    (flat*deltat/(fb_dt_k+(flat+flvag)*deltat).gt.0.96).or.
c7     &    (flvag*deltat/(fb_dt_k+(flat+flvag)*deltat).gt.0.96) ) then
c7            if(fb_dt_k.gt.10..or.flat*deltat.gt.10..or.
c7   &         flvag*deltat.gt.10.) then
c7            write(99,*)'WARNING icebtraj: ',
c7     &      'relative melt>0.96 berg(',l,') layer(',k,'of',ksubm,'):'
c7            write(99,*)'basal=',fb_dt_k/(fb_dt_k+(flat+flvag)*deltat),
c7     &      'lateral=',flat*deltat/(fb_dt_k+(flat+flvag)*deltat),
c7     &      'waves=',flvag*deltat/(fb_dt_k+(flat+flvag)*deltat),
c7     &       fb_dt_k,flat*deltat,flvag*deltat,'Twater:',tocean,'k',k
c7            endif
c7	  endif

c----------------------------------------------------------------------
c----------------------------------------------------------------------
c--- dVol_icb(i,j,k) is sum of basal and lateral melt minus overlap
c----------------------------------------------------------------------
c-
!      print*,'vor dmrmabversion 1'
#if ( dmrmabversion == 0 )
          dVol_icb(i,j,k)=dVol_icb(i,j,k)+pond_icb(l)*(
c- daily basal melt is bottom area times thickness fb_dt_k:
     &                fb_dt_k*1.5*wiceb(l)*wiceb(l)
c- side melt area is wiceb*1.5*fldt + 1.5*wiceb*fl_dt minus overlap,
c- height is dziceb times extrapolation correction (1+remsim):
c7 get rid  of xpcorr; put remsim part in toplayer after k-loop
     &               +dziceb*fl_dt*(3*wiceb(l)-1.5*fl_dt)
c- minus overlap with bottom melt:
     &               -fb_dt_k*fl_dt*(3*wiceb(l)-1.5*fl_dt) )

          volmeltday(l,k) = volmeltday(l,k)+pond_icb(l)*(
     &                fb_dt_k*1.5*wiceb(l)*wiceb(l)
     &               +dziceb*fl_dt*(3*wiceb(l)-1.5*fl_dt)
     &               -fb_dt_k*fl_dt*(3*wiceb(l)-1.5*fl_dt) )

#else

!      print*,'IN dmrmabversion 1'
cdmr&mab --- basal melting first
          dVol_icb(i,j,k) = dVol_icb(i,j,k) + pond_icb(l) * (
c- daily basal melt is bottom area times thickness fb_dt_k:
     &  fb_dt_k*1.5*wiceb(l)*wiceb(l)
cdmr&mab --- lateral melt on side 1.5*wiceb
     & + 1.5*wiceb(l)*MAX(dziceb-fb_dt_k,0.0)*fl_dt
cdmr&mab --- lateral melt on side (wiceb-fl_dt)
     &+ MAX((wiceb(l)-fl_dt),0.0)*MAX(dziceb-fb_dt_k,0.0)*fl_dt
     &                                                    )
          volmeltday(l,k) = volmeltday(l,k) + pond_icb(l) * (
     &  fb_dt_k*1.5*wiceb(l)*wiceb(l)
     & + 1.5*wiceb(l)*MAX(dziceb-fb_dt_k,0.0)*fl_dt
     &+ MAX((wiceb(l)-fl_dt),0.0)*MAX(dziceb-fb_dt_k,0.0)*fl_dt
     &                                                    )


cdmr&mab ---
cdmr&mab --- So far, we are left with an iceberg that has a shape
cdmr&mab ---  at depth k of:
cdmr&mab --- shape(k) = (wiceb(l)-fl_dt)*(1.5*wiceb(l)-fl_dt)*MAX(dziceb-fb_dt_k,0.0)
cdmr&mab --- shape(k) = ( 1.5*wiceb**2 - 2.5*wiceb*fl_dt + fl_dt**2 )*(dziceb-fb_dt)
cdmr&mab ---

#endif
cdbr commenting redirected we didnt help so its not below here
c- we2top
c- redirect wave erosion to top layer (ks2) from lower layers
#if ( dmrmabversion == 0 )
          if(k.ne.ks2) then
            dVol_icb(i,j,ks2) = dVol_icb(i,j,ks2) + pond_icb(l)*(
c- side melt area is wiceb*1.5*fldt + 1.5*wiceb*fl_dt minus overlap,
c- height is dziceb
     &        dziceb*flvag_dt*(3*wiceb(l)-1.5*flvag_dt) )
            dVol_icb(i,j,k)   = dVol_icb(i,j,k) - pond_icb(l)*(
     &        dziceb*flvag_dt*(3*wiceb(l)-1.5*flvag_dt) )

            volmeltday(l,ks2) = volmeltday(l,ks2)+pond_icb(l)*(
     &        dziceb*flvag_dt*(3*wiceb(l)-1.5*flvag_dt) )

            volmeltday(l,k) = volmeltday(l,k)-pond_icb(l)*(
     &        dziceb*flvag_dt*(3*wiceb(l)-1.5*flvag_dt) )
          endif
c- we2top end
#else


cdmr&mab --- Wave erosion on two sides ...
          if(k.ne.ks2) then
            dVol_icb(i,j,ks2) = dVol_icb(i,j,ks2) + pond_icb(l)*(
cdmr&mab --- lateral melt on side (1.5*wiceb-fl_dt)
     >      +(1.5*wiceb(l)-fl_dt)*MAX(dziceb-fb_dt_k,0.0)*flvag_dt
!     >      + MAX(0.0d0,(1.5*wiceb(l)-fl_dt))*MAX(dziceb-fb_dt_k,0.0)*flvag_dt
cdmr&mab --- lateral melt on side (wiceb-fl_dt)-flvag_dt
     >      + MAX(0.0d0,(wiceb(l)-fl_dt)-flvag_dt)*MAX(dziceb-fb_dt_k,0.0)
     > *flvag_dt
     >                                                          )
            volmeltday(l,ks2) = volmeltday(l,ks2) + pond_icb(l)*(
     >      +(1.5*wiceb(l)-fl_dt)*MAX(dziceb-fb_dt_k,0.0)*flvag_dt
     >      + MAX(0.0d0,(wiceb(l)-fl_dt)-flvag_dt)*MAX(dziceb-fb_dt_k,0.0)
     > *flvag_dt
     >                                                          )
cdmr&mab --- suppress the water added to the top layer from the lower levels
            dVol_icb(i,j,k) = dVol_icb(i,j,k) - pond_icb(l)*(
cdmr&mab --- lateral melt on side (1.5*wiceb-fl_dt)
     >      +  (1.5*wiceb(l)-fl_dt)*MAX(dziceb-fb_dt_k,0.0)
     >         *flvag_dt
!     >      +  MAX(0.0d0,(1.5*wiceb(l)-fl_dt))*MAX(dziceb-fb_dt_k,0.0)
!     >         *flvag_dt
cdmr&mab --- lateral melt on side (wiceb-fl_dt)-flvag_dt
     >      + MAX(0.0d0,(wiceb(l)-fl_dt)-flvag_dt)*MAX(dziceb-fb_dt_k,0.0)
     >        *flvag_dt
     >                                                          )
                volmeltday(l,k) =volmeltday(l,k) - pond_icb(l)*(
     >      +  (1.5*wiceb(l)-fl_dt)*MAX(dziceb-fb_dt_k,0.0)
     >         *flvag_dt
     >      + MAX(0.0d0,(wiceb(l)-fl_dt)-flvag_dt)*MAX(dziceb-fb_dt_k,0.0)
     >        *flvag_dt
     >                                                          )

          endif
cdmr&mab ---
cdmr&mab --- So far, we are left with an iceberg that has a shape
cdmr&mab ---  at depth k of:
cdmr&mab --- shape(k) = MAX(0.0,(wiceb-fl_dt-flvag_dt))*MAX(0.0,(1.5*wiceb-fl_dt-flvag_dt))*MAX(dziceb-fb_dt_k,0.0)
cdmr&mab ---
#endif


c7
c7          write(99,'(A,I6,A,I5,2I4,A,F32.12)')
c7     &      'it:',n,', dVol_icb(',i,j,k,')=',dVol_icb(i,j,k)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
 700    continue
c--- end k-loop over layers

cdbr
 710    continue
cdbr

#if ( dmrmabversion == 0 )

c--- average lateral melt is dwiceb = weighted SUM/hiceb {m/day}
        dwiceb=dwiceb/hiceb(l)

c----------------------------------------------------------------------
c7 put remsim part in surface layer instead of using xpcorr
c bars have width dwiceb and height remsim*hiceb
        tip_of_the_iceberg=pond_icb(l)*(
     &        remsim*hiceb(l)*dwiceb*(3*wiceb(l)-1.5*dwiceb)
c7 add leftover basal melt (non-zero if it exceeded submarine hiceb)
     &        +fb_dt*wiceb(l)*1.5*wiceb(l)
c subtract overlap
     &        -fb_dt*dwiceb*(3*wiceb(l)-1.5*dwiceb) )
        dVol_icb(i,j,ks2) = dVol_icb(i,j,ks2) + tip_of_the_iceberg
        volmeltday(l,ks2) = volmeltday(l,ks2) + tip_of_the_iceberg

c7        write(99,*)'tip_of_the_iceberg:',tip_of_the_iceberg
c7        write(99,*)'divide by dVol_icb_ks2=',
c7   &      tip_of_the_iceberg/dVol_icb(i,j,ks2)

#else
cdmr&mab ---
cdmr&mab --- Deal with the "tip of the icebergs"
cdmr&mab ---  The emerged part of the iceberg is still in original size
cdmr&mab ---   (it has not been eroded or melted in the previous loop)
cdmr&mab ---

       if (fb_dt.GT.0.0d0) then
          tip_of_the_iceberg=1.5*wiceb(l)**2*fb_dt*pond_icb(l)
          dVol_icb(i,j,ks2) = dVol_icb(i,j,ks2) + tip_of_the_iceberg
          volmeltday(l,ks2) = volmeltday(l,ks2) + tip_of_the_iceberg
       endif

cdmr&mab ---
cdmr&mab --- Add all melting to the dVol
cmab-not to be done (otherwise some melt counted twice!!!
cmab dVol_icb(i,j,:)=dVol_icb(i,j,:) + volmeltday(l,:)
cmab-not to be done (otherwise some melt counted twice!!!
cdmr&mab ---
#endif
c----------------------------------------------------------------------


#if ( dmrmabversion == 0 )

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- 3.4) NEW SIZE; eliminate small bergs; volume; -difference=melt fonte_icb
c----------------------------------------------------------------------

        wiceb(l)=wiceb(l)-dwiceb
c NOTE: berg sinks a little bit to update the ratio hnew to remsim
c- dbugging: changed calculation of hnew (see above),
c- to prevent melt exceeding iceberg height
c CHECKed        hiceb(l)=(hnew+remsim*hiceb(l))/(1+remsim)
        hiceb(l)=hnew/(1+remsim)

cdbr
        if (hiceb(l).lt.0. .or. wiceb(l).lt.0) then
cmab          write(99,*)'icebtraj: negative new dimensions set to zero',
cmab &      'hiceb(',l,')=',hiceb(l),'wiceb:',wiceb(l)
c        endif

          hiceb(l)=max(hiceb(l),zero)
          wiceb(l)=max(wiceb(l),zero)
cmab          write(99,*)'checking if max(..,zero) worked:',
cmab &      'hiceb(',l,')=',hiceb(l),'wiceb:',wiceb(l)
        endif
!        if(mod(nit,10).eq.0) then
!        write(250,'(1I8,4F40.15)')l, wiceb(l) ,hiceb(l),hnew/(1+remsim)
!     &    ,fb_dt
!        endif
cdbr
#else
        if (kiceb(l) .gt. -1) then
cdmr&mab --- New total height of the icebergs is hnew
!        if(mod(nit,10).eq.0) then
!        write(250,'(1I8,5F40.15)')l, wiceb(l) ,hiceb(l),hnew
!     &    ,fb_dt,(1+remsim)
!        endif

           if ( hnew.lt.0d0) print*,'ALERT!!!hnew',hnew

           V_new(l) = MAX(0.0d0,(vol_in(l) - SUM(volmeltday(l,:))))

          if (V_new(l) .gt. 0d0 .and.hnew .gt.0d0 .and. pond_icb(l).gt.0d0)
     >    then
            wiceb(l) = SQRT(V_new(l) / (hnew * 1.5*pond_icb(l)))
            hiceb(l)=hnew/(1+remsim)

            if (hiceb(l).lt.0. .or. wiceb(l).lt.0) then
cmab          write(99,*)'icebtraj: negative new dimensions set to zero',
cmab &      'hiceb(',l,')=',hiceb(l),'wiceb:',wiceb(l)
c        endif

              hiceb(l)=max(hiceb(l),zero)
              wiceb(l)=max(wiceb(l),zero)
              kiceb(l)=ks2+1
              goto 790
cmab          write(99,*)'checking if max(..,zero) worked:',
cmab &      'hiceb(',l,')=',hiceb(l),'wiceb:',wiceb(l)
            endif

          else
            hiceb(l) = 0d0
            wiceb(l)=0d0
            kiceb(l)=ks2+1
            goto 790
          endif
!mab          if(mod(nit,360).eq.0) then
!mab            write(250,'(1I8,3F40.15)')l,V_new(l),vol_in(l),SUM(volmeltday(l,:))
!mab             write(250,'(1I8,3F40.15)')l, wiceb(l) ,hiceb(l),hnew/(1+remsim)
!mab           endif
        else
          goto 800
        endif
cdmr&mab ---
#endif
c----------------------------------------------------------------------
c bergs smaller than 5 are reduced to zero (instant melt)...
c BEFORE the volume difference is calculated, so we have to
c add this volume loss to dVol_icb(i,j,surfacelayer) to remain consistent
c these bergs are ELIMINATED, along with other melted bergs (warn??),
c in 'mise a jour', below

        if (wiceb(l).le.5.) then
          dVol_icb(i,j,ks2)=dVol_icb(i,j,ks2)+pond_icb(l)*
     &      hiceb(l)*(1+remsim)*wiceb(l)*1.5*wiceb(l)
          volmeltday(l,ks2) = volmeltday(l,ks2)+ pond_icb(l)*
     &      hiceb(l)*(1+remsim)*wiceb(l)*1.5*wiceb(l)
         wiceb(l)=0.
	     hiceb(l)=0.
	     kiceb(l)=ks2+1
        else
          if (hiceb(l).le.5.) then
            dVol_icb(i,j,ks2)=dVol_icb(i,j,ks2)+pond_icb(l)*
     &        hiceb(l)*(1+remsim)*wiceb(l)*1.5*wiceb(l)
            volmeltday(l,ks2) = volmeltday(l,ks2)+ pond_icb(l)*
     &        hiceb(l)*(1+remsim)*wiceb(l)*1.5*wiceb(l)

            wiceb(l)=0.
	        hiceb(l)=0.
            kiceb(l)=ks2+1
          endif
        endif
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c 	NEW VOLUME and "filling the icb-array"??
c- re-calcul du volume de l'iceberg "l" et remplissage des tableaux
c----------------------------------------------------------------------
c due to the (1+remsim) factor,
c the new underwater berg geometry is extrapolated to the above water part!!
c therefore we had to give the lateral melt per layer an
c (1+remsim) bonus to obtain the same total volume change...

        V_new(l)=wiceb(l)*wiceb(l)*1.5*hiceb(l)*(remsim+1)*pond_icb(l)
        if (V_new(l).lt.0.) V_new(l)=0.

c----------------------------------------------------------------------
c putting the volume difference into fonte_icb(i,j), which is
c only used when melt-fluxes are not separated per layer,
c (USED TO BE fwf in ec_co2oc.f) and lhf into phiss1 in thersf.f)
C now (2011) it's per layer in phiss2 in thersf

        fonte_icb(i,j)=fonte_icb(i,j)+((V_old(l)-V_new(l))*
     &    pond_icb(l))

c---------------------------------------------------------------------
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- skip to here when avoiding melt:

 790  continue

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c 4)	Updating iceberg layer-depth kiceb(l)
c-	mise a jour de "kiceb" :
c----------------------------------------------------------------------
c kiceb(l) is the number of the iceberg bottoms ocean layer,
c  	counting down from the surface (ks2=20).
c 	[:78: ksubm=max(kiceb(l),kfs(i,j)) (kfs=bottomdepth)]
c kiceb=0 means: iceberg grounded in the water, see line 160,updated here(~:561:)
c kiceb=-1 means: iceberg was ELIMINATEd. Volume not counted if
c	*) terminated on land or off the grid and no valid old gridcell
c kiceb=-2 means free array spot
c kiceb=ks2+1 means fully melted; ready to be ELIMINATEd

c- calculate volume of grounded bergs after melt to determine grounded melt
        if (kiceb(l).eq.0) then
          vol_grounded_after=vol_grounded_after+
     &      pond_icb(l)*1.5*wiceb(l)*wiceb(l)*hiceb(l)*(1+remsim)
          if (yn(l).lt.0) then
	     vol_grounded_after_south=vol_grounded_after_south +
     &      pond_icb(l)*1.5*wiceb(l)*wiceb(l)*hiceb(l)*(1+remsim)
          endif
	     endif
c-

cdbr cut kiceb update from here and moved it to the end of the loop
cdbr (after rollover and elimination of bergs that melted thru surface)

c--- eliminate bergs that melted thru surface layer
        if (kiceb(l).ge.(ks2+1).or.hiceb(l).le.0.) then
c- ELIMINATEd:
          if(hiceb(l).ne.0d0) print*,'eliminated!1 ', kiceb(l),hiceb(l)
          kiceb(l)=-1
          uiceb(l) = 0.
          viceb(l) = 0.
          nbricb_moins=nbricb_moins+(1*pond_icb(l))
          if(j.le.int(jmax/2))nbricb_moins_s=nbricb_moins_s
     &                     +(1*pond_icb(l))
          if(j.gt.int(jmax/2))nbricb_moins_n=nbricb_moins_n
     &                     +(1*pond_icb(l))
	  goto 800
        endif

c--- eliminate bergs that melted to < 5m. width
        if (wiceb(l).le.0.) then
!          print*,'eliminated!2 '
c- ELIMINATEd:
          kiceb(l)=-1
          uiceb(l) = 0.
          viceb(l) = 0.
          nbricb_moins=nbricb_moins+(1*pond_icb(l))
          if(j.le.int(jmax/2))nbricb_moins_s=nbricb_moins_s
     &                +(1*pond_icb(l))
          if(j.gt.int(jmax/2))nbricb_moins_n=nbricb_moins_n
     &           +(1*pond_icb(l))
          goto 800
        endif
c---

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c 5) 	ROLL_OVER
c-- 	Icebergs qui tournent sur eux-memes:
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        if(kiceb(l).gt.-1) then
          wwiceb=wiceb(l)*1.5
          htot=(1+remsim)*hiceb(l)
          ratio = wwiceb / htot
          cfstab=(sqrt(cfst1+(cfst2/htot)))

C          if (mod(n,nfricb).eq.0) then
c- this works... writing once per month if nfricb = 30
C            if (mod(l,100).eq.3)then
C              if((ratio/cfstab).lt.(1.4)) then
C                write (99,*)'icebtraj ROLL-OVER',
C     &      ' test mod(n,nfricb).eq.0..,l,ww,htot,ww/htot,cfstab:(',
C     &       n,nfricb,')',l,wwiceb,htot,wwiceb/htot,cfstab
C              endif
C            endif
C          endif
c---
          if (ratio.lt.cfstab) then
!          print*,'Berg rolls!'
	    i_rolled=i_rolled+1
            old_volume=(1+remsim)*hiceb(l)*wiceb(l)*1.5*wiceb(l)
c-
cwa            write(99,*)'icebtraj: INSTABLE iceberg(',l,') ',
cwa     &        '(ww/htot < cfstab )=(',
cwa     &         ratio,'<',cfstab,') ROLLING at ',
cwa     &        'it_icb=',n,'H,W=', hiceb(l), wiceb(l),
cwa     &        'old_volume=',old_volume
c-
c-- the geometry of this trafo is wrong: small bug!
            hiceb(l)=wiceb(l)/(1+remsim)
            wiceb(l)=sqrt(htot*wwiceb/1.5)
            wwiceb=wiceb(l)*1.5
            htot=(1+remsim)*hiceb(l)
c--
            new_volume=(1+remsim)*hiceb(l)*wiceb(l)*1.5*wiceb(l)
c-
cwa            write(99,*)'icebtraj ROLLED: ',
cwa     &        'New H,W=', hiceb(l),wiceb(l),
cwa     &        'new_volume=',new_volume
c-
            if(new_volume/old_volume.ne.1.) then
               vol_lost_roll=vol_lost_roll+
     &            old_volume-new_volume

!              print*,'vol_lost:',vol_lost_roll
c-
cwa              write(99,*)'ALERT icebtraj: Volume changed by roll-over',
cwa     &          'new/old volume=',new_volume/old_volume,
cwa     &          'new-old volume=',new_volume-old_volume
c-
            endif
          endif
c---
        endif

c------ recalculate depth kiceb of all live -incl. grounded- icebergs
        if(kiceb(l).ge.0) then
c-
          kiceb(l)=ks2
          do 795 k=ks2,2,-1
            if (-zw(kiceb(l)).le.hiceb(l)) then
              kiceb(l)=k-1
              if(kiceb(l) .eq. -1) print*,'kiceb -1!!'
            else
              goto 795
            endif
 795      continue
c-
        endif
c------


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- end of the main loop over all bergs 'l'

!mab        if(mod(nit,360).eq.0) then
!mab        write(250,'(2I8,3F40.15)') nit,l
!mab      & ,sum(volmeltday(l,:)),sum(dVol_icb(:,:,:)),vol_orig(l)
!mab        endif
 800  continue
!mab       close(250)


#if (dmrmabversion == 0)
         vol_rest(:,:) = 0d0
      do l=1,lmx
        if(kiceb(l).gt.-1) then
        vol_diff(l) = vol_in(l)-(V_new(l)+sum(volmeltday(l,:)))

        vol_melt(l) = vol_melt(l)+ sum(volmeltday(l,:)) + vol_diff(l)

        if(abs(vol_diff(l)) .gt. 200d0) then
          open(250,file='outputdata/ism/vol_diff.dat',position='append')
          write(250,'(2I6,4F20.8)') l, kiceb(l),vol_diff(l),vol_in(l),V_new(l)
     &    ,sum(volmeltday(l,:))
!mab          write(250,'(2I6,9F20.8)') l, kiceb(l),vol_diff(l),vol_in(l),V_new(l)
!mab     &   ,sum(volmeltday(l,:)),(V_new(l)+sum(volmeltday(l,:))),vol_melt(l)
!mab     &   ,vol_orig(l),vol_orig(l)-vol_melt(l),sum(dVol_icb(:,:,:))
          close(250)
        endif

          i=1+nint((xn(l)-xwi1)/dxwi)
          j=1+nint((yn(l)-ywj1)/dywj)


          dVol_icb(i,j,ks2)=dVol_icb(i,j,ks2)+vol_diff(l)

          else

           vol_melt(l) = vol_melt(l)+ sum(volmeltday(l,:)) + vol_diff(l)

        endif

       if(kiceb(l).le.-1) then
          i=1+nint((xn(l)-xwi1)/dxwi)
          j=1+nint((yn(l)-ywj1)/dywj)
          vol_rest(i,j)=vol_orig(l)-vol_melt(l)
          vol_orig(l)=0d0
          vol_melt(l)=0d0

!mab       open(25,file='outputdata/icebergs/checkvol.dat',position='append',form='formatted')
          dVol_icb(i,j,ks2)=dVol_icb(i,j,ks2)+vol_rest(i,j)
!mab          if (abs(vol_rest(i,j)) .gt. 10d0) then
!mab            write(25,'(3I6,1F20.8)') i,j,l,dVol_icb(i,j,ks2),vol_rest(i,j)
!mab          endif
!mab          close(25)

        endif
       enddo
#endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c 6) WRITING VOLUMES
c	write warnings and grand total

c--- total melt per layer today ---------------------------------------
c note: local array apparently automatically reset each call
      do kl1=1,ks2
        dvol_k(kl1)=0.
      enddo
      dvol=0.
      sum_jono3=0.
c---
      do kl2=1,ks2
cdb: YES ; are these i,j limits cool? dVol make sense theer?
        do i=1,imax
          do j=1,jmax
            dvol_k(kl2)= dvol_k(kl2)+dVol_icb(i,j,kl2)
c-
            if(dVol_icb(i,j,kl2).lt.0.) then
              write(99,*)'RED ALERT icebtraj:',
     &         'at it',n,'melted icb volume negative (i,j,layer):',
     &         'dvol_icb(',i,',',j,',',kl2,')=',dVol_icb(i,j,kl2)
            endif
c-
cmab            if(dVol_icb(i,j,kl2).gt.10000000000) then
cmab              write(99,*)'ALERT icebtraj:',
cmab     &         'at it',n,'melted icb volume>10000000000 (i,j,layer):',
cmab     &         'dvol_icb(',i,',',j,',',kl2,')=',dVol_icb(i,j,kl2)
cmab            endif
c-
          enddo
        enddo
        dvol=dvol+dvol_k(kl2)

cdb        write(99,*)dvol,kl2,'dvol cumulative'
c-
cdb        if(dvol_k(kl2).lt.0.) then
cdbR        write(99,*)'RED ALERT icebtraj:',
cdbr     &      'at it',n,'layer-total melted icb volume negative',
cdbr     &      'dvol_k(',kl2,')=',dvol_k(kl2)
cdb        endif
c-
      enddo
c7      write(99,*)'todays dvol=',dvol



c--- total melt per year ------------------------------
      yr_vmg=yr_vmg+(vol_grounded-vol_grounded_after)
      yr_vmgs=yr_vmgs+(vol_grounded_south-vol_grounded_after_south)
      yr_va=yr_va+vol_all
      yr_ioob=yr_ioob+i_OOB
      yr_ir=yr_ir+i_rolled
      yr_vlOOB=yr_vlOOB+vol_lost_OOB
      yr_vlOOBs=yr_vlOOBs+vol_lost_OOB_south
      yr_vlr=yr_vlr+vol_lost_roll
      yr_vmt=yr_vmt+dvol

      do i=1,imax
        do j=1,jmax
          fonte_icbtr_dp=fonte_icbtr_dp+fonte_icb(i,j)
        enddo
      enddo
c--- cdbR --- monthly output ------------------------
      if(mod(n,30).eq.0)then
!mab: commtend, otherwise too much output in mouchard
!        write(99,*)
!        do kl9=1,20
!          write(99,'(A,I6,A,I4,A,F32.12)')
!     &      'it:',n,', dvol_k(',kl9,')=',dvol_k(kl9)
!        enddo
cdbR debug RATIO dVol/fonte
        fonte_icbtr=0.
        do i=1,imax
          do j=1,jmax
            fonte_icbtr=fonte_icbtr+fonte_icb(i,j)
          enddo
        enddo
!        write(99,*)'ratio todays fonte/todays dvol; ',
!     &    'fonte_icbtr/dvol',fonte_icbtr/dvol
!        write(99,*)'todays cum FONTE in icebtraj:',fonte_icbtr
!        write(99,*)'total cum FONTE in icebtraj:',fonte_icbtr_dp
! mab: end of commented part
c7 sum_jono4 gives 'same' result as
c sum-yr_vmtk_lyr{cum dvol_k} and cum-time-sum[todays_dvol(sumKI)] yr_vmt
c they can differ at 5th decimal... apparently a rounding error
c7      do kl4=10,ks2
c7        sum_jono4=sum_jono4+dvol_k(kl4)
c7      enddo
c7      write(99,*)'yr_cum (not reset)jono4:Ksum{dvol_k}=',sum_jono4
cdb: apparantly yr_vmtk(ks2) array got reset automatically at each call..
c      write(99,*)'(yr_vmtk_20 b4 adding dvol_k(20):)',yr_vmtk_20
c7      yr_vmtk_20=yr_vmtk_20+dvol_k(20)
c7      yr_vmtk_19=yr_vmtk_19+dvol_k(19)
c7      yr_vmtk_18=yr_vmtk_18+dvol_k(18)
c7      yr_vmtk_17=yr_vmtk_17+dvol_k(17)
c7      yr_vmtk_16=yr_vmtk_16+dvol_k(16)
c7      yr_vmtk_15=yr_vmtk_15+dvol_k(15)
c7      yr_vmtk_14=yr_vmtk_14+dvol_k(14)
c7      yr_vmtk_13=yr_vmtk_13+dvol_k(13)
c7      yr_vmtk_12=yr_vmtk_12+dvol_k(12)
c7      yr_vmtk_11=yr_vmtk_11+dvol_k(11)
c7      yr_vmtk_10=yr_vmtk_10+dvol_k(10)
c7      write(99,*)'sum of yr_vmtk_lyr{cum dvol_k} numbers=',
c7     &                 yr_vmtk_20+yr_vmtk_19+yr_vmtk_18+yr_vmtk_17+
c7     &                 yr_vmtk_16+yr_vmtk_15+yr_vmtk_14+yr_vmtk_13+
c7     &                 yr_vmtk_12+yr_vmtk_11+yr_vmtk_10

cdbR
        write(99,*)'cum-time-sum[todays_dvol(sumKI)] yr_vmt=',yr_vmt
        write(99,*)'ratio f/dv _dp/yr_vmt=',fonte_icbtr_dp/yr_vmt
      endif
c---

cdb-
c------------ yearly icebtraj output -------------
      if (mod(n,360).eq.0) then
        write(99,*)
        write(99,*)'icebtraj: YEAR_TOTAL iceberg volumes and melts:',
     &    ' at it_icb',n
c- lost volume ratio
        write(99,*)yr_vlOOB,': volume lost by',yr_ioob,
     &    ' Out Off Bounds icebergs;',
     &    yr_vlOOBs,': of which in the south (eliminated all)'
        write(99,*)yr_vlr,': volume lost by',yr_ir,'rolling bergs'
	write(99,*)yr_va,': total volume all bergs before melt',
     &    ' lmx=',lmx
        write(99,*)(yr_vlOOB+yr_vlr)/yr_va,
     &    ': RATIO (total lost volume)/(total iceberg volume)'
c- todays grounded volume ratio
!mab: commtend, otherwise too much output in mouchard
!        write(99,*)vol_grounded,': grounded volume TODAY',
!     &    ' of', i_grounded,' icebergs'
!        write(99,*)vol_all,': volume of all bergs today'
!	write(99,*)vol_grounded/vol_all,
!     &    ': todays-RATIO (grounded volume)/(total volume)'
!        write(99,*)vol_grounded_south/vol_all,
!     &    ': todays-RATIO (grounded volume SOUTH)/(total volume)'
c- grounded melt ratio
!        write(99,*)yr_vmg,
!     &    ': year-total grounded melt {m3}..SOUTH:',yr_vmgs
!	write(99,*)yr_vmt,': year-total all iceberg melt'
!        write(99,*)yr_vmt+vol_all,'(:yt melt+vol bergs today)',
!     &    'minus tv allb4melt:',yr_va,'=',yr_vmt+vol_all-yr_va
!	write(99,*)yr_vmg/yr_vmt,': RATIO (grounded melt)/(total melt)'
!	write(99,*)yr_vmgs/yr_vmt,
!     &    ': RATIO (grounded melt SOUTH)/(total melt)'
!        write(99,*)
c- total melt per layer
!        write(99,*)yr_vmtk_20, 'year-total vol flux to surface'
!        write(99,*)yr_vmtk_19, 'year-total vol flux to lyr 19'
!        write(99,*)yr_vmtk_18, 'year-total vol flux to lyr 18'
!        write(99,*)yr_vmtk_17, 'year-total vol flux to lyr 17'
!        write(99,*)yr_vmtk_16, 'year-total vol flux to lyr 16'
!        write(99,*)yr_vmtk_15, 'year-total vol flux to lyr 15'
!        write(99,*)yr_vmtk_14, 'year-total vol flux to lyr 14'
!        write(99,*)yr_vmtk_13, 'year-total vol flux to lyr 13'
!        write(99,*)yr_vmtk_12, 'year-total vol flux to lyr 12'
!        write(99,*)yr_vmtk_11, 'year-total vol flux to lyr 11'
!        write(99,*)yr_vmtk_10, 'year-total vol flux to lyr 10'
!        write(99,*)'sum of these yr_vmtk_lyr numbers=',
!     &                 yr_vmtk_20+yr_vmtk_19+yr_vmtk_18+yr_vmtk_17+
!     &                 yr_vmtk_16+yr_vmtk_15+yr_vmtk_14+yr_vmtk_13+
!     &                 yr_vmtk_12+yr_vmtk_11+yr_vmtk_10

!        write(99,*)
!mab: end of commented part

c-
	yr_vlOOB=0.
	yr_ioob=0
	yr_vlOOBs=0.
	yr_vlr=0.
	yr_ir=0
	yr_va=0.
	yr_vmg=0.
	yr_vmgs=0.
	yr_vmt=0.
        yr_vmtk_20=0.
        yr_vmtk_19=0.
        yr_vmtk_18=0.
        yr_vmtk_17=0.
        yr_vmtk_16=0.
        yr_vmtk_15=0.
        yr_vmtk_14=0.
        yr_vmtk_13=0.
        yr_vmtk_12=0.
        yr_vmtk_11=0.
        yr_vmtk_10=0.

      endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- update amounts{eliminated}
      nbricb=nbricb-nbricb_moins+nbricb_plus
      nbricb_n=nbricb_n-nbricb_moins_n+nbricb_plus_n
      nbricb_s=nbricb_s-nbricb_moins_s+nbricb_plus_s
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      end subroutine icebtraj
