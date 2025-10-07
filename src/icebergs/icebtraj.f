#include "choixcomposantes.h"

      subroutine icebtraj(l)
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
c	*) melted volume {m3/day} put in dVol_icb(i,j,k),
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
c c7: march07 getting rid off xpcorr. remsim part of melt piped to surface layer
c thus make sure you dont produce grounded bergs
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!! START_OF_USE_SECTION
      use global_constants_mod, only: dp, ip, tK_zero_C, deg_to_rad            
      use const_mod, only: zero
      use bloc0_mod, only: tms, kfs, u, v, unsdx, unsdy, zw,dz,ks2,scal,
     >    xslon, yslat
      use bloc_mod, only: smx, smy
      use ice_mod, only: imax, jmax, kmax
      use iceberg_mod, only: dticeb,
     >    dVol_icb,
     >    lmx,kiceb,wiceb,hiceb,remsim, pond_icb, uiceb,viceb,
     >    xn, yn, repuls, wind10_u, wind10_v,
     >    icbmax,iceberg_info_out_id,
     >    d_edge,ti,fbcf,flcf1,flcf2,flcf3,cfst1,cfst2,
     >    iceb_circbounc, iceb_findocean,
     >    uiceb_reg,viceb_reg
      use reper_mod, only: xwi1, dxwi, ywj1, dywj
      use to_and_from_CLIO, only: get_indexes_C, get_lonlat_C
      use delta_latlon, only: get_Nlatlon
      use united_Grid_mod, only: grid_elmnt, getUnifiedGrid_indx   
      
!! END_OF_USE_SECTION

      implicit none

c--variables locales :
      real(kind=dp), dimension(kmax) :: un, vn

      real(kind=dp) :: V_new(icbmax)
      real(kind=dp) :: vol_in(icbmax)
      real(kind=dp) :: volmeltday(icbmax,kmax)

c- melting variables
      real(kind=dp) :: temp, vit, hnew, tocean
      real(kind=dp) :: flat, flvag, flvag_dt, fl_sum, fl_dt
      real(kind=dp) :: dziceb, fb, fb_dt, fb_dt_k

c- check on mass conservation after rolling bergs
      real(kind=dp) :: old_volume, new_volume

      real(kind=dp), dimension(2)   :: lon_lat_C_c, lon_lat_C_i, lon_lat_C_j
      integer(kind=ip):: n,l,i,j,ksubm,k,inew,jnew, indexes_C(2)
      real(kind=dp) :: x,xx,yy,xx1,yy1,htot,xs, xnew,
     >                   ynew,xfron,yfron,wwiceb,tip_of_the_iceberg,
     >                   cfstab, V_tot

      real(kind=dp) :: uiceb_temp, viceb_temp
      type(grid_elmnt) :: my_GE
      
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
     
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|      
! Initialize some variables
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      volmeltday(:,:) = 0d0
      V_new(:) = 0d0
      vol_in(:) = 0d0

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-
c- main loop over all icebergs
c-
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      do 800 l=1,lmx
      
c- Negative kiceb means there is no berg to be moved or melted
c- thus TERMINATEd icebergs do not enter this routine again.
        if (kiceb(l).le.-1) goto 800
        vol_in(l) = wiceb(l)*wiceb(l)*1.5*hiceb(l)*(remsim+1)*pond_icb(l)
        if (vol_in(l).le.0d0) then !--- ELIMINATEd:
          kiceb(l)= -1
          uiceb(l) = 0.
          viceb(l) = 0.
          goto 800
        endif
       
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- position iceberg au pas de temps n+1 initial position (and depth) of iceberg "l"

c- Here the (int) grid numbers of the icebergs position (xn,yn) are determined
c- feedomL.f: 1er point scalaire (1,1) lat,long en degre: xwi1=25.5, ywj1=-79.5
c- (xi,yj) are the left bottom corners of the grids (i,j),
c- (presuming xwi1,ywi1 are in the center.)
c- (xx,yy) are distance (0 to 1) to the berg; inverse 'interpolation weights'
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Position of icebergs on clio grid
        indexes_C = get_indexes_C(xn(l), yn(l))
        i=indexes_C(1)
        j=indexes_C(2)

!- Distance interpolation weights
! PB: I had to use min/max in order to keep the weights between 0 and 1; should be checked.
! PB: for now don't use these weights, just use center-point values
!        lon_lat_C_c = get_lonlat_C(i,j)
!        xx=(xn(l)-lon_lat_C_c(1))
!        if (xx.ge.0) then
!          lon_lat_C_i = get_lonlat_C(i+1,j)
!          xx = min(xx / (lon_lat_C_i(1)-lon_lat_C_c(1)),0.5)+(1./2.)
!        else
!          lon_lat_C_i = get_lonlat_C(i-1,j)
!          xx = max(xx / (lon_lat_C_c(1)-lon_lat_C_i(1)),-0.5)+(1./2.)       
!        endif
!        yy=(yn(l)-lon_lat_C_c(2))
!        if (yy.ge.0) then
!          lon_lat_C_j = get_lonlat_C(i,j+1)
!          yy = min(yy / (lon_lat_C_j(2)-lon_lat_C_c(2)),0.5)+(1./2.) 
!        else
!          lon_lat_C_j = get_lonlat_C(i,j-1)
!          yy = max(yy / (lon_lat_C_c(2)-lon_lat_C_j(2)),-0.5)+(1./2.)
!        endif

!        if (xx.lt.0.or.xx.gt.1.or.yy.lt.0.or.yy.gt.1) then
!            write(iceberg_info_out_id,*) 'Error xx or yy',l,xx,yy,
!     &          xn(l),yn(l),lon_lat_C_c, lon_lat_C_i, lon_lat_C_j
!        endif
                        
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Calculate depth of the bottom of the iceberg in k-levels; either iceberg bottom or bottom of the ocean (meaning it is grounded)
        ksubm=max(kiceb(l),kfs(i,j))

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  Movement of iceberg "l"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Ocean speed per layer interpolated at icb position
c- ps.  u(.,.,.) (look like corner values: vectors!)
c- Calcul des vit. oc. a la posit del iceb aux diff niveaux:
c-----------------------------------------------------------------------
! PB Using xx=0 and yy=0, yy1=0 and xx1=0 will result in no weighting of all iceberg dynamics forcings (in icebtraj and icebdyn)
        xx = 0
        yy = 0
        xx1 = 0
        yy1 = 0
        do  k=ksubm,ks2
!          yy1=1.+tms(i,j-1,k)*(tms(i,j+1,k)*yy-1.)
!          xx1=1.+tms(i-1,j,k)*(tms(i+1,j,k)*xx-1.)
          un(k) = (1.-yy1)*( xx*u(i+1,j,k)+(1.-xx)*u(i,j,k) )
     &          +  yy1*( xx*u(i+1,j+1,k)+(1.-xx)*u(i,j+1,k) )
          vn(k) = (1.-yy)*( xx1*v(i+1,j,k)+(1.-xx1)*v(i,j,k) )
     &          +  yy*( xx1*v(i+1,j+1,k)+(1.-xx1)*v(i,j+1,k) )
        end do
            

        
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
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
          kiceb(l)=0
          uiceb(l)=0.
          viceb(l)=0.
        endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- calculate acceleration of (non-stranded) iceberg l by calling icebdyn
c-----------------------------------------------------------------------

        call icebdyn(i,j,l,un,vn,xx,yy)

!PB Change uiceb(l) and viceb(l) from the rotated CLIO grid to a regular grid if on the Northern Hemisphere
!PB use uiceb_temp and viceb_temp to store these variables on the regular grid, use these to calculate the new iceberg location and as iceberg velocity output to netcdf
        my_GE = getUnifiedGrid_indx(xn(l),yn(l),is_degree=.true.)
        
        if (.not. my_GE%is_inWW) then ! my_GE%is_inWW is false if you are on the rotated grid
          uiceb_reg(l) = viceb(l)*-1
          viceb_reg(l) = uiceb(l)
        else
          uiceb_reg(l) = uiceb(l)
          viceb_reg(l) = viceb(l)        
        endif
        
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  New position of iceberg "l":
c-----------------------------------------------------------------------

        call get_Nlatlon(uiceb_reg(l),viceb_reg(l),dticeb,yn(l),xn(l),ynew,xnew)
        
        call iceb_circbounc(xnew,ynew)
        
        indexes_C = get_indexes_C(xnew, ynew)
        inew=indexes_C(1)
        jnew=indexes_C(2)


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c test if bergs moved onto LAND-cells, then put at (fron= is) water edge,
c at intersection of move-direction and 'fron' and
c REPULSEd orthogonnally from fron.
c same for icebergs moving into waters that are too shallow
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!        if ((tms(inew,jnew,ks2).eq.0).or.
!     &      (kiceb(l).lt.kfs(inew,jnew)) ) then
!           
!          if (inew.lt.i) then
!!            xfron=xwi1+dxwi*(i-1)-0.5*dxwi
!            lon_lat_C_c = get_lonlat_C(inew,jnew)
!            lon_lat_C_i = get_lonlat_C(inew-1,jnew)
!            xfron = lon_lat_C_c(1)-(lon_lat_C_c(1)-lon_lat_C_i(1))*0.5
!            ynew=yn(l)+
!     &        ((ynew-yn(l))/(xnew-xn(l)-d_edge))*(xfron-xn(l))
!            xnew=xfron+d_edge
!            indexes_C = get_indexes_C(xnew, ynew)
!            inew=indexes_C(1)
!            jnew=indexes_C(2)
        
!            uiceb(l) = repuls
!          else
!            if (inew.gt.i) then
!!              xfron=xwi1+dxwi*(i-1)+0.5*dxwi
!              lon_lat_C_c = get_lonlat_C(inew,jnew)
!              lon_lat_C_i = get_lonlat_C(inew+1,jnew)
!              xfron = lon_lat_C_c(1)+(lon_lat_C_i(1)-lon_lat_C_c(1))*0.5
!              ynew=yn(l)+
!     &          ((ynew-yn(l))/(xnew-xn(l)+d_edge))*(xfron-xn(l))
!              xnew=xfron-d_edge
!              indexes_C = get_indexes_C(xnew, ynew)
!              inew=indexes_C(1)
!              jnew=indexes_C(2)

!              uiceb(l) = -repuls
!            endif
!          endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----| 
c- xnew and ynew could be changed, so have to check again
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----| 
      
!          if ((tms(inew,jnew,ks2).eq.0).or.
!     &        (kiceb(l).lt.kfs(inew,jnew)) ) then
      
!            if (jnew.lt.j) then
!!              yfron=ywj1+dywj*(j-1)-0.5*dywj
!              lon_lat_C_c = get_lonlat_C(i,j)
!              lon_lat_C_j = get_lonlat_C(i,j-1)
!              yfron = lon_lat_C_c(2)-(lon_lat_C_c(2)-lon_lat_C_j(2))*0.5
!              xnew=xn(l)+
!     &          ((xnew-xn(l))/(ynew-yn(l)-d_edge))*(yfron-yn(l))
!              ynew=yfron+d_edge
!              indexes_C = get_indexes_C(xnew, ynew)
!              inew=indexes_C(1)
!              jnew=indexes_C(2)
        
!              viceb(l) = repuls
!            else
!              if (jnew.gt.j) then
!!                yfron=ywj1+dywj*(j-1)+0.5*dywj
!                lon_lat_C_c = get_lonlat_C(i,j)
!                lon_lat_C_j = get_lonlat_C(i,j+1)
!                yfron = lon_lat_C_c(2)+(lon_lat_C_j(2)-lon_lat_C_c(2))*0.5
!                xnew=xn(l)+
!     &            ((xnew-xn(l))/(ynew-yn(l)+d_edge))*(yfron-yn(l))
!                ynew=yfron-d_edge
!                indexes_C = get_indexes_C(xnew, ynew)
!                inew=indexes_C(1)
!                jnew=indexes_C(2)
  
!                viceb(l) = -repuls
!              endif
!            endif
!          endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- if the icebergs are still on land after repulsion,
c- adding their volume to the surface layer of their old gridcell,
c- or adding their volume to the closest ocean cell using nearest neighbour
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!PB as a fix because missing repulsion, put icebergs in the closest ocean cell if on land
          if (tms(inew,jnew,ks2).eq.0.) then
c- If the new and original cells are land, find closest ocean cell using nearest neighbour
            if ( (tms(i,j,ks2).eq.0.)) then
              write(iceberg_info_out_id,*)  
     &          'Moved to closest ocean cell',l,xnew,ynew             
	      call iceb_findocean(l, inew,jnew, kiceb(l))
!            else
!              write(iceberg_info_out_id,*)  
!     &          'Melted in original location',l,xnew,ynew         
!            endif
c- putting melt volume in old grid-cell or nearest ocean cell
!            dVol_icb(i,j,(ks2))=dVol_icb(i,j,(ks2))+
!     &         (wiceb(l)*wiceb(l)*1.5*hiceb(l)*(remsim+1)*pond_icb(l))
!            volmeltday(l,ks2)=volmeltday(l,ks2)+
!     &         (wiceb(l)*wiceb(l)*1.5*hiceb(l)*(remsim+1)*pond_icb(l))

!            kiceb(l) = -1 ! Melted at the spot
!            uiceb(l) = 0.
!            viceb(l) = 0.
!            goto 800
            
c- stranding/grounding berg stopped at front and repulsed succesfully
!          else
	    xn(l)=xnew
            yn(l)=ynew
          endif

c- warning for chronicly grounded bergs
          if (kiceb(l).lt.kfs(inew,jnew)) then
            kiceb(l)=0
	    uiceb(l)=0.
	    viceb(l)=0.
          endif

c------ else bergs were not stranding on land-cells or grounding in shallows
        else ! ((tms(inew,jnew,ks2).eq.0).or.(kiceb(l).lt.kfs(inew,jnew)) )
          xn(l)=xnew
          yn(l)=ynew
        endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Get new in and jn positions on CLIO grid corresponding to new lat-lon position.
c----------------------------------------------------------------------
        indexes_C = get_indexes_C(xn(l), yn(l))
        i=indexes_C(1)
        j=indexes_C(2)               
               
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Check iceberg location if iceberg is not melted
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|  
        if (i.le.0.or.i.gt.imax.or.j.le.0.or.j.gt.jmax) then
            write(iceberg_info_out_id,*)
     &           'Iceberg still of the grid, ERROR', 
     &           l,i,j,xn(l),yn(l),tms(i,j,ks2),kiceb(l)
        endif
          
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

 200  continue


c from here on chronicly grounded bergs (kiceb(l)=0) are also treated
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c Iceberg MELT and DiSINTEGRATION
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c check for zero hiceb, because of divisions:

        if(hiceb(l).le.0.) then
	  write(iceberg_info_out_id,*)'ALERT icebtraj; hiceb<=0'
	  goto 790
	endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- BASAL MELT
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c basal melt based on ice-shelf parameterisation (3.17: Weeks and Campbell 1973)
c Rb = 0.58 |Vw-Vi|^0.8 (Tw-Ti) / L^0.2 (m/day)
cmab: why **0.4 if above the equation is with 0.8???

        temp=(scal(i,j,ksubm,1)-ti)
        vit=((un(ksubm)-uiceb(l))*(un(ksubm)-uiceb(l))+
     &       (vn(ksubm)-viceb(l))*(vn(ksubm)-viceb(l)))**0.4
        wwiceb=wiceb(l)*1.5

cdmr&mab --- I do not understand why the following equation pertains to wwiceb only (dmr)
        if(wwiceb .gt. 0d0) then
          fb=fbcf*vit*temp/(wwiceb**0.2)
        else
          fb=0d0
        endif

c-- no negative melt:
        if (fb.lt.0.) then
          fb=0.
        endif

        fb_dt=fb*dticeb

c correct for melt exceeding iceberg-size
        if ((fb_dt).gt.((1+remsim)*hiceb(l))) then
          fb_dt=(1+remsim)*hiceb(l)
        endif

c- hnew is used to find the new, buoyancy adjusted, height. see NOTE below
         hnew=(1+remsim)*hiceb(l)-fb_dt

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--- Lateral wave erosion in surface layer is a function of sea-state
c-   It is extrapolated over total iceberg height (1+remsim)hiceb
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c Rvague = 0.5 Ss    (m/day!!) , with
c Ss = -5 + (32+2|Va|)^0.5 (=1 if Va=0)
c ("Both above and below waterline")
c Bigg et al '97, quoting El-Tahan et al '87

        xs = -5 + sqrt( 32+2*
     &  sqrt(wind10_u(i,j)*wind10_u(i,j)
     &  +wind10_v(i,j)*wind10_v(i,j)))
        flvag=flcf3*xs

c- we2top wave erosion to top layer
c- separate flvag_dt (local var in this k-loop; dont touch flvag),
c- but keep fl_sum->fl_dt for compensation purposes
        flvag_dt=flvag*dticeb

        if ((flvag_dt).lt.0.) then
          write(iceberg_info_out_id,*) 
     &     'cdbr WARNING icebtraj.f; flvag_dt negative'
           flvag_dt = 0.0d0
        endif

        if ((flvag_dt).gt.wiceb(l)) then
          write(iceberg_info_out_id,*) 'cdbr WARNING icebtraj.f; flvag_dt>wiceb'
           flvag_dt = wiceb(l)
        endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--- Calculate Lateral melt and total melted volume per layer.
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        do k=ksubm,ks2
          dziceb=dz(k)
c- iceberg height per layer
          if (k.eq.kiceb(l)) dziceb=hiceb(l)+zw(k+1) ! dmr&mab: works only if zw is negative!! checked!
c- tip of the iceberg in bottom layer

          if(dziceb.lt.0) then  !should set dziceb to zero??? --- dmr&mab
            write(iceberg_info_out_id,*)'ALERT! dziceb(',l,'):',dziceb,
     &       '(k:',k,') at',n, 'in',i,j,'kiceb:',kiceb(l)
          endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Lateral melt ("Buoyant convection")
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c Lateral melt with a first and second order term according to El-Tahan 87:
c Rlatt=7.62e-3 Tw +1.29e-3 Tw^2 (m/day!!)

          tocean=scal(i,j,k,1)-tK_zero_C

          if(tocean.gt.0d0) then
            flat=flcf1*tocean+flcf2*tocean*tocean
          else
            flat=0.
          endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-- Total lateral melt per layer fl_dt into dVol_icb(ijk)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--- total lateral melt is fl_dt {m/day} (<wiceb):
          fl_sum=flat+flvag

          fl_dt=fl_sum*dticeb

          if(fl_dt.gt.wiceb(l)) then
            fl_dt=wiceb(l)
          endif

c--- chop basal melt into bottom layers ----------------------
          fb_dt_k=fb_dt

c- fb_dt_k is daily basal melt for layer k
          if (fb_dt.gt.dziceb.and.k.lt.ks2) then
c- limit it if melting into next layer (surface layer excluded)
            fb_dt_k=dziceb
	  endif

c- calculate leftovers for next layer:
          fb_dt=fb_dt-fb_dt_k

          if (fb_dt .lt. 0.) then
            write(iceberg_info_out_id,*) 
     &       'cdbr WARNING icebtraj.f; fb_dt leftover<0!?', fb_dt
            fb_dt = 0d0
          endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--- dVol_icb(i,j,k) is sum of basal and lateral melt minus overlap
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
cdmr&mab --- basal melting first
          dVol_icb(i,j,k) = dVol_icb(i,j,k)
     &      + pond_icb(l) * (
c- daily basal melt is bottom area times thickness fb_dt_k:
     &  fb_dt_k*1.5*wiceb(l)*wiceb(l)
cdmr&mab --- lateral melt on side 1.5*wiceb
     & + 1.5*wiceb(l)*MAX(dziceb-fb_dt_k,0.0)*fl_dt
cdmr&mab --- lateral melt on side (wiceb-fl_dt)
     & + MAX((wiceb(l)-fl_dt),0.0)*MAX(dziceb-fb_dt_k,0.0)*fl_dt
     &                                                    )
	   volmeltday(l,k) = volmeltday(l,k) + pond_icb(l) * (
     &  fb_dt_k*1.5*wiceb(l)*wiceb(l)
     & + 1.5*wiceb(l)*MAX(dziceb-fb_dt_k,0.0)*fl_dt
     & + MAX((wiceb(l)-fl_dt),0.0)*MAX(dziceb-fb_dt_k,0.0)*fl_dt
     &                                                    )

cdmr&mab ---
cdmr&mab --- So far, we are left with an iceberg that has a shape
cdmr&mab ---  at depth k of:
cdmr&mab --- shape(k) = (wiceb(l)-fl_dt)*(1.5*wiceb(l)-fl_dt)*MAX(dziceb-fb_dt_k,0.0)
cdmr&mab --- shape(k) = ( 1.5*wiceb**2 - 2.5*wiceb*fl_dt + fl_dt**2 )*(dziceb-fb_dt)
cdmr&mab ---

c- redirect wave erosion to top layer (ks2) from lower layers
cdmr&mab --- Wave erosion on two sides ...
          if(k.ne.ks2) then
            dVol_icb(i,j,ks2) = dVol_icb(i,j,ks2) 
     >      + pond_icb(l)*(
cdmr&mab --- lateral melt on side (1.5*wiceb-fl_dt)
     >      +  (1.5*wiceb(l)-fl_dt)*MAX(dziceb-fb_dt_k,0.0)*flvag_dt
cdmr&mab --- lateral melt on side (wiceb-fl_dt)-flvag_dt
     >      +  MAX(0.0d0,(wiceb(l)-fl_dt)-flvag_dt)*MAX(dziceb-fb_dt_k,0.0)
     >      *flvag_dt
     >                                                          )
	            volmeltday(l,ks2) = volmeltday(l,ks2) + pond_icb(l)*(
     >      +  (1.5*wiceb(l)-fl_dt)*MAX(dziceb-fb_dt_k,0.0)*flvag_dt
     >      +  MAX(0.0d0,(wiceb(l)-fl_dt)-flvag_dt)*MAX(dziceb-fb_dt_k,0.0)
     >      *flvag_dt
     >                                                          )
cdmr&mab --- suppress the water added to the top layer from the lower levels
            dVol_icb(i,j,k) = dVol_icb(i,j,k) 
     >      - pond_icb(l)*(
cdmr&mab --- lateral melt on side (1.5*wiceb-fl_dt)
     >      +  (1.5*wiceb(l)-fl_dt)*MAX(dziceb-fb_dt_k,0.0)
     >         *flvag_dt
cdmr&mab --- lateral melt on side (wiceb-fl_dt)-flvag_dt
     >      +  MAX(0.0d0,(wiceb(l)-fl_dt)-flvag_dt)*MAX(dziceb-fb_dt_k,0.0)
     >         *flvag_dt)
	            volmeltday(l,k) =volmeltday(l,k) - pond_icb(l)*(
     >      +  (1.5*wiceb(l)-fl_dt)*MAX(dziceb-fb_dt_k,0.0)
     >         *flvag_dt
     >      +  MAX(0.0d0,(wiceb(l)-fl_dt)-flvag_dt)*MAX(dziceb-fb_dt_k,0.0)
     >         *flvag_dt
     >                                                          )

          endif
cdmr&mab ---
cdmr&mab --- So far, we are left with an iceberg that has a shape
cdmr&mab ---  at depth k of:
cdmr&mab --- shape(k) = MAX(0.0,(wiceb-fl_dt-flvag_dt))*MAX(0.0,(1.5*wiceb-fl_dt-flvag_dt))*MAX(dziceb-fb_dt_k,0.0)
cdmr&mab ---

c----------------------------------------------------------------------
        end do ! end k-loop over layers

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c Deal with the "tip of the icebergs"
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
cdmr&mab ---  The emerged part of the iceberg is still in original size
cdmr&mab ---   (it has not been eroded or melted in the previous loop)

       if (fb_dt.GT.0.0d0) then
          tip_of_the_iceberg=1.5*wiceb(l)**2*fb_dt*pond_icb(l)
          dVol_icb(i,j,ks2) = dVol_icb(i,j,ks2)
     >    + tip_of_the_iceberg
          volmeltday(l,ks2) = volmeltday(l,ks2) + tip_of_the_iceberg
       endif

cdmr&mab ---
cdmr&mab --- Add all melting to the dVol_icb
cmab-not to be done (otherwise some melt counted twice!!!
cmab dVol_icb(i,j,:)=dVol_icb(i,j,:) + volmeltday(l,:)
cmab-not to be done (otherwise some melt counted twice!!!

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- REMOVE MELTED ICE VOLUME AND CALCULATE NEW SIZE
c----------------------------------------------------------------------

        if (kiceb(l) .gt. -1) then
           if ( hnew.lt.0d0) print*,'ALERT!!!hnew',hnew

           V_new(l) = MAX(0.0d0,(vol_in(l) - SUM(volmeltday(l,:))))
           if (V_new(l) .gt. 0d0 .and.hnew .gt.0d0 .and. pond_icb(l).gt.0d0)
     >       then
             wiceb(l) = SQRT(V_new(l) / (hnew * 1.5*pond_icb(l)))
             hiceb(l)=hnew/(1+remsim)

             if (hiceb(l).lt.0. .or. wiceb(l).lt.0) then
               hiceb(l)=max(hiceb(l),zero)
               wiceb(l)=max(wiceb(l),zero)
               kiceb(l)=ks2+1
               goto 790
            endif

          else
            hiceb(l) = 0d0
            wiceb(l)=0d0
            kiceb(l)=ks2+1
            goto 790
          endif
        else
          goto 800
        endif

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
          wiceb(l)=0d0
	  hiceb(l)=0d0
	  kiceb(l)=ks2+1
        else
          if (hiceb(l).le.5.) then
            dVol_icb(i,j,ks2)=dVol_icb(i,j,ks2)+pond_icb(l)*
     &        hiceb(l)*(1+remsim)*wiceb(l)*1.5*wiceb(l)
            volmeltday(l,ks2) = volmeltday(l,ks2)+ pond_icb(l)*
     &        hiceb(l)*(1+remsim)*wiceb(l)*1.5*wiceb(l)

            wiceb(l)=0d0
	    hiceb(l)=0d0
            kiceb(l)=ks2+1
          endif
        endif

c----------------------------------------------------------------------
c 	NEW VOLUME and "filling the icb-array"??
c- re-calcul du volume de l'iceberg "l" et remplissage des tableaux
c----------------------------------------------------------------------
c due to the (1+remsim) factor,
c the new underwater berg geometry is extrapolated to the above water part!!
c therefore we had to give the lateral melt per layer an
c (1+remsim) bonus to obtain the same total volume change...

        V_new(l)=wiceb(l)*wiceb(l)*1.5*hiceb(l)*(remsim+1)*pond_icb(l)
        if (V_new(l).lt.0.) V_new(l)=0d0

c---------------------------------------------------------------------
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- skip to here when avoiding melt:

 790  continue

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c Updating iceberg layer-depth kiceb(l)
c-	mise a jour de "kiceb" :
cc---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--- eliminate bergs that melted thru surface layer
        if (kiceb(l).ge.(ks2+1).or.hiceb(l).le.0.) then
          if(hiceb(l).ne.0d0) write(iceberg_info_out_id,*),
     &     'berg melted through the surface
     &                        eliminated! ', kiceb(l),hiceb(l)
          kiceb(l)=-1 ! ELIMINATED
          uiceb(l) = 0d0
          viceb(l) = 0d0
        endif

c--- eliminate bergs that melted to < 5m. width
        if (kiceb(l).gt.-1.and.wiceb(l).le.0.) then
          kiceb(l)=-1 ! ELIMINATED
          uiceb(l) = 0d0
          viceb(l) = 0d0
        endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  ROLL_OVER
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        if(kiceb(l).gt.-1) then
          wwiceb=wiceb(l)*1.5
          htot=(1+remsim)*hiceb(l)
          cfstab=(sqrt(cfst1+(cfst2/htot)))

          if ((wwiceb / htot).lt.cfstab) then
!          Berg rolls
            old_volume=(1+remsim)*hiceb(l)*wiceb(l)*1.5*wiceb(l)

c-- the geometry of this trafo is wrong: small bug!
            hiceb(l)=wiceb(l)/(1+remsim)
            wiceb(l)=sqrt(htot*wwiceb/1.5)
            wwiceb=wiceb(l)*1.5

            new_volume=(1+remsim)*hiceb(l)*wiceb(l)*1.5*wiceb(l)

            if ((new_volume-old_volume)/old_volume*100.gt.1E-10) then         
              write(iceberg_info_out_id,*) 
     &         'Iceberg volume lost by rolling',l,
     &              (new_volume-old_volume)/old_volume*100,
     &              'percentage of total volume'
            endif
          endif
        endif

c------ recalculate depth kiceb of all live -incl. grounded- icebergs
        if(kiceb(l).ge.0) then

          kiceb(l)=ks2
          do 795 k=ks2,2,-1
            if (-zw(kiceb(l)).le.hiceb(l)) then
              kiceb(l)=k-1
              if(kiceb(l) .eq. -1) print*,'kiceb -1!!'
            else
              goto 795
            endif
 795      continue
        endif


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- end of the main loop over all bergs 'l'
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        
 800  continue   
      
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      end subroutine icebtraj
