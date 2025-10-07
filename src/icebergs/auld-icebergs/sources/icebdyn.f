      subroutine icebdyn(i,j,l,un,vn,xx,yy)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  modif : 17/02/00
c- modif 24/10/19: cleaned use declaration sections              

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip             
      use const_mod, only: gpes, rho0
      use para0_mod, only: imax, jmax, kmax
      use para_mod, only:
      use bloc0_mod, only: kfs, ks1, dz, ks2, zw, tmu
     >  , uns2dx, uns2dy, eta
      use bloc_mod, only: smx, q, smy, fs2cor
      use ice_mod, only: uwind, vwind, hgbq, rhog, albq
      use iceberg_mod, only: seaicevel, wind10v, it_icb
     >   ,roiceb, ndbug, kiceb, wiceb, remsim, hiceb, dticeb
     >   ,wind10_u, wind10u, wind10_v, uiceb, viceb, xn, yn
     >   ,urep, vrep  
      use dynami_mod, only: ug, vg
      use reper_mod, only:
!! END_OF_USE_SECTION
     
      implicit none

!! START_OF_INCLUDE_SECTION

!PB#include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "reper.com"
! [SCRPTCOM] #include "iceberg.com"
! [SCRPTCOM] #include "ice.com"
! [SCRPTCOM] #include "dynami.com"

!! END_OF_INCLUDE_SECTION

c--variables locales:
      logical :: flag
C     double precision olduiceb,oldviceb

c--terme r=wave radiation, f=coriolis, w=water drag, b=barotrop, p=barocline, c=

c- dummy variables :
      real(kind=dblp), dimension(kmax) :: un, vn

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c JONO_out   --------------------------------------------------------------
c wrt_zero=1 prints empty grids involved in interpolation to 2058=fort.2058
c e=0.0001   is epsilon for determination of 'zero'; label is loop-var
c listi,j=   list of grids with wind='zero' (<e), involved in interpolation
c listi,jSI= list of grids with zero sea-ice speed
c aantal, aantalSI = counter for number of zero's in array list
      integer(kind=ip):: wrt_zero, listi(imax*jmax), listj(imax*jmax),
     >        listiSI(imax*jmax),listjSI(imax*jmax), aantal, aantalSI,
     >        i_xists(imax,jmax), i2, j2

!PB variables added after imposing implicit none
      integer(kind=ip) ::  l,n,ksubm,i,j,k,ii,jj,label,kountj,kounti
      real(kind=dblp) :: roair,cdragw,cdraga,Awater,Aair,hair,dziceb,
     >   unsmas,prclnx,prclny,xx,yy,pressx,pressy,prxdt, prydt,ff2cor,
     >   f2cdt,e,windu,windv,uvwd,wave, ccwave,frdtx,frdty,uvemd,cvar,
     >   cvelx,cvely,bvar, wvelx,wvely,uocean,vocean,velmd,bk,bb,
     >   dvar,dvelx,dvely,ui,vi,cdragi,dd,cdiag,unsdet, uu,vv,uunew,
     >   vvnew,olduiceb,oldviceb

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      aantal=1
      aantalSI = 0
      seaicevel = 0
       wind10v=0
      wrt_zero=1
      if (((aantal.gt.l).and.(it_icb.eq.1)).or.aantal.lt.0) then
        aantal = 0
        aantalSI = 0
      endif
c JONO_out end   ----------------------------------------------------------
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c--constantes utilisees:
      roiceb=910d0
      roair=1.275d0

c---nouveau cdrag air (doublage du vent)
c     cdraga=1.3d0
C      cdraga=5.2d0
c       cdraga=2.6d0
      cdragw=0.9d0
c JONO Bigg says
c JONO_pw3
      cdraga=1.3d0
c     cdragw=0.9d0
c JONO 11-6-04 for consistency, let's take the same dragcoefficient as the wind-water interaction: 2.1:
c      cdraga=2.1d0

c JONO 11-6 Bigg suggested for the effective area A:
c Awater= factor 1 (head-on) for water and
c Aair=   factor 1.77 = |1.5sin(theta)| + |cos(theta)| , with theta=45% "in accordance with Eckham theory"

C JONO_4allCR 29-9 having corrected the error (factor rho0) now switching the A's on again
C JONO 24-9 5noSeadrag using Awater=0 to shut off the seawaterdrag
C JONO_4onWR 28-9 setting effective areas to zero except for the waveradiation, to obtain' only WR'
      Awater=1.0d0
C      Awater=0.0
C JONO 24-9 5noW using Aair=0 to shut off the winddrag this affects the waveradiation and the winddrag (and resetting Awater to1)
      Aair=1.768d0
C      Aair=0.0
C JONO 27-9 4noWR A... reset to Spw3 values, now shutting off ONLY the wave radiation (see below)(OUTDATED SEE 4onWR 28-9)

      n=it_icb
      flag = (n.eq.ndbug).and.(mod(l,10000).lt.1)

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Calcul des forces:
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c JONO_says ksubm is submerged part of the berg? (in layer# or meters?)
c in layers because it's the maximum of kiceb(l) and bottomdepth kfs
c kiceb is iteratively set to match the height hiceb? in deficeb:376
      ksubm=max(kiceb(l),kfs(i,j))

c-calcul de la hauteur non submergee
c JONO_says hair=the part above the water?
c if ksubm<2, there is no part above the water?
c (maybe should change this to also use hair for stranded bergs??)
c i don't understand why remsim is not used. and what does the layer thickness dz have to do with it?
c maybe this is only for stranded bergs...?
      hair=0.
      dziceb=0.
      if ((ksubm-1).ge.1) then
        do 100 k=ks1,(ksubm-1)
          hair=hair+dz(k)
 100    continue
      endif
      if (hair.le.(0.)) hair=0.

c- masse icebergs: SIDE=1.5*FRONT...
      unsmas = 1./(roiceb*1.5*wiceb(l)*wiceb(l)*(remsim+1.)*hiceb(l))
c      write(*,*) 'mass',unsmas

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Pression barocline :

      prclnx = 0.
      prclny = 0.
      do 130 k=ksubm,ks2
        dziceb=dz(k)
        if (k.eq.kiceb(l)) dziceb=hiceb(l)+zw(k+1)
c-
        ii=i
        jj=j
        prclnx = prclnx + dziceb*(1.-xx)*(1.-yy)*
     &     tmu(ii,jj,k)*smx(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               +(q(ii-1,jj,k)-q(ii,jj-1,k)))
        prclny = prclny + dziceb*(1.-xx)*(1.-yy)*
     &     tmu(ii,jj,k)*smy(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               -(q(ii-1,jj,k)-q(ii,jj-1,k)))
c-
        ii=i+1
        jj=j
        prclnx = prclnx + dziceb*xx*(1.-yy)*
     &     tmu(ii,jj,k)*smx(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               +(q(ii-1,jj,k)-q(ii,jj-1,k)))
        prclny = prclny + dziceb*xx*(1.-yy)*
     &     tmu(ii,jj,k)*smy(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               -(q(ii-1,jj,k)-q(ii,jj-1,k)))
c-
        ii=i
        jj=j+1
        prclnx = prclnx + dziceb*(1.-xx)*yy*
     &     tmu(ii,jj,k)*smx(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               +(q(ii-1,jj,k)-q(ii,jj-1,k)))
        prclny = prclny + dziceb*(1.-xx)*yy*
     &     tmu(ii,jj,k)*smy(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               -(q(ii-1,jj,k)-q(ii,jj-1,k)))
c-
        ii=i+1
        jj=j+1
        prclnx = prclnx + dziceb*xx*yy*
     &     tmu(ii,jj,k)*smx(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               +(q(ii-1,jj,k)-q(ii,jj-1,k)))
        prclny = prclny + dziceb*xx*yy*
     &     tmu(ii,jj,k)*smy(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               -(q(ii-1,jj,k)-q(ii,jj-1,k)))
 130  continue

C      if ((hiceb(l).gt.0).AND.(unsmas.gt.0)) then
C             termexp(l)=prclnx*uns2dx/(hiceb(l)*unsmas*rho0)
C             termeyp(l)=prclny*uns2dy/(hiceb(l)*unsmas*rho0)
C      endif

      prclnx = prclnx*uns2dx*dticeb/hiceb(l)
      prclny = prclny*uns2dy*dticeb/hiceb(l)
c      write(*,*) 'pressx',prclnx,prclnx,uns2dx,dticeb,xx,yy
c      write(*,*) 'pressy',prclny

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Pression barotrope :
      pressx=0.
      pressy=0.
      k = ks2

c JONO_dbg JONO_CHECK REDUNDANT? these are not commons nor used...
c they are calculated again and used in icebtraj, in an interpolation i don't understand...
c 14 april 2004: commenting them out
c      yy1=1.+tms(i,j-1,k)*(tms(i,j+1,k)*yy-1.)
c      xx1=1.+tms(i-1,j,k)*(tms(i+1,j,k)*xx-1.)
c JONO_dbg JONO_CHECK REDUNDANT JONO_end
c
c JONO in these 'vector' interpolations the q -OR eta-values from a shift to the lower left
c JONO are used with the mask of a shift to the upper right...
c JONO also note the sign differences before (ii-1,jj)...
c JONO might not be a shift but a corner interpolation!

      ii=i
      jj=j
      pressx=pressx+smx(ii,jj,3)*tmu(ii,jj,k)*(1.-xx)*(1.-yy)*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         + (eta(ii-1,jj)-eta(ii,jj-1)) )
      pressy=pressy+smy(ii,jj,3)*tmu(ii,jj,k)*(1.-xx)*(1.-yy)*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         - (eta(ii-1,jj)-eta(ii,jj-1)) )
C     if (flag) write (6,'(A,I8,I4,2F11.4)')
C    &           ' dyn: pressx,y 1: ', n, l , pressx,pressy

      ii=i+1
      jj=j
      pressx=pressx+smx(ii,jj,3)*tmu(ii,jj,k)*xx*(1.-yy)*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         + (eta(ii-1,jj)-eta(ii,jj-1)) )
      pressy=pressy+smy(ii,jj,3)*tmu(ii,jj,k)*xx*(1.-yy)*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         - (eta(ii-1,jj)-eta(ii,jj-1)) )
C     if (flag) write (6,'(A,I8,I4,2F11.4)')
C    &           ' dyn: pressx,y 2: ', n, l , pressx,pressy

      ii=i
      jj=j+1
      pressx=pressx+smx(ii,jj,3)*tmu(ii,jj,k)*(1.-xx)*yy*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         + (eta(ii-1,jj)-eta(ii,jj-1)) )
      pressy=pressy+smy(ii,jj,3)*tmu(ii,jj,k)*(1.-xx)*yy*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         - (eta(ii-1,jj)-eta(ii,jj-1)) )
C     if (flag) write (6,'(A,I8,I4,2F11.4)')
C    &           ' dyn: pressx,y 3: ', n, l , pressx,pressy

      ii=i+1
      jj=j+1
      pressx=pressx+smx(ii,jj,3)*tmu(ii,jj,k)*xx*yy*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         + (eta(ii-1,jj)-eta(ii,jj-1)) )
      pressy=pressy+smy(ii,jj,3)*tmu(ii,jj,k)*xx*yy*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         - (eta(ii-1,jj)-eta(ii,jj-1)) )
C     if (flag) write (6,'(A,I8,I4,2F11.4)')
C    &           ' dyn: pressx,y 4: ', n, l , pressx,pressy

      prxdt = pressx*gpes*uns2dx*dticeb
      prydt = pressy*gpes*uns2dy*dticeb

C     if (flag) write (6,'(2(A,2F10.6))') ' dyn: press.eta,q _X:',
C    &          prxdt, prclnx, ' ; _Y:', prydt, prclny

C      if (unsmas.gt.0) then
C             termexb(l)=pressx*uns2dx*gpes/(unsmas*rho0)
C             termeyb(l)=pressy*uns2dy*gpes/(unsmas*rho0)
C      endif

      prxdt = prxdt + prclnx
      prydt = prydt + prclny
c      write(*,*) 'pressxcomp',prxdt
c      write(*,*) 'pressycomp',prydt

c---pression equilibre hydrostatique:le calcul s'effectue apres le calcul de la
c JONO_says ???
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Coriolis:
c      ff2cor = ( (1.-yy)*((1.-xx)*fs2cor(i,j)+xx*fs2cor(i+1,j))
c     &         +  yy*((1.-xx)*fs2cor(i,j+1)+xx*fs2cor(i+1,j+1)) )
      ff2cor = ( (1.-yy)*( (1.-xx)*tmu(i,j,ks2)*fs2cor(i,j)
     &                   +  xx*tmu(i+1,j,ks2)*fs2cor(i+1,j) )
     &         +  yy*( (1.-xx)*tmu(i,j+1,ks2)*fs2cor(i,j+1)
     &               +  xx*tmu(i+1,j+1,ks2)*fs2cor(i+1,j+1) ) )
      f2cdt = ff2cor*dticeb
c      write(*,*) 'coriolis,f2cdt',f2cdt
c JONO_dbg JONO_CHECK says: I think this interpolation interprets empty/land cells as ZERO instead of ignoring them

c JONO_dbg 14 april 2004 adding write corr-map:
       if(n.eq.60.and.l.eq.10) then
          write (1059, '(A)') 'CORIOLIS-map,fs2cor(ij):'
          do 133 j2=1,jmax
            do 132 i2=1,imax
              write (1059,'(128F6.4)') fs2cor(i2,j2)
 132  continue
 133  continue
       endif
       if(n.eq.60.and.l.eq.10) then
          write (1059, '(A)') 'CORIOLIS-map,f2cdt:'
          do j2=1,jmax
            do i2=1,imax
              write (1059,'(128E8.2)') fs2cor(i2,j2)
            enddo
          enddo
       endif
c JONO_end

C     if (flag) write (6,'(A,I8,I4,3F10.6)')
C    &     ' dyn: pr_x,y, f_cor/2 : ', n, l , prxdt,prydt, f2cdt

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- masse:
c modif 2/03/00 la masse est calculee plus haut
c      unsmas = 1./(roiceb*1.5*wiceb(l)*wiceb(l)*(remsim+1.)*hiceb(l))

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Air drag + vagues:
c******************************************************************************
c JONO_icb_cpl_wind
c introduced new wind vector [wind10_u(i,j),wind10_v]
c will be interpolated to windu,windv at icb-position here:
c p.s.
c The JONO_out lines in this section list the gridcells with zero wind,
c  used (erronously) in the interpolation. These are probably all land-cells.
c  The interpolation is therefor limited to wind-containing cells.
c  PROBLEM: 3 point interpolation?? .

c***** JONO_dbg JONO_CHECK
c*****   maybe this example from forcing.f helps:
c*****
c*****  Does not use coastal points for determining wind
c*****  at oceanic grid.
c*****
c*****      do j=js1,js2
c*****        jp1 = j+1
c*****        do i=is1(j),is2(j)
c*****          ip1         = (i+1)
c*****          usp         = 1.0/
c*****     &    max(1.0e-12*one,tmu(i,j,ks2)+tmu(ip1,j,ks2)
c*****     &                   +tmu(i,jp1,ks2)+tmu(ip1,jp1,ks2))
c*****          uwind(i,j)  = (ugw(i,j)*tmu(i,j,ks2)+
c*****     &                   ugw(ip1,j)*tmu(ip1,j,ks2)+
c*****     &                   ugw(i,jp1)*tmu(i,jp1,ks2)+
c*****     &                   ugw(ip1,jp1)*tmu(ip1,jp1,ks2))*usp
c*****          vwind(i,j)  = (vgw(i,j)*tmu(i,j,ks2)+
c*****     &                   vgw(ip1,j)*tmu(ip1,j,ks2)+
c*****     &                   vgw(i,jp1)*tmu(i,jp1,ks2)+
c*****     &                   vgw(ip1,jp1)*tmu(ip1,jp1,ks2))*usp
c*****        enddo
c*****      enddo
c***********************************************************************
c***** with ugw, vgw: geostrophic wind at the corners of the grid squares
c*****--3.1. Calculate geostrophic wind at the corners of the grid squares.
c---------------------------------------------------------------------
c*****
c*****      rhlim=2.0*omega*sin(15.0*radian)
c*****      do 330 j=2,jmax
c*****          jm1 = j-1
c*****        do 320 i=1,imax
c*****          im1      =  (i-1)+(imax-2)*(1/i)
c*****          rhoa     =  psbq(i,j)/(287.04*tabq(i,j))
C*****         rhoafn   =  rhoa*zfn(i,j)+1d-18
c*****          rhoafn   =  rhoa*sign(one,zfn(i,j))*max(abs(zfn(i,j)),rhlim)
C*****         write(89,*) i,j,rhoa,zfn(i,j)
c*****          ugw(i,j) = -(-(alambd(i,j,1,1,2,1)+
c*****     &                       alambd(i,j,1,2,2,1))*psl(i,jm1)-
c*****     &                 (alambd(i,j,1,1,1,1)
c*****     &                      +alambd(i,j,1,2,1,1))*psl(im1,jm1)+
c*****     &                 (alambd(i,j,1,1,2,2)
c*****     &                      -alambd(i,j,1,2,2,2))*psl(i,j)+
c*****     &                 (alambd(i,j,1,1,1,2)
c*****     &                      -alambd(i,j,1,2,1,2))*psl(im1,j))/
c*****     &                rhoafn
c*****
c*****          vgw(i,j) =  ((alambd(i,j,2,2,2,1)
c*****     &                     -alambd(i,j,2,1,2,1))*psl(i,jm1)+
c*****     &                 (alambd(i,j,2,2,2,2)
c*****     &                     -alambd(i,j,2,1,2,2))*psl(i,j)-
c*****     &                 (alambd(i,j,2,2,1,1)
c*****     &                     +alambd(i,j,2,1,1,1))*psl(im1,jm1)-
c*****     &                 (alambd(i,j,2,2,1,2)
c*****     &                     +alambd(i,j,2,1,1,2))*psl(im1,j))/
c*****     &                rhoafn
c*****320     continue
c*****330   continue
c*****
c***** ,with
c***** alambd { usden,d1d2p(i,j),dx2(i,j),dx1(i,j),h1(i,j),h2(i,j) }, see defgrid
c***** forcing.f:
c*****      psl(i,j)=((rpr(i,j,1,mft)-rpr(i,j,1,mit))/dtt-
c*****     &          ((3.0*ttat*ttat-1.0)*rpr(i,j,2,mit)-
c*****     &          (3.0*ttbt*ttbt-1.0)*rpr(i,j,2,mft))*dtts6+
c*****     &          ampr(i,j))*100.0
c*****, with        read(71)(/72) (rpr(i,j,1 (/2),1),i=1,imax)
c*****
c**********************************************************************
c***** JONO_dbg JONO_CHECK end

Cpw    windu=((1.-yy)*((1.-xx)*wind10_u(i,j)+xx*wind10_u(i+1,j))
Cpw   &     + yy*((1.-xx)*wind10_u(i,j+1)+xx*wind10_u(i+1,j+1)))
Cpw    windv=((1.-yy)*((1.-xx)*wind10_v(i,j)+xx*wind10_v(i+1,j))
Cpw   &     + yy*((1.-xx)*wind10_v(i,j+1)+xx*wind10_v(i+1,j+1)))
c JONO_pure_wind 14 april 2004 COMMENTED THESE LINES and CHANGED windu,v TO wind10_u,v(i,j)
c ALSO commenting all write statements that write windu
c To solve the 'interpolating from empty land-grids' problem,
c and to save cpu-time, gonna use the not-interpolated wind instead!
c I think i can just use the wind10_u,v at (i,j), won't ever be land? CHECK!
c
c_pw  windu = wind10_u(i,j)
c_pw  windv = wind10_v(i,j)
c (still COMMENTED because refering straight to wind10_u,v(i,j))

c JONO_dbg
c      if (flag) write (6,'(2(A,I6),2(A,2F7.3))')
c     &      'Wind: n:',n,' l:',l,' windu,v: (',windu,windv,
c     &      ') wind10_u,v:', wind10_u(i,j), wind10_v(i,j)

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c ------------------------------------------------------------------------------
c JONO_out ---------------------------------------------------------------------
c5 JONO5 cutting the writing of zero's in neighbouring gridcells, adding wind10_v(i,j):
      if (wrt_zero.eq.1) then
      e=0.0001
c _u; i+1
c5      if  (abs(wind10_u(i+1,j)).lt.e) then
c5        do label=0,(1+aantal)
c5cCheck if allready in the list    write(1059,*) 'doing',n,l,label,aantal
c5           if( ( listi(label).eq.(i+1) ).and.
c5     &         ( listj(label).eq.(j  ) ) ) then
c5c else list it (note difference)  write(1059,*) 'dupl in',n,l,label,aantal,i+1,j
c5             goto 135
c5           endif
c5        enddo
c5        aantal = aantal + 1
c5        listi(aantal) = i+1
c5        listj(aantal) = j
c5        write (2058,'(A,2I4,A,6I7)')
c5     &  'wind zero in gridcell (',listi(aantal),listj(aantal),
c5     &  ') n,l,label,aantal,i,j:', n, l, label,aantal,i,j
c5      endif
c5c j+1
c5      if  (abs(wind10_u(i,j+1)).lt.e) then
c5        do label=0,(1+aantal)
c5           if(( listi(label).eq.(i  ) ).and.
c5     &        ( listj(label).eq.(j+1) )) goto 135
c5        enddo
c5        aantal = aantal + 1
c5        listi(aantal) = i
c5        listj(aantal) = j+1
c5        write (2058,'(A,2I4,A,6I7)')
c5     &  'wind zero in gridcell (',listi(aantal),listj(aantal),
c5     &  ') n,l,label,aantal,i,j:', n, l, label,aantal,i,j
c5      endif
c5c i+1, j+1
c5      if  (abs(wind10_u(i+1,j+1)).lt.e) then
c5        do label=0,(1+aantal)
c5           if(( listi(label).eq.(i+1) ).and.
c5     &        ( listj(label).eq.(j+1) )) goto 135
c5        enddo
c5        aantal = aantal + 1
c5        listi(aantal) = i+1
c5        listj(aantal) = j+1
c5        write (2058,'(A,2I4,A,6I7)')
c5     &  'wind zero in gridcell (',listi(aantal),listj(aantal),
c5     &  ') n,l,label,aantal,i,j:', n, l, label,aantal,i,j
c5      endif
c i,j
      if  (abs(wind10_u(i,j)).lt.e) then
        do label=1,(1+aantal)
           if(( listi(label).eq.(i  ) ).and.
     &        ( listj(label).eq.(j  ) )) goto 135
        enddo
        aantal = aantal + 1
        listi(aantal) = i
        listj(aantal) = j
        wind10u = wind10u + 1
!        write (2058,'(A,2I4,A,6I7)')
!     &  'wind10_u zero in gridcell (',listi(aantal),listj(aantal),
!     &  ') n,l,label,aantal,i,j:', n, l, label,aantal,i,j
      endif

      if  (abs(wind10_v(i,j)).lt.e) then
        do label=1,(1+aantal)
           if(( listi(label).eq.(i  ) ).and.
     &        ( listj(label).eq.(j  ) )) goto 135
        enddo
        aantal = aantal + 1
        listi(aantal) = i
        listj(aantal) = j
        wind10v = wind10v + 1
!        write (2058,'(A,2I4,A,6I7)')
!     &  'wind10_v zero in gridcell (',listi(aantal),listj(aantal),
!     &  ') n,l,label,aantal,i,j:', n, l, label,aantal,i,j
      endif

 135  continue
c     &     (abs(wind10_v(i+1,j)).lt.e).or.
c     &     (abs(wind10_v(i,j+1)).lt.e).or.
c     &     (abs(wind10_v(i+1,j+1)).lt.e).or.
c     &     (abs(wind10_v(i,j)).lt.e) ) then
      endif
c end "if(wrt_zero=T)"  ---- end JONO_out -----------------------------
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      if (flag) then
      e=0.0001
      if ( (abs(wind10_u(i+1,j)).lt.e).or.
     &     (abs(wind10_u(i,j+1)).lt.e).or.
     &     (abs(wind10_u(i+1,j+1)).lt.e).or.
     &     (abs(wind10_u(i,j)).lt.e).or.
     &     (abs(wind10_v(i+1,j)).lt.e).or.
     &     (abs(wind10_v(i,j+1)).lt.e).or.
     &     (abs(wind10_v(i+1,j+1)).lt.e).or.
     &     (abs(wind10_v(i,j)).lt.e) ) then
       write (1059,'(A,E10.2,A,I7,A,I7,A,2I7,A)')
     &  'WIND smaller than e=',e,'   icb',l,' n:',n,' grid(',i,j,')'
       write (1059, '(4(A,E10.2))')
     &     '  wind10_u(i+1,j): ', wind10_u(i+1,j),
     &     '  wind10_u(i,j+1): ', wind10_u(i,j+1),
     &     ' wind10_u(i+1,j+1): ', wind10_u(i+1,j+1),
     &     '    wind10_u(i,j): ', wind10_u(i,j),
     &     '  wind10_v(i+1,j): ', wind10_v(i+1,j),
     &     '  wind10_v(i,j+1): ', wind10_v(i,j+1),
     &     ' wind10_v(i+1,j+1): ', wind10_v(i+1,j+1),
     &     '    wind10_v(i,j): ', wind10_v(i,j)
      endif
      endif

      if (flag)
     &   write (1059,'(/,A,I5,A,2I4,A)')
     &  'id:',l,' cpl MINUS fix : u,v-WINDS at grid i,j(',i,j,')'
      if (flag) then
        write (1059,'(A,F6.2,A,F5.2,A1,F6.3,A,F7.3,A)')
     &  ' w10u: [',
     & wind10_u(i,j),'] - ',uwind(i,j),'=',(wind10_u(i,j)-uwind(i,j)),
     &  '  and local interpol.: [' ,windu , ']'
       write (1059,'(A,F6.2,A,F5.2,A1,F6.3,A,F7.3,A)')
     &  ' w10v: [',
     & wind10_v(i,j),'] - ',vwind(i,j),'=',(wind10_v(i,j)-vwind(i,j)),
     &  '                       [' ,windv , ']'
      endif

c JONO_fix_wind This was the original fixed wind:
c JONO_fix_wind      windu=((1.-yy)*((1.-xx)*uwind(i,j)+xx*uwind(i+1,j))
c JONO_fix_wind     &     + yy*((1.-xx)*uwind(i,j+1)+xx*uwind(i+1,j+1)))
c JONO_fix_wind      windv=((1.-yy)*((1.-xx)*vwind(i,j)+xx*vwind(i+1,j))
c JONO_fix_wind     &     + yy*((1.-xx)*vwind(i,j+1)+xx*vwind(i+1,j+1)))
c JONOend
c******************************************************************************

c JONO_REDUNDANT leftovers i think
c     write(*,*) 'windu et windv',windu,windv
c--essai de doublage du vent
C      windu=2*windu
C      windv=2*windv

c---wave radiation:

c JONO Smith;"for waves absorbed at a deep vertical wall ":
c wave radiation (3.15):
c        Fr = 1/4 * RHOw g a^2 L (Va/|Va|)
c, with Va wind speed (vector)
c, a = 01010125*|Va|^2 wave amplitude
c, L length of the face normal to incoming waves ( CHECK: is factor 1.5 0r 1 reasonable)
C       JONO 11-6 using Aair-1.768
c CHECK: shouldn't I us the vector difference between wind and iceb speed???
c (i think in this case might not be needed???)
c JONO 29-9
c dticb=86400=#sec/day
c unsmas=1/icbmass=w*1.5*w*h*(1+remsim) (how about the density!=> it will cancel against rho0, so I will remove rho0 below(=factor 1000!) )
c rho0=1030 (kg/m3)
c gpes=9.8

c JONO_pure_wind: change to
      uvwd = wind10_u(i,j)*wind10_u(i,j) + wind10_v(i,j)*wind10_v(i,j)
c      uvwd = windu*windu + windv*windv

      wave = 0.010125d0*uvwd
c JONO_REDUNDANT,.5*0.02=0.01                           wave = 0.5*uvwd*0.02025d0
c JONO_dbg CHECK REDUNDANT, not mentioned anywhere     	uvcoef = uvcoef*uvcoef

cJONO 11-6       ccwave = 0.25*dticeb*unsmas*rho0*gpes*1.5*wiceb(l)*
cJONO 11-6     &                            wave*wave/sqrt(uvwd)

C JONO_4onWR setting Aair=1.678 only for the waveradiation
C      ccwave = 0.25*dticeb*unsmas*rho0*gpes*Aair*wiceb(l)*
C     &                            wave*wave/sqrt(uvwd)
C JONO_4onWRcr JONO_dbg 29-9 rho0 cancels from the mass
C      ccwave = 0.25*dticeb*unsmas*rho0*gpes*1.768*wiceb(l)*
C     &                            wave*wave/sqrt(uvwd)
C JONO_4allCR canceld rho0 and reset to Aair to obtain CORRECT RUN
      ccwave = 0.25*dticeb*unsmas*gpes*Aair*wiceb(l)*
     &                            wave*wave/sqrt(uvwd)

C JONO_4noWR setting frdtx,y to zero to shut off the waveradiation:
C      frdtx = 0.0
C      frdty = 0.0
c JONO_pure_wind: change to
      frdtx = ccwave*wind10_u(i,j)
      frdty = ccwave*wind10_v(i,j)
c      frdtx = ccwave*windu
c      frdty = ccwave*windv
c JONO explains frdtx,y = dt*MASS*Frx,y = acceleration

c JONO_REDUNDANT leftovers i think
c      termexr(l)=0.25*rho0*gpes*1.5*wiceb(l)*windu*
c     &                            wave*wave/sqrt(uvwd)
c      termeyr(l)=0.25*rho0*gpes*1.5*wiceb(l)*windv*
c     &                            wave*wave/sqrt(uvwd)

c---air drag:

c JONO dynamic equation: (airdrag)
c        Fa=1/2 Rhoa Ca Aa |va-vi| (va-vi)
c, with Ca=5.2 {=1.3 (Smith93) *4 "we'r doubling the windstrength"}
c, Aa is the "section transversale a la direction de la tension exercee
c par la vent"  (CHECK: why 1/4 in waverad and 1/2 in other drag?)
c BIGGs/Smith:
c "cross-sectional area (above the water line) in a vertical plane normal to the wind"

c JONO_pure_wind: change windu,v to wind10_u,v(i,j)
      uvemd=sqrt((wind10_u(i,j)-uiceb(l))*(wind10_u(i,j)-uiceb(l))
     &            +(wind10_v(i,j)-viceb(l))*(wind10_v(i,j)-viceb(l)))
c      uvemd=sqrt((windu-uiceb(l))*(windu-uiceb(l))
c     &            +(windv-viceb(l))*(windv-viceb(l)))

c JONO 11-6      cvar = 0.5*dticeb*unsmas*cdraga*roair*1.5*wiceb(l)*
c JONO 11-6     &                         remsim*hiceb(l)*uvemd
      cvar = 0.5*dticeb*unsmas*cdraga*roair*Aair*wiceb(l)*
     &                         remsim*hiceb(l)*uvemd

c JONO_pure_wind: change windu,v to wind10_u,v(i,j)
      cvelx=cvar*( wind10_u(i,j)-uiceb(l) )
      cvely=cvar*( wind10_v(i,j)-viceb(l) )
c      cvelx=cvar*windu
c      cvely=cvar*windv
c JONO_CHECK why not windu-uiceb??? same goes for frdtx. doesn't this cause flying icebergs?

c JONO_REDUNDANT leftovers i think
c     termexc(l)=0.5*cdraga*roair*1.5*wiceb(l)*remsim*hiceb(l)*
c    &              uvemd*(windu-uiceb(l))/10e3
c     termeyc(l)=0.5*cdraga*roair*1.5*wiceb(l)*remsim*hiceb(l)*
c    &              uvemd*(windv-viceb(l))/10e3

C     if (flag) write (6,'(A,I8,I4,1P3E11.3)')
C    &      'Wind: n,l, windu,v, uvemd:', n,l, windu,windv, uvemd
      if (flag) write(6,'(A,1P5E11.3)')
     & 'A.drag+vagues: C,C_x,y, frdt_x,y=',
     &                  cvar,cvelx,cvely,frdtx,frdty


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- water drag:
      bvar=0.
      wvelx=0.
      wvely=0.

c      write(*,*) 'uocean',uocean

      uocean = 0.
      vocean = 0.

c      write(*,*) 'puoceanp',uocean

c     termexw(l)=0.
c     termeyw(l)=0.

      do 140 k=ksubm,ks2
        dziceb=dz(k)
        if (k.eq.kiceb(l)) dziceb=hiceb(l)+zw(k+1)
        velmd=sqrt((un(k)-uiceb(l))*(un(k)-uiceb(l))
     &            +(vn(k)-viceb(l))*(vn(k)-viceb(l)))
C       bk=cdragw*rho0*1.5*wiceb(l)*dziceb*velmd*dticeb*unsmas
        bk = dziceb*velmd
        bvar=bvar+bk
        wvelx=wvelx+bk*un(k)
        wvely=wvely+bk*vn(k)
c       termexw(l)=termexw(l)+(bk*(un(k)-uiceb(l)))
c       termeyw(l)=termeyw(l)+(bk*(vn(k)-viceb(l)))

        uocean = uocean + dziceb*un(k)
        vocean = vocean + dziceb*vn(k)
 140  continue
c JONO 11-6      bb = 0.5*dticeb*unsmas*cdragw*rho0*1.5*wiceb(l)
      bb = 0.5*dticeb*unsmas*cdragw*rho0*Awater*wiceb(l)
      bvar = bb*bvar
      wvelx = bb*wvelx
      wvely = bb*wvely
c     termexw(l)=termexw(l)*0.5*cdragw*rho0*1.5*wiceb(l)/10e3
c     termeyw(l)=termeyw(l)*0.5*cdragw*rho0*1.5*wiceb(l)/10e3

      uocean = uocean/hiceb(l)
      vocean = vocean/hiceb(l)
c      write(*,*) 'bvar',bvar

c--calcul pression equilibre hydrostat(le 2 parceque coriolis est predivise par
c---------------------------------------------
c      prxdt=-dticeb*ff2cor*2*vocean
c      prydt=(dticeb*ff2cor*2*uocean)

C       if (flag) write (6,'(A,4F13.6)')
C    &   ' dyn: velmd, B, wx,y=', velmd, bvar,wvelx,wvely

c---+-S--1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- sea ice drag:

      dvar=0.
      dvelx=0.
      dvely=0.

c---determination du niveau de glace de mer
c---ui=vitesse de la glace
c---vi=vitesse de la glace
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      ui= ( (1.-yy)*((1.-xx)*ug(i,j)+xx*ug(i+1,j))
     &         +  yy*((1.-xx)*ug(i,j+1)+xx*ug(i+1,j+1)) )
      vi= ( (1.-yy)*((1.-xx)*vg(i,j)+xx*vg(i+1,j))
     &         +  yy*((1.-xx)*vg(i,j+1)+xx*vg(i+1,j+1)) )
      cdragi=.9

c*******************************************************************************
c JONO_out ---------------------------------------------------------------------
      if (wrt_zero.eq.1) then
      e=0.000001
c _u; i+1
      if  (abs(ug(i+1,j)).lt.e) then
        do label=1,(1+aantalSI)
cCheck if allready in the list    write(1059,*) 'doing',n,l,label,aantalSI
           if( ( listiSI(label).eq.(i+1) ).and.
     &         ( listjSI(label).eq.(j  ) ) ) then
c else list it (note difference)  write(1059,*) 'dupl in',n,l,label,aantalSI,i+1,j
             goto 150
           endif
        enddo
        aantalSI = aantalSI + 1
        listiSI(aantalSI) = i+1
        listjSI(aantalSI) = j
        seaicevel = seaicevel + 1
!        write (2058,'(A,2I4,A,6I7)')
!     &'SEA-ICE speed zero in gridcell (',listiSI(aantalSI),listjSI(aantalSI),
!     &  ') n,l,label,aantalSI,i,j:', n, l, label,aantalSI,i,j
      endif
c j+1
      if  (abs(ug(i,j+1)).lt.e) then
        do label=1,(1+aantalSI)
           if(( listiSI(label).eq.(i  ) ).and.
     &        ( listjSI(label).eq.(j+1) )) goto 150
        enddo
        aantalSI = aantalSI + 1
        listiSI(aantalSI) = i
        listjSI(aantalSI) = j+1
        seaicevel = seaicevel + 1
!        write (2058,'(A,2I4,A,6I7)')
!     &'SEA-ICE speed zero in gridcell (',listiSI(aantalSI),listjSI(aantalSI),
!     &  ') n,l,label,aantalSI,i,j:', n, l, label,aantalSI,i,j
      endif
c i+1, j+1
      if  (abs(ug(i+1,j+1)).lt.e) then
        do label=1,(1+aantalSI)
           if(( listiSI(label).eq.(i+1) ).and.
     &        ( listjSI(label).eq.(j+1) )) goto 150
        enddo
        aantalSI = aantalSI + 1
        listiSI(aantalSI) = i+1
        listjSI(aantalSI) = j+1
        seaicevel = seaicevel + 1
!        write (2058,'(A,2I4,A,6I7)')
!     &'SEA_ICE speed zero in gridcell (',listiSI(aantalSI),listjSI(aantalSI),
!     &  ') n,l,label,aantalSI,i,j:', n, l, label,aantalSI,i,j
      endif
c i,j
      if  (abs(ug(i,j)).lt.e) then
        do label=1,(1+aantalSI)
           if(( listiSI(label).eq.(i  ) ).and.
     &        ( listjSI(label).eq.(j  ) )) goto 150
        enddo
        aantalSI = aantalSI + 1
        listiSI(aantalSI) = i
        listjSI(aantalSI) = j
        seaicevel = seaicevel + 1
!        write (2058,'(A,2I4,A,6I7)')
!     &'SEA-ICE speed zero in this grd (',listiSI(aantalSI),listjSI(aantalSI),
!     &  ') n,l,label,aantalSI,i,j:', n, l, label,aantalSI,i,j
      endif
 150  continue

c     &     (abs(wind10_v(i+1,j)).lt.e).or.
c     &     (abs(wind10_v(i,j+1)).lt.e).or.
c     &     (abs(wind10_v(i+1,j+1)).lt.e).or.
c     &     (abs(wind10_v(i,j)).lt.e) ) then

c JONO_out print SORTED map
      if (it_icb.eq.359.and.l.eq.1) then
       do kountj=1,jmax
        do kounti=1,imax
         do label=1,(1+aantalSI)
           i_xists(kounti,kountj)=0
           if(( listiSI(label).eq.kounti ).and.
     &        ( listjSI(label).eq.kountj )) then
             i_xists(kounti,kountj)=1
             goto 160
           endif
         enddo
160   continue
        enddo
       enddo
!       write (2058,'(3I15))') wind10u,wind10v,seaicevel
!     & ((i_xists(kounti,kountj),kountj=0,jmax),kounti=0,imax)
      endif

!       write (2058,'(3I15))') wind10u,wind10v,seaicevel
      endif
c end "if(wrt_zero=T)"  ---- end JONO_out -----------------------------
c***********************************************************************

c---cdragi=Cdrag de glace de mer

c---hgbq=epaisseur moyenne
c---hinbq=epaisseur de neige
c---rhog=densite sea ice
c---albq=fraction de glace


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---boucle de calcul du drag de la glace de mer

        dziceb=hgbq(i,j)
        velmd=sqrt((ui-uiceb(l))*(ui-uiceb(l))
     &            +(vi-viceb(l))*(vi-viceb(l)))

        dvar = dziceb*velmd
        dvelx=dvelx+dvar*ui
        dvely=dvely+dvar*vi

c JONO 11-6      dd = 0.5*dticeb*unsmas*cdragi*rhog*1.5*wiceb(l)*(1-albq(i,j))
      dd = 0.5*dticeb*unsmas*cdragi*rhog*Awater*wiceb(l)*(1-albq(i,j))
      dvar = dd*dvar
      dvelx = dd*dvelx
      dvely = dd*dvely

c      write(*,*) 'dvar',dvar

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--nouvelles vitesses de l'iceberg:


      cdiag = 1.+bvar+cvar+dvar
      unsdet = 1./(cdiag*cdiag+f2cdt*f2cdt)
      uu = uiceb(l)*(cdiag-f2cdt*f2cdt) + viceb(l)*f2cdt*(1.+cdiag)
      uu = uu + cdiag*(wvelx+cvelx+dvelx+prxdt+frdtx)
     &        + f2cdt*(wvely+cvely+dvely+prydt+frdty)

      vv = viceb(l)*(cdiag-f2cdt*f2cdt) - uiceb(l)*f2cdt*(1.+cdiag)
      vv = vv + cdiag*(wvely+cvely+dvely+prydt+frdty)
     &        - f2cdt*(wvelx+cvelx+dvely+prxdt+frdtx)
c      write(*,*) 'uu',uu
c      write(*,*) 'vv',vv

      if (flag) then
c-- bilan de tous les termes :
        uunew = uu*unsdet
        vvnew = vv*unsdet

        write(6,'(A,I6,A,2F8.3,A,I8)')
     &  ' dyn: Bilan  ------------+ Iceb. no=', l,
     &  ' ; X,Y=', xn(l),yn(l), ' ; iter=', n
        write(6,'(2A)') '  Iceb n |  n+1  | Ocean | Coriolis',
     &                 '| -grad.P |  W.drag |  A.drag |  Wave '
        write(6,'(A,3F8.4,5F10.6)') 'u',
     &    uiceb(l), uunew, uocean, f2cdt*(viceb(l)+vvnew),
     &    prxdt, wvelx-bvar*uunew, cvelx-cvar*uunew, frdtx
        write(6,'(A,3F8.4,5F10.6)') 'v',
     &    viceb(l), vvnew, vocean, -f2cdt*(uiceb(l)+uunew),
     &    prydt, wvely-bvar*vvnew, cvely-cvar*vvnew, frdty
        write(6,'(A)') '      -------------------+'
      endif

      olduiceb=uiceb(l)
      oldviceb=viceb(l)

      uiceb(l) = uu*unsdet
     &         + (1.-yy) * ( (1.-xx)*urep(i,j)+xx*urep(i+1,j) )
     &         +  yy * ( (1.-xx)*urep(i,j+1)+xx*urep(i+1,j+1) )
      viceb(l) = vv*unsdet
     &         + (1.-yy) * ( (1.-xx)*vrep(i,j)+xx*vrep(i+1,j) )
     &         +  yy * ( (1.-xx)*vrep(i,j+1)+xx*vrep(i+1,j+1) )

c     termexf(l)=f2cdt*(oldviceb+viceb(l))/dticeb
c     termeyf(l)=f2cdt*(olduiceb+uiceb(l))/dticeb

c------norme des differents terme
c      tnormew=((termexw(l)**2)+(termeyw(l)**2))**0.5
c      tnormec=((termexc(l)**2)+(termeyc(l)**2))**0.5
c      tnormef=((termexf(l)**2)+(termeyf(l)**2))**0.5
c      tnormeb=((termexb(l)**2)+(termeyb(l)**2))**0.5
c      tnormep=((termexp(l)**2)+(termeyp(l)**2))**0.5

c      total=tnormew+tnormec+tnormef+tnormeb+tnormep

c      termexw(l)=(tnormew/total)*100
c      termexc(l)=(tnormec/total)*100
c      termexf(l)=(tnormef/total)*100
c      termexb(l)=(tnormeb/total)*100
c      termexp(l)=(tnormep/total)*100

c---comparaison avec bigg
c---/&
c     tnormew=((termexw(l)**2)+(termeyw(l)**2))**0.5
c     tnormec=((termexc(l)**2)+(termeyc(l)**2))**0.5
c     tnormef=((termexf(l)**2)+(termeyf(l)**2))**0.5
c     tnormeb=(((termexb(l)+termexp(l))**2)+
c    &         ((termeyb(l)+termeyp(l))**2))**0.5
c     tnormep=0.

c     total=tnormew+tnormec+tnormef+tnormeb+tnormep

c     termexw(l)=(tnormew/total)*100
c     termexc(l)=(tnormec/total)*100
c     termexf(l)=(tnormef/total)*100
c     termexb(l)=(tnormeb/total)*100
c     termexp(l)=(tnormep/total)*100

c--la norme de la vitesse
c     termeyp(l)=((uiceb(l)**2)+(viceb(l)**2))**.5

c      write(*,'(A6,I1,f10.2)') 'normev ',l,termeyw(l)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      if (flag) write(6,'(A,I5,1P5E11.3)')
     &  ' dyn: New_u,v, uu,vv,unsdet:l=', l, uiceb(l),viceb(l),
     &                                    uu,vv,unsdet
C      write(6,'(A,I8,I4,3F11.6)') 'Iceb 51: ',n
C    &                  ,l,velmd,uiceb(l),viceb(l)

c-fin dyn
c      write (6,*) 'findyn'

      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      end subroutine icebdyn

