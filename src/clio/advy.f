!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:42 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:42 CET 2009

!
      SUBROUTINE advy(dt,vt,crh,sm,s0,sx,sxx,sy,syy,sxy)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!
!  PRATHER S (1986) advection scheme (periodic conditions).
!

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip
      
      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use dynami_mod
!! END_OF_USE_SECTION

!

      real(kind=dblp), dimension(imax,jmax) :: vt, sm, s0, sx, sxx, sy
     &               , syy, sxy, f0, fx, fxx, fy, fyy, fxy, fm
     
      real(kind=dblp)::  alg, alg1, alg1q, bet, dt, crh, slpmax, s1max
     >               , s1new, s2new, zin0, alf, alfq, alf1, alf1q, temp
     
      integer(kind=ip):: i, j, jj
      

c~       dimension vt(imax,jmax)
c~       dimension sm(imax,jmax)
c~       dimension s0(imax,jmax),sx(imax,jmax),sxx(imax,jmax),
c~      &          sy(imax,jmax),syy(imax,jmax),sxy(imax,jmax)
c~       dimension f0(imax,jmax),fx(imax,jmax),fxx(imax,jmax),
c~      &          fy(imax,jmax),fyy(imax,jmax),fxy(imax,jmax),
c~      &                fm(imax,jmax)
!
!  bet: Coefficient with values 0/1 to avoid over-writing
!        of arrays.
!
      dimension alg(imax,jmax),alg1(imax,jmax),alg1q(imax,jmax),
     &                bet(imax,jmax)
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  1) Limitation of moments.                                           |
!-----------------------------------------------------------------------
!
      do 5 j=1,jmax
         do 5 i=1,imax
         bet(i,j)=0.0
          alg(i,j)=0.0
          alg1(i,j)=1.0-alg(i,j)
          alg1q(i,j)=alg1(i,j)*alg1(i,j)
         f0(i,j)=0.0
         fx(i,j)=0.0
         fxx(i,j)=0.0
         fy(i,j)=0.0
         fyy(i,j)=0.0
         fxy(i,j)=0.0
         fm(i,j)=0.0
 5    continue

      do 15 j=1,jmax
        do 10 i=1,imax
          slpmax   = max(zero,s0(i,j))
          s1max    = 1.5*slpmax
          s1new    = min(s1max,max(-s1max,sy(i,j)))
          s2new    = min((2.0*slpmax-0.3334*abs(s1new)),
     &               max(abs(s1new)-slpmax,syy(i,j)))
!         s0(i,j)  = slpmax
          sy(i,j)  = s1new
          syy(i,j) = s2new
          sxy(i,j) = min(slpmax,max(-slpmax,sxy(i,j)))
10      continue
15    continue
!
!--1.1 CASE OF EMPTY GRIDS.
!--------------------------
!
      do 29 j=1,jmax
        do 27 i=1,imax
          slpmax   = max(zero,s0(i,j))
          zin0     = 1.0-max(zero,sign(one,-slpmax))
          sx(i,j)  = sx(i,j)*zin0
          sxx(i,j) = sxx(i,j)*zin0
          sy(i,j)  = sy(i,j)*zin0
          syy(i,j) = syy(i,j)*zin0
          sxy(i,j) = sxy(i,j)*zin0
27      continue
29    continue
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  2) Apply mask.                                                      |
!-----------------------------------------------------------------------
!
      do 40 j=1,jmax
        do 30 i=1,imax
          sx(i,j)  = sx(i,j)*zindfa(i,j)
          sxx(i,j) = sxx(i,j)*zindfa(i,j)
          sy(i,j)  = sy(i,j)*zindfa(i,j)
          syy(i,j) = syy(i,j)*zindfa(i,j)
          sxy(i,j) = sxy(i,j)*zindfa(i,j)
30      continue
40    continue
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  3) Initialize volumes of boxes.                                     |
!-----------------------------------------------------------------------
!
      do j=1,jmax
        do i=1,imax
            sm(i,j) = crh*area(i,j)+(1.0-crh)*sm(i,j)
        enddo
      enddo
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  4) Calculate fluxes and moments between boxes j<-->j+1              |
!-----------------------------------------------------------------------
!
      do 80 j=js1,js2
        do 70 i=is1(j),is2(j)
!
!--4.1. Flux from j to j+1 if v.gt.0.0.
!--------------------------------------
!     
          bet(i,j) = max(zero,sign(one,vt(i,j)))
          alf      = (max(zero,vt(i,j))*dt*dxs1(i,j))/sm(i,j)
          alfq     = alf*alf
          alf1     = 1.0-alf
          alf1q    = alf1*alf1
          fm(i,j)  = alf*sm(i,j)
          f0(i,j)  = alf*(s0(i,j)+alf1*(sy(i,j)+(alf1-alf)*syy(i,j)))
          fy(i,j)  = alfq*(sy(i,j)+3.0*alf1*syy(i,j))
          fyy(i,j) = alf*alfq*syy(i,j)
          fx(i,j)  = alf*(sx(i,j)+alf1*sxy(i,j))
          fxy(i,j) = alfq*sxy(i,j)
          fxx(i,j) = alf*sxx(i,j)
!
!  4.1.1. Readjustment of moments remaining in the box.
!
          sm(i,j)  = sm(i,j)-fm(i,j)
          s0(i,j)  = s0(i,j)-f0(i,j)
          sy(i,j)  = alf1q*(sy(i,j)-3.0*alf*syy(i,j))
          syy(i,j) = alf1*alf1q*syy(i,j)
          sx(i,j)  = sx(i,j)-fx(i,j)
          sxx(i,j) = sxx(i,j)-fxx(i,j)
          sxy(i,j) = alf1q*sxy(i,j)
70        continue
80    continue
!
      do 100 j=js1,js2
        jj = j+1
        do 90 i=is1(j),is2(j)
!
!--4.2. Flux from j+1 to j if v.lt.0.0.
!--------------------------------------
!
          alf        = (max(zero,-vt(i,j))*dt*dxs1(i,j))/sm(i,jj)
          alg(i,j)   = alf
          alfq       = alf*alf
          alf1       = 1.0-alf
          alg1(i,j)  = alf1
          alf1q      = alf1*alf1
          alg1q(i,j) = alf1q
          fm(i,j)    = fm(i,j)+alf*sm(i,jj)
          f0(i,j)    = f0(i,j)+
     &                   alf*(s0(i,jj)-alf1*(sy(i,jj)-(alf1-alf)*syy(i,jj)))
          fy(i,j)    = fy(i,j)+alfq*(sy(i,jj)-3.0*alf1*syy(i,jj))
          fyy(i,j)   = fyy(i,j)+alf*alfq*syy(i,jj)
          fx(i,j)    = fx(i,j)+alf*(sx(i,jj)-alf1*sxy(i,jj))
          fxy(i,j)   = fxy(i,j)+alfq*sxy(i,jj)
          fxx(i,j)   = fxx(i,j)+alf*sxx(i,jj)
90        continue
100   continue
      do 120 j=js1,js2
        jj = j-1
        do 110 i=is1(j),is2(j)
!
!  4.2.1. Readjustment of moments remaining in the box.
!
          sm(i,j)  = bet(i,jj)*sm(i,j)
     &               +(1.0-bet(i,jj))*(sm(i,j)-fm(i,jj))
          s0(i,j)  = bet(i,jj)*s0(i,j)
     &               +(1.0-bet(i,jj))*(s0(i,j)-f0(i,jj))
          sy(i,j)  = alg1q(i,jj)*(sy(i,j)
     &               +3.0*alg(i,jj)*syy(i,j))
          syy(i,j) = alg1(i,jj)*alg1q(i,jj)*syy(i,j)
          sx(i,j)  = bet(i,jj)*sx(i,j)
     &               +(1.0-bet(i,jj))*(sx(i,j)-fx(i,jj))
          sxx(i,j) = bet(i,jj)*sxx(i,j)
     &               +(1.0-bet(i,jj))*(sxx(i,j)-fxx(i,jj))
          sxy(i,j) = alg1q(i,jj)*sxy(i,j)
110     continue
120   continue
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  5) Put the temporary moments into appropriate neighboring boxes.    |
!-----------------------------------------------------------------------
!
      do 140 j=js1,js2
        jj = j-1
        do 130 i=is1(j),is2(j)
!
!
!--5.1. Flux from j to j+1 if v.gt.0.0.
!--------------------------------------
!
          sm(i,j)  = bet(i,jj)*(sm(i,j)+fm(i,jj))
     &                +(1.0-bet(i,jj))*sm(i,j)
          alf      = bet(i,jj)*fm(i,jj)/sm(i,j)
          alf1     = 1.0-alf
          temp     = alf*s0(i,j)-alf1*f0(i,jj)
          s0(i,j)  = bet(i,jj)*(s0(i,j)+f0(i,jj))
     &               +(1.0-bet(i,jj))*s0(i,j)
          sy(i,j)  = bet(i,jj)*(alf*fy(i,jj)+alf1*sy(i,j)+3.0*temp)
     &                     +(1.0-bet(i,jj))*sy(i,j)
          syy(i,j) = bet(i,jj)*(alf*alf*fyy(i,jj)+alf1*alf1*syy(i,j)
     &               +5.0*(alf*alf1*(sy(i,j)-fy(i,jj))
     &                   -(alf1-alf)*temp))+
     &               (1.0-bet(i,jj))*syy(i,j)
          sxy(i,j) = bet(i,jj)*(alf*fxy(i,jj)+alf1*sxy(i,j)+
     &               3.0*(-alf1*fx(i,jj)+alf*sx(i,j)))+
     &                     (1.0-bet(i,jj))*sxy(i,j)
          sx(i,j)  = bet(i,jj)*(sx(i,j)+fx(i,jj))
     &               +(1.0-bet(i,jj))*sx(i,j)
          sxx(i,j) = bet(i,jj)*(sxx(i,j)+fxx(i,jj))
     &               +(1.0-bet(i,jj))*sxx(i,j)
130        continue
140   continue
!
      do 160 j=js1,js2
        do 150 i=is1(j),is2(j)
!
!--5.2. Flux from j+1 to j if v.lt.0.0.
!--------------------------------------
!
          sm(i,j)  = bet(i,j)*sm(i,j)
     &               +(1.0-bet(i,j))*(sm(i,j)+fm(i,j))
          alf      = (1.0-bet(i,j))*fm(i,j)/sm(i,j)
          alf1     = 1.0-alf
          temp     = -alf*s0(i,j)+alf1*f0(i,j)
          s0(i,j)  = bet(i,j)*s0(i,j)+
     &                    (1.0-bet(i,j))*(s0(i,j)+f0(i,j))
          sy(i,j)  = bet(i,j)*sy(i,j)+
     &               (1.0-bet(i,j))*(alf*fy(i,j)+alf1*sy(i,j)+3.0*temp)
          syy(i,j) = bet(i,j)*syy(i,j)+
     &               (1.0-bet(i,j))*(alf*alf*fyy(i,j)+alf1*alf1*syy(i,j)+
     &          5.0*(alf*alf1*(-sy(i,j)+fy(i,j))+(alf1-alf)*temp))
          sxy(i,j) = bet(i,j)*sxy(i,j)+
     &                   (1.0-bet(i,j))*(alf*fxy(i,j)+alf1*sxy(i,j)+
     &              3.0*(alf1*fx(i,j)-alf*sx(i,j)))
          sx(i,j)  = bet(i,j)*sx(i,j)+
     &               (1.0-bet(i,j))*(sx(i,j)+fx(i,j))
          sxx(i,j) = bet(i,j)*sxx(i,j)
     &               +(1.0-bet(i,j))*(sxx(i,j)+fxx(i,j))
150     continue
160   continue
!
!-- raccord cyclique
!
        call raccord(sm(1,1),zero,1,8)
        call raccord(s0(1,1),zero,1,8)
        call raccord(sx(1,1),zero,1,8)
        call raccord(sxx(1,1),zero,1,8)
        call raccord(sy(1,1),zero,1,8)
        call raccord(syy(1,1),zero,1,8)
        call raccord(sxy(1,1),zero,1,8)
      return
      end
