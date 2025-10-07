      subroutine iceb_dist(niceb,i,j,xiceb,yiceb)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--ancien (28/03/00) nom = dcotes(.f)
c-controle de la distance par rapport a la cote
c-pour les iceberg de distribution fixe
!- modif 24/10/19: cleaned use declaration sections              

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip             
      use const_mod, only:
      use para0_mod, only:
      use para_mod, only:
      use bloc0_mod, only: ks2, tms
      use bloc_mod, only: cmx, cmy
      use iceberg_mod, only: xn, yn, d_cotes
      use reper_mod, only: dxwi, xwi1,dywj,ywj1
!! END_OF_USE_SECTION

      implicit none
      
!! START_OF_INCLUDE_SECTION

!PB#include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "iceberg.com"
! [SCRPTCOM] #include "reper.com"

!! END_OF_INCLUDE_SECTION

!PB variables added after imposing implicit none
      integer(kind=ip):: i,j,niceb
      real(kind=dblp) :: distance,xold,yold,xiceb,yiceb

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|     
c--boucle sur les icebergs

c-calcul de i et j de chaque icebergs


c-test de terre ou eau
c---EST----------
c        i=1+nint((xn(l)-xwi1)/dxwi)
c        j=1+nint((yn(l)-ywj1)/dywj)

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c      write(*,*) 'iceb,bordEST:',xiceb,((i-0.5)*dxwi+xwi1),xn(niceb)
       if (tms((i+1),j,ks2).eq.0)  then
c        write (*,*) 'terre bord EST:numiceb=',niceb

c-----on est en terre du cote est
c-------->test de la distance avec la cote

c--distance en degres d'angle
         distance=abs((i-0.5)*dxwi+xwi1-xiceb)
c-distance en m (smx*dwxi=distance d'une maille/nbre de degre par maille)
         distance=distance*cmx(i,j,0)/dxwi
c         write(*,*) 'distance cotes:', distance,d_cotes

c--------si la distance est plus petite qu'un min
c--------alors on regle xn pour que ca le soit
          if (distance.lt.d_cotes) then
c            write (*,*) 'coor EST:distance,d_cotes',distance,d_cotes
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
            xold=xiceb
            xn(niceb)=xn(niceb)-
     &              ((d_cotes-distance)/cmx(i,j,0)*dxwi)

c--nouveau calcul de la distance
            distance=abs((i-0.5)*dxwi+xwi1-xn(niceb))
            distance=distance*cmx(i,j,0)/dxwi

c           write (*,*) 'EST:niceb,xold,xnnew,distancenew,d_cotes',
c    &                  niceb,xold,xn(niceb),distance,d_cotes
          endif
        endif

c---OUEST
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c      write(*,*) 'iceb,bordOUEST:',xiceb,((i-1)*dxwi+xwi1)
c      write(*,*) 'tmsOUEST:',tms((i-1.5),j,ks2)
      if (tms((i-1),j,ks2).eq.0)  then
c        write (*,*) 'terre bord OUEST:numiceb=',niceb
         distance=abs((i-1.5)*dxwi+xwi1-xiceb)
         distance=distance*cmx(i,j,0)/dxwi
c         write(*,*) 'distance cotes:', distance,d_cotes

         if (distance.lt.d_cotes) then
            xold=xiceb
            xn(niceb)=xn(niceb)+
     &              ((d_cotes-distance)/cmx(i,j,0)*dxwi)
            distance=abs((i-1.5)*dxwi+xwi1-xn(niceb))
            distance=distance*cmx(i,j,0)/dxwi
c           write (*,*) 'OUEST:niceb,xold,xnnew,distancenew,d_cotes',
c    &                  niceb,xold,xn(niceb),distance,d_cotes
         endif
      endif

c---NORD
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c      write(*,*) 'iceb,bordNORD:',yiceb,((j-2)*dywj+ywj1)
      if (tms(i,(j+1),ks2).eq.0)  then
c        write (*,*) 'terre bord NORD:numiceb=',niceb
         distance=abs((j-2)*dywj+ywj1-yiceb)
         distance=distance*cmy(i,j,0)/dywj
c         write(*,*) 'distance cotes:', distance,d_cotes

         if (distance.lt.d_cotes) then
            yold=yiceb
            yn(niceb)=yn(niceb)-
     &              ((d_cotes-distance)/cmy(i,j,0)*dywj)
            distance=abs((j-2)*dywj+ywj1-yn(niceb))
            distance=distance*cmy(i,j,0)/dywj
c           write (*,*) 'NORD:niceb,xold,xnnew,distancenew,d_cotes',
c    &                  niceb,yold,yn(niceb),distance,d_cotes
         endif
      endif

c---SUD
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c      write(*,*) 'iceb,bordSUD:',yiceb,((j-1)*dywj+ywj1)
      if (tms(i,(j-1),ks2).eq.0)  then
c        write (*,*) 'terre bord SUD:numiceb=',niceb
         distance=abs((j-1)*dywj+ywj1-yiceb)
         distance=distance*cmy(i,j,0)/dywj
c         write(*,*) 'distance cotes:', distance,d_cotes

         if (distance.lt.d_cotes) then
            yold=yiceb
            yn(niceb)=yn(niceb)+
     &              ((d_cotes-distance)/cmy(i,j,0)*dywj)
            distance=abs((j-1)*dywj+ywj1-yn(niceb))
            distance=distance*cmy(i,j,0)/dywj
c           write (*,*) 'SUD:niceb,xold,xnnew,distancenew,d_cotes',
c    &                  niceb,yold,yn(niceb),distance,d_cotes
         endif
      endif

c      write(*,*) 'fin iceberg',niceb
c      read(*,*)

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---end

      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine iceb_dist -
      end subroutine iceb_dist
