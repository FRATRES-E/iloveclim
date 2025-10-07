!dmr -- Ajout du choix optionnel des composantes - Tue Dec 15 16:33:22 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Tue Dec 15 16:33:22 CET 2009

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  modif : 17/02/00
c- modif 24/10/19: cleaned use declaration sections              
! iceberg dynamics is called per iceberg (incluing its i,j,n etc values)

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      subroutine icebdyn(i,j,l,un,vn,xx,yy)
      
!! START_OF_USE_SECTION

      use global_constants_mod, only: dp, ip, vol_mass_dens_air
     >    , vol_mass_dens_ice
      use const_mod, only: gpes, rho0
      use para0_mod, only: imax, jmax, kmax
      use para_mod, only:
      use bloc0_mod, only: kfs, dz, ks2, zw, tmu
     >  , uns2dx, uns2dy, eta
      use bloc_mod, only: smx, q, smy, fs2cor
      use ice_mod, only: hgbq, rhog, albq
      use iceberg_mod, only: it_icb,lmx
     >   , kiceb, wiceb, remsim, hiceb, dticeb
     >   ,wind10_u, wind10_v, uiceb, viceb
     >   ,urep, vrep, cdragw, cdraga, Awater, Aair
     >   ,cdragi, iceberg_info_out_id
      use dynami_mod, only: ug, vg
            
!! END_OF_USE_SECTION
     

      implicit none
  
      integer(kind=ip):: l,ksubm,i,j,k,ii,jj

      real(kind=dp), dimension(kmax) :: un, vn

      real(kind=dp) :: dziceb,v_icb_max,v_fastest,velocity_icb,
     >   unsmas,prclnx,prclny,xx,yy,pressx,pressy,prxdt, prydt,ff2cor,
     >   f2cdt,uvwd,ccwave,frdtx,frdty,uvemd,cvar,
     >   cvelx,cvely,bvar, wvelx,wvely,velmd,bk,bb,
     >   dvar,dvelx,dvely,ui,vi,dd,cdiag,unsdet, uu,vv
         
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|      
! Initialize some variables
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      bvar=0.
      wvelx=0.
      wvely=0.

      dvar=0.
      dvelx=0.
      dvely=0.
      
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Calcul des forces:
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      ksubm=max(kiceb(l),kfs(i,j)) ! submurged part of the iceberg (# k-levels)

c- masse iceberg: SIDE=1.5*FRONT...
      unsmas = 1./(vol_mass_dens_ice*1.5*wiceb(l)*wiceb(l)*(remsim+1.)*hiceb(l))


c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Pression barocline :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      prclnx = 0.
      prclny = 0.

       
        
      do k=ksubm,ks2
        dziceb=dz(k)
        if (k.eq.kiceb(l)) dziceb=hiceb(l)+zw(k+1)

        ii=i
        jj=j
        prclnx = prclnx + dziceb*(1.-xx)*(1.-yy)*
     &     tmu(ii,jj,k)*smx(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               +(q(ii-1,jj,k)-q(ii,jj-1,k)))
        prclny = prclny + dziceb*(1.-xx)*(1.-yy)*
     &     tmu(ii,jj,k)*smy(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               -(q(ii-1,jj,k)-q(ii,jj-1,k)))

        ii=i+1
        jj=j
        prclnx = prclnx + dziceb*xx*(1.-yy)*
     &     tmu(ii,jj,k)*smx(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               +(q(ii-1,jj,k)-q(ii,jj-1,k)))
        prclny = prclny + dziceb*xx*(1.-yy)*
     &     tmu(ii,jj,k)*smy(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               -(q(ii-1,jj,k)-q(ii,jj-1,k)))

        ii=i
        jj=j+1
        prclnx = prclnx + dziceb*(1.-xx)*yy*
     &     tmu(ii,jj,k)*smx(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               +(q(ii-1,jj,k)-q(ii,jj-1,k)))
        prclny = prclny + dziceb*(1.-xx)*yy*
     &     tmu(ii,jj,k)*smy(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               -(q(ii-1,jj,k)-q(ii,jj-1,k)))

        ii=i+1
        jj=j+1
        prclnx = prclnx + dziceb*xx*yy*
     &     tmu(ii,jj,k)*smx(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               +(q(ii-1,jj,k)-q(ii,jj-1,k)))
        prclny = prclny + dziceb*xx*yy*
     &     tmu(ii,jj,k)*smy(ii,jj,3)*((q(ii-1,jj-1,k)-q(ii,jj,k))
     &                               -(q(ii-1,jj,k)-q(ii,jj-1,k)))
      end do

      prclnx = prclnx*uns2dx*dticeb/hiceb(l)
      prclny = prclny*uns2dy*dticeb/hiceb(l)

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Pression barotrope :
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      pressx=0.
      pressy=0.

      ii=i
      jj=j
      pressx=pressx+smx(ii,jj,3)*tmu(ii,jj,ks2)*(1.-xx)*(1.-yy)*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         + (eta(ii-1,jj)-eta(ii,jj-1)) )
      pressy=pressy+smy(ii,jj,3)*tmu(ii,jj,ks2)*(1.-xx)*(1.-yy)*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         - (eta(ii-1,jj)-eta(ii,jj-1)) )

      ii=i+1
      jj=j
      pressx=pressx+smx(ii,jj,3)*tmu(ii,jj,ks2)*xx*(1.-yy)*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         + (eta(ii-1,jj)-eta(ii,jj-1)) )
      pressy=pressy+smy(ii,jj,3)*tmu(ii,jj,ks2)*xx*(1.-yy)*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         - (eta(ii-1,jj)-eta(ii,jj-1)) )

      ii=i
      jj=j+1
      pressx=pressx+smx(ii,jj,3)*tmu(ii,jj,ks2)*(1.-xx)*yy*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         + (eta(ii-1,jj)-eta(ii,jj-1)) )
      pressy=pressy+smy(ii,jj,3)*tmu(ii,jj,ks2)*(1.-xx)*yy*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         - (eta(ii-1,jj)-eta(ii,jj-1)) )

      ii=i+1
      jj=j+1
      pressx=pressx+smx(ii,jj,3)*tmu(ii,jj,ks2)*xx*yy*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         + (eta(ii-1,jj)-eta(ii,jj-1)) )
      pressy=pressy+smy(ii,jj,3)*tmu(ii,jj,ks2)*xx*yy*
     &         ( (eta(ii-1,jj-1)-eta(ii,jj))
     &         - (eta(ii-1,jj)-eta(ii,jj-1)) )

      prxdt = pressx*gpes*uns2dx*dticeb
      prydt = pressy*gpes*uns2dy*dticeb

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Coriolis:
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      ff2cor = ( (1.-yy)*( (1.-xx)*tmu(i,j,ks2)*fs2cor(i,j)
     &                   +  xx*tmu(i+1,j,ks2)*fs2cor(i+1,j) )
     &         +  yy*( (1.-xx)*tmu(i,j+1,ks2)*fs2cor(i,j+1)
     &               +  xx*tmu(i+1,j+1,ks2)*fs2cor(i+1,j+1) ) )
      f2cdt = ff2cor*dticeb

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Air drag
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c JONO dynamic equation: (airdrag)
c        Fa=1/2 Rhoa Ca Aa |va-vi| (va-vi)
c, with Ca=5.2 {=1.3 (Smith93) *4 "we'r doubling the windstrength"}
c, Aa is the "section transversale a la direction de la tension exercee
c par la vent"  (CHECK: why 1/4 in waverad and 1/2 in other drag?)
c BIGGs/Smith:
c "cross-sectional area (above the water line) in a vertical plane normal to the wind"

          
      uvemd=sqrt((wind10_u(i,j)-uiceb(l))*(wind10_u(i,j)-uiceb(l))
     &            +(wind10_v(i,j)-viceb(l))*(wind10_v(i,j)-viceb(l)))

      cvar = 0.5*dticeb*unsmas*cdraga*vol_mass_dens_air*Aair*wiceb(l)*
     &                         remsim*hiceb(l)*uvemd

      cvelx=cvar*( wind10_u(i,j)-uiceb(l) )
      cvely=cvar*( wind10_v(i,j)-viceb(l) )

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---wave radiation
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c for waves absorbed at a deep vertical wall ":
c wave radiation (3.15):
c        Fr = 1/4 * RHOw g a^2 L (Va/|Va|)
c with Va wind speed (vector)
c a = 01010125*|Va|^2 wave amplitude
c L length of the face normal to incoming waves 

      uvwd = wind10_u(i,j)*wind10_u(i,j) + wind10_v(i,j)*wind10_v(i,j)

      ccwave = 0.25*dticeb*unsmas*gpes*Aair*wiceb(l)*
     &        (0.5*0.02025d0)**2*uvwd**(3/2)

      frdtx = ccwave*wind10_u(i,j)
      frdty = ccwave*wind10_v(i,j)

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- water drag:
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      do k=ksubm,ks2
        dziceb=dz(k)
        if (k.eq.kiceb(l)) dziceb=hiceb(l)+zw(k+1)
        velmd=sqrt((un(k)-uiceb(l))*(un(k)-uiceb(l))
     &            +(vn(k)-viceb(l))*(vn(k)-viceb(l)))
        bk = dziceb*velmd
        bvar=bvar+bk
        wvelx=wvelx+bk*un(k)
        wvely=wvely+bk*vn(k)

      end do

      bb = 0.5*dticeb*unsmas*cdragw*rho0*Awater*wiceb(l)
      bvar = bb*bvar
      wvelx = bb*wvelx
      wvely = bb*wvely

c---+---1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- sea ice drag:
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c---ui=vitesse de la glace
c---vi=vitesse de la glace
c---cdragi=Cdrag de glace de mer
c---hgbq=epaisseur moyenne
c---rhog=densite sea ice
c---albq=fraction de glace

      ui= ( (1.-yy)*((1.-xx)*ug(i,j)+xx*ug(i+1,j))
     &         +  yy*((1.-xx)*ug(i,j+1)+xx*ug(i+1,j+1)) )
      vi= ( (1.-yy)*((1.-xx)*vg(i,j)+xx*vg(i+1,j))
     &         +  yy*((1.-xx)*vg(i,j+1)+xx*vg(i+1,j+1)) )

      dziceb=hgbq(i,j)
      velmd=sqrt((ui-uiceb(l))*(ui-uiceb(l))
     &            +(vi-viceb(l))*(vi-viceb(l)))

      dvar = dziceb*velmd
      dvelx=dvelx+dvar*ui
      dvely=dvely+dvar*vi

      dd = 0.5*dticeb*unsmas*cdragi*rhog*Awater*wiceb(l)*(1-albq(i,j))
      dvar = dd*dvar
      dvelx = dd*dvelx
      dvely = dd*dvely

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--New iceberg velocity
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      cdiag = 1.+bvar+cvar+dvar
      unsdet = 1./(cdiag*cdiag+f2cdt*f2cdt)

      uu = uiceb(l)*(cdiag-f2cdt*f2cdt) + viceb(l)*f2cdt*(1.+cdiag)
      uu = uu + cdiag*(wvelx+cvelx+dvelx+prxdt+prclnx+frdtx)
     &        + f2cdt*(wvely+cvely+dvely+prydt+prclny+frdty)

      vv = viceb(l)*(cdiag-f2cdt*f2cdt) - uiceb(l)*f2cdt*(1.+cdiag)
      vv = vv + cdiag*(wvely+cvely+dvely+prxdt+prclnx+frdty)
     &        - f2cdt*(wvelx+cvelx+dvely+prydt+prclny+frdtx)
      

      uiceb(l) = uu*unsdet
     &         + (1.-yy) * ( (1.-xx)*urep(i,j)+xx*urep(i+1,j) )
     &         +  yy * ( (1.-xx)*urep(i,j+1)+xx*urep(i+1,j+1) )
      viceb(l) = vv*unsdet
     &         + (1.-yy) * ( (1.-xx)*vrep(i,j)+xx*vrep(i+1,j) )
     &         +  yy * ( (1.-xx)*vrep(i,j+1)+xx*vrep(i+1,j+1) )

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- this is a warning for high speed bergs {m/s}
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        v_icb_max = 5.
        if (kiceb(l).gt.-1) then
          velocity_icb=sqrt(uiceb(l)*uiceb(l)+viceb(l)*viceb(l))
        endif

        if (velocity_icb.ge.v_icb_max) then
          write(iceberg_info_out_id,*) 'Speeding iceberg',
     &     l,velocity_icb,'{m/s}?','indexes i,j: ',i,j
        endif

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      end subroutine icebdyn

