!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:46 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:46 CET 2009

      SUBROUTINE forcat(iyear,xjour)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!-This routine calculates the heat fluxes from atmospheric data
! read by forcng.
! Derniere modification 04/06/95
! Differences entre forcat.f et formel.f tcn(i,j,1) - scal(i,j,ks2,1)
! et le runoff, oceacn.com

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip      

#if ( ISOOCN >= 1 ) 
      USE iso_param_mod, ONLY : ieau
#endif

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use ice_mod
      
      use newunit_clio_mod, only: clio3_out_id      
      
!! END_OF_USE_SECTION

      integer(kind=ip) :: iyear
      real(kind=dblp)  :: xjour

      real(kind=dblp), dimension(96) :: ws, zmue, zalcnp
      
      integer(kind=ip) :: jj, j, i, indaet, ijour, ih, k, ii1, ii2, jj1, jj2
      

      real(kind=dblp)  :: zind1, zind2, zind3, dday, dcor, dec, pdecli
     >                  , sdec, cdec, slat, zps, zpc, zljour, dws
     >                  , zlmidi, zalcnq, zmudum, zalb, zalcn, zalbp
     >                  , zaldum, es, evg, evo, e, albg, frsdtg, frsdto
     >                  , albo, frsdrg, frsdfg, frsdro, frsdfo, zinda
     >                  , rhoa, zqsat, dteta, deltaq, ztvmoy, obouks
     >                  , psims, psihs, psils, obouki, xins, psimi
     >                  , psihi, psili, stab, psim, psih, psil, zzero
     >                  , cmn, chn, cln, cmcmn, ch, ce, drghce, drghch


c~       dimension ws(96),zmue(96),zalcnp(96)

!     write(clio3_out_id,*) 'debut de forcat.f'
      integer(kind=ip) :: iflgaa, iflgab, iflago, nintsr
      
      data iflgaa /0/
      data iflgab /0/
      data iflago /0/
      data nintsr /24/

      real(kind=dblp) :: zeps, zeps0, zeps1, dpi
!
      zeps  = 1d-20
      zeps0 = 1.d-13
      zeps1 = 1.d-06
      dpi   = 2.*pi
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  1) READ ATMOSPHERIC FORCING.                                        |
!-----------------------------------------------------------------------
!
          call forcng(iyear,xjour)
!
            do jj=jcl1,jcl2
               tenagx(ims1-1,jj) = tenagx(ims2,jj)    
               tenagx(ims2+1,jj) = tenagx(ims1,jj)
               tenagy(ims1-1,jj) = tenagy(ims2,jj)
               tenagy(ims2+1,jj) = tenagy(ims1,jj)
            enddo
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  2) COMPUTATION OF SNOW PRECIPITATION                                |
!-----------------------------------------------------------------------
!
!--2.1. FRACTION OF SNOW PRECIPITATION (LEDLEY, 1985).
!-----------------------------------------------------
!
          do j=js1,js2
            do i=is1(j),is2(j)
              zind1       = max(zero,sign(one,253.0-tabq(i,j)))
              zind2       = max(zero,sign(one,272.0-tabq(i,j)))
              zind3       = max(zero,sign(one,281.0-tabq(i,j)))
#if ( ISOOCN >= 1 )
              hnplbq(i,j,ieau) = (zind1+(1.0-zind1)*
#else
              hnplbq(i,j) = (zind1+(1.0-zind1)*
#endif
     &                      (zind2*(0.5+(272.0-tabq(i,j))/38.0)+
     &                            (1.0-zind2)
     &                       *zind3*((281.0-tabq(i,j))/18.0)))*
#if ( ISOOCN >= 1 )
     &                      fwat(i,j,ieau)
#else
     &                      fwat(i,j)
#endif


#if ( ISOOCN >= 1 )
               WRITE(*,*) "forcat : héhé ", fwat(i,j,ieau), ieau
               fwat(i,j,ieau)=fwat(i,j,ieau)+runoff(i,j)*86400*1000
#else
               WRITE(*,*) "forcat : héhé ", fwat(i,j), runoff(i,j)
               fwat(i,j)=fwat(i,j)+runoff(i,j)*86400*1000
#endif
            enddo
          enddo
!
!--2.1. FRACTION OF SNOW PRECIPITATION (ROSS & WALSH, 1987).
!-----------------------------------------------------------
!
!          do j=js1,js2
!            do i=is1(j),is2(j)
!              if (tabq(i,j).lt.268.15) then
!                hnplbq(i,j) = fwat(i,j)
!              else 
!                if (tabq(i,j).lt.276.15) then
!                  hnplbq(i,j) = ((276.15-tabq(i,j))/8.0)*fwat(i,j)
!                else 
!                  hnplbq(i,j) = 0.0
!                endif
!              endif
!            enddo
!          enddo
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  3) COMPUTATION OF SOLAR IRRADIANCE.                                 |
!-----------------------------------------------------------------------
!
          indaet = 1
          ijour  = int(xjour)
          dday = xjour*2.*pi/yeaday
          dcor=1.000d0+0.0013d0*sin(dday)+0.0342d0*cos(dday)

!          write(clio3_out_id,*) 'dcor: ',dcor,' jour : ', ijour
          dec  = pdecli(indaet,ijour)
          dec  = dec*radian
          sdec = sin(dec)
          cdec = cos(dec)
          do j=js1,js2
            do i=is1(j),is2(j)
            slat   = covrai(i,j)
            zps    = slat*sdec
            zpc    = cos(asin(slat))*cdec
            zljour = acos(-sign(one,zps)
     &               *min(one,sign(one,zps)*(zps/zpc)))
            dws    = (2.0*zljour)/real(nintsr)
            zlmidi = asin(( zps +zpc ))/radian
            zalcnq = 0.0
            do k=1,nintsr
              ws(k)     = zljour-(real(k)-0.5)*dws
              zmue(k)   = max(zero,zps+zpc*cos(ws(k)))
              zalcnp(k) = 0.05/(1.1*zmue(k)**1.4+0.15)
              zalcnq    = zalcnq+zalcnp(k)*dws
            enddo
            zalcnq = zalcnq/max(2.0*zljour,zeps0)
            zmudum = 0.4
!
!
              ih=max0(0,isign(1,j-jeq))
!
!--ALBEDO.
!
                    call shine(ih,zmudum,tfsn,tfsg,ts(i,j),hgbq(i,j)
     &                  ,hnbq(i,j),zalb,zalcn,zalbp,zaldum)
!
!--SHORTWAVE RADIATION.
!
!             evg=611.0*10.0**(9.5*(tdew(i,j)-273.16)
!    &            /(tdew(i,j)-7.660))
!             evo=611.0*10.0**(7.5*(tdew(i,j)-273.16)
!    &            /(tdew(i,j)-35.86))
!             es1=611.0*10.0**(7.5*(tabq(i,j)-273.16)
!    &            /(tabq(i,j)-7.660))
!             es2=611.0*10.0**(7.5*(tabq(i,j)-273.16)
!    &           /(tabq(i,j)-35.86))
!             evg=tdew(i,j)*es1
!             evo=tdew(i,j)*es2
              es =611.0*exp(min(sign(17.269*one,
     &            tabq(i,j)-too),sign(21.875*one,
     &            tabq(i,j)-too))*abs(tabq(i,j)-too)/
     &            (tabq(i,j)-35.86+max(zero,sign(28.2*one,
     &            too-tabq(i,j)))))
              evg=tdew(i,j)*es
              evo=tdew(i,j)*es
              e  =tdew(i,j)*es
              qabq(i,j) = (0.622*e)/(psbq(i,j)-(1.0-0.622)*e)
!
!--ZILLMAN.
!
               albg =  (1.0-cloud(i,j))*zalbp+cloud(i,j)*zalb
              frsdtg = 0.0
              frsdto = 0.0
              do k=1,nintsr
                albo   = (1.0-cloud(i,j))*zalcnp(k)+cloud(i,j)*zalcn
!               frsdtg = frsdtg+dws*(1.0-albg)*
!    &                   (1368.0*zmue(k)*zmue(k))/
!    &                         ((zmue(k)+2.7)*evg*1.0e-05+1.085*zmue(k)+0.10)
                 frsdto = frsdto+dws*(1.0-albo)*
     &                   (1368.0*zmue(k)*zmue(k))/
     &                         ((zmue(k)+2.7)*evo*1.0e-05+1.085*zmue(k)+0.10)
               enddo
!
!--Computation of the solar heat flux absorbed at
!  the snow/ice and ocean surfaces.
!
!             fsolg(i,j) =(1.0-0.6*cloud(i,j)
!    &                     *cloud(i,j)*cloud(i,j))*frsdtg/dpi
!             fsolcn(i,j)=(1.0-0.6*cloud(i,j)
!    &                    *cloud(i,j)*cloud(i,j))*frsdto/dpi
!             fsolg(i,j) = 0.9*min(one,(1-.62*cloud(i,j)+.0019*zlmidi))*
!    &                     frsdtg/dpi
              fsolcn(i,j)= 0.9*min(one,(1-.62*cloud(i,j)+.0019*zlmidi))*
     &                     frsdto/dpi

!--Computation of the effective sea-water albedo.
!
!             albege(i,j) = 1.0-fsolg(i,j)/
!    &                      max(fsolg(i,j)/(1.0-albg),zeps0)
              albecn(i,j) = 1.0-fsolcn(i,j)/
     &                      max(fsolcn(i,j)/(1.0-((1.0-cloud(i,j))
     &                      *zalcnq+cloud(i,j)*zalcn)),zeps0)
!
!--SHINE AND CRANE.
!
              frsdrg = 0.0
              frsdfg = 0.0
              frsdro = 0.0
              frsdfo = 0.0
              do k=1,nintsr
                 frsdrg = frsdrg+dws*
     &                         (1368.0*zmue(k)*zmue(k)*(1.0-zalbp))/
     &                         (1.2*zmue(k)+(1.0+zmue(k))*evg*1.0e-05+0.0455)
                frsdfg  = frsdfg+dws*
     &                         ((53.5+1274.5*zmue(k))*sqrt(zmue(k))
     &                   *(1.0-0.996*zalb))/
     &                   (1.0+0.139*(1.0-0.9435*zalb)*tauc(i,j))
!               frsdro = frsdro+dws*
!    &                         (1368.0*zmue(k)*zmue(k)*(1.0-zalcnp(k)))/
!    &                         (1.2*zmue(k)+(1.0+zmue(k))*evo*1.0e-05+0.0455)
!               frsdfo = frsdfo+dws*
!    &                         ((53.5+1274.5*zmue(k))*sqrt(zmue(k))
!    &                   *(1.0-0.996*zalcn))/
!    &                   (1.0+0.139*(1.0-0.9435*zalcn)*tauc(i,j))
              enddo
!
!--Computation of the solar heat flux absorbed at
!  the snow/ice and ocean surfaces.
!
              zinda       = 1.0-max(zero,sign(one,-(-0.5-albq(i,j))))
              fsolg(i,j)  = ((1.0-cloud(i,j))*frsdrg
     &                      +cloud(i,j)*frsdfg)/dpi
              fsolcn(i,j) = (((1.0-cloud(i,j))*frsdro
     &                      +cloud(i,j)*frsdfo)/dpi)*zinda
     &                      +(1.-zinda)*fsolcn(i,j)

!--Computation of the effective snow/ice albedo.
!
              albege(i,j) = (((frsdrg/(1.0-zalbp))*zalbp+
     &                      (frsdfg/(1.0-zalb))*zalb)/
     &                      max(frsdrg/(1.0-zalbp)
     &                      +frsdfg/(1.0-zalb),zeps0))
     &                      *(1-max(zero,sign(one,-hgbq(i,j))))

!--Taking into account the ellipsity of the earth orbit
              fsolg(i,j)  =fsolg(i,j)*dcor
              fsolcn(i,j) =fsolcn(i,j)*dcor
              ii1=85
              ii2=55
              jj1=2
              jj2=2

!
!--Computation of the effective sea-water albedo.
!
!             albecn(i,j) = ((frsdro/(1.0-zalcnq))*zalcnq+
!    &                        (frsdfo/(1.0-zalcn))*zalcn)/
!    &                      max(frsdro/(1.0-zalcnq)
!    &                      +frsdfo/(1.0-zalcn),zeps0)
!
            enddo
          enddo
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  4) PARTIAL COMPUTATION OF HEAT, WATER AND MOMENTUM FLUXES           |
!-----------------------------------------------------------------------
!
          do 170 j=js1,js2
            do 160 i=is1(j),is2(j)
!
! rhoa: Air density, obtained from the equation of
!        state for dry air.
!
              rhoa = psbq(i,j)/(287.04*tabq(i,j))
!
!--4.1. Calculate turbulent heat fluxes over water.
!--------------------------------------------------
!
!NCAR: bulk sensible and latent heat fluxes.
!
!  es: Saturation vapor pressure.
!  e : Vapor pressure.
!
               es         =  611.0*10.0**(7.5*(scal(i,j,ks2,1)
     &                      -273.16)/(scal(i,j,ks2,1)-35.86))
!              es         =  611.0*10.0**(7.5*(tcn(i,j,1)
!    &                      -273.16)/(tcn(i,j,1)-35.86))
!             es2        =  611.0*10.0**(7.5*(tabq(i,j)
!    &                      -273.16)/(tabq(i,j)-35.86))
!              e          =  611.0*10.0**(7.5*(tdew(i,j)
!    &                      -273.16)/(tdew(i,j)-35.86))
!             e          = tdew(i,j)*es2
              zqsat      = (0.622*es)/(psbq(i,j)-(1.0-0.622)*es)
!             qabq(i,j)  = (0.622*e)/(psbq(i,j)-(1.0-0.622)*e)
!

!--computation of drag coefficients following Bunker
!
!              tvirt      = tabq(i,j)*(1+.608*qabq(i,j))
!              tvirtsu    = scal(i,j,ks2,1)*(1+.608*zqsat)
!              ce         = ceb(vabq(i,j),tabq(i,j)-scal(i,j,ks2,1))
!              ceorg      = ceb(vabq(i,j),tvirt-tvirtsu)
!              ce         = .92*ceorg
!              ch         = .87*ceorg

!              ce         = ceb(vabq(i,j),tabq(i,j)-tcn(i,j,1))
!              ce         = 1.75e-03
!             ch         = ce

!--drag coefficients from Large and Pond

!  Stability parameters
              dteta  = scal(i,j,ks2,1)-tabq(i,j)
!             dteta  = tcn(i,j,1)-tabq(i,j)
              deltaq = qabq(i,j)-zqsat
              ztvmoy = tabq(i,j)*(1.+2.2E-3*tabq(i,j)*qabq(i,j))
              Obouks = -70.0*10.*(dteta+3.2E-3*ztvmoy*ztvmoy*deltaq)/
     &                 (vabq(i,j)*vabq(i,j)*ztvmoy)
              Obouks =max(zero,Obouks)
              psims  = -7.0*Obouks
              psihs  =  psims
              psils  =  psims
              Obouki =  -100.0*10.0*(dteta
     &                    +2.2E-3*ztvmoy*ztvmoy*deltaq)/
     &                 (vabq(i,j)*vabq(i,j)*ztvmoy)
              Obouki = min(zero,Obouki)
              xins   = (1-16.*Obouki)**.25
              psimi  = 2*log((1+xins)/2)+
     &                 log((1+xins**2)/2)-2*atan(xins)+pi/2
              psihi  = 2*log((1+xins**2)/2)
              psili  = psihi
              stab   = max(zero,sign(one,dteta))
              psim   = stab*psimi+(1.0-stab)*psims
              psih   = stab*psihi+(1.0-stab)*psihs
              psil   = stab*psili+(1.0-stab)*psils

!  computation of intermediate values
              zzero  = .032*1.5E-3*vabq(i,j)*vabq(i,j)/gpes
              cmn    = (vkarmn/log(10./zzero))**2
              chn    = 0.0327*vkarmn/log(10./zzero)
              cln    = 0.0346*vkarmn/log(10./zzero)
              cmcmn  = 1/(1-sqrt(cmn)*psim/vkarmn)

!  ch and ce
              ch     = chn*cmcmn/(1-chn*psih/(vkarmn*sqrt(cmn)))
              ce     = cln*cmcmn/(1-cln*psil/(vkarmn*sqrt(cmn)))
!
!End compuation of ch and ce
!
              drghce     = rhoa*ce*vabq(i,j)
              drghch     = rhoa*ch*vabq(i,j)
!              fcscn(i,j) = drghch*1004.0*(tcn(i,j,1)-tabq(i,j))
               fcscn(i,j) = drghch*1004.0*(scal(i,j,ks2,1)-tabq(i,j))
              flecn(i,j) = drghce*2.5e+06*(zqsat-qabq(i,j))
!
!
!--4.3. CALCULATE EVAPORATION OVER WATER.
!----------------------------------------
!
              fevabq(i,j) = flecn(i,j)/cevap
!
160         continue
170       continue
!
!
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!- fin de la routine forcat -
      end
