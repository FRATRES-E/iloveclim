!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009

!dmr --- 2015-10-23
!dmr --- Added the computations of freshwater for the Arctic by Frazer Davies
!dmr --- Following flag set to one activates that in the output
!dmr --- FRAZER_ARCTIC 1 -> moved to choixcomposantes
!dmr ---

      SUBROUTINE inforun(nn99)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!    calcule et ecrit  quelques variables globales, sur le fichier "evolu",
!    frequence "ninfo", serie chronologique.
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  modif : 09/08/02


      use comemic_mod,  only:
      use const_mod,    only: cpo, gpes, one, rho0, svrdrp, yeaday, zero

      use para0_mod,    only: imax, jmax, kmax, nsmax
      use para_mod,     only: nbsmax
      use bloc0_mod   , only: dx, dy, iberp, jberp, jeq, js1, js2, ju1
     &  , ju2, ks1, ks2, ku2, nsav, spvr, tpstot, unsdy, fqajc, daeta
     &  , is1, dzw, eta, ub, vb, u, dz, v, z, dts, alphmi, scal, scalr
     &  , fss, b, is2, tmu, iu1, scal0, ahs, tms, iu2, w, unsdx
      use bloc_mod    , only: aire, ctmi, cmxy, zsurfo, cmx, smy,smx,cmy
     &  , zvols, zvolv, zvolw, koutpu, ninfo, nstart, ntmoy, numit
     &  , zsurf, zvolo
      use isoslope_mod, only: viso, uiso
      use ice_mod,      only: ficebergn, ficebergs, toticesm, xlg, vwx
     &  , albq, hgbq, hnbq, ddtb
      use dynami_mod,   only: bound, vg, dxs1, ug

      use reper_mod,    only: dlat, fmtw, ndhsf, nferme, nocean, nvhsf
     &  , nvinfo, ylat1, iehsf, ishsf, jshsf, iszon, titvar, vinfor
     &  , zmdeta, jehsf, iezon, scalwr, ktsum, vinfom

      use global_constants_mod, only: dblp=>dp, ip

      use ipcc_output_mod, only: moc, thex, tmc
      use newunit_clio_mod, only: clio3_out_id, mouchard_id, evolu_id
     &  , testevolu_id      


      implicit none


      real(dblp), dimension(imax,jmax) :: heaticbism
      real(dblp)                       :: vicbismn,vicbisms

      common /icbism/ heaticbism,vicbismn,vicbisms

!--variables locales :
      integer(ip), dimension(nsmax) :: nnvk
      real(dblp), dimension(3)      :: vvk
      real(dblp), dimension(0:1)    :: zitsum
      logical                       :: flgout

!--variables locales a conserver d'un appel a l'autre -> dans "reper .com" .


! --- dmr Local variables following implicit none declaration

      integer(ip) :: i, ii, ii1, ii2, ip1, j, jj, jj1, jj2, k, kk, n
     &             , nhsf, nm1n, nm2n, nn, nninfo, nniter, nnp, ns, nv
     &             , nn99
!PB
     &             , countjj
      real(dblp)  :: aalpha, ccydif, cny, convfx, ctmobs, factke, phiy
     &             , ssc2, sumk, therma, v2cd2, vber, vbord, vfram, vv
     &             , vv2, vvn, vvp, vvpsm, yy, zold, zz
!PB
     &             , sszonm, ss, vvpa, vvpasm, vvpd, vvpdsm
     &             , vvtfa, vvtfd, vvtba, vvtbd, vvtca, vvtcd
     &             , vva, vvd, flux, ccxdif, uu2, u2cd2, phix, cnx
     &             , saltatlantic, volatlantic,saltatlantic_old
     &             , saltatlantic_new, saltatlantic_dt
!PB

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Prepare & debute le remplissage de "vinfor".                    |
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  => acces pendant le run, pour calculer et/ou ecrire Var.evolution
!-----------------------------------------------------------------------

!- flag = TRUE <=> write on "evolu" file at this iter. :
      flgout = mod(numit,ninfo).eq.0 .or. numit.eq.1

!- Initialisation (a zero) de "vinfor":
      do nv=1,nvinfo
        vinfor(nv) = 0.
      enddo 

      vinfor(1) = DFLOAT(numit)
      vinfor(2) = tpstot / ( 86400. * yeaday )
      nv = 2

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Calcul de quelques variables globales :                         |
!-----------------------------------------------------------------------

      if (flgout) then

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--2.1 Variables deja cumulees dans une autre routine :

!-----
!- evaluation de la variation de eta au cours de l it. par W(kfond) :
!  -->  transferre le 07/12/93 dans "flucor"
!-----

!--Moyenne sur tout le bassin de : 1) l energie dissipee par Ajust.Conv.
!                                  2) variation de l elevat. au cours de l it.
      do j=js1,js2
       do i=is1(j),is2(j)
        vinfor(nv+1) = vinfor(nv+1) + ctmi(i,j,ks2,0) * fqajc(i,j,1)
        vinfor(nv+3) = vinfor(nv+3) + ctmi(i,j,ks2,0) * daeta(i,j)
        daeta(i,j) = 0.
       enddo 
      enddo 
      vinfor(nv+1) = vinfor(nv+1) * zsurf
      vinfor(nv+3) = vinfor(nv+3) * zmdeta

!--Integre la Freq.d AjC sur le volume oceanique global :
      do k=ks1+1,ks2
        sumk = 0.
        kk = k - 1
        do j=js1,js2
         do i=is1(j),is2(j)
           sumk = sumk + ctmi(i,j,kk,0) * fqajc(i,j,k)
         enddo 
        enddo 
        vinfor(nv+2) = vinfor(nv+2) + sumk * dzw(k)
      enddo

!- fraction du volume total, en o/oo :
      vinfor(nv+2) = vinfor(nv+2) * zvolw * 1000.
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
      nv = nv + 4

      if (flgout .or. ntmoy.eq.2) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--2.2 Variables non-cumulees dans une autre routine :

!--Moyenne sur tout le bassin de : 3) de l elevation
!--            (ancienne version : 3) de l energ. due a l elevation)
      do j=js1,js2
       do i=is1(j),is2(j)
        vinfor(nv) = vinfor(nv)
     &             + ctmi(i,j,ks2,0) * eta(i,j)
!    &             + ctmi(i,j,ks2,0) * eta(i,j) * eta(i,j)
       enddo 
      enddo 
!     vinfor(nv) = 0.5 * gpes * rho0 * vinfor(nv) * zsurf
      vinfor(nv) = vinfor(nv) * zsurf
!-----
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--2.3 Horizontal Transport (Sverdrup) , entre qq points particuliers :
!- passage de Drake (nhsf=1), Indonesie (nhsf=2)
      do nhsf=1,nvhsf
        nn = nv + nhsf
        if (iehsf(nhsf).eq.0) then
!- Integre (-ub).(-dy) :
          ii = ishsf(nhsf)
          do j=jshsf(nhsf),jehsf(nhsf)
            vinfor(nn) = vinfor(nn) + cmy(ii,j,3) * ub(ii,j)
          enddo 
          vinfor(nn) = vinfor(nn) * dy * svrdrp
        else
!- Integre vb.dx :
          jj = jshsf(nhsf)
          do i=ishsf(nhsf),iehsf(nhsf)
            vinfor(nn) = vinfor(nn) + cmx(i,jj,3) * vb(i,jj)
          enddo 
          vinfor(nn) = vinfor(nn) * dx * svrdrp
        endif
      enddo 
      nv = nv + nvhsf
!-- Florida Strait (i=177, j=79,83)
!L30  ii=177
!L30  jj1=79
!L30  jj2=83
!L30  nv=nv+1
!L30  do j=jj1,jj2
!L30    vinfor(nv) = vinfor(nv) + cmy(ii,j,3) * ub(ii,j)
!L30  enddo
!L30  vinfor(nv) = vinfor(nv) * dy * svrdrp

!--2.4 Horizontal Transport (Sverdrup) , entre qq points particuliers :
!- calcul separe Flux Sud / Flux Nord ou Flux Ouest / Flux Est
!-----(nhsf definitions)
!  4 = Dan.Strait   , 5 = Icel.-Scotl , 6 = Fram Strait , 7 = N.W.Passage
!  8 = Spitz-Norway , 9 = Gibraltar ,  10 = Austr.-NewZe
!-----
!ic0  do 520 nhsf=4,min(66,ndhsf)

#if ( FRAZER_ARCTIC == 0 )
      do nhsf=4,min(7,ndhsf)
#else
!dmr --- Following line is for computation defined by Frazer Davies
      do nhsf=4,min(8,ndhsf)
!dmr ---
#endif
       if (iehsf(nhsf).eq.0) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Integre (u.dy.dz) :
        i = ishsf(nhsf)
        do k=ks1,ks2
          vvn = 0.
          vvp = 0.
          do j=jshsf(nhsf),jehsf(nhsf)
            vvn = vvn + tmu(i,j,k) * cmy(i,j,3) * min(zero, u(i,j,k))
            vvp = vvp + tmu(i,j,k) * cmy(i,j,3) * max(zero, u(i,j,k))
          enddo 
          vinfor(nv+1) = vinfor(nv+1) + vvn * dz(k)
          vinfor(nv+2) = vinfor(nv+2) + vvp * dz(k)
        enddo 
        vinfor(nv+1) = vinfor(nv+1) * dy * svrdrp
        vinfor(nv+2) = vinfor(nv+2) * dy * svrdrp
       else
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Integre (v.dx.dz) :
        j = jshsf(nhsf)
        do k=ks1,ks2
          vvn = 0.
          vvp = 0.
          do i=ishsf(nhsf),iehsf(nhsf)
            vvn = vvn + tmu(i,j,k) * cmx(i,j,3) * min(zero, v(i,j,k))
            vvp = vvp + tmu(i,j,k) * cmx(i,j,3) * max(zero, v(i,j,k))
          enddo
          vinfor(nv+1) = vinfor(nv+1) + vvn * dz(k)
          vinfor(nv+2) = vinfor(nv+2) + vvp * dz(k)
        enddo 
        vinfor(nv+1) = vinfor(nv+1) * dx * svrdrp
        vinfor(nv+2) = vinfor(nv+2) * dx * svrdrp
       endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
       nv = nv + 2
      enddo 
       
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--2.5 diagnostics of the THC circulation
!- Atlantique Prod GIN (= le Max entre 68 et 75 N)
      zold=0
      if (zold.eq.1) then
      yy = 68.
      jj1= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      yy = 75.
      jj2= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      do jj=jj1,jj2
       vvpsm = 0.0
       do k=ks1,ks2
         if (z(k).gt.zz) cycle
         vv = 0.0
         do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
           vv = vv + cmx(i,jj,3)*v(i,jj,k)
         enddo
         vvpsm = vvpsm - vv * dz(k)
         vinfor(nn) = max(vinfor(nn),vvpsm)
       enddo
      enddo
      vinfor(nn) = vinfor(nn) * dx * svrdrp
      nv = nv + 1
!- Atlantique Prod (= le Max entre 45 et 65 N)
      yy = 45.
      jj1= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      yy = 65.
      jj2= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      do jj=jj1,jj2
       vvpsm = 0.0
       do k=ks1,ks2
         if (z(k).gt.zz) cycle
         vv = 0.0
         do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
          vv = vv + cmx(i,jj,3)*v(i,jj,k)
         enddo
         vvpsm = vvpsm - vv * dz(k)
         vinfor(nn) = max(vinfor(nn),vvpsm)
       enddo
      enddo  
      vinfor(nn) = vinfor(nn) * dx * svrdrp
      nv = nv + 1
!- NADW Outflow and AABW inflow into the Atlantic
!    (= le Max et le min a 20 S)
      yy = -20.0
      jj = 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vvpsm = 0.0
      do k=ks1,ks2
        if (z(k).gt.zz) cycle
        vv = 0.0
        do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
          vv = vv + cmx(i,jj,3)*v(i,jj,k)
        enddo
        vvpsm = vvpsm - vv * dz(k)
        vinfor(nn)  = max(vinfor(nn),vvpsm)
        vinfor(nn+3)= min(vinfor(nn+3),vvpsm)
      enddo
      vinfor(nn) = vinfor(nn) * dx * svrdrp
      vinfor(nn+3) = vinfor(nn+3) * dx * svrdrp
      nv = nv + 1
!- AABW Prod (= le Min Sud de 60S )
      yy = -70.0
      jj1= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      yy = -60.0
      jj2= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      do jj=jj1,jj2
       vvpsm = 0.0
       do k=ks1,ks2
         if (z(k).gt.zz) cycle
         vv = 0.0
         do i=iu1(jj),iu2(jj)
           vv = vv + cmx(i,jj,3)*v(i,jj,k)
         enddo
         vvpsm = vvpsm - vv * dz(k)
         vinfor(nn) = min(vinfor(nn),vvpsm)
       enddo
      enddo 
      vinfor(nn) = vinfor(nn) * dx * svrdrp
      nv = nv + 1
!- AABW export (= le Min a 30S)
      yy = -30.0
      jj= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vvpsm = 0.0
      do k=ks1,ks2
        if (z(k).gt.zz) cycle
        vv = 0.0
         do i=iu1(jj),iu2(jj)
          vv = vv + cmx(i,jj,3)*v(i,jj,k)
         enddo
        vvpsm = vvpsm - vv * dz(k)
        vinfor(nn)  = min(vinfor(nn),vvpsm)
      enddo
      vinfor(nn) = vinfor(nn) * dx * svrdrp
      nv = nv + 2
      else
!     write(mouchard_id,*) start new diag
      nm1n=nv
!- Atlantic Prod GIN (= Max between 68 and 75 N)
      yy = 68
      jj1= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      yy = 75.
      jj2= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vinfor(nn)=0.0
      do jj=jj1,jj2
        do k=ks1,ks2
         if (z(k).gt.zz) cycle
         vinfor(nn) = max(vinfor(nn),(vwx(jj,k,3)))
        enddo
      enddo
      vinfor(nn) = vinfor(nn) * dx * svrdrp
      nv = nv + 1
!- Atlantique Prod (= Max between 45 and 75 N)
      yy = 45.
      jj1= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      yy = 75.
      jj2= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vinfor(nn)=0.0
      do jj=jj1,jj2
        do k=ks1,ks2
         if (z(k).gt.zz) cycle
!        write(mouchard_id,*) jj,k,vinfor(nn),vwx(jj,k,3)
         vinfor(nn) = max(vinfor(nn),(vwx(jj,k,3)))
        enddo
      enddo
      moc=vinfor(nn)
      nv = nv + 1   ! =20 ?
!- NADW Outflow and AABW inflow into the Atlantic
!    (= Max and  min at 20 S)
      yy = -20.0
      jj = 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vinfor(nn)=0.0
      vinfor(nn+3)=0.0
      do k=ks1,ks2
        if (z(k).gt.zz) cycle
        vinfor(nn)  = max(vinfor(nn),(vwx(jj,k,3)))
        vinfor(nn+3)= min(vinfor(nn+3),(vwx(jj,k,3)))
      enddo
!     write(mouchard_id,*) 'NADW,AABW',vinfor(nn),vinfor(nn+3)
      nv = nv + 1
!- AABW Prod (= Min southward of 60S )
      yy = -70.0
      jj1= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      yy = -60.0
      jj2= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vinfor(nn)=0.0
      do jj=jj1,jj2
        do k=ks1,ks2
         if (z(k).gt.zz) cycle
         vinfor(nn) = min(vinfor(nn),(vwx(jj,k,0)))
        enddo
      enddo
      nv = nv + 1
!- AABW export (= Min at 30S)
      yy = -30.0
      jj= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      zz = -500.0
      nn = nv + 1
      vinfor(nn)=0.0
      do k=ks1,ks2
       if (z(k).gt.zz) cycle
       vinfor(nn) = min(vinfor(nn),(vwx(jj,k,0)))
      enddo
      nv = nv + 2
      endif
      nm2n=nv

!- Heat flux at 30 S
      yy = -30.0
      jj= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      vvpsm = 0.0
      ssc2 = 2.0 * scal0(ks1,1)
!     ssc2 = 2.0 * 273.15
      do k=ks1,ks2
        vv = 0.0
        ccydif = ahs(k) * unsdy
        do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
          vv2 = 0.5 * (v(i,jj,k)+v(i+1,jj,k))
     &        + viso(i,jj,k)
          v2cd2 = 0.5 * cmx(i,jj,2) * vv2
          cny = smy(i,jj,2) * unsdy * dts(k) * vv2
          aalpha = min( one, abs(cny) + alphmi(k) )
          phiy = tms(i,jj,k) * tms(i,jj-1,k) * (
     &         v2cd2 * (scal(i,jj-1,k,1) + scal(i,jj,k,1) - ssc2)
!    &         + ( alphay(i,jj,k) + cmxy(i,jj,2)*ccydif )
     &         + (aalpha*abs(v2cd2) + cmxy(i,jj,2)*ccydif)
     &               * (scal(i,jj-1,k,1) - scal(i,jj,k,1))  )
!-----
!         phi1=tms(i,jj,k) * tms(i,jj-1,k) * (
!    &          v2cd2 * (scal(i,jj-1,k,1) + scal(i,jj,k,1) - ssc2 ))
!         phi2=tms(i,jj,k) * tms(i,jj-1,k) * (
!    &        cmxy(i,jj,2) * ccydif
!    &               * (scal(i,jj-1,k,1) - scal(i,jj,k,1))  )
!         write(95,*) phi1,phi2,phiy
          vv = vv + phiy
        enddo
        vvpsm = vvpsm + vv*dz(k)
      enddo
      nv = nv + 1
      vinfor(nv) = vvpsm * dx * rho0*cpo * 1.D-15

!- Calculate Atlantic average salinity for Atlantic salt budget
!calculations
      saltatlantic=0.0
      volatlantic=0.0
       do j=1,imax
        do i=iszon(j,nbsmax),iezon(j,nbsmax)
         do k=1,kmax
          saltatlantic=saltatlantic+aire(i,j)*dz(k)*
     &     tms(i,j,k)*scal(i,j,k,2) !m2*m*g/kg -> m3*g/kg
          volatlantic=volatlantic+aire(i,j)*dz(k)*tms(i,j,k)
          enddo
        enddo
      enddo
      saltatlantic=saltatlantic/volatlantic !m3*g/kg/m3 -> g/kg
!      write(mouchard_id,*) 'Atlantic mean
!     &                       salinity [g/kg]',saltatlantic

!- Salt flux at 30 S (Fs30A)
      yy = -30.0
      jj= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      vvpsm = 0.0
      ssc2 = 2.0 * scal0(ks1,2) ! should we use saltatlantic as reference salinity?
      do k=ks1,ks2
        vv = 0.0
        ccydif = ahs(k) * unsdy
        do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
         vv2 = 0.5 * (v(i,jj,k)+v(i+1,jj,k))
     &         + viso(i,jj,k)
         v2cd2 = 0.5 * cmx(i,jj,2) * vv2
         cny = smy(i,jj,2) * unsdy * dts(k) * vv2
         aalpha = min( one, abs(cny) + alphmi(k) )
         phiy = tms(i,jj,k) * tms(i,jj-1,k) * (
     &       v2cd2 * (scal(i,jj-1,k,2) + scal(i,jj,k,2) - ssc2)
     &       + (aalpha*abs(v2cd2) + cmxy(i,jj,2)*ccydif)
     &       * (scal(i,jj-1,k,2) - scal(i,jj,k,2))  )
         vv = vv + phiy
       enddo
       vvpsm = vvpsm + vv*dz(k)
      enddo
      nv = nv + 1
      vinfor(nv) = vvpsm * dx
      
!- Advective salt flux at 30 S (Fsa30A)
      yy = -30.0
      jj= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      vvpsm = 0.0
      do k=ks1,ks2
        vv = 0.0
        ccydif = ahs(k) * unsdy
        do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
         vv2 = 0.5 * (v(i,jj,k)+v(i+1,jj,k))
     &         + viso(i,jj,k)
         v2cd2 = 0.5 * cmx(i,jj,2) * vv2 ! why 'times 0.5?
         phiy = tms(i,jj,k) * tms(i,jj-1,k) * (v2cd2 * 
     &       0.5*(scal(i,jj-1,k,2) + scal(i,jj,k,2)) )
         vv = vv + phiy
       enddo
       vvpsm = vvpsm + vv*dz(k)
      enddo
      nv = nv + 1
      vinfor(nv) = vvpsm * dx
!      write(mouchard_id,*) 'Advective salt flux at 30 S
!     &                   [1e9 kg/s]', vinfor(nv)/1e9
       
!- Diffisive salt flux at 30 S (Fsd30A)
      vvpsm = 0.0
      do k=ks1,ks2
        vv = 0.0
        ccydif = ahs(k) * unsdy
        do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
         vv2 = 0.5 * (v(i,jj,k)+v(i+1,jj,k))
     &         + viso(i,jj,k)
         v2cd2 = 0.5 * cmx(i,jj,2) * vv2
         cny = smy(i,jj,2) * unsdy * dts(k) * vv2
         aalpha = min( one, abs(cny) + alphmi(k) )
         phiy = tms(i,jj,k) * tms(i,jj-1,k) * (
     &       (aalpha*abs(v2cd2) + cmxy(i,jj,2)*ccydif)
     &       * (scal(i,jj-1,k,2) - scal(i,jj,k,2))  )
         vv = vv + phiy
       enddo
       vvpsm = vvpsm + vv*dz(k)
      enddo
      nv = nv + 1
      vinfor(nv) = vvpsm * dx     
!      write(mouchard_id,*) 'Diffisive salt flux at 30 S 
!     &                   [1e9 kg/s]', vinfor(nv)/1e9

! Advective and diffusive salt fluxes at northern boundary of the Atlantic basin (Fram
! Strait + Barents Strait + Canadian Archipelico)
      vvtfa = 0.0
      vvtfd = 0.0
      vvtba = 0.0
      vvtbd = 0.0
      vvtca = 0.0
      vvtcd = 0.0

      do k=ks1,ks2
         vva = 0.0
         vvd = 0.0
         !Fram Strait 
         do ii=106,107
              jj=55   
              ccydif = ahs(k) * unsdy
              vv2 = 0.5 * (v(ii,jj,k)+v(ii+1,jj,k))
     &                + viso(ii,jj,k)
              v2cd2 = 0.5 * cmx(ii,jj,2) * vv2
              ! advective part
              phiy = tms(ii,jj,k) * tms(ii,jj-1,k) * (
     &             v2cd2 * 0.5*(scal(ii,jj-1,k,2) + scal(ii,jj,k,2)))
              vva = vva + phiy
              ! diffusive part
              cny = smy(ii,jj,2) * unsdy * dts(k) * vv2
              aalpha = min( one, abs(cny) + alphmi(k) )
              phiy = tms(ii,jj,k) * tms(ii,jj-1,k) * (
     &                 (aalpha*abs(v2cd2) + cmxy(ii,jj,2)*ccydif)
     &                * (scal(ii,jj-1,k,2) - scal(ii,jj,k,2))  )
              vvd = vvd + phiy
         enddo
         vvtfa = vvtfa + vva*dz(k)
         vvtfd = vvtfd + vvd*dz(k)*(-1) ! *-1 according to eqation 9 in Ragani & Dijkstra 2025

         ! Barents Strait
         vva = 0.0
         vvd = 0.0
         do jj=54,56
              ii=110
              ccxdif = ahs(k) * unsdx
              uu2 = 0.5 * (u(ii,jj,k)+u(ii,jj+1,k))
     &                + uiso(ii,jj,k)
              u2cd2 = 0.5 * cmy(ii,jj,2) * uu2
              ! advective part
              phix = tms(ii,jj,k) * tms(ii-1,jj,k) * (
     &               u2cd2 * 0.5*(scal(ii-1,jj,k,2)+scal(ii,jj,k,2)))
              vva = vva + phix
              ! diffusive part
              cnx = smx(ii,jj,2) * unsdx * dts(k) * uu2
              aalpha = min( one, abs(cnx) + alphmi(k) )
              phix = tms(ii,jj,k) * tms(ii-1,jj,k) * (
     &                 (aalpha*abs(u2cd2) + cmxy(ii,jj,2)*ccxdif)
     &                * (scal(ii-1,jj,k,2) - scal(ii,jj,k,2))  )
              vvd = vvd + phix
         enddo
         vvtba = vvtba + vva*dz(k)
         vvtbd = vvtbd + vvd*dz(k)*(-1) ! *-1 according to eqation 9 in Ragani & Dijkstra 2025

         ! Canadian Archipelico
         vva = 0.0
         vvd = 0.0
         do ii=101,102
               jj=57
               ccydif = ahs(k) * unsdy
               vv2 = 0.5 * (v(ii,jj,k)+v(ii+1,jj,k))
     &                 + viso(ii,jj,k)
               v2cd2 = 0.5 * cmx(ii,jj,2) * vv2
               ! advective part         
               phiy = tms(ii,jj,k) * tms(ii,jj-1,k) * (
     &             v2cd2 * 0.5*(scal(ii,jj-1,k,2) + scal(ii,jj,k,2)))
               vva = vva + phiy
               ! diffusive part
               cny = smy(ii,jj,2) * unsdy * dts(k) * vv2
               aalpha = min( one, abs(cny) + alphmi(k) )
               phiy = tms(ii,jj,k) * tms(ii,jj-1,k) * (
     &                   (aalpha*abs(v2cd2) + cmxy(ii,jj,2)*ccydif)
     &                   * (scal(ii,jj-1,k,2) - scal(ii,jj,k,2))  )
               vvd = vvd + phiy
         enddo
         vvtca = vvtca + vva*dz(k)
         vvtcd = vvtcd + vvd*dz(k)*(-1) ! *-1 according to eqation 9 in Ragani & Dijkstra 2025

      enddo ! k=ks1,ks2

      nv = nv + 1
      vinfor(nv) = ((vvtfa * dx)
     &             +  (vvtba * dy)
     &             +  (vvtca * dx))!/1000*scal0(ks1,2))
!      write(mouchard_id,*) 
!     &            'Advective salt flux north. bound. [1e9 kg/s]',
!     & vinfor(nv)/1e9
      nv = nv + 1
      vinfor(nv) = ((vvtfd * dx)
     &             +  (vvtbd * dy)
     &             +  (vvtcd * dx))!/1000*scal0(ks1,2))
!      write(mouchard_id,*) 
!     &            'Diffusive salt flux north. bound. [1e9 kg/s]',
!     &            vinfor(nv)/1e9

! Atlantic salt content tendency
      saltatlantic_new = 0.0
      saltatlantic_dt = 0.0
      do j=1,imax
        do i=iszon(j,nbsmax),iezon(j,nbsmax)
         do k=1,kmax
          saltatlantic_new = saltatlantic_new+aire(i,j)*dz(k)*
     &     tms(i,j,k)*scal(i,j,k,2) !m2*m*g/kg -> m3*g/kg
          enddo
        enddo
      enddo
      saltatlantic_dt = (saltatlantic_new-saltatlantic_old)/ddtb
      saltatlantic_old = saltatlantic_new
!      write(mouchard_id,*) 'Atlantic salt tendency [1e9 kg/s]'
!     &                    ,saltatlantic_dt/1e9
      nv = nv + 1
      vinfor(nv) = saltatlantic_dt


!- Mov parameter at 30 S (meridional freshwater flux at 30S related to overturning circulation)
      yy = -30.0 ! Set latitude
      jj= 1 + nint( (yy - ylat1) / dlat + 0.5 )
      vvpsm = 0.0

      do k=ks1,ks2
        ! Get zonal mean salinity
        ss = 0.0
        do i=iszon(jj,nbsmax),iezon(jj,nbsmax) ! nbsmax==3==Atlantic
             ss = ss + 0.5 * (scal(i,jj-1,k,2) + scal(i,jj,k,2) )
        enddo
        sszonm = ss / (iezon(jj,nbsmax)-iszon(jj,nbsmax)) ! sszonm = zonal mean salinity at depth k

        ! Calculate Mov at 30S
        vv = 0.0 ! horizontal velocity?
        do i=iszon(jj,nbsmax),iezon(jj,nbsmax)
             vv2    = 0.5 * (v(i,jj,k)+v(i+1,jj,k))
     &                 + viso(i,jj,k)
             v2cd2 = 0.5 * cmx(i,jj,2) * vv2
             phiy = tms(i,jj,k) * tms(i,jj-1,k)
     &              * (v2cd2 * (sszonm - saltatlantic) )
             vv   = vv + phiy
         enddo
      vvpsm = vvpsm + vv*dz(k)
      enddo

      nv = nv + 1
      vinfor(nv) = vvpsm * dx * (-1/saltatlantic) / 1e6   ! Movs freshwater flux in Sv
!      vinfor(nv) = vvpsm * (-1/scal0(ks1,2)) * dx / 1e6   ! Movs freshwater flux in Sv
!      write(mouchard_id,*) 'Mov 30S [Sv]', vinfor(nv)
      
            
!- Salt flux at Bering Strait (fsber)
      jj=65
      ii=102
      vvpsm = 0.0
      ssc2 = 2.0 * scal0(ks1,2)
      vv = 0.0
      ccydif = ahs(ks2) * unsdy
        vv2 = 0.5 * (v(ii,jj,ks2)+v(ii+1,jj,ks2))
     &         + viso(ii,jj,ks2)
        v2cd2 = 0.5 * cmx(ii,jj,2) * vv2
        cny = smy(ii,jj,2) * unsdy * dts(ks2) * vv2
        aalpha = min( one, abs(cny) + alphmi(ks2) )
        phiy = tms(ii,jj,ks2) * tms(ii,jj-1,ks2) * (
     &       v2cd2 * (scal(ii,jj-1,ks2,2) + scal(ii,jj,ks2,2) - ssc2)
     &       + (aalpha*abs(v2cd2) + cmxy(ii,jj,2)*ccydif)
     &       * (scal(ii,jj-1,ks2,2) - scal(ii,jj,ks2,2))  )
         vv = vv + phiy
         vvpsm = vvpsm + vv*dz(ks2)

      nv = nv + 1  ! 26 ?
      vinfor(nv) = vvpsm * dx

   
      
#if ( FRAZER_ARCTIC == 1 )
!----------+----------+----------+----------+----------+----------+----------+----------+----------+
!---Added on 16th September 2011 by Frazer.J.Davies

!---CALCULATION OF FRESHWATER FLUX (FF) in (mSv) at various locations in the model.

!---NB: These are calculated using a reference salintiy of 34.8. This can altered if required.
!---NB: 1Sv = 10^6 m3s-1

!---FF = Intergral of velocity * ((SALref-salinity)/SALref)

!---For a more detailed explanantion for calculating the FF see(amongst others):

!---'Steele et al., (1996), A Simple model study of the Arctic Ocean freshwater balance, 1979-1985,
!---Journal of Geophyscial Research, Vol.101,NO.C9, pg 20833-20848'.
!---
!---'Holland et al., (2006), Simulated Arctic Ocean Freshwater Budgets in the 20th and 21st Centuries
!---Journal of Climate, Vol 19, Pg 6221-6242'.


!---These locations where the FF is calculated are as follows:

!---1)Bering Strait "BSFF"
!---2)Fram Strait "FSFF"
!---3)Barents Strait "BAFF" (Defined as the boudary between Svalbard and Norway).
!---4)CAA "CAAFF"
!---5)Kara sea "KSFF"

!---These variable save been added to the informe.f file.

!--The sign (+/-) of the "flux" term represenst a source (+) or sink (-) of freshwater
!----------+----------+----------+----------+----------+----------+----------+--
!---1) Bering Strait "BSFP & BSFN"

      jj=jberp
      ii=iberp
      vvpsm = 0.0
      ssc2 = 2.0 * 34.8
      vvn = 0.0
      vvp = 0.0
      vvt = 0.0
      vv = 0.0
      vvppsm = 0.0
      vvnpsm = 0.0
      vvtpsm = 0.0
      do k=ks1,ks2
      ccydif = ahs(k) * unsdy
        vv2 = 0.5 * (v(ii,jj,k)+v(ii+1,jj,k))
     &         + viso(ii,jj,k)
        v2cd2 = 0.5 * cmx(ii,jj,2) * vv2
        cny = smy(ii,jj,2) * unsdy * dts(k) * vv2
        aalpha = min( one, abs(cny) + alphmi(k) )
        phiy = tms(ii,jj,k) * tms(ii,jj-1,k) * (
     &       v2cd2 * (scal(ii,jj-1,k,2)+scal(ii,jj,k,2) - ssc2)
     &       + (aalpha*abs(v2cd2) + cmxy(ii,jj,2)*ccydif)
     &       * (scal(ii,jj-1,k,2) - scal(ii,jj,k,2))  )
        flux = (phiy/34.8)*(-1)

        vvt = vvt + flux
        vvp = vvp + max(0.0,flux)
        vvn = vvn + min(0.0,flux)

      vvtpsm = vvtpsm + vvt*dz(k)
      vvppsm = vvppsm + vvp*dz(k)
      vvnpsm = vvnpsm + vvn*dz(k)

      enddo
      nv = nv + 1
      vinfor(nv) = (((vvtpsm*dx)/1000)*31.104)
!      write(*,*), "bsff", vinfor(nv)
      nv = nv + 1
      vinfor(nv) = (((vvppsm * dx)/1000)*31.104)
!      write(*,*), "bsfp", vinfor(nv)
      nv = nv + 1
      vinfor(nv) = (((vvnpsm * dx)/1000)*31.104)
!      write(*,*), "bsfn", vinfor(nv)

!---2) Fram Strait  FSFP & FSFN----+----------+----------+----------+----------+
! PB advective and diffusive flux (sum == total)
      jj=55
      ssc2 = 2.0 * 34.8
      vvn = 0.0
      vvp = 0.0
      vvt = 0.0
      vvppsm = 0.0
      vvnpsm = 0.0
      vvtpsm = 0.0
      do k=ks1,ks2
      vvp = 0.0
      vvn = 0.0
      vvt = 0.0
        do 607 ii=106,107
          ccydif = ahs(k) * unsdy
          vv2 = 0.5 * (v(ii,jj,k)+v(ii+1,jj,k))
     &           + viso(ii,jj,k)
          v2cd2 = 0.5 * cmx(ii,jj,2) * vv2
          cny = smy(ii,jj,2) * unsdy * dts(k) * vv2
          aalpha = min( one, abs(cny) + alphmi(k) )
          phiy = tms(ii,jj,k) * tms(ii,jj-1,k) * (
     &         v2cd2 * (scal(ii,jj-1,k,2) + scal(ii,jj,k,2) - ssc2)
     &         + (aalpha*abs(v2cd2) + cmxy(ii,jj,2)*ccydif)
     &         * (scal(ii,jj-1,k,2) - scal(ii,jj,k,2))  )
          flux = (phiy/34.8)*(-1)
          vvt = vvt + flux
          vvp = vvp + max(0.0,flux)
          vvn = vvn + min(0.0,flux)

        enddo
      vvtpsm = vvtpsm + vvt*dz(k)
      vvppsm = vvppsm + vvp*dz(k)
      vvnpsm = vvnpsm + vvn*dz(k)
      enddo
      nv = nv + 1
      vinfor(nv) = (((vvtpsm * dx)/1000)*31.104)
!      write(*,*), "fsff", vinfor(nv)
      nv = nv + 1
      vinfor(nv) = (((vvppsm * dx)/1000)*31.104)
!      write(*,*), "fsfp", vinfor(nv)
      nv = nv + 1
      vinfor(nv) = (((vvnpsm * dx)/1000)*31.104)
!      write(*,*), "fsfn", vinfor(nv)

!---3) Barents Sea "BAFP & BAFN"----------+---------+----------+----------+----+

!-original start
      ii=110
      uupsm = 0.0
      ssc2 = 2.0 * 34.8
      uu = 0.0
      uun = 0.0
      uup = 0.0
      uutpsm = 0.0
      uuppsm = 0.0
      uunpsm = 0.0
      do k=ks1,ks2
      uun = 0.0
      uup = 0.0
      uut = 0.0
        do jj=54,56
          ccxdif = ahs(k) * unsdx
          uu2 = 0.5 * (u(ii,jj,k)+u(ii,jj+1,k))
     &           + uiso(ii,jj,k)
          u2cd2 = 0.5 * cmy(ii,jj,2) * uu2
          cnx = smx(ii,jj,2) * unsdx * dts(k) * uu2
          aalpha = min( one, abs(cnx) + alphmi(k) )
          phix = tms(ii,jj,k) * tms(ii-1,jj,k) * (
     &         u2cd2 * (scal(ii-1,jj,k,2)+scal(ii,jj,k,2) - ssc2)
     &         + (aalpha*abs(u2cd2) + cmxy(ii,jj,2)*ccxdif)
     &         * (scal(ii-1,jj,k,2) - scal(ii,jj,k,2))  )
          flux = (phix/34.8)*(-1)
          uut = uut + flux
          uup = uup + max(0.0,flux)
          uun = uun + min(0.0,flux)
        enddo
      uutpsm = uutpsm + uut*dz(k)
      uuppsm = uuppsm + uup*dz(k)
      uunpsm = uunpsm + uun*dz(k)
      enddo
      nv = nv + 1
      vinfor(nv) = (((uutpsm * dy)/1000)*31.104)
!      write(*,*), "baff", vinfor(nv)
      nv = nv + 1
      vinfor(nv) = (((uuppsm * dy)/1000)*31.104)
!      write(*,*), "bafp", vinfor(nv)
      nv = nv + 1
      vinfor(nv) = (((uunpsm * dy)/1000)*31.104)
!      write(*,*), "bafn", vinfor(nv)
!-------original end
!--test
!      ii=108
!      uupsm = 0.0
!      uu2 = 0.0
!      do k=ks1,ks2
!      do jj=53,55
!      ccxdif = ahs(k) * unsdx
!        uu2 = 0.5 * (u(ii,jj,k)+u(ii,jj+1,k))
!     &         + uiso(ii,jj,k)
!        u2cd2 = 0.5 * cmy(ii,jj,2) * uu2
!        cnx = smx(ii,jj,2) * unsdx * dts(k) * uu2
!        aalpha = min( one, abs(cnx) + alphmi(k) )
!        phix = tms(ii,jj,k) * tms(ii-1,jj,k) * (
!     &       u2cd2 * (scal(ii-1,jj,k,2)+scal(ii,jj,k,2) - ssc2)
!     &       + (aalpha*abs(u2cd2) + cmxy(ii,jj,2)*ccxdif)
!     &       * (scal(ii-1,jj,k,2) - scal(ii,jj,k,2))  )
!        sal = ((scal(ii-1,jj,k,2)+scal(ii,jj,k,2))*0.5)
!        flux = uu2*((34.8-sal)/34.8)
!        uup = uup + flux
!        uun = uun + min(0.0,flux)
!      enddo
!      uuppsm = uuppsm + uup*dz(k)
!      uunpsm = uunpsm + uun*dz(k)
!      enddo
!      nv = nv + 1
!      vinfor(nv) = (((uuppsm * dy)/1000)*31.104)
!      write(*,*),"flux", vinfor(nv)
!----advection (adve)
      ii=108
      uupsm = 0.0
      ssc2 = 2.0 * 34.8
      uu = 0.0
      uun = 0.0
      uup = 0.0
      adv = 0.0
      div = 0.0
      advt =0.0
      do k=ks1,ks2
      adv = 0.0
        do jj=53,55
          ccxdif = ahs(k) * unsdx
          uu2 = 0.5 * (u(ii,jj,k)+u(ii,jj+1,k))
     &           + uiso(ii,jj,k)
          u2cd2 = 0.5 * cmy(ii,jj,2) * uu2
          cnx = smx(ii,jj,2) * unsdx * dts(k) * uu2
          aalpha = min( one, abs(cnx) + alphmi(k) )
          ad = (u2cd2 * (scal(ii-1,jj,k,2)+scal(ii,jj,k,2) - ssc2))
          advflux = (ad/34.8)*(-1)
          adv = adv + advflux
        enddo
       advt = advt + adv * dz(k)
      enddo
      nv = nv + 1
      vinfor(nv) = (((advt*dy)/1000)*31.104)
!----diffusion (diff)
      ii=108
      uupsm = 0.0
      ssc2 = 2.0 * 34.8
      uu = 0.0
      uun = 0.0
      uup = 0.0
      adv = 0.0
      div = 0.0
      divt= 0.0
      do k=ks1,ks2
      div = 0.0
        do jj=53,55
          ccxdif = ahs(k) * unsdx
          uu2 = 0.5 * (u(ii,jj,k)+u(ii,jj+1,k))
     &           + uiso(ii,jj,k)
          u2cd2 = 0.5 * cmy(ii,jj,2) * uu2
          cnx = smx(ii,jj,2) * unsdx * dts(k) * uu2
          aalpha = min( one, abs(cnx) + alphmi(k) )
          di = ((aalpha*abs(u2cd2) + cmxy(ii,jj,2)*ccxdif)
     &         * (scal(ii-1,jj,k,2) - scal(ii,jj,k,2))  )
          divflux = (di/34.8)*(-1)
          div = div + divflux
        enddo
       divt = divt + div * dz(k)
      enddo
      nv = nv + 1
      vinfor(nv) = (((divt*dy)/1000)*31.104)

!---west Arctic Ocean Freshwater Content (AOW)
      jj=58
      ii=100
      aofc = 0.0
      aofc1 = 0.0
      aofc2 = 0.0
      aofc3 = 0.0
      vol = 0.0
      do k=ks1,ks2
      aofc1 = 0.0
        do ii=100,105
          aofc2 = 0.0
           do jj=55,65
             vol = (dx * dy * dz(k) * tms(ii,jj,k) * 1e-9)
             aofc = ((1-(scal(ii,jj,k,2)/34.8)) * vol)
             aofc1 = aofc1 + aofc
!            aofc1 = aofc1 + max(0.0,aofc)
           enddo
        aofc2 = aofc2 + aofc1
        enddo
      aofc3 = aofc3 + aofc2
      enddo
      nv = nv + 1
      vinfor(nv) = aofc3
*      write(*,*), "aow", vinfor(nv)

*      ii=108
*      sal = 0.0
*      salt = 0.0
*      tz = 0.0
*      do jj=53,54
*      tz = 0.0
*      sal = 0.0
*      do  k=ks1,ks2
*      sali =  (tms(ii,jj,k) * scal(ii,jj,k,2)) * dz(k)
*        sal = sal + sali
*        tz = tz + dz(k) * tms(ii,jj,k)
*      enddo
*       salt = salt + (sal/tz)
*      enddo
*      nv = nv + 1
*      vinfor(nv) = salt/2

!---East Arctic Ocean Freshwater Content (AOE)
      jj=58
      ii=100
      aofc = 0.0
      aofc1 = 0.0
      aofc2 = 0.0
      aofc3 = 0.0
      vol = 0.0
      do k=ks1,ks2
        aofc1 = 0.0
        do ii=106,112
          aofc2 = 0.0
           do jj=55,65
             vol = (dx * dy * dz(k) * tms(ii,jj,k) * 1e-9)
             aofc = ((1-(scal(ii,jj,k,2)/34.8)) * vol)
             aofc1 = aofc1 + aofc
!            aofc1 = aofc1 + max(0.0,aofc)
           enddo
        aofc2 = aofc2 + aofc1
        enddo
      aofc3 = aofc3 + aofc2
      enddo
      nv = nv + 1
      vinfor(nv) = aofc3
*      write(*,*), "aoe", vinfor(nv)

*      ii=110
*      ssc2 = 2.0 * 34.8
*      uu = 0.0
*      velo = 0.0
*      vel = 0.0
*      do jj=54,56
*      tz = 0.0
*      velo = 0.0
*      do k=ks1,ks2
*      uu2 = 0.5 * (u(ii,jj,k) + u(ii,jj+1,k))
*      vel = tms(ii,jj,k) * (0.5 * u(ii,jj,k) + u(ii,jj+1,k)) * dz(k)
*      velo = velo + vel
*      tz = tz +dz(k) * tms(ii,jj,k)

*      enddo
*       ve = ve + (velo/tz)
*      enddo
*      nv = nv + 1
*      vinfor(nv) = ve/3

!
      jj=58
      ii=100
*these values for ii and jj make no difference
      aofc = 0.0
      aofc1 = 0.0
      aofc2 = 0.0
      aofc3 = 0.0
      vol = 0.0
      do k=ks1,ks2
      aofc1 = 0.0
        do ii=100,112
          aofc2 = 0.0
           do jj=55,65
             vol = (dx * dy * dz(k) * tms(ii,jj,k) * 1e-9)
             aofc = ((1-(scal(ii,jj,k,2)/34.8)) * vol)
             aofc1 = aofc1 + aofc
!            aofc1 = aofc1 + max(0.0,aofc)
           enddo
        aofc2 = aofc2 + aofc1
        enddo
      aofc3 = aofc3 + aofc2
      enddo
      nv = nv + 1
      vinfor(nv) = aofc3
*      write(*,*), "aofc", vinfor(nv)

!*

!---4) CAA "CAFP & CAFN"----------+----------+----------+----------+----------+

      jj=57
      vvpsm = 0.0
      ssc2 = 2.0 * 34.8
      vv = 0.0
      vvp = 0.0
      vvn = 0.0
      vvt = 0.0
      vvtpsm = 0.0
      vvppsm = 0.0
      vvnpsm = 0.0
      do k=ks1,ks2
        vvp = 0.0
        vvn = 0.0
        vvt = 0.0
        do ii=101,102
          ccydif = ahs(k) * unsdy
          vv2 = 0.5 * (v(ii,jj,k)+v(ii+1,jj,k))
     &           + viso(ii,jj,k)
          v2cd2 = 0.5 * cmx(ii,jj,2) * vv2
          cny = smy(ii,jj,2) * unsdy * dts(k) * vv2
          aalpha = min( one, abs(cny) + alphmi(k) )
          phiy = tms(ii,jj,k) * tms(ii,jj-1,k) * (
     &         v2cd2 * (scal(ii,jj-1,k,2) + scal(ii,jj,k,2) - ssc2)
     &         + (aalpha*abs(v2cd2) + cmxy(ii,jj,2)*ccydif)
     &         * (scal(ii,jj-1,k,2) - scal(ii,jj,k,2))  )
          flux = (phiy/34.8)*(-1)
          vvt = vvt + flux
          vvp = vvp + max(0.0,flux)
          vvn = vvn + min(0.0,flux)
        enddo
        vvtpsm = vvtpsm + vvt*dz(k)
        vvppsm = vvppsm + vvp*dz(k)
        vvnpsm = vvnpsm + vvn*dz(k)
      enddo
      nv = nv + 1
      vinfor(nv) = (((vvtpsm * dx)/1000)*31.104)
!      write(*,*), "caff", vinfor(nv)
      nv = nv + 1
      vinfor(nv) = (((vvppsm * dx)/1000)*31.104)
!      write(*,*), "cafp", vinfor(nv)
      nv = nv + 1
      vinfor(nv) = (((vvnpsm * dx)/1000)*31.104)
!      write(*,*), "cafn", vinfor(nv)

!---5) Kara Sea "KSFP & KSFN"----------+----------+----------+----------+

      jj=55
      vvpsm = 0.0
      ssc2 = 2.0 * 34.8
      vv = 0.0
      vvp = 0.0
      vvn = 0.0
      vvt = 0.0
      vvtpsm = 0.0
      vvppsm = 0.0
      vvnpsm = 0.0
      do k=ks1,ks2
        vvn = 0.0
        vvp = 0.0
        vvt = 0.0
        do ii=109,112
          ccydif = ahs(k) * unsdy
          vv2 = 0.5 * (v(ii,jj,k)+v(ii+1,jj,k))
     &           + viso(ii,jj,k)
          v2cd2 = 0.5 * cmx(ii,jj,2) * vv2
          cny = smy(ii,jj,2) * unsdy * dts(k) * vv2
          aalpha = min( one, abs(cny) + alphmi(k) )
          phiy = tms(ii,jj,k) * tms(ii,jj-1,k) * (
     &         v2cd2 * (scal(ii,jj-1,k,2) + scal(ii,jj,k,2) - ssc2)
     &         + (aalpha*abs(v2cd2) + cmxy(ii,jj,2)*ccydif)
     &         * (scal(ii,jj-1,k,2) - scal(ii,jj,k,2))  )
          flux = (phiy/34.8)*(-1)
          vvt = vvt + flux
          vvp = vvp + max(0.0,flux)
          vvn = vvn + min(0.0,flux)
        enddo
        vvtpsm = vvtpsm + vvt*dz(k)
        vvppsm = vvppsm + vvp*dz(k)
        vvnpsm = vvnpsm + vvn*dz(k)
      enddo
      nv = nv + 1
      vinfor(nv) = (((vvtpsm * dx)/1000)*31.104)
!      write(*,*), "ksff", vinfor(nv)
      nv = nv + 1
      vinfor(nv) = (((vvppsm * dx)/1000)*31.104)
!      write(*,*), "ksfp", vinfor(nv)
      nv = nv + 1
      vinfor(nv) = (((vvnpsm * dx)/1000)*31.104)
!      write(*,*), "ksfn", vinfor(nv)

#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- downslpoing out of the antarctic shelf ( from j=1 to jmsud ) :
!d     nv = nv + 1
!d    vinfor(nv) = 0.0
!d    ijlim = jmsud*imax
!d    do nnt=1,nslpdw
!d     nl=nslp(nnt)
!d     if ( kslp(nl).ge.ksud .and. ijslp(nl).le.ijlim )
!d   &   vinfor(nv) = vinfor(nv) + abs(uvcslp(nnt))
!d    enddo
!d    vinfor(nv) = vinfor(nv) * dx * svrdrp
!- downslpoing out of the arctic ( from j=jmnor to jmax ) :
!d    nv = nv + 1
!d    vinfor(nv) = 0.0
!d    ijlim = 1+(jmnor-1)*imax
!d    do nnt=1,nslpdw
!d     nl=nslp(nnt)
!d     if ( kslp(nl).ge.knor .and. ijslp(nl).ge.ijlim )
!d   &   vinfor(nv) = vinfor(nv) + abs(uvcslp(nnt))
!d    enddo
!d    vinfor(nv) = vinfor(nv) * dx * svrdrp
!d
      if (flgout .or. ntmoy.eq.2) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--2.6 Integrale globale a faible variabilite temporelle :

!--Moyenne sur tout le bassin des Scal., |Scal-Obs| et Scal.-Obs (Niv 1) :
      nnp = 3
      if (ninfo.lt.0) nnp = 2
      nnvk(2) = nocean + ks1   ! 76 + 1 = 77
      nnvk(1) = nnvk(2) - ks2  ! 77 - 20 = 57
! nv = 26
      do ns=1,nsmax   ! 1=temp, 2=salt
!- 1er niveau (surface) :
        k = ks2    ! =20
          vvk(1) = 0.
          vvk(2) = 0.
          vvk(3) = 0.
          do j=js1,js2
           do i=is1(j),is2(j)
            ctmobs = min( ctmi(i,j,k,0), (scalr(i,j,k,ns)-spvr) )
            vvk(1) = vvk(1) + ctmi(i,j,k,0) * scal(i,j,k,ns)
            vv2 = ctmobs * (scal(i,j,k,ns) - scalr(i,j,k,ns))
            vvk(2) = vvk(2) + vv2
            vvk(3) = vvk(3) + abs(vv2)
           enddo
          enddo
          vinfor(nv+1) = vvk(1) * dz(k)   !27/30
          vinfor(nv+3) = vvk(3) * dz(k)   !29/32
          vinfor(nv+2) = vvk(2) * zsurfo(k) !28/31
          if (ns.le.2) vinfor(nnvk(ns)-k) = vvk(3) * zsurfo(k) ! 37/57
!- autres niveaux :
        do k=ks2-1,ks1,-1
          vvk(1) = 0.
          vvk(2) = 0.
          vvk(3) = 0.
          do j=js1,js2
           do i=is1(j),is2(j)
            ctmobs = min( ctmi(i,j,k,0), (scalr(i,j,k,ns)-spvr) )
            vvk(1) = vvk(1) + ctmi(i,j,k,0) * scal(i,j,k,ns)
            vv2 = ctmobs * (scal(i,j,k,ns) - scalr(i,j,k,ns))
            vvk(2) = vvk(2) + vv2
            vvk(3) = vvk(3) + abs(vv2)
           enddo
          enddo
          vinfor(nv+1) = vinfor(nv+1) + vvk(1) * dz(k) !27/30
          vinfor(nv+3) = vinfor(nv+3) + vvk(3) * dz(k) !29/32
          if (ns.le.2) vinfor(nnvk(ns)-k) = vvk(nnp) * zsurfo(k) ! 37/57
        enddo
        vinfor(nv+1) = vinfor(nv+1) * zvols + scalwr(ns)  !27/30
        if (ns.eq.1) tmc=vinfor(nv+1)*cpo*rho0  ! TMC?
        vinfor(nv+3) = vinfor(nv+3) * zvolo  !29/32
        nv = nv + 3
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--Moyenne sur tout le bassin des |w| (m/s) :
      nv = nv + 1 ! 33
      do k=ks1+1,ks2
        vvk(1)=0.
        do j=js1,js2
         do i=is1(j),is2(j)
          vvk(1) = vvk(1) + ctmi(i,j,k-1,0) * abs(w(i,j,k))
         enddo
        enddo  
        vinfor(nv) = vinfor(nv) + vvk(1) * dzw(k)
      enddo
      vinfor(nv) = vinfor(nv) * zvolw

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--Moyenne sur tout le bassin des Vitesses (m/s) :
!     nvkek = nnvk(1) - ks2
      factke = rho0 * 0.5
      do k=ks1,ks2
        do n=1,3
          vvk(n) = 0.
        enddo
        do j=ju1,ju2
         do i=iu1(j),iu2(j)
          vvk(1) = vvk(1) + cmxy(i,j,3) * abs(u(i,j,k))
          vvk(2) = vvk(2) + cmxy(i,j,3) * abs(v(i,j,k))
          vvk(3) = vvk(3) + cmxy(i,j,3) *
     &           ( u(i,j,k) * u(i,j,k) + v(i,j,k) * v(i,j,k) )
         enddo
        enddo
        vvk(3) = vvk(3) * factke
!       vinfor(nvkek-k) = vvk(3) * zsurfv(k)
        do n=1,3
          vinfor(nv+n) = vinfor(nv+n) + vvk(n) * dz(k)
        enddo
      enddo
      do n=1,3
        vinfor(nv+n) = vinfor(nv+n) * zvolv
      enddo
      nv = nv + 3
!-----
      endif
      nv = nocean

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--2.7 Integrale globale (ou partielle) a forte variabilite temporelle :

      do j=jeq,js2
        do i=is1(j),is2(j)
          if (tms(i,j,ks2).eq.1) then
            vinfor(nv+1)=vinfor(nv+1)+(1.0-albq(i,j))*aire(i,j)
            if (albq(i,j).lt.0.15) vinfor(nv+3)=vinfor(nv+3) + aire(i,j)
            if (albq(i,j).lt.0.85) vinfor(nv+5)=vinfor(nv+5) + aire(i,j)
            vinfor(nv+7)=vinfor(nv+7)+albq(i,j)*aire(i,j)*(1-albq(i,j))
     &                        /max((1-albq(i,j)),1e-6*one)
            vinfor(nv+9)=vinfor(nv+9)
     &                   +(1.0-albq(i,j))*aire(i,j)*hgbq(i,j)
            vinfor(nv+11)=vinfor(nv+11)
     &                   +(1.0-albq(i,j))*aire(i,j)*hnbq(i,j)
            vinfor(nv+13)=vinfor(nv+13)
     &                   +(1.0-albq(i,j))*aire(i,j)*hgbq(i,j)*
     &                    (ug(i,j)*ug(i,j)+vg(i,j)*vg(i,j))
          endif
        enddo
      enddo
      vinfor(nv+13)=sqrt(vinfor(nv+13)/max(vinfor(nv+9),1e-5*one))
      nv=nv+1

      do j=js1,jeq-1
        do i=is1(j),is2(j)
          if (tms(i,j,ks2).eq.1) then
            vinfor(nv+1)=vinfor(nv+1)+(1.0-albq(i,j))*aire(i,j)
            if (albq(i,j).lt.0.15) vinfor(nv+3)=vinfor(nv+3) + aire(i,j)
            if (albq(i,j).lt.0.85) vinfor(nv+5)=vinfor(nv+5) + aire(i,j)
            vinfor(nv+7)=vinfor(nv+7)+albq(i,j)*aire(i,j)*(1-albq(i,j))
     &                        /max((1-albq(i,j)),1e-6*one)
            vinfor(nv+9)=vinfor(nv+9)+(1.0-albq(i,j))
     &                   *aire(i,j)*hgbq(i,j)
            vinfor(nv+11)=vinfor(nv+11)+(1.0-albq(i,j))
     &                   *aire(i,j)*hnbq(i,j)
            vinfor(nv+13)=vinfor(nv+13)
     &                   +(1.0-albq(i,j))*aire(i,j)*hgbq(i,j)*
     &                    (ug(i,j)*ug(i,j)+vg(i,j)*vg(i,j))
          endif
        enddo
      enddo
      vinfor(nv+13)=sqrt(vinfor(nv+13)/max(vinfor(nv+9),1e-5*one))
      nv=nv+13
      nv=nv+1

!     write(mouchard_id,*) nv, nvinfo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Calcul separe  Flux vers Nord / Flux vers Sud .
!----------+----------+----------+----------+----------+----------+----------+

!----------+----------+----------+----------+----------+----------+----------+
!---CALCULATING THE FRESHWATER CONTAINED IN ICE
!----------+----------+----------+----------+----------+----------+----------+
!---dxs1:The 'x' dimension of (ii,jj).
!---hnbq:The snowthickness of (ii,jj).
!---albq:The ice fraction of (ii,jj). % of the (ii,jj) containing ice.
!---hgbq:The ice thickness of (ii,jj).
!---0.9:Density of ice(900kgm-3)/density of seawater(1000kgm-3)

!---i) Fram Strait
!---'FRAG'
#if ( HRCLIO == 0 )
      j = 56
      ii1 = 107
      ii1 = 106
      ii2 = 108
#elif ( HRCLIO == 1 )
      j = 109
      ii1 = 209
      ii2 = 215
#endif
      vfram = 0.0
      vbord = 1.0+(1.0-bound)
      do i=ii1-1,ii2
          ip1    = i+1
          vfram  = (vg(i,j)*tmu(i,j,ks2)+vg(ip1,j)*tmu(ip1,j,ks2))/
     &             (max(tmu(i,j,ks2)+tmu(ip1,j,ks2),vbord))
          jj     = j-max(0,int(sign(one,vfram)))
          vinfor(nv) = vinfor(nv)+vfram*dxs1(i,j-1)*1e-6*
     &              (hnbq(i,jj)*0.33+hgbq(i,jj)*0.9)*(1.0-albq(i,jj))
      enddo
!-----
       nv=nv+1

#if ( FRAZER_ARCTIC == 1 )
!--- FSIF
!f      j = 55
!f      ii1 = 106
!f      ii2 = 107
!f       vbord = 1.0+(1.0-bound)
!f       do i=ii1,ii2
!f          ip1    = i+1
!f          vfram  = (vg(i,j)*tmu(i,j,ks2)+vg(ip1,j)*tmu(ip1,j,ks2))/
!f     &             (max(tmu(i,j,ks2)+tmu(ip1,j,ks2),vbord))
!f          jj     = j-max(0,int(sign(one,vfram)))
!f          vinfor(nv) = vinfor(nv)+vfram*dxs1(i,j-1)*1e-3*
!f     &              (hnbq(i,jj)*0.33+hgbq(i,jj)*0.9*(28.8/34.8))
!f     &              *(1.0-albq(i,jj))*31.104
!f       enddo
!f       write (*,*) "fsif1" , vinfor(nv)
!f       nv = nv + 1

      j = 55
      vfram = 0.0
      vbord = 1.0+(1.0-bound)
      vv = 0.0
      flux = 0.0
      do i=106,107
        vv  = (vg(i,j)*tmu(i,j,ks2))+(vg(i+1,j)*tmu(i+1,j,ks2))/
     &             (max(tmu(i,j,ks2)+tmu(i+1,j,ks2),vbord))
       flux = vv * dxs1(i,j-1) * 0.031536 *
     &           (((hnbq(i,j)*0.33) + (hgbq(i,j)*0.9)) * (28.8/34.8))
     &           * (1.0-albq(i,j))
          fsif = fsif + flux
      enddo
!       write(*,*), "fsif",fsif
       vinfor(nv) = fsif
       nv = nv + 1
#endif

!----------+----------+----------+----------+----------+----------+----------+
!---ii) The Barents Sea/Kara Boundary

!---'SPNG'

#if ( HRCLIO == 0 )
#if ( FRAZER_ARCTIC == 1 )
      ii1 = 109
      ii2 = 112
       vbord = 1.0+(1.0-bound)
      vfram = 0.0
      vv = 0.0
      flux = 0.0
#else
      j = 55
      ii1 = 109
      ii2 = 114
#endif
#elif ( HRCLIO == 1 )
      j = 108
      ii1 = 216
      ii2 = 225
#endif
#if ( FRAZER_ARCTIC == 1 )
       do i=ii1,ii2
#else
       do i=ii1-1,ii2
#endif
          ip1    = i+1
          vfram  = (vg(i,j)*tmu(i,j,ks2)+vg(ip1,j)*tmu(ip1,j,ks2))/
     &             (max(tmu(i,j,ks2)+tmu(ip1,j,ks2),vbord))
          jj     = j-max(0,int(sign(one,vfram)))
          vinfor(nv) = vinfor(nv)+vfram*dxs1(i,j-1)*1e-6*
     &              (hnbq(i,jj)*0.33+hgbq(i,jj)*0.9)*(1.0-albq(i,jj))
       enddo
!-----
       nv=nv+1
#if ( FRAZER_ARCTIC == 1 )
!---'KSIF'
!       j = 55
!       ii1 = 109
!       ii2 = 113
!      ii1 = 109
!      ii2 = 112
!       vfram = 0.0
!       vbord = 1.0+(1.0-bound)
!       do i=ii1,ii2
!        do i=109,112
!          ip1    = i+1
!          vfram  = (vg(i,j)*tmu(i,j,ks2)+vg(ip1,j)*tmu(ip1,j,ks2))/
!     &             (max(tmu(i,j,ks2)+tmu(ip1,j,ks2),vbord))
!          jj     = j-max(0,int(sign(one,vfram)))
!          vinfor(nv) = vinfor(nv)+vfram*dxs1(i,j-1)*1e-3*
!     &              (hnbq(i,jj)*0.33+hgbq(i,jj)*0.9*(28.8/34.8))
!     &              *(1.0-albq(i,jj))*31.104
!        enddo
!       enddo
!       nv = nv + 1

      j = 55
      vfram = 0.0
      vbord = 1.0+(1.0-bound)
      vv = 0.0
      flux = 0.0
      do i=109,113
        vv  = (vg(i,j)*tmu(i,j,ks2))+(vg(i+1,j)*tmu(i+1,j,ks2))/
     &             (max(tmu(i,j,ks2)+tmu(i+1,j,ks2),vbord))
       flux = vv * dxs1(i,j-1) * 0.031536 *
     &           (((hnbq(i,j)*0.33) + (hgbq(i,j)*0.9)) * (28.8/34.8))
     &           * (1.0-albq(i,j))
          ksif = fsif + flux
      enddo
       vinfor(nv) = ksif
       nv = nv + 1

#endif

!----------+----------+----------+----------+----------+----------+----------+
!---iii)The Bering Strait 'BERG'

       vber=vg(iberp,jberp)*tmu(iberp,jberp,ku2)/vbord*1e-6
       jj     = jberp-max(0,int(sign(one,vber)))
#if ( FRAZER_ARCTIC == 1 )
       vinfor(nv)=vber*dxs1(iberp,jberp-1)*(hgbq(iberp,jj)+
     &                           hnbq(iberp,jj))*(1.0-albq(iberp,jj))+
     &            vber*dxs1(iberp-1,jberp-1)*hgbq(iberp-1,jj)+
     &                          hnbq(iberp-1,jj)*(1.0-albq(iberp-1,jj))
#else
       vinfor(nv)=vber*dxs1(iberp,jberp-1)*hgbq(iberp,jj)
     &                                *(1.0-albq(iberp,jj))+
     &            vber*dxs1(iberp-1,jberp-1)*hgbq(iberp-1,jj)
     &                                *(1.0-albq(iberp-1,jj))
#endif
#if ( FRAZER_ARCTIC == 1 )
        nv = nv + 1
!---'BSIF'
      j = 65
      ii1 = 102
      ii2 = 102
      vfram = 0.0
       vbord = 1.0+(1.0-bound)
       do i=ii1,ii2
          ip1    = i+1
          vfram  = (vg(i,j)*tmu(i,j,ks2)+vg(ip1,j)*tmu(ip1,j,ks2))/
     &             (max(tmu(i,j,ks2)+tmu(ip1,j,ks2),vbord))
          jj     = j-max(0,int(sign(one,vfram)))
          vinfor(nv) = vinfor(nv)+vfram*dxs1(i,j-1)*1e-3*
     &              (hnbq(i,jj)*0.33+hgbq(i,jj)*0.9*(28.8/34.8))
     &              *(1.0-albq(i,jj))*31.104
       enddo
       nv = nv + 1
!
!      j = 65
!      vfram = 0.0
!      vbord = 1.0+(1.0-bound)
!      vv = 0.0
!      flux = 0.0
!      do i=102,103
!        vv  = (vg(i,j)*tmu(i,j,ks2))+(vg(i+1,j)*tmu(i+1,j,ks2))/
!     &             (max(tmu(i,j,ks2)+tmu(i+1,j,ks2),vbord))
!       flux = vv * dxs1(i,j-1) * 0.031536 *
!     &           (((hnbq(i,j)*0.33) + (hgbq(i,j)*0.9)) * (28.8/34.8))
!     &           * (1.0-albq(i,j))
!          bsif = bsif + flux
!       enddo
!       vinfor(nv) = bsif
!       nv = nv + 1


!----------+----------+----------+----------+----------+----------+----------+
!---iv) CAA Boundary 'CAIF'

!      j = 57
!      j = 53
!      ii1 = 101
!      ii2 = 103
!      vfram = 0.0
!       vbord = 1.0+(1.0-bound)
!       do i=ii1,ii2
!          ip1    = i+1
!          vfram  = (vg(i,j)*tmu(i,j,ks2)+vg(ip1,j)*tmu(ip1,j,ks2))/
!     &             (max(tmu(i,j,ks2)+tmu(ip1,j,ks2),vbord))
!          jj     = j-max(0,int(sign(one,vfram)))
!          vinfor(nv) = vinfor(nv)+vfram*dxs1(i,j-1)*1e-3*
!     &              (hnbq(i,jj)*0.33+hgbq(i,jj)*0.9*(28.8/34.8))
!     &              *(1.0-albq(i,jj))*31.104
!       enddo
!       nv = nv + 1

      j = 57
      vfram = 0.0
      vbord = 1.0+(1.0-bound)
      vv = 0.0
      flux = 0.0
      do i=101,103
        vv  = (vg(i,j)*tmu(i,j,ks2))+(vg(i+1,j)*tmu(i+1,j,ks2))/
     &             (max(tmu(i,j,ks2)+tmu(i+1,j,ks2),vbord))
       flux = vv * dxs1(i,j-1) * 0.031536 *
     &           (((hnbq(i,j)*0.33) + (hgbq(i,j)*0.9)) * (28.8/34.8))
     &           * (1.0-albq(i,j))
          caif = caif + flux
      enddo
       vinfor(nv) = caif
       nv = nv + 1

!----------+----------+----------+----------+----------+----------+----------+
!---v) Barents Sea Boundary 'BAIF'
!      i = 110
!      jj1 = 54
!      jj2 = 56
!      vfram = 0.0
!       vbord = 1.0+(1.0-bound)
!       do j=jj1,jj2
!          jp1    = j+1
!          vfram  = (vg(j,i)*tmu(j,i,ks2)+vg(jp1,i)*tmu(jp1,i,ks2))/
!     &             (max(tmu(j,i,ks2)+tmu(jp1,i,ks2),vbord))
!          ii     = i-max(0,int(sign(one,vfram)))
!          vinfor(nv) = vinfor(nv)+vfram*dxs1(i-1,j)*1e-3*
!     &              (hnbq(ii,j)*0.33+hgbq(ii,j)*0.9*(28.8/34.8))
!     &              *(1.0-albq(ii,j))*31.104
!       write(*,*), "baif", vinfor(nv)
!       enddo
!       nv = nv + 1

      i = 110
      vfram = 0.0
      vbord = 1.0+(1.0-bound)
      uu = 0.0
      vv = 0.0
      flux = 0.0
      do j=54,56
        uu  = (ug(i,j)*tmu(i,j,ks2))+(ug(i,j+1)*tmu(i,j+1,ks2))/
     &             (max(tmu(i,j,ks2)+tmu(i,j+1,ks2),vbord))
!       write(*,*), "uu", ug(i,j),i,j
       flux = uu * dxs1(i,j) * 0.031536 *
     &           (((hnbq(i,j)*0.33) + (hgbq(i,j)*0.9)) * (28.8/34.8))
     &           * (1.0-albq(i,j))
          baif = baif + flux
      enddo
!       write(*,*), "baif", baif
       vinfor(nv) = baif
       nv = nv + 1

!-----
!      j = 53
!      vfram = 0.0
!      vbord = 1.0+(1.0-bound)
!     vv = 0.0
!      uu = 0.0
!      flux = 0.0
!      do i=109,110
!        vv  = (vg(i,j)*tmu(i,j,ks2))+(vg(i+1,j)*tmu(i+1,j,ks2))/
!     &             (max(tmu(i,j,ks2)+tmu(i+1,j,ks2),vbord))
!        write(*,*), "vv", vg(i,j),i,j
!       flux = vv * dxs1(i,j) * 0.031536 *
!     &           (((hnbq(i,j)*0.33) + (hgbq(i,j)*0.9)) * (28.8/34.8))
!     &           * (1.0-albq(i,j))
!          baif2 = baif2 + flux
!       enddo
!       vinfor(nv) = baif2
!       nv = nv + 1


!----------+----------+----------+----------+----------+----------+----------+
#else
       nv=nv+1
#endif

!---Computation of thermal expansion of the ocean
       therma=0.0
       do k=ks1,ks2
        do j=js1,js2
         do i=is1(j),is2(j)
          therma=therma+dz(k)*ctmi(i,j,k,0)*b(i,j,k)
         enddo
        enddo
       enddo
       vinfor(nv)=(-1)*therma*zsurf/gpes
       thex=vinfor(nv)

!---iceshelf and icebergs
!-driess:LOVECLIM includes icebergs from ISM, excess of snow and
!-correction linked to the temperature anomalies
       nv=nv+1
       vinfor(nv)=toticesm
!mab       write(mouchard_id,*) 'informe',ficebergn,ficebergs
!mab       write(mouchard_id,*) 'informe2',ficebergn/xlg*360.,ficebergs/xlg*360.
       nv=nv+1
       vinfor(nv)=ficebergn/xlg*360.
#if ( ISM == 1 )
       if(flgism) vinfor(nv)=vinfor(nv)+vicbismn
#endif
       nv=nv+1
       vinfor(nv)=ficebergs/xlg*360.
#if ( ISM == 1 )
       if(flgism) vinfor(nv)=vinfor(nv)+vicbisms
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        if (nv.ne.nvinfo) then
          write(clio3_out_id,'(2(A,I4))') 
     &     'STOP in "informe" : compute nv=', nv,
     &     ' var. <> nvinfo=', nvinfo
          stop
        endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Accumulation before averaging                                   |
!-----------------------------------------------------------------------

      do nv=2,nvinfo
        vinfom(nv) = vinfom(nv) + vinfor(nv)
      enddo

      if (.not.flgout) return

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Ecriture sur le fichier type ".evol" :                          |
!-----------------------------------------------------------------------

!--Comptabilise le Nb d iter. cumulees dans "vinfom" :
      nninfo = abs(ninfo)
      nniter = numit - nstart + 1
!- facteur utilise pour Var. Non moyennee : zitsum(0)
      zitsum(0) = 1.
!- facteur utilise pour Var. Moyennee sur tts iter. : zitsum(1)
      if (nniter.eq.0) then
        zitsum(1) = 0.
      elseif (nniter.lt.nninfo) then
!- cas de moyenne sur moins de "ninfo" iter. (en debut de run) :
        zitsum(1) = 1. /  DFLOAT(nniter)
      else
!- cas de moyenne sur "ninfo" iter. (cas standard) :
        zitsum(1) = 1. / DFLOAT(nninfo)
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--Moyenne et transfert : vinfom -> vinfor ; reset "vinfom" :
!     do nv=2,nvinfo
!       vinfor(nv) = vinfom(nv) * zitsum(ktsum(nv))
!       vinfom(nv) = 0.0
!     enddo
      do nv=2,nm1n-1
        vinfor(nv) = vinfom(nv) * zitsum(ktsum(nv))
        vinfom(nv) = 0.0
      enddo
      do nv=nm2n,nvinfo
        vinfor(nv) = vinfom(nv) * zitsum(ktsum(nv))
        vinfom(nv) = 0.0
      enddo
      do nv=nm1n,nm2n-1
        vinfom(nv) = 0.0
      enddo


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!- Ecriture sur fichier :
      write(evolu_id,fmtw) (titvar(nv),vinfor(nv),nv=1,nvinfo)
#if ( ISM >= 2 )
      write(testevolu_id,fmtw) (titvar(nv),vinfor(nv),nv=1,nvinfo)
#endif
!- Fermeture du fichier :
      if(numit.eq.nferme) then
        write(*,*) ' evolu (in inforun.f), numit,nferme "',numit,nferme
        write(evolu_id,*)
        close(evolu_id)
#if ( ISM >= 2 )
        close(testevolu_id)
#endif
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  6 ) Traitement pour d'autres sorties et/ou d'autres routines.       |
!-----------------------------------------------------------------------

      if (koutpu.le.1 .or. mod(numit,nsav).ne.0) then
!--Remise a zero de fss, Energ. Aj.Conv & Freq.Aj.Conv :
        do k=1,kmax
         do j=1,jmax
          do i=1,imax
            fqajc(i,j,k) = 0.
          enddo
         enddo
        enddo

        do ns=0,nsmax
         do j=1,jmax
          do i=1,imax
           fss(i,j,ns) = 0.
          enddo
         enddo
        enddo
      else
!--Preparation des tableaux avant ecriture sur fichier de resultats :
        do k=ks1,ks2
         do j=js1,js2
          do i=is1(j),is2(j)
           fqajc(i,j,k) = zitsum(1) * fqajc(i,j,k)
          enddo
         enddo
        enddo
!- compute time average surface Fluxes :
        convfx = zitsum(1)
        do ns=0,nsmax
         if (ns.eq.1) convfx = zitsum(1) * ( dz(ks2) / dts(ks2) )
         do j=js1,js2
          do i=is1(j),is2(j)
           fss(i,j,ns) = convfx * fss(i,j,ns)
          enddo
         enddo
        enddo
        call raccord(fqajc(1,1,1), zero, kmax, 4)
        call raccord(fss(1,1,0), zero, 3, 0)
      endif

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine inforun -
      end subroutine inforun
