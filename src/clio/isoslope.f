!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:49 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:49 CET 2009

      SUBROUTINE isoslope
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!======================================================================
!  Computation of isopycnal slope : SLOPE CUT-OFF Limiter (Cox, 1987)
!======================================================================
!- Attention au signe !  drho,y,z = - d(rho) / dx,y,z ; idem pour drox,y,z1,2,4
!- modif : "epsil2" est ajoute dans "drhoz"
!  modif : 02/06/99

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip           

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use isoslope_mod
      use reper_mod

!!     USE OMP_LIB
      use mchd99_mod, only: nn99

      use newunit_clio_mod, only: clio3_out_id, mouchard_id
       
!! END_OF_USE_SECTION


      integer(kind=ip), parameter :: ntabmx = 100 
      real(kind=dblp), dimension(imax,jmax,kmax), save :: drhox, drhoy
     &, drhoz
      real(kind=dblp), dimension(imax,jmax), save :: rossby
      real(kind=dblp), dimension(-ntabmx:ntabmx), save :: ztabf1
      real(kind=dblp), dimension(-1:ntabmx), save :: ztabf2
      real(kind=dblp), save :: slopc, unsdds, unsddr
      integer(kind=ip), save:: ndd1mx, ndd1mn, ndd2mx
      
!    & slopc, slopd, refcur, radmn, radmx, domdds, ddslop, ddrati,
!- slopd ... ddrati => facultatif              

!--- loads of locales
      real(kind=dblp) :: ddrati, ddslop, domdds, drox4, droy4, droz1
     &                 , droz2, epsil2, radmn, radmx, refcur, slopd
     &                 , xx, zz, ccxy, drr, dss, rrn, slp4x, slp4y
     &                 , slp4z, ssn, ztap1, ztap2
     
      integer(kind=ip):: i, ie, is, iuo, j, k, km0, km1, kp0, kp1, nn


!- n99=2 => ecritures auxiliaires sur fichier "mouchard", unit=99
c~       common / mchd99 / nn99

      epsil2 = epsil * epsil

      if (numit.eq.nstart) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- debut du traitement specifique 1ere it.
!-----------------------------------------------------------------------

!- initialisation des variables "isopyc + GM" en common :
      do k=1,kmax
       do j=1,jmax
        do i=1,imax
          c1x(i,j,k) = zero
          c2y(i,j,k) = zero
          c4x(i,j,k) = zero
          c4y(i,j,k) = zero
          c4xgm(i,j,k) = zero
          c4ygm(i,j,k) = zero
          drhox(i,j,k) = zero
          drhoy(i,j,k) = zero
          drhoz(i,j,k) = zero
          uiso(i,j,k) = zero
          viso(i,j,k) = zero
!dmr --- [DELETED] unused in isoslope         dscalz(i,j,k) = zero
        enddo
       enddo
      enddo

      do k=1,kmax+1
       do j=1,jmax
        do i=1,imax
          wiso(i,j,k) = zero
!dmr --- [DELETED] unused in isoslope          fisoz(i,j,k) = zero
        enddo
       enddo
      enddo

      if(ai(kmax).eq.zero.and.aitd(kmax).eq.zero) return

       write(clio3_out_id,'(2A,2F6.4)')
     &   ' *** isoslope ; iso,gm90 *** :',
     &   ' SlopMx_ISO,GM =', slopemax(kmax), slopmgm(kmax)

!- computation of new masks
      do k=1,kmax
       do j=1,jmax
        do i=1,imax
          if (i.ne.1) ttm1(i,j,k) = tms(i,j,k)*tms(i-1,j,k)
          if (j.ne.1) then
            ttm2(i,j,k) = tms(i,j,k)*tms(i,j-1,k)
          else
            ttm2(i,j,k) = zero
          endif
        enddo
        ttm1(1,j,k) = ttm1(ims2,j,k)
       enddo
      enddo

!- setup index jsdom1,jsdom2(i) to cover the computation domain :
      do i=1,imax
        jsdom1(i) = jmax
        jsdom2(i) = 1
      enddo
      do j=js1,js2
       do i=is1(j),is2(j)
         if (tms(i,j,ks2).eq.one) then
           jsdom1(i) = min(jsdom1(i),j)
           jsdom2(i) = max(jsdom2(i),j)
         endif
       enddo
      enddo
!- write index jsdom1,jsdom2 on file "mouchard"
      if (nn99.eq.2) then
        write(mouchard_id,'(A)') ' isoslope ; i / jsdom1 / jsdom2 :'
        do is=1,imax,18
          ie = min(imax,is+17)
          write(mouchard_id,'(A,18I4)') '    i = ', (i,i=is,ie)
          write(mouchard_id,'(A,18I4)') ' jsdom1=', (jsdom1(i),i=is,ie)
          write(mouchard_id,'(A,18I4)') ' jsdom2=', (jsdom2(i),i=is,ie)
        enddo
        write(mouchard_id,*)
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- limiteur F1 & F2 : calcul de "Rossby_Radius" & tabulation de F1 & F2 :
!-----------------------------------------------------------------------

!- mise en place des coeffs (cf papier Large et al) :
!  slopc=Sc ; slopd=Sd ; refcur="c"=1rst.Barocl.Wave Speed ;
!  [radmn,radmx] = interval for Rossby radius (R)
!-  discretis.F1 : domaine : -domdds*Sd < S < domdds*Sd ; Pas=ddslop*Sd
!-  discretis.F2 : domaine : 0 < ratio=|z|/R*S < 1 ; Pas=ddrati
! NB: slopc, slopd, ... /Smax in file "run.param"
        slopc = slopmgm(kmax-1)
        slopd = slopmgm(kmax-2)
        refcur = slopmgm(kmax-3)
        radmn  = slopmgm(kmax-4)
        radmx  = slopmgm(kmax-5)
        domdds = slopmgm(kmax-6)
        ddslop = slopmgm(kmax-7)
        ddrati = slopmgm(kmax-8)
!- si Increment(ddslop or ddrati) Nul => pas de limiteur (F1=1 or F2=1)
        if (ddslop.gt.zero) then
          slopd  = max(slopd,epsil)
          domdds = max(domdds,ddslop)
          unsdds = 1.0 / (slopd*ddslop)
          ndd1mx = nint(domdds/ddslop)
          ndd1mn = -ndd1mx-1
        else
          slopd  = 1.
          domdds = 1.
          ddslop = 0.
          unsdds = 0.
          ndd1mx = 0
          ndd1mn = 0
          ztabf1(0) = 1.
        endif
        if (ddrati.gt.zero) then
          unsddr = 1.0 / ddrati
          ndd2mx = nint(unsddr)
        else
          unsddr = 0.
          ndd2mx = -1
          ztabf2(0) = 1.
        endif
!---------
        if (ndd1mx+1.gt.ntabmx .or. ndd2mx+1.gt.ntabmx) then
          write(clio3_out_id,'(2A)') 'STOP in "isoslope" :',
     &               ' array ztabf1,ztabf2 sous-dimensione !'
          write(clio3_out_id,'(2(A,2I6))') ' Max Index=', ndd1mx+1, ndd2mx+1,
     &                          ' > Dim=', ntabmx
          stop
        endif
!--------
        zz = 0.5 / tanh(domdds)
        do nn=1+ndd1mn,ndd1mx
          xx = ddslop*nn
          ztabf1(nn) = 0.5 + zz*tanh(-xx)
        enddo
        ztabf1(-ndd1mx-1) = ztabf1(-ndd1mx)
        ztabf1( ndd1mx+1) = ztabf1( ndd1mx)
!--------
        do nn=0,ndd2mx
          xx = ddrati*nn
          ztabf2(nn) = 0.5 + 0.5*sin( pi*(xx-0.5) )
        enddo
        ztabf2(-1) = ztabf2(0)
        ztabf2(ndd2mx+1) = ztabf2(ndd2mx)
!--------
        do j=1,jmax
         do i=1,imax
           rossby(i,j) = 0.5*refcur/max(abs(fs2cor(i,j)),epsil)
           rossby(i,j) = min(radmx,max(radmn,rossby(i,j)))
         enddo
        enddo

!--------
!- Sortie sur fichier "iso_tab.chk" des fonctions F1 & F2 :
!        include 'iso_slp_tab.inc'
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin du traitement specifique 1ere it.
      endif

      if (ai(kmax).eq.zero.and.aitd(kmax).eq.zero) return

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Compute (minus)GRADIENTS OF BUYANCY B .
!-----------------------------------------------------------------------

      do k=ks1,ks2
        do j=js1,js2+1
         do i=isf1(j),isf2(j)+1
          drhox(i,j,k) = ttm1(i,j,k) *
     &                 (b(i-1,j,k)-b(i,j,k))*unsdx*smx(i,j,1)
          drhoy(i,j,k) = ttm2(i,j,k) *
     &                 (b(i,j-1,k)-b(i,j,k))*unsdy*smy(i,j,2)
         enddo
        enddo
        do j=1,jmax
         do i=1,imax
          drhoz(i,j,k) = max( bvf(i,j,k), epsil2 )
         enddo
        enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Start  external loop on all levels, index "k".
c~ !$DIR SPP LOOP_PARALLEL
c~ !$DIR SPP LOOP_PRIVATE(km1,km0,kp1,kp0,j,i)
c~ !$DIR SPP LOOP_PRIVATE(droz1,droz2,drox4,droy4, slp4x,slp4y,slp4z)
c~ !$DIR SPP LOOP_PRIVATE(nn,ssn,dss,ztap1, rrn,drr,ztap2, ccxy)
c~ !$OMP PARALLEL
c~ !$OMP DO PRIVATE (k,j,i,km1,km0,kp1,kp0)
c~ !$OMP&   PRIVATE (droz1,droz2,drox4,droy4,slp4x,slp4y,slp4z)
c~ !$OMP&   PRIVATE (nn,ssn,dss,ztap1,rrn,drr,ztap2, ccxy)
c~ !$OMP&  SCHEDULE(static)
      do k=ks1,ks2
!-----
      km1 = max(k-1, ks1)
      km0 = km1 + 1
      kp1 = min(k+1, ks2)
      kp0 = kp1 - 1

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Compute Isopycnal Slope at point 1 (for Diffus. Isopyc).        |
!-----------------------------------------------------------------------

      do j=js1,js2
       do i=is1(j),is2(j)+1
         droz1 = ( drhoz(i,j,km0) + drhoz(i-1,j,kp1)
     &           + drhoz(i,j,kp1) + drhoz(i-1,j,km0) ) /
     &           ( tms(i,j,km1) + tms(i-1,j,kp0)
     &           + tms(i,j,kp0) + tms(i-1,j,km1) + epsil )
         c1x(i,j,k) = ttm1(i,j,k) * drhox(i,j,k)
     &              / max( droz1, abs(drhox(i,j,k)/slopemax(k)) )
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Compute Isopycnal Slope at point 2 (for Diffus. Isopyc).        |
!-----------------------------------------------------------------------

      do j=js1,js2+1
       do i=isf1(j),isf2(j)
         droz2 = ( drhoz(i,j,km0) + drhoz(i,j-1,kp1)
     &           + drhoz(i,j,kp1) + drhoz(i,j-1,km0) ) /
     &           ( tms(i,j,km1) + tms(i,j-1,kp0)
     &           + tms(i,j,kp0) + tms(i,j-1,km1) + epsil )
         c2y(i,j,k) = ttm2(i,j,k) * drhoy(i,j,k)
     &              / max( droz2, abs(drhoy(i,j,k)/slopemax(k)) )
       enddo
      enddo

      if (k.ne.ks1) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Compute Isopycnal Slope at point 4 (for Diffus. Isopyc & GM scheme).
!-----------------------------------------------------------------------

      do j=js1,js2
       do i=is1(j),is2(j)
         drox4 = ( drhox(i,j,k) + drhox(i+1,j,km1)
     &           + drhox(i,j,km1) + drhox(i+1,j,k) ) /
     &           ( ttm1(i,j,k) + ttm1(i+1,j,km1)
     &           + ttm1(i,j,km1) + ttm1(i+1,j,k) + epsil )
         droy4 = ( drhoy(i,j,k) + drhoy(i,j+1,km1)
     &           + drhoy(i,j,km1) + drhoy(i,j+1,k) ) /
     &           ( ttm2(i,j,k) + ttm2(i,j+1,km1)
     &           + ttm2(i,j,km1) + ttm2(i,j+1,k) + epsil )

!- slope for Diffus. Isopyc. :
         c4x(i,j,k) = tms(i,j,k-1) * drox4
     &              / max( drhoz(i,j,k), abs(drox4/slopmgm(ks2)) )
         c4y(i,j,k) = tms(i,j,k-1) * droy4
     &              / max( drhoz(i,j,k), abs(droy4/slopmgm(ks2)) )

!- slope for GM scheme :
         slp4x = drox4 / max( drhoz(i,j,k), abs(drox4/slopmgm(ks2)) )
         slp4y = droy4 / max( drhoz(i,j,k), abs(droy4/slopmgm(ks2)) )

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! limit slope*Ai <- Large et al, 1997, JPO, Vol.27, p 2445+2446
! test : only applied on A_ITD
! --> ATTENTION : 1) put limit of A_ITD on slope_x,y  --> c4x & c4y
! and 2) use only slopmgm(kmax), for other k<kmax :
!       => slopc,slopd,refcur,Radmn,Radmx, domdds, ddslop, ddrati
!----------------------------------------------------------------------
!     unsdds = 1.0 / (slopd*ddslop)
!     ndd1mx = nint(domdds/ddslop)
!     ndd1mn = -ndd1mx-1
!     unsddr = 1.0 / ddrati
!     ndd2mx = nint(unsddr)
!---------
!- limiteur F1 (slp4z discretise, par pas de "ddslop") :
          slp4z = abs(slp4x) + abs(slp4y) + epsil2
          ssn = (slp4z - slopc) * unsdds
          nn = nint(ssn-0.5)
          dss = ssn - nn
          nn = min(ndd1mx,max(ndd1mn,nn))
          ztap1 = ztabf1(nn) + dss * (ztabf1(nn+1) - ztabf1(nn))
!- limiteur F2 (ratio discretise, par pas de "ddrati") :
!         ratio = -zw(k) / (rossby(i,j) * slp4z)
!         rrn = ratio * unsddr
          rrn = -zw(k)*unsddr / (rossby(i,j) * slp4z)
          rrn=min(rrn,1D9)
          nn = nint(rrn-0.5)
          drr = rrn - nn
          nn = min(ndd2mx,max(-1,nn))
          ztap2 = ztabf2(nn) + drr * (ztabf2(nn+1) - ztabf2(nn))
!---------
          ccxy = tms(i,j,k-1)*ztap1*ztap2
          c4xgm(i,j,k) = ccxy*slp4x
          c4ygm(i,j,k) = ccxy*slp4y
!---------
!- Sortie sur fichier "iso_slp.chk" de la colonne "icheck,jcheck" :
!         include 'iso_slp_out.inc'
!---------
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
       end if ! k.ne.ks1
!- End of external loop on all levels, index "k".
      enddo

c~ !$OMP END DO
c~ !$OMP END PARALLEL

      do k=ks1,ks2 
       if (k.ne.ks1) then

!- CYCLIC Raccord [c4x] for GM90
        if (ltest.ge.1) then
          do j=jcl1,jcl2
            c4xgm(1,j,k) = c4xgm(ims2,j,k)
            c4xgm(imax,j,k) = c4xgm(ims1,j,k)
          enddo
        endif
!- BERING Raccord [c4y] for GM90 (attention minus sign)
        if (ltest.eq.3) then
          c4ygm(iberpm,jberp,k) = -c4ygm(ibera, jberam,k)
          c4ygm(iberp, jberp,k) = -c4ygm(iberam,jberam,k)
          c4ygm(iberam,jbera,k) = -c4ygm(iberp, jberpm,k)
          c4ygm(ibera, jbera,k) = -c4ygm(iberpm,jberpm,k)
        endif

      end if ! k.ne.ks1
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- End of external loop on all levels, index "k".
      enddo

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine isoslope -
      end
