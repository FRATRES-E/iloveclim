!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:52 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:52 CET 2009

      SUBROUTINE scadhv(ns,ngroup)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! calcul des differentes contributions (Avection, Diffusion, Horiz., Vert.)
! de l'Eq. d'Evolution du scalaire "ns" .
! ngroup : type de regroupement : 1,3,5 Advec.Centr (Decentr=0) / Diffus. seule
!                                 2,4,6 Advec.Centr+Decentr / Advec.Decentr
!--
! "Cclp" <-> en accord avec "scali.f.clp"
! "Cstd" <-> en accord avec "scali.f.std"
! "Ccsp" <-> en accord avec "scali.f.csp"
!--
! alphah(1) intervient : terme de decentrement du a correction maillage z
!--
!  modif : 16/05/96

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip           

      use const_mod

      use para0_mod
      use para_mod
      use vareq_mod, only: haterm, hdterm, alphxo, alphyo

      use bloc0_mod
      use bloc_mod
      
      use newunit_clio_mod, only: clio3_out_id      
!! END_OF_USE_SECTION


!--parameter local :
      integer(kind=ip), parameter :: mtest = 0

!--variables locales :

!dmr [NOEQUI]      dimension scavat(imax,jmax,kmax+1), scavdt(imax,jmax,kmax+1)
!dmr [NOEQUI]      dimension z2alph(imax,jmax,2:kmax)
!dmr [NOEQUI]      equivalence ( scavat(1,1,1), w(1,1,1) ) !dmr  w(imax,jmax,kmax+1)
!dmr [NOEQUI]      equivalence ( scavdt(1,1,1), q2turb(1,1,1) )
!dmr [NOEQUI]      equivalence ( z2alph(1,1,2), fqajc(1,1,2) )
!dmr [NOEQUI]      dimension phixd(imax,jmax) , phiyd(imax,jmax)
!dmr [NOEQUI]      dimension phixa(imax,jmax) , phiya(imax,jmax)
!dmr [NOEQUI]      equivalence ( phixd(1,1) , phihhh(1,1,1) )
!dmr [NOEQUI]      equivalence ( phiyd(1,1) , phihhh(1,1,2) )
!dmr [NOEQUI]      equivalence ( phixa(1,1) , phihhh(1,1,3) )
!dmr [NOEQUI]      equivalence ( phiya(1,1) , phihhh(1,1,4) )
!dmr [NOEQUI]      dimension alphax(imax,jmax,kmax), alphay(imax,jmax,kmax)
!dmr [NOEQUI]      equivalence ( alphax(1,1,1) , alphxo(1,1,1) )
!dmr [NOEQUI]      equivalence ( alphay(1,1,1) , alphyo(1,1,1) )

!--variables locales :
      real(kind=dblp), dimension(ixjmax)   :: scaij
      integer(kind=ip), dimension(ixjmax)  :: nij
      real(kind=dblp), dimension(kmax)     :: ccwcor, scaref
      real(kind=dblp), dimension(imax,kmax):: wliss, varloc

!--- more locales

      real(kind=dblp) :: cc0w, cc1w, cc2w, ccstep, ccwid2, ssc2rf, sss
     &                 , sstm, sum2s, sumw, sumw1, sumw2, ww2, wwliss
     &                 , cc2x, cc2y, ccx, ccxdif, ccy, ccydif, u1cd2
     &                 , v2cd2
      integer(kind=ip):: i, j, jj, k, kk, l, n, ngroup, nnn, ns, nsum

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--choix de type de regroupement "ngroup" :
      if (ngroup.eq.0) then
        write(clio3_out_id,'(A)')
     &    ' Choix du type de regroupement (ngroup) ?'
        write(clio3_out_id,'(A)') ' 1,3 : Decentr. = 0 /'
     &               //' Advect.Centr / Diffus.Only'
        write(clio3_out_id,'(A)') ' 2,4 : Decentrement /'
     &               //' Advect.Centr+Decentr / Advec.Decentr'
        write(clio3_out_id,'(A)') ' Advect. Verticale :'
     &     //' 1,2 -> Contribution = W.grad(S) ; 3,4 -> Flux = W.S-S0'
        write(clio3_out_id,'(20X,A)')
     &    ' 5,6 -> Flux = W.(S-Sk) et calcule Sk'
        read(*,*) ngroup
      endif

      write(clio3_out_id,*) 'scadhv : debut , ns,ngroup =', ns, ngroup

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  0 ) preparation des coeff. intervenant dans les flux d'advection .  |
!-----------------------------------------------------------------------

!--initialisation (indispensable) :
      do n=1,4
       do j=1,jmax
        do i=1,imax
         phihhh(i,j,n) = 0.0
        enddo
       enddo
      enddo
      do k=1,kmax+1
       do j=1,jmax
        do i=1,imax
         phizzz(i,j,k,1) = 0.0
         phizzz(i,j,k,2) = 0.0
        enddo
       enddo
      enddo
      do k=2,kmax
       do j=1,jmax
        do i=1,imax
         fqajc(i,j,k) = 0.0
        enddo
       enddo
      enddo
!     do k=1,kmax+1
!      do i=1,imax
!       phiz(i,k) = 0.0
!      enddo
!     enddo

!--Mise en place des coeff (Decentrement Vertical) : (extrait de scali)
      cc0w = 0.5 * alphaz(ks1)
      do k=ks1+1,ks2
       cc1w = 0.5 * alphaz(k)
       cc2w = 0.5 * (1.0 - alphaz(k))
       do j=js1,js2
        do i=is1(j),is2(j)
          ccstep = ( tms(i-1,j,k) - tms(i-1,j,k-1) )
     &           + ( tms(i+1,j,k) - tms(i+1,j,k-1) )
     &           + ( tms(i,j-1,k) - tms(i,j-1,k-1) )
     &           + ( tms(i,j+1,k) - tms(i,j+1,k-1) )
          fqajc(i,j,k) = cc1w + min(cc2w , cc0w * ccstep)
        enddo
       enddo
      enddo

!--preparation des coeffs dependant de la verticale :
!     ccwid2 = 0.25
      ccwid2 = 0.25 * alphah(1)
      ssc2rf = 2.0 * scal0(ks1,ns)
      do k=1+ks1,ks2
        scaref(k) = scal0(k-1,ns) + scal0(k,ns)
!       ccwabs(k) = 0.5 * alphaz(k)
        ccwcor(k) = ccwid2 * unsdzw(k) * (dz(k) - dz(k-1))
     &                     * (1. - alphaz(k))
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Traitement des differents regroupements - Calcul de Sk .
!-----------------------------------------------------------------------

      if (mod(ngroup,2).eq.0) then
!--mise a zero de la Diffusivite :
        do k=1,kmax
         ahs(k) = 0.0
!        do j=1,jmax
!         do i=1,imax
!          avsdz(i,j,k) = 0.0
!         enddo
!        enddo
        enddo

      else
!-----
!--Supression du decentrement :
!       do k=ks1,ks2
!         ccwabs(k) = 0.0
!         ccwcor(k) = 0.0
!       enddo
      write(clio3_out_id,*)
     & 'scadhv : ngroup=1,3,5 => Annule Dec.H (alphax,y <- 0)'
        do k=1,kmax
         do j=1,jmax
          do i=1,imax
           alphxo(i,j,k) = 0.0
           alphyo(i,j,k) = 0.0
          enddo
         enddo
        enddo
!-----
      endif

      if (ngroup.ge.5) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Calcul de scaref(k) / Somme[(Flx)^2] mimimale :
      loop_170: do k=1+ks1,ks2
!- mtest = parameter : 0,1,2 selon le type de calcul de "scaref".
        if (mtest.ge.1) then
        sumw2 = 0.
        sum2s = 0.
        do j=js1,js2
         do i=is1(j),is2(j)
          if (mtest.eq.2) then
            ww2 = ctmi(i,j,k-1,0) * w(i,j,k) * w(i,j,k)
          else
            ww2 = ctmi(i,j,k-1,0) * abs(w(i,j,k))
          endif
          sumw2 = sumw2 + ww2
          sum2s = sum2s + ww2 * (scal(i,j,k-1,ns) + scal(i,j,k,ns))
         enddo
        enddo
        if (sumw2.gt.epsil) scaref(k) = sum2s / sumw2

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        else
!--Calcul de scaref(k) / Somme[ |Flx| ] mimimale :
!- classement de scal(k) par ordre croissant :
        n = 0
        sumw = 0.
!       loop_150outer: do j=js1,js2
        loop_150outer: do jj=0,js2-js1
         j = (js1 + jj/2)*mod(jj+1,2) + (js2 - (jj-1)/2 )*mod(jj,2)
         loop_150inner: do i=is1(j),is2(j)
          if (k.le.kfs(i,j)) cycle loop_150inner
          sumw = sumw + ctmi(i,j,k-1,0) * abs(w(i,j,k))
          n = n + 1
          nnn = i + imax*(j-1)
          sss  = scal(i,j,k-1,ns) + scal(i,j,k,ns)
          nij(n) = nnn
          scaij(n) = sss
          if (n.eq.1) cycle loop_150inner
          loop_145: do l=n-1,1,-1
            if (sss.ge.scaij(l)) then
              nij(l+1) = nnn
              scaij(l+1) = sss
              cycle loop_150inner
            endif
            nij(l+1) = nij(l)
            scaij(l+1) = scaij(l)
          enddo loop_145
          nij(1) = nnn
          scaij(1) = sss
         enddo loop_150inner
        enddo loop_150outer
        nsum = n
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- verification :
!       ss0 = -10.
!       do n=1,nsum
!         i = 1 + mod(nij(n)-1,imax)
!         j = 1 + (nij(n)-1)/imax
!         if (n.eq.1.or.n.eq.nsum) write(66,'(A,2I6,3I4,F8.3)')
!    &       'n,nij,i,j,k,Sn :', n,nij(n),i,j,k,0.5*scaij(n)
!         if (scaij(n).lt.ss0) then
!           write(clio3_out_id,*) 'scadhv : Pb1 dans le classement ! => arret'
!           write(clio3_out_id,*) 'n,nij(n),i,j,k :',n,nij(n),i,j,k
!           write(clio3_out_id,*) 'scai, scal(k-1), scal(k) :', 0.5*scaij(n),
!    &                  scal(i,j,k-1,ns), scal(i,j,k,ns)
!           stop
!         endif
!         ss0 = scaij(n)
!         sss = scal(i,j,k-1,ns) + scal(i,j,k,ns)
!         if (abs(sss-ss0).ge.epsil) then
!           write(clio3_out_id,*) 'scadhv : Pb2 dans le classement ! => arret'
!           write(clio3_out_id,*) 'n,nij(n),i,j,k :',n,nij(n),i,j,k
!           write(clio3_out_id,*) 'scai, scal(k-1), scal(k) :', 0.5*scaij(n),
!    &                  scal(i,j,k-1,ns), scal(i,j,k,ns)
!           stop
!         endif
!       enddo
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        sumw1 = 0.
        sumw2 = sumw * 0.5
!- cherche scaref(k) tel que Som.|w|(S<Sref) = Som.|w|(S>Sref)
        do n=1,nsum
          i = 1 + mod(nij(n)-1,imax)
          j = 1 + (nij(n)-1)/imax
          sumw1 = sumw1 + ctmi(i,j,k-1,0) * abs(w(i,j,k))
          if (sumw1.ge.sumw2) then
            scaref(k) = scaij(n)
            exit
          endif
        enddo
        endif
      enddo loop_170
      write( *,'(2A,6F8.3,10(/,10F8.3))') ' scadhv, Avect.Vert.,',
     &  ' scaref(k)=', (0.5*scaref(k),k=ks2,1+ks1,-1)
!-----
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Flux Verticaux advect. & diffus. -
!-----------------------------------------------------------------------

      loop_280: do j=js1,js2
!--debut de la boucle sur l'indice de latitude "j" .

!- 2.1) Lissage de la vitesse W :
!-----
      do k=ks1+1,ks2
       kk = k-1
       do i=is1(j),is2(j)
        wwliss = w(i-1,j,k) + w(i+1,j,k) + w(i,j-1,k) + w(i,j+1,k)
        sstm = epsil +
     &   tms(i-1,j,kk) + tms(i+1,j,kk) + tms(i,j-1,kk) + tms(i,j+1,kk)
        wwliss = 0.5 * ( wwliss / sstm + w(i,j,k) )
        wliss(i,k) = min(zero, max(w(i,j,k),wwliss) )
     &             + max(zero, min(w(i,j,k),wwliss) )
!       wliss(i,k) = min(zero, w(i,j,k) )
!    &               max(zero, min(w(i,j,k),wwliss) )
       enddo
      enddo

!- 2.2) Flux et terme d'Advection Verticale :
!-----
      if (ngroup.le.2) then
!--calcul termes verticaux advection centree :
      do k=ks1+1,ks2
       do i=is1(j),is2(j)
        phizzz(i,j,k,1) = tms(i,j,k-1) * w(i,j,k)
     &                * ( scal(i,j,k-1,ns) - scal(i,j,k,ns) )
       enddo
      enddo
!--Bilan termes verticaux advection centree :
      do k=ks1,ks2
       do i=is1(j),is2(j)
!--Attention : scavat & w equivalent ! => calcul de scavat ds varloc :
        varloc(i,k) = 0.5 * (phizzz(i,j,k,1) + phizzz(i,j,k+1,1))
       enddo
      enddo

      else
!-----
!--calcul Flux verticaux advection centree, avec Scal - Scal0(1)
      do k=ks1+1,ks2
       do i=is1(j),is2(j)
        phizzz(i,j,k,1) = tms(i,j,k-1) * w(i,j,k) * 0.5 *
     &        ( (scal(i,j,k-1,ns) + scal(i,j,k,ns)) - ssc2rf )
       enddo
      enddo
!--Bilan Flux verticaux advection :
      do k=ks1,ks2
       do i=is1(j),is2(j)
        varloc(i,k) = phizzz(i,j,k,1) - phizzz(i,j,k+1,1)
       enddo
      enddo
!--calcul Flux verticaux advection centree, avec Scal - Scalref(k)
      do k=ks1+1,ks2
       do i=is1(j),is2(j)
        phizzz(i,j,k,1) = tms(i,j,k-1) * w(i,j,k) * 0.5 *
     &        ( (scal(i,j,k-1,ns) + scal(i,j,k,ns)) - scaref(k) )
       enddo
      enddo
!-----
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!- 2.3) Flux et terme d'Advect.Decentr. & Diffusion Verticale :
!-----
!--calcul des flux verticaux advection decentree et diffusion.
      if (mod(ngroup,2).eq.1) then
        do k=ks1+1,ks2
         do i=is1(j),is2(j)
           phizzz(i,j,k,2) = tms(i,j,k-1) * avsdz(i,j,k)
     &             * ( scal(i,j,k-1,ns) - scal(i,j,k,ns) )
         enddo
        enddo
      else
!- calcul du decentrement :
        do k=ks1+1,ks2
         do i=is1(j),is2(j)
          phizzz(i,j,k,2) = tms(i,j,k-1) *
     &           ( fqajc(i,j,k)*abs(w(i,j,k)) + ccwcor(k)*wliss(i,k) )
!std & (avsdz(i,j,k) + ccwabs(k)*abs(w(i,j,k)) + ccwcor(k)*wliss(i,k) )
!clp &           max( zero , avsdz(i,j,k) + ccwcor(k)*wliss(i,k) )
!csp &           max( avsdz(i,j,k), ccwcor(k)*w(i,j,k) )
        enddo
       enddo
      if(ngroup.eq.2) then
!--Diffusivite equivalente -> mise dans phizzz(1) :
        do k=ks1+1,ks2
         do i=is1(j),is2(j)
          phizzz(i,j,k,1) = phizzz(i,j,k,2) + tms(i,j,k-1)*avsdz(i,j,k)
          phizzz(i,j,k,2) = phizzz(i,j,k,2)
     &                    * ( scal(i,j,k-1,ns) - scal(i,j,k,ns) )
        enddo
       enddo
      else
        do k=ks1+1,ks2
         do i=is1(j),is2(j)
          phizzz(i,j,k,2) = phizzz(i,j,k,2)
     &                    * ( scal(i,j,k-1,ns) - scal(i,j,k,ns) )
          phizzz(i,j,k,1) = phizzz(i,j,k,1) + phizzz(i,j,k,2)
         enddo
        enddo
      endif
      endif
!--fin du calcul des flux verticaux.

!--Bilan des Flux Verticaux Diffusion / Decentrement :
      do k=ks1,ks2
       do i=is1(j),is2(j)
!dmr [NOEQUI]
        w(i,j,k) = varloc(i,k)
        q2turb(i,j,k) = phizzz(i,j,k,2) - phizzz(i,j,k+1,2)
       enddo
      enddo

!--fin de la boucle sur l'indice de latitude "j".
      enddo loop_280

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Flux Horizontaux advect. & diffus -
!-----------------------------------------------------------------------

      loop_480: do k=ks1,ks2
!--Debut de la boucle sur l'indice de niveau "k" :

      ccx  = unsdx * dts(k)
      cc2x = ccx + ccx
      ccxdif = ahs(k) * unsdx
!--calcul des flux dans la direction x :
      do j=js1,js2
        do i=is1(j),is2(j)+1
          u1cd2 = 0.25 * cmy(i,j,1) * ( u(i,j,k) + u(i,j+1,k) )
          phihhh(i,j,1) = tms(i,j,k) * tms(i-1,j,k)
     &             * (alphxo(i,j,k) + smxy(i,j,1) * ccxdif)
     &             * (scal(i-1,j,k,ns) - scal(i,j,k,ns))
          phihhh(i,j,3) = u1cd2 * (scal(i-1,j,k,ns) - scal(i,j,k,ns))
        enddo
      enddo

!--calcul des flux dans la direction y :
      ccy  = unsdy * dts(k)
      cc2y = ccy + ccy
      ccydif = ahs(k) * unsdy
      do j=js1,1+js2
        do i=isf1(j),isf2(j)
          v2cd2 = 0.25 * cmx(i,j,2) * ( v(i,j,k) + v(i+1,j,k) )
          phihhh(i,j,2) = tms(i,j,k) * tms(i,j-1,k)
     &             * (alphyo(i,j,k) + cmxy(i,j,2) * ccydif)
     &             * (scal(i,j-1,k,ns) - scal(i,j,k,ns))
          phihhh(i,j,4) = v2cd2 * (scal(i,j-1,k,ns) - scal(i,j,k,ns))
        enddo
      enddo

!--bilan des Flux Diff & Advect. Decentree / Advect. Centree :
      ccx  = ccx * dz(k) / dts(k)
      ccy  = ccy * dz(k) / dts(k)
      do j=js1,js2
        do i=is1(j),is2(j)
          hdterm(i,j,k) = smxy(i,j,0) *
     &                  ( ccx * (phihhh(i,j,1) - phihhh(i+1,j,1))
     &                  + ccy * (phihhh(i,j,2) - phihhh(i,j+1,2)) )
          haterm(i,j,k) = smxy(i,j,0) *
     &                  ( ccx * (phihhh(i,j,3) + phihhh(i+1,j,3))
     &                  + ccy * (phihhh(i,j,4) + phihhh(i,j+1,4)) )
        enddo
      enddo

!--Fin de la boucle sur l'indice de niveau "k".
      enddo loop_480

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Regroupement de termes ; restaure les tx. Decentrement .        |
!-----------------------------------------------------------------------

      if (mod(ngroup,2).eq.0) then

!--Somme Advect.Centre + Advect.Decentree :
        do k=1,kmax
         do j=1,jmax
          do i=1,imax
!dmr [NOEQUI]
           w(i,j,k) = w(i,j,k) + q2turb(i,j,k)
           haterm(i,j,k) = haterm(i,j,k) + hdterm(i,j,k)
          enddo
         enddo
        enddo

!--Restaure les Taux de decentrement alphax & alphay :
      do k=ks1,ks2
        do j=js1,js2
         do i=is1(j),is2(j)+1
          u1cd2 = 0.25 * cmy(i,j,1) * ( u(i,j,k) + u(i,j+1,k) )
!         alphax(i,j,k) = tms(i,j,k) * tms(i-1,j,k)
!    &             * alphax(i,j,k) / max(epsil, abs(u1cd2) )
          alphxo(i,j,k) = alphxo(i,j,k) / max(epsil, abs(u1cd2) )
         enddo
        enddo
        do j=js1,1+js2
         do i=isf1(j),isf2(j)
          v2cd2 = 0.25 * cmx(i,j,2) * ( v(i,j,k) + v(i+1,j,k) )
!         alphay(i,j,k) = tms(i,j,k) * tms(i,j-1,k)
!    &             * alphay(i,j,k) / max(epsil, abs(v2cd2) )
          alphyo(i,j,k) = alphyo(i,j,k) / max(epsil, abs(v2cd2) )
         enddo
        enddo
      enddo

      endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      write(clio3_out_id,*) 'scadhv : fin'

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine scadhv -
      end
