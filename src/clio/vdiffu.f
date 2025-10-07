!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:57 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:57 CET 2009

      SUBROUTINE vdiffu(nn99)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Calcul des coefficients de DIFFUsion Verticale (divise par dz) : avsdz, avudz
!  Viscosite & Diffusivite Verticale : fonction de la Freq. de Br.Vais. (bvf)
!   suivant la formulation de Pacanowski et Philander (1981)
! Modification : Diffusivite V. Simplifie = avkb + avk0 / (1+5Ri)^3
!   Les seuils portent sur les termes en 1/(1+5Ri)
!--
! Diffusivite Ocean Profond : AvsN2 * (N2)^q  (q = -0.6 environ)
!  Inactif si kavsmx=0 ; Interfaces concernes : k=2,kavsmx
!  bornes (sur Avs) globales (k=1) et par niveau : avsmn & avsmx (lu sur fich)
!  => equivalent sur N2 : bvfmn & bvfmx
!  bvfmn(1) & bvfmx(1) definissent la discretication de la fct Avs(N2)
!   par pas de "dn2" suivant N2 : avsn(nav),nav=0,navmax
!--
!  modif : 11/05/99

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod

      use newunit_clio_mod, only: clio3_out_id, mouchard_id

!! END_OF_USE_SECTION


      integer(kind=ip), parameter :: navmax = 10000

!--variables locales conservees d'un appel a l'autre :
      real(kind=dblp), save       :: epsil2, ref0n2, ref0m2, uv2bnd
     &               , ccudm, cc1zm, kk0ri, unsdn2, bvfmix

      real(kind=dblp), dimension(kmax), save:: avu0ri, avs0ri, avu1ri
     &               , avs1ri, avs2mx, bvfmn, bvfmx

      real(kind=dblp), dimension(imax, jmax), save:: unszek


      real(kind=dblp), dimension(0:navmax+1), save :: avsn

!--variables locales :
      real(kind=dblp), dimension(imax, jmax):: riums, bvfloc
      real(kind=dblp), dimension(kmax)      :: avsmn, avsmx

!--valeurs intervenant dans la formule de Pac&Phi (certaines deja ds le common)
      real(kind=dblp) :: alpha
      data alpha / 5.0d0     /


      real(kind=dblp) :: avsmix, cc2zm, ccdz, ccm2, ccm2n, cczzm, deno
     &                 , dn2, tm4u, ud, unsm2, unsqav, unszz4, uvzbnd
     &                 , uz, vd, vz, xx, zzekm, bvfavs, demi, rifac
     &                 , unsden, xav

      integer(kind=ip):: i, j, k, nav, nn, nn99

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  0 ) 1ere Iter, Determine 1er niv. a calculer ; Initialisation .     |
!-----------------------------------------------------------------------

      if (numit.le.nstart) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Formulation Avs(N2) : mise en place des Coeff.

      do k=1,kmax
        avsmn(k) = 1.
        avsmx(k) = 1.
      enddo

!- Lecture des parametres pour la formulation Avs(N2) <- dans defcst
      if_kavsmx_positive: if (kavsmx > 0) then
      kavsmx = min(kavsmx,kmax)

!- Mise en place des coeff. & Verification :
      unsqav = one / qavs
      do k=1,kavsmx
        avsmn(k) = ccfmn * avsmn(k)
        avsmx(k) = ccfmx * avsmx(k)
        if ( avsmn(k).lt.avsmn(1) .or. avsmn(k).gt.avsmx(k) .or.
     &       avsmx(k).gt.avsmx(1) ) then
          write(clio3_out_id,*)
     &      'ARRET : vdiffu, coeff. avsmn,avsmx False', k
          stop
        endif
        bvfmn(k) = (avsmx(k)/avsn2)**unsqav
        bvfmx(k) = (avsmn(k)/avsn2)**unsqav
      enddo
      dn2 = (bvfmx(1) - bvfmn(1)) / DFLOAT(navmax)
      if ( dn2.lt.epsil*epsil ) then
          write(clio3_out_id,*)
     &      'ARRET : vdiffu, N2 Discr. step too small', dn2
          stop
      endif
      unsdn2 = one / dn2

!-  Mise en place de la Fonct. Tabulee : Avs(N2)
      do nav=0,navmax+1
        avsn(nav) = avsn2 * ( bvfmn(1) + dn2 * DFLOAT(nav) )**qavs
      enddo

      if (nn99.eq.2) then
!--Impression sur fichier mouchard :
        write(mouchard_id,'(2A,I4,1P2E16.8)') 'vdiffu, Avs(N2) :',
     &     ' kavsmx,qavs,avsn2 =', kavsmx, qavs, avsn2
        write(99,*) 'navmax(=N), dn2 =', navmax, dn2
        nn = navmax / 2
        write(mouchard_id,'(A,1P4E14.6)') ' avsn(0,1,N/2,N) :',
     &       avsn(0), avsn(1), avsn(nn), avsn(navmax)
        write(mouchard_id,
     &   '(A,(1P6E12.5))') 'avsmn :', (avsmn(k),k=1,kavsmx)
        write(mouchard_id,
     &   '(A,(1P6E12.5))') 'avsmx :', (avsmx(k),k=1,kavsmx)
        write(mouchard_id,
     &   '(A,(1P6E12.5))') 'bvfmn :', (bvfmn(k),k=1,kavsmx)
        write(mouchard_id,
     &   '(A,(1P6E12.5))') 'bvfmx :', (bvfmx(k),k=1,kavsmx)
        write(mouchard_id,*)
      endif

      endif if_kavsmx_positive
!- fin du bloc de definition des Coeff. pour formulation Avs(N2).
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      epsil2 = epsil * epsil
!- valeur typique de N2 : 1.E-5 (s-2)
      ref0n2 = 1.d-5 * epsil

      avsmix = avnub(1)
      bvfmix = avnu0(1)
      ref0m2 = avnub(1)
      uvzbnd = avnu0(1)
      ccudm = avkb(1)
      zzekm = avk0(1)

      xx = uvzbnd * 100.0
        write(clio3_out_id,'(2A,F6.2,A1,I4,A1)')
!    &   ' vdiffu : Pac.&Phil. Ris,Riu separes'
!    &   ' vdiffu : Pac.&Phil. avu=f(min{N2}) '
!    &   ' vdiffu : Pac.&Phil. Ris = Moy.Riu  '
     &   ' vdiffu : Pac.&Phil. avu=f(min.N2) & Ris = Moy.Riu'
     &  ,', MinM2 :', ccudm, ',',  nint(zzekm), 'm'
!    &  ,'<-MinM2 :', ccudm, ',',  nint(zzekm), 'm'
!       write(clio3_out_id,'(2A,1PE8.1,A,I2,2A,0PF4.1,A,I4,A)')
!    &   ' vdiffu : Pac&Phil, avu=f(min.N2), Ris=4.Riu'
!    &   ' vdiffu : Pac&Phil, avu=f(min.N2), Ris=4.RiU'
!    &  ,', MinM2:', ref0m2, ',', nint(xx), 'cm/s'
!    &  ,',', ccudm, ',', nint(zzekm), 'm'
        if (lstab.eq.0) write(66,'(2(A,1PE10.3))')
     &   '          AVS +', avsmix, ' si bvf <', bvfmix
        if (kavsmx.ne.0) write(66,'(3X,A,I4,1P2E11.3)')
     &   ' Avs = AvsN2*(N2)^q ; kmx,AvsN2,q =', kavsmx, avsn2, qavs

!--1er Niveau a calculer :
      kk0ri = ks2
      avu0ri(ks1) = 0.0
      avs0ri(ks1) = 0.0
      do k=ks1+1,ks2
        if (kk0ri.eq.ks2 .and. (avk0(k).ne.zero .or. avnu0(k).ne.zero))
     &      kk0ri = k - 1
        avu0ri(k) = unsdzw(k) * avnub(k)
        avu1ri(k) = unsdzw(k) * avnu0(k)
        avs0ri(k) = unsdzw(k) * avkb(k)
        avs1ri(k) = unsdzw(k) * avk0(k)
        avs2mx(k) = unsdzw(k) * avsmix * 0.5
        if (k.le.kavsmx) avs0ri(k) = 0.
      enddo

      uv2bnd = uvzbnd * uvzbnd
      if (ccudm.ge.epsil) then
!--profondeur pour le calcul du cisaillement minimun :
        cc1zm = ccudm * ccudm * 4.
        cc2zm = cc1zm / ( epsil + zzekm * zzekm )
        cc2zm = cc2zm / (omega * omega)
        do j=1,jmax
         do i=1,imax
          unszek(i,j) = cc2zm * fs2cor(i,j) * fs2cor(i,j)
         enddo
        enddo
      endif

!- Fin du traitement specifique de la 1ere Iter.
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

!--Debut de la boucle externe sur l'indice de niveau k :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,j,unszz4,uz,vz,unsm2,cczzm,ud,vd,ccm2,ccm2n)
!$DIR SPP LOOP_PRIVATE(riums,ccdz,tm4u,deno,rifac,unsden)
!$DIR SPP LOOP_PRIVATE(nav,bvfloc,bvfavs,xav)
      do k=ks1+1,ks2
!-----

      if (k.le.kk0ri) then
!--Niveaux profond : Avu,Avs independent de Ri :
        do j=js1,js2
         do i=ims1,ims2
          avudz(i,j,k) = avu0ri(k)
          avsdz(i,j,k) = avs0ri(k)
         enddo
        enddo

      else
!--Niveaux peu profond : Calcul de Avu,Avs en fonction de Ri :
        unszz4 = unsdzw(k) * unsdzw(k) * 4.

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Ridchardson Number at velocity point .                          |
!-----------------------------------------------------------------------
!  avudz <- Ri used for momentum,  [=Min(N2)/M2] or [=Moy(N2)/M2]
!  riums <- Ri used for scalar  ,  [=Moy(N2)/M2]

!--Initialisation :
      do j=1,jmax
       do i=1,imax
        riums(i,j) = ref0n2
       enddo
      enddo

      if (ccudm.lt.epsil) then
!--Evaluation directe du cissaillement :
        do j=ju1,ju2
         do i=iu1(j),iu2(j)
          uz = u(i,j,k) - u(i,j,k-1)
          vz = v(i,j,k) - v(i,j,k-1)
          unsm2 = (tmu(i,j,k-1) + epsil2)
     &          / (unszz4 * (uz*uz + vz*vz) + epsil)
!- Ridchardson Number at velocity point, for scalar   eq. (riums) :
          riums(i,j) = ( (bvf(i-1,j-1,k) + bvf(i,j,k))
     &                 + (bvf(i-1,j,k) + bvf(i,j-1,k)) ) * unsm2
!- Ridchardson Number at velocity point, for momentum eq. (avudz) :
          avudz(i,j,k) = min(bvf(i-1,j-1,k), bvf(i,j,k),
     &                       bvf(i-1,j,k), bvf(i,j-1,k)) * 4. * unsm2
!         avudz(i,j,k) = riums(i,j)
!-
         enddo
        enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      else
!--Evaluation du cissaillement : Minore par (|V| / max(Dekman,Z))
        cczzm = cc1zm / (zw(k) * zw(k))
        do j=ju1,ju2
         do i=iu1(j),iu2(j)
          uz = u(i,j,k) - u(i,j,k-1)
          vz = v(i,j,k) - v(i,j,k-1)
          ud = u(i,j,k) + u(i,j,k-1)
          vd = v(i,j,k) + v(i,j,k-1)
          ccm2  = (uz*uz + vz*vz) * unszz4 + epsil
          ccm2n = (ud*ud + vd*vd) * min(cczzm, unszek(i,j))
          unsm2 = (tmu(i,j,k-1) + epsil2) / max(ccm2, ccm2n)
!- Ridchardson Number at velocity point, for scalar   eq. (riums) :
          riums(i,j) = ( (bvf(i-1,j-1,k) + bvf(i,j,k))
     &               + (bvf(i-1,j,k) + bvf(i,j-1,k)) ) * unsm2
!- Ridchardson Number at velocity point, for momentum eq. (avudz) :
!         avudz(i,j,k) = ( (bvf(i-1,j-1,k) + bvf(i,j,k))
!    &                   + (bvf(i-1,j,k) + bvf(i,j-1,k)) )
          avudz(i,j,k) = min( bvf(i-1,j-1,k), bvf(i,j,k),
     &                        bvf(i-1,j,k), bvf(i,j-1,k) ) * 4.0
!    &                 * (tmu(i,j,k-1) + epsil2) / ccm2
     &                 * unsm2
!         avudz(i,j,k) = riums(i,j)
!-
         enddo
        enddo
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Ridchardson Number at scalar point :                            |
!-----------------------------------------------------------------------

!--calcul (dans avsdz) du nombre de Ridchardson au point scalaire :
!     unszz  = unsdzw(k) * unsdzw(k)
!C    unsz16 = unsdzw(k) * unsdzw(k) * 0.0625
!-  Weat Point Only :
!     do j=js1,js2
!      do i=is1(j),is2(j)
!       uz = ( (u(i,j,k)-u(i,j,k-1)) + (u(i+1,j+1,k)-u(i+1,j+1,k-1)) )
!    &     + ( (u(i,j+1,k)-u(i,j+1,k-1)) + (u(i+1,j,k)-u(i+1,j,k-1)) )
!       vz = ( (v(i,j,k)-v(i,j,k-1)) + (v(i+1,j+1,k)-v(i+1,j+1,k-1)) )
!    &     + ( (v(i,j+1,k)-v(i,j+1,k-1)) + (v(i+1,j,k)-v(i+1,j,k-1)) )
!       tm4u = tmu(i,j,k)+tmu(i+1,j+1,k)+tmu(i,j+1,k)+tmu(i+1,j,k)
!       cm2 = epsil + unszz*(uz*uz + vz*vz)
!C      avsdz(i,j,k) = bvf(i,j,k) / cm2
!       avsdz(i,j,k) = bvf(i,j,k) * tm4u * tm4u / cm2
!      enddo
!     enddo

!--raccord cyclique et autre (riums) :
      if (ltest.ge.1) then
        do j=jcl1,jcl2
          riums(ims1-1,j) = riums(ims2,j)
          riums(ims2+1,j) = riums(ims1,j)
        enddo
      endif
      if (ltest.eq.3) then
        riums(ibera, jbera) = riums(iberp, jberp)
      endif
!     call raccord(riums(1,1), 0., 1, 7)

!-  Moyenne des 4 Ri(pt.Vitesse) - Weat Point Only :
      ccdz = unsdzw(k) * unsdzw(k)
      do j=js1,js2
       do i=is1(j),is2(j)
        tm4u = tmu(i,j,k-1) + tmu(i+1,j+1,k-1)
     &       + tmu(i,j+1,k-1) + tmu(i+1,j,k-1) + epsil2
        avsdz(i,j,k) = ( (riums(i+1,j+1) + riums(i,j))
     &                 + (riums(i+1,j) + riums(i,j+1)) ) / tm4u
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Coefficients de diffusion Vert. pour les Vitesses .             |
!-----------------------------------------------------------------------

!--Coeff. Diffus. Vert. pour les Scalaires : avudz = Viscosite Vert. / dz
      do j=ju1,ju2
       do i=iu1(j),iu2(j)
        deno =  1.0 + alpha*avudz(i,j,k)
        rifac =  1. / (epsil + deno*deno )
        avudz(i,j,k) = avu0ri(k) + avu1ri(k) * min(rifumx, rifac)
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Coefficients de diffusion Vert. pour les Scalaires .            |
!-----------------------------------------------------------------------

!--Coeff. Diffus. Vert. pour les Scalaires : avsdz = Diffusivite Vert. / dz
!- (Pac.Phil. simplifie) :
      do j=js1,js2
       do i=is1(j),is2(j)
        unsden = 1.0 / ( abs(1.0 + alpha*avsdz(i,j,k)) + epsil )
        rifac = unsden * unsden * unsden
        avsdz(i,j,k) = avs0ri(k) + avs1ri(k) * min( rifsmx, rifac )
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Fin du traitement separe selon k, Avu,Avs independant/dependant de Ri
      endif

      if (k.le.kavsmx) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  6 ) Diffusion Verticale (Deep Ocean) Fonct. de (N2)^q .
!-----------------------------------------------------------------------

!--Initialisation :
      do j=1,jmax
       do i=1,imax
        bvfloc(i,j) = bvfmn(1)
       enddo
      enddo

!- Applique le minorant/Majorant :
      do j=js1,js2
       do i=is1(j)-1,is2(j)+1
        bvfloc(i,j) = max( bvfmn(k), min( bvfmx(k), bvf(i,j,k) ))
       enddo
      enddo

!- Raccord cyclique & autres <- inutile.

!--Calcul de Avs(N2) et incorpore dans avsdz :
      do j=js1,js2
       do i=is1(j),is2(j)
!- Max de N2 sur 9 points :
        bvfavs = max( bvfloc(i,j),
     & bvfloc(i-1,j-1),bvfloc(i+1,j-1),bvfloc(i-1,j+1),bvfloc(i+1,j+1),
     &   bvfloc(i-1,j), bvfloc(i+1,j), bvfloc(i,j-1), bvfloc(i,j+1) )
!- calcul approche de Avs(N2) :
        xav = (bvfavs - bvfmn(1)) * unsdn2
        nav = int(xav)
        xav = xav - DFLOAT(nav)
        avsdz(i,j,k) = avsdz(i,j,k)
     &    + unsdzw(k)*( (1.-xav)*avsn(nav) + xav*avsn(nav+1) )
       enddo
      enddo

!--Fin du traitement Niveaux profonds (k=<kavsmx) : Avs Fonct de (N2)^q
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Fin de la boucle externe sur l'indice de niveau k .
      enddo

      if (lstab.eq.0) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  8 ) Ajust.Conv. par augmentation de la Diffusion Verticale .
!-----------------------------------------------------------------------

!- Ajoute avsmix(/dzw) si bvf < bvfmix :
      demi = 0.5
      do k=1+ks1,ks2
       do j=js1,js2
        do i=is1(j),is2(j)
         avsdz(i,j,k) = avsdz(i,j,k)
     &         + ( avs2mx(k) + sign(avs2mx(k), bvfmix-bvf(i,j,k)) )
         fqajc(i,j,k) = fqajc(i,j,k)
     &         + ( demi + sign(demi, bvfmix-bvf(i,j,k)) )
        enddo
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine vdiffu -
      end
