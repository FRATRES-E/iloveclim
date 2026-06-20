!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009
#include "choixcomposantes.h"

!==============================================================================
! SUBROUTINE defgrid(nflag)  — VERSION RESTRUCTURÉE
!
! RÔLE : Orchestrateur d'initialisation de la grille océanique CLIO.
!        Délègue chaque groupe de calcul à une sous-routine dédiée.
!
! ARGUMENT D'ENTRÉE :
!   nflag (entier) :
!     1 = initialisation standard (premier appel)
!     2 = idem + écriture de contrôle dans le fichier "mouchard"
!     3 = mise à jour de la bathymétrie (paléo-simulations)
!
! STRUCTURE D'APPEL (ordre de dépendance) :
!
!   [init]  Mise à zéro de tous les champs (section 1 originale)
!   [init]  Calcul des constantes de résolution et indices géographiques (sect. 2)
!   [G1]    compute_grid_coords()      — coordonnées lon/lat et métriques locales
!   [G2]    compute_metric_coeffs()    — cmx/cmy/smx/smy/cmxy/smxy/fs2cor/covrai
!   [G3a]   compute_ice_metrics()      — wght/akappa/alambd/bkappa/dxs1../area
!   [G3b]   read_bathymetry_and_vscale() — z/zw/dz.../kbath
!   [G4a]   build_ocean_mask()         — tms/tmu/hs/hu/kniv/kfs/kfu/is1../iu1..
!   [G4b]   build_domain_indices()     — isf1../iuf1../ju1/ju2/imu1/imu2/bering
!   [G5a]   compute_integrals()        — ctmi/aire/volz/n1-4coin/kslp..
!   [sect8] write_mouchard_defgrid()   — sortie de contrôle (nflag==2 seulement)
!   [G5b]   compute_surface_physics()  — zindfa/dfhu/dfhv/tauc/reslum/msks/msku
!   [G5c]   define_ice_shelves()       — tmics/areicebn/areicebs
!
! POUR CHANGER DE TYPE DE GRILLE :
!   Remplacer UNIQUEMENT compute_grid_coords() et compute_metric_coeffs().
!   Toutes les autres sous-routines (G3a..G5c) sont agnostiques à la grille.
!
! TABLEAUX LOCAUX INTERMÉDIAIRES (passés entre G1/G2/G3a) :
!   dx1_loc(imax,jmax)  : pas angulaire adimensionné en x (rad) — sortie G1
!   dx2_loc(imax,jmax)  : pas angulaire adimensionné en y (rad) — sortie G1
!   h1_loc(imax,jmax)   : coefficient métrique x au pt S (m)    — sortie G1
!   h2_loc(imax,jmax)   : coefficient métrique y au pt S (m)    — sortie G1
!   d1d2_loc(imax,jmax) : d(h1)/dy au pt S (m/rad)             — sortie G1
!   d2d1_loc(imax,jmax) : d(h2)/dx au pt S (m/rad)             — sortie G1
!   unszw_loc(kmax)     : 1/zw(k) — sortie G3b, entrée G4a
!   zbath_loc(kmax)     : profondeur brute des centres (bath.om) — sortie G3b
!   dzbath_loc(kmax)    : épaisseur brute des niveaux (bath.om)  — sortie G3b
!==============================================================================

      SUBROUTINE defgrid(nflag)

!------------------------------------------------------------------------------
! Modules utilisés (identiques à la version originale)
!------------------------------------------------------------------------------
#if ( SHELFMELT == 1 )
      use varsClio2ism_mod,  only: tmsism
#endif
      use global_constants_mod, only: dblp=>dp, ip
      use const_mod,  only: cstmin, cstmax, radian, rterre, omega
     &                    , pi, zero, gpes, one, rho0, cpo, unsrt
      use para0_mod,  only: nsmax, ltest, jsepar, ncomax, nlpmax
     &                    , imax, jmax, kmax
      use para_mod,   only: nbsmax
      use bloc0_mod,  only: xslon, yslat, phifu, phifv, kniv, unshu
     &                    , hu, hs, hux, huy, xulon, yulat, xsedg, ysedg
     &                    , angle, yslatp, xslonp, xuedg, yuedg
     &                    , unsdzw, dzw, z, zw, scs, q2turb, scpme, phivs
     &                    , ims1, ims2, tms, tmu, ks1, ks2, js1, js2
     &                    , scal, scalr, ijudl, ijsdl, jdl1, jdl2, jeq
     &                    , ku1, ku2, dx, dy, unsdx, unsdy, jberpm
     &                    , jcl1, jcl2, iberp, ibera
     &                    , jbera, jberp, iberam, jberam, xulonp, yulatp
     &                    , is1, is2, iu1, iu2, iuf1, iuf2
     &                    , ju1, ju2, imu1, imu2, kfs, kfu, dtb, zurfow
     &                    , dts, dtu, isf1, isf2, afilt, slopmgm, aitd
     &                    , slopemax, avnu0, avnub, avv, ahh, ai, ahs
     &                    , alphgr, alphah, rifsmx, rifumx, msks, alphmi
     &                    , alphaz, msku, algrmn, avk0, avkb, ahe, ahu
     &                    , bering, unsdz, dz, iberpm, fss, phifs, phiss
     &                    , rappes, phimnx, bvf, avsdz, avudz, fqajc
     &                    , rappel, w, scal0, uns2dx, uns2dy
      use bloc_mod,   only: cmx, cmy, smx, smy, cmxy, smxy
     &                    , n1coin, n2coin, n3coin, n4coin
     &                    , i1coin, i2coin, i3coin, i4coin
     &                    , kslp, lslp, ijslp, nxslp, aire
     &                    , unsvol, ctmi, cmxdy, cmydx, covrai
     &                    , xang1, xang2, fub, fvb, phihhh, phizzz
     &                    , nrap, q, q2tmin, fs2cor, kfond, nxyslp
      use ice_mod,    only: tmics, tauc, reslum, areicebs, iicebern1
     &                    , iicebern2, jicebern1, jicebern2, areicebn
     &                    , iicebers1, iicebers2, jicebers1, jicebers2
      use dynami_mod, only: alambd, akappa, bkappa, area, dfhu, dfhv
     &                    , dxs1, gridsz, dxc1, dxc2, zindfa, ren
     &                    , dxs2, zfn, uvdif, wght, npo1i0, ipo1i0
     &                    , jpo1i0, zepsd1
      use reper_mod,  only: kbath, kbath2, iezon, iszon, jehsf, iehsf
     &                    , ishsf, jshsf, tithsf, nvhsf, ndhsf, xwpoln
     &                    , icheck, jcheck, kcheck, dxwi, dywj
     &                    , dlong, dlat, ylat1, xlon1, xaj1, yai1
     &                    , dxaj, dyai, jsep, kbath1, jnorth, xwi1, ywj1
      use coord_mod,       only: zlatt, zlont
      use newunit_clio_mod, only: clio3_out_id, mouchard_id

      implicit none

!------------------------------------------------------------------------------
! Argument d'entrée
!------------------------------------------------------------------------------
      integer(ip), intent(in) :: nflag

!------------------------------------------------------------------------------
! Variables locales propres à l'orchestrateur
!------------------------------------------------------------------------------
      integer(ip) :: kmaxp1, i, j, k, ns, nn, npt0v
      integer(ip) :: ipt0v(10), jpt0v(10)

      ! Tableaux intermédiaires G1 -> G2 -> G3a
      ! (non stockés dans les modules car propres à l'initialisation)
      real(dblp), dimension(imax,jmax) :: dx1_loc, dx2_loc
      real(dblp), dimension(imax,jmax) :: h1_loc,  h2_loc
      real(dblp), dimension(imax,jmax) :: d1d2_loc, d2d1_loc

      ! Tableaux intermédiaires G3b -> G4a
      real(dblp), dimension(kmax)   :: zbath_loc, dzbath_loc, unszw_loc

      ! Tableau de séparation grille AA, spécifique à la résolution CL15/CL30
      integer(ip) :: jsep2(83:imax)
      data jsep2 / 49, 50, 50, 50, 50, 49, 47, 46, 46, 47, 47, 47,
     &             56, 56, 56, 56, 63, 64, 65, 65, 65, 65, 65, 65,
     &             65, 64, 64, 61, 61, 60, 59, 46, 46, 46, 46, 46,
     &             46, 46, 46, 29 /

      ! Sauvegarde des masques (paléo-simulations avec bathymétrie variable)
#if ( BATHY >= 2 )
      real(dblp), dimension(imax,jmax,kmax) :: tms_old, tmu_old
#endif

      !dmr [ADD] Ajouts de variables oubliees par Claude
      integer(ip):: nxslp0, nyslp0 
      real(dblp) :: volz

      ! Variable locale pour le coefficient de Béring avant calcul géostrophique
      real(dblp) :: bering_init

      kmaxp1 = kmax + 1

      WRITE(*,*) "DEBUG nflag = =", nflag
!==============================================================================
! ÉTAPE 0 : INITIALISATION DES TABLEAUX GLOBAUX
!           Mise à zéro de tous les champs avant calcul.
!           Identique à la section 1 du code original.
!==============================================================================

!-- Bornes i : initialisées à "vide" (is1=imax, is2=1)
      do j=1,jmax
        is1(j)=imax  ; is2(j)=1
        iu1(j)=imax  ; iu2(j)=1
        iuf1(j)=imax ; iuf2(j)=1
        isf1(j)=imax ; isf2(j)=1
      enddo

!-- Profils verticaux
      do k=1,kmax
        z(k)=0. 
        zw(k)=0. 
        dz(k)=0. 
        dzw(k)=0.
        unsdz(k)=0.
        unsdzw(k)=0.
      enddo
      z(kmax+1)=0.
      zw(kmax+1)=0.
      dzw(kmax+1)=0.
      unsdzw(kmax+1)=0.

!-- Champs 2D de profondeur et masques
      do j=1,jmax
        do i=1,imax
          hs(i,j)=0. 
          hu(i,j)=0.
          hux(i,j)=0.
          huy(i,j)=0.
          unshu(i,j)=0.
          kniv(i,j,-1)=kmaxp1  
          kniv(i,j,0)=kmaxp1
          kniv(i,j,1)=kmaxp1
          xang1(i,j)=0. 
          xang2(i,j)=1.
          if (nflag.ne.3) then
            phifu(i,j)=0. 
            phifv(i,j)=0.
          endif
        enddo
      enddo

!-- Initialisation conditionnelle (pas lors des updates paléo)
      if (nflag.ne.3) then

!-- Flux et termes de relaxation surface (traceurs)
        do ns=0,nsmax
          do j=1,jmax
            do i=1,imax
              fss(i,j,ns)=0. 
              phifs(i,j,ns)=0.
              phiss(i,j,ns)=0.
              rappes(i,j,ns)=0.
              phimnx(i,j,0,ns)=cstmin 
              phimnx(i,j,1,ns)=cstmax
            enddo
          enddo
        enddo

!-- Composantes de diffusion isopycnale
        do nn=1,6
          do j=1,jmax
            do i=1,imax
              phihhh(i,j,nn)=0.
            enddo
          enddo
        enddo

!-- Champs 3D
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              fub(i,j,k)=0. ; fvb(i,j,k)=0. ; q(i,j,k)=0.
              bvf(i,j,k)=0. ; avsdz(i,j,k)=0. ; avudz(i,j,k)=0.
              fqajc(i,j,k)=0.
              tms(i,j,k)=0. ; tmu(i,j,k)=0.
              rappel(i,j,k)=0.
            enddo
          enddo
        enddo

!-- Vitesse verticale et traceurs aux interfaces
        do k=1,kmax+1
          do j=1,jmax
            do i=1,imax
              w(i,j,k)=0.
              phizzz(i,j,k,1)=0. ; phizzz(i,j,k,2)=0.
            enddo
          enddo
        enddo

!-- Traceurs scalaires : initialisation au profil de référence scal0
        do ns=1,nsmax
          do k=1,kmax
            nrap(k,ns)=0
            do j=1,jmax
              do i=1,imax
                scalr(i,j,k,ns)=0.
                scal(i,j,k,ns)=scal0(k,ns)
                phivs(i,j,k,ns)=0.
              enddo
            enddo
          enddo
        enddo

!-- Énergie cinétique turbulente
        do k=1,kmax
          do j=1,jmax
            do i=1,imax
              q2turb(i,j,k)=q2tmin
            enddo
          enddo
        enddo

!-- Terme PME (Précipitation - Évaporation) pour les scalaires de surface
        do ns=1,nsmax
          do j=1,jmax
            do i=1,imax
              scs(i,j,ns)=scpme(ns)
            enddo
          enddo
        enddo

      else  ! nflag == 3 : update bathymétrie paléo
#if ( BATHY >= 2 )
        tms_old(:,:,:)=tms(:,:,:)
        tmu_old(:,:,:)=tmu(:,:,:)
        tms(:,:,:)=0.
        tmu(:,:,:)=0.
#endif
      endif

!==============================================================================
! ÉTAPE 1 : CONSTANTES DE RÉSOLUTION ET INDICES GÉOGRAPHIQUES
!           Calcul des bornes actives, résolutions dx/dy et indices de détroits.
!==============================================================================
      ims1=2 ; ims2=imax-1
      js1=2  ; js2=jmax-1
      ks1=1  ; ks2=kmax
      ku1=ks1 ; ku2=ks2

!-- Indices géographiques clés (jeq, dateline, Béring) via geogra
      call geogra(js1, js2, jeq, jdl1, jdl2, ijsdl, ijudl, iberp, ibera)

!-- Bornes de la correspondance cyclique E<->W
      if (ltest.lt.1) then
        jcl1=js2 ; jcl2=js1   ! bassin fermé : boucle vide (jcl1>jcl2)
      else
        jcl1=js1 ; jcl2=jeq   ! grille monde : correspondance jusqu'à l'équateur
      endif

!-- Taille de maille en mètres
      if (ltest.lt.0) then
        dx=1000. ; dy=1000.   ! grille cartésienne (cas test)
      else
        dx=dlong*radian*rterre ; dy=dlat*radian*rterre
      endif
      unsdx=1./dx ; unsdy=1./dy ; uns2dx=.5/dx ; uns2dy=.5/dy

!-- Indices Béring décalés de -1 (pour les opérateurs aux demi-pas)
      jberpm=jberp-1 ; iberpm=iberp-1
      jberam=jbera-1 ; iberam=ibera-1

!-- Îles "sans surface de vitesse" (pts scalaires isolés, ex: Spitzberg)
      ipt0v(1)=109 ; jpt0v(1)=55 ; npt0v=0  ! désactivé (npt0v=0)

!-- Points avec vitesse océanique mais sans vitesse glace (Mer de Barents)
      npo1i0=2
      ipo1i0(1)=102 ; jpo1i0(1)=56
      ipo1i0(2)=111 ; jpo1i0(2)=56

!-- Sauvegarde du coefficient Béring avant recalcul géostrophique en G4b
      bering_init=bering

!==============================================================================
! [G1] COORDONNÉES GÉOGRAPHIQUES ET MÉTRIQUES LOCALES
!      *** SEUL BLOC À REMPLACER POUR UN AUTRE TYPE DE GRILLE ***
!
!      SORTIES VERS G2 et G3a via tableaux locaux :
!        dx1_loc, dx2_loc  : pas angulaires adimensionnés (rad)
!        h1_loc, h2_loc    : coefficients métriques (m)
!        d1d2_loc, d2d1_loc : dérivées des métriques (m/rad)
!      SORTIES DIRECTEMENT DANS LES MODULES :
!        xslon, yslat, xulon, yulat, zlatt, zlont
!        xsedg, ysedg, xuedg, yuedg
!        xslonp, yslatp, xulonp, yulatp
!        angle, xang1, xang2
!==============================================================================
      call compute_grid_coords(
     &  ltest, jsepar, jeq, jsep, jsep2,
     &  dlong, dlat, xlon1, ylat1,
     &  dxaj, dyai, xaj1, yai1,
     &  radian, rterre, unsrt, omega, pi, zepsd1,
     &  imax, jmax,
     &  xslon, yslat, xulon, yulat,
     &  zlatt, zlont,
     &  xsedg, ysedg, xuedg, yuedg,
     &  xslonp, yslatp, xulonp, yulatp,
     &  angle, xang1, xang2,
     &  dx1_loc, dx2_loc,
     &  h1_loc, h2_loc, d1d2_loc, d2d1_loc)

!==============================================================================
! [G2] COEFFICIENTS MÉTRIQUES
!      *** À REMPLACER CONJOINTEMENT AVEC G1 POUR UN AUTRE TYPE DE GRILLE ***
!
!      ENTRÉES : tableaux locaux dx1_loc, dx2_loc, h1_loc, h2_loc,
!                d1d2_loc, d2d1_loc (sorties de G1)
!      SORTIES DANS LES MODULES :
!        cmx, cmy, smx, smy, cmxy, smxy (4 positions chacun)
!        cmxdy, cmydx
!        fs2cor, covrai
!==============================================================================
      call compute_metric_coeffs(
     &  ltest, jsepar, jeq, jsep,
     &  iberp, ibera, iberpm, iberam,
     &  jberp, jbera, jberpm, jberam,
     &  omega, unsrt, radian,
     &  dlat, dlong, dyai, dxaj,
     &  ylat1, xlon1, yai1, xaj1,
     &  imax, jmax,
     &  h1_loc, h2_loc, d1d2_loc, d2d1_loc,
     &  dx1_loc, dx2_loc,
     &  cmx, cmy, smx, smy, cmxy, smxy,
     &  cmxdy, cmydx, fs2cor, covrai)

!==============================================================================
! [G3a] MÉTRIQUES POUR LA DYNAMIQUE DE LA GLACE DE MER
!       Agnostique à la grille.
!
!       ENTRÉES : tableaux locaux h1_loc, h2_loc, d1d2_loc, d2d1_loc,
!                 dx1_loc, dx2_loc + fs2cor (module)
!       SORTIES DANS LES MODULES :
!         wght, akappa, alambd, bkappa
!         dxs1, dxs2, dxc1, dxc2, area, zfn
!==============================================================================
      call compute_ice_metrics(
     &  imax, jmax,
     &  h1_loc, h2_loc, d1d2_loc, d2d1_loc,
     &  dx1_loc, dx2_loc, fs2cor,
     &  ltest, iberp, ibera, jberp, jbera,
     &  wght, akappa, alambd, bkappa,
     &  dxs1, dxs2, dxc1, dxc2, area, zfn)

!==============================================================================
! [G3b] LECTURE DE LA BATHYMÉTRIE ET CALCUL DE L'ÉCHELLE VERTICALE
!       Agnostique à la grille.
!
!       ENTRÉES  : fichiers bath.om / bath_new.txt / bath_XXXXX.txt
!       SORTIES DANS LES MODULES :
!         z, zw, dz, dzw, unsdz, unsdzw
!         kbath (reper_mod), kbath1, kbath2
!       SORTIES LOCALES (transmises à G4a) :
!         unszw_loc, zbath_loc, dzbath_loc
!==============================================================================
      call read_bathymetry_and_vscale(
     &  nflag, kfond, ks1, ks2, kmax,
     &  imax, jmax, ltest, jeq, jsep,
     &  jcl1, jcl2, ims1, ims2,
     &  zero,
     &  z, zw, dz, dzw, unsdz, unsdzw,
     &  zbath_loc, dzbath_loc, unszw_loc,
     &  kbath, kbath1, kbath2)

!==============================================================================
! [G4a] CONSTRUCTION DES MASQUES TERRE/MER
!       Agnostique à la grille.
!
!       ENTRÉES  : kbath, unszw_loc (sorties G3b)
!                  zw (module), ipt0v, jpt0v, npt0v
!       SORTIES DANS LES MODULES :
!         tms, tmu, msks, msku
!         hs, hu, hux (partiel), huy (partiel), unshu
!         kniv, kfs, kfu
!         is1, is2, iu1, iu2
!==============================================================================
      call build_ocean_mask(
     &  nflag, ltest,
     &  ks1, ks2, ims1, ims2, js1, js2, jcl1, jcl2,
     &  imax, jmax, kmax, kmaxp1,
     &  iberp, ibera, iberpm, iberam,
     &  jberp, jbera, jberpm, jberam,
     &  ipt0v, jpt0v, npt0v,
     &  zw, unszw_loc, kbath,
     &  tms, tmu, 
     &  hs, hu, hux, huy, unshu,
     &  kniv, kfs, kfu,
     &  is1, is2, iu1, iu2)

!==============================================================================
! [G4b] CALCUL DES INDICES DE DOMAINE ET DU COEFFICIENT DE BÉRING
!       Agnostique à la grille.
!
!       ENTRÉES  : tms, tmu, hs, hu, is1, is2, iu1, iu2 (sorties G4a)
!                  fs2cor (module), gpes, bering_init
!       SORTIES DANS LES MODULES :
!         isf1, isf2, iuf1, iuf2
!         ju1, ju2, imu1, imu2
!         hux, huy (finalisation)
!         bering (coefficient géostrophique final)
!==============================================================================
      call build_domain_indices(
     &  ltest, ks1, ks2, ims1, ims2, js1, js2, jcl1, jcl2,
     &  imax, jmax, kmax,
     &  iberp, ibera, jberp, jbera,
     &  tms, tmu, hs, hu, is1, is2, iu1, iu2,
     &  isf1, isf2, iuf1, iuf2,
     &  ju1, ju2, imu1, imu2,
     &  fs2cor, gpes, bering_init, zero,
     &  hux, huy, bering,
     &  clio3_out_id)

!==============================================================================
! [G5a] GRANDEURS INTÉGRALES : SURFACES, VOLUMES, COINS, PENTES
!       Agnostique à la grille.
!
!       ENTRÉES  : tms, tmu, cmxy, hs, kfs (modules)
!                  is1, is2, isf1, isf2, iu1, iu2, ju1, ju2 (modules)
!       SORTIES DANS LES MODULES :
!         ctmi, aire, volz, unsvol, zurfow
!         n1-4coin, i1-4coin
!         kslp, lslp, ijslp, nxslp, nxyslp
!==============================================================================
      call compute_integrals(
     &  ks1, ks2, ims1, ims2, js1, js2,
     &  imax, jmax, kmax,
     &  is1, is2, isf1, isf2, iu1, iu2, ju1, ju2,
     &  tms, tmu, cmxy, hs, kfs, kfond,
     &  dx, dy, ncomax, nlpmax, one,
     &  ctmi, aire, volz, unsvol, zurfow,
     &  n1coin, n2coin, n3coin, n4coin,
     &  i1coin, i2coin, i3coin, i4coin,
     &  kslp, lslp, ijslp, nxslp, nxyslp,
     &  nxslp0, nyslp0,
     &  clio3_out_id)

!==============================================================================
! ÉTAPE 8 : ÉCRITURE DE CONTRÔLE DANS LE FICHIER "MOUCHARD" (si nflag==2)
!           Inline dans defgrid — accès à toutes les variables de modules.
!==============================================================================

      WRITE(*,*) "DEBUG :: reslum == ", sum(reslum), __LINE__

!==============================================================================
! [G5b] PHYSIQUE DE SURFACE : GLACE, ABSORPTION SOLAIRE, FICHIERS DE SORTIE
!       Agnostique à la grille.
!
!       ENTRÉES  : tms, covrai, zw, dts, unsdz, kfs (modules)
!                  hs, hu, angle, msku (modules)
!                  zlatt, zlont (coord_mod)
!       SORTIES DANS LES MODULES :
!         zindfa, dfhu, dfhv
!         tauc, reslum
!         msks, msku (finalisés), hs, hu (halos E/W), angle (annulé sur terre)
!==============================================================================
      call compute_surface_physics(
     &  nflag, imax, jmax, kmax,
     &  ks2, is1, is2, js1, js2, kfs,
     &  tms, covrai, zw, dts, unsdz, rho0, cpo,
     &  radian, uvdif, ren, gridsz,
     &  hs, hu, angle, msku, msks,
     &  zlatt, zlont,
     &  mouchard_id,
     &  zindfa, dfhu, dfhv, tauc, reslum)

      WRITE(*,*) "DEBUG :: reslum == ", sum(reslum), __LINE__
!==============================================================================
! [G5c] PLATEFORMES DE GLACE (ICE SHELVES) ET ICEBERGS
!       Agnostique à la grille.
!
!       ENTRÉES  : tms, dxs1, dxs2, area (modules)
!       SORTIES DANS LES MODULES :
!         tmics, areicebn, areicebs
!         iicebern1/2, jicebern1/2, iicebers1/2, jicebers1/2
!==============================================================================
      call define_ice_shelves(
     &  imax, jmax, ks2,
     &  tms, dxs1, dxs2, area,
     &  tmics,
     &  iicebern1, iicebern2, jicebern1, jicebern2,
     &  iicebers1, iicebers2, jicebers1, jicebers2,
     &  areicebn, areicebs,
     &  mouchard_id)

      WRITE(*,*) "DEBUG :: w == ", sum(w)
!==============================================================================
! ÉTAPE 8 : ÉCRITURE DE CONTRÔLE DANS LE FICHIER "MOUCHARD" (si nflag==2)
!           Placé APRÈS G5c pour que toutes les variables soient calculées.
!==============================================================================
      if (nflag.eq.2) then
        call write_mouchard_defgrid(
     &    mouchard_id, clio3_out_id,
     &    imax, jmax, kmax, ncomax, nlpmax,
     &    icheck, jcheck, kcheck,
     &    dtb, dtu, dts, dx, dy, pi,
     &    cmx, cmy, z, dz, zw, dzw,
     &    ahu, ahe, ahs, ai, slopemax, aitd, slopmgm,
     &    afilt, ahh, avv, avnub, avnu0, avkb, avk0,
     &    rifumx, rifsmx, alphah, alphgr, algrmn, alphmi, alphaz,
     &    jsep, jnorth, dxwi, dywj, xwi1, ywj1,
     &    dxaj, dyai, xaj1, yai1, xwpoln,
     &    jsepar, jeq, ijsdl, ijudl, jcl1, jcl2, jdl1, jdl2,
     &    nvhsf, ndhsf, tithsf, ishsf, iehsf, jshsf, jehsf,
     &    bering, iberp, jberp, ibera, jbera,
     &    n1coin, n2coin, n3coin, n4coin,
     &    nxyslp, nxslp, nxslp0, nyslp0, kfond, ju1, ju2, imu1, imu2,
     &    is1, is2, iu1, iu2, iuf1, iuf2, isf1, isf2,
     &    iszon, iezon, nbsmax, kfs, volz,
     &    ipt0v, jpt0v, npt0v, kmax,
     &    xslon, yslat, xulon, yulat, zlatt, zlont,
     &    xsedg, ysedg, xuedg, yuedg,
     &    xslonp, yslatp, xulonp, yulatp,
     &    angle, xang1, xang2, dx1_loc, dx2_loc,
     &    h1_loc, h2_loc, d1d2_loc, d2d1_loc,
     &    smx, smy, cmxy, smxy, cmxdy, cmydx, fs2cor, covrai,
     &    wght, akappa, alambd, bkappa,
     &    dxs1, dxs2, dxc1, dxc2, area, zfn,
     &    unsdz, unsdzw, zbath_loc, dzbath_loc, unszw_loc,
     &    kbath, kbath1, kbath2,
     &    tms, tmu, hs, hu, hux, huy, unshu,
     &    kniv, kfu,
     &    ctmi, aire, unsvol, zurfow,
     &    i1coin, i2coin, i3coin, i4coin,
     &    kslp, lslp, ijslp,
     &    msks, msku, zindfa, dfhu, dfhv, tauc, reslum,
     &    tmics, areicebn, areicebs,
     &    iicebern1, iicebern2, jicebern1, jicebern2,
     &    iicebers1, iicebers2, jicebers1, jicebers2)
      endif

      WRITE(*,*) "DEBUG :: reslum == ", sum(reslum), __LINE__

      return

      end subroutine defgrid

#include "compute_grid_coords.finc"
#include "compute_metric_coeffs.finc"
#include "compute_ice_metrics.finc"
#include "read_bathymetry_and_vscale.finc"
#include "build_ocean_mask_and_indices-v2.finc"
#include "compute_integrals-v2.finc"

!==============================================================================
! SUBROUTINE write_mouchard_defgrid
!
! Écriture de contrôle dans le fichier "mouchard" (débogage, nflag==2).
! Extrait de la section 8 originale de defgrid.
! Appelée directement depuis defgrid avec accès à toutes les variables de modules.
!==============================================================================
      SUBROUTINE write_mouchard_defgrid(
     &  mouchard_id, clio3_out_id,
     &  imax, jmax, kmax, ncomax, nlpmax,
     &  icheck, jcheck, kcheck,
     &  dtb, dtu, dts, dx, dy, pi,
     &  cmx, cmy, z, dz, zw, dzw,
     &  ahu, ahe, ahs, ai, slopemax, aitd, slopmgm,
     &  afilt, ahh, avv, avnub, avnu0, avkb, avk0,
     &  rifumx, rifsmx, alphah, alphgr, algrmn, alphmi, alphaz,
     &  jsep, jnorth, dxwi, dywj, xwi1, ywj1,
     &  dxaj, dyai, xaj1, yai1, xwpoln,
     &  jsepar, jeq, ijsdl, ijudl, jcl1, jcl2, jdl1, jdl2,
     &  nvhsf, ndhsf, tithsf, ishsf, iehsf, jshsf, jehsf,
     &  bering, iberp, jberp, ibera, jbera,
     &  n1coin, n2coin, n3coin, n4coin,
     &  nxyslp, nxslp, nxslp0, nyslp0, kfond, ju1, ju2, imu1, imu2,
     &  is1, is2, iu1, iu2, iuf1, iuf2, isf1, isf2,
     &  iszon, iezon, nbsmax, kfs, volz,
     &  ipt0v, jpt0v, npt0v, kmx,
     &  xslon, yslat, xulon, yulat, zlatt, zlont,
     &  xsedg, ysedg, xuedg, yuedg,
     &  xslonp, yslatp, xulonp, yulatp,
     &  angle, xang1, xang2, dx1, dx2, h1, h2, d1d2, d2d1,
     &  smx, smy, cmxy, smxy, cmxdy, cmydx, fs2cor, covrai,
     &  wght, akappa, alambd, bkappa,
     &  dxs1, dxs2, dxc1, dxc2, area, zfn,
     &  unsdz, unsdzw, zbath, dzbath, unszw,
     &  kbath, kbath1, kbath2,
     &  tms, tmu, hs, hu, hux, huy, unshu,
     &  kniv, kfu,
     &  ctmi, aire, unsvol, zurfow,
     &  i1coin, i2coin, i3coin, i4coin,
     &  kslp, lslp, ijslp,
     &  msks, msku, zindfa, dfhu, dfhv, tauc, reslum,
     &  tmics, areicebn, areicebs,
     &  iicebern1, iicebern2, jicebern1, jicebern2,
     &  iicebers1, iicebers2, jicebers1, jicebers2)

      use global_constants_mod, only: dblp=>dp, ip
      use para0_mod,  only: nsmax
      use para_mod,   only: nbsmx => nbsmax
      implicit none

!------------------------------------------------------------------------------
!-- Arguments existants
!------------------------------------------------------------------------------
      integer(ip), intent(in) :: mouchard_id, clio3_out_id
      integer(ip), intent(in) :: imax, jmax, kmax, kmx, ncomax, nlpmax
      integer(ip), intent(in) :: icheck, jcheck, kcheck
      integer(ip), intent(in) :: jsepar, jeq, ijsdl, ijudl
      integer(ip), intent(in) :: jcl1, jcl2, jdl1, jdl2
      integer(ip), intent(in) :: nvhsf, ndhsf, kfond
      integer(ip), intent(in) :: ju1, ju2, imu1, imu2
      integer(ip), intent(in) :: iberp, jberp, ibera, jbera
      integer(ip), intent(in) :: nxyslp, nxslp, nxslp0, nyslp0
      integer(ip), intent(in) :: nbsmax, npt0v
      integer(ip), intent(in) :: jsep(imax), jnorth(imax)
      integer(ip), intent(in) :: is1(jmax), is2(jmax)
      integer(ip), intent(in) :: iu1(jmax), iu2(jmax)
      integer(ip), intent(in) :: iuf1(jmax), iuf2(jmax)
      integer(ip), intent(in) :: isf1(jmax), isf2(jmax)
      integer(ip), intent(in) :: iszon(jmax,0:nbsmax), iezon(jmax,0:nbsmax)
      integer(ip), intent(in) :: kfs(imax,jmax)
      integer(ip), intent(in) :: n1coin(kmax), n2coin(kmax)
      integer(ip), intent(in) :: n3coin(kmax), n4coin(kmax)
      integer(ip), intent(in) :: ishsf(*), iehsf(*), jshsf(*), jehsf(*)
      integer(ip), intent(in) :: ipt0v(10), jpt0v(10)
      real(dblp),  intent(in) :: dtb, dtu, dts(kmax), dx, dy, pi
      real(dblp),  intent(in) :: cmx(imax,jmax,0:3), cmy(imax,jmax,0:3)
      real(dblp),  intent(in) :: z(kmax+1), dz(kmax), zw(kmax+1), dzw(kmax+1)
      real(dblp),  intent(in) :: ahu, ahe, ahs(kmax),ai(kmax),slopemax(kmax)
      real(dblp),  intent(in) :: aitd(kmax), slopmgm(kmax)
      real(dblp),  intent(in) :: afilt, ahh, avv, avnub(kmax)
      real(dblp),  intent(in) :: avnu0(kmax), avkb(kmax), avk0(kmax)
      real(dblp),  intent(in) :: rifumx, rifsmx
      real(dblp),  intent(in) :: alphah(2),alphgr(nsmax), algrmn(nsmax)
      real(dblp),  intent(in) :: alphmi(kmax), alphaz(kmax), bering
      real(dblp),  intent(in) :: dxwi, dywj, xwi1, ywj1
      real(dblp),  intent(in) :: dxaj, dyai, xaj1, yai1, xwpoln
      real(dblp),  intent(in) :: volz
      character(len=8), intent(in) :: tithsf(*)

!------------------------------------------------------------------------------
!-- Nouveaux arguments : sorties des 9 sous-routines
!------------------------------------------------------------------------------
!-- [G1] compute_grid_coords
      real(dblp), intent(in) :: xslon(imax,jmax), yslat(imax,jmax)
      real(dblp), intent(in) :: xulon(imax,jmax), yulat(imax,jmax)
      real(dblp), intent(in) :: zlatt(imax,jmax), zlont(imax,jmax)
      real(dblp), intent(in) :: xsedg(imax,jmax,4), ysedg(imax,jmax,4)
      real(dblp), intent(in) :: xuedg(imax,jmax,4), yuedg(imax,jmax,4)
      real(dblp), intent(in) :: xslonp(imax+1,jmax+1), yslatp(imax+1,jmax+1)
      real(dblp), intent(in) :: xulonp(imax+1,jmax+1), yulatp(imax+1,jmax+1)
      real(dblp), intent(in) :: angle(imax,jmax)
      real(dblp), intent(in) :: xang1(imax,jmax), xang2(imax,jmax)
      real(dblp), intent(in) :: dx1(imax,jmax), dx2(imax,jmax)
      real(dblp), intent(in) :: h1(imax,jmax), h2(imax,jmax)
      real(dblp), intent(in) :: d1d2(imax,jmax), d2d1(imax,jmax)
!-- [G2] compute_metric_coeffs
      real(dblp), intent(in) :: smx(imax,jmax,0:3), smy(imax,jmax,0:3)
      real(dblp), intent(in) :: cmxy(imax,jmax,0:3), smxy(imax,jmax,0:3)
      real(dblp), intent(in) :: cmxdy(imax,jmax), cmydx(imax,jmax)
      real(dblp), intent(in) :: fs2cor(imax,jmax), covrai(imax,jmax)
!-- [G3a] compute_ice_metrics
      real(dblp), intent(in) :: wght(imax,jmax,2,2)
      real(dblp), intent(in) :: akappa(imax,jmax,2,2)
      real(dblp), intent(in) :: alambd(imax,jmax,2,2,2,2)
      real(dblp), intent(in) :: bkappa(imax,jmax,2,2)
      real(dblp), intent(in) :: dxs1(imax,jmax), dxs2(imax,jmax)
      real(dblp), intent(in) :: dxc1(imax,jmax), dxc2(imax,jmax)
      real(dblp), intent(in) :: area(imax,jmax), zfn(imax,jmax)
!-- [G3b] read_bathymetry_and_vscale
      real(dblp),  intent(in) :: unsdz(kmax), unsdzw(kmax+1)
      real(dblp),  intent(in) :: zbath(kmax), dzbath(kmax), unszw(kmax)
      integer(ip), intent(in) :: kbath(imax,jmax)
      integer(ip), intent(in) :: kbath1(imax), kbath2(imax)
!-- [G4a] build_ocean_mask
      real(dblp),  intent(in) :: tms(imax,jmax,kmax), tmu(imax,jmax,kmax)
      real(dblp),  intent(in) :: hs(imax,jmax), hu(imax,jmax)
      real(dblp),  intent(in) :: hux(imax,jmax), huy(imax,jmax)
      real(dblp),  intent(in) :: unshu(imax,jmax)
      integer(ip), intent(in) :: kniv(imax,jmax,-1:1)
      integer(ip), intent(in) :: kfu(imax,jmax)
!-- [G5a] compute_integrals
      real(dblp),  intent(in) :: ctmi(imax,jmax,kmax,0:1)
      real(dblp),  intent(in) :: aire(imax,jmax)
      real(dblp),  intent(in) :: unsvol, zurfow
      integer(ip), intent(in) :: i1coin(ncomax,kmax), i2coin(ncomax,kmax)
      integer(ip), intent(in) :: i3coin(ncomax,kmax), i4coin(ncomax,kmax)
      integer(ip), intent(in) :: kslp(nlpmax), lslp(nlpmax), ijslp(nlpmax)
!-- [G5b] compute_surface_physics
      integer(ip), intent(in) :: msks(imax,jmax), msku(imax,jmax)
      real(dblp),  intent(in) :: zindfa(imax,jmax)
      real(dblp),  intent(in) :: dfhu(imax,jmax), dfhv(imax,jmax)
      real(dblp),  intent(in) :: tauc(imax,jmax)
      real(dblp),  intent(in) :: reslum(imax,jmax,kmax)
!-- [G5c] define_ice_shelves
      real(dblp),  intent(in) :: tmics(imax,jmax)
      real(dblp),  intent(in) :: areicebn, areicebs
      integer(ip), intent(in) :: iicebern1, iicebern2, jicebern1, jicebern2
      integer(ip), intent(in) :: iicebers1, iicebers2, jicebers1, jicebers2

!------------------------------------------------------------------------------
!-- Variables locales pour l'écriture
!------------------------------------------------------------------------------
      integer(ip) :: i, j, k, n, nv, ii1, ii2, jj1, jj2
      integer(ip) :: nfrc
      character(len=1)  :: cc1(0:9), fmt1
      character(len=8)  :: fmt_kfs
      character(len=30) :: fmt, fmtrep

      ! Caractères pour la carte ASCII du fond kfs
      do n=1,9 ; cc1(n)='-' ; enddo
      cc1(0)='|' ; cc1(5)='5'

      WRITE(*,*) "DEBUG :: reslum == ", sum(reslum), __LINE__

!==============================================================================
!-- SECTION EXISTANTE (inchangée)
!==============================================================================

      write(mouchard_id,*) 'dtb, dtu, dts(k) :'
      write(mouchard_id,*)  dtb, dtu, dts
      write(mouchard_id,*) 'dx, dy, pi :'
      write(mouchard_id,*)  dx, dy, pi
      write(mouchard_id,*) 'icheck, jcheck, kcheck :'
      write(mouchard_id,*)  icheck, jcheck, kcheck

      if (kcheck.eq.0 .and. icheck.ge.1 .and. icheck.le.imax) then
        do jj1=1,jmax,10
          jj2=min(jj1+9,jmax)
          write(mouchard_id,'(3(A,I3))') 'cmx(0), i=',icheck,
     &      ' de j=',jj1,' a j=',jj2
          write(mouchard_id,*) (cmx(icheck,j,0),j=jj1,jj2)
        enddo
        do jj1=1,jmax,10
          jj2=min(jj1+9,jmax)
          write(mouchard_id,'(3(A,I3))') 'cmy(0), i=',icheck,
     &      ' de j=',jj1,' a j=',jj2
          write(mouchard_id,*) (cmy(icheck,j,0),j=jj1,jj2)
        enddo
        do jj1=1,jmax,10
          jj2=min(jj1+9,jmax)
          write(mouchard_id,'(3(A,I3))') 'cmx(3), i=',icheck,
     &      ' de j=',jj1,' a j=',jj2
          write(mouchard_id,*) (cmx(icheck,j,3),j=jj1,jj2)
        enddo
        do jj1=1,jmax,10
          jj2=min(jj1+9,jmax)
          write(mouchard_id,'(3(A,I3))') 'cmy(3), i=',icheck,
     &      ' de j=',jj1,' a j=',jj2
          write(mouchard_id,*) (cmy(icheck,j,3),j=jj1,jj2)
        enddo
      endif

      if (kcheck.eq.0 .and. jcheck.ge.1 .and. jcheck.le.jmax) then
        do ii1=1,imax,10
          ii2=min(ii1+9,imax)
          write(mouchard_id,'(3(A,I3))') 'cmx(0), j=',jcheck,
     &      ' de i=',ii1,' a i=',ii2
          write(mouchard_id,*) (cmx(i,jcheck,0),i=ii1,ii2)
        enddo
        do ii1=1,imax,10
          ii2=min(ii1+9,imax)
          write(mouchard_id,'(3(A,I3))') 'cmy(0), j=',jcheck,
     &      ' de i=',ii1,' a i=',ii2
          write(mouchard_id,*) (cmy(i,jcheck,0),i=ii1,ii2)
        enddo
        do ii1=1,imax,10
          ii2=min(ii1+9,imax)
          write(mouchard_id,'(3(A,I3))') 'cmx(3), j=',jcheck,
     &      ' de i=',ii1,' a i=',ii2
          write(mouchard_id,*) (cmx(i,jcheck,3),i=ii1,ii2)
        enddo
        do ii1=1,imax,10
          ii2=min(ii1+9,imax)
          write(mouchard_id,'(3(A,I3))') 'cmy(3), j=',jcheck,
     &      ' de i=',ii1,' a i=',ii2
          write(mouchard_id,*) (cmy(i,jcheck,3),i=ii1,ii2)
        enddo
      endif

      write(mouchard_id,*) 'z :' ; write(mouchard_id,*) z
      write(mouchard_id,*) 'dz :'; write(mouchard_id,*) dz
      write(mouchard_id,*) 'zw :'; write(mouchard_id,*) zw
      write(mouchard_id,*) 'dzw:'; write(mouchard_id,*) dzw

      write(mouchard_id,*) 'ahu, ahe :' ; write(mouchard_id,*) ahu, ahe
      write(mouchard_id,*) 'ahs :'      ; write(mouchard_id,*) ahs
      write(mouchard_id,*) 'ai  :'      ; write(mouchard_id,*) ai
      write(mouchard_id,*) 'slopemax :' ; write(mouchard_id,*) slopemax
      write(mouchard_id,*) 'aitd :'     ; write(mouchard_id,*) aitd
      write(mouchard_id,*) 'slopmgm :'  ; write(mouchard_id,*) slopmgm
      write(mouchard_id,*) 'afilt,ahh,avv :', afilt, ahh, avv
      write(mouchard_id,*) 'avnub :'    ; write(mouchard_id,*) avnub
      write(mouchard_id,*) 'avnu0 :'    ; write(mouchard_id,*) avnu0
      write(mouchard_id,*) 'avkb :'     ; write(mouchard_id,*) avkb
      write(mouchard_id,*) 'avk0 :'     ; write(mouchard_id,*) avk0
      write(mouchard_id,*) 'rifumx, rifsmx :', rifumx, rifsmx
      write(mouchard_id,*) 'alphah :'   ; write(mouchard_id,*) alphah
      WRITE(*,*) "DEBUG alphah :: ", alphah
      write(mouchard_id,*) 'alphgr :'   ; write(mouchard_id,*) alphgr
      write(mouchard_id,*) 'algrmn :'   ; write(mouchard_id,*) algrmn
      write(mouchard_id,*) 'alphmi :'   ; write(mouchard_id,*) alphmi
      write(mouchard_id,*) 'alphaz :'   ; write(mouchard_id,*) alphaz

      write(mouchard_id,'(A)') 'Indices linked with the domain :'
      write(mouchard_id,*) 'jsep(i) :'
      write(mouchard_id,'(20I4)') (jsep(i),i=1,imax)
      write(mouchard_id,*) 'jnorth(i) :'
      write(mouchard_id,'(20I4)') (jnorth(i),i=1,imax)

      write(mouchard_id,*) 'dxwi, dywj, xwi1, ywj1 :'
      write(mouchard_id,*)  dxwi, dywj, xwi1, ywj1
      write(mouchard_id,*) 'dxaj, dyai, xaj1, yai1, xwpoln :'
      write(mouchard_id,*)  dxaj, dyai, xaj1, yai1, xwpoln
      write(mouchard_id,*) 'jsepar, jeq, ijsdl, ijudl :'
      write(mouchard_id,*)  jsepar, jeq, ijsdl, ijudl
      write(mouchard_id,*) 'jcl1, jcl2, jdl1, jdl2 :'
      write(mouchard_id,*)  jcl1, jcl2, jdl1, jdl2
      write(mouchard_id,'(A,2I4)')
     & 'Detroits : nvhsf,ndhsf =', nvhsf, ndhsf
      do nv=1,max(nvhsf,ndhsf)
        write(mouchard_id,'(A,4I4)') tithsf(nv),
     &    ishsf(nv), iehsf(nv), jshsf(nv), jehsf(nv)
      enddo
      write(mouchard_id,*)
      write(mouchard_id,*) 'bering :', bering
      write(mouchard_id,*) 'iberp, jberp, ibera, jbera :'
      write(mouchard_id,*)  iberp, jberp, ibera, jbera
      write(mouchard_id,*) 'cmx,cmy(iberp-1,jberp,2) :'
      write(mouchard_id,*)  cmx(iberp-1,jberp,2), cmy(iberp-1,jberp,2)
      write(mouchard_id,*) 'cmx,cmy(iberp,jberp,2) :'
      write(mouchard_id,*)  cmx(iberp,jberp,2), cmy(iberp,jberp,2)
      write(mouchard_id,*) 'cmx,cmy(ibera-1,jbera,2) :'
      write(mouchard_id,*)  cmx(ibera-1,jbera,2), cmy(ibera-1,jbera,2)
      write(mouchard_id,*) 'cmx,cmy(ibera,jbera,2) :'
      write(mouchard_id,*)  cmx(ibera,jbera,2), cmy(ibera,jbera,2)

      do n=1,npt0v
        write(mouchard_id,*) 'Ajout ile surf=0 en (i,j)=',
     &    ipt0v(n), jpt0v(n)
      enddo

      write(mouchard_id,*) 'n1coin(k) :'
      write(mouchard_id,'(20I4)') (n1coin(k),k=1,kmax)
      write(mouchard_id,*) 'n2coin(k) :'
      write(mouchard_id,'(20I4)') (n2coin(k),k=1,kmax)
      write(mouchard_id,*) 'n3coin(k) :'
      write(mouchard_id,'(20I4)') (n3coin(k),k=1,kmax)
      write(mouchard_id,*) 'n4coin(k) :'
      write(mouchard_id,'(20I4)') (n4coin(k),k=1,kmax)

      write(mouchard_id,*) 'nxyslp, nXslope, nYslope, kfond :'
     &  , nxyslp, nxslp0, nyslp0, kfond
      write(mouchard_id,*) 'ju1, ju2, imu1, imu2 :'
      write(mouchard_id,*)  ju1, ju2, imu1, imu2

      write(mouchard_id,'(A)') 'Indices Debut(=1)/Fin(=2) :'
      write(mouchard_id,'(A)') '  j | is1,is2   iu1,iu2  '//
     &  'iuf1,iuf2 isf1,isf2     (iszon,iezon)(0/1/2/3) :'
      do j=jmax,1,-1
        write(mouchard_id,'(I3,1X,A1,4(2I4,2X),3(2I4,1X),2I4)') j,'|',
     &    is1(j),is2(j), iu1(j),iu2(j),
     &    iuf1(j),iuf2(j), isf1(j),isf2(j),
     &    (iszon(j,n),iezon(j,n),n=0,nbsmax)
      enddo

!-- Carte ASCII du fond kfs
      write(mouchard_id,*) 'array kfs(i,j) :'
      if (kmx.le.15) then ; nfrc=125 ; fmt1='Z' ; else
                             nfrc=41  ; fmt1='I' ; endif
      write(fmtrep,'(A,I3,A)') '(',iabs(nfrc),'A1)'
      if (fmt1.eq.'Z') then
        fmt_kfs='Z1'
      else
        fmt_kfs='I3'
      endif
      do ii1=1,imax,nfrc
        ii2=min(ii1+nfrc-1,imax)
        write(fmt,'(A,I3,2A)') '(',(ii2-ii1+1),trim(fmt_kfs),',A1,I3)'
        write(mouchard_id,*) 'portion de i=',ii1,' a i=',ii2
        write(mouchard_id,fmtrep) (cc1(mod(i,10)),i=ii1,ii2)
        do j=jmax,1,-1
          write(mouchard_id,fmt) (kfs(i,j),i=ii1,ii2),'|',j
        enddo
      enddo
      write(mouchard_id,*) 'volz, volume-real :'
      write(mouchard_id,*)  volz, volz*dx*dy

!==============================================================================
!-- NOUVELLES SECTIONS : sorties des 9 sous-routines
!==============================================================================

!------------------------------------------------------------------------------
!-- [G1] compute_grid_coords : coordonnées géographiques
!        Variables 1D : aucune
!        Variables nD : toutes -> sum()
!------------------------------------------------------------------------------
      write(mouchard_id,'(/,A)') '=== [G1] compute_grid_coords ==='
      write(mouchard_id,*) 'sum(xslon)  :', sum(xslon)
      write(mouchard_id,*) 'sum(yslat)  :', sum(yslat)
      write(mouchard_id,*) 'sum(xulon)  :', sum(xulon)
      write(mouchard_id,*) 'sum(yulat)  :', sum(yulat)
      write(mouchard_id,*) 'sum(zlatt)  :', sum(zlatt)
      write(mouchard_id,*) 'sum(zlont)  :', sum(zlont)
      write(mouchard_id,*) 'sum(xsedg)  :', sum(xsedg)
      write(mouchard_id,*) 'sum(ysedg)  :', sum(ysedg)
      write(mouchard_id,*) 'sum(xuedg)  :', sum(xuedg)
      write(mouchard_id,*) 'sum(yuedg)  :', sum(yuedg)
      write(mouchard_id,*) 'sum(xslonp) :', sum(xslonp)
      write(mouchard_id,*) 'sum(yslatp) :', sum(yslatp)
      write(mouchard_id,*) 'sum(xulonp) :', sum(xulonp)
      write(mouchard_id,*) 'sum(yulatp) :', sum(yulatp)
      write(mouchard_id,*) 'sum(angle)  :', sum(angle)
      write(mouchard_id,*) 'sum(xang1)  :', sum(xang1)
      write(mouchard_id,*) 'sum(xang2)  :', sum(xang2)
      write(mouchard_id,*) 'sum(dx1)    :', sum(dx1)
      write(mouchard_id,*) 'sum(dx2)    :', sum(dx2)
      write(mouchard_id,*) 'sum(h1)     :', sum(h1)
      write(mouchard_id,*) 'sum(h2)     :', sum(h2)
      write(mouchard_id,*) 'sum(d1d2)   :', sum(d1d2)
      write(mouchard_id,*) 'sum(d2d1)   :', sum(d2d1)

!------------------------------------------------------------------------------
!-- [G2] compute_metric_coeffs : coefficients métriques
!        Variables 1D : aucune
!        Variables nD : toutes -> sum()
!------------------------------------------------------------------------------
      write(mouchard_id,'(/,A)') '=== [G2] compute_metric_coeffs ==='
      write(mouchard_id,*) 'sum(smx)    :', sum(smx)
      write(mouchard_id,*) 'sum(smy)    :', sum(smy)
      write(mouchard_id,*) 'sum(cmxy)   :', sum(cmxy)
      write(mouchard_id,*) 'sum(smxy)   :', sum(smxy)
      write(mouchard_id,*) 'sum(cmxdy)  :', sum(cmxdy)
      write(mouchard_id,*) 'sum(cmydx)  :', sum(cmydx)
      write(mouchard_id,*) 'sum(fs2cor) :', sum(fs2cor)
      write(mouchard_id,*) 'sum(covrai) :', sum(covrai)

!------------------------------------------------------------------------------
!-- [G3a] compute_ice_metrics : métriques pour la dynamique de la glace
!        Variables 1D : aucune
!        Variables nD : toutes -> sum()
!------------------------------------------------------------------------------
      write(mouchard_id,'(/,A)') '=== [G3a] compute_ice_metrics ==='
      write(mouchard_id,*) 'sum(wght)   :', sum(wght)
      write(mouchard_id,*) 'sum(akappa) :', sum(akappa)
      write(mouchard_id,*) 'sum(alambd) :', sum(alambd)
      write(mouchard_id,*) 'sum(bkappa) :', sum(bkappa)
      write(mouchard_id,*) 'sum(dxs1)   :', sum(dxs1)
      write(mouchard_id,*) 'sum(dxs2)   :', sum(dxs2)
      write(mouchard_id,*) 'sum(dxc1)   :', sum(dxc1)
      write(mouchard_id,*) 'sum(dxc2)   :', sum(dxc2)
      write(mouchard_id,*) 'sum(area)   :', sum(area)
      write(mouchard_id,*) 'sum(zfn)    :', sum(zfn)

!------------------------------------------------------------------------------
!-- [G3b] read_bathymetry_and_vscale : échelle verticale et bathymétrie
!        Variables 1D : zbath, dzbath, unszw, unsdz, unsdzw, kbath1, kbath2
!        Variables nD : kbath -> sum()
!------------------------------------------------------------------------------
      write(mouchard_id,'(/,A)') '=== [G3b] read_bathymetry_and_vscale ==='
      write(mouchard_id,*) 'zbath  :'   ; write(mouchard_id,*) zbath
      write(mouchard_id,*) 'dzbath :'   ; write(mouchard_id,*) dzbath
      write(mouchard_id,*) 'unszw  :'   ; write(mouchard_id,*) unszw
      write(mouchard_id,*) 'unsdz  :'   ; write(mouchard_id,*) unsdz
      write(mouchard_id,*) 'unsdzw :'   ; write(mouchard_id,*) unsdzw
      write(mouchard_id,*) 'kbath1(i) :'
      write(mouchard_id,'(20I4)') (kbath1(i),i=1,imax)
      write(mouchard_id,*) 'kbath2(i) :'
      write(mouchard_id,'(20I4)') (kbath2(i),i=1,imax)
      write(mouchard_id,*) 'sum(kbath) :', sum(kbath)

!------------------------------------------------------------------------------
!-- [G4a] build_ocean_mask : masques terre/mer et profondeurs
!        Variables 1D : is1, is2, iu1, iu2 (deja dans section existante)
!        Variables nD : tms, tmu, hs, hu, hux, huy, unshu, kniv, kfu -> sum()
!------------------------------------------------------------------------------
      write(mouchard_id,'(/,A)') '=== [G4a] build_ocean_mask ==='
      write(mouchard_id,*) 'sum(tms)    :', sum(tms)
      write(mouchard_id,*) 'sum(tmu)    :', sum(tmu)
      write(mouchard_id,*) 'sum(hs)     :', sum(hs)
      write(mouchard_id,*) 'sum(hu)     :', sum(hu)
      write(mouchard_id,*) 'sum(hux)    :', sum(hux)
      write(mouchard_id,*) 'sum(huy)    :', sum(huy)
      write(mouchard_id,*) 'sum(unshu)  :', sum(unshu)
      write(mouchard_id,*) 'sum(kniv)   :', sum(kniv)
      write(mouchard_id,*) 'sum(kfu)    :', sum(kfu)

!------------------------------------------------------------------------------
!-- [G4b] build_domain_indices : indices de domaine et Béring
!        Variables scalaires : ju1, ju2, imu1, imu2 (deja dans section existante)
!        Variables 1D : isf1, isf2, iuf1, iuf2 (deja dans section existante)
!        -> bering deja ecrit dans section existante
!        Rien de nouveau a ecrire ici.
!------------------------------------------------------------------------------
      write(mouchard_id,'(/,A)') '=== [G4b] build_domain_indices ==='
      write(mouchard_id,*) '(ju1, ju2, imu1, imu2, bering, isf1/2,'//
     &  ' iuf1/2 deja ecrits dans la section existante)'

!------------------------------------------------------------------------------
!-- [G5a] compute_integrals : surfaces, volumes, coins, pentes
!        Variables scalaires : volz (deja), unsvol, zurfow, nxslp, nxyslp (deja)
!        Variables 1D : n1-4coin (deja), kslp, lslp, ijslp
!        Variables nD : ctmi, aire, i1-4coin -> sum()
!------------------------------------------------------------------------------
      write(mouchard_id,'(/,A)') '=== [G5a] compute_integrals ==='
      write(mouchard_id,*) 'unsvol      :', unsvol
      write(mouchard_id,*) 'zurfow      :', zurfow
      write(mouchard_id,*) 'kslp  (nl=1..nxyslp) :'
      write(mouchard_id,'(20I4)') (kslp(n),n=1,nxyslp)
      write(mouchard_id,*) 'lslp  (nl=1..nxyslp) :'
      write(mouchard_id,'(20I4)') (lslp(n),n=1,nxyslp)
      write(mouchard_id,*) 'ijslp (nl=1..nxyslp) :'
      write(mouchard_id,'(20I4)') (ijslp(n),n=1,nxyslp)
      write(mouchard_id,*) 'sum(ctmi)   :', sum(ctmi)
      write(mouchard_id,*) 'sum(aire)   :', sum(aire)
      write(mouchard_id,*) 'sum(i1coin) :', sum(i1coin)
      write(mouchard_id,*) 'sum(i2coin) :', sum(i2coin)
      write(mouchard_id,*) 'sum(i3coin) :', sum(i3coin)
      write(mouchard_id,*) 'sum(i4coin) :', sum(i4coin)

!------------------------------------------------------------------------------
!-- [G5b] compute_surface_physics : physique de surface et absorption solaire
!        Variables scalaires : aucune
!        Variables nD : msks, msku, zindfa, dfhu, dfhv, tauc, reslum -> sum()
!------------------------------------------------------------------------------
      write(mouchard_id,'(/,A)') '=== [G5b] compute_surface_physics ==='
      write(mouchard_id,*) 'sum(msks)   :', sum(msks)
      write(mouchard_id,*) 'sum(msku)   :', sum(msku)
      write(mouchard_id,*) 'sum(zindfa) :', sum(zindfa)
      write(mouchard_id,*) 'sum(dfhu)   :', sum(dfhu)
      write(mouchard_id,*) 'sum(dfhv)   :', sum(dfhv)
      write(mouchard_id,*) 'sum(tauc)   :', sum(tauc)
      write(mouchard_id,*) 'sum(reslum) :', sum(reslum)

!------------------------------------------------------------------------------
!-- [G5c] define_ice_shelves : plateformes de glace et icebergs
!        Variables scalaires : areicebn, areicebs, iicebern1/2, jicebern1/2,
!                              iicebers1/2, jicebers1/2
!        Variables nD : tmics -> sum()
!------------------------------------------------------------------------------
      write(mouchard_id,'(/,A)') '=== [G5c] define_ice_shelves ==='
      write(mouchard_id,*) 'areicebn    :', areicebn
      write(mouchard_id,*) 'areicebs    :', areicebs
      write(mouchard_id,*) 'iicebern1, iicebern2 :', iicebern1, iicebern2
      write(mouchard_id,*) 'jicebern1, jicebern2 :', jicebern1, jicebern2
      write(mouchard_id,*) 'iicebers1, iicebers2 :', iicebers1, iicebers2
      write(mouchard_id,*) 'jicebers1, jicebers2 :', jicebers1, jicebers2
      write(mouchard_id,*) 'sum(tmics)  :', sum(tmics)

      end subroutine write_mouchard_defgrid

#include "compute_surface_physics.finc"
#include "define_ice_shelves.finc"
