!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2026 iLOVECLIM / LUDUS coding group

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [defgrid_kernel_mod]
!
!>     @brief Noyau AGNOSTIQUE A LA GRILLE de l'initialisation de defgrid (etapes G3a a G5c).
!
!      DESCRIPTION:
!>     Regroupe les sept sous-routines de defgrid qui ne dependent PAS du type de grille : elles ne consomment que des
!!     metriques deja calculees (h1/h2, dx1/dx2, fs2cor, cmxy...) et des masques, transmis exclusivement par arguments. Ce
!!     noyau est reutilisable tel quel quelle que soit la grille produite par defgrid_compute_mod.
!!
!!       [G3a] compute_ice_metrics        -- poids d'interpolation, tenseurs rheologiques, longueurs de cotes, surfaces
!!       [G3b] read_bathymetry_and_vscale -- lecture bathymetrie, profils verticaux z/zw/dz/dzw, fermetures de bordures
!!       [G4a] build_ocean_mask           -- masques terre/mer tms/tmu, profondeurs hs/hu, indices de fond kniv/kfs/kfu
!!       [G4b] build_domain_indices       -- bornes de flux, bornes globales, coefficient geostrophique de Bering
!!       [G5a] compute_integrals          -- surfaces ctmi/aire, volume, coins bathymetriques, listes de pentes
!!       [G5b] compute_surface_physics    -- masques de surface, absorption solaire, fichiers lat/long/mask
!!       [G5c] define_ice_shelves         -- plateformes de glace tmics, zones de fonte des icebergs N/S
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module defgrid_kernel_mod

       use global_constants_mod, only: dblp => dp, ip

       implicit none

       private

       public :: compute_ice_metrics
       public :: read_bathymetry_and_vscale
       public :: build_ocean_mask
       public :: build_domain_indices
       public :: compute_integrals
       public :: compute_surface_physics
       public :: define_ice_shelves

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: [compute_ice_metrics]
!
!>     @brief [G3a] Metriques pour la dynamique de la glace de mer (agnostique a la grille).
!
!      DESCRIPTION:
!>     Calcule les poids d'interpolation bilineaire (wght), les coefficients metriques interpoles aux coins et milieux de
!!     cotes, les derivees aux coins, le double-Coriolis (zfn), les tenseurs rheologiques (akappa, alambd, bkappa) et les
!!     longueurs de cotes (dxs1, dxs2, dxc1, dxc2, area). Consomme uniquement h1, h2, d1d2, d2d1, dx1, dx2 et fs2cor.
!!
!!     ENTREES :
!!       imax, jmax                 : dimensions de la grille
!!       h1, h2, d1d2, d2d1, dx1, dx2 : metriques locales (sorties G1)
!!       fs2cor                     : facteur de Coriolis = 2*Omega*sin(lat_U) (sortie G2)
!!       ltest, iberp, ibera, jberp, jbera : type de bassin et detroit de Bering (correction dxs1)
!!
!!     SORTIES :
!!       wght(i,j,2,2)       : poids d'interpolation bilineaire au coin SW (pt U)
!!       akappa(i,j,2,2)     : coefficients pour le gradient de deformation rheologique
!!       alambd(i,j,2,2,2,2) : coefficients pour la divergence du tenseur des contraintes
!!       bkappa(i,j,2,2)     : coefficients pour la divergence (advection/diffusion scalaire)
!!       dxs1, dxs2(i,j)     : longueurs des cotes N et E de la maille S (m)
!!       dxc1, dxc2(i,j)     : longueurs effectives x et y au centre (m)
!!       area(i,j)           : surface de la maille (m2)
!!       zfn(i,j)            : 2 * facteur de Coriolis au pt U
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine compute_ice_metrics(                                   &
        imax, jmax,                                                     &
        h1, h2, d1d2, d2d1, dx1, dx2, fs2cor,                           &
        ltest, iberp, ibera, jberp, jbera,                              &
        wght, akappa, alambd, bkappa,                                   &
        dxs1, dxs2, dxc1, dxc2, area, zfn)

       implicit none

!-- Declaration des arguments
       integer(ip), intent(in) :: imax, jmax, ltest, iberp, ibera, jberp, jbera
       real(dblp),  intent(in) :: h1(imax,jmax), h2(imax,jmax)
       real(dblp),  intent(in) :: d1d2(imax,jmax), d2d1(imax,jmax)
       real(dblp),  intent(in) :: dx1(imax,jmax), dx2(imax,jmax)
       real(dblp),  intent(in) :: fs2cor(imax,jmax)

       real(dblp), intent(out) :: wght(imax,jmax,2,2)
       real(dblp), intent(out) :: akappa(imax,jmax,2,2)
       real(dblp), intent(out) :: alambd(imax,jmax,2,2,2,2)
       real(dblp), intent(out) :: bkappa(imax,jmax,2,2)
       real(dblp), intent(out) :: dxs1(imax,jmax), dxs2(imax,jmax)
       real(dblp), intent(out) :: dxc1(imax,jmax), dxc2(imax,jmax)
       real(dblp), intent(out) :: area(imax,jmax)
       real(dblp), intent(out) :: zfn(imax,jmax)

!-- Variables locales (tableaux internes non exposes en sortie)
       real(dblp) :: h1p(imax,jmax), h2p(imax,jmax)
       real(dblp) :: h1pp(imax,jmax), h2pp(imax,jmax)
       real(dblp) :: d1d2p(imax,jmax), d2d1p(imax,jmax)
       integer(ip) :: i, j, im1, ip1, jm1
       real(dblp)  :: usden

!==============================================================================
!-- 4.1 : Poids d'interpolation bilineaire wght aux coins des mailles
!         wght(i,j,n,m) : poids pour interpoler du coin SW a partir des 4 pts S voisins
!         pondere par les aires relatives des sous-mailles
!==============================================================================
       do j=2,jmax
         jm1=j-1
         do i=1,imax
           im1 = (i-1)+(imax-2)*(1/i)  ! i-1 avec periodicite zonale
           usden = 1.0/((dx1(i,j)+dx1(im1,j))*(dx2(i,j)+dx2(i,jm1)))
           wght(i,j,1,1) = usden*dx1(i,j)  *dx2(i,j)    ! NE du coin = (i,j)
           wght(i,j,1,2) = usden*dx1(i,j)  *dx2(i,jm1)  ! SE du coin = (i,j-1)
           wght(i,j,2,2) = usden*dx1(im1,j)*dx2(i,jm1)  ! SW du coin = (i-1,j-1)
           wght(i,j,2,1) = usden*dx1(im1,j)*dx2(i,j)    ! NW du coin = (i-1,j)
         enddo
       enddo

!==============================================================================
!-- 4.2 : Interpolation de h1, h2 et de leurs derivees aux coins et milieux de cotes
!==============================================================================
       do i=1,imax
         h1p(i,1)   = 0.0
         h2p(i,1)   = 1.0e+20  ! valeur sentinelle (bord j=1, jamais utilise)
         d1d2p(i,1) = 1.0e+20
         d2d1p(i,1) = 0.0
         zfn(i,1)   = fs2cor(i,1)*2
       enddo

       do j=1,jmax-1
         do i=1,imax
           im1 = (i-1)+(imax-2)*(1/i)
           ip1 = (i+1)-(imax-2)*(i/imax)
!-- h1p, h2p : interpolation bilineaire au coin SW (position j+1)
           h1p(i,j+1) = h1(i,j+1)*wght(i,j+1,2,2)                       &
                       +h1(im1,j+1)*wght(i,j+1,1,2)                     &
                       +h1(i,j)  *wght(i,j+1,2,1)                       &
                       +h1(im1,j)*wght(i,j+1,1,1)
           h2p(i,j+1) = h2(i,j+1)*wght(i,j+1,2,2)                       &
                       +h2(im1,j+1)*wght(i,j+1,1,2)                     &
                       +h2(i,j)  *wght(i,j+1,2,1)                       &
                       +h2(im1,j)*wght(i,j+1,1,1)
!-- d1d2p, d2d1p : derivees aux coins par differences finies ponderees
           d1d2p(i,j+1) = 2.0*                                         &
               (dx1(i,j)*(-h1(im1,j)+h1(im1,j+1))                       &
               +dx1(im1,j)*(-h1(i,j)+h1(i,j+1)))/                       &
               ((dx1(i,j)+dx1(im1,j))*(dx2(i,j)+dx2(i,j+1)))
           d2d1p(i,j+1) = 2.0*                                         &
               (dx2(i,j+1)*(h2(i,j)-h2(im1,j))                          &
               +dx2(i,j)*(h2(i,j+1)-h2(im1,j+1)))/                      &
               ((dx1(i,j)+dx1(im1,j))*(dx2(i,j)+dx2(i,j+1)))
!-- h1pp, h2pp : interpoles aux milieux des cotes (pour dxs1, dxs2)
           h1pp(i,j) = (dx2(i,j)*h1(i,j+1)+dx2(i,j+1)*h1(i,j))/         &
                       (dx2(i,j)+dx2(i,j+1))
           h2pp(i,j) = (dx1(i,j)*h2(ip1,j)+dx1(ip1,j)*h2(i,j))/         &
                       (dx1(i,j)+dx1(ip1,j))
           zfn(i,j+1) = 2*fs2cor(i,j+1)
         enddo
       enddo
       do i=1,imax
         ip1 = (i+1)-(imax-2)*(i/imax)
         h1pp(i,jmax) = 0.0  ! bord nord : pas de maille au-dessus
         h2pp(i,jmax) = (dx1(i,jmax)*h2(ip1,jmax)+dx1(ip1,jmax)*h2(i,jmax))/ &
                        (dx1(i,jmax)+dx1(ip1,jmax))
       enddo

!==============================================================================
!-- akappa : coefficients pour le gradient de deformation rheologique
!==============================================================================
       do j=1,jmax
         do i=1,imax
           akappa(i,j,1,1) = 1.0/(2.0*h1(i,j)*dx1(i,j))
           akappa(i,j,1,2) = d1d2(i,j)/(4.0*h1(i,j)*h2(i,j))
           akappa(i,j,2,1) = d2d1(i,j)/(4.0*h1(i,j)*h2(i,j))
           akappa(i,j,2,2) = 1.0/(2.0*h2(i,j)*dx2(i,j))
         enddo
       enddo

!==============================================================================
!-- 4.3 : alambd : divergence du tenseur des contraintes (rheologie Hibler)
!==============================================================================
       do j=2,jmax
         jm1=j-1
         do i=1,imax
           im1=(i-1)+imax*(1/i)
           usden = 1.0/(h1p(i,j)*h2p(i,j)                               &
                  *(dx1(i,j)+dx1(im1,j))                                &
                  *(dx2(i,j)+dx2(i,jm1)))
           alambd(i,j,2,2,2,1) = usden*2.0*dx2(i,j)  *h2(i,jm1)
           alambd(i,j,2,2,2,2) = usden*2.0*dx2(i,jm1)*h2(i,j)
           alambd(i,j,2,2,1,1) = usden*2.0*dx2(i,j)  *h2(im1,jm1)
           alambd(i,j,2,2,1,2) = usden*2.0*dx2(i,jm1)*h2(im1,j)
           alambd(i,j,1,1,2,1) = usden*2.0*dx1(im1,j)*h1(i,jm1)
           alambd(i,j,1,1,1,1) = usden*2.0*dx1(i,j)  *h1(im1,jm1)
           alambd(i,j,1,1,2,2) = usden*2.0*dx1(im1,j)*h1(i,j)
           alambd(i,j,1,1,1,2) = usden*2.0*dx1(i,j)  *h1(im1,j)
           alambd(i,j,1,2,1,1) = usden*d1d2p(i,j)*dx2(i,j)  *dx1(i,j)
           alambd(i,j,1,2,2,1) = usden*d1d2p(i,j)*dx2(i,j)  *dx1(im1,j)
           alambd(i,j,1,2,1,2) = usden*d1d2p(i,j)*dx2(i,jm1)*dx1(i,j)
           alambd(i,j,1,2,2,2) = usden*d1d2p(i,j)*dx2(i,jm1)*dx1(im1,j)
           alambd(i,j,2,1,1,1) = usden*d2d1p(i,j)*dx2(i,j)  *dx1(i,j)
           alambd(i,j,2,1,2,1) = usden*d2d1p(i,j)*dx2(i,j)  *dx1(im1,j)
           alambd(i,j,2,1,1,2) = usden*d2d1p(i,j)*dx2(i,jm1)*dx1(i,j)
           alambd(i,j,2,1,2,2) = usden*d2d1p(i,j)*dx2(i,jm1)*dx1(im1,j)
         enddo
       enddo

!==============================================================================
!-- 4.4 : Longueurs de cotes et surface des mailles
!         dxs1, dxs2 : longueurs aux milieux des cotes N et E (m)
!         dxc1, dxc2 : longueurs effectives au centre (m)
!         area        : surface de la maille (m2)
!==============================================================================
       do j=1,jmax
         do i=1,imax
           dxs1(i,j) = h1pp(i,j)*dx1(i,j)
           dxs2(i,j) = h2pp(i,j)*dx2(i,j)
           dxc1(i,j) = h1(i,j)  *dx1(i,j)
           dxc2(i,j) = h2(i,j)  *dx2(i,j)
           area(i,j)  = dxc1(i,j)*dxc2(i,j)
         enddo
       enddo

!-- bkappa : coefficients pour la divergence (advection/diffusion scalaire)
       do j=2,jmax
         do i=1,imax
           im1=(i-1)+(imax-2)*(1/i)
           bkappa(i,j,1,1) = dxs2(i,j)  /area(i,j)
           bkappa(i,j,1,2) = dxs2(im1,j)/area(i,j)
           bkappa(i,j,2,2) = dxs1(i,j)  /area(i,j)
           bkappa(i,j,2,1) = dxs1(i,j-1)/area(i,j)
         enddo
       enddo
       do i=1,imax
         im1=(i-1)+(imax-2)*(1/i)
         bkappa(i,1,1,1) = dxs2(i,1)  /area(i,1)
         bkappa(i,1,1,2) = dxs2(im1,1)/area(i,1)
         bkappa(i,1,2,2) = dxs1(i,1)  /area(i,1)
         bkappa(i,1,2,1) = dxs1(i,1)  /area(i,1)
       enddo

!-- Correction de dxs1 au detroit de Bering (jonction AA<->WW)
       if (ltest.eq.3 .and. iberp.ne.0) then
         dxs1(ibera-1,jbera-1) = dxs1(iberp,  jberp-1)
         dxs1(ibera,  jbera-1) = dxs1(iberp-1,jberp-1)
       endif

      end subroutine compute_ice_metrics

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: [read_bathymetry_and_vscale]
!
!>     @brief [G3b] Lecture de la bathymetrie et calcul de l'echelle verticale (agnostique a la grille).
!
!      DESCRIPTION:
!>     Lit les fichiers de bathymetrie (bath.om, puis eventuellement un fichier texte pour BATHY>=1), construit les profils
!!     verticaux z, zw, dz, dzw, unsdz, unsdzw, et applique les fermetures aux bordures du bassin ainsi que la verification
!!     de consistance bathymetrique.
!!
!!     ENTREES :
!!       nflag, kfond, ks1, ks2, kmax, imax, jmax : mode, fond fixe, bornes et dimensions
!!       ltest, jeq, jsep(imax)                   : type de bassin et jonction WW/AA
!!       jcl1, jcl2, ims1, ims2                   : bornes cycliques et horizontales
!!       zero                                     : constante 0.0 en double precision
!!
!!     SORTIES :
!!       z, zw(kmax+1)            : profondeurs des centres et interfaces (m)
!!       dz(kmax), dzw(kmax+1)    : epaisseurs et distances entre centres (m)
!!       unsdz(kmax), unsdzw(kmax+1) : inverses
!!       zbath, dzbath(kmax)      : profondeur et epaisseur brutes (bath.om)
!!       unszw(kmax)              : 1/zw(k) (utilise pour hs, hu dans G4a)
!!       kbath(imax,jmax)         : nombre de niveaux actifs (0=terre)
!!       kbath1, kbath2(imax)     : sauvegardes des bordures S et N de kbath
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine read_bathymetry_and_vscale(                            &
        nflag, kfond, ks1, ks2, kmax,                                   &
        imax, jmax, ltest, jeq, jsep,                                   &
        jcl1, jcl2, ims1, ims2,                                         &
        zero,                                                           &
        z, zw, dz, dzw, unsdz, unsdzw,                                  &
        zbath, dzbath, unszw, kbath,                                    &
        kbath1, kbath2)

#if ( BATHY >= 2 || NC_BERG == 2 )
       use comemic_mod,          only: iyear
       use palaeo_timer_mod,     only: palaeo_year
       use update_clio_bathy_tools, only: la_date
#endif
#if ( BRINES >= 3 )
       use comemic_mod,          only: iyear
       use palaeo_timer_mod,     only: palaeo_year
       use brines_mod,           only: la_date_brines
#endif
       implicit none

!-- Declaration des arguments
       integer(ip), intent(in)  :: nflag, kfond, ks1, ks2, kmax
       integer(ip), intent(in)  :: imax, jmax, ltest, jeq, jcl1, jcl2
       integer(ip), intent(in)  :: ims1, ims2
       integer(ip), intent(in)  :: jsep(imax)
       real(dblp),  intent(in)  :: zero

       real(dblp), intent(out)  :: z(kmax+1), zw(kmax+1)
       real(dblp), intent(out)  :: dz(kmax), dzw(kmax+1)
       real(dblp), intent(out)  :: unsdz(kmax), unsdzw(kmax+1)
       real(dblp), intent(out)  :: zbath(kmax), dzbath(kmax), unszw(kmax)
       integer(ip), intent(out) :: kbath(imax,jmax)
       integer(ip), intent(out) :: kbath1(imax), kbath2(imax)

!-- Variables locales
       integer(ip) :: i, j, k
       integer(ip) :: bath_om_id, bath_id
#if ( BATHY >= 2 )
       character(len=30) :: name_file
       character(len=5)  :: charI
#endif
#if ( BATHY >= 1 )
       real        :: test_wet(4), kbath_bas
#endif

!==============================================================================
!-- Lecture du fichier binaire bath.om (profils verticaux + kbath si BATHY==0)
!==============================================================================
       open(newunit=bath_om_id,                                        &
         file='inputdata/clio/bath.om', status='old', form='UNFORMATTED')
       read(bath_om_id) zbath    ! profondeurs des centres de niveaux (m)
       read(bath_om_id) dzbath   ! epaisseurs des niveaux (m)
#if ( BATHY == 0 )
       read(bath_om_id) kbath    ! bathymetrie standard
#endif
       close(bath_om_id)

!==============================================================================
!-- Lecture eventuelle d'un fichier de bathymetrie texte (BATHY >= 1)
!==============================================================================
#if ( NC_BERG == 2 )
       la_date = palaeo_year - iyear
       write(*,*) 'The date for ice sheet update is ', la_date
#endif

#if ( BATHY == 1 )
       open(newunit=bath_id, file='inputdata/clio/bath_new.txt',       &
            status='old', form='FORMATTED')
#elif ( BATHY >= 2 )
       la_date = palaeo_year - iyear
       write(*,*) 'The date for bathy update is ', la_date
       write(charI,'(I5.5)') la_date
       name_file = 'inputdata/clio/bath_'//trim(charI)//'.txt'
       write(*,*) name_file
       open(newunit=bath_id, file=name_file, status='old', form='FORMATTED')
#endif

#if ( BRINES >= 3 )
       la_date_brines = palaeo_year - iyear
#endif

#if ( BATHY >= 1 )
       do i=1,imax
         read(bath_id,*) kbath(i,:)
       enddo
       close(bath_id)
#endif

!==============================================================================
!-- Calcul de l'echelle verticale depuis bath.om
!   Convention : zbath/dzbath sont numerotes du fond (k=1=ks1) vers la surface
!   (k=kmax=ks2), mais z,zw,dz sont stockes fond=1, surface=kmax.
!==============================================================================
       zw(ks2+1) = 0.0   ! interface superieure de la couche de surface = 0 m
       do k=ks2,ks1,-1
         z(k)     = -zbath(ks1+ks2-k)   ! profondeur du centre du niveau k (m, <0)
         dz(k)    = dzbath(ks1+ks2-k)   ! epaisseur du niveau k (m)
         zw(k)    = zw(k+1) - dz(k)     ! profondeur de l'interface inferieure
         unsdz(k) = 1.0 / dz(k)
         if (zw(k).ne.zero) unszw(k) = 1.0 / zw(k)
       enddo

       do k=ks1+1,ks2
         dzw(k)    = z(k) - z(k-1)    ! distance entre centres des niveaux k et k-1
         unsdzw(k) = 1.0 / dzw(k)
       enddo

!==============================================================================
!-- Limitation du fond si kfond > 0 (fond plat / bathymetrie tronquee)
!==============================================================================
       if (kfond.gt.0) then
         do j=1,jmax
           do i=1,imax
             kbath(i,j) = min(kfond, kbath(i,j)*kmax)
           enddo
         enddo
       endif

!==============================================================================
!-- Fermeture des bordures du bassin (rangees et colonnes fantomes)
!==============================================================================
!-- Sauvegarde et fermeture des bordures S et N
       do i=1,imax
         kbath1(i)    = kbath(i,1)
         kbath(i,1)   = 0     ! bordure SUD solide
         kbath2(i)    = kbath(i,jmax)
         kbath(i,jmax)= 0     ! bordure NORD solide
       enddo
!-- Fermeture des bordures E et W
       do j=1,jmax
         kbath(1,j)   = 0     ! bordure OUEST solide
         kbath(imax,j)= 0     ! bordure EST solide
       enddo
!-- Fermeture a la jonction WW/AA (hors equateur)
       if (ltest.ge.2) then
         do i=1,imax
           if (jsep(i).ne.jeq) kbath(i,jsep(i)-1) = 0
         enddo
       endif

!==============================================================================
!-- Verification de la consistance bathymetrique (BATHY >= 1)
!   Un point ne peut pas etre plus profond que tous ses voisins
!   (serait isole sans echange horizontal) -> on reduit kbath au max des voisins.
!==============================================================================
#if ( BATHY >= 1 )
       do i=2,imax-1
         do j=2,jmax-1
           test_wet(1) = min(kbath(i-1,j+1),kbath(i,j+1),kbath(i-1,j))   ! NW
           test_wet(2) = min(kbath(i,j+1),kbath(i+1,j+1),kbath(i+1,j))   ! NE
           test_wet(3) = min(kbath(i+1,j),kbath(i+1,j-1),kbath(i,j-1))   ! SE
           test_wet(4) = min(kbath(i,j-1),kbath(i-1,j-1),kbath(i-1,j))   ! SW
           if (kbath(i,j).gt.test_wet(1) .and.                          &
               kbath(i,j).gt.test_wet(2) .and.                          &
               kbath(i,j).gt.test_wet(3) .and.                          &
               kbath(i,j).gt.test_wet(4)) then
             kbath_bas    = maxval(test_wet(:))
             kbath(i,j)   = kbath_bas
           endif
         enddo
       enddo
#endif

      end subroutine read_bathymetry_and_vscale

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: [build_ocean_mask]
!
!>     @brief [G4a] Construction des masques terre/mer tms, tmu et des profondeurs hs, hu (agnostique a la grille).
!
!      DESCRIPTION:
!>     Pour chaque point (i,j), determine si le point est oceanique (kbath>0) et construit les masques tms (pts scalaires S)
!!     et tmu (pts vitesse U), les profondeurs hs, hu, et les indices de fond kniv, kfs, kfu. Applique les halos de
!!     periodicite et les correspondances au detroit de Bering.
!!
!!     ENTREES :
!!       nflag, ltest                       : mode et type de bassin
!!       ks1, ks2, ims1, ims2, js1, js2     : bornes actives verticales et horizontales
!!       jcl1, jcl2, imax, jmax, kmax, kmaxp1 : bornes cycliques, dimensions, indice fictif sous le fond
!!       iberp, ibera, iberpm, iberam, jberp, jbera, jberpm, jberam : indices du detroit de Bering
!!       ipt0v(10), jpt0v(10), npt0v        : coordonnees et nombre d'iles sans vitesse
!!       zw(kmax+1), unszw(kmax)            : profondeurs des interfaces et 1/zw
!!       kbath(imax,jmax)                   : nombre de niveaux actifs (inout : halos remplis)
!!
!!     SORTIES :
!!       tms, tmu(i,j,k)  : masques scalaire S et vitesse U (1.0=mer, 0.0=terre)
!!       hs, hu(i,j)      : profondeurs d'eau aux pts S et U (m)
!!       hux, huy(i,j)    : profondeurs inter-U en x et y (min des voisins)
!!       unshu(i,j)       : 1/hu(i,j)
!!       kniv(i,j,-1:1)   : indice vertical du fond (-1=pt U, 0=pt S, +1=max 4 coins)
!!       kfs, kfu(i,j)    : niveaux du fond aux pts S et U
!!       is1, is2, iu1, iu2(j) : bornes i actives (scalaire et vitesse) par ligne j
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine build_ocean_mask(                                      &
        nflag, ltest,                                                   &
        ks1, ks2, ims1, ims2, js1, js2, jcl1, jcl2,                     &
        imax, jmax, kmax, kmaxp1,                                       &
        iberp, ibera, iberpm, iberam,                                   &
        jberp, jbera, jberpm, jberam,                                   &
        ipt0v, jpt0v, npt0v,                                            &
        zw, unszw, kbath,                                               &
        tms, tmu,                                                       &
        hs, hu, hux, huy, unshu,                                        &
        kniv, kfs, kfu,                                                 &
        is1, is2, iu1, iu2)

#if ( SHELFMELT == 1 )
       use varsClio2ism_mod, only: tmsism
#endif
       implicit none

!-- Arguments
       integer(ip), intent(in)  :: nflag, ltest
       integer(ip), intent(in)  :: ks1, ks2, ims1, ims2, js1, js2
       integer(ip), intent(in)  :: jcl1, jcl2, imax, jmax, kmax, kmaxp1
       integer(ip), intent(in)  :: iberp, ibera, iberpm, iberam
       integer(ip), intent(in)  :: jberp, jbera, jberpm, jberam
       integer(ip), intent(in)  :: ipt0v(10), jpt0v(10), npt0v
       real(dblp),  intent(in)  :: zw(kmax+1), unszw(kmax)
       integer(ip), intent(inout)  :: kbath(imax,jmax)

       real(dblp), intent(inout) :: tms(imax,jmax,kmax)
       real(dblp), intent(inout) :: tmu(imax,jmax,kmax)
       real(dblp), intent(out)   :: hs(imax,jmax), hu(imax,jmax)
       real(dblp), intent(out)   :: hux(imax,jmax), huy(imax,jmax)
       real(dblp), intent(out)   :: unshu(imax,jmax)
       integer(ip), intent(out)  :: kniv(imax,jmax,-1:1)
       integer(ip), intent(out)  :: kfs(imax,jmax), kfu(imax,jmax)
       integer(ip), intent(out)  :: is1(jmax), is2(jmax)
       integer(ip), intent(out)  :: iu1(jmax), iu2(jmax)

!-- Variables locales
       integer(ip) :: i, j, k, n, km, kkm, km3, kkm3, km1, nn0vit, ii

!==============================================================================
!-- 6.1 : Points scalaires (centre de maille) -- masque tms
!==============================================================================
       do j=js1,js2
         is1(j)=imax
         is2(j)=1
         do i=ims1,ims2
           km=kbath(i,j)
           if (km.gt.0) then
             kkm         = ks1+ks2-km
             kniv(i,j,0) = kkm
             hs(i,j)     = -zw(kkm)
             if (is1(j).eq.imax) is1(j)=i
             is2(j) = i
             do k=kkm,kmax
               tms(i,j,k)=1.0
             enddo
           endif
         enddo
       enddo

!==============================================================================
!-- 6.2 : Correspondance cyclique E<->W pour kbath (ouverture des bords)
!==============================================================================
       if (ltest.ge.1) then
         do j=jcl1,jcl2
           kbath(ims1-1,j) = kbath(ims2,j)  ! halo W
         enddo
       endif
!-- Correspondance Bering
       if (ltest.eq.3 .and. iberp.ne.0) then
         kbath(iberp, jberp) = kbath(iberam,jberam)
         kbath(iberpm,jberp) = kbath(ibera, jberam)
       endif

!==============================================================================
!-- 6.3 : Points de vitesse (coin SW) -- masque tmu
!==============================================================================
       do j=js1,js2
         iu1(j)=imax ; iu2(j)=1
         do i=ims1,ims2
           km3 = min(kbath(i-1,j-1),kbath(i,j),kbath(i-1,j),kbath(i,j-1))
!-- Exclusion des iles sans vitesse
           nn0vit=1
           do n=1,npt0v
             nn0vit=min(nn0vit, abs(i-ipt0v(n))+abs(j-jpt0v(n)))
           enddo
           if (km3.gt.0 .and. nn0vit.eq.1) then
             kkm3        = ks1+ks2-km3
             kniv(i,j,-1)= kkm3
             hu(i,j)     = -zw(kkm3)
             unshu(i,j)  = -unszw(kkm3)
             if (iu1(j).eq.imax) iu1(j)=i
             iu2(j)=i
             do k=kkm3,kmax
               tmu(i,j,k)=1.0
             enddo
           endif
         enddo
       enddo

!==============================================================================
!-- 6.4 : Halos de periodicite (colonnes fantomes W et E)
!==============================================================================
       if (ltest.ge.1) then
         do j=jcl1,jcl2
           kbath(ims2+1,j)   = kbath(ims1,j)
           kniv(ims1-1,j,0)  = kniv(ims2,j,0)
           kniv(ims2+1,j,0)  = kniv(ims1,j,0)
           kniv(ims1-1,j,-1) = kniv(ims2,j,-1)
           kniv(ims2+1,j,-1) = kniv(ims1,j,-1)
           hu(ims1-1,j)      = hu(ims2,j)
           hu(ims2+1,j)      = hu(ims1,j)
           unshu(ims1-1,j)   = unshu(ims2,j)
           unshu(ims2+1,j)   = unshu(ims1,j)
           do k=ks1,ks2
             tms(ims1-1,j,k) = tms(ims2,j,k)
             tms(ims2+1,j,k) = tms(ims1,j,k)
             tmu(ims1-1,j,k) = tmu(ims2,j,k)
             tmu(ims2+1,j,k) = tmu(ims1,j,k)
           enddo
         enddo
       endif

!==============================================================================
!-- 6.5 : Correspondances finales au detroit de Bering (iberp != ibera)
!==============================================================================
       if (ltest.eq.3 .and. iberp.ne.ibera) then
         kbath(ibera, jbera) = kbath(iberpm,jberpm)
         kbath(iberam,jbera) = kbath(iberp, jberpm)
         kniv(iberp, jberp,0)  = kniv(iberam,jberam,0)
         kniv(iberpm,jberp,0)  = kniv(ibera, jberam,0)
         kniv(ibera, jbera,0)  = kniv(iberpm,jberpm,0)
         kniv(iberam,jbera,0)  = kniv(iberp, jberpm,0)
         kniv(ibera, jbera,-1) = kniv(iberp, jberp,-1)
         hu(ibera,jbera)       = hu(iberp,jberp)
         unshu(ibera,jbera)    = unshu(iberp,jberp)
         do k=ks1,ks2
           tms(iberp, jberp,k) = tms(iberam,jberam,k)
           tms(iberpm,jberp,k) = tms(ibera, jberam,k)
           tms(ibera, jbera,k) = tms(iberpm,jberpm,k)
           tms(iberam,jbera,k) = tms(iberp, jberpm,k)
           tmu(ibera, jbera,k) = tmu(iberp, jberp,k)
         enddo
       endif

!==============================================================================
!-- 6.6 : Masque "maximum" aux coins kniv(i,j,+1)
!         = profondeur maximale des 4 pts S voisins ; utilise par la glace
!==============================================================================
       do j=js1,jmax
         do i=ims1,imax
           km1 = max(kbath(i-1,j-1),kbath(i,j),kbath(i-1,j),kbath(i,j-1))
           kniv(i,j,1) = kmaxp1 - km1
         enddo
       enddo
       if (ltest.ge.1) then
         do j=jcl1,jcl2
           kniv(ims1-1,j,1) = kniv(ims2,j,1)
           kniv(ims2+1,j,1) = kniv(ims1,j,1)
         enddo
       endif
       if (ltest.eq.3 .and. iberp.ne.ibera) then
         kniv(ibera,  jbera,1) = kniv(iberp,  jberp,1)
         kniv(ibera+1,jbera,1) = kniv(iberp-1,jberp,1)
         kniv(ibera-1,jbera,1) = kniv(iberp+1,jberp,1)
         do ii=-1,1
           kniv(iberp+ii,jberp+1,1) = kmaxp1  ! fermeture au nord du Bering
         enddo
       endif

!==============================================================================
!-- 6.9 : Indices du fond kfs et kfu
!==============================================================================
       do j=1,jmax
         do i=1,imax
           kfs(i,j) = min(kmax, kniv(i,j,0))
           kfu(i,j) = min(kmax, kniv(i,j,-1))
         enddo
       enddo

#if ( SHELFMELT == 1 )
       tmsism(:,:,:) = tms(:,:,:)  ! copie du masque S pour le modele de calotte
#endif

      end subroutine build_ocean_mask

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: [build_domain_indices]
!
!>     @brief [G4b] Indices de domaine, bornes de flux et coefficient du detroit de Bering (agnostique a la grille).
!
!      DESCRIPTION:
!>     A partir des masques tms, tmu et des bornes is1/is2/iu1/iu2 (sorties de G4a), calcule les bornes de flux isf1/isf2,
!!     iuf1/iuf2, les bornes globales ju1/ju2, imu1/imu2, les profondeurs hux/huy, et le coefficient geostrophique du
!!     detroit de Bering. Corrige egalement les bornes vides pour eviter les boucles nulles.
!!
!!     ENTREES :
!!       ltest, ks1, ks2, ims1, ims2, js1, js2, jcl1, jcl2, imax, jmax, kmax : bornes et dimensions
!!       iberp, ibera, jberp, jbera : indices du detroit de Bering
!!       tms, tmu, hs, hu           : masques et profondeurs (sorties G4a)
!!       is1, is2, iu1, iu2(j)      : bornes actives (inout : correction des bornes vides)
!!       fs2cor, gpes, bering_in, zero : Coriolis, pesanteur, coefficient de controle, 0.0
!!
!!     SORTIES :
!!       isf1, isf2(j)    : bornes i pour les flux scalaires
!!       iuf1, iuf2(j)    : bornes i pour les flux de vitesse
!!       ju1, ju2         : bornes j avec au moins un pt U actif
!!       imu1, imu2       : bornes i globales pour la glace
!!       hux, huy(i,j)    : profondeurs inter-U min en x et y
!!       bering           : coefficient geostrophique du detroit de Bering (m/s)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine build_domain_indices(                                  &
        ltest, ks1, ks2, ims1, ims2, js1, js2, jcl1, jcl2,              &
        imax, jmax, kmax,                                               &
        iberp, ibera, jberp, jbera,                                     &
        tms, tmu, hs, hu, is1, is2, iu1, iu2,                           &
        isf1, isf2, iuf1, iuf2,                                         &
        ju1, ju2, imu1, imu2,                                           &
        fs2cor, gpes, bering_in, zero,                                  &
        hux, huy, bering,                                               &
        clio3_out_id)

       implicit none

!-- Arguments
       integer(ip), intent(in)  :: ltest, ks1, ks2, ims1, ims2, js1, js2
       integer(ip), intent(in)  :: jcl1, jcl2, imax, jmax, kmax
       integer(ip), intent(in)  :: iberp, ibera, jberp, jbera
       real(dblp),  intent(in)  :: tms(imax,jmax,kmax), tmu(imax,jmax,kmax)
       real(dblp),  intent(in)  :: hs(imax,jmax), hu(imax,jmax)
       integer(ip), intent(inout)  :: is1(jmax), is2(jmax)
       integer(ip), intent(inout) :: iu1(jmax), iu2(jmax)
       integer(ip), intent(out) :: isf1(jmax), isf2(jmax)
       integer(ip), intent(out) :: iuf1(jmax), iuf2(jmax)
       integer(ip), intent(out) :: ju1, ju2, imu1, imu2
       real(dblp),  intent(in)  :: fs2cor(imax,jmax), gpes, bering_in, zero
       real(dblp), intent(out)  :: hux(imax,jmax), huy(imax,jmax)
       real(dblp), intent(inout):: bering
       integer(ip), intent(in)  :: clio3_out_id

!-- Variables locales
       integer(ip) :: i, j

!==============================================================================
!-- 6.7 : Bornes isf1/isf2 et iuf1/iuf2
!==============================================================================
       do j=js1,1+js2
         isf1(j) = min(is1(j-1),is1(j))
         isf2(j) = max(is2(j-1),is2(j))
       enddo

       do j=1,jmax-1
         iuf1(j) = min(iu1(j),iu1(j+1))
         iuf2(j) = max(iu2(j),iu2(j+1))
         do i=1,imax-1
           huy(i,j) = min(hu(i,j),hu(i,j+1))
           hux(i,j) = min(hu(i,j),hu(i+1,j))
         enddo
       enddo

!==============================================================================
!-- 6.8 : Bornes globales imu1/imu2, ju1/ju2 ; correction des bornes vides
!==============================================================================
       imu1=imax ; imu2=1 ; ju1=jmax ; ju2=1
       do j=1,jmax
         imu1 = min(imu1,iu1(j))
         imu2 = max(imu2,iu2(j))
         if (iu1(j).gt.iu2(j)) then
           iu2(j)=ims2             ! ligne vide : boucle nulle garantie
         else
           if (ju1.eq.jmax) ju1=j
           ju2=j
         endif
         if (is1(j).gt.is2(j))   is2(j)=ims2
         if (isf1(j).gt.isf2(j)) isf2(j)=ims2
         if (iuf1(j).gt.iuf2(j)) iuf2(j)=ims2
       enddo

!==============================================================================
!-- Coefficient geostrophique du detroit de Bering
!   bering = alpha * g * H / (2*f) ; regule le flux Pacifique->Arctique
!==============================================================================
       if (ltest.eq.3 .and. iberp.ne.ibera) then
         bering = bering_in * gpes * hu(iberp,jberp)                    &
                  * 0.5 / fs2cor(iberp,jberp)
         bering = max(zero, bering)
       endif

      end subroutine build_domain_indices

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: [compute_integrals]
!
!>     @brief [G5a] Grandeurs integrales : surfaces, volumes, coins bathymetriques, pentes (agnostique a la grille).
!
!      DESCRIPTION:
!>     Calcule les integrales de surface ctmi (element de surface pondere par le masque), les aires des mailles, le volume
!!     total oceanique, la surface totale pour w (zurfow), les tableaux de coins bathymetriques (n1-4coin, i1-4coin), et les
!!     listes de pentes pour la diffusion isopycnale de Gent-McWilliams (kslp, lslp, ijslp, nxslp, nxyslp).
!!
!!     ENTREES :
!!       ks1, ks2, ims1, ims2, js1, js2, imax, jmax, kmax : bornes et dimensions
!!       is1, is2, isf1, isf2, iu1, iu2(j), ju1, ju2       : bornes actives
!!       tms, tmu(i,j,k)   : masques terre/mer S et U
!!       cmxy(i,j,0:3)     : elements de surface normalises (sortie G2)
!!       hs(i,j), kfs(i,j) : profondeur et niveau du fond au pt S
!!       kfond             : niveau du fond fixe (<-1 = pentes actives)
!!       dx, dy, ncomax, nlpmax, one : taille de maille, tailles max, constante 1.0
!!       clio3_out_id      : unite du fichier de sortie principal
!!
!!     SORTIES :
!!       ctmi(i,j,k,0:1)  : integrale tms*cmxy (0=pt S, 1=pt U)
!!       aire(i,j)         : surface de la maille en 10^12 m2
!!       volz, unsvol      : volume total oceanique normalise et son inverse
!!       zurfow            : surface totale de l'ocean (10^12 m2)
!!       n1-4coin(k), i1-4coin(n,k) : comptes et indices des coins bathymetriques
!!       kslp, lslp, ijslp(nl) : niveau, direction et indice lineaire des pentes
!!       nxslp, nxyslp, nxslp0, nyslp0 : comptes de pentes (x, total, sauvegardes)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine compute_integrals(                                     &
        ks1, ks2, ims1, ims2, js1, js2,                                 &
        imax, jmax, kmax,                                               &
        is1, is2, isf1, isf2, iu1, iu2, ju1, ju2,                       &
        tms, tmu, cmxy, hs, kfs, kfond,                                 &
        dx, dy, ncomax, nlpmax, one,                                    &
        ctmi, aire, volz, unsvol, zurfow,                               &
        n1coin, n2coin, n3coin, n4coin,                                 &
        i1coin, i2coin, i3coin, i4coin,                                 &
        kslp, lslp, ijslp, nxslp, nxyslp,                               &
        nxslp0, nyslp0,                                                 &
        clio3_out_id)

       implicit none

!-- Arguments
       integer(ip), intent(in)  :: ks1, ks2, ims1, ims2, js1, js2
       integer(ip), intent(in)  :: imax, jmax, kmax
       integer(ip), intent(in)  :: is1(jmax), is2(jmax)
       integer(ip), intent(in)  :: isf1(jmax), isf2(jmax)
       integer(ip), intent(in)  :: iu1(jmax), iu2(jmax)
       integer(ip), intent(in)  :: ju1, ju2
       real(dblp),  intent(in)  :: tms(imax,jmax,kmax), tmu(imax,jmax,kmax)
       real(dblp),  intent(in)  :: cmxy(imax,jmax,0:3)
       real(dblp),  intent(in)  :: hs(imax,jmax)
       integer(ip), intent(in)  :: kfs(imax,jmax), kfond, ncomax, nlpmax
       real(dblp),  intent(in)  :: dx, dy, one
       integer(ip), intent(in)  :: clio3_out_id

       real(dblp), intent(out)  :: ctmi(imax,jmax,kmax,0:1)
       real(dblp), intent(out)  :: aire(imax,jmax)
       real(dblp), intent(out)  :: volz, unsvol, zurfow
       integer(ip), intent(out) :: n1coin(kmax), n2coin(kmax)
       integer(ip), intent(out) :: n3coin(kmax), n4coin(kmax)
       integer(ip), intent(out) :: i1coin(ncomax,kmax), i2coin(ncomax,kmax)
       integer(ip), intent(out) :: i3coin(ncomax,kmax), i4coin(ncomax,kmax)
       integer(ip), intent(out) :: kslp(nlpmax), lslp(nlpmax), ijslp(nlpmax)
       integer(ip), intent(out) :: nxslp, nxyslp, nxslp0, nyslp0

!-- Variables locales
       integer(ip) :: i, j, k, n, nl, ij, llm
       real(dblp)  :: xx4tms, ccdxdy

!==============================================================================
!-- 7.1 : Tableaux des coins bathymetriques (n1-4coin, i1-4coin)
!         Un coin = pt S avec exactement 2 voisins actifs parmi les 4 N/S/E/W.
!         Classifie selon la position du pt U voisin le plus proche.
!==============================================================================
       do k=ks1,ks2
         n1coin(k)=0 ; n2coin(k)=0 ; n3coin(k)=0 ; n4coin(k)=0
         do j=js1,js2
           ij=(j-1)*imax
           do i=is1(j),is2(j)
             xx4tms = tms(i-1,j,k)+tms(i+1,j,k)                         &
                     +tms(i,j-1,k)+tms(i,j+1,k)
             if (xx4tms.eq.2.0d0 .and. tms(i,j,k).eq.one) then
               if (tmu(i,j,k).eq.one) then
                 n=1+n1coin(k)
                 if (n.gt.ncomax) goto 1710
                 i1coin(n,k)=i+ij ; n1coin(k)=n
               elseif (tmu(i+1,j,k).eq.one) then
                 n=1+n2coin(k)
                 if (n.gt.ncomax) goto 1710
                 i2coin(n,k)=i+ij ; n2coin(k)=n
               elseif (tmu(i,j+1,k).eq.one) then
                 n=1+n3coin(k)
                 if (n.gt.ncomax) goto 1710
                 i3coin(n,k)=i+ij ; n3coin(k)=n
               elseif (tmu(i+1,j+1,k).eq.one) then
                 n=1+n4coin(k)
                 if (n.gt.ncomax) goto 1710
                 i4coin(n,k)=i+ij ; n4coin(k)=n
               endif
             endif
           enddo
         enddo
       enddo

!==============================================================================
!-- 7.3 : Listes des pentes bathymetriques (kslp, lslp, ijslp)
!         Transition verticale dans tms : pt actif a k, voisin inactif a k-1.
!==============================================================================
       nl=0
!-- Direction X
       do k=ks1+1,ks2
         do j=js1,js2
           do i=is1(j),is2(j)+1
             if (tms(i-1,j,k).ne.one .or. tms(i,j,k).ne.one) cycle
             if (tms(i-1,j,k-1).eq.one .and. tms(i,j,k-1).eq.0.) then
               nl=nl+1 ; if (nl.gt.nlpmax) goto 1730
               ijslp(nl)=i-1+(j-1)*imax ; kslp(nl)=k ; lslp(nl)=1
             endif
             if (tms(i-1,j,k-1).eq.0. .and. tms(i,j,k-1).eq.one) then
               nl=nl+1 ; if (nl.gt.nlpmax) goto 1730
               ijslp(nl)=i+(j-1)*imax   ; kslp(nl)=k ; lslp(nl)=-1
             endif
           enddo
         enddo
       enddo
       nxslp=nl
!-- Direction Y
       do k=ks1+1,ks2
         do j=js1,js2+1
           do i=isf1(j),isf2(j)
             if (tms(i,j-1,k).ne.one .or. tms(i,j,k).ne.one) cycle
             if (tms(i,j-1,k-1).eq.one .and. tms(i,j,k-1).eq.0.) then
               nl=nl+1 ; if (nl.gt.nlpmax) goto 1730
               ijslp(nl)=i+(j-2)*imax ; kslp(nl)=k ; lslp(nl)=imax
             endif
             if (tms(i,j-1,k-1).eq.0. .and. tms(i,j,k-1).eq.one) then
               nl=nl+1 ; if (nl.gt.nlpmax) goto 1730
               ijslp(nl)=i+(j-1)*imax ; kslp(nl)=k ; lslp(nl)=-imax
             endif
           enddo
         enddo
       enddo
       nxyslp=nl
       nxslp0=nxslp ; nyslp0=nl-nxslp
!-- Desactivation des pentes si fond plat ou bathymetrie simplifiee
       if (kfond.gt.-2) then
         nxslp=0 ; nxyslp=0
       endif

!==============================================================================
!-- 7.5 : Surfaces ctmi, aire, volz, unsvol, zurfow
!==============================================================================
       do llm=0,1
         do k=1,kmax
           ctmi(:,:,k,llm)=0.
         enddo
       enddo

       do k=1,kmax
         do j=js1,js2
           do i=is1(j),is2(j)
             ctmi(i,j,k,0)=tms(i,j,k)*cmxy(i,j,0)  ! surface S normalisee
           enddo
         enddo
       enddo

       do k=1,kmax
         do j=ju1,ju2
           do i=iu1(j),iu2(j)
             ctmi(i,j,k,1)=tmu(i,j,k)*cmxy(i,j,3)  ! surface U normalisee (coin SW)
           enddo
         enddo
       enddo

       volz=0.0
       do j=js1,js2
         do i=is1(j),is2(j)
           volz=volz+cmxy(i,j,0)*hs(i,j)
         enddo
       enddo
       unsvol=1.0/volz

       ccdxdy=dx*dy/1.0d+12
       do j=1,jmax
         do i=1,imax
           aire(i,j)=ctmi(i,j,ks2,0)*ccdxdy  ! surface en 10^12 m2
         enddo
       enddo

       zurfow=0.0
       do j=js1,js2
         do i=is1(j),is2(j)
           zurfow=zurfow+aire(i,j)*tms(i,j,ks2)
         enddo
       enddo

       return

!-- Gestion des erreurs
 1710 continue
       write(clio3_out_id,'(A,2I5)')                                    &
         'Stop in compute_integrals : ncomax too small ! (k,ncomax)= ', &
         k, ncomax
       stop
 1730 continue
       write(clio3_out_id,'(A,4I5)')                                    &
         'Stop in compute_integrals : nlpmax too small ! (i,j,k,nlpmax)= ', &
         i, j, k, nlpmax
       stop

      end subroutine compute_integrals

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: [compute_surface_physics]
!
!>     @brief [G5b] Physique de surface : masques glace, absorption solaire, fichiers de sortie (agnostique a la grille).
!
!      DESCRIPTION:
!>     Calcule le masque de surface pour la glace (zindfa), les coefficients de diffusion horizontale (dfhu, dfhv),
!!     l'epaisseur optique nuageuse (tauc), l'absorption solaire niveau par niveau (reslum), finalise les masques msks/msku
!!     et angle, et ecrit les fichiers lat.dat/long.dat/mask.dat.
!!
!!     ENTREES :
!!       nflag, imax, jmax, kmax, ks2 : mode et dimensions
!!       is1, is2(j), js1, js2, kfs(i,j) : bornes actives et niveau du fond
!!       tms(i,j,k), covrai(i,j)      : masque S et sin(lat geographique)
!!       zw(kmax+1), dts(kmax), unsdz(kmax) : profils verticaux
!!       rho0, cpo, radian, uvdif, ren, gridsz : constantes physiques et de diffusion
!!       zlatt, zlont(i,j)            : latitude/longitude U en radians
!!       mouchard_id                 : unite du fichier de debogage
!!
!!     SORTIES :
!!       zindfa(i,j)       : masque de surface pour la glace (= tms(i,j,ks2))
!!       dfhu, dfhv(i,j)   : coefficients de diffusion horizontale glace en x et y
!!       tauc(i,j)         : epaisseur optique nuageuse (Chou et al. 1981)
!!       reslum(i,j,k)     : fraction du rayonnement solaire restant au niveau k
!!       msks, msku(i,j)   : masques de surface finalises (inout pour msku)
!!       angle(i,j)        : angle mis a zero sur les pts terre (inout)
!!       hs, hu(i,j)       : profondeurs avec halos E/W remplis (inout)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine compute_surface_physics(                               &
        nflag, imax, jmax, kmax,                                        &
        ks2, is1, is2, js1, js2, kfs,                                   &
        tms, covrai, zw, dts, unsdz, rho0, cpo,                         &
        radian, uvdif, ren, gridsz,                                     &
        hs, hu, angle, msku, msks,                                      &
        zlatt, zlont,                                                   &
        mouchard_id,                                                    &
        zindfa, dfhu, dfhv, tauc, reslum)

       use newunit_clio_mod, only: clio3_out_id
       implicit none

!-- Arguments
       integer(ip), intent(in)  :: nflag, imax, jmax, kmax, ks2, js1, js2
       integer(ip), intent(in)  :: is1(jmax), is2(jmax), kfs(imax,jmax)
       real(dblp),  intent(in)  :: tms(imax,jmax,kmax), covrai(imax,jmax)
       real(dblp),  intent(in)  :: zw(kmax+1), dts(kmax), unsdz(kmax)
       real(dblp),  intent(in)  :: rho0, cpo, radian, uvdif, ren, gridsz
       real(dblp),  intent(in)  :: zlatt(imax,jmax), zlont(imax,jmax)
       integer(ip), intent(in)  :: mouchard_id

       real(dblp), intent(inout)  :: hs(imax,jmax), hu(imax,jmax)
       real(dblp), intent(inout)  :: angle(imax,jmax)
       integer(ip), intent(inout) :: msku(imax,jmax)
       integer(ip), intent(out)   :: msks(imax,jmax)
       real(dblp),  intent(out)   :: zindfa(imax,jmax)
       real(dblp),  intent(out)   :: dfhu(imax,jmax), dfhv(imax,jmax)
       real(dblp),  intent(out)   :: tauc(imax,jmax)
       real(dblp),  intent(out)   :: reslum(imax,jmax,0:kmax+1)

!-- Variables locales
       integer(ip) :: i, j, k, ip1, jp1, indx, indxp1, klim
       integer(ip) :: typeaux_dat_id, lat_dat_id, lon_dat_id, mask_id
       integer(ip) :: ntypcn(imax,jmax)
       real(dblp)  :: rrro(imax,jmax), dd1o(imax,jmax), dd2o(imax,jmax)
       real(dblp)  :: ah, alat, clat, dl, dr, znivsp, znivin

       real(dblp), parameter :: tauco(20) =                             &
         (/6.6d0,6.6d0,7.0d0,7.2d0,7.1d0,6.8d0,6.5d0,6.6d0,7.1d0,7.6d0, &
           6.6d0,6.1d0,5.6d0,5.5d0,5.8d0,5.8d0,5.6d0,5.6d0,5.6d0,5.6d0/)
       real(dblp), parameter :: zrr(5) = (/0.58d0,0.62d0,0.67d0,0.77d0,0.78d0/)
       real(dblp), parameter :: zd1(5) = (/0.35d0,0.60d0,1.00d0,1.50d0,1.40d0/)
       real(dblp), parameter :: zd2(5) = (/23.0d0,20.0d0,17.0d0,14.0d0,7.90d0/)

!==============================================================================
!-- 9.1 : Masque de surface pour la glace et coefficients de diffusion horizontale
!==============================================================================
       do j=1,jmax
         do i=1,imax
           zindfa(i,j)=tms(i,j,ks2)
         enddo
       enddo

       ah=(uvdif/ren)*gridsz   ! coefficient de diffusion horizontale glace (m2/s)
       do j=1,jmax-1
         jp1=j+1
         do i=1,imax
           ip1=(i+1)-(imax-2)*(i/imax)
           dfhu(i,j)=tms(i,j,ks2)*tms(ip1,j, ks2)*ah  ! diffusion en x
           dfhv(i,j)=tms(i,j,ks2)*tms(i,  jp1,ks2)*ah  ! diffusion en y
         enddo
       enddo

!==============================================================================
!-- 9.2 : Epaisseur optique nuageuse tauc (Chou et al. 1981)
!         Interpolation lineaire dans tauco selon la latitude (bandes de 5 deg)
!==============================================================================
       if (nflag.ne.3) then
         do j=js1,js2
           do i=is1(j),is2(j)
             alat    = asin(covrai(i,j))/radian
             clat    = (95.0d0-alat)/10.0d0
             indx    = 1+int(clat)
             indxp1  = indx+1
             dl      = clat-int(clat)
             dr      = 1.0d0-dl
             tauc(i,j) = dr*tauco(indx)+dl*tauco(indxp1)
           enddo
         enddo
       endif

!==============================================================================
!-- 10 : Absorption solaire par niveau (modele a 2 longueurs d'onde,
!        Paulson & Simpson 1977)
!==============================================================================
!-- Lecture du fichier des types d'eau
       open(newunit=typeaux_dat_id,                                    &
            file='inputdata/clio/typeaux.dat', status='old')
       do j=1,jmax
         read(typeaux_dat_id,*) (ntypcn(i,j),i=1,imax)
       enddo
       close(typeaux_dat_id)

!-- Affectation des proprietes optiques par type d'eau
       do j=1,jmax
         do i=1,imax
           rrro(i,j)=zrr(ntypcn(i,j))
           dd1o(i,j)=zd1(ntypcn(i,j))
           dd2o(i,j)=zd2(ntypcn(i,j))
         enddo
       enddo

!-- Calcul de reslum(i,j,k) : fraction du rayonnement solaire a l'interface k
       do j=1,jmax
         do i=1,imax
           do k=1,kmax
             reslum(i,j,k)=0.0
           enddo
           znivin=0.0
           klim=max(11,kfs(i,j))
           do k=klim,kmax
             znivsp = rrro(i,j)*exp(max(-500.d0,zw(k+1)/dd1o(i,j)))     &
                     +(1.0d0-rrro(i,j))*exp(max(-500.d0,zw(k+1)/dd2o(i,j)))
             reslum(i,j,k) = (znivsp-znivin)*tms(i,j,k)                 &
                             *dts(k)/dts(ks2)*unsdz(k)/(rho0*cpo)
             znivin=znivsp
           enddo
         enddo
       enddo

!==============================================================================
!-- Finalisation des masques de surface et de l'angle de rotation
!==============================================================================
       do j=1,jmax
         hs(1,j)   =hs(imax-1,j) ; hs(imax,j)   =hs(2,j)
         hu(1,j)   =hu(imax-1,j) ; hu(imax,j)   =hu(2,j)
         do i=1,imax
           msks(i,j)=1 ; if (hs(i,j).eq.0) msks(i,j)=0
           msku(i,j)=1 ; if (hu(i,j).eq.0) msku(i,j)=0
           angle(i,j)=angle(i,j)*float(msku(i,j))  ! angle nul sur les pts terre
         enddo
       enddo

!==============================================================================
!-- Ecriture des fichiers de coordonnees pour le post-traitement
!==============================================================================
       open(newunit=lat_dat_id,  file='outputdata/ocean/lat.dat')
       open(newunit=lon_dat_id,  file='outputdata/ocean/long.dat')
       open(newunit=mask_id,     file='outputdata/ocean/mask.dat')
       do j=1,jmax
#if ( HRCLIO == 0 )
         write(lat_dat_id, '(122(F10.5))') (zlatt(i,j),i=1,imax)
         write(lon_dat_id, '(122(F10.5))') (zlont(i,j),i=1,imax)
         write(mask_id,    '(122(I1))'  ) (msks(i,j), i=1,imax)
#elif ( HRCLIO == 1 )
         write(lat_dat_id, '(242(F10.5))') (zlatt(i,j),i=1,imax)
         write(lon_dat_id, '(242(F10.5))') (zlont(i,j),i=1,imax)
         write(mask_id,    '(242(I1))'  ) (msks(i,j), i=1,imax)
#endif
       enddo
       write(lat_dat_id,*) ; write(lon_dat_id,*) ; write(mask_id,*)
       close(lat_dat_id)   ; close(lon_dat_id)   ; close(mask_id)

      end subroutine compute_surface_physics

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: [define_ice_shelves]
!
!>     @brief [G5c] Plateformes de glace (ice shelves) et zones icebergs (agnostique a la grille).
!
!      DESCRIPTION:
!>     Initialise le masque tmics des ice shelves (surface effective de contact ice shelf/ocean en m2) pour les 8
!!     plateformes antarctiques nommees, et calcule les surfaces des zones de fonte des icebergs Nord (areicebn) et Sud
!!     (areicebs).
!!
!!     NOTA : Les indices geographiques codes en dur (i=93..103, j=2, etc.) sont specifiques a la resolution CL15/CL30.
!!            Pour une autre resolution, remplacer ces indices ou les lire depuis un fichier de configuration.
!!
!!     ENTREES :
!!       imax, jmax, ks2   : dimensions et indice de surface
!!       tms(i,j,ks2)      : masque de surface (1=mer, 0=terre)
!!       dxs1, dxs2(i,j)   : longueurs des cotes N et E de la maille S (m)
!!       area(i,j)         : surface de la maille (m2)
!!       mouchard_id       : unite du fichier de debogage
!!
!!     SORTIES :
!!       tmics(i,j)        : surface effective de contact ice shelf (m2 ; 0 si absent)
!!       iicebern1/2, jicebern1/2 : bornes de la zone icebergs Nord
!!       iicebers1/2, jicebers1/2 : bornes de la zone icebergs Sud
!!       areicebn, areicebs : surfaces oceaniques des zones icebergs Nord et Sud (m2)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine define_ice_shelves(                                    &
        imax, jmax, ks2,                                                &
        tms, dxs1, dxs2, area,                                          &
        tmics,                                                          &
        iicebern1, iicebern2, jicebern1, jicebern2,                     &
        iicebers1, iicebers2, jicebers1, jicebers2,                     &
        areicebn, areicebs,                                             &
        mouchard_id)

       implicit none

!-- Arguments
       integer(ip), intent(in)  :: imax, jmax, ks2, mouchard_id
       real(dblp),  intent(in)  :: tms(imax,jmax,ks2)
       real(dblp),  intent(in)  :: dxs1(imax,jmax), dxs2(imax,jmax)
       real(dblp),  intent(in)  :: area(imax,jmax)

       real(dblp),  intent(out) :: tmics(imax,jmax)
       integer(ip), intent(out) :: iicebern1, iicebern2, jicebern1, jicebern2
       integer(ip), intent(out) :: iicebers1, iicebers2, jicebers1, jicebers2
       real(dblp),  intent(out) :: areicebn, areicebs

!-- Variables locales
       integer(ip) :: i, j
       real(dblp)  :: effecta

!==============================================================================
!-- Zones de fonte des icebergs
!   areicebn/s = surface oceanique active dans la zone Nord/Sud (m2)
!==============================================================================
!-- Zone icebergs NORD (Atlantique Nord, ~50-55N)
       iicebern1=98 ; iicebern2=115 ; jicebern1=46 ; jicebern2=55
       areicebn=0.0
       do i=iicebern1,iicebern2
         do j=jicebern1,jicebern2
           areicebn=areicebn+area(i,j)*tms(i,j,ks2)
         enddo
       enddo

!-- Zone icebergs SUD (Antarctique circumpolaire)
       iicebers1=2 ; iicebers2=121 ; jicebers1=1 ; jicebers2=9
       areicebs=0.0
       do i=iicebers1,iicebers2
         do j=jicebers1,jicebers2
           areicebs=areicebs+area(i,j)*tms(i,j,ks2)
         enddo
       enddo

       write(mouchard_id,*) 'surface iceberg melting N,S',              &
                             areicebn, areicebs

!==============================================================================
!-- Plateformes de glace (ice shelves)
!   effecta = epaisseur effective de la face verticale de contact (m)
!   tmics(i,j) = tms(i,j,ks2) * dxs1_ou_dxs2(i,j) * effecta
!
!   Plateformes nommees (indices specifiques CL15/CL30) :
!     1. Ross Ice Shelf (i=93..103, j=2)     2. Est Weddell 1 (i=104..109, j=3..4)
!     3. Est Weddell 2  (i=110..117, j=4..5) 4. Larsen        (i=93, j=4..5, dxs2)
!     5. Amery          (i=14..17, j=5)      6. Ross (suite)  (i=48..61, j=2)
!     7. Getz           (i=73..77, j=3)      8. George VI     (i=84..87, j=4)
!==============================================================================
       effecta=5.0d3   ! epaisseur effective de face verticale (m)

       tmics(:,:)=0.0  ! initialisation -- pas de plateforme par defaut

!-- 1. Ross Ice Shelf
       do i=93,103
         tmics(i,2)=tms(i,2,ks2)*dxs1(i,2)*effecta
       enddo
       do i=93,96
         tmics(i,3)=tms(i,3,ks2)*dxs1(i,2)*effecta
       enddo

!-- 2. Est Weddell 1
       do i=104,106
         tmics(i,3)=tms(i,3,ks2)*dxs1(i,3)*effecta
       enddo
       do i=107,109
         tmics(i,4)=tms(i,4,ks2)*dxs1(i,4)*effecta
       enddo

!-- 3. Est Weddell 2
       do i=110,114
         tmics(i,4)=tms(i,4,ks2)*dxs1(i,4)*effecta
       enddo
       do i=115,117
         tmics(i,5)=tms(i,5,ks2)*dxs1(i,5)*effecta
       enddo

!-- 4. Larsen (cote est de la maille, dxs2)
       do j=4,5
         tmics(93,j)=tms(93,j,ks2)*dxs2(93,j)*effecta
       enddo

!-- 5. Amery
       do i=14,17
         tmics(i,5)=tms(i,5,ks2)*dxs1(i,5)*effecta
       enddo

!-- 6. Ross (suite)
       do i=48,61
         tmics(i,2)=tms(i,2,ks2)*dxs1(i,2)*effecta
       enddo

!-- 7. Getz
       do i=73,77
         tmics(i,3)=tms(i,3,ks2)*dxs1(i,3)*effecta
       enddo

!-- 8. George VI
       do i=84,87
         tmics(i,4)=tms(i,4,ks2)*dxs1(i,4)*effecta
       enddo

      end subroutine define_ice_shelves

      end module defgrid_kernel_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
