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
!      MODULE: [defgrid_compute_mod]
!
!>     @brief Bloc GRILLE-DEPENDANT de l'initialisation : coordonnees et coefficients metriques.
!
!      DESCRIPTION:
!>     Regroupe les deux seules sous-routines de defgrid qui dependent du TYPE de grille (bi-domaine spherique WW + spherique
!!     tournee AA) : compute_grid_coords (G1) et compute_metric_coeffs (G2). Pour porter le modele sur une autre grille, on
!!     remplace CE module en bloc, sans toucher au noyau agnostique (defgrid_kernel_mod) ni a l'orchestrateur (defgrid_mod).
!!
!!     G1 produit les coordonnees lon/lat, les coins de mailles, l'angle de rotation grille->geo, et les metriques locales
!!     h1/h2/d1d2/d2d1/dx1/dx2 consommees par G2 et par le noyau. G2 en derive les coefficients metriques cmx/cmy/smx/smy/
!!     cmxy/smxy aux 4 positions, les derivees cmxdy/cmydx, le facteur de Coriolis fs2cor et le sinus de latitude covrai.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module defgrid_compute_mod

       use global_constants_mod, only: dblp => dp, ip

       implicit none

       private

       public :: compute_grid_coords
       public :: compute_metric_coeffs

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: [compute_grid_coords]                                                          *** DEPEND DE LA GRILLE ***
!
!>     @brief [G1] Coordonnees geographiques de la grille bi-domaine WW+AA et metriques locales.
!
!      DESCRIPTION:
!>     Calcule, pour chaque point (i,j), les coordonnees geographiques (lon/lat en degres et radians), les coins des mailles,
!!     l'angle de rotation grille->geo, et les metriques locales h1/h2 necessaires aux sous-routines G2 et G3a.
!!
!!     ENTREES :
!!       ltest          : type de bassin (0=rect, <0=cartesien, >=3=monde reel)
!!       jsepar         : derniere ligne j de la grille WW (avant transition AA)
!!       jeq            : indice j de jonction WW/AA
!!       jsep(imax)     : pour chaque colonne i, indice j de debut de grille AA
!!       jsep2(83:imax) : indices secondaires de separation (specifiques CL15/CL30)
!!       dlong, dlat    : resolution en longitude/latitude grille WW (degres)
!!       xlon1, ylat1   : origine de la grille WW (degres)
!!       dxaj, dyai     : resolution de la grille AA (degres)
!!       xaj1, yai1     : origine de la grille AA (degres)
!!       radian         : pi/180
!!       rterre         : rayon de la Terre (m)
!!       unsrt          : 1/rterre
!!       omega          : vitesse angulaire de la Terre (rad/s)
!!       pi             : 3.14159...
!!       zepsd1         : epsilon pour eviter les divisions par zero
!!       imax, jmax     : dimensions de la grille
!!
!!     SORTIES :
!!       xslon/yslat(i,j)   : longitude/latitude du pt scalaire S (degres)
!!       xulon/yulat(i,j)   : longitude/latitude du pt vitesse U (degres)
!!       zlatt/zlont(i,j)   : latitude/longitude du pt vitesse U (radians)
!!       xsedg/ysedg(i,j,4) : longitudes/latitudes des 4 coins de la maille S (SW,SE,NE,NW)
!!       xuedg/yuedg(i,j,4) : longitudes/latitudes des 4 coins de la maille U
!!       xslonp/yslatp      : coin SW de la maille S, etendu a (imax+1,jmax+1) pour Ferret
!!       xulonp/yulatp      : coin SW de la maille U etendu
!!       angle(i,j)         : angle de rotation grille->geo (radians)
!!       xang1/xang2(i,j)   : composantes cos/sin de la rotation (grille AA uniquement)
!!       dx1/dx2(i,j)       : pas angulaires adimensionnes en x/y (rad)
!!       h1/h2(i,j)         : coefficients metriques x/y au pt S (m)
!!       d1d2/d2d1(i,j)     : d(h1)/dy, d(h2)/dx au pt S (m/rad)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine compute_grid_coords(                                   &
        ltest, jsepar, jeq, jsep, jsep2,                                &
        dlong, dlat, xlon1, ylat1,                                      &
        dxaj, dyai, xaj1, yai1,                                         &
        radian, rterre, unsrt, omega, pi, zepsd1,                       &
        imax, jmax,                                                     &
        xslon, yslat, xulon, yulat,                                     &
        zlatt, zlont,                                                   &
        xsedg, ysedg, xuedg, yuedg,                                     &
        xslonp, yslatp, xulonp, yulatp,                                 &
        angle, xang1, xang2,                                            &
        dx1, dx2, h1, h2, d1d2, d2d1)

       implicit none

!-- Declaration des arguments
       integer(ip), intent(in) :: ltest, jsepar, jeq, imax, jmax
       integer(ip), intent(in) :: jsep(imax)
       integer(ip), intent(in) :: jsep2(83:imax)
       real(dblp),  intent(in) :: dlong, dlat, xlon1, ylat1
       real(dblp),  intent(in) :: dxaj, dyai, xaj1, yai1
       real(dblp),  intent(in) :: radian, rterre, unsrt, omega, pi, zepsd1

       real(dblp), intent(out) :: xslon(imax,jmax), yslat(imax,jmax)
       real(dblp), intent(out) :: xulon(imax,jmax), yulat(imax,jmax)
       real(dblp), intent(out) :: zlatt(imax,jmax), zlont(imax,jmax)
       real(dblp), intent(out) :: xsedg(imax,jmax,4), ysedg(imax,jmax,4)
       real(dblp), intent(out) :: xuedg(imax,jmax,4), yuedg(imax,jmax,4)
       real(dblp), intent(out) :: xslonp(imax+1,jmax+1), yslatp(imax+1,jmax+1)
       real(dblp), intent(out) :: xulonp(imax+1,jmax+1), yulatp(imax+1,jmax+1)
       real(dblp), intent(out) :: angle(imax,jmax)
       real(dblp), intent(out) :: xang1(imax,jmax), xang2(imax,jmax)
       real(dblp), intent(out) :: dx1(imax,jmax), dx2(imax,jmax)
       real(dblp), intent(out) :: h1(imax,jmax),  h2(imax,jmax)
       real(dblp), intent(out) :: d1d2(imax,jmax), d2d1(imax,jmax)

!-- Variables locales
       integer(ip) :: i, j, ipm, jpm
       real(dblp)  :: ysdeg, yudeg, yurad, yu0rad, yu2rad, ysrad, xsrad
       real(dblp)  :: ccsysw, sscysw, ccsyuw, sscyuw, ssnyuw, ssnysw
       real(dblp)  :: ccsysa, sscysa, ccsyua, sscyua, ssnyua
       real(dblp)  :: ccmxdy, ccfcor, ccfcu8, ccmydx
       real(dblp)  :: xsdeg, xudeg, xurad
       real(dblp)  :: xu0m90, xudeg_aa, yyalim, yfac
       real(dblp)  :: szphi, zlam
       real(dblp)  :: csxua(jmax), csxsa(jmax)

!==============================================================================
!-- Initialisation a 90 deg (valeur sentinelle ; sera ecrasee par les vraies valeurs)
!==============================================================================
       do j=1,jmax
         do i=1,imax
           xslon(i,j)=90.0 ; yslat(i,j)=90.0
           xulon(i,j)=90.0 ; yulat(i,j)=90.0
           xsedg(i,j,1)=90.0 ; xsedg(i,j,2)=90.0
           xsedg(i,j,3)=90.0 ; xsedg(i,j,4)=90.0
           ysedg(i,j,1)=90.0 ; ysedg(i,j,2)=90.0
           ysedg(i,j,3)=90.0 ; ysedg(i,j,4)=90.0
           xuedg(i,j,1)=90.0 ; xuedg(i,j,2)=90.0
           xuedg(i,j,3)=90.0 ; xuedg(i,j,4)=90.0
           yuedg(i,j,1)=90.0 ; yuedg(i,j,2)=90.0
           yuedg(i,j,3)=90.0 ; yuedg(i,j,4)=90.0
           angle(i,j)=0.0
           xang1(i,j)=0.0
           xang2(i,j)=1.0
         enddo
       enddo
       do j=1,jmax+1
         do i=1,imax+1
           xslonp(i,j)=90.0 ; yslatp(i,j)=90.0
           xulonp(i,j)=90.0 ; yulatp(i,j)=90.0
         enddo
       enddo

!==============================================================================
!-- 3.1 : GRILLE SPHERIQUE CLASSIQUE (WW) : lon/lat standard
!==============================================================================

       do j=1,jsepar
         ysdeg = ylat1 + dlat * DFLOAT(j-1)
         yudeg = ysdeg - 0.5 * dlat
         yurad = (ysdeg - 0.5 * dlat) * radian

         if (ltest.lt.0) then
!-- Grille cartesienne (cas test) : metriques constantes a 45N
           yu0rad = 45.0 * radian
           ccsysw = 1.0  ; sscysw = 1.0
           ccsyuw = 1.0  ; sscyuw = 1.0
           ssnysw = sin(yu0rad) ; ssnyuw = sin(yu0rad)
           ccmxdy = 0.0
           ccfcor = omega * sin(yu0rad)
           ccfcu8 = -0.25 * omega * cos(yu0rad)
         else
!-- Grille spherique reelle
           ccsysw = cos(ysdeg * radian)
           sscysw = 1.0 / ccsysw
           ccsyuw = cos(yurad)
           sscyuw = 1.0 / ccsyuw
           ssnyuw = sin(yurad)
           ssnysw = sin(ysdeg*radian)
           ccmxdy = -ssnyuw * unsrt
           ccfcor = omega * ssnyuw
           ccfcu8 = -0.25 * omega * ccsyuw
         endif

         do i=1,imax
           zlatt(i,j)  = yurad
           yslat(i,j)  = ysdeg
           ysedg(i,j,1)= yudeg
           ysedg(i,j,2)= yudeg
           ysedg(i,j,3)= yudeg + dlat
           ysedg(i,j,4)= yudeg + dlat
           yulat(i,j)  = yudeg
           yuedg(i,j,1)= ysdeg - dlat
           yuedg(i,j,2)= ysdeg - dlat
           yuedg(i,j,3)= ysdeg
           yuedg(i,j,4)= ysdeg
           xsdeg = xlon1 + dlong * DFLOAT(i-1)
           xurad = (xsdeg - 0.5 * dlong) * radian
           xudeg = xsdeg - 0.5 * dlong
           xslon(i,j)  = xsdeg
           xsedg(i,j,1)= xudeg
           xsedg(i,j,2)= xudeg + dlong
           xsedg(i,j,3)= xudeg + dlong
           xsedg(i,j,4)= xudeg
           xulon(i,j)  = xudeg
           if ((xsdeg-dlong).le.0.0d0) then
             xuedg(i,j,1)= xsdeg - dlong + 360.0d0
           else
             xuedg(i,j,1)= xsdeg - dlong
           endif
           xuedg(i,j,2)= xsdeg
           xuedg(i,j,3)= xsdeg
           if ((xsdeg-dlong).le.0.0d0) then
             xuedg(i,j,4)= xsdeg - dlong + 360.0d0
           else
             xuedg(i,j,4)= xsdeg - dlong
           endif
           zlont(i,j) = xurad
!-- Metriques pour la dynamique de la glace (section 4)
           dx1(i,j)  = dlong * radian
           yu2rad    = (ysdeg + 0.5*dlat) * radian
           dx2(i,j)  = sin(yu2rad) - sin(yurad)
           h1(i,j)   = rterre * ccsysw
           h2(i,j)   = rterre / ccsysw
           d1d2(i,j) = -rterre * sin(ysdeg*radian) * sscysw
           d2d1(i,j) = 0.0
         enddo

!-- Angle de rotation (quasi nul en grille WW)
         do i=1,imax-1
           angle(i,j)=atan2((yulat(i+1,j)-yulat(i,j)),                  &
                            (xulon(i+1,j)-xulon(i,j)))
         enddo
         angle(imax,j)=angle(imax-1,j)

       enddo  ! fin boucle j grille WW

!==============================================================================
!-- 3.3 : GRILLE SPHERIQUE TOURNEE (AA) : Atlantique Nord + Arctique
!==============================================================================

       if (ltest.ge.3) then

         xu0m90 = xaj1 - 1.5 * dxaj - 90.
         do j=jeq,jmax
           xudeg_aa  = xu0m90 + dxaj * DFLOAT(j)
           csxua(j)  = sin(xudeg_aa * radian)
           csxsa(j)  = sin((xudeg_aa + 0.5*dxaj) * radian)
         enddo
         yyalim = 90.0 - abs(dyai)

         do i=1,imax
           ysdeg = yai1 + dyai * DFLOAT(i-1)
           if (abs(ysdeg).ge.yyalim) then
             ccsysa=1. ; yurad=0. ; ysrad=0. ; ccsyua=1. ; ssnyua=0.
           else
             ccsysa = cos(ysdeg * radian)
             yurad  = (ysdeg - 0.5 * dyai) * radian
             ysrad  = ysdeg * radian
             ccsyua = cos(yurad)
           endif
           sscysa = 1.0 / ccsysa
           sscyua = 1.0 / ccsyua
           ccmydx = sin(yurad) * unsrt

           do j=jsep(i),jmax
             dx2(i,j)  = dlat * radian
             yu2rad    = (ysdeg + 0.5*dyai) * radian
             dx1(i,j)  = -sin(yu2rad) + sin(yurad)
             h1(i,j)   = rterre / ccsysa
             h2(i,j)   = rterre * ccsysa
             d1d2(i,j) = 0.0
             d2d1(i,j) = +rterre * sin(ysdeg*radian) * sscysa

!-- Transformation rotation spherique AA -> lon/lat geographique
             xsrad = (xaj1 + dxaj*DFLOAT(j-1)) * radian
             xurad = (xaj1 - 1.5*dxaj + dxaj*DFLOAT(j)) * radian
             yslat(i,j) = -asin(cos(ysrad)*cos(xsrad)) / radian
             yulat(i,j) = -asin(cos(yurad)*cos(xurad)) / radian
             zlatt(i,j) = -asin(cos(yurad)*cos(xurad))
             szphi      =  cos(yurad)*cos(xurad)
             zlam = atan2(cos(ysrad)*sin(xsrad), sin(ysrad))
             xslon(i,j) = (zlam + pi) / radian + 69.0
             zlam = atan2(cos(yurad)*sin(xurad), sin(yurad))
             xulon(i,j) = (zlam + pi) / radian + 69.0
             zlont(i,j) = zlam + pi + 69.0*pi/180.0
!-- Angles de rotation grille AA -> geographique
             xang1(i,j) = (szphi*cos(zlam)) / max(zepsd1,cos(yurad))
             xang2(i,j) = (sin(zlam))       / max(zepsd1,cos(yurad))
           enddo
         enddo ! fin boucle i grille AA

!-- Affectation des coins des mailles en grille AA
         do i=1,imax
           do j=jsep(i),jmax
             ipm=min(i+1,imax) ; jpm=min(j+1,jmax)
             xsedg(i,j,1)=xulon(i  ,j  ) ; xsedg(i,j,2)=xulon(ipm,j  )
             xsedg(i,j,3)=xulon(ipm,jpm) ; xsedg(i,j,4)=xulon(i  ,jpm)
             ysedg(i,j,1)=yulat(i  ,j  ) ; ysedg(i,j,2)=yulat(ipm,j  )
             ysedg(i,j,3)=yulat(ipm,jpm) ; ysedg(i,j,4)=yulat(i  ,jpm)
             ipm=max(i-1,1) ; jpm=max(j-1,1)
             xuedg(i,j,1)=xslon(ipm,jpm) ; xuedg(i,j,2)=xslon(i  ,jpm)
             xuedg(i,j,3)=xslon(i  ,j  ) ; xuedg(i,j,4)=xslon(ipm,j  )
             yuedg(i,j,1)=yslat(ipm,jpm) ; yuedg(i,j,2)=yslat(i  ,jpm)
             yuedg(i,j,3)=yslat(i  ,j  ) ; yuedg(i,j,4)=yslat(ipm,j  )
           enddo
         enddo

!-- Lissage des coordonnees dans la zone de transition (colonnes i=1..82)
         do i=1,82
           yfac = 0.25 / float(jmax - jsep(i))
           jpm  = jsep(i) - 1
           do j=jsep(i),jmax
             xslon(i,j)   = xslon(i,jpm)
             yslat(i,j)   = yslat(i,jpm) + yfac*(j-jpm)
             xsedg(i,j,1) = xsedg(i,jpm,1) ; xsedg(i,j,2) = xsedg(i,jpm,2)
             xsedg(i,j,3) = xsedg(i,jpm,3) ; xsedg(i,j,4) = xsedg(i,jpm,4)
             ysedg(i,j,1) = yslat(i,j)-0.5*yfac ; ysedg(i,j,2) = yslat(i,j)-0.5*yfac
             ysedg(i,j,3) = yslat(i,j)+0.5*yfac ; ysedg(i,j,4) = yslat(i,j)+0.5*yfac
             xulon(i,j)   = xulon(i,jpm)
             yulat(i,j)   = yulat(i,jpm) + yfac*(float(j-jpm)-0.5)
             xuedg(i,j,1) = xuedg(i,jpm,1) ; xuedg(i,j,2) = xuedg(i,jpm,2)
             xuedg(i,j,3) = xuedg(i,jpm,3) ; xuedg(i,j,4) = xuedg(i,jpm,4)
             yuedg(i,j,1) = yulat(i,j)-0.5*yfac ; yuedg(i,j,2) = yulat(i,j)-0.5*yfac
             yuedg(i,j,3) = yulat(i,j)+0.5*yfac ; yuedg(i,j,4) = yulat(i,j)+0.5*yfac
           enddo
         enddo

!-- Meme lissage pour la region Arctique (colonnes i=83..imax, specifique CL15/CL30)
#if ( HRCLIO == 0 )
         do i=83,imax
           if (jsep2(i).ne.jmax) then
#elif ( HRCLIO == 1 )
         jsep2(imax)=jmax
         do i=imax,imax
           if (jsep2(i).ne.jmax) then
#endif
             yfac = 0.25 / float(jmax - jsep2(i))
             jpm  = jsep2(i) - 1
             do j=jsep2(i),jmax
               xslon(i,j)   = xslon(i,jpm)
               yslat(i,j)   = yslat(i,jpm) + yfac*(j-jpm)
               xsedg(i,j,1) = xsedg(i,jpm,1) ; xsedg(i,j,2) = xsedg(i,jpm,2)
               xsedg(i,j,3) = xsedg(i,jpm,3) ; xsedg(i,j,4) = xsedg(i,jpm,4)
               ysedg(i,j,1) = yslat(i,j)-0.5*yfac ; ysedg(i,j,2) = yslat(i,j)-0.5*yfac
               ysedg(i,j,3) = yslat(i,j)+0.5*yfac ; ysedg(i,j,4) = yslat(i,j)+0.5*yfac
               xulon(i,j)   = xulon(i,jpm)
               yulat(i,j)   = yulat(i,jpm) + yfac*(float(j-jpm)-0.5)
               xuedg(i,j,1) = xuedg(i,jpm,1) ; xuedg(i,j,2) = xuedg(i,jpm,2)
               xuedg(i,j,3) = xuedg(i,jpm,3) ; xuedg(i,j,4) = xuedg(i,jpm,4)
               yuedg(i,j,1) = yulat(i,j)-0.5*yfac ; yuedg(i,j,2) = yulat(i,j)-0.5*yfac
               yuedg(i,j,3) = yulat(i,j)+0.5*yfac ; yuedg(i,j,4) = yulat(i,j)+0.5*yfac
             enddo
           endif
         enddo

!-- Angle de rotation pour la grille AA
         do i=1,imax-1
           do j=jsep(i),jmax
             angle(i,j)=atan2((yulat(i+1,j)-yulat(i,j)),                &
                              (xulon(i+1,j)-xulon(i,j)))
           enddo
         enddo
         i=imax
         do j=jsep(i),jmax
           angle(i,j)=angle(i-1,j)
         enddo

       endif  ! fin section grille AA

!-- Construction des tableaux de coins etendus (imax+1,jmax+1) pour Ferret
       do j=1,jmax
         do i=1,imax
           xslonp(i,j) = xsedg(i,j,1) ; yslatp(i,j) = ysedg(i,j,1)
           xulonp(i,j) = xuedg(i,j,1) ; yulatp(i,j) = yuedg(i,j,1)
         enddo
         xslonp(imax+1,j) = xsedg(imax,j,2)
         yslatp(imax+1,j) = ysedg(imax,j,2)
         xulonp(imax+1,j) = xuedg(imax,j,2)
         yulatp(imax+1,j) = yuedg(imax,j,2)
       enddo
       do i=1,imax
         xslonp(i,jmax+1) = xsedg(i,jmax,4)
         yslatp(i,jmax+1) = ysedg(i,jmax,4)
         xulonp(i,jmax+1) = xuedg(i,jmax,4)
         yulatp(i,jmax+1) = yuedg(i,jmax,4)
       enddo
       xslonp(imax+1,jmax+1) = xsedg(imax,jmax,3)
       yslatp(imax+1,jmax+1) = ysedg(imax,jmax,3)
       xulonp(imax+1,jmax+1) = xuedg(imax,jmax,3)
       yulatp(imax+1,jmax+1) = yuedg(imax,jmax,3)

      end subroutine compute_grid_coords

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: [compute_metric_coeffs]                                                        *** DEPEND DE LA GRILLE ***
!
!>     @brief [G2] Coefficients metriques de la grille (4 positions par point).
!
!      DESCRIPTION:
!>     A partir des metriques locales h1/h2 et de leurs derivees (sorties de G1), remplit les tableaux cmx, cmy, smx, smy,
!!     cmxy, smxy (4 positions chacun), les derivees cmxdy, cmydx, le facteur de Coriolis fs2cor et le sinus de latitude
!!     geographique covrai. Applique les corrections de jonction WW/AA et du detroit de Bering.
!!
!!     CONVENTION DES 4 POSITIONS (3eme indice) :
!!       0 = centre de la maille (pt S)   1 = face ouest (i-1/2)
!!       2 = face sud (j-1/2)             3 = coin SW (i-1/2, j-1/2) = pt U
!!
!!     ENTREES :
!!       ltest, jsepar, jeq, jsep(imax)               : geometrie de la grille
!!       iberp, ibera, iberpm, iberam                 : indices i du detroit de Bering (WW, AA) et decales
!!       jberp, jbera, jberpm, jberam                 : indices j du detroit de Bering (WW, AA) et decales
!!       omega, unsrt, radian                         : constantes physiques
!!       dlat, dlong, dyai, dxaj, ylat1, xlon1, yai1, xaj1 : resolutions et origines WW/AA
!!       imax, jmax                                   : dimensions de la grille
!!       h1, h2, d1d2, d2d1, dx1, dx2                  : metriques locales (sorties de G1)
!!
!!     SORTIES :
!!       cmx/cmy(i,j,0:3)   : h1/rterre, h2/rterre aux 4 positions
!!       smx/smy(i,j,0:3)   : 1/cmx, 1/cmy
!!       cmxy/smxy(i,j,0:3) : cmx*cmy (element de surface normalise) et son inverse
!!       cmxdy/cmydx(i,j)   : d(cmx)/dy, d(cmy)/dx au coin SW
!!       fs2cor(i,j)        : facteur de Coriolis = 2*Omega*sin(lat_U)
!!       covrai(i,j)        : sin(lat_S geographique)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine compute_metric_coeffs(                                 &
        ltest, jsepar, jeq, jsep,                                       &
        iberp, ibera, iberpm, iberam,                                   &
        jberp, jbera, jberpm, jberam,                                   &
        omega, unsrt, radian,                                           &
        dlat, dlong, dyai, dxaj,                                        &
        ylat1, xlon1, yai1, xaj1,                                       &
        imax, jmax,                                                     &
        h1, h2, d1d2, d2d1, dx1, dx2,                                   &
        cmx, cmy, smx, smy, cmxy, smxy,                                 &
        cmxdy, cmydx, fs2cor, covrai)

       implicit none

!-- Declaration des arguments
       integer(ip), intent(in) :: ltest, jsepar, jeq, imax, jmax
       integer(ip), intent(in) :: jsep(imax)
       integer(ip), intent(in) :: iberp, ibera, iberpm, iberam
       integer(ip), intent(in) :: jberp, jbera, jberpm, jberam
       real(dblp),  intent(in) :: omega, unsrt, radian
       real(dblp),  intent(in) :: dlat, dlong, dyai, dxaj
       real(dblp),  intent(in) :: ylat1, xlon1, yai1, xaj1
       real(dblp),  intent(in) :: h1(imax,jmax), h2(imax,jmax)
       real(dblp),  intent(in) :: d1d2(imax,jmax), d2d1(imax,jmax)
       real(dblp),  intent(in) :: dx1(imax,jmax), dx2(imax,jmax)

       real(dblp), intent(out) :: cmx(imax,jmax,0:3), cmy(imax,jmax,0:3)
       real(dblp), intent(out) :: smx(imax,jmax,0:3), smy(imax,jmax,0:3)
       real(dblp), intent(out) :: cmxy(imax,jmax,0:3), smxy(imax,jmax,0:3)
       real(dblp), intent(out) :: cmxdy(imax,jmax), cmydx(imax,jmax)
       real(dblp), intent(out) :: fs2cor(imax,jmax), covrai(imax,jmax)

!-- Variables locales
       integer(ip) :: i, j
       real(dblp)  :: ysdeg, yudeg, yurad, yu0rad
       real(dblp)  :: ccsysw, sscysw, ccsyuw, sscyuw, ssnyuw, ssnysw
       real(dblp)  :: ccsysa, sscysa, ccsyua, sscyua
       real(dblp)  :: ccmxdy, ccmydx, ccfcor
       real(dblp)  :: xu0m90, xudeg_aa, yyalim, ysrad
       real(dblp)  :: csxua(jmax), csxsa(jmax)

!==============================================================================
!-- 3.1 : Coefficients metriques pour la grille WW
!==============================================================================

       do j=1,jsepar
         ysdeg = ylat1 + dlat * DFLOAT(j-1)
         yurad = (ysdeg - 0.5 * dlat) * radian

         if (ltest.lt.0) then
!-- Grille cartesienne (cas test)
           ccsysw=1.0 ; sscysw=1.0 ; ccsyuw=1.0 ; sscyuw=1.0
           ssnyuw=sin(45.0*radian) ; ssnysw=sin(45.0*radian)
           ccmxdy=0.0
           ccfcor=omega*sin(45.0*radian)
         else
           ccsysw = cos(ysdeg*radian)  ; sscysw = 1.0/ccsysw
           ccsyuw = cos(yurad)         ; sscyuw = 1.0/ccsyuw
           ssnyuw = sin(yurad)         ; ssnysw = sin(ysdeg*radian)
           ccmxdy = -ssnyuw * unsrt
           ccfcor = omega * ssnyuw
         endif

         do i=1,imax
!-- cmx = h1/R = cos(lat) en grille WW
           cmx(i,j,0)=ccsysw ; cmx(i,j,1)=ccsysw
           cmx(i,j,2)=ccsyuw ; cmx(i,j,3)=ccsyuw
           smx(i,j,0)=sscysw ; smx(i,j,1)=sscysw
           smx(i,j,2)=sscyuw ; smx(i,j,3)=sscyuw
!-- cmy = h2/R = 1 en grille WW (h2=R uniforme en latitude)
           cmy(i,j,0)=1. ; cmy(i,j,1)=1. ; cmy(i,j,2)=1. ; cmy(i,j,3)=1.
           smy(i,j,0)=1. ; smy(i,j,1)=1. ; smy(i,j,2)=1. ; smy(i,j,3)=1.
!-- cmxy = cmx*cmy
           cmxy(i,j,0)=ccsysw ; cmxy(i,j,1)=ccsysw
           cmxy(i,j,2)=ccsyuw ; cmxy(i,j,3)=ccsyuw
           smxy(i,j,0)=sscysw ; smxy(i,j,1)=sscysw
           smxy(i,j,2)=sscyuw ; smxy(i,j,3)=sscyuw
!-- Derivees
           cmxdy(i,j) = ccmxdy  ! d(h1)/dy = -sin(lat_U)/R
           cmydx(i,j) = 0.      ! d(h2)/dx = 0 (h2 uniforme en x)
!-- Coriolis et sin(lat geo)
           fs2cor(i,j) = ccfcor
           covrai(i,j) = ssnysw
         enddo
       enddo  ! fin boucle j grille WW

!==============================================================================
!-- 3.3 : Coefficients metriques pour la grille AA
!==============================================================================

       if (ltest.ge.3) then

         xu0m90 = xaj1 - 1.5*dxaj - 90.
         do j=jeq,jmax
           xudeg_aa = xu0m90 + dxaj * DFLOAT(j)
           csxua(j) = sin(xudeg_aa * radian)
           csxsa(j) = sin((xudeg_aa + 0.5*dxaj) * radian)
         enddo
         yyalim = 90.0 - abs(dyai)

         do i=1,imax
           ysdeg = yai1 + dyai * DFLOAT(i-1)
           if (abs(ysdeg).ge.yyalim) then
             ccsysa=1. ; yurad=0. ; ysrad=0. ; ccsyua=1.
           else
             ccsysa = cos(ysdeg*radian)
             yurad  = (ysdeg - 0.5*dyai)*radian
             ysrad  = ysdeg*radian
             ccsyua = cos(yurad)
           endif
           sscysa = 1.0/ccsysa ; sscyua = 1.0/ccsyua
           ccmydx = sin(yurad)*unsrt

           do j=jsep(i),jmax
!-- cmx = h1/R = 1 en grille AA (direction x = ex-longitude)
             cmx(i,j,0)=1. ; cmx(i,j,1)=1. ; cmx(i,j,2)=1. ; cmx(i,j,3)=1.
             smx(i,j,0)=1. ; smx(i,j,1)=1. ; smx(i,j,2)=1. ; smx(i,j,3)=1.
!-- cmy = h2/R = cos(lat_AA) en grille AA (direction y = ex-latitude)
             cmy(i,j,0)=ccsysa ; cmy(i,j,1)=ccsyua
             cmy(i,j,2)=ccsysa ; cmy(i,j,3)=ccsyua
             smy(i,j,0)=sscysa ; smy(i,j,1)=sscyua
             smy(i,j,2)=sscysa ; smy(i,j,3)=sscyua
!-- cmxy (note : positions 1 et 2 sont des ratios inverses, termes croises)
             cmxy(i,j,0)=ccsysa ; cmxy(i,j,1)=sscyua
             cmxy(i,j,2)=sscysa ; cmxy(i,j,3)=ccsyua
             smxy(i,j,0)=sscysa ; smxy(i,j,1)=ccsyua
             smxy(i,j,2)=ccsysa ; smxy(i,j,3)=sscyua
!-- Derivees
             cmxdy(i,j) = 0.       ! d(h1)/dy = 0 en grille AA
             cmydx(i,j) = ccmydx   ! d(h2)/dx = sin(lat_AA_U)/R
!-- Coriolis en grille AA : f = 2*Omega*cos(lat_AA)*sin(lon_AA)
             fs2cor(i,j) = omega * ccsyua * csxua(j)
             covrai(i,j) = ccsysa * csxsa(j)
           enddo

!-- Correction a la jonction WW/AA (j=jeq) : moyenne des metriques des 2 grilles
           if (jsep(i).eq.jeq) then
             cmy(i,jeq,2)  = 0.5*(1. + cmy(i,jeq,0))
             cmy(i,jeq,3)  = 0.5*(1. + cmy(i,jeq,1))
             smy(i,jeq,2)  = 1./cmy(i,jeq,2)
             smy(i,jeq,3)  = 1./cmy(i,jeq,3)
             cmxy(i,jeq,2) = smy(i,jeq,2)
             cmxy(i,jeq,3) = cmy(i,jeq,3)
             smxy(i,jeq,2) = cmy(i,jeq,2)
             smxy(i,jeq,3) = smy(i,jeq,3)
             cmydx(i,jeq)  = 0.5 * ccmydx
           endif
         enddo  ! fin boucle i grille AA

!-- Corrections metriques au detroit de Bering (jonction WW<->AA)
         if (iberp.ne.0) then
           cmx(ibera, jbera,2) = cmx(iberpm,jberp,2)
           smx(ibera, jbera,2) = smx(iberpm,jberp,2)
           cmy(ibera, jbera,2) = cmy(iberpm,jberp,2)
           smy(ibera, jbera,2) = smy(iberpm,jberp,2)
           cmxy(ibera, jbera,2) = cmxy(iberpm,jberp,2)
           smxy(ibera, jbera,2) = smxy(iberpm,jberp,2)
           cmx(iberam,jbera,2) = cmx(iberp, jberp,2)
           smx(iberam,jbera,2) = smx(iberp, jberp,2)
           cmy(iberam,jbera,2) = cmy(iberp, jberp,2)
           smy(iberam,jbera,2) = smy(iberp, jberp,2)
           cmxy(iberam,jbera,2) = cmxy(iberp, jberp,2)
           smxy(iberam,jbera,2) = smxy(iberp, jberp,2)
!-- Corrections diagnostiques (cmxy au coin SW)
           cmxy(ibera, jbera,3) = cmxy(iberp, jberp,3)
           smxy(ibera, jbera,3) = smxy(iberp, jberp,3)
           cmxy(iberp, jberp,0) = cmxy(iberam,jberam,0)
           cmxy(iberpm,jberp,0) = cmxy(ibera, jberam,0)
           cmxy(ibera, jbera,0) = cmxy(iberpm,jberpm,0)
           cmxy(iberam,jbera,0) = cmxy(iberp, jberpm,0)
           smxy(iberp, jberp,0) = smxy(iberam,jberam,0)
           smxy(iberpm,jberp,0) = smxy(ibera, jberam,0)
           smxy(ibera, jbera,0) = smxy(iberpm,jberpm,0)
           smxy(iberam,jbera,0) = smxy(iberp, jberpm,0)
         endif

       endif  ! fin section grille AA

      end subroutine compute_metric_coeffs

      end module defgrid_compute_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
