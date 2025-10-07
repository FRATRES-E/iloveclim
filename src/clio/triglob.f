!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:56 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:56 CET 2009

      SUBROUTINE triglob(kvar,kspv,im,jm)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
! Modification des valeurs kspv (-> ksp+1) du masque (entier) "kmasq" :
!  supprime les "kspv" "hors domaine" qui ne correspondent pas a du continent
!
!  modif : 09/08/94

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use reper_mod
      
      use newunit_clio_mod, only: clio3_out_id, mouchard_id
      
!! END_OF_USE_SECTION


!--dummy variables :
      real(kind=dblp), dimension(imax,jmax) ::  kvar


!--- loads of locales
      real(kind=dblp) :: ccsya, degre, dxamin, dxwmin, dyamin, dywmin
     &                 , ssnyw, unmin, xa, xx, xx0, xx0w, xxard, ya, yy
     &                 , yyard, ccsyw, ssnya, xw, xxwrd, yw, yywrd

      integer(kind=ip):: i, iea, im, isa, j, jm, kspv, nx, ny, ia, iw
     &                 , ja, jw

     
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|

      degre = 180. / pi
      unmin =  1. - epsil

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  1 ) Determination des Indices Debut/Fin de la grille AA :           |
!-----------------------------------------------------------------------

      isa = 0
      iea = 0
      do i=imax,1,-1
       do j=jsep(i),jmax
        if (kvar(i,j).ne.kspv) then
          isa = i
          iea = max(i,iea)
        endif
       enddo
      enddo
      isa = max(isa-1, 1)
      iea = min(iea+1, imax)

!- pour la cyclicite de la grille WW :
      xx0w = xwi1 + 0.5 * dxwi

!- inferieur a la maille :
      dxamin = unmin * dxaj
      dyamin = unmin * dyai
      dxwmin = unmin * dxwi
      dywmin = unmin * dywj

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  2 ) Traitement separe des 3 domaines : WW, AA, ni AA ni WW :
!-----------------------------------------------------------------------

      jloop: do j=1,jm
       iloop: do i=1,im
       
        if (kvar(i,j).ne.kspv) cycle iloop

        if (j.ge.jsep(i) .and. i.ge.isa .and. i.le.iea) then
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  3 ) Traitement d'un point de la grille AA .                         |
!-----------------------------------------------------------------------

      xa = xaj1 + dxaj * DFLOAT(j-1)
      ya = yai1 + dyai * DFLOAT(i-1)

      xx0 = xa - 0.5 * dxamin
      xx  = xx0
      yy  = ya - 0.5 * dyamin
      do ny=1,2
       do nx=1,2
!--Examine chaque coin de la maille (i,j) , grille AA :

!--Calcul de la position sur la grille WW :
        xxard = xx * radian
        yyard = yy * radian
        ccsya = cos(yyard)
        ssnyw = -ccsya * cos(xxard)
        if ( abs(ssnyw).ge.unmin ) then
          write(clio3_out_id,*) 'ARRET dans "triglob" , (i,j) = ', i, j
          write(clio3_out_id,*)
     &      'Pt Continental de AA trop pre du Pole de WW !'
          stop
        endif
        yw = asin(ssnyw) * degre
        xw = atan2(ccsya * sin(xxard), sin(yyard)) * degre
!- changement d'intervalle pour xw  ]xx0w,xx0w+360]  :
        xw = xx0w + mod(xwpoln + xw - xx0w + untour, untour)

!--Recherche de la maille de WW correspondant :
        iw = 1 + nint( (xw-xwi1) / dxwi )
        jw = 1 + nint( (yw-ywj1) / dywj )
!--Si (i,j) coincide avec maille oceanique de WW => Exclusion
        if ( iw.ge.1 .and. iw.le.imax .and.
     &       jw.ge.1 .and. jw.lt.jsep(iw) .and.
     &       kvar(iw,jw).lt.kspv ) then
          kvar(i,j) = kspv+1        !--Exclusion du point (i,j)
          cycle iloop               !  (repeated the line preceding enddo iloop)
        endif
!-
         xx = xx + dxamin
        enddo
        xx = xx0
        yy = yy + dyamin
       enddo

      cycle iloop

!--fin du traitement du pt de AA .
        endif

        if (j.lt.jsep(i)) then
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  5 ) Traitement des donnees de la grille WW .                        |
!-----------------------------------------------------------------------

      xw = xwi1 + dxwi * DFLOAT(i-1)
      yw = ywj1 + dywj * DFLOAT(j-1)

      xx0 = xw - 0.5 * dxwmin
      xx  = xx0
      yy  = yw - 0.5 * dywmin
      do ny=1,2
       do nx=1,2
!--Examine chaque coin de la maille (i,j) , grille WW :

!--Calcul de la position sur la grille AA :
!- changement d'origine des longitudes :
        xxwrd = mod((untour + xx - xwpoln), untour) * radian
        yywrd = yy * radian
        ccsyw = cos(yywrd)
        ssnya = ccsyw * cos(xxwrd)
        if ( abs(ssnya).ge.unmin ) then
          write(clio3_out_id,*) 'ARRET dans "triglob" , (i,j) = ', i, j
          write(clio3_out_id,*)
     &      'Pt Continental de WW trop pre du Pole de AA !'
          stop
        endif
        ya = asin(ssnya) * degre
        xa = atan2(ccsyw * sin(xxwrd), -sin(yywrd)) * degre
!- changement d intervalle pour xa  ]+0,+360]  :
        xa = mod(xa + untour, untour)

!--Recherche de la maille de AA correspondant :
        ia = 1 + nint( (ya-yai1) / dyai )
        ja = 1 + nint( (xa-xaj1) / dxaj )
!------
!     if (i.eq.icheck .and. j.eq.jcheck) then
!       write(mouchard_id,*) 'nx, ny, xx, yy :'
!       write(mouchard_id,*)  nx, ny, xx, yy
!       write(mouchard_id,*) 'ia, ja, xa, ya , kvar(ia,ja) :'
!       write(mouchard_id,*)  ia, ja, xa, ya , kvar(ia,ja)
!     endif
!------
!--Si (i,j) coincide avec maille oceanique de AA => Exclusion
        if ( ia.ge.isa .and. ia.le.iea .and.
     &       ja.ge.jsep(ia) .and. ja.le.jmax .and.
     &       kvar(ia,ja).lt.kspv ) then
          kvar(i,j) = kspv+1        !--Exclusion du point (i,j)
          cycle iloop               !  (repeated the line preceding enddo iloop)
        endif
!--
        xx = xx + dxwmin
        enddo
       xx = xx0
       yy = yy + dywmin
       enddo
       
      cycle iloop

!--fin du traitement du pt de WW .
        endif

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|

!--Exclusion du point (i,j) :
        kvar(i,j) = kspv+1

       enddo iloop
      enddo jloop

!       write(mouchard_id,*) 'icheck, jcheck, kvar(icheck,jcheck) :'
!       write(mouchard_id,*)  icheck, jcheck, kvar(icheck,jcheck)

      return
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!- fin de la routine triglob -
      end
