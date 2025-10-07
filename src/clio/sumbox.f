!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:55 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:55 CET 2009

      SUBROUTINE sumbox(phix,phiy, u1flow,v2flow, flxbox,
     &                  scabox, surfbx, cmfx,cmfy,
     &                  ibond1,ibond2, jlim1,jlim2, nbox, k, ns, nn99)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! calcule (dans flxbox) le bilan des flux phix & phiy (Entrant et Sortant)
! pour les differentes boites ; ainsi que
! la Somme (-> scabox) du scalaire sur la paroi (et surface eq. -> surfbx)
! N.B: Sans initialisation de flxbox,scabox ni surfbx (=> doit etre fait avant)
!  modif : 13/06/96


!! START_OF_USE_SECTION

      use const_mod,            only: one
      use para0_mod,            only: imax, jmax, nsmax
      use para_mod,             only:
      use bloc0_mod,            only: tms, scal
      use bloc_mod,             only:
      use global_constants_mod, only: dblp=>dp, ip

!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"

!! END_OF_INCLUDE_SECTION

!--dummy variables :
      real(dblp),  dimension(imax,jmax) :: phix, phiy, u1flow, v2flow
      real(dblp),  dimension(0:14,*)    :: flxbox
      real(dblp),  dimension(*)         :: scabox, surfbx
      integer(ip), dimension(jmax,*)    :: ibond1, ibond2
      integer(ip), dimension(*)         :: jlim1, jlim2

!- variables locales :
      logical sflag


      real(dblp) :: cmfx, cmfy
      integer(ip):: nbox, k, ns, nn99, icr, icr1, ii, ii1, ii2, isigne
     &            , jcr, jj, m, nb
      real(dblp) :: ttm1,ttm2,xx0,xxs
      real(dblp) :: alpha,c2alph


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initialisation.
!-----------------------------------------------------------------------

!- taux de decentrement :
      alpha = 0.0

      sflag = ns.ge.1 .and. ns.le.nsmax

      do 600 nb=1,nbox
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Bilan sur chacun des murs de chacune des boites "nb" :
!-----------------------------------------------------------------------

!--Somme le flux phiy sur le mur Sud (m=3,4,5) :
      jcr = jlim1(nb)
      do 210 ii=ibond1(jcr,nb),ibond2(jcr,nb)
        xx0 = 0.5 * cmfy * phiy(ii,jcr)
        xxs = xx0 * sign(one,v2flow(ii,jcr))
        flxbox( 4,nb) = flxbox(4,nb) + xxs + xx0
        flxbox( 5,nb) = flxbox(5,nb) + xxs - xx0
 210  continue
      if (sflag) then
!- moyenne du scalaire correspondant sur le mur :
        do 215 ii=ibond1(jcr,nb),ibond2(jcr,nb)
          ttm2 = tms(ii,jcr-1,k)*tms(ii,jcr,k)
          surfbx(nb) = surfbx(nb) + ttm2
          c2alph = sign(alpha,v2flow(ii,jcr))
          scabox(nb) = scabox(nb) + ttm2 * 0.5 * (
     &           (1.+c2alph)*scal(ii,jcr-1,k,ns)
     &         + (1.-c2alph)*scal(ii, jcr, k,ns)  )
 215    continue
      endif

!--Somme le flux phiy sur le mur Nord (m=6,7,8) :
      jcr = jlim2(nb)
      do 230 ii=ibond1(jcr,nb),ibond2(jcr,nb)
        xx0 = 0.5 * cmfy * phiy(ii,jcr+1)
        xxs = xx0 * sign(one,v2flow(ii,jcr+1))
        flxbox( 7,nb) = flxbox(7,nb) + xxs - xx0
        flxbox( 8,nb) = flxbox(8,nb) + xxs + xx0
 230  continue
      if (sflag) then
!- moyenne du scalaire correspondant sur le mur :
        do 235 ii=ibond1(jcr,nb),ibond2(jcr,nb)
          ttm2 = tms(ii,jcr,k)*tms(ii,jcr+1,k)
          surfbx(nb) = surfbx(nb) + ttm2
          c2alph = sign(alpha,v2flow(ii,jcr+1))
          scabox(nb) = scabox(nb) + ttm2 * 0.5 * (
     &           (1.+c2alph)*scal(ii, jcr, k,ns)
     &         + (1.-c2alph)*scal(ii,jcr+1,k,ns)  )
 235    continue
      endif

      do 290 jj=jlim1(nb),jlim2(nb)
!--Somme le flux phix sur le mur West (m=9,10,11) :
        icr = ibond1(jj,nb)
        xx0 = 0.5 * cmfx * phix(icr,jj)
        xxs = xx0 * sign(one,u1flow(icr,jj))
        flxbox(10,nb) = flxbox(10,nb) + xxs + xx0
        flxbox(11,nb) = flxbox(11,nb) + xxs - xx0
        if (sflag) then
!- moyenne du scalaire correspondant sur le mur :
          ttm1 = tms(icr-1,jj,k)*tms(icr,jj,k)
          surfbx(nb) = surfbx(nb) + ttm1
          c2alph = sign(alpha,u1flow(icr,jj))
          scabox(nb) = scabox(nb) + ttm1 * 0.5 * (
     &           (1.+c2alph)*scal(icr-1,jj,k,ns)
     &         + (1.-c2alph)*scal( icr, jj,k,ns)  )
        endif

!--Contribution du flux phiy sur le mur West (m=9,10,11) :
        icr1 = ibond1(jj+1,nb)
        isigne = sign(1,icr-icr1)
        ii1 = min(icr,icr1)
        ii2 = max(icr,icr1) - 1
        do 250 ii=ii1,ii2
          xx0 = 0.5 * cmfy * phiy(ii,jj+1)
          xxs = xx0 * sign(one,v2flow(ii,jj+1)) * isigne
          flxbox(10,nb) = flxbox(10,nb) + xxs + xx0
          flxbox(11,nb) = flxbox(11,nb) + xxs - xx0
 250    continue
        if (sflag) then
          do 255 ii=ii1,ii2
!- moyenne du scalaire correspondant sur le mur :
            ttm2 = tms(ii,jj,k)*tms(ii,jj+1,k)
            surfbx(nb) = surfbx(nb) + ttm2
            c2alph = sign(alpha,v2flow(ii,jj+1))
            scabox(nb) = scabox(nb) + ttm2 * 0.5 * (
     &              (1.+c2alph)*scal(ii, jj, k,ns)
     &            + (1.-c2alph)*scal(ii,jj+1,k,ns)  )
 255      continue
        endif

!--Somme le flux phix sur le mur Est (m=12,13,14) :
        icr = ibond2(jj,nb)
        xx0 = 0.5 * cmfx * phix(icr+1,jj)
        xxs = xx0 * sign(one,u1flow(icr+1,jj))
        flxbox(13,nb) = flxbox(13,nb) + xxs - xx0
        flxbox(14,nb) = flxbox(14,nb) + xxs + xx0
        if (sflag) then
!- moyenne du scalaire correspondant sur le mur :
          ttm1 = tms(icr,jj,k)*tms(icr+1,jj,k)
          surfbx(nb) = surfbx(nb) + ttm1
          c2alph = sign(alpha,u1flow(icr+1,jj))
          scabox(nb) = scabox(nb) + ttm1 * 0.5 * (
     &           (1.+c2alph)*scal( icr, jj,k,ns)
     &         + (1.-c2alph)*scal(icr+1,jj,k,ns)  )
        endif

!--Contribution du flux phiy sur le mur Est (m=12,13,14) :
        icr1 = ibond2(jj+1,nb)
        isigne = sign(1,icr-icr1)
        ii1 = min(icr,icr1) + 1
        ii2 = max(icr,icr1)
        do 270 ii=ii1,ii2
          xx0 = 0.5 * cmfy * phiy(ii,jj+1)
          xxs = xx0 * sign(one,v2flow(ii,jj+1)) * isigne
          flxbox(13,nb) = flxbox(13,nb) + xxs - xx0
          flxbox(14,nb) = flxbox(14,nb) + xxs + xx0
 270    continue
        if (sflag) then
          do 275 ii=ii1,ii2
!- moyenne du scalaire correspondant sur le mur :
            ttm2 = tms(ii,jj,k)*tms(ii,jj+1,k)
            surfbx(nb) = surfbx(nb) + ttm2
            c2alph = sign(alpha,v2flow(ii,jj+1))
            scabox(nb) = scabox(nb) + ttm2 * 0.5 * (
     &              (1.+c2alph)*scal(ii, jj, k,ns)
     &            + (1.-c2alph)*scal(ii,jj+1,k,ns)  )
 275      continue
        endif
!----
 290  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Somme des bilans partiel .
!-----------------------------------------------------------------------

      do 310 m=3,14,3
!- Total(0:3) = Entrant(1:3) - Sortant(2:3) :
        flxbox(m,nb) = flxbox(m+1,nb) - flxbox(m+2,nb)
 310  continue

      do 330 m=0,2
!- Total(<3) = Somme des Flx sur les 4 murs :
        flxbox(m,nb) = flxbox(m+3,nb) + flxbox(m+6,nb)
     &               + flxbox(m+9,nb) + flxbox(m+12,nb)
 330  continue

!--fin de la boucle sur la boite "nb".
 600  continue

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine sumbox -
      end
