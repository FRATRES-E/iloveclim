!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:56 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:56 CET 2009

      SUBROUTINE uvbfet(ccsplt)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  Fusion des routines conti2d + uvb (= uvbe + uvbi) :
!  Resolution complete du mode barotrope, Avec Filtre sur "eta" (ahe <> 0)
!- Detroit de Bering : Flux = coeff*(EtaPac-EtaArc), evalue sur 2 x 4 points.
!  modif : 23/09/99


!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
!! END_OF_USE_SECTION


!--variables locales equivalentes :
!dmr [NOEQUI]      dimension tm1x2(imax,jmax), tm2x2(imax,jmax)
!dmr [NOEQUI]      equivalence ( tm1x2(1,1)  , phihhh(1,1,1) )
!dmr [NOEQUI]      equivalence ( tm2x2(1,1)  , phihhh(1,1,2) )


!dmr [NOEQUI]      dimension unsdet(imax,jmax), phiypy(imax,jmax)
!dmr [NOEQUI]      dimension phiypx(imax,jmax), phiymx(imax,jmax)

!dmr [NOEQUI]      equivalence ( unsdet(1,1) , phihhh(1,1,3) )
!dmr [NOEQUI]      equivalence ( phiypy(1,1) , phihhh(1,1,4) )
!dmr [NOEQUI]      equivalence ( phiypx(1,1) , phihhh(1,1,5) )
!dmr [NOEQUI]      equivalence ( phiymx(1,1) , phihhh(1,1,6) )

!--pour raison de place memoire :

      real(kind=dblp), dimension(imax)     :: phixj, phjxpx, phixju
     &                                      , phixjv
      real(kind=dblp), dimension(imax,jmax):: phiyu, phiyv
     &                                      , detadx, detady


      real(kind=dblp) :: afdtb, ccde, ccdif, ccfexp, ccfimp, ccgdx
     &                 , ccgdy, ccsplt, ccxdif, ccydif, cu, cv, etaarc
     &                 , etapac,  pshber, vvsber
      integer(kind=ip):: i, j

c~       dimension phixj(imax), phjxpx(imax)
c~       dimension phixju(imax), phixjv(imax)
c~       dimension phiyu(imax,jmax), phiyv(imax,jmax)
c~       dimension detadx(imax,jmax), detady(imax,jmax)

!--pour la sortie sur fichier "splout" des variations de eta pendant
!   la derniere iteration  -- sauvegarde de "eta" dans "etasub" :
!     include 'split2.com'

!--initialisation de phi_x,y,xpx,ypy,ypx,umx : (effectuee dans "barot")

      ccde = 0.5 * dtb * ahe * unsdx * unsdy
      ccxdif = ahu * unsdx
      ccydif = ahu * unsdy
      ccgdx = gpes * uns2dx
      ccgdy = gpes * uns2dy
      vvsber = vbspl(iberp,jberp)
      ccfexp = 2. - 2.*txifcb
      ccfimp = 2.*txifcb*dtb

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) calcul de l'elevation, cas avec Filtre (ahe non nul) :          |
!-----------------------------------------------------------------------

!- Calcul de phiypx & phiymx : transfere a la fin de "uvbfet" et dans "barot".
      j=js1
      do 210 i=isf1(j),isf2(j)
!dmr [NOEQUI]
        phihhh(i,j,4) = (eta(i,j-1)-eta(i,j)) * phihhh(i,j,2)
 210  continue
      j=1+js2
      do 220 i=isf1(j),isf2(j)
!dmr [NOEQUI]
        phihhh(i,j,4) = (eta(i,j-1)-eta(i,j)) * phihhh(i,j,2)
 220  continue

!--raccord Bering : Flux croises nuls :
      if (ltest.eq.3) then
        phihhh(iberp,jberp,5) = 0.
        phihhh(iberp,jberp,6) = 0.
      endif

!--Debut de la 1ere boucle externe sur l'indice de latitude j :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,ccdif,phixj,phjxpx)
      do 300 j=js1,js2
!-----
!--Calcul des Flux :
      do 230 i=is1(j),is2(j)+1
        phixj(i) = cmy(i,j,1) * (ub(i,j)+ub(i,j+1))
!dmr [NOEQUI]
        phjxpx(i) = (eta(i-1,j)-eta(i,j)) * phihhh(i,j,1)
 230  continue

!--Bilan des Flux :
!dmr [NOEQUI]
      do 250 i=is1(j),is2(j)
        eta(i,j) = eta(i,j)
     &         - dtb * w(i,j,ks2+1)
     &         + dtb * smxy(i,j,0) *
     &         ( uns2dx*(phixj(i)-phixj(i+1))
!    &         + uns2dy*(phiy(i,j)-phiy(i,j+1)) )
     &         + uns2dy*( cmx(i, j, 2)*(vb(i, j )+vb(i+1, j ))
     &                  - cmx(i,j+1,2)*(vb(i,j+1)+vb(i+1,j+1)) ))
     &  + ccde*smxy(i,j,0)*( (phjxpx(i)-phjxpx(i+1))
     &                     + (phihhh(i,j,4)-phihhh(i,j+1,4))
     &                     - (phihhh(i,j,5)-phihhh(i+1,j+1,5))
     &                     - (phihhh(i+1,j,6)-phihhh(i,j+1,6))  )
 250  continue
      do 255 i=is1(j),is2(j)
        eta(i,j) =  tms(i,j,ks2) * eta(i,j)
        etaspl(i,j) = etaspl(i,j) + ccsplt*eta(i,j)
 255  continue

!--calcul des flux "visqueux" dans la direction y. (issu de "uvb")
      do 270 i=iuf1(j),iuf2(j)
        ccdif = cmxy(i,j,1) * ccydif
        phiyu(i,j) = ccdif * (ub(i,j)-ub(i,j+1))
        phiyv(i,j) = ccdif * (vb(i,j)-vb(i,j+1))
 270  continue

!--Fin de la 1ere boucle externe sur l'indice de latitude j .
 300  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) raccords cycliques et autres  pour "eta" .                      |
!-----------------------------------------------------------------------

      call raccord(eta, zero, 1, 8)

!--pour la sortie sur fichier "splout" des variations de eta pendant
!   la derniere iteration  -- calcul de "etasub" et ecriture sur fichier :
!     include 'split3.com'

!--Fin de l'ancienne routine "conti2d" - Debut de l'ancienne routine "uvb".

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) calcul des flux "visqueux" dans les directions x et y .     |
!-----------------------------------------------------------------------

!--Debut de la 2nd  boucle externe sur l'indice de latitude j :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,ccdif,phixju,phixjv,cu,cv,afdtb)
      do 700 j=ju1,ju2
!-----

!--calcul des flux "visqueux" dans la direction x.
      do 450 i=iu1(j)-1,iu2(j)
        ccdif = smxy(i,j,2) * ccxdif
        phixju(i) = ccdif * (ub(i,j)-ub(i+1,j))
        phixjv(i) = ccdif * (vb(i,j)-vb(i+1,j))
 450  continue

!--fin de calcul des flux horizontaux.

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Calcul du Forcage lie a l'elevation et Bilan des flux .         |
!-----------------------------------------------------------------------

      do 530 i=iu1(j),iu2(j)
        detadx(i,j) = ccgdx * smx(i,j,3) *
     &              ((eta(i-1,j-1)-eta(i,j))+(eta(i-1,j)-eta(i,j-1)))
        detady(i,j) = ccgdy * smy(i,j,3) *
     &              ((eta(i-1,j-1)-eta(i,j))-(eta(i-1,j)-eta(i,j-1)))
        cu = fub(i,j,ks2) + ccfexp*fs2cor(i,j)*vb(i,j)
     &     + hu(i,j)*detadx(i,j)
     &     + smxy(i,j,3) * ( unsdx*(phixju(i-1)-phixju(i))
     &                     + unsdy*(phiyu(i,j-1)-phiyu(i,j)) )
        cv = fvb(i,j,ks2) - ccfexp*fs2cor(i,j)*ub(i,j)
     &     + hu(i,j)*detady(i,j)
     &     + smxy(i,j,3) * ( unsdx*(phixjv(i-1)-phixjv(i))
     &                     + unsdy*(phiyv(i,j-1)-phiyv(i,j)) )
        ub(i,j) = ub(i,j) + dtb * cu
        vb(i,j) = vb(i,j) + dtb * cv
 530  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  6 ) Coriolis : partie implicite .                                   |
!-----------------------------------------------------------------------
!dmr [NOEQUI]
      do 630 i=iu1(j),iu2(j)
        afdtb = fs2cor(i,j)*ccfimp
        cu = ub(i,j) + afdtb*vb(i,j)
        cv = vb(i,j) - afdtb*ub(i,j)
        ub(i,j) = cu * phihhh(i,j,3)
        vb(i,j) = cv * phihhh(i,j,3)
        ubspl(i,j) = ubspl(i,j) + ccsplt*ub(i,j)
        vbspl(i,j) = vbspl(i,j) + ccsplt*vb(i,j)
 630  continue

!--Filtre : calcul des flux diagonaux et N-S (issu de "conti2d").
      do 660 i=isf1(j),isf2(j)
!dmr [NOEQUI]
        phihhh(i,j,4) = (eta(i,j-1)-eta(i,j)) * phihhh(i,j,2)
 660  continue
      do 670 i=iu1(j),iu2(j)+1
        phihhh(i,j,5) = (eta(i-1,j-1)-eta(i,j))*tmu(i,j,ks2)*cmxy(i,j,3)
        phihhh(i,j,6) = (eta(i,j-1)-eta(i-1,j))*tmu(i,j,ks2)*cmxy(i,j,3)
 670  continue

!--Fin de la 2nd boucle externe sur l'indice de latitude j .
 700  continue

!--Parametrisation du Flux dans le Detroit de Bering :
!-  EtaPac & Eta.Arc => evalues sur 4 points chacune, separes d'au moins
!-  1 maille de la position du vrai detroit.

      if (ltest.eq.3 .and. iberp.ne.ibera) then
        etapac = eta(iberp-1,jberp-3) + eta(iberp,jberp-3)
     &         + eta(iberp-1,jberp-2) + eta(iberp,jberp-2)
        etaarc = eta(ibera,jbera-3) + eta(ibera+1,jbera-3)
     &         + eta(ibera,jbera-2) + eta(ibera+1,jbera-2)
        pshber = bering * 0.25 * (etapac - etaarc)
        vb(iberp,jberp) = smx(iberp,jberp,2) * unsdx * pshber
        ub(iberp,jberp) = 0.
        vbspl(iberp,jberp) = vvsber + ccsplt*vb(iberp,jberp)
        ubspl(iberp,jberp) = 0.
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  8 ) Raccord pour ub et vb .                                         |
!-----------------------------------------------------------------------

      call raccord( ub, zero, 1, 23)
      call raccord( vb, zero, 1, 35)

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine uvbfet -
      end
