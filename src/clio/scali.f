!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009

      SUBROUTINE scali
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!**AVANCEMENT IMPLICITE D'UN PAS DE TEMPS SCALAIRE DES SCALAIRES**
!  Traitement implicite des termes de rappel, diffusion et advection .
!    advect & diffus : taux d'implicitete au choix.
!  advection suivant la verticale :  w(k) * [ S(k-1) + S(k) ] / 2
!- resolution separee ou non, selon "ns"
!  modif : 23/03/98

!! START_OF_USE_SECTION

      use const_mod,            only:

      use para0_mod,            only: imax, jmax, kmax, nsmax
      use para_mod,             only:
      use bloc0_mod,            only: js1, js2, ks1, ks2, txiads, tms
     &                              , is1, is2, unsdz, rappel, scalr, w
     &                              , alphaz, dts, scal, unsdzw, ai, avsdz
      use bloc_mod,             only: nstart, numit
      use isoslope_mod,         only: wiso, c4y, c4x

      use global_constants_mod, only: dblp=>dp, ip

      use newunit_clio_mod, only: clio3_out_id

!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "isoslope.com"

!! END_OF_INCLUDE_SECTION

!--variables a conserver d'un appel a l'autre :
      logical, dimension(nsmax), save :: bcalc
c~       common / bscali /
c~      &  bcalc(nsmax)

      integer(ip), dimension(nsmax), save :: nscom
c~       common / nscali /
c~      &  nscom(nsmax)

      real(dblp), dimension(imax,jmax,2:kmax,nsmax), save :: z2alph

c~       common / xscali /
c~      &  z2alph(imax,jmax,2:kmax,nsmax)

!--variables locales :
      real(dblp), dimension(imax,kmax,nsmax) :: aa, bb, cc
      real(dblp), dimension(imax,kmax)       :: ff
!     dimension ccwabs(kmax), ccwcor(kmax), ccwdmy(kmax)
!     dimension ccwddn(kmax), ccwdup(kmax)
      real(dblp), dimension(kmax)            :: cczdt
!     dimension ccvois(-kmax:kmax), alpha(imax,kmax)

      real(dblp)  :: cc0w, cc1w, cc2w, ccadv, ccdif, ccdt, ccstep, ccwi
      integer(ip) :: i, j, k, nn, nnsep, ns



!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  0 ) preparation des coeff. intervenant dans les flux d'advection .  |
!-----------------------------------------------------------------------

      if(numit.eq.nstart) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--Resolution commune / separee :
      nnsep = nsmax
      do ns=1,nsmax
        bcalc(ns) = .TRUE.
        nscom(ns) = ns
        do nn=1,ns-1
         if (bcalc(ns).and.(alphaz(1-nn+ks2).eq.alphaz(1-ns+ks2))) then
           bcalc(ns) = .FALSE.
           nscom(ns) = nn
           nnsep = nnsep - 1
         endif
        enddo
      enddo

!     write(clio3_out_id,'(A,) ' scali : alphaz = run.par + Dec. near a step'
      write(clio3_out_id,
     & '(A,I2,A,(6F5.2))') ' scali : cal=', nnsep,
     &    ' ; alphaz = Dec.near_step + F(ns):',
     &    (alphaz(1-ns+ks2),ns=1,nsmax)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- initialisation :

!dmr [###] Simpler versions ...
!      do ns=1,nsmax
!       do k=2,kmax
!        do j=1,jmax
!         do i=1,imax
!          z2alph(i,j,k,ns) = 0.0
!         enddo
!        enddo
!       enddo
!      enddo
!dmr [###] Simpler versions ...

      z2alph(1:imax,1:jmax,2:kmax,1:nsmax) = 0.0d0

!--Mise en place des coeff (Decentrement Vertical) :
      cc0w = 0.5 * alphaz(ks1)
      do ns=1,nsmax
       if (bcalc(ns)) then
         cc1w = 0.5 * txiads * alphaz(1-ns+ks2)
         cc2w = 0.5 * (1.0 - alphaz(1-ns+ks2))
         do k=ks1+1,ks2
!         cc1w = 0.5 * txiads * alphaz(k)
!         cc2w = 0.5 * (1.0 - alphaz(k))
          do j=js1,js2
           do i=is1(j),is2(j)
             ccstep = ( tms(i-1,j,k) - tms(i-1,j,k-1) )
     &              + ( tms(i+1,j,k) - tms(i+1,j,k-1) )
     &              + ( tms(i,j-1,k) - tms(i,j-1,k-1) )
     &              + ( tms(i,j+1,k) - tms(i,j+1,k-1) )
             z2alph(i,j,k,ns) = cc1w + min(cc2w , cc0w * ccstep)
           enddo
          enddo
	 enddo
       endif
      enddo
      
!--Fin du traitement specifique de la 1ere itt.
      endif

! convention :
!    CCW* coeff relatif a l'interface (= place des W).
!    CCZ* coeff relatif au centre des boites (= place de T,S)
!    CCWI : tient compte du taux Implicite.
!    CCWABS(k) intervient avec |w| ;
!              CCWUP(k) avec scal(k) ; CCWDowN(k) avec scal(k-1)
!  pour le taux de Decentrement calcule a partir du Nb de courant : CCWDUP(k)
!    avec dts(k), CCWDDowN(k) avec dts(k-1) ey CCWDMY(k) avec la moyenne des 2
!    CCZDT**(k) intervient pour le bilan des flux V, boite "k" ;

      ccwi   = 0.5 * txiads
      do k=ks1,ks2
        cczdt(k) = dts(k) * unsdz(k)
!       ccwabs(k) = ccwi * alphaz(k)
      enddo
!     ccwid2 = 0.25 * txiads
!     do k=ks1+1,ks2
!       ccwcor(k) = ccwid2 * unsdzw(k) * (dz(k) - dz(k-1))
!    &                     * (1. - alphaz(k))
!       ccwdup(k) = ccwi * unsdzw(k) * dts(k)
!       ccwddn(k) = ccwi * unsdzw(k) * dts(k-1)
!       ccwdmy(k) = ccwi * unsdzw(k) * (dts(k) + dts(k-1)) * 0.5
!       cczdif(k) = txidfs * unsdzw(k)
![###]    enddo

!--Debut de la boucle externe sur l'indice de latitude j :
c~ !$DIR SPP LOOP_PARALLEL
c~ !$DIR SPP LOOP_PRIVATE(i,k,ns,nn,ccdif,ccadv,ccdt,aa,bb,cc,ff)

c~ !$OMP PARALLEL
c~ !$OMP DO PRIVATE(i,k,ns,nn,ccdif,ccadv,ccdt,aa,bb,cc,ff)
      jloop: do j=js1,js2

!-----
!--debut de la boucle sur "ns".
      nsloop: do ns=1,nsmax

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) construction de la matrice (commune a certains scalaires).      |
!-----------------------------------------------------------------------

!--Ininialisation de la diagonale (aa) et Mise en place du 2nd membre (ff) :
!- transfert de scal(*,j,*,ns) dans ff(*,*) & Traite le terme rappel implic.
      do k=ks1,ks2
       do i=is1(j),is2(j)
        aa(i,k,ns) = 1.0 + rappel(i,j,k)
        ff(i,k) = scal(i,j,k,ns) + rappel(i,j,k) * scalr(i,j,k,ns)
       enddo
      enddo
      
      if (bcalc(ns)) then
!--Traitement commun a certains scalaires :

!-N.B.-: avsdz = txidfs * Dif.Vert. / dzw !!
      do k=ks1+1,ks2
       do i=is1(j),is2(j)
        ccdif = avsdz(i,j,k)
     &        + ai(k)*( c4x(i,j,k)*c4x(i,j,k)
     &                 +c4y(i,j,k)*c4y(i,j,k) )*unsdzw(k)
     &        + z2alph(i,j,k,ns) * abs(w(i,j,k)+wiso(i,j,k))
!is0 &        + z2alph(i,j,k,ns) * abs(w(i,j,k))
!    &        + ccwabs(k) * abs(w(i,j,k))
!    &        + ccwcor(k) * w(i,j,k)
!    &        + ccwdmy(k) * w(i,j,k) * w(i,j,k)
!    &        + ccwddn(k) * w(i,j,k) * w(i,j,k)
!    &        + ccwdup(k) * w(i,j,k) * w(i,j,k)
        ccadv = ccwi * w(i,j,k)
     &         +ccwi * wiso(i,j,k)

!- effet du flux phiz(k) sur S(k) : aa(k)*S(k) + bb(k)*S(k-1)
        ccdt = cczdt(k) * tms(i,j,k-1)
        bb(i,k,ns) = ccdt * (-ccdif - ccadv)
        aa(i,k,ns) = aa(i,k,ns) + ccdt * (ccdif - ccadv)
!- effet du flux phiz(k) sur S(k-1) : aa(k-1)*S(k-1) + cc(k-1)*S(k)
        ccdt = cczdt(k-1) * tms(i,j,k-1)
        aa(i,k-1,ns) = aa(i,k-1,ns) + ccdt * (ccdif + ccadv)
        cc(i,k-1,ns) = ccdt * (-ccdif + ccadv)
       enddo
      enddo
 !**fin construction de la matrice du systeme commun**


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) decomposition LU - (methode et notation : cf Linear Algebra pp165-167)
!-----------------------------------------------------------------------

!--calcul de 1/alpha(k) dans aa(i,k) et beta(k) dans bb(i,k)
      do i=is1(j),is2(j)
        aa(i,ks1,ns) = 1.0 / aa(i,ks1,ns)
      enddo
      do k=ks1+1,ks2
       do i=is1(j),is2(j)
        bb(i,k,ns) = bb(i,k,ns) * aa(i,k-1,ns)
        aa(i,k,ns) = 1.0 / ( aa(i,k,ns) - bb(i,k,ns) * cc(i,k-1,ns) )
       enddo
      enddo
!**fin decomposition LU commune**

!--Fin de la partie commune a certains scalaires / debut du traitement separe
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) substitutions avant et arriere pour chaque scalaire :           |
!-----------------------------------------------------------------------

      nn = nscom(ns)
!--calcul de g(k) dans ff(i,k) :
      do k=ks1+1,ks2
       do i=is1(j),is2(j)
        ff(i,k) = ff(i,k) - bb(i,k,nn) * ff(i,k-1)
       enddo
      enddo
      
!--calcul de x(k) dans scal(i,j,k,ns) :
      do i=is1(j),is2(j)
       scal(i,j,ks2,ns) = ff(i,ks2) * aa(i,ks2,nn)
      enddo
      do k=ks2-1,ks1,-1
       do i=is1(j),is2(j)
        scal(i,j,k,ns) = (ff(i,k) - cc(i,k,nn)*scal(i,j,k+1,ns))
     &                 * aa(i,k,nn)
       enddo
      enddo
!--fin substitutions avant et arriere pour chaque scalaire.
!-----

      enddo nsloop 

      enddo jloop
c~ !$OMP END DO
c~ !$OMP END PARALLEL

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine scali -
      end
