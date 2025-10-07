!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009

      SUBROUTINE barot
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  prepare l integration du mode barotrope.
!  remise a zero des variables eta/ub/vb-spl (accumulateur&moyennes).
!  modif : 26/03/99
!  [UPDATED] to avoid the use of equivalence statements, 2018-05-31
!            -> suppressed all local variables

#define SPLOUT 0

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip

      use const_mod, only: zero

      use para0_mod, only: imax, jmax
      use bloc0_mod
      use bloc_mod

#if ( SPLOUT == 1 )
!--pour la sortie sur fichier "splout" des variations de eta pendant
!   la derniere iteration :
      use reper_mod
#endif

!! END_OF_USE_SECTION



      real(kind=dblp) :: corfac, afdtb
      integer(kind=ip):: i, j, k

!dmr [NOEQUI] Code below was the equivalence statement
!dmr           2018-05-31 proposal to use phihhh directly without equivalence
!dmr [NOTA] phihhh a pour dimensions imax,jmax,6 ...


! #if ( SPLOUT == 1 )
!--pour la sortie sur fichier "splout" des variations de eta pendant
!   la derniere iteration :
! [SCRPTCOM] #include "reper.com"
!     include 'split1.com'
! #endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initialisation.                                                 |
!-----------------------------------------------------------------------

      do j=1,jmax
       do i=1,imax
        etaspl(i,j)   = 0.0
        ubspl(i,j)    = 0.0
        vbspl(i,j)    = 0.0
        phihhh(i,j,1) = 0.0 ! tm1x2
        phihhh(i,j,2) = 0.0 ! tm2x2
        phihhh(i,j,5) = 0.0 ! phiypx
        phihhh(i,j,6) = 0.0 ! phiymx
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Preparation de l'integration du mode barotrope .                |
!-----------------------------------------------------------------------

      corfac = 2.*txifcb*dtb

!--Mise en place des masques "interface" :
      j = 1+js2
      do i=ims1,ims2
        phihhh(i,j,2) = tms(i,j-1,ks2) * tms(i,j,ks2)
     &             * 2.0 * cmx(i,j,2) * cmy(i,j,2)
      enddo

!--Calcul de la Tension de fond : <- transfere dans "flucor".

!--Debut de la boucle externe sur l'indice de latitude j :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,k,afdtb)
      do j=js1,js2
!-----

      do i=ims1,ims2+1
        phihhh(i,j,1) = tms(i-1,j,ks2) * tms(i,j,ks2)
     &             * 2.0 * cmx(i,j,1) * cmy(i,j,1)
      enddo

      do i=ims1,ims2
        phihhh(i,j,2) = tms(i,j-1,ks2) * tms(i,j,ks2)
     &             * 2.0 * cmx(i,j,2) * cmy(i,j,2)
      enddo

!--Somme les termes de forcage barocline dans fub(ks2) & fvb(ks2) :
!- surface et fond :
      do i=iu1(j),iu2(j)
        fub(i,j,ks2) = fub(i,j,ks2) + phifu(i,j) - phisu(i,j)
        fvb(i,j,ks2) = fvb(i,j,ks2) + phifv(i,j) - phisv(i,j)
!       fub(i,j,ks2) = fub(i,j,ks2) - phisu(i,j)
!       fvb(i,j,ks2) = fvb(i,j,ks2) - phisv(i,j)
      enddo
      do k=ks1,ks2-1
       do i=iu1(j),iu2(j)
        fub(i,j,ks2) = fub(i,j,ks2) + fub(i,j,k)
        fvb(i,j,ks2) = fvb(i,j,ks2) + fvb(i,j,k)
       enddo
      enddo

!--determinant, resolution semi Implic. Corriolis :
      do i=iu1(j),iu2(j)
        afdtb = fs2cor(i,j)*corfac
        phihhh(i,j,3) = tmu(i,j,ku2) / (1.0 + afdtb*afdtb ) ! unsdet
      enddo

      if (ahe.ne.zero) then
!--Filtre : calcul des flux diagonaux et N-S (issu de "conti2d"), 1ere ss_iter.
        do i=isf1(j),isf2(j)
          phihhh(i,j,4) = (eta(i,j-1)-eta(i,j)) * phihhh(i,j,2) ! phiypy
        enddo
        do i=iu1(j),iu2(j)+1
          phihhh(i,j,5)=(eta(i-1,j-1)-eta(i,j))*tmu(i,j,ks2)*cmxy(i,j,3)
          phihhh(i,j,6)=(eta(i,j-1)-eta(i-1,j))*tmu(i,j,ks2)*cmxy(i,j,3)
        enddo
      endif

!--Fin de la 1ere boucle externe sur l'indice de latitude j .
      enddo

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine barot -
      end subroutine barot
