!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009

      subroutine alphdkc(alphkx, alphky, grdmin, k, ns)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  calcule les taux de decentrement alphax & alphay dans les 2 directions,
!  pour le niveau "k", et les passent par argument (alphkx & alphky).
! -> uniquement cas des Coeff. "alpha" visant a supprimer le Depassement .
!   en sortie : multiplie par |u1.cmy/2|, |v2.cmx/2|
!  modif : 25/05/99
!dmr --- Updated: 2018-11-15

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp
      use const_mod, only: one

      use para0_mod, only: imax, jmax, ixjmax, kmax, nlpmax, nrpmax
     &   , ncomax, ltest
      use para_mod , only: nltmax, nbsmax, ntocn, ntatm

      use bloc0_mod   , only: ks2, js1, js2, is1, is2, isf1, isf2, jcl1, jcl2
     &                , ims1, ims2, dts, unsdx, unsdy, alphmi, alphgr
     &                , u, v, tms, scal
     &                , iberp, iberpm, jberp, jberpm, ibera, jbera, iberam
     &                , jberam
      use bloc_mod    , only: n1coin, n2coin, n3coin, n4coin, i1coin, i2coin
     &                , i3coin, i4coin, cmx, cmy, smx, smy, numit, nstart
      use isoslope_mod, only: uiso, viso
!! END_OF_USE_SECTION

      implicit none

!--By reference variables :

      real(dblp), dimension(imax,jmax), intent(inout) :: alphkx, alphky
      integer,                          intent(in)    :: ns, k
      real(dblp),                       intent(in)    :: grdmin

!--variables locales :
!      dimension alphxx(imax,jmax) , alphyy(imax,jmax)
      real(dblp), dimension(imax,jmax) :: alphxx = 0.0_dblp, alphyy = 0.0_dblp

! dmr Those two are a former equivalence statement
      real(dblp), dimension(imax,jmax) :: phiax, phiay
      real(dblp), dimension(ixjmax)    :: phicx, phicy

!-- variables declarees suite implicit none
      integer :: i,j,n
      real    :: cc2x, cc2y, u1pd2, v2pd2, alpha


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-!
! dmr  Start of main code
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-!

#if ( CTRL_FIRST_ITER == 1 )
      if (numit.eq.nstart) then
       if (k.eq.ks2) then
        write(66,'(A,I3,F10.6,1PE14.6)') ' alphdkc : sans alphah, |A| ;'
     &              //' ns,alphgr,algrmn =', ns, -alphgr(ns), grdmin
       endif
      endif
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Calcul des Coeff. "alpha" visant a supprimer le Depassement .   |
!-----------------------------------------------------------------------

!--calcul des gradients dans les 2 directions :
      do j=js1,js2
       do i=is1(j),1+is2(j)
        phiax(i,j) = (scal(i,j,k,ns) - scal(i-1,j,k,ns))
     &             * tms(i,j,k) * tms(i-1,j,k)
       enddo
      enddo

      do j=js1,1+js2
       do i=isf1(j),isf2(j)
        phiay(i,j) = (scal(i,j,k,ns) - scal(i,j-1,k,ns))
     &             * tms(i,j,k) * tms(i,j-1,k)
       enddo
      enddo

      phicx(:) = RESHAPE(phiax,(/ixjmax/))
      phicy(:) = RESHAPE(phiay,(/ixjmax/))

!--virage dans les coins (4 types de coins) :
!- coin 1 : tmu(i,j) = 1.
      do n=1,n1coin(k)
        phicx( i1coin(n,k)+1 )    = -phicy( i1coin(n,k) )
        phicy( i1coin(n,k)+imax ) = -phicx( i1coin(n,k) )
      enddo
!- coin 2 : tmu(i+1,j) = 1.
      do n=1,n2coin(k)
        phicx( i2coin(n,k) )      =  phicy( i2coin(n,k) )
        phicy( i2coin(n,k)+imax ) =  phicx( i2coin(n,k)+1 )
      enddo
!- coin 3 : tmu(i,j+1) = 1.
      do n=1,n3coin(k)
        phicx( i3coin(n,k)+1 )    =  phicy( i3coin(n,k)+imax )
        phicy( i3coin(n,k) )      =  phicx( i3coin(n,k) )
      enddo
!- coin 4 : tmu(i+1,j+1) = 1.
      do n=1,n4coin(k)
        phicx( i4coin(n,k) )      = -phicy( i4coin(n,k)+imax )
        phicy( i4coin(n,k) )      = -phicx( i4coin(n,k)+1 )
      enddo

      phiax(:,:) = RESHAPE(phicx,(/imax,jmax/))
      phiay(:,:) = RESHAPE(phicy,(/imax,jmax/))

!--calcul des coeff alphxx & alphyy :
      do j=js1,js2
       do i=is1(j),is2(j)
        alphxx(i,j) = -alphgr(ns)
     &              * abs(phiax(i+1,j) - phiax(i,j))
     &            / ( abs(phiax(i+1,j) + phiax(i,j)) + grdmin )

        alphyy(i,j) = -alphgr(ns)
     &              * abs(phiay(i,j+1) - phiay(i,j))
     &            / ( abs(phiay(i,j+1) + phiay(i,j)) + grdmin )
       enddo
      enddo

!--raccord cyclique et autre (alphxx, alphyy) :
#if ( L_TEST > 0 )
!dmr [CLIOpH]      if (ltest.ge.1) then
        do j=jcl1,jcl2
          alphxx(ims1-1,j) = alphxx(ims2,j)
          alphxx(ims2+1,j) = alphxx(ims1,j)
        enddo
!dmr [CLIOpH]      endif
#endif

#if ( L_TEST == 3 )
!dmr [CLIOpH]      if (ltest.eq.3) then
        alphyy(iberpm,jberp) = alphyy(ibera, jberam)
        alphyy(iberp, jberp) = alphyy(iberam,jberam)
        alphyy(iberam,jbera) = alphyy(iberp, jberpm)
        alphyy(ibera, jbera) = alphyy(iberpm,jberpm)
!dmr [CLIOpH]      endif
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Calcul des Taux de decentrement a partir des Coeff "alpha" :    |
!-----------------------------------------------------------------------

      cc2x = 2. * unsdx * dts(k)
!--calcul du taux de decentrement dans la direction x :
      do j=js1,js2
        do i=is1(j),is2(j)+1

          u1pd2 = abs(0.25 * (u(i,j,k) + u(i,j+1,k))
     &               +0.5  *  uiso(i,j,k)            )
          alpha = smx(i,j,1) * cc2x * u1pd2 + alphmi(k)
          alphkx(i,j) = cmy(i,j,1) * u1pd2 *
     &       min(one, max(alpha, alphxx(i,j), alphxx(i-1,j)))
        enddo
      enddo

!--calcul du taux de decentrement dans la direction y :
      cc2y  = 2. * unsdy * dts(k)
      do j=js1,1+js2
        do i=isf1(j),isf2(j)
          v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k))
     &                +0.5  *  viso(i,j,k)            )
          alpha = smy(i,j,2) * cc2y * v2pd2 + alphmi(k)
          alphky(i,j) = cmx(i,j,2) * v2pd2 *
     &       min(one, max(alpha, alphyy(i,j), alphyy(i,j-1)))
        enddo
      enddo

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine alphdkc -
      end subroutine alphdkc
