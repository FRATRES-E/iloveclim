!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:42 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:42 CET 2009

      SUBROUTINE alphdec(alphax,alphay,ns)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  calcule les taux de decentrement alphax & alphay dans les 2 directions .
!   en sortie : multiplie par |u1.cmy/2|, |v2.cmx/2|
!  modif : 25/05/99

!! START_OF_USE_SECTION

      use const_mod   , only: epsil, zero, one

      use para0_mod   , only: imax, jmax, kmax  ! ltest
      use para_mod    , only: 
      use bloc0_mod   , only: is1, is2, js1, js2, isf1, isf2, ks1, ks2
     &                , u, v, scal, alphgr, algrmn, alphmi, unsdx, unsdy
     &                , tms, dts, jcl1, jcl2, ims1, ims2, iberpm, jberpm
     &                , iberp, jberp, ibera, jbera, iberam, jberam
      use bloc_mod    , only: cmx, cmy, smx, smy, numit, nstart
      use isoslope_mod, only: uiso, viso
!! END_OF_USE_SECTION

      implicit none

!--By reference variables :
      integer                        , intent(in)  :: ns
      real, dimension(imax,jmax,kmax), intent(out) :: alphax, alphay

!--variables locales :
      real, dimension(imax,jmax) :: phiax, phiay, alphxx, alphyy

!-- variables declarees suite implicit none
      integer :: i,j,k
      real    :: cc2x, cc2y, alpha, u1pd2, v2pd2, alpha1, alpha2, grdmin

#if ( CTRL_FIRST_ITER == 1 )     
      if (numit.eq.nstart) then
       if (alphgr(ns).gt.zero) then
        write(66,'(A,I3,F10.6,1PE14.6)') ' alphdec : sans alphah, |A| ;'
     &            //' ns,alphgr,algrmn =', ns, alphgr(ns), algrmn(ns)
       elseif (alphgr(ns).eq.zero .and. ns.eq.1) then
        write(66,'(A,I3,F10.6,1PE14.6)') ' alphdec : sans alphah ;'
       endif
      endif
#endif      

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  0 ) preparation des coeff. intervenant dans les flux d'advection .  |
!-----------------------------------------------------------------------

      if ( alphgr(ns).eq.zero ) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1a) Decentrement Minimum pour stabilite = Nb.de courant :
!-----------------------------------------------------------------------

!--Debut de la boucle externe sur l'indice de niveau k :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,j,cc2x,u1pd2,alpha,cc2y,v2pd2)
      do k=1,kmax
!-----

      cc2x  = 2. * unsdx * dts(k)
!--calcul du taux de decentrement dans la direction x :
      do j=js1,js2
        do i=is1(j),is2(j)+1
          u1pd2 = abs( 0.25 * (u(i,j,k) + u(i,j+1,k))
     &                + 0.5 *  uiso(i,j,k)            )
          alpha = smx(i,j,1) * cc2x * u1pd2 + alphmi(k)
          alphax(i,j,k) = cmy(i,j,1) * u1pd2 * min(one, alpha)
        enddo
       enddo

!--calcul du taux de decentrement dans la direction y :
      cc2y  = 2. * unsdy * dts(k)
      do j=js1,1+js2
        do i=isf1(j),isf2(j)
          v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k))
     &                + 0.5 *  viso(i,j,k)            )
          alpha = smy(i,j,2) * cc2y * v2pd2 + alphmi(k)
          alphay(i,j,k) = cmx(i,j,2) * v2pd2 * min(one, alpha)
       enddo
      enddo

!--Fin de la boucle externe sur l'indice de niveau k .
      enddo ! on k

      return

      elseif ( alphgr(ns).eq.2.0d0 ) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1b) Schema Upwind :
!-----------------------------------------------------------------------

!--Debut de la boucle externe sur l'indice de niveau k :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,j,cc2x,u1pd2,alpha,cc2y,v2pd2)
      do k=1,kmax

!--calcul du taux de decentrement dans la direction x :
      do j=js1,js2
        do i=is1(j),is2(j)+1
          u1pd2 = abs( 0.25 * (u(i,j,k) + u(i,j+1,k))
     &                + 0.5 *  uiso(i,j,k)            )
          alphax(i,j,k) = cmy(i,j,1) * u1pd2
       enddo
      enddo

!--calcul du taux de decentrement dans la direction y :
      do j=js1,1+js2
        do i=isf1(j),isf2(j)
          v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k))
     &                + 0.5 *  viso(i,j,k)            )
          alphay(i,j,k) = cmx(i,j,2) * v2pd2
       enddo
      enddo

!--Fin de la boucle externe sur l'indice de niveau k .
      enddo ! on k
!-----
      return

      elseif ( alphgr(ns).gt.zero) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Calcul des Coeff. "alpha" : formule "TVM" .                     |
!-----------------------------------------------------------------------

      grdmin = max( epsil, algrmn(ns) )

!--Debut de la boucle externe sur l'indice de niveau k :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,j,phiax,phiay,alpha1,alpha2,alphxx,alphyy)
!$DIR SPP LOOP_PRIVATE(cc2x,u1pd2,alpha,cc2y,v2pd2)
      do k=ks1,ks2
!----

!--calcul des gradients dans les 2 directions :
      do j=js1,js2
       do i=is1(j),1+is2(j)
        phiax(i,j) = abs(scal(i,j,k,ns) - scal(i-1,j,k,ns))
       enddo
      enddo
       
      do j=js1,1+js2
       do i=isf1(j),isf2(j)
        phiay(i,j) = abs(scal(i,j,k,ns) - scal(i,j-1,k,ns))
        enddo
      enddo

!--calcul des coeff alphax & alphyy :
      do j=js1,js2
       do i=is1(j),is2(j)
        alpha1 = ( phiax(i+1,j) - phiax(i,j) )
     &         / ( phiax(i+1,j) + phiax(i,j) + grdmin )
        alpha2 = ( phiay(i,j+1) - phiay(i,j) )
     &         / ( phiay(i,j+1) + phiay(i,j) + grdmin )
        alphxx(i,j) = alphgr(ns) * abs(alpha1)
     &              * tms(i+1,j,k) * tms(i,j,k) * tms(i-1,j,k)
        alphyy(i,j) = alphgr(ns) * abs(alpha2)
     &              * tms(i,j+1,k) * tms(i,j,k) * tms(i,j-1,k)
        enddo
      enddo

!     do 210 j=js1,js2
!      do 210 i=is1(j),is2(j)
!       dp = abs( scal(i+1,j,k,ns) - scal(i,j,k,ns) )
!       dm = abs( scal(i,j,k,ns) - scal(i-1,j,k,ns) )
!       alpha = ( dp - dm ) / ( dp + dm + epsil )
!       alpha = alphgr(ns) * alpha * alpha
!C      alpha = alphgr(ns) * abs( alpha )
!       alphxx(i,j) = alpha * tms(i+1,j,k)*tms(i,j,k)*tms(i-1,j,k)
!210  continue
!     do 230 j=js1,js2
!      do 230 i=is1(j),is2(j)
!       dp = abs( scal(i,j+1,k,ns) - scal(i,j,k,ns) )
!       dm = abs( scal(i,j,k,ns) - scal(i,j-1,k,ns) )
!       alpha = (dp - dm ) / ( dp + dm + epsil )
!       alpha = alphgr(ns) * alpha * alpha
!C      alpha = alphgr(ns) * abs( alpha )
!       alphyy(i,j) = alpha * tms(i,j+1,k)*tms(i,j,k)*tms(i,j-1,k)
!230  continue

!--raccord cyclique et autre (alphxx, alphyy) :
#if ( L_TEST > 1 )
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
!is0      u1pd2 = abs( 0.25 * (u(i,j,k) + u(i,j+1,k)) )
          u1pd2 = abs( 0.25 * (u(i,j,k) + u(i,j+1,k))
     &                +0.5  *  uiso(i,j,k)            )
          alpha = smx(i,j,1) * cc2x * u1pd2 + alphmi(k)
!         cnx = smx(i,j,1) * cc2x * u1d2
!         alpha = abs(cnx) + alphmi(k)
!    &          + alphah(1) * (2.0 - tmu(i,j,k) - tmu(i,j+1,k))
          alphax(i,j,k) = cmy(i,j,1) * u1pd2 *
     &       min(one, max(alpha, alphxx(i,j), alphxx(i-1,j)))
        enddo
      enddo

!--calcul du taux de decentrement dans la direction y :
      cc2y  = 2. * unsdy * dts(k)
      do j=js1,1+js2
        do i=isf1(j),isf2(j)
!is0      v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k)) )
          v2pd2 = abs( 0.25 * (v(i,j,k) + v(i+1,j,k))
     &                + 0.5 *  viso(i,j,k)            )
          alpha = smy(i,j,2) * cc2y * v2pd2 + alphmi(k)
!         cny = smy(i,j,2) * cc2y * v2d2
!         alpha = abs(cny) + alphmi(k)
!    &          + alphah(2) * (2.0 - tmu(i,j,k) - tmu(i+1,j,k))
          alphay(i,j,k) = cmx(i,j,2) * v2pd2 *
     &       min(one, max(alpha, alphyy(i,j), alphyy(i,j-1)))
        enddo
      enddo

!--Fin de la boucle externe sur l'indice de niveau k .
      enddo ! on k
      else
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Calcul des Coeff. "alpha" visant a supprimer le Depassement .   |
!-----------------------------------------------------------------------

      grdmin = max( epsil, algrmn(ns) )

!--Debut la de boucle externe sur l'indice de niveau k :
!$DIR SPP LOOP_PARALLEL
      do k=ks1,ks2
!----
        call alphdkc(alphax(1,1,k),alphay(1,1,k),grdmin,k,ns)
!--Fin de la boucle externe sur l'indice de niveau k .
       enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine alphdec -
      end subroutine alphdec
