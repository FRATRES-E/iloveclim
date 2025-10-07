!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:42 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:42 CET 2009

      SUBROUTINE alph2dc(alphax,alphay,ns)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  calcule les taux de decentrement alphax & alphay dans les 2 directions .
!   en sortie : multiplie par |u1.cmy/2|, |v2.cmx/2|
!  generalisation 2 D du calcul "suppression du depassement".
!  modif : 30/12/97
!  ppmodif: 19-03-97: gm90
!dmr --- update: 2018-11-14

!! START_OF_USE_SECTION

      use const_mod   , only: epsil, zero

      use para0_mod   , only: imax, jmax, kmax ! ltest
      use para_mod    , only: 
      use bloc0_mod   , only: alphgr, algrmn, alphmi, unsdx, unsdy, dts
     &                , ks1, ks2, js1, js2, is1, is2
     &                , u, v, scal, tms, isf1, isf2, jcl1, jcl2, ims1, ims2
     &                , iberp, jberp, ibera, jbera
      use bloc_mod    , only: numit, nstart, cmx, cmy, smx, smy
      use isoslope_mod, only: uiso, viso

      use newunit_clio_mod, only: clio3_out_id
      
!! END_OF_USE_SECTION

      implicit none

!--By reference variables :
      integer                        , intent(in)  :: ns
      real, dimension(imax,jmax,kmax), intent(out) :: alphax, alphay

!--variables locales :
      real, dimension(imax,jmax) :: grdx1, grdy2, u1cd2, v2cd2, unsdiv
     &    , u0ct, v0ct, phx0e, phx0w, phy0s, phy0n


!-- variables declarees suite implicit none
      integer :: i,j,k
      real    :: alpgrs, grdmin, u0cdiv, v0cdiv, grdx0d, grdy0d, phx0t
     &         , phy0t, alpmng, cc2xg, cc2yg, u1pcd2, v2pcd2, alpha

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  0 ) preparation des coeff. intervenant dans les flux d'advection .  |
!-----------------------------------------------------------------------

      alpgrs = -1. - alphgr(ns)
      grdmin = max( epsil, algrmn(ns) )
      if(numit.eq.nstart) then
       if ( alpgrs.lt.epsil ) then
        write(clio3_out_id,*) 'ARRET dans "alph2dc" ; param alphgr mauvais !'
        stop
       endif

        write(clio3_out_id,
     &             '(A,I3,F10.6,1PE14.6)') ' alph2dc : sans alphah, |A| ;'
     &             //' ns,alphgr,algrmn =', ns, alpgrs, grdmin
      endif

!--Debut de la boucle externe sur l'indice de niveau k :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,j,phx0w,phx0e,phy0s,phy0n)
!$DIR SPP LOOP_PRIVATE(u1cd2,v2cd2,grdx1,grdy2,u0ct,v0ct)
!$DIR SPP LOOP_PRIVATE(u0cdiv,v0cdiv,unsdiv,grdx0d,grdy0d,phx0t,phy0t)
!$DIR SPP LOOP_PRIVATE(alpmng,cc2xg,cc2yg,u1pcd2,v2pcd2,alpha)
!$DIR SPP LOOP_PRIVATE(alphyn,alphys,alphxe,alphxw)
      do k=ks1,ks2
!-----

!--Initialisation :
      do j=1,jmax
       do i=1,imax
        phx0w(i,j) = 0.
        phx0e(i,j) = 0.
        phy0s(i,j) = 0.
        phy0n(i,j) = 0.
       enddo
      enddo
       
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Calcul des Coeff. "alpha" visant a supprimer le Depassement .   |
!-----------------------------------------------------------------------

!--calcul des vitesses(flux) u1 & v2 et des gradients dans les 2 directions :
!- Attention : suppose dx = dy (si non, doivent etre pris en compte)

      do j=js1,js2
       do i=is1(j),1+is2(j)
        u1cd2(i,j) = 0.25 * cmy(i,j,1) * (u(i,j,k) + u(i,j+1,k))
     &              +0.5  * cmy(i,j,1) *  uiso(i,j,k)
        grdx1(i,j) = (scal(i,j,k,ns) - scal(i-1,j,k,ns))
     &             * tms(i,j,k) * tms(i-1,j,k) * smx(i,j,1)
       enddo
      enddo
      
      do j=js1,1+js2
       do i=isf1(j),isf2(j)
        v2cd2(i,j) = 0.25 * cmx(i,j,2) * (v(i,j,k) + v(i+1,j,k))
     &              +0.5  * cmx(i,j,2) *  viso(i,j,k)
        grdy2(i,j) = (scal(i,j,k,ns) - scal(i,j-1,k,ns))
     &             * tms(i,j,k) * tms(i,j-1,k) * smy(i,j,2)
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--repartition des flux :
      do j=js1,js2
       do i=is1(j),is2(j)

        u0ct(i,j) = min(zero, max(u1cd2(i,j), u1cd2(i+1,j)) )
     &            + max(zero, min(u1cd2(i,j), u1cd2(i+1,j)) )
        v0ct(i,j) = min(zero, max(v2cd2(i,j), v2cd2(i,j+1)) )
     &            + max(zero, min(v2cd2(i,j), v2cd2(i,j+1)) )
        u0cdiv = u1cd2(i,j) - u1cd2(i+1,j)
        v0cdiv = v2cd2(i,j) - v2cd2(i,j+1)
        unsdiv(i,j) = 1. / sign( max(epsil,abs(u0cdiv),abs(v0cdiv)),
     &                           (u0cdiv - v0cdiv) )

        phx0w(i,j) = u1cd2( i, j) - u0ct(i,j)
        phx0e(i,j) = u1cd2(i+1,j) - u0ct(i,j)
        phy0s(i,j) = v2cd2(i, j ) - v0ct(i,j)
        phy0n(i,j) = v2cd2(i,j+1) - v0ct(i,j)
        
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--calcul des coeff phx/y-e/w/n/s :
      do j=js1,js2
       do i=is1(j),is2(j)
       
!--grad_diag :
        grdx0d = unsdiv(i,j) *
     &    ( phx0w(i,j) * grdx1(i,j) + phx0e(i,j) * grdx1(i+1,j) )
        grdy0d = unsdiv(i,j) *
     &    ( phy0s(i,j) * grdy2(i,j) + phy0n(i,j) * grdy2(i,j+1) )
!--coef alpha transverse, |A| :

        phx0t = abs( u0ct(i,j)*(grdx1(i+1,j) - grdx1(i,j)) )
     &         / ( grdmin + abs(grdx1(i+1,j) + grdx1(i,j)) )
        phy0t = abs( v0ct(i,j)*(grdy2(i,j+1) - grdy2(i,j)) )
     &         / ( grdmin + abs(grdy2(i,j+1) + grdy2(i,j)) )

        phx0w(i,j) = phx0t + abs( phx0w(i,j)*(grdx1( i, j) - grdy0d) )
     &                       / ( grdmin + abs(grdx1( i, j) + grdy0d) )
        phx0e(i,j) = phx0t + abs( phx0e(i,j)*(grdx1(i+1,j) + grdy0d) )
     &                       / ( grdmin + abs(grdx1(i+1,j) - grdy0d) )
        phy0s(i,j) = phy0t + abs( phy0s(i,j)*(grdy2(i, j ) + grdx0d) )
     &                       / ( grdmin + abs(grdy2(i, j ) - grdx0d) )
        phy0n(i,j) = phy0t + abs( phy0n(i,j)*(grdy2(i,j+1) - grdx0d) )
     &                       / ( grdmin + abs(grdy2(i,j+1) + grdx0d) )
     
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--(demi) raccord cyclique pour phx_e/w :
      do j=jcl1,jcl2
        phx0e(ims1-1,j) = phx0e(ims2,j)
        phx0w(ims2+1,j) = phx0w(ims1,j)
      enddo

#if ( L_TEST == 3 )
!dmr [CLIOpH]      if (ltest.eq.3) then
!--Bering : (demi) raccord pour phy_n/s :
        phy0s(iberp-1,jberp) = phy0n( ibera, jbera-1)
        phy0s( iberp, jberp) = phy0n(ibera-1,jbera-1)
        phy0s(ibera-1,jbera) = phy0n( iberp, jberp-1)
        phy0s( ibera, jbera) = phy0n(iberp-1,jberp-1)
!dmr [CLIOpH]      endif
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Calcul des Taux de decentrement a partir des Coeff "alpha" :    |
!-----------------------------------------------------------------------

      alpmng = alphmi(k) / alpgrs
!--calcul du taux de decentrement dans la direction x :
      cc2xg = 2. * unsdx * dts(k) / alpgrs
      do j=js1,js2
       do i=is1(j),is2(j)+1
         u1pcd2 = abs(u1cd2(i,j))
         alpha =  u1pcd2 *
     &     ( alpmng + smx(i,j,1) * smy(i,j,1) * cc2xg * u1pcd2 )

         alphax(i,j,k) = min( u1pcd2,
     &          alpgrs * max(alpha, phx0e(i-1,j), phx0w(i,j)) )

       enddo
      enddo

!--calcul du taux de decentrement dans la direction y :
      cc2yg  = 2. * unsdy * dts(k) / alpgrs
      do j=js1,1+js2
       do i=isf1(j),isf2(j)
         v2pcd2 = abs(v2cd2(i,j))
         alpha =  v2pcd2 *
     &      ( alpmng + smx(i,j,2) * smy(i,j,2) * cc2yg * v2pcd2 )

         alphay(i,j,k) = min( v2pcd2,
     &          alpgrs * max(alpha, phy0n(i,j-1), phy0s(i,j)) )

       enddo
      enddo

!--Fin de la boucle externe sur l'indice de niveau k .
      enddo ! on k

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine alph2dc -
      end subroutine alph2dc
