!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:52 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:52 CET 2009

      SUBROUTINE scadew(scalat,scaldt,alphax,alphay,ns)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  resolution de l'equation d'evolution du scalaire "ns",
!  dans les 2 directions Horizontales, pendant 1 pas de temps.
!--------------------------------------------------
! - Direction X (E-W) d abord, Y (N-S) ensuite -- |
!--------------------------------------------------
!  modif : 31/05/99
!-dmr updated: 2018-11-20

!! START_OF_USE_SECTION

      use const_mod,    only:

      use para0_mod,    only: imax, jmax, kmax 
      use para_mod,     only: 
      use bloc0_mod,    only: u, v, tms, scal, unsdx, ahs, dts, is1, is2, unsdy
     &             , isf1, isf2, ks1, ks2, ims1, ims2, jcl1, jcl2, js1, js2
     &             , iberpm, jberpm, iberp, jberp, ibera, jbera, iberam, jberam
      use bloc_mod,     only: cmy, cmx, smxy, cmxy
      use isoslope_mod, only: uiso, viso
!! END_OF_USE_SECTION

      implicit none

!--By reference variables :
      real, dimension(imax,jmax,kmax) :: scalat, scaldt, alphax, alphay
      integer, intent(in)             :: ns

!--variables locales :

      integer                    :: i, j, k
      real                       :: ccx, cc2x, ccxdif, u1cd2, divx, ccy
     &                            , cc2y, ccydif,v2cd2
      real, dimension(imax,jmax) :: phix, phiy, phidiv


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  0 ) initialisation .
!-----------------------------------------------------------------------

!--Debut de la boucle externe sur l'indice de niveau k :
c~ !$DIR SPP LOOP_PARALLEL
c~ !$DIR SPP LOOP_PRIVATE(i,j,divx,ccx,cc2x,ccxdif,u1cd2,phix)
c~ !$DIR SPP LOOP_PRIVATE( phidiv, ccy,cc2y,ccydif,v2cd2,phiy)
      do k=ks1,ks2
!-----

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Calcul des taux de decentrement dans la direction X  (E-W) .    |
! => transfere dans la routine "alphdec" appelee avant "scadew" ou "scadns"
!--raccord cyclique et autre (alphax,alphay) :
!     call raccord (alphax, 0., 1, 12)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Calcul des Flux advect. & diffus., dans la direction X (E-W) .  |
!-----------------------------------------------------------------------

      ccx  = unsdx * dts(k)
      cc2x = ccx + ccx
      ccxdif = ahs(k) * unsdx
!--calcul des flux dans la direction x :
      do j=js1,js2
        do i=is1(j),is2(j)+1
          u1cd2 = 0.25 * cmy(i,j,1) * ( u(i,j,k) + u(i,j+1,k) )
     &           +0.5  * cmy(i,j,1) *   uiso(i,j,k)

          phidiv(i,j) = cc2x * u1cd2
          phix(i,j) = tms(i,j,k) * tms(i-1,j,k) * (
     &          u1cd2 * (scalat(i,j,k) + scalat(i-1,j,k))
     &        + (alphax(i,j,k) + smxy(i,j,1) * ccxdif )
     &               * (scalat(i-1,j,k) - scalat(i,j,k)) )
        enddo
      enddo
      
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Bilan des Flux advect. & diffus., dans la direction X (E-W) .   |
!-----------------------------------------------------------------------

!--Prise en compte du bilan des flux selon X : Mise a jour de scalat() :
      do j=js1,js2
        do i=is1(j),is2(j)
          divx = scal(i,j,k,ns) *
     &              smxy(i,j,0) * (phidiv(i,j) - phidiv(i+1,j))
          scaldt(i,j,k) = scaldt(i,j,k)  + divx
          scalat(i,j,k) = scalat(i,j,k)
     &            + smxy(i,j,0) * ccx * (phix(i,j) - phix(i+1,j))
     &            - divx
        enddo
      enddo

!--Fin de la 1ere boucle externe sur l'indice de niveau k .
!490  continue

!--raccord cyclique et autre (scalat) :

#if ( L_TEST >= 1 )
!dmr [CLIOpH]      if (ltest.ge.1) then
        do j=jcl1,jcl2
          scalat(ims1-1,j,k) = scalat(ims2,j,k)
          scalat(ims2+1,j,k) = scalat(ims1,j,k)
        enddo
!dmr [CLIOpH]      endif
#endif
#if ( L_TEST == 3 )
!dmr [CLIOpH]      if (ltest.eq.3) then
        scalat(iberpm,jberp,k) = scalat(ibera, jberam,k)
        scalat(iberp, jberp,k) = scalat(iberam,jberam,k)
        scalat(iberam,jbera,k) = scalat(iberp, jberpm,k)
        scalat(ibera, jbera,k) = scalat(iberpm,jberpm,k)
!dmr [CLIOpH]      endif
#endif      

!--raccord cyclique et autre (scalat) :
!     call raccord(scalat(1,1,k), 0., 1, 0)

!--Debut de la 2nd  boucle externe sur l'indice de niveau k :
!C$DIR SPP LOOP_PARALLEL
!C$DIR SPP LOOP_PRIVATE(i,j,ccy,cc2y,ccydif,v2cd2,phiy)
!     do 790 k=ks1,ks2
!-----

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Calcul des taux de decentrement dans la direction Y (N-S) .     |
! => transfere dans la routine "alphdec" appelee avant "scadew" ou "scadns"
!--raccord cyclique et autre (alphax,alphay) :
!     call raccord (alphay, 0., 1, 24)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  6 ) Calcul des Flux advect. & diffus., dans la direction Y (N-S) .  |
!-----------------------------------------------------------------------

      ccy  = unsdy * dts(k)
      cc2y = ccy + ccy
      ccydif = ahs(k) * unsdy
!--calcul des flux dans la direction y :
      do j=js1,1+js2
        do i=isf1(j),isf2(j)
          v2cd2 = 0.25 * cmx(i,j,2) * ( v(i,j,k) + v(i+1,j,k) )
     &           +0.5  * cmx(i,j,2) *   viso(i,j,k)

          phiy(i,j) = tms(i,j,k) * tms(i,j-1,k) * (
     &          v2cd2 * (scalat(i,j-1,k) + scalat(i,j,k))
     &        + (alphay(i,j,k) + cmxy(i,j,2) * ccydif )
     &               * (scalat(i,j-1,k) - scalat(i,j,k))  )
        enddo
      enddo


!--fin du calcul des flux.

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  7 ) Bilan des Flux advect. & diffus., dans la direction Y (N-S) .   |
!-----------------------------------------------------------------------

!--Prise en compte du bilan des flux selon Y : Mise a jour de scalat() :
      do j=js1,js2
        do i=is1(j),is2(j)
          scalat(i,j,k) = scalat(i,j,k)
     &            + smxy(i,j,0) * ccy * (phiy(i,j) - phiy(i,j+1))
        enddo
      enddo


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Fin de la boucle externe sur l'indice de niveau k .
      enddo ! on k

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine scadew -
      end subroutine scadew
