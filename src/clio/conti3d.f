!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:43 CET 2009

      SUBROUTINE conti3d

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- calcul de la vitesse verticale w depuis la surface jusqu au fond.
!  modif : 06/10/98


#if ( COMATM == 1 )
      use comemic_mod, only:
#endif

      use const_mod, only:

      use para0_mod, only: imax, jmax, kmax
      use para_mod,  only:
      use bloc0_mod, only: is1, is2, js1, js2, ks1, ks2, zurfow, zfluxm,
     &                     u, v, w, dz, tms, zflux0, unsdy, unsdx, dts,
     &                     phiss
#if( APPLY_UNCORFWF == 1)
     &                     , w_uncor
#endif

      use bloc_mod,  only: aire, cmx, cmy, smxy

      use global_constants_mod, only: dblp=>dp, ip

      use newunit_clio_mod, only: mouchard_id !clio3_out_id ! used in commented dbg
      
      implicit none

!-AM - pour LOCH: (u,v,w) au meme instant que T & S
!      dynam3 en common avec loch.com
      real*8 uloch,vloch,wloch
      common /dynam3/uloch(imax,jmax,kmax),vloch(imax,jmax,kmax),
     &               wloch(imax,jmax,kmax)
!-AM - pour dilution icesheets


       real(kind=dblp) :: dzd2, unsdts
       integer(kind=ip) :: i, j, k


!--variables locales equivalentes (pour raison de place memoire) :
      real(dblp), dimension(imax) :: phixj
!     dimension phix(imax,jmax), phiy(imax,jmax)

!-AM - terme rappel sur la sal. pour FW Flx -> voir thersf.f


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!     [step/s] = 1 / [s/step]
      unsdts = 1.0 / dts(ks2)

!--Debut de la boucle externe sur l'indice de latitude j :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,k,phixj,dzd2)
      do j=js1,js2
!-----

!--initialisation :
      do i=is1(j),is2(j)
!       w(i,j,ks1) = 0.0
!dbug     if (abs(scal(i,j,ks2,2)).lt.epsil) then
!dbug      write(clio3_out_id,*) 'ARRET : conti3d, scal(i,j,ks2,2) too small'
!dbug      write(clio3_out_id,*) 'scal(i=',i,',j=',j,')=',scal(i,j,ks2,2)
!dbug      stop
!dbug     endif
!-
!-AM
! vertical velocity: [m/s] = [1]*[step/s]*[m/step]
        w(i,j,ks2+1) = tms(i,j,ks2) * unsdts * phiss(i,j,0)
!-
!- phimnx(0,1) => limits Fresh Water Flux : Unit = m/s !!
!ic0    w(i,j,ks2+1) =
!ic0 &    max( phimnx(i,j,0,0), min(phimnx(i,j,1,0),w(i,j,ks2+1)) )
!fd0    w(i,j,ks2+1) = 0.0
        enddo
        enddo

!     mean value of w
      zflux0=0.0
      do j=js1,js2
         do i=is1(j),is2(j)
            zflux0=zflux0+aire(i,j)*tms(i,j,ks2)*w(i,j,ks2+1)
         enddo
      enddo
      zflux0=zflux0/zurfow
      write(mouchard_id,*) 'zflux',zflux0,w(15,15,ks2+1),zfluxm
!
!---> Modifications AM ---> see flucor
!
#if( APPLY_UNCORFWF == 1)
!dmr&mb --- Apply uncorrected freshwater flux  > changes volume of the ocean
!     [m/s] = [m/s] + [1] * [step/s] * [m/step]
      WHERE(abs(w_uncor).GT.0.0d0)
         w(:,:,ks2+1) = w(:,:,ks2+1)
     &        + tms(:,:,ks2) * unsdts * w_uncor(:,:)
      ENDWHERE
      write(mouchard_id,'(A,F12.8,A)') 'UNCORFLOW2: ',
     &     zflux0-(sum(w(:,:,ks2)*tms(:,:,ks2)*aire(:,:))/zurfow),'[m]'
#endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Bilan de l'Eq. Continuite et Somme directement dans w .         |
!-----------------------------------------------------------------------

      do j=js1,js2
!--boucle sur l'indice de niveau, "k" :
      do k=ks2,ks1,-1

      do i=is1(j),is2(j)+1
        phixj(i) = cmy(i,j,1) * (u(i,j+1,k)+u(i,j,k))
      enddo
!       phiy(i,j) = cmx(i,j,2) * (v(i+1,j,k)+v(i,j,k))

      dzd2 = 0.5 * dz(k)
      do i=is1(j),is2(j)
        w(i,j,k) = w(i,j,k+1) + dzd2 * smxy(i,j,0) *
!    &           ( unsdx * (phix(i+1,j)-phix(i,j))
!    &           + unsdy * (phiy(i,j+1)-phiy(i,j)) )
     &           ( unsdx * (phixj(i+1)-phixj(i))
     &           + unsdy * ( cmx(i,j+1,2)*(v(i+1,j+1,k)+v(i,j+1,k))
     &                     - cmx(i, j, 2)*(v(i+1, j, k)+v(i, j, k)) ))
!       w(i,j,k+1) = w(i,j,k+1) * tms(i,j,k)
        w(i,j,k) = w(i,j,k) * tms(i,j,k)
       enddo

!-----
      enddo
!--fin boucle sur "k".


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Fin de la boucle externe sur l'indice de latitude j .
      enddo


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Prepare Flux en surface .                                       |
!-----------------------------------------------------------------------

!--Valeur du scalaire associee au Flux de masse en surf. (= Pluie / Evap.)
!ic0  do ns=1,nsmax
!ic0    if (scpme(ns).eq.spvr) then
!ic0      do j=js1,js2
!ic0       do i=is1(j),is2(j)
!ic0         scs(i,j,ns) = scal(i,j,ks2,ns)
!ic0       enddo
!ic0      enddo
!ic0    endif
!icO  enddo


#if ( CYCC == 1 )
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!-AM - pour LOCH: (u,v,w) au meme instant que T & S
!-----------------------------------------------------------------------
      do k=1,kmax
       do j=1,jmax
        do i=1,imax
         uloch(i,j,k)=u(i,j,k)
         vloch(i,j,k)=v(i,j,k)
         wloch(i,j,k)=w(i,j,k)
        enddo
       enddo
      enddo

#endif

      return

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine conti3d -
      end
