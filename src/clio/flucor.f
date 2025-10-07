!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:46 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:46 CET 2009

      SUBROUTINE flucor

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! routine de calcul des corrections de flux (pour les scalaires)
!  qui interviennent au fond de l ocean.
! ( ATTENTION : Pas encore multiplie par le pas de temps ! )
!  modif : 03/07/98

!! START_OF_USE_SECTION

      use const_mod, only:

      use para0_mod, only: imax, jmax, nsmax
      use para_mod,  only:
      use bloc0_mod, only: js1, js2, ks2, spvr, vcor, zfluxm, kfs, is1
     &             , scssv, kniv, scpme, is2, tms, zw, scal, phifs, w
     &             , fss, daeta

      use bloc_mod,  only: deriv, cmxy, unsvol

      use global_constants_mod, only: dblp=>dp, ip

!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"

!! END_OF_INCLUDE_SECTION

!--variables locales :
      real dum,wdum(imax,jmax)
!-dmr
      real(dblp) :: kvol
      integer(ip):: k1000, j50S, k2000

      integer(ip):: i,j,ns,kfd
!-dmr



!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Transfert dans PhiFond de Wfond , puis mise a  zero .           |
!-----------------------------------------------------------------------

!--Stockage dans phifs(-,-,1) de W(kfd) ; Mise a zero de W(kfd) ;
!- et pour la routine informe : Somme de |W(kfd)| stockee dans daeta .
      do j=js1,js2
       do i=is1(j),is2(j)
!dmr&mb --- kfs = local bottom of the ocean (layer index)
        kfd = kfs(i,j)
        phifs(i,j,1) = w(i,j,kfd)
        daeta(i,j) = daeta(i,j) + abs( w(i,j,kfd) )
        w(i,j,kfd) = 0.0
        fss(i,j,0) = fss(i,j,0) + w(i,j,ks2+1)
       enddo
      enddo 
!
!-AM
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1b ) Modification of sea surface w for mass change .           |
!-----------------------------------------------------------------------
!-- stockage dans wdum de w(ks2+1) modifiee
!dmr&mb --- zfluxm is the yearly average of w
      dum=zfluxm*vcor
!     dum=zflux0*vcor <- moyenne instantanee
      do j=js1,js2
        do i=is1(j),is2(j)
          wdum(i,j)=w(i,j,ks2+1)-dum*tms(i,j,ks2)
       enddo
      enddo
!-AM

      do ns=nsmax,1,-1
!--boucle sur les scalaires (ns) :


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Calcul des Flux au fond .
!-----------------------------------------------------------------------

!--calcul de la correction sur le flux entrant par le fond :
      do j=js1,js2
       do i=is1(j),is2(j)
         phifs(i,j,ns) = phifs(i,j,1) * scal(i,j,kfs(i,j),ns)
       enddo 
      enddo 

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Calcul de la correction globale (=la derive) .                  |
!-----------------------------------------------------------------------

!--calcul du gain total entrant (=la derive), suite a la correction :
      deriv(ns) = 0.0
!-AM
!--derive associee au Flux de masse en surf. (= Pluie / Evap.)
!--   -> depend du type de scalaire

!-dmr --- Modify the computation of deriv so that the changes are
!-dmr applied only close to the surface
!-dmr      k1000 = 1 ! 1 is all depths
      k1000 = 10 ! 8 is 1000 meters
      k2000 = 7 ! 5 is 2300 meters, 7 is 1225
      j50S = 13 ! 13 is 64°S
      kvol = 0.0d0
!-dmr
      if (scpme(ns).eq.spvr) then
       do j=js1,js2
         do i=is1(j),is2(j)
           deriv(ns) = deriv(ns) + cmxy(i,j,0) *
     &          ( phifs(i,j,ns) - w(i,j,ks2+1)*scal(i,j,ks2,ns) )
         enddo 
       enddo 
      else
       do j=js1,js2
         do i=is1(j),is2(j)
           deriv(ns) = deriv(ns) + cmxy(i,j,0) *
     &          ( phifs(i,j,ns) - wdum(i,j)*scssv(ns) )
!-dmr
!dmr --- Volume "in metres" on which I limit the computation, not unsvol
!dmr --- The following is equal to unsvol if k1000 = kniv(i,j,0)

           IF (j.LE.j50S) THEN ! only south of 50°S in all basins
             kvol = kvol + zw(max(kniv(i,j,0),k2000)) * cmxy(i,j,0) * (-1.0d0)
           ELSE
             kvol = kvol + zw(max(kniv(i,j,0),k1000)) * cmxy(i,j,0) * (-1.0d0)
           ENDIF
!-dmr
         enddo 
       enddo 
      endif

       if (scpme(ns).eq.spvr) then

!-AM
         deriv(ns) = deriv(ns) * unsvol

!-dmr
      else

        deriv(ns) = deriv(ns) / kvol

      endif
!-dmr


!--fin de la boucle sur les scalaire.
      enddo 

      return

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine flucor -
      end
