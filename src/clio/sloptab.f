!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:54 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:54 CET 2009

      SUBROUTINE sloptab(uslpfx, vslpfx, valfil, kslpdw, nn99)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Appele par "process".
! Retablit dans 2 tableaux imax X jmax X (2+nsmax)
!    la hauteur de chute et les Transports de masse et de scalaire
!    impliques dans le courant de Dowsloping.
!-----
!  modif : 30/01/98

!! START_OF_USE_SECTION

      use const_mod, only:

      use para0_mod, only: imax, nlpmax, ixjmax, jmax, kmax, nsmax
      use para_mod,  only:
      use bloc0_mod, only: dx, scal, z
      use bloc_mod,  only: nslpdw, nslp, ijslp, kslp, lslp, uvcslp

      use global_constants_mod, only: dblp=>dp, ip

!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION


!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"

!! END_OF_INCLUDE_SECTION

!--dummy variables :
      integer(ip), dimension(nlpmax)            :: kslpdw
      real(dblp),  dimension(imax,jmax,-1:nsmax):: uslpfx, vslpfx

!--variables locales equivalentes :
      real(dblp),  dimension(ixjmax,kmax,nsmax) :: scalhk
!dmr [NOEQUI]      equivalence ( scalhk(1,1,1) , scal(1,1,1,1) )

!--variables locales :

       real(dblp)  :: valfil, uvflow
       integer(ip) :: nn99, dzdw, i, iijj, ij, j, k, kdw, l, n
     &              , nl, ns

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initialisation par valfil .
!-----------------------------------------------------------------------

!--Debut du remplissage de u/vslpfx :

      do  n=-1,nsmax
       do  j=1,jmax
        do  i=1,imax
          uslpfx(i,j,n) = valfil
          vslpfx(i,j,n) = valfil
        enddo
       enddo
      enddo

!dmr [NOEQUI]
       scalhk(:,:,:) = RESHAPE(scal,(/ixjmax,kmax,nsmax/))
!dmr [NOEQUI]

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Transfert uvslp -> u/v_slpfx(1) .
!-----------------------------------------------------------------------

      do 250 n=1,nslpdw
        nl = nslp(n)
        kdw = kslpdw(n)
        ij = ijslp(nl)
        k = kslp(nl)
        l = lslp(nl)
        iijj = ij + max(0,l)
        i = 1 + mod(iijj-1,imax)
        j = 1 + (iijj-1) / imax
        dzdw = z(k) - z(kdw)
        uvflow = -l
        uvflow = uvcslp(n) * sign(dx,uvflow)
        if (abs(l).eq.1) then
          uslpfx(i,j,-1) = dzdw
          uslpfx(i,j,0) = uvflow
          do 230 ns=1,nsmax
            uslpfx(i,j,ns) = uvflow
     &                * (scalhk(ij+l,k,ns) - scalhk(ij,k,ns))
 230      continue
        else
          vslpfx(i,j,-1) = dzdw
          vslpfx(i,j,0) = uvflow
          do 240 ns=1,nsmax
            vslpfx(i,j,ns) = uvflow
     &                * (scalhk(ij+l,k,ns) - scalhk(ij,k,ns))
 240      continue
        endif
 250  continue

!--Fin du remplissage de u/vslpfx . ---

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine sloptab -
      end
