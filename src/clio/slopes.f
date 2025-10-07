!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009

      SUBROUTINE slopes(scathk, alphhk, cdtsxz, kslpdw,  ns)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  appelee par "scale" ; pour le scalaire "ns" ,
! Calcule et integre ds scadt le courant de DownSloping .
!  modif : 16/05/96

!! START_OF_USE_SECTION

      use const_mod

      use para0_mod, only: imax, ixjmax, kmax, nlpmax
      use bloc0_mod, only: scal, xslop
      use bloc_mod,  only: smxy, nslp, nslpdw, ijslp, kslp, lslp
     &            , alpslp, uvcslp
!! END_OF_USE_SECTION

      implicit none

!--variables locales equivalentes
!-- (!dmr only intent(in) && scalhk is used only at ns given in the input call)

      real, dimension(ixjmax) :: smxyhk
      real, dimension(ixjmax,kmax) :: scalhk

!--By reference variables :
      real, dimension(ixjmax,kmax),  intent(inout) :: scathk
      real, dimension(ixjmax,kmax,2),intent(inout) :: alphhk
      real, dimension(kmax),         intent(in)    :: cdtsxz
      integer, dimension(nlpmax),    intent(in)    :: kslpdw
      integer                      , intent(in)    :: ns 

!-- variables declarees suite implicit none
      integer :: k, l, nldw, ij, kk, lk, ijc, nnc, lx, imax1p
      real    :: zzslp, sscal

!- initialise ds scale : cdtsxz(k) =  dts(k) * unsdx * unsdz(k)
!   ATTENTION :  dx = dy  indispensable !

!dmr [NOEQUI]
      smxyhk(:) = RESHAPE(smxy(:,:,0),(/ ixjmax /))
      do k=1,kmax
        scalhk(:,k) = RESHAPE(scal(:,:,k,ns),(/ ixjmax /))
      enddo
!dmr [NOEQUI]

      if (xslop.lt.epsil) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Modification du decentrement (xslop = 0) :
!-----------------------------------------------------------------------

      imax1p = imax + 1
      do  nldw=1,nslpdw
        ij = ijslp(nslp(nldw))
        kk = kslp(nslp(nldw))
        l  = lslp(nslp(nldw))
!- decentrement = 1 :
        lx = mod((imax1p+l),imax) - 1
        ijc = ij + max(-lx,0) + max(l,lx)
        nnc = 2 - abs(lx)
        alphhk(ijc,kk,nnc) = alpslp(nldw)
      enddo

      else
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  6 ) Prise en compte du Flux (xslop > 0 ) :
!-----------------------------------------------------------------------

      imax1p = imax + 1
      do nldw=1,nslpdw
        ij = ijslp(nslp(nldw))
        kk = kslp(nslp(nldw))
        l  = lslp(nslp(nldw))
        
        if (kslpdw(nldw).le.kk) then
!- permutation de la quantite : uvcslp(nl)*dt/dx - de kslpdw a kslp :
          zzslp = smxyhk(ij) * uvcslp(nldw)
          sscal = scalhk(ij+l,kk)

          do k=kslpdw(nldw),kk
            scathk(ij,k) = scathk(ij,k)
     &         + zzslp * cdtsxz(k) * (sscal-scalhk(ij,k))
            sscal = scalhk(ij,k)
          enddo
          
          scathk(ij+l,kk) = scathk(ij+l,kk) + smxyhk(ij+l) *
     &      uvcslp(nldw) * cdtsxz(kk) * (sscal-scalhk(ij+l,kk))

        endif
        
      enddo ! on nldw

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine slopes -
      end subroutine slopes
