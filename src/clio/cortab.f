!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009

      SUBROUTINE cortab(ytab,ymodif,yspv,ymdspv,
     &                  im,is,ie,js,je,kmodif,kspv)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  Appele par "correct",
!  Modifie les valeurs du tableau ytab, entre les indice [is,ie]x[js,je]
!  suivant kmodif : 0,3 -> substituion, 1,4 -> addition ; 2,5 -> multiplication
!     0,1,2 par val.unique ymodif(1,1) ou 3,4,5 val.correspondante ymodif(i,j)
!          kspv   : 0 sans spv, 1 uniquement spv, 2 sauf spv.
!--------
!  modif : 22/09/97



!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

!! END_OF_USE_SECTION


!- dummy variables :

      integer(kind=ip) :: im
      real(kind=dblp),dimension(im,*) :: ytab, ymodif     
      
c~       dimension ytab(im,*), ymodif(im,*)
!- variables locales


      real(kind=dblp)  :: yspv, ymdspv
      integer(kind=ip) :: is, ie, js, je, kmodif, kspv, i, j
      
      logical :: flmod0, flmod1, flspv1, flspv2

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initialisation :
!-----------------------------------------------------------------------

      flmod0 = mod(kmodif,3).eq.0
      flmod1 = mod(kmodif,3).eq.1
      flspv1 = kspv.eq.1
      flspv2 = kspv.eq.2

      if (kmodif.lt.3) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Traitement avec Valeur Unique :
!-----------------------------------------------------------------------

      do j=js,je
       do i=is,ie
         if (ytab(i,j).ne.yspv) then
!--Modifie que si = spv :
           if (flspv1) cycle
         else
!--Ne modifie pas si = spv
           if (flspv2) cycle
         endif
         if (flmod0) then
!--Substitution :
           ytab(i,j) = ymodif(1,1)
         elseif (flmod1) then
!--Addition :
           ytab(i,j) = ytab(i,j) + ymodif(1,1)
         else
!--Multiplication :
           ytab(i,j) = ytab(i,j) * ymodif(1,1)
         endif
       enddo
      enddo

      else
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Traitement avec Valeur correspondante :
!-----------------------------------------------------------------------

      if (kmodif.eq.0) then
!--Substitution :
        do j=js,je
         do i=is,ie
          if (ymodif(i,j).eq.ymdspv) cycle
          if (ytab(i,j).ne.yspv) then
!- Modifie que si = spv :
            if(kspv.eq.1) cycle
          else
!- Ne modifie pas si = spv
            if(kspv.eq.2) cycle
          endif
          ytab(i,j) = ymodif(i,j)
         enddo
        enddo
!-----

      elseif (kmodif.eq.1) then
!--Addition :
        do j=js,je
         do i=is,ie
          if (ymodif(i,j).eq.ymdspv) cycle
          if (ytab(i,j).ne.yspv) then
!- Modifie que si = spv :
            if(kspv.eq.1) cycle
          else
!- Ne modifie pas si = spv
            if(kspv.eq.2) cycle
          endif
          ytab(i,j) = ytab(i,j) + ymodif(i,j)
         enddo
        enddo
!-----

      else
!--Multiplication :
        do j=js,je
         do i=is,ie
          if (ymodif(i,j).eq.ymdspv) cycle
          if (ytab(i,j).ne.yspv) then
!- Modifie que si = spv :
            if(kspv.eq.1) cycle
          else
!- Ne modifie pas si = spv
            if(kspv.eq.2) cycle
          endif
          ytab(i,j) = ytab(i,j) * ymodif(i,j)
         enddo
        enddo
!-----
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

      return

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine cortab -
      end
