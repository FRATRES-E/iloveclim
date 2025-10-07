!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine sert a  rearranger les tableaux icebergs pour ne 
!       pas utiliser trop d'espace memoire (supprime les cases inutiles
!       d'un grand tableau)
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche  & Marianne Buegelmayer
!      Date   : 23 mai 2012
!      Derniere modification :  31 mai 2012
!-----|--1--------2---------3---------4---------5---------6---------7-|

      SUBROUTINE rearrange_lmx_arrays(livingicb)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : livingicb, nomre d'icebergs vivants
!       Variables de sortie : RAS
!-----|--1--------2---------3---------4---------5---------6---------7-|

      use const_mod

#include "type.com"
#include "para.com"
#include "bloc.com"
#include "reper.com"
#include "iceberg.com"

       INTEGER, INTENT(INOUT) :: livingicb

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       INTEGER :: licb,licb1
       INTEGER :: i,j,aa,bb,lmx_new


!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Business starts ...
!-----|--1--------2---------3---------4---------5---------6---------7-|

       licb = 1
       licb1 = 0

       DO WHILE ((kiceb(licb).LE.-1).AND.(licb.LT.lmx))
           licb = licb + 1  
           licb1 = licb1 + 1  
       ENDDO
      
       IF (licb1.GE.1) THEN

          startxn(1:lmx-licb1) = startxn(licb:lmx)
          startyn(1:lmx-licb1) = startyn(licb:lmx)
          xn(1:lmx-licb1) = xn(licb:lmx)
          yn(1:lmx-licb1) = yn(licb:lmx)
          starthiceb(1:lmx-licb1) = starthiceb(licb:lmx)
          startwiceb(1:lmx-licb1) = startwiceb(licb:lmx)
          hiceb(1:lmx-licb1) = hiceb(licb:lmx)
          wiceb(1:lmx-licb1) = wiceb(licb:lmx)
          uiceb(1:lmx-licb1) = uiceb(licb:lmx)
          viceb(1:lmx-licb1) = viceb(licb:lmx)
          pond_icb(1:lmx-licb1) = pond_icb(licb:lmx)
          termBIGx(1:lmx-licb1) = termBIGx(licb:lmx)
          termBIGy(1:lmx-licb1) = termBIGy(licb:lmx)
          kiceb(1:lmx-licb1) = kiceb(licb:lmx)
          siceb(1:lmx-licb1) = siceb(licb:lmx)
          vol_orig(1:lmx-licb1) = vol_orig(licb:lmx)
!mab: vol_melt used in old way of calculating iceberg melt
!mab: in icebtraj.f
          vol_melt(1:lmx-licb1) = vol_melt(licb:lmx)
          id_abs_table(1:lmx-licb1) = id_abs_table(licb:lmx)
          flagstartwritten(1:lmx-licb1) = flagstartwritten(licb:lmx)

          lmx = lmx-licb1

       ENDIF
            

       END SUBROUTINE rearrange_lmx_arrays
