!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:51 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:51 CET 2009

      SUBROUTINE redtab(ytab,nmax,numfil,kkr)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Appele par "redforc", "redrunb".
! Lecture sur fich. numfil, en un seul bloc, des nmax elements du tableau ytab.
!  modif : 16/15/96

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip           

!! END_OF_USE_SECTION

!- dummy variables :
      integer(kind=ip) :: nmax, kkr, numfil
      
      real(kind=dblp), dimension(nmax) :: ytab

      kkr = 1
      read(numfil,end=920,err=910) ytab

      return
 910  continue
      kkr = -1
      return
 920  continue
      kkr = 0
      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine redtab -
      end
