!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:52 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:52 CET 2009

      SUBROUTINE savtab(ytab,nmax,numfil,kkw)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Appele par "savrunb".
! Ecriture sur fich. numfil, en un seul bloc, des nmax elements du tableau ytab.
!  modif : 16/15/96

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip           

!! END_OF_USE_SECTION

    
!- dummy variables :
      integer(kind=ip) :: nmax, numfil, kkw
      real(kind=dblp), dimension(nmax):: ytab

      kkw = 1
      write(numfil,err=910) ytab

      return

 910  continue
      kkw = -1
      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine savtab -
      end
