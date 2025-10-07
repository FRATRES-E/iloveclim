!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009

      SUBROUTINE sepgl(vinp1, vinp2, voutw, vouta,
     &                 spvbin, spvout, ioutw, joutw, iouta, jouta,
     &                 kmniv, kref, kexcl, ijdl, jeqm, ksgn)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!   Appele par "ncdfout".
!   Separation des 2 grilles dans 2 tableaux voutw & vouta
!-------
!  modif : 09/08/94

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip, snlp=>sp
      
      use para0_mod
      use para_mod
      use reper_mod
!! END_OF_USE_SECTION

!--dummy variables :
      real*4 spvbin, spvout
      
      real(kind=dblp), dimension(imax,jmax)  :: vinp1, vinp2, kmniv
      real(kind=snlp), dimension(ioutw,joutw)::voutw
      real(kind=snlp), dimension(iouta,jouta)::vouta

      integer(kind=ip):: iouta, ioutw, jouta, joutw

!--- loads of locales
      integer(kind=ip):: i, ia, iim, ijdl, j, ja, jeq, jeqm, jj, jjm
     &                 , kexcl, kref, ksgn


!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  0 ) Initialisation des tableaux avec "spvout" :
!-----------------------------------------------------------------------

      do j=1,joutw
       do i=1,ioutw
        voutw(i,j) = spvout
       enddo
      enddo

      do j=1,jouta
       do i=1,iouta
        vouta(i,j) = spvout
       enddo
      enddo

      if (ltest.le.2) then
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  1 ) traitement du tableau grille WW , simple transfert :
!-----------------------------------------------------------------------

      if (ksgn.eq.2) then
!--Traitement particulier pour w :
        do j=1,joutw
         do i=1,ioutw
          if(kmniv(i,j).le.kref) then
            voutw(i,j) = vinp1(i,j) - vinp2(i,j)
          elseif(kmniv(i,j).le.kexcl) then
            voutw(i,j) = spvbin
          endif
         enddo
        enddo
      else
        do j=1,joutw
         do i=1,ioutw
          if(kmniv(i,j).le.kref) then
            voutw(i,j) = vinp1(i,j)
          elseif(kmniv(i,j).le.kexcl) then
            voutw(i,j) = spvbin
          endif
         enddo
        enddo
      endif

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  2 ) traitement du tableau grille AA , simple transfert :
!-----------------------------------------------------------------------

      if (ksgn.eq.2) then
!--Traitement particulier pour w :
        do ja=1,jouta
         j = ja + jsepar - 1
         do i=1,iouta
          if(kmniv(i,j).le.kref) then
            vouta(i,ja) = vinp1(i,j) - vinp2(i,j)
          elseif(kmniv(i,j).le.kexcl) then
            vouta(i,ja) = spvbin
          endif
         enddo
        enddo
      else
        do ja=1,jouta
         j = ja + jsepar - 1
         do i=1,iouta
          if(kmniv(i,j).le.kref) then
            vouta(i,ja) = vinp1(i,j)
          elseif(kmniv(i,j).le.kexcl) then
            vouta(i,ja) = spvbin
          endif
         enddo
        enddo
      endif

!--
      else
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  3 ) traitement du tableau grille WW , transfert si j =< jsep(i) :
!-----------------------------------------------------------------------

      jeq   = jeqm + 2

      if (ksgn.eq.2) then
!--Traitement particulier pour w :
        do i=1,ioutw
!        jjm = min(jsep(i), joutw)
         jjm = min(max(jsep(i)-1, jeq), joutw)
         do j=1,jjm
          if(kmniv(i,j).le.kref) then
            voutw(i,j) = vinp1(i,j) - vinp2(i,j)
          elseif(kmniv(i,j).le.kexcl) then
            voutw(i,j) = spvbin
          endif
         enddo
        enddo
      else
        do i=1,ioutw
         jjm = min(max(jsep(i)-1, jeq), joutw)
         do j=1,jjm
          if(kmniv(i,j).le.kref) then
            voutw(i,j) = vinp1(i,j)
          elseif(kmniv(i,j).le.kexcl) then
            voutw(i,j) = spvbin
          endif
         enddo
        enddo
      endif

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  4 ) traitement du tableau grille AA , retournement :
!-----------------------------------------------------------------------

      if (ksgn.eq.2) then
!--Traitement particulier pour w :
        do ja=1,jouta
         i = ijdl - ja
         jj = jsep(i) - jeq
         iim = jj + min( 2 , jj + 1 )
         do ia=iim,iouta
           j = jeqm + ia
           if (kmniv(i,j).le.kref) then
             vouta(ia,ja) = vinp1(i,j) - vinp2(i,j)
           elseif(kmniv(i,j).le.kexcl) then
             vouta(ia,ja) = spvbin
           endif
         enddo
        enddo
      elseif (ksgn.eq.-1) then
!--Avec changement de signe :
        do ja=1,jouta
         i = ijdl - ja
         jj = jsep(i) - jeq
         iim = jj + min( 2 , jj + 1 )
         do ia=iim,iouta
           j = jeqm + ia
           if (kmniv(i,j).le.kref) then
             vouta(ia,ja) = -vinp2(i,j)
           elseif(kmniv(i,j).le.kexcl) then
             vouta(ia,ja) = spvbin
           endif
         enddo
        enddo
      else
!--Sans changement de signe :
        do ja=1,jouta
         i = ijdl - ja
         jj = jsep(i) - jeq
         iim = jj + min( 2 , jj + 1 )
         do ia=iim,iouta
           j = jeqm + ia
           if (kmniv(i,j).le.kref) then
             vouta(ia,ja) = vinp2(i,j)
           elseif(kmniv(i,j).le.kexcl) then
             vouta(ia,ja) = spvbin
           endif
         enddo
        enddo
      endif

!--
      endif

      return
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!- fin de la routine sepgl -
      end
