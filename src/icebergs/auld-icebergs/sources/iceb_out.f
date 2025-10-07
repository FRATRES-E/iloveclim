      subroutine iceb_out
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  routine qui transforme ce qui faut pour tenir compte du tableau VARIABLE
c  des icebergs
c  terme12345=p=barocline,b=barotrope,w=waterdrag,c=airdrag,f=coriolis
c--ancien (08/06/00) nom (file) = output(.f)
c  modif : 06/04/99
c- modif 24/10/19: cleaned use declaration sections              

c--definition

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip             
      use const_mod, only:
      use para0_mod, only:
      use para_mod, only:
      use bloc0_mod, only:
      use bloc_mod, only:
      use datadc_mod, only:
      use ice_mod, only: rhog
      use iceberg_mod, only: nftemp1,nftemp2,nftemp3,necriture,xn,yn,
     >    hiceb, wiceb, uiceb, viceb, lmx
      
!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

!PB#include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "iceberg.com"
! [SCRPTCOM] #include "ice.com"

!! END_OF_INCLUDE_SECTION

c--variables locales :
      character*8 fmtinf
      character*30 ficbtraj,fmtmas,fterre
      character*50 ccfile
      character*20 ftemp,fname
      integer(kind=ip) :: kiteration,kflag,kflag2,kflag3(5)

!PB variables added after imposing implicit none
      integer(kind=ip) :: nfmasse,nfterre,nfouticeb, nftermy,
     >                    i,nligne,l,iboucle
      real(kind=dblp) :: daysec,tmassetotale,tempo

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--instructions "data" :
! [SCRPTCOM] #include "datadc.com"

 1000 format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
 1111 format(3(F7.1,1X,F7.3,1X),I3,A)
 1222 format (4(F7.2,1x),1x,2(I7,1x))

c--constante
      daysec=86400
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c--variable en plus :
c      remsim=0.1494

C     write(nftemp1)
C     write(nftemp2)
C     write(nftemp3)
C     close(nftemp1)
C     close(nftemp2)
C     close(nftemp3)

c--definition des numeros de fichiers
      nfmasse=50
      nfterre=51
      nftemp1=52
      nftemp2=53
      nftemp3=54
C     nficeb=55
      nfouticeb=56
      nfmasse=57
      nfterre=58
C     nfterm=59
      nftermy=60

c--ouverture des fichiers
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C     open(nftemp1,file='tempxyhw.out',status='unknown',form=
C    &    'unformatted')
c     open(nftemp2,file='temptaille.out',status='unknown',form=
c    &    'unformatted')
C     open(nftemp3,file='tempterme.out',status='unknown',form=
C    &    'unformatted')



c---ecriture de la masse totale de glace d'iceberg sur un run
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C       open(nfmasse,file='icebmass.out',status='unknown')
C       fmtmas='(1PE14.6)'
C       write(nfmasse,1000) fmtmas, 999.,1, jm, 0,1
C       write(nfmasse,1111) 3., 0., 0., ddtt, 0., 0., 0
C       write(nfmasse,'(A,I4,A)') 'Masse totale de glace par iteration'
C       write(nfmasse,'(10000(A,I3))')
C    &        ('I', l,l=1,min(220,lmx))

c---ecriture du fichier de presence de terre
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C       open(nfterre,file='terre.out',status='unknown')
C       fterre='(1000(f5.2,1x))'
C
c-ecriture
C       do 50 j=jmax,1,-1
C         write (nfterre,fterre) (tms(i,j,ks2),i=1,120)
C         write (nfterre,'(1000I1)') (int(tms(i,j,ks2)),i=1,imax)
C  50   continue

c-fermeture
C       close(nfterre)

c----boucle principale
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c      write (*,*) 'nombre d iteration:',nit,' lmx',lmx
      do 2000 i=1,(necriture+1)
c     write(*,*) 'iteration enregistree: ',i
        nligne=0
C        read(nftemp1) nligne
C        read(nftemp1) (xn(l),l=1,nligne)
C        read(nftemp1) (yn(l),l=1,nligne)
C        read(nftemp1) (hiceb(l),l=1,nligne)
C        read(nftemp1) (wiceb(l),l=1,nligne)
C        read(nftemp1) (uiceb(l),l=1,nligne)
C        read(nftemp1) (viceb(l),l=1,nligne)

c--pour les terme il faut prendre des precautions car il y a une ligne de moins
c-(une iteration de moins que pour les xn,yn,etc.
c- MaB this if-loop is useless as everything was commented when i first saw it
c therefore i comment the if and endif - line as well (117, 131)
c        if (i.le.necriture) then
c          read(nftemp3) ndefault
c          read(nftemp3) (termexp(l),l=1,nligne)
c          read(nftemp3) (termexb(l),l=1,nligne)
c          read(nftemp3) (termexw(l),l=1,nligne)
c          read(nftemp3) (termexc(l),l=1,nligne)
c          read(nftemp3) (termexf(l),l=1,nligne)
c          read(nftemp3) (termeyp(l),l=1,nligne)
c          read(nftemp3) (termeyb(l),l=1,nligne)
c          read(nftemp3) (termeyw(l),l=1,nligne)
c          read(nftemp3) (termeyc(l),l=1,nligne)
c          read(nftemp3) (termeyf(l),l=1,nligne)
c          read(nftemp3) (termBIGx(l),l=1,nligne)
c          read(nftemp3) (termBIGy(l),l=1,nligne)
c        endif

c--correction des yn non egaux a 999
c--ces yn contiennent le num d'iteration de regeneration de l'iceb

      do 1400 l=1,nligne
         if (xn(l).eq.(999.)) then
           yn(l)=999.
           hiceb(l)=999.
           wiceb(l)=999.
           uiceb(l)=999.
           viceb(l)=999.
c          termexp(l)=999.
c          termexb(l)=999.
c          termexw(l)=999.
c          termexc(l)=999.
c          termexf(l)=999.
c          termeyp(l)=999.
c          termeyb(l)=999.
c          termeyw(l)=999.
c          termeyc(l)=999.
c          termeyf(l)=999.
c          termBIGx(l)=999.
c          termBIGy(l)=999.

          endif
1400  continue

c--completion des vides
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- MaB this if-loop is useless as everything was commented when i first saw it
ctherefore i comment the if, do, continue and endif - lines as well(163,165,185,186)
c      if ((lmx-nligne).ne.0) then
c        write (*,*) 'completion',nit,nligne
c        do 1500 j=1,(lmx-nligne)
C         yn(j+nligne)=999.
C         hiceb(j+nligne)=999.
C         wiceb(j+nligne)=999.
C         uiceb(j+nligne)=999.
C         viceb(j+nligne)=999.
c         termexp(j+nligne)=999.
c         termexb(j+nligne)=999.
c         termexw(j+nligne)=999.
c         termexc(j+nligne)=999.
c         termexf(j+nligne)=999.
c         termeyp(j+nligne)=999.
c         termeyb(j+nligne)=999.
c         termeyw(j+nligne)=999.
c         termeyc(j+nligne)=999.
c         termeyf(j+nligne)=999.
c         termBIGx(j+nligne)=999.
c         termBIGy(j+nligne)=999.

c          write (*,*) 'completion'
c1500    continue
c      endif



c--impression de mon fichier trajectoire pour l'antarctique
C     write(56,ficbtraj) (xn(l),yn(l),l,l=1,lmx)

c---impression de la masse totale
      tmassetotale=0.
      do 1501 iboucle=1,lmx
         if (xn(iboucle).ne.(999.)) then
c           if ((wiceb(iboucle).ge.(0.)).AND.(hiceb(iboucle).ge.(0.))
c     &          .AND.(rhog.ge.(0.))) then
           tempo=(1.5*wiceb(iboucle)*wiceb(iboucle)
     &                           *hiceb(iboucle)*rhog)

c           write(*,*) 'masse',tempo

           if (tempo.ge.(0.)) then
           tmassetotale=tmassetotale+((1.5*wiceb(iboucle)*
     &      wiceb(iboucle)*hiceb(iboucle)*rhog))

           endif
c          massetotale=massetotale+1
         endif
1501  continue
C     write(nfmasse,fmtmas) (tmassetotale*1e-9)
c      write(nfmasse,fmtmas) (nsr_frq(1))

c----terme des equations
c      write(nfterm,fmttrm)
c     &      (termexp(l),termexb(l),termexw(l),
c     &       termexc(l),termexf(l),termBIGx(l),l=1,lmx)

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c      write(nftermy,fmttrm)
c     &      (termexp(l),termexb(l),termexw(l),termexc(l),termexf(l),
c     &      termeyp(l),termeyb(l),termeyw(l),termeyc(l),termeyf(l),
c     &      l=1,lmx)

c      write (*,*) 'iteration:',i

2000  continue

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c-fermeture des fichiers

c--fermeture du fichier 'terme.out'
C     write(nfterm,*)
C     close(nfterm)

c--fermeture du fichier 'tempxyhw.out'
C     close(nftemp1)

c--fermeture du fichier 'temptaille.ou
C     close(nftemp2)

c--fermeture du fichier 'tempterme.out'
C     close(nftemp3)

c--fermeture du fichier 'termetotal.out'
C     write(nftermy,*)
C     close(nftermy)

c-fermeture trajectoires antarctique
      close(58)

c-fermeture ecriture masse icebger
C     write(nfmasse,*)
C     close(nfmasse)

c-------------------------------
      close(55)
c--fin

      return

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine iceb_out -
      end subroutine iceb_out
