C     subroutine rdcour(var,spval,xi1,dxi,nszmax,nbcmax,im,jm,nturn,
C    &                  nccv,numf,fildat,tit1,titvar)
      program evo2gra
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c   Read data (usual format )(jm > 0 => j increasing) file = fildat .
c   work on several variables (im) , 1 per column
c  modif : 24/05/95
 
      parameter (nbcmax = 130)
      parameter (nchmax = 90)
cdmr ---
      parameter (maxime = 15000, maximien = 120)
      parameter (nszmax = maxime*maximien)
c--dummy variables :
      real var(nszmax)
      real*4 sort(maximien,maxime)
      character*5 fildat
      character*80 tit1
      character*40 titvar(nbcmax)
 
c--variables locales
      character*2  cc2
      character*10 cc10
      character*30 fmt , fmtitr , fmtw
      character*40 ccline
      character*70 line
      character*(nchmax) tit2
 
 1000 format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
 1111 format(3(F7.1,1X,F7.3,1X),I3,A)
 
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  1 ) Ouverture et lecture de l'entete du fichier :                   |
c-----------------------------------------------------------------------
 
      numf=5
C     write(6,*) len(tit1), len(titvar(1))
cdmr --- Version sans questions ...

c      write(6,*) 'Name of the data file ?'
c      call rdstdin(numf,fildat)
 
      fildat='evolu'
cdmr ---

      open(unit=10, file = fildat, status='old')
 
      ltest = 2
      kprec = 0
      km = 0
      jprec = 0
 110  continue
      if(km.eq.0) then
        read(10,'(2A)') fmt, line
        read(line,*)  ovsc, im, jm, km, nfrc
        read(10,*) xi1, dxi, yj1, dyj, zk1, dzk, ntyp
        jmp = abs(jm)
      endif
 
c entree au clavier du numero du tableau (ou bloc) a traiter
      if (km.eq.1) then
        ks = 0
      else
 120    continue
cdmr --- Versions sans questions ...
c        write(6,'(A)') 'Number of the arrays to vizualise ? '
c     &     // '(-num < 0 for prolongation ; num > 0 if last)'
c        call rdstdin(numf,ccline)
c        read(ccline,*) ks
        ks = 1
cdmr ---
        if(km.ne.0.and.iabs(ks).gt.iabs(km)) then
          write(6,*) 
     &    'Greater than the number of arrays in the file = ',km
          goto 120
        endif
      endif
      if(ks.eq.0) then
        write(6,*) 'num = 0 ! par defaut lecture du tableau qui suit'
        ks = 1
        ltest = 1
      elseif (ks.gt.0) then
        ks = max0(1,ks - kprec)
        ltest = 0
      else
        ks = -ks
        ks = max0(1,ks - kprec)
        ltest = 2
      endif
 
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  2 ) lecture standard de plusieurs tableaux [ou blocs si km=0]       |
c-----------------------------------------------------------------------
 
c--positionnement au debut du tableau (km # 0) a traiter :
      if  (km.ne.0.and.ks.ge.2) then
        if (nfrc.ge.1) then
          njj = (im-1)/nfrc + 1
          njj = (ks - 1) * (njj*jmp + 3)
        else
          njj = (im-1)/(-nfrc) + 1
          njj = (ks - 1) * (njj*(jmp+1) + 3)
        endif
        do 210 j=1,njj
           read(10,*)
 210    continue
      endif
 
c--positionnement au debut du bloc (km = 0) a traiter :
      if (km.eq.0.and.ks.ge.2) then
        ktot = 1
        do 230 k=1,ks-1
c- lecture d'un bloc :
          if (nfrc.ge.1) then
            njj = (im-1)/nfrc + 1
            njj = ktot * (njj*jmp + 3)
          else
            njj = (im-1)/(-nfrc) + 1
            njj = ktot * (njj*(jmp+1) + 3)
          endif
          do 220 j=1,njj
            read(10,*)
 220      continue
          read(10,'(2A)') fmt, line
          read(line,*)  ovsc, im, jm, ktot, nfrc
          read(10,*) xi1, dxi, yj1, dyj, zk1, dzk, ntyp
          ktot = max0(1,iabs(ktot))
          jmp = abs(jm)
 230    continue
      endif
 
      if(jprec.eq.0) spval = ovsc
 
c--verification des dimensions :
      nsz = (jmp+jprec) * im
      if(nsz.gt.nszmax) then
        write(6,*) 'Arret dans rdcour :'
     a     //' Pb de sous-dimensionnement des tableaux !'
        stop
      endif
 
c--lecture du 1er titre du tableau :
      cc10 = tit1(47:56)
      read(10,'(A)') tit1
      if (km.ne.0) then
        write(6,*) 'Titre associe au tableau traite (1 lignes) :'
        write(6,'(A79)') tit1
      else
        write(6,*) 'Titre du 1er bloc (1 lignes) :'
        write(6,'(A79)') tit1
        if(len(tit1).ge.56.and.kprec.ne.0) tit1(47:56) = cc10
      endif
 
c--lecture du 2nd titre du tableau :
      do 250 i=1,nbcmax
        titvar(i) = ' '
 250  continue
      if(dxi.eq.0.) then
        nccv = nint(xi1)
        nbc = im
      else
        nccv = nint(yj1)
        nbc = jmp
      endif
      nbc = min0(nbc,nbcmax)
      if(nccv.ge.1.and.nccv.le.nchmax) then
c- Morcellement du 2nd titre :
cccccccccccccc
        write(cc2,'(I2)') nccv
        if(nccv.le.9) cc2 = cc2(2:)//' '
        write(fmtitr,'(A,I3,A)') '(',nbcmax,'A'//cc2//')'
        read(10,fmtitr) (titvar(i),i=1,nbc)
      else
        nccv = 0
        read(10,'(A)') titvar(1)
      endif
 
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c  3 ) lecture standard du tableau selectione :                        |
c-----------------------------------------------------------------------
 
c--lecture directe (j crois.) ou indirecte (j decrois.) selon le signe de jm :
      jj1 = max(1,-jm)
      jj2 = max(jm,1)
      jj3 = sign(1,jm)
      if (nfrc.ge.0) then
        do 310 j=jj1,jj2,jj3
          nj = (jprec + j -1) * im
          read(10,fmt) (var(i+nj),i=1,im)
 310    continue
        read(10,*)
      else
        do 320 iis = 1,im,-nfrc
          iie = min(im,iis-nfrc-1)
          do 315 j=jj1,jj2,jj3
            nj = (jprec + j -1) * im
            read(10,fmt) (var(i+nj),i=iis,iie)
 315      continue
          read(10,*)
 320    continue
      endif
 
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
 
c--Uniformisation de la valeur speciale :
      if(ovsc.ne.spval) then
        do 410 j=1,jmp
         nj = (jprec + j - 1) * im
         do 410 i=1,im
           if(var(i+nj).eq.ovsc) var(i+nj) = spval
 410    continue
      endif
 
c--Prolongation du set de donnees par celles d'un autre tableau :
      jprec = jprec + jmp
      kprec = kprec + ks
      if(ltest.eq.1) then
        write(6,*) 'Prolonge avec la lecture d un autre tableau ? '
     &           //'(0 non, > 0 oui) :'
        call rdstdin(numf,ccline)
        read(ccline,*) ltest
      endif
      if(ltest.ge.1) goto 110
 
      jm = jprec
      close(10)
 
c--mise en place de "nturn" :
      nturn = 0
      if (km.le.0) nturn = nint(dzk)
 
C     if (fact.ne.1.0) then
C       write(6,*) 'Facteur Multiplicatif ? '
C       read(6,'(E13.6)') fact
C       read(5,*) fact
C       write(6,*) 'FACTEUR MULTIPLICATIF : ',fact
C       write(6,*) 'Taper RETURN '
C       read(5,*)
C     endif
 
C     do 450 n=1,nsz
C       yy = var(n)
C       if(yy.eq.ovsc) then
C         var(n) = 111111.
C       else
C         var(n) = yy * fact
C       endif
C450  continue
 
C ecriture du tableau :
C     do 500 j=jm,1,-1
C       nnj = (j-1)*im
C       write(11,fmtw) (var(i+nnj),i=1,im)
C500  continue
C     close(11)
       
      titvar(2)='Tyr  '
      titvar(27)='Tmc '
      titvar(28)='T1mo '
      titvar(29)='aTmo'
      titvar(30)='Sm30 '
      titvar(31)='S1mo '
      titvar(32)='aSmo'
      titvar(33)='a_w'
      titvar(34)='a_u'
      titvar(35)='a_v'
      titvar(37)='Tmo1'
      titvar(38)='Tmo2'
      titvar(39)='Tmo3'
      titvar(40)='Tmo4'
      titvar(41)='Tmo5'
      titvar(42)='Tmo6'
      titvar(43)='Tmo7'
      titvar(44)='Tmo8'
      titvar(45)='Tmo9'
      titvar(46)='Tmo10'
      titvar(47)='Tmo11'
      titvar(48)='Tmo12'
      titvar(49)='Tmo13'
      titvar(50)='Tmo14'
      titvar(51)='Tmo15'
      titvar(52)='Tmo16'
      titvar(53)='Tmo17'
      titvar(54)='Tmo18'
      titvar(55)='Tmo19'
      titvar(56)='Tmo20'
      titvar(57)='Smo1'
      titvar(58)='Smo2'
      titvar(59)='Smo3'
      titvar(60)='Smo4'
      titvar(61)='Smo5'
      titvar(62)='Smo6'
      titvar(63)='Smo7'
      titvar(64)='Smo8'
      titvar(65)='Smo9'
      titvar(66)='Smo10'
      titvar(67)='Smo11'
      titvar(68)='Smo12'
      titvar(69)='Smo13'
      titvar(70)='Smo14'
      titvar(71)='Smo15'
      titvar(72)='Smo16'
      titvar(73)='Smo17'
      titvar(74)='Smo18'
      titvar(75)='Smo19'
      titvar(76)='Smo20'
 
c--sorties
C     open(unit=21,file='outev.dat')
C     write(fmtw,'(A,I3,A7,A)') '(',jm,'(E13.5)',')'
C     do 600 i=1,im
C      do 600 j=1,jm
C        nnj = (j-1)*im
C       write(21,*) titvar(i)
C       write(21,fmtw) (var(i+(j-1)*im),j=1,jm)
C       write(21,*) var(i+(j-1)*im)
C600  continue
       do 610 i=1,im
          do 610 j=1,jm
             sort(i,j) = var(i+(j-1)*im)
 610   continue
C     write(6,*) '1',sort(1,1),sort(85,1)
C     write(6,*) '2',sort(1,2),sort(85,2)

C      jpxy=1*4
      jpxy=1*1

      PRINT*, "This is Didier M.'s output 0 : im,jm,jpxy",im,jm,jpxy

      open (8, FILE='gra'//fildat//'-mon.dat', FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=jpxy)
      write(6,*) 'open file OK'

      PRINT*, "This is Didier M.'s output  1 : im,jm,jpxy",im,jm,jpxy

      IREC=1
      do j=1,jm
      do i=1,im
C       write(8,REC=IREC) (sort(i,j),J=1,jm)
        write(8,REC=IREC) sort(i,j)
       IREC=IREC+1
      enddo
      enddo
      close(8)
      write(6,*) 'write Z OK'
      open (9, FILE='gra'//fildat//'-mon.ctl')
      write (9,'(A6,A16)') 'dset ^','gra'//fildat//'-mon.dat'
      write (9,'(A6,E12.4)') 'undef ',spval
cdmr      write (9,'(A18)') 'options big_endian'
      write (9,'(A)') 'title test of grads: evolu'
      write (9,'(A)') 'xdef 1 linear 1 1'
      write (9,'(A)') 'ydef 1 linear 1 1'
      write (9,'(A)') 'zdef  1 linear 5  5'
      write (9,'(A5,I5,A8,A9,A4)') 'tdef ',jm,' linear ',
     &                   '1jan0001 ','30dy'
      write (9,'(A5,I3)') 'vars ',im
      do i=1,im
        write (9,'(A7,A12)') titvar(i),'   1  9   99'
      enddo
      write (9,'(A)') 'endvars'
      close (9)

c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      end
      subroutine rdstdin(numf,line)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  ReaD STanDard INput : Lecture de la chaine de caractere "line"
c    sur l'unite "numf' si la fin du fichier n'est pas atteinte,
c    sur l'unite 5 (standard) si non.
c  modif :  20/07/93

      character*(*) line

      read(numf,'(A)',end=100) line
      write(6,'(A,I3,2A)') 'Unit', numf, ' Lu : ', line
      return
 100  continue
      numf = 5
      read(numf,'(A)') line
      write(6,'(A,I3,2A)') 'Unit', numf, ' Lu : ', line
      return
c---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
c- fin de la routine rdcour -
      end
