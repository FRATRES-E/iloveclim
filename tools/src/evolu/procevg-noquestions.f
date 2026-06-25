      PROGRAM PROCEV
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c   lit un fichier de type "evolu" et le modifie : ajout/supress colonnes ;
c ou Moyenne/Min/Max sur plusieurs lignes.
c  modif : 22/11/99
 
c--declaration implicite de type (standard fortran) :
      implicit double precision (a-h,o-z)
 
      parameter (imax = 200, jlismx = 99)
      parameter (nchmx1 = 120 , nchmx2 = 1000)
 
c--variables locales
      logical bflag
      dimension idel(imax), iadd(0:imax), ipnt(imax), ikeep(imax)
      dimension var(imax), varnit(imax)
      dimension varinp(imax,0:jlismx-1), wliss(jlismx)
 
      character*16 cc16
      character*20 titv(imax), tit00
      character*30 fmt, fmtfc, fmtrw
      character*40 fichr, fichw
      character*80 line
      character*(nchmx1) titr1, titr2
      character*(nchmx2) titrb
 
 1000 format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
 1111 format(3(F7.1,1X,F7.3,1X),I3,A)
 
      epsil = 1.d-10
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  1 ) Initialisation - Noms des fichiers I/O  -  Ouverture et Lecture |
c-----------------------------------------------------------------------
 

cdmr --- version sans questions : lit evolu-orig ecrit evolu
c      write(6,*) 'Nom du fichier a traiter (max 40 cc) ?'
c      read(5,'(A)') fichr
c      write(6,*) 'Nom du fichier resultat  (max 40 cc) ?'
c      read(5,'(A)') fichw
 
      fichr='evolu-orig'
      fichw='evolu'


      nncr = 0
      nncw = 0
      do 110 n=1,len(fichr)
        if (fichr(n:n).ne.' ') nncr = n
        if (fichw(n:n).ne.' ') nncw = n
 110  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Ouverture, Lecture de l'en-tete et des titres du fichier a traiter :
 
      open(30, file = fichr(:nncr), status='old')
      read(30,'(2A)') fmtfc , line
      read(line,*) spv, im0, jm, km, nfrc
      read(30,*) xi1, dxi, yj1, dyj, zk1, dzk, ntyp
      read(30,'(A)') titr1
      read(30,'(A)') titrb
 
      if (im0.gt.imax) then
        write(6,*) 'Pb. Sous Dim. im,imax=', im0, imax
        close(30)
        stop
      endif
 
      nct1 = 0
      do 120 n=1,len(titr1)
        if (titr1(n:n).ne.' ') nct1 = n
 120  continue
      write(6,'(2A)') titr1(:nct1)
 
c- extraction du format utilisee pour R/W :
      fmtrw = fmtfc
      bflag = 1
      do 130 n=2,len(fmtrw)
        if (bflag .and. fmtrw(n:n).eq.'X') then
          fmtrw(n-1:n) = 'A'//fmtrw(n-1:n-1)
          bflag = 0
        endif
 130  continue
      write(6,'(2A)') ' fmtrw = ', fmtrw
 
      ncc = nint(xi1)
      ncctb = ncc * im0
      ncctb = min( ncctb, nchmx2 )
 
      im1 = im0
      do 150 i=1,im0
        ikeep(i) = 1
        ipnt(i) = i
        titv(i)(:ncc) = titrb(1+(i-1)*ncc:i*ncc)
 150  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2 ) Definition de la modification :
c----------------------------------------------------------------------
 
c--Impression a l'ecran (pour faciliter le choix) :
      write(6,'(A,I4,A,I4,A)') 'Total : ',im0,' series , de ',
     &      jm,' (= jm) enregistrements'
      if (ncc.ne.0) then
        fmt = '(4(A,I3,1X,A))'
        write(6,fmt) (' No = ',i,titv(i)(:ncc),i=1,im0)
        write(6,*)
      endif
 
cmdr --- version sans questions : pas de modifications de colonnes
c      write(6,'(A)') ' Input Nb de colonne a Ajouter (+Nb) '
c     &             //'/ a Suprimer (-Nb) ? (0 => autre)'
c      read(5,*) modcol
      modcol=0
    
c-----
      liscol = 0
      if (modcol.eq.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Lissage par moyenne glissante :
cdmr --- Version sans questions : pas de moyenne glissante !
!        write(6,'(2A)') ' Lissage par moyenne glissante :',
!     &                  ' Input Nb de colonnes concernnees ?'
!        read(5,*) liscol
        liscol = 0  
cdmr --- 
        liscol = max(0,liscol)
c------
       if (liscol.ge.1) then
        write(6,'(2A)') ' No 1ere colonne a Lisser : lis1c= ?'
        read(5,*) lis1c
        lis2c = lis1c+liscol-1
        if (lis1c.lt.1 .or. lis2c.gt.im0) then
          write(6,*) 'Not valid input : lis1c,lis2c,im0=',
     &               lis1c,lis2c,im0
          liscol = 0
        endif
       endif
c------
       if (liscol.ge.1) then
        write(6,'(3(A,I4))') ' Lissage des colonnes',
     &     lis1c, ' -->', lis2c, ' , par moyenne glissante'
        write(6,'(2A)') '  sur jliss valeurs (jliss impair) :',
     &                  ' Input jliss = ?'
        read(5,*) jliss
        jliss = min(jlismx,jliss)
        jliss = (jliss - 1) / 2
        jliss = max(0,jliss)
        jliss = 1 + 2*jliss
        write(6,'(A,I3,A)') ' Poids des jliss(=', jliss, ') valeurs ?'
        read(5,*) (wliss(j),j=1,jliss)
       endif
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
      inter1 = 0
      if (modcol.eq.0 .and. liscol.eq.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Condense le fichier par Moy/Min/Max, par groupe de "inter1/inter0" lignes
        inter0 = nint(1000.*dyj)
        write(6,'(2A)') 'Regroupement de lignes :',
     &   ' Applique Moyenne ou Min/Max a chaque groupe de lignes.'
cdmr --- Version sans questions
c        write(6,'(2A,I5,A)')' Nb_d''iter par groupe de lignes ?',
c     &   ' (Nb_actuel=', inter0, ') ; Input inter1= ?'
c        read(5,*) inter1
        inter1=360
c        write(6,'(2A)') ' Nb_Variables a traiter par Min/Max ?',
c     &   ' (defaut=Moy.) ; Input nbmnx = ?'
c        read(5,*) nbmnx
        nbmnx = 0
cdmr ---
    
        nbmnx = max(0,nbmnx)
        if (inter1.eq.0) then
          close(30)
          stop
        elseif (inter1.lt.inter0) then
          write(6,'(3(A,I6))') 'STOP <== Bloc trop petit : inter1=',
     &        inter1, ' <', inter0, ' = Initial !'
          close(30)
          stop
        elseif (im0+nbmnx.gt.imax) then
          write(6,*) 'Pb. Sous Dim. im,imax=', im0+modcol, imax
          close(30)
          stop
        elseif (nbmnx.gt.0) then
          write(6,'(2A)') ' Input liste (1 line) des No_Var.',
     &      '(<->input file) a traiter par Min/Max :'
          read(5,*) (iadd(n),n=1,nbmnx)
          kok = 1
          do 205 n=2,nbmnx
            if (iadd(n).le.iadd(n-1)) kok = 0
 205      continue
          if (iadd(1).le.1 .or. iadd(nbmnx).gt.im0 .or. kok.eq.0) then
            write(6,'(A,(20I4))') 'STOP <= Pb in Liste No_Var :',
     &                    (iadd(n),n=1,nbmnx)
            close(30)
            stop
          endif
        endif
c-----
      elseif (modcol.gt.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Ajout de colonnes :
        if (im0+modcol.gt.imax) then
          write(6,*) 'Pb. Sous Dim. im,imax=', im0+modcol, imax
          close(30)
          stop
        endif
        write(6,'(A)') 'Ajout de Colonnes (Dans l''ordre !) :'
     &               //' Num_Var = celui du fichier final'
        write(6,'(A,I2,A)') ' Input (for each Var+) : tit_Var (',
     &                      ncc, ' cc), Num_Var, Value(Var) ?'
        do 210 n=1,modcol
          read(5,'(2A)') titv(im0+n)(:ncc), line
          read(line,*) iadd(n), var(im0+n)
 210    continue
c-----
      elseif (modcol.lt.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c--Supression de colonnes :
        write(6,*) 'Input : Numero des colonnes a suprimer ?'
        read(5,*) (idel(n),n=1,-modcol)
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  2.b) Prepare la modif, titres et autres :
c-----------------------------------------------------------------------
 
      if (modcol.gt.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        iadd(0) = 0
        do 230 n=1,modcol
          ii = im0+n
          if (iadd(n).le.im1+1 .and. iadd(n).gt.iadd(n-1)) then
            do 220 i=im1,iadd(n),-1
              ipnt(i+1) = ipnt(i)
 220        continue
            im1 = im1 + 1
            ipnt(iadd(n)) = ii
            write(6,'(3A,I4,A,1PE13.5)') ' +Var : ', titv(ii)(:ncc),
     &                    ', No=', iadd(n), ', Val=', var(ii)
          else
            write(6,'(3A,I4,A)') '=>Var : ', titv(ii)(:ncc),
     &                    ', No=', iadd(n), ', Cannot be added !'
            iadd(n) = iadd(n-1)
          endif
 230    continue
c- prepare 2nd ligne de titre et decale de 2cc chaque titre_Var
        n2cc = min(ncc+2,20)
        n0cc = n2cc - 2
        do 240 i=1,im1
          titrb(1+(i-1)*ncc:i*ncc) = titv(ipnt(i))(:ncc)
          titv(ipnt(i))(:n2cc) = '  '//titv(ipnt(i))(:n0cc)
 240    continue
        ncctb = ncc * im1
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
      if (modcol.lt.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        do 260 n=1,-modcol
          if (idel(n).ge.1 .and. idel(n).le.im0 .and.
     &        ikeep(idel(n)).eq.1) then
            im1 = im1 - 1
            ikeep(idel(n)) = 0
            write(6,'(A,I4,2X,A)') ' remove Var = ', idel(n),
     &                               titv(idel(n))(:ncc)
          else
            write(6,'(A,I4)') ' Error : Cannot remove Var.No=',
     &                          idel(n)
          endif
 260    continue
c- prepare 2nd ligne de titre (supression des titres correspondant) :
        ii = 0
        do 270 i=1,im0
          if (ikeep(i).eq.1) then
            ii = ii + 1
            ipnt(ii) = i
            titrb(1+(ii-1)*ncc:ii*ncc) = titv(i)(:ncc)
          endif
 270    continue
        ncctb = ncc * im1
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
      if (inter1.ne.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        iadd(nbmnx+1) = 0
        n = 0
        nn1 = n + 1
        do 280 i=1,im0
          if (i.eq.iadd(nn1)) then
            ipnt(i+n) = i
            n = nn1
            nn1 = n + 1
            ipnt(i+n) = im0 + n
c- Modif titre [call chguplw, Conversion -> Lower_case / Upper_case] :
            call chguplw(titv(i)(:ncc), titv(i),     1)
            call chguplw(titv(i)(:ncc), titv(im0+n), 2)
            write(6,'(A,I4,3A,2I4,A)')  ' Var i=', i,
     &       ' => Min/Max : ', titv(i)(:ncc)//' '//titv(im0+n)(:ncc),
     &       ' (no=', i+n-1, i+n, ')'
          else
            ipnt(i+n) = i
          endif
 280    continue
        im1 = im1 + n
C       write(6,'(A,3I4,A)') ' im0, nbmnx, im1 =', im0, nbmnx, im1,
C    &         ' ; ipnt(i) :'
C       write(6,'(20I4)') (ipnt(i),i=1,im1)
C       write(6,*)
c- prepare 2nd ligne de titre et decale de 2cc chaque titre_Var
        n2cc = min(ncc+2,20)
        n0cc = n2cc - 2
        titv(1)(:n2cc) = titv(1)(:n0cc)//'  '
        do 290 i=2,im1
          titrb(1+(i-1)*ncc:i*ncc) = titv(ipnt(i))(:ncc)
          titv(ipnt(i))(:n2cc) = '  '//titv(ipnt(i))(:n0cc)
 290    continue
        ncctb = ncc * im1
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3.a ) Lecture et reecriture sur fichier : Cas Supress/Ajout Sans Modif.
c-----------------------------------------------------------------------
 
      write(6,'(2A)') 'Debut de Lecture / Reecriture sur ',
     &                fichw(:nncw)
      open(36, file = fichw(:nncw), status='unknown')
 
      if (modcol.ne.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Conserve le meme Nb de lignes , change les titres des var. :
 
      write(36,1000) fmtfc, spv, im1, jm, km, nfrc
      write(36,1111) xi1, dxi, yj1, dyj, zk1, dzk, ntyp
      write(36,'(A)') titr1(:nct1)
      write(36,'(A)') titrb(:ncctb)
 
c--Supression/Ajout Sans Modif :
      do 310 j=1,jm
        read(30,fmtrw,end=500)  (titv(i),var(i),i=1,im0)
        write(36,fmtrw) (titv(ipnt(i)),var(ipnt(i)),i=1,im1)
 310  continue
 
      elseif (liscol.ne.0) then
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  3.b ) Lecture et reecriture sur fichier : Cas avec Lissage :
c-----------------------------------------------------------------------
 
c- Conserve la meme entete, les memes titres, le meme Nb de lignes :
      write(cc16,'(2(A,I4))') ' Lis.C', lis1c, ' -', lis2c
      write(36,1000) fmtfc, spv, im1, jm, km, nfrc
      write(36,1111) xi1, dxi, yj1, dyj, zk1, dzk, ntyp, cc16
      write(36,'(A)') titr1(:nct1)
      write(36,'(A)') titrb(:ncctb)
 
      jlisd2 = (jliss - 1) / 2
      do 390 j=1,jm+jlisd2
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        jj = mod(j,jliss)
c- jcw = indice du centre de l'interval [j+1-jliss,j] = [jcw-jlisd2,jcw+jlisd2]
        jcw = j - jlisd2
        jjw = mod(jcw,jliss)
c- lecture :
        if (j.le.jm) read(30,fmtrw,end=500)
     &                   (titv(i),varinp(i,jj),i=1,im0)
        if (j.gt.jlisd2) then
c- Sans lissage : remplissage de var() par le record au centre de l'interval :
          do 350 i=1,im0
            var(i) = varinp(i,jjw)
 350      continue
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- Avec lissage, Somme de jl1 a jl2 :
          jl1 = max(1,jcw-jlisd2)
          jl2 = min(jm,j)
          do 370 i=lis1c,lis2c
            wwsum = 0.
            vvsum = 0.
            do 360 jl=jl1,jl2
              jjl = mod(jl,jliss)
              jwl = 1 + jl - jcw + jlisd2
              vvv = varinp(i,jjl)
              if (vvv.ne.spv) then
                wwsum = wwsum + wliss(jwl)
                vvsum = vvsum + wliss(jwl)*vvv
              endif
 360        continue
c- moyenne glogale conservee (serie etandue sur les bords a l'origine spv) :
            if (wwsum.gt.epsil) var(i) = vvsum / wwsum
c- conserve les memes valeurs speciales :
C           if (var(i).ne.spv)  var(i) = vvsum / wwsum
 370      continue
c- ecriture de var() :
          write(36,fmtrw) (titv(i),var(i),i=1,im0)
        endif
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 390  continue
 
c-----
 
      else
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  4 ) Condense (Moy/min/Max) par groupe de "inter1/inter0" lignes :
c-----------------------------------------------------------------------
 
c--Lit la 1ere ligne ; Calcule le Nb de lignes en sortie :
      read(30,fmtrw,end=500)  (tit00,var(i),i=1,im0)
 
c- elimine Pb troncature du au format :
      xxit1 = var(1)/DFLOAT(inter0)
      itr1 = nint(xxit1) * inter0
      itr0 = itr1 - inter0
      itend = itr0 + inter0 * jm
 
      idphas = mod(itr0+inter1,inter1)
      if (idphas.ne.0) then
        write(6,'(2A,I4,A)')' Synchoniser le debut de fichier :',
     &   ' Yes=1, No=0 ? (idphas=', idphas, ' )'
        read(5,*) kdphas
        if (kdphas.eq.0) idphas = 0
      endif
 
      midj = (1 + inter0) / 2
      jj1   = (inter1 + itr1  - midj - idphas) / inter1
      jjend = (inter1 + itend - midj - idphas) / inter1
      jm1 = jjend - jj1 + 1
      yj11 = (idphas + jj1 * inter1) / 1000.
      dyj1 = inter1 / 1000.
 
      xjm = DFLOAT(jm * inter0) / DFLOAT(inter1)
      write(6,'(A,2I8,2I6,F8.3,F12.3)')
     & ' itr1,itend, j1,jend, yj1,xjm=', itr1,itend,jj1,jjend,yj11,xjm
 
c--Ecriture de l'en-tete et titre du fichier de sortie :
      write(36,1000) fmtfc, spv, im1, jm1, km, nfrc
      write(36,1111) xi1, dxi, yj11, dyj1, zk1, dzk, ntyp
      write(36,'(A)') titr1(:nct1)
      write(36,'(A)') titrb(:ncctb)
 
c--Initialisation :
      itr = itr1
      jjp = jj1
      jsum = 1
      do 405 n=1,nbmnx
        var(n+im0) = var(iadd(n))
 405  continue
      do 410 i=1,im1
        varnit(i) = var(ipnt(i))
 410  continue
 
      do 480 j=2,jm
        read(30,fmtrw,end=500)  (tit00,var(i),i=1,im0)
        itr = itr + inter0
        jjn = (inter1 + itr - midj - idphas) / inter1
c-----
        if (jjn.eq.jjp) then
c--Reste dans le meme groupe de lignes :
        n = 1
        ii = 1
        varnit(1) = var(1)
        do 430 i=2,im0
          ii = ii + 1
          if (i.eq.iadd(n)) then
            varnit(ii) = min(varnit(ii),var(i))
            ii = ii + 1
            varnit(ii) = max(varnit(ii),var(i))
            n = n + 1
          else
            varnit(ii) = varnit(ii) + var(i)
          endif
 430    continue
        jsum = jsum + 1
c-----
        endif
 
c--Changement de groupe de lignes.
        if (jjn.gt.jjp.or.j.eq.jm) then
c-----
c- calcul des moyennes & ecriture sur fichier, ligne(jjp) :
        zsum = 0.
        if (jsum.ne.0) zsum = 1.0 / DFLOAT(jsum)
        n = 1
        do 440 i=2,im0
          if (i.ne.iadd(n)) then
            ii = i + n - 1
            varnit(ii) = varnit(ii) * zsum
          else
            n = n + 1
          endif
 440    continue
        write(titv(1),'(A,I2)') titv(1)(:n0cc), jsum
        write(36,fmtrw) (titv(ipnt(i)),varnit(i),i=1,im1)
c-----
        if (jjp.eq.jj1) write(6,'(A,4I8)')
     &      '1rst line : jj,jsum,j,itr :', jjp,jsum,j-1,itr-inter0
        if (j.eq.jm) write(6,'(A,4I8)')
     &      'last line : jj,jsum,j,itr :', jjp,jsum,j-1,itr-inter0
c-----
c- debute le groupe de lignes suivant :
        jsum = 1
        do 445 n=1,nbmnx
          var(n+im0) = var(iadd(n))
 445    continue
        do 450 i=1,im1
          varnit(i) = var(ipnt(i))
 450    continue
c-----
c- Ecriture d'une ligne isolee (=la derniere du fichier initial) :
        if (jjn.gt.jjp.and.j.eq.jm) then
          write(titv(1),'(A,I2)') titv(1)(:n0cc), jsum
          write(36,fmtrw) (titv(ipnt(i)),varnit(i),i=1,im1)
          write(6,'(A,4I8)')
     &      'Last line : jj,jsum,j,itr :', jjn,jsum,j,itr
        endif
        jjp = jjn
c-----
        endif
 480  continue
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  5 ) Fin de fichier ; Fermeture & Arret ; Cas Fichier incomplet .
c-----------------------------------------------------------------------
 
c- termine ecriture & ferme les fichiers.
      read(30,*)
      write(36,*)
 
      close(30)
      close(36)
      stop
 
c--Traitement des cas d'erreur de lecture :
 500  continue
      write(6,*) 'Probleme de LECTURE, fichier = ', fichr(:nncr)
      stop
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine PROCEV -
      end
      subroutine chguplw(linein, linout, k)
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c  exchange Upper_case <-> Lower_case character  (selon "k")
c Input = linein ; Output = linout   (called with same var. in/out => it works)
c-    (Majuscule: 65 -> 90 ; minuscule: 97 -> 122)
c  modif : 22/02/99
 
c- dummy variables :
      character*(*) linein, linout
 
c- local variables :
      dimension icor(0:7)
      data icor / 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 /
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c- Definit la transformation :
      if (k.eq.0) then
c-- Lower_case <--> Upper_case :
        icor(2) = 3
        icor(3) = 2
      elseif(k.eq.1) then
c-- All char --> lower_case :
        icor(2) = 3
        icor(3) = 3
      elseif(k.eq.2) then
c-- All char --> Upper_case :
        icor(2) = 2
        icor(3) = 2
      else
        icor(2) = 2
        icor(3) = 3
      endif
 
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
 
c- Applique la transformation :
      nc = len(linein)
      do 200 n=1,nc
        i = ichar(linein(n:n))
        ic = i / 32
        ir = mod(i,32)
        if (ir.ge.1 .and. ir.le.26) ic = icor(ic)
        j = 32*ic + ir
        linout(n:n) = char(j)
 200  continue
 
      return
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
c- fin de la routine chguplw -
      end
