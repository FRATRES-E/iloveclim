!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:56 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:56 CET 2009

      SUBROUTINE traitkd(qw,qlim,kmp,idmax,jdmax,kdmax,im,jm,
     &                   kinter,kcase,ijexcl,numf,titf,fmt)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Traitement de la bathymetrie : Amenagement des points enclaves en profondeur
!  si kmp > k(voisin)=Max{4kmu} > kinter : intervient.
! Selon kcase : 1 : demande le(s) niveau(x) de remplacement.
!               2 : approfondit automatiquement si QWmoyen > qlim.
!               3 : idem + suppression des enclaves restantes (si k > kinter)
!  modif : 31/10/96

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use para0_mod
      use para_mod
!! END_OF_USE_SECTION


!--dummy variables :
      
      integer(kind=ip)                             :: idmax,jdmax,kdmax
      
      integer(kind=ip), dimension(idmax, jdmax)    :: kmp
      integer(kind=ip), dimension(*)               :: ijexcl
      real(kind=dblp), dimension(kdmax,idmax,jdmax):: qw
      integer(kind=ip), dimension(4)               :: ksolu, lsolu
     &                             , isolu, jsolu
      real(kind=dblp), dimension(kmax,4)           :: qwmoy
      integer(kind=ip), dimension(imax,jmax)       :: kmu, lmu, kmv
     &                            , modifk

      

      character*(*) titf
      character*(*) fmt
!--local variables :

       integer(kind=ip) :: i, ii, ii1, ii2, im, j, jj, jj1, jj2, jm, kb
     &                   , kcase, kinter, kmm, l, ll, modif, n, nbsolu
     &                   , nexcl, nmodif, nn, npass, numf, iii, jjj, k
     &                   , kk1, kk2, knew, kq, kredon, ksol, nsol
     &                   , nsupp
     
       real(kind=dblp)  :: qlim

      character(len=20) :: fmtwr
      character(len=70) :: line

      numf = abs(numf)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  0 ) Initialisation (pour plusieurs passages) .
!-----------------------------------------------------------------------

!- initialisation :
      npass = 0
      nmodif = 0
      fmtwr = '(A7,I1,2(A,3I3))'
      do j=1,jmax
       do i=1,imax
         modifk(i,j) = 0
       enddo
      enddo

!--Reperage des zones a ne pas traiter :
      nexcl = 0
      do
        nn = 4 * nexcl + 1
        if (ijexcl(nn).ne.0) then
          nexcl = nexcl + 1
          cycle
        else
          exit
        endif
      enddo

!-----------------------------------------------------------------------
!--Ici debute le passage en revue de tous les points / modifk = 0
      loop_modif: do

      modif = 0
      npass = npass + 1

      write(6   ,'(1x,2A,I6,A,I3,A)') titf, 'Deep Enclosed pts :',
     &   nmodif, ' pts Deeper, Begin New checking (No=', npass, ' )'
      write(numf,'(1x,2A,I6,A,I3,A)') titf, 'Deep Enclosed pts :',
     &   nmodif, ' pts Deeper, Begin New checking (No=', npass, ' )'
!     write(numf,'(1x,3A,I3,A)') titf, 'Traite Enclave Profonde :',
!    &  ' Approfondissement (passage No :', npass, ' )'


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Detection des points enclaves .                                 |
!-----------------------------------------------------------------------

!--Calcul des Niv. Prof. Pts. Vitesse (grille B) :

!- initialisation :
      do j=1,jmax
       do i=1,imax
         kmu(i,j) = kdmax
         kmv(i,j) = kdmax
         lmu(i,j) = 0
         modifk(i,j) = modifk(i,j) + 1
       enddo
      enddo

!--Traitement des exclusions :
      do n=1,nexcl
        nn = 4 * (n-1)
        ii1 = ijexcl(nn + 1)
        ii2 = ijexcl(nn + 2)
        jj1 = ijexcl(nn + 3)
        jj2 = ijexcl(nn + 4)
        do j=jj1,jj2
         do i=ii1,ii2
           modifk(i,j) =  0
         enddo
       enddo
      enddo

      do j=2,jm
       do i=2,im
        do l=0,3
          ii = i - mod(l,2)
          jj = j - l/2
          if (kmp(ii,jj).lt.kmu(i,j)) then
            kmv(i,j) = kmu(i,j)
            kmu(i,j) = kmp(ii,jj)
            lmu(i,j) = l
          else
            kmv(i,j) = min(kmv(i,j),kmp(ii,jj))
          endif
        enddo
       enddo
      enddo

!----------------------------------------------------------
!- Boucle sur l'ensemble des points Interieurs au domaine -
!----------------------------------------------------------
      loop_390outer: do j=2,jm-1
       loop_390inner: do i=2,im-1

        kb = max( kmu(i,j), kmu(i+1,j), kmu(i,j+1), kmu(i+1,j+1) )
        if (kb.le.kinter) cycle loop_390inner
        if (kb.eq.kmp(i,j)) cycle loop_390inner
        if (modifk(i,j).ne.1) cycle loop_390inner

!--Detection des points enclaves - avec amelioration possible :
        nbsolu = 0
        kmm = kdmax
        do l=0,3
          ii = i + mod(l,2)
          jj = j + l/2
          if (kmu(ii,jj).gt.kinter.and.kmv(ii,jj).gt.kb) then
            ll = lmu(ii,jj)
            iii = ii - mod(ll,2)
            jjj = jj - ll/2
!- supression des redondances :
            kredon = 0
            do n=1,nbsolu
             if (iii.eq.isolu(n).and.jjj.eq.jsolu(n) ) kredon = 1
            enddo
            if (kredon.eq.0) then
              nbsolu = nbsolu + 1
              lsolu(nbsolu) = l
              kmm = min(kmm,kmu(ii,jj))
              isolu(nbsolu) = iii
              jsolu(nbsolu) = jjj
            endif
          endif
        enddo
        if (nbsolu.eq.0) cycle loop_390inner
        kmm = kmm + 1


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Traitement d'un point "enclave en profondeur" .                 |
!-----------------------------------------------------------------------

        if (kcase.eq.1) then
          write(6   ,'(2(I4,A),2I4)') nbsolu,
     &     ' Solutions pour desenclaver kmp=',kmp(i,j),' en (i,j)= ',i,j
          write(6   ,fmt) (qw(k,i,j),k=kmm,kmp(i,j))
        endif
        write(numf,'(2(I4,A),2I4)') nbsolu,
     &   ' Solutions pour desenclaver kmp=',kmp(i,j),' en (i,j)= ',i,j
        write(numf,fmt) (qw(k,i,j),k=kmm,kmp(i,j))

        do n=1,4
         do k=1,kmax
           qwmoy(k,n) = -1.
         enddo
        enddo



        do n=1,nbsolu

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Traitement de chaque solution possible :
        l = lsolu(n)
        ii = i + mod(l,2)
        jj = j + l/2
        iii = isolu(n)
        jjj = jsolu(n)

        kk1 = kmu(ii,jj) + 1
        kk2 = kmv(ii,jj)

        ksolu(n) = kmp(iii,jjj)
        do k=kk1,kk2
          qwmoy(k,n) = 0.5 * ( qw(k,iii,jjj) + qw(k,i,j) )
          if (qwmoy(k,n).ge.qlim) ksolu(n) = max(ksolu(n),k)
        enddo
        if (kcase.eq.1) then
          write(6   ,'(2I4,A,I4,A)')  iii, jjj,' =(i,j) ; kmp=',
     &                    kmp(iii,jjj), '  puis QW et QWmoy :'
          write(6   ,fmt) (qw(k,iii,jjj),k=kmm,kmp(i,j))
          write(6   ,fmt) (qwmoy(k,n),k=kmm,kmp(i,j))
        endif

        write(numf,'(2I4,A,I4,A)')  iii, jjj,' =(i,j) ; kmp=',
     &                   kmp(iii,jjj), '  puis QW et QWmoy :'
        write(numf,fmt) (qw(k,iii,jjj),k=kmm,kmp(i,j))
        write(numf,fmt) (qwmoy(k,n),k=kmm,kmp(i,j))

       enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Definition et Application de la solution .                      |
!-----------------------------------------------------------------------

!--Choix de la solution adoptee :

      if (kcase.eq.1) then
!--Choix sur mesure defini en interactif :
      write(66,*) 'Input Remplacement ? Centre / sol(n) (sur 1 ligne)'
     &         //'(-1 : aucune modif)'
      read(5,'(A)') line
      read(line,*) knew
      if (knew.ne.-1) then
        read(line,*) knew, (ksolu(n),n=1,nbsolu)

        if (kmp(i,j).ne.knew ) then
          modif =  1
          nmodif = nmodif + 1
          if (nmodif.eq.10)  fmtwr = '(A6,I2,2(A,3I3))'
          if (nmodif.eq.100) fmtwr = '(A5,I3,2(A,3I3))'
          write(6   ,fmtwr) 'AppEncl', nmodif, ' :',
     &                  i, j, kmp(i,j), ' -->', knew
          write(numf,fmtwr) 'AppEncl', nmodif, ' :',
     &                  i, j, kmp(i,j), ' -->', knew
          kmp(i,j) = knew
          do jj=j-1,j+1
           do ii=i-1,i+1
             modifk(ii,jj) = 0
           enddo
          enddo
        endif
      else
        nbsolu = 0
      endif

      else
!--Aprofondissement automatique :
      nsol = 1
      ksol = ksolu(1)
      do n=2,nbsolu
        if (ksolu(n).gt.ksol) then
          ksolu(nsol) = kmp(isolu(nsol),jsolu(nsol))
          nsol = n
          ksol = ksolu(n)
        elseif (ksolu(n).eq.ksol .and.
     &          qwmoy(ksolu(n),n).gt.qwmoy(ksol,nsol)) then
          ksolu(nsol) = kmp(isolu(nsol),jsolu(nsol))
          nsol = n
        else
          ksolu(n) = kmp(isolu(n),jsolu(n))
        endif
      enddo
       if (ksol.le.kb) then
!--Annulle la modif si elle n'entraine pas de reelle amelioration :
          ksolu(nsol) = kmp(isolu(nsol),jsolu(nsol))
       endif
!--fin traitement specifique selon kcase.
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Prise en compte de la modification :
      do n=1,nbsolu
        iii = isolu(n)
        jjj = jsolu(n)
        if (kmp(iii,jjj).ne.ksolu(n) ) then
          modif =  1
          nmodif = nmodif + 1
          if (nmodif.eq.10)  fmtwr = '(A6,I2,2(A,3I3))'
          if (nmodif.eq.100) fmtwr = '(A5,I3,2(A,3I3))'
!         write(6   ,fmtwr) 'AppEncl', nmodif, ' :',
!    &          iii, jjj, kmp(iii,jjj), ' -->', ksolu(n)
          write(numf,fmtwr) 'AppEncl', nmodif, ' :',
     &          iii, jjj, kmp(iii,jjj), ' -->', ksolu(n)
          kmp(iii,jjj) = ksolu(n)
          ii1 = max(iii-1,1)
          ii2 = min(iii+1,im)
          jj1 = max(jjj-1,1)
          jj2 = min(jjj+1,jm)
          do jj=jj1,jj2
           do ii=ii1,ii2
             modifk(ii,jj) = 0
           enddo
          enddo
        endif
        enddo
      write(numf,*)

!----------------------------------------------------------
!--Fin du Traitement du point kmp(i,j)=1 .
        enddo loop_390inner
       enddo loop_390outer


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--Un autre passage :
      if (modif < 1) exit loop_modif

      enddo loop_modif

      write(numf,*) 'Total ', titf, ' :',nmodif,' Approfondissements.'
      write(numf,*)

      if (kcase.lt.3) return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Suppression des Points Enclaves .
!----------------------------------------------------------

      write(6   ,*) titf, 'Traite Enclave Profonde : Suppression'
      write(numf,*) titf, 'Traite Enclave Profonde : Suppression'

!--Met en place d'un nouveau tableau "modifk" :
      do j=1,jmax
       do i=1,imax
         modifk(i,j) = 1
       enddo
      enddo

!- Traitement des exclusions :
      do n=1,nexcl
        nn = 4 * (n-1)
        ii1 = ijexcl(nn + 1)
        ii2 = ijexcl(nn + 2)
        jj1 = ijexcl(nn + 3)
        jj2 = ijexcl(nn + 4)
        do j=jj1,jj2
         do i=ii1,ii2
           modifk(i,j) =  0
         enddo
        enddo
      enddo

      nsupp = 0
      fmtwr = '(A7,I1,2(A,3I3))'
      do j=2,jm-1
       do i=2,im-1
!- if ( kmp > kb > kinter ) => suppress
        kb = max( kmu(i,j), kmu(i+1,j), kmu(i,j+1), kmu(i+1,j+1) )
        if (kmp(i,j).gt.kb.and.kb.gt.kinter.and.modifk(i,j).eq.1) then
          nsupp = nsupp + 1
          if (nsupp.eq.10)  fmtwr = '(A6,I2,2(A,3I3))'
          if (nsupp.eq.100) fmtwr = '(A5,I3,2(A,3I3))'
          kq = min(kmp(i,j)+1,kmax)
!         write(6   ,fmtwr) 'SupEncl', nsupp, ' :',
!    &                  i, j, kmp(i,j), ' -->', kb
          write(numf,fmtwr) 'SupEncl', nsupp, ' :',
     &                  i, j, kmp(i,j), ' -->', kb
          write(numf,fmt) (qw(k,i,j),k=1,kq)
          kmp(i,j) = kb
        endif
       enddo
      enddo
      write(6   ,*) 'Total ', titf, ' :',nsupp,' Suppressions.'
      write(numf,*) 'Total ', titf, ' :',nsupp,' Suppressions.'
      write(numf,*)

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine traitkd -
      end
