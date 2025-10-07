!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:56 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:56 CET 2009

      SUBROUTINE traitkm(qw,kmp,idmax,jdmax,kdmax,im,jm,numf,titf,fmt)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
! Traitement de la bathymetrie : Suppression des zones peu profondes
!    (kmp = 1 --> 0 ou --> 2 selon ses voisins).
! Supression des enclaves (= pts entourer de 3 pts terre) :
!  modif : 28/07/94

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip
      
      use newunit_clio_mod, only: clio3_out_id
      

!! END_OF_USE_SECTION


!+    include 'const.com'
!-    include 'para.com'

!--dummy variables :

      integer(kind=ip)                             :: idmax,jdmax,kdmax
      
      integer(kind=ip), dimension(idmax,jdmax)     :: kmp
      real(kind=dblp), dimension(kdmax,idmax,jdmax)::qw
    
      character*(*) titf
      character*(*) fmt

!--local variables :
      integer(kind=ip), dimension(4)               :: kordon
      integer(kind=ip) :: i, ie, ii, im, is, isav, j, je, jj, jm, js
     &                  , jsav, k, kk, kprtqw, kq, n, nbcas, nbprec
     &                  , npas, npermu, numf


      kprtqw = 0
      if (numf.lt.0) kprtqw = 1
      numf = abs(numf)

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  Remplacement des "kmp=1" , par plusieurs passages (npas) successifs |
!-----------------------------------------------------------------------

      npas  = 0
      nbcas = -1
      loop_replace_kmp1: do
      npas = npas + 1
      write(numf,*) titf, 'Suppress kmp=1, passage :', npas
      nbprec = nbcas
      nbcas = 0
      do j=1,jm
       do i=1,im
        if (kmp(i,j).eq.1) then

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!--Traitement d'un point de niveau 1 .

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  2 ) Classement, par ordre croissant, des 4 Nb de Niv. entourant i,j |
!-----------------------------------------------------------------------

!--initialisation .
      do n=1,4
        kordon(n) = kmp(i,j)
      enddo
      if (i.gt.1)  kordon(1) = kmp(i-1,j)
      if (j.gt.1)  kordon(2) = kmp(i,j-1)
      if (i.lt.im) kordon(3) = kmp(i+1,j)
      if (j.lt.jm) kordon(4) = kmp(i,j+1)

!--Classement par ordre croissant :
      do
        npermu = 0
        do n=1,3
          kk = kordon(n+1)
          if (kordon(n).gt.kk) then
            npermu = 1
            kordon(n+1) = kordon(n)
            kordon(n) = kk
          endif
        enddo
        if (npermu == 0) exit
      enddo

!--Fin du classement .

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  3 ) Remplacement d'un pt de niveau 1 selon les 4 Nb de Niv voisin . |
!-----------------------------------------------------------------------

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
       if (kordon(2).ge.2) then
!- au moins 3 voisins > 1 :  1 --> 2
        write(numf,'(2A,2I4)') titf, ' k=1 --> 2 en i,j :', i, j
        if (kprtqw.eq.1) write(numf,fmt) (qw(k,i,j),k=1,3)
        kmp(i,j) = 2
      elseif (kordon(3).ge.2) then
!- au moins 2 voisins > 1 :
        if (kordon(3).le.3) then
!- 2 ds [0,1], 1 3eme ds [2,3] : Zone peu profonde etandue (2pts) conservee :
          write(numf,'(2A,2I4)') titf, ' k=1 --> 2 en i,j :', i, j
          if (kprtqw.eq.1) write(numf,fmt) (qw(k,i,j),k=1,3)
          kmp(i,j) = 2
        elseif (kordon(2).eq.0) then
!- 2 pts = 0, 2 autres > 3 : Zone peu profonde isolee supprimee :
          write(numf,'(2A,2I4)') titf, ' k=1 --> 0 en i,j :', i, j
          if (kprtqw.eq.1) write(numf,fmt) (qw(k,i,j),k=1,3)
          kmp(i,j) = 0
        else
!- indecis :
          write(numf,'(2A,2I4)') titf, ' k=1 --> ? en i,j :', i, j
          nbcas = nbcas + 1
        endif
      elseif (kordon(3).eq.0) then
!- au moins 3 voisins = 0 :  1 --> 0
        write(numf,'(2A,2I4)') titf, ' k=1 --> 0 en i,j :', i, j
        if (kprtqw.eq.1) write(numf,fmt) (qw(k,i,j),k=1,3)
        kmp(i,j) = 0
      else
!- indecis :
        write(numf,'(2A,2I4)') titf, ' k=1 --> ? en i,j :', i, j
        nbcas = nbcas + 1
      endif
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
          isav = i
          jsav = j
!--Fin du Traitement du point kmp(i,j)=1 .

        endif
       enddo
      enddo
!--fin du passage "npas" .

      if (nbcas.eq.nbprec) then
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  4 ) Intervention dans les cas "delicat" (=traitement stationaire).  |
!-----------------------------------------------------------------------
      write(clio3_out_id,*) ' Intervention manuelle necesssaire'
     &         //' pour remplacer les pts a 1 niv. de Prof.'
      write(clio3_out_id,
     & '(2A,2I4)') titf, 'Localisation : (i,j) = ', isav, jsav
      is = max(1,isav-5)
      ie = min(im,isav+5)
      js = max(1,jsav-5)
      je = min(jm,jsav+5)
      write(clio3_out_id,'(12I3)') (mod(i,100),i=is,ie)
      do j=je,js,-1
        write(clio3_out_id,'(12I3)') (kmp(i,j),i=is,ie), j
      enddo
!     write(clio3_out_id,*) 'Intervention sur un voisin : Entrer i et j ?'
!     read(5,*) i, j
      i = isav
      j = jsav
      write(clio3_out_id,*) 'Profondeur de Remplacement ( 0 ou 2 ) ? '
      read(5,*) kk
      write(numf,'(A)') 'Intervention Manuelle :'
      write(numf,'(2A,I3,A,2I4)') titf, ' k=1 -->', kk,
     &          ' en i,j :', i, j
      if (kprtqw.eq.1) write(numf,fmt) (qw(k,i,j),k=1,3)
      kmp(i,j) = kk

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
      endif

      if (nbcas < 1) exit
      enddo loop_replace_kmp1

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  5 ) Souligne les Enclaves (= pts avec 3 voisins Terre ) :
!-----------------------------------------------------------------------

      write(numf,*)
      write(numf,'(A)') 'Supression des  pts Isoles :'

!--Supression des points Isoles :
      do j=2,jm-1
       do i=2,im-1
        if (kmp(i+1,j).ne.0.and.kmp(i-1,j).eq.0.and.kmp(i+1,j).eq.0
     &     .and.kmp(i,j-1).eq.0.and.kmp(i,j+1).eq.0) kmp(i,j) = 0
       enddo
      enddo

      write(6   ,*) 'Pointage des points Enclaves :'
      write(numf,*) 'Pointage des points Enclaves :'

      do j=1,jm
       do i=1,im
        kk = 6
        if (kmp(i,j).ne.0) kk = 4
!REFACTORING: Note that the following tests may lead to out-of-bounds addressing of kmp
!REFACTORING: if the two subconditions are evaluated in parallel (the first one does not
!REFACTORING: a priori prvent the second one from being carried out)
!REFACTORING: To be reformulated as shown below.
        if (i.gt.1 .and.kmp(i-1,j).eq.0) kk = kk - 1
        if (i.lt.im.and.kmp(i+1,j).eq.0) kk = kk - 1
        if (j.gt.1 .and.kmp(i,j-1).eq.0) kk = kk - 1
        if (j.lt.jm.and.kmp(i,j+1).eq.0) kk = kk - 1
!REFACTORING:        if (i > 1) then
!REFACTORING:          if (kmp(i-1,j) == 0) kk = kk - 1
!REFACTORING:        endif
!REFACTORING:        if (i < im) then
!REFACTORING:          if (kmp(i+1,j) == 0) kk = kk - 1
!REFACTORING:        endif
!REFACTORING:        if (j > 1) then
!REFACTORING:          if (kmp(i,j-1) == 0) kk = kk - 1
!REFACTORING:        endif
!REFACTORING:        if (j < jm) then
!REFACTORING:          if (kmp(i,j+1) == 0) kk = kk - 1
!REFACTORING:        endif


!--Cas d'un point Enclave :
        if (kk.eq.1) then

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
      write(6   ,'(A,I4,A,2I4)') 'Enclave kmp=', kmp(i,j),
     &                           ' en (i,j)= ', i, j
      write(numf,'(A,I4,A,2I4)') 'Enclave kmp=', kmp(i,j),
     &                           ' en (i,j)= ', i, j
      kq = min(1+kmp(i,j),kdmax)
      is = max(1,i-1)
      ie = min(im,i+1)
      js = max(1,j-1)
      je = min(jm,j+1)
      if (kprtqw.eq.1) then
        do jj=je,js,-1
         do ii=is,ie
          write(numf,'(2I4,A,I4,A)')  ii, jj,' =(i,j) ; kmp=',
     &                   kmp(ii,jj), '  et QW :'
          write(numf,fmt) (qw(k,ii,jj),k=1,kq)
         enddo
        enddo
        write(numf,*)
      endif
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
        endif
       enddo
      enddo

      return

!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!- fin de la routine traitkm -
      end
