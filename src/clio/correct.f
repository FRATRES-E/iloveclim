!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009

      SUBROUTINE correct(nn99,filcor)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  Appelee par "redforc" : Correction du forcage selon le fichier filcor.
!  Convention : ktyp < 3 remplace , 2 < ktyp < 6 additione , 5 < ktyp multiplie
!    mod(ktyp,3) = 0 spv sans effet, = 1 modif si spv, = 2 modif sauf si spv
!  modification locale ([is,ie]x[js,je] -> 2nd ligne) definie par valeur unique
!    ou globale par valeur correspondante lue sur fichier(-> 2nd ligne).
!--------
!  modif : 04/08/98

!! START_OF_USE_SECTION

      use const_mod, only:

      use para0_mod, only: imax, jmax, kmax, nrpmax, nsmax
      use para_mod,  only:
      use bloc0_mod, only: phisu, phisv, scalr, phiss, rappel, rappes, phimnx
     &                   , rapint
      use bloc_mod,  only: nrap, ijrap

      use global_constants_mod, only: dblp=>dp, ip
      use newunit_clio_mod, only: clio3_out_id, mouchard_id

!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

! [IMPLCTNONE] #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"

!! END_OF_INCLUDE_SECTION

!- dummy variables :
      character*(*) filcor

!- variables locales equivalentes :
      real(dblp), dimension(imax,jmax,2) :: phisuv
c~       equivalence ( phisuv(1,1,1) , phisu(1,1) )
c~       equivalence ( phisuv(1,1,2) , phisv(1,1) )

!- variables locales :
      real(dblp), dimension(imax,jmax) :: ycor

      character*40 filtab
      character*70 line

      real(dblp)  :: spv, spvcor
      integer(ip) :: nx, ns, ntabc, nvcor, nnadd, nn, nbcor, nt, nrapmx, ktyp
     &             , nc, kspv, kmodif, kcor1, kcor, kcor2, k, jcor1, jcor2, j
     &             , nn99, i, icor2, icor1

      integer(ip) :: filcor_id, filtab_id

 1200 format(A,3(2I4,A),1PE11.3)
 1300 format(A,2(2I4,A),I3,2I2,A,1PE11.3)
 1400 format(A,2(2I4,A),2I3,A,2I2,A,1PE11.3)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initialisation ; Ouverture du fichier "forc.corr" .
!-----------------------------------------------------------------------

      nrapmx = 0

      phisuv(:,:,1) = phisu(:,:)
      phisuv(:,:,2) = phisv(:,:)

      if (nn99.eq.2) write(mouchard_id,'(2A)')
     &  'Forcage Modifie par le fichier ', filcor

      open(newunit = filcor_id, file = filcor, status = 'old')
      read(filcor_id,*)
      read(filcor_id,*)
      read(filcor_id,*)
      read(filcor_id,*)
      read(filcor_id,*)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      read(filcor_id,*) ntabc
      do nt=1,ntabc
!--Definition du type de modification :
!  Nbcor = Nb de lignes (a lire) pour definir les modifs ;
!  nvcor = tableau modifie ; kcor = 3eme indice (ex:niveau k) du Tab modifie
!  ktyp : cf en-tete ; spv = Special-Value a prendre en compte.
!-----
      read(filcor_id,'(A)') line
      if (nn99.eq.2) write(mouchard_id,'(A)') line
      read(filcor_id,*) nbcor, nvcor, kcor, ktyp, spv
      kspv = mod(ktyp,3)
      kmodif = ktyp / 3

      do nc=1,max(nbcor,min(1,-nbcor))
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Modif par portion ou en totalite d'un tableau (2D) de forcage :
!-----------------------------------------------------------------------

      if (nbcor.lt.0) then
!- modif de la totalite d'un Tab.2D (Nouvelles valeurs lues sur fichier) :
        read(filcor_id,'(A)') filtab
        icor1 = 1
        icor2 = imax
        jcor1 = 1
        jcor2 = jmax
        kmodif = kmodif + 3
!- lecture(1 niveau complet) du fichier "filtab" :
        open(newunit=filtab_id,file=filtab,status='OLD',form='UNFORMATTED')
        read(filtab_id) spvcor
        read(filtab_id) ycor
        close(filtab_id)
        if (nn99.eq.2) write(mouchard_id,'(2A)')
     &                ' -> modif globale, fichier = ', filtab
      elseif (nvcor.eq.6) then
!-----
!- ajoute (dans liste specifique) un forcage local :
!- NB : read i1,i2 + j1,j2 + k1,k2 + Y_new
        read(filcor_id,*) icor1, icor2, jcor1, jcor2, kcor1, kcor2, ycor(1,1)
        icor2 = max(icor1,icor2,-imax*icor2)
        icor1 = max(icor1,1)
        jcor2 = max(jcor1,jcor2,-jmax*jcor2)
        jcor1 = max(jcor1,1)
        kcor2 = max(kcor1,kcor2,-kmax*kcor2)
        kcor1 = max(kcor1,1)
        spvcor = spv
!-----
      else
!- modif d'un Tab.2D, par portion :
        read(filcor_id,*) icor1, icor2, jcor1, jcor2, ycor(1,1)
        icor2 = max(icor1,icor2,-imax*icor2)
        icor1 = max(icor1,1)
        jcor2 = max(jcor1,jcor2,-jmax*jcor2)
        jcor1 = max(jcor1,1)
        spvcor = spv
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Traitements des differents cas de tableaux a modifier :
!-----------------------------------------------------------------------

      if (nvcor.eq.1) then
          if (kcor.lt.1 .or. kcor.gt.2) goto 910
          call cortab(phisuv(1,1,kcor),ycor,spv,spvcor,
     &                imax,icor1,icor2,jcor1,jcor2,kmodif,kspv)
          if (nn99.eq.2) write(mouchard_id,1200) ' modif phiss(ns) : i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ; u/v,ktyp', kcor,ktyp,
     &      ' ; V.cor=', ycor(1,1)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      elseif (nvcor.eq.2) then
          if (kcor.lt.0 .or. kcor.gt.nsmax) goto 910
          call cortab(phiss(1,1,kcor),ycor,spv,spvcor,
     &                imax,icor1,icor2,jcor1,jcor2,kmodif,kspv)
          if (nn99.eq.2) write(mouchard_id,1200) ' modif phiss(ns) : i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ;  ns,ktyp', kcor,ktyp,
     &      ' ; V.cor=', ycor(1,1)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      elseif (nvcor.eq.3) then
!-- kcor = ns * k
          if (kcor.lt.1 .or. kcor.gt.nsmax*kmax) goto 910
          k = kcor - 1
          ns = k / kmax + 1
          k = mod(k,kmax) + 1
          call cortab(scalr(1,1,k,ns),ycor,spv,spvcor,
     &                imax,icor1,icor2,jcor1,jcor2,kmodif,kspv)
          if (nn99.eq.2) write(mouchard_id,1300) ' modif scalr(ns) : i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ; k,ns,ktyp', k,ns,ktyp,
     &      ' ; V.cor=', ycor(1,1)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      elseif (nvcor.eq.4) then
          if (kcor.lt.1 .or. kcor.gt.kmax) goto 910
          call cortab(rappel(1,1,kcor),ycor,spv,spvcor,
     &                imax,icor1,icor2,jcor1,jcor2,kmodif,kspv)
          if (nn99.eq.2) write(mouchard_id,1200) ' modif rappel(k) : i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ;  k, ktyp', kcor,ktyp,
     &      ' ; V.cor=', ycor(1,1)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      elseif (nvcor.eq.5) then
          if (kcor.lt.0 .or. kcor.gt.nsmax) goto 910
          call cortab(rappes(1,1,kcor),ycor,spv,spvcor,
     &                imax,icor1,icor2,jcor1,jcor2,kmodif,kspv)
          if (nn99.eq.2) write(mouchard_id,1200) ' modif rappes(ns): i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ;  ns,ktyp', kcor,ktyp,
     &      ' ; V.cor=', ycor(1,1)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      elseif (nvcor.eq.6) then
!---Rappel (explic.) autre qu'en surface (--> tableau specifique) :
          if (kcor.lt.1 .or. kcor.gt.nsmax) goto 910
          ns = kcor
          nnadd = (icor2 - icor1 + 1) * (jcor2 - jcor1 + 1)
          do k=kcor1,kcor2
            nn = nrap(k,ns)
            nrap(k,ns) = nrap(k,ns) + nnadd
            nrapmx = max(nrapmx, nrap(k,ns))
            if (nrapmx.le.nrpmax) then
              do j=jcor1,jcor2
               do i=icor1,icor2
                 nn = nn + 1
                 ijrap(nn,k,ns) = (j - 1)*imax + i - 1
                 rapint(nn,k,ns) = ycor(1,1)
               end do
              end do
            endif
          end do
          if (nn99.eq.2) write(mouchard_id,1400) ' rapint(k,ns): i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ; k', kcor1,kcor2,
     &      ' ; ns,typ', ns,ktyp, ' ; V.cor=', ycor(1,1)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      elseif (nvcor.eq.7) then
!-- kcor = 0,1 min/Max Frsh-W-Flx ; 2,3 min/Max Heat-Flx ; 4,5 min/Max Salt-Flx
          if (kcor.lt.0 .or. kcor.ge.3*nsmax) goto 910
          ns = kcor/2
          nx = mod(kcor,2)
          call cortab(phimnx(1,1,nx,ns),ycor,spv,spvcor,
     &                imax,icor1,icor2,jcor1,jcor2,kmodif,kspv)
          if (nn99.eq.2) write(mouchard_id,1300) ' modif phimnx(ns): i',
     &      icor1,icor2,' ; j', jcor1,jcor2, ' ; n,ns,ktyp',nx,ns,ktyp,
     &      ' ; V.cor=', ycor(1,1)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      else
          goto 900
      endif

!--fin du traitement d'un tableau 2 D.
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      end do

!--fin du traitement de la modificcation.
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      end do

      close(filcor_id)

      phisu(:,:) = phisuv(:,:,1)
      phisv(:,:) = phisuv(:,:,2)

      if (nrapmx.le.nrpmax) then
        write(clio3_out_id,'(2A)') 'Forcage Modifie par le fichier ', filcor
        return
      else
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  9 ) Traitements des cas d'erreurs .
!-----------------------------------------------------------------------

        write(clio3_out_id,'(A,I10,A)')'STOP in "correct", dimension nrpmax=',
     &              nrpmax, '= TOO SMALL !'
        write(clio3_out_id,'(2A)')     '=> Change file : ', filcor
        write(clio3_out_id,
     &           '(A,I10,A)')'or set nrpmax (in "para.Com") to at least'
     &            , nrapmx, ' and compile the code again.'
        stop
      endif

 900  continue
      write(clio3_out_id,'(2A)') 'STOP in "correct", ERROR in : ',filcor
      write(clio3_out_id,'(2A)') ' modif : ', line
      write(clio3_out_id,'(A,I8,A)') ' var. nvcor=', nvcor, '  <-- Not found !'
      stop

 910  continue
      write(clio3_out_id,'(2A)') 'STOP in "correct", ERROR in : ',filcor
      write(clio3_out_id,'(2A)') ' modif : ', line
      write(clio3_out_id,'(A,I3,A,I8,A)') ' var. nvcor=', nvcor,
     &   ' index kcor=', kcor, '  <-- out of range !'
      stop

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine correct -
      end
