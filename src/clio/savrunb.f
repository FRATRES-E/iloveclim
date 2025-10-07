!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:52 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:52 CET 2009

      SUBROUTINE savrunb(nnt,nn99,ccfile)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  sortie des resultats sur fichier binaire
!     1ere partie utilisable pour faire redemarrer le programme ;
!  fichier de sortie 'res/n/.om' avec "n" decroissant jusqu a zero
!-----
!  modif : 02/07/98 ; cleanup, dmr, 2021-05-24
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!! START_OF_USE_SECTION

      use para0_mod,  only: imax, jmax, kmax, nsmax, ixjmax
      use para_mod,   only: nvmax
      use moment_mod, only: vicmom, copy_to_vicmon, copy_from_vicmon
      use vareq_mod,  only: uslpfx,vslpfx

      use bloc0_mod,  only: tpstot, scal, fss, eta, ub, vb, fqajc, b, u
     &  , v, q2turb, w, bvf, avsdz, avudz

      use bloc_mod  , only: nstart, numit, refexp, umoy, vmoy

!--- dmr [UNUSED] nvreac, nvrfc, nvrfs, nvrhac, nvrs, nvrt

      use datadc_mod, only: nvrajc, nvral, nvras, nvrau, nvrb
     &  , nvret, nvrfq, nvrfw, nvrhg, nvrhn
     &  , nvrmom, nvrn2, nvrqs, nvrtbq, nvrtgx, nvrtgy
     &  , nvrtke, nvrts, nvru, nvrub, nvrug, nvrum, nvrusl, nvrv, nvrvb
     &  , nvrvg, nvrvm, nvrvsl, nvrw, nvrxzo

      use ice_mod,    only: albq, fsbbq, hgbq, hnbq, qstobq, tbq, tenagx
     &  , tenagy, ts, xzo
      use dynami_mod, only: ug, vg
      use varno_mod,  only: nvrl, titcv, ltyp, krlm

      use global_constants_mod, only: dblp=>dp, ip

      use newunit_clio_mod, only: clio3_out_id, mouchard_id

!! END_OF_USE_SECTION

      implicit none


!- dummy variables :
      integer(ip),   intent(in) :: nnt, nn99
      character*(*), intent(in) :: ccfile

!- variables locales :
      integer(ip), dimension(nvmax) :: nnvv
      character*120 ccline

      integer(ip) :: kkerr, kksize, n, nbvar, nnc, nnc1, nnct, nnfbin 
     >             , ns, nv

      integer(ip) :: ccfile4_id
      
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Preparation de l'Ecriture ; definition et ecriture de l'entete  |
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      nnct = len(ccline)
      ccline(:6) = 'Wbin :'
      nnc = 6

!- calcul du nb de tableaux a ecrire :
      nbvar = 0
      do nv=1,nvmax
        if (mod(nvrl(nv),2).eq.1) then
          nbvar = nbvar + 1
          nnvv(nbvar) = nv
        endif
      enddo

      if (nbvar.eq.0) then
        write(clio3_out_id,*)
     &    'savrunb : aucun tableau a ecrire sur fichier !'
        return
      endif

!--Ouverture et ecriture de l'en-tete

      open(newunit=ccfile4_id,file=ccfile,
     &     status='unknown',form='UNFORMATTED')
     
      nnfbin = ccfile4_id
      
      write(ccfile4_id) numit, tpstot, refexp
      write(ccfile4_id) nbvar

      do n=1,nbvar
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Ecriture sur fichier variable par variable .
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


      nv = nnvv(n)
!     write(clio3_out_id,'(A,2I4)') titcv(nv), nv, nvrl(nv)

!- taille du tableau (=Nb d'elements) a ecrire :
      kksize = krlm(nv)
      if (ltyp(nv).ge.0) kksize = kksize * ixjmax

!- Ecriture du No="nv", de la taille et du nom (3cc) de la variable "nv"
      write(ccfile4_id) nv, kksize, titcv(nv)

!- Recherche du vrai tableau et ecriture par appel a la routine "savtab" :
      if (nv.le.nsmax) then
           call savtab(scal(1,1,1,nv),kksize,nnfbin,kkerr)
      elseif (nv.ge.nvrfw .and. nv.le.nvrfw+nsmax) then
           ns = nv-nvrfw
           call savtab(fss(1,1,ns),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvret) then
           call savtab(eta(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrub) then
           call savtab(ub(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrvb) then
           call savtab(vb(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrajc) then
           call savtab(fqajc(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrb) then
           call savtab(b(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvru) then
           call savtab(u(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrv) then
           call savtab(v(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrtke) then
           call savtab(q2turb(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrw) then
           call savtab(w(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrn2) then
           call savtab(bvf(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvras) then
           call savtab(avsdz(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrau) then
           call savtab(avudz(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrusl) then
           call savtab(uslpfx(1,1,-1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrvsl) then
           call savtab(vslpfx(1,1,-1),kksize,nnfbin,kkerr)
!-----
      elseif (nv.eq.nvrum) then
           call savtab(umoy(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrvm) then
           call savtab(vmoy(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrhg) then
           call savtab(hgbq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrfq) then
           call savtab(fsbbq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrqs) then
           call savtab(qstobq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvral) then
           call savtab(albq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrhn) then
           call savtab(hnbq(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrts) then
           call savtab(ts(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrug) then
           call savtab(ug(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrvg) then
           call savtab(vg(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrtbq) then
           call savtab(tbq(1,1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrxzo) then
           call savtab(xzo(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrtgx) then
           call savtab(tenagx(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrtgy) then
           call savtab(tenagy(1,1),kksize,nnfbin,kkerr)
      elseif (nv.eq.nvrmom) then
           call copy_to_vicmon()
           call savtab(vicmom(1,1,1),kksize,nnfbin,kkerr)
!----
      else
        write(clio3_out_id,'(3A,I3,A)') 'ARRET, savrunb : Ecriture de ',
     &         titcv(nv), ' nv=', nv, ' Pas Prevue !'
        goto 900
      endif
      if (kkerr.eq.-1) then
        write(clio3_out_id,'(3A,I3,A)') 'ARRET, savrunb : Error writing ',
     &         titcv(nv), ' nv=', nv
        goto 900
      endif

!--Consigne les noms des variables ecrites sur "ccline" :
      nnc1 = 1 + nnc
      nnc = min(3+nnc1,nnct)
      if(nnc1.le.nnc) ccline(nnc1:nnc) = ' '//titcv(nv)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--fin de l'ecriture de la nariable "nv".
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      enddo 
      
      close(ccfile4_id)

      write(clio3_out_id,
     & '(A,I11,A,I3,2A)') 'Exp '//refexp//', Iter', numit,
     &  ' :', nbvar, ' variables written  on  File ', ccfile

      if (nn99.eq.2) then
        write(clio3_out_id,'(A)') ccline(:nnc)
        write(clio3_out_id,'(A,I11,A,I3,2A)')
     &    'Exp '//refexp//', Iter', numit,
     &  ' :', nbvar, ' variables written  on  File ', ccfile
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  7 ) Remise a zero des variables moyennnes et extremes .
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      if (nnt.ge.2 .and. numit.ge.nstart) then

        fqajc(:,:,:) = 0.0
        fss(:,:,:) = 0.0

      endif

      return

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  9 ) Traitement des cas "problematiques" .                           |
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


 900  continue
      write(clio3_out_id,'(A,I11,A,I3,2A)')
     &  'Exp '//refexp//', Iter', numit,
     &  ' :', n-1, ' var. written & STOP,   File ', ccfile
      write(clio3_out_id,'(A)') ccline(:nnc)

      if (nn99.eq.2) close(mouchard_id)
      close(ccfile4_id)

      stop

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine savrunb -
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      end subroutine savrunb
