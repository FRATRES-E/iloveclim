!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:51 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:51 CET 2009

      SUBROUTINE redrunb(nnt,nn99,ccfile)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  redemarrage a partir de l'etat definit par le fichier binaire "rest.om" .
! Entree : nnt = 0,1  => lecture du fichier de resultat Nouveau(>05/96) format
!          nnt > 3    => lecture du fichier de resultat Ancien( <05/96) format
!              = 4,5  => lit  ancien fichier de resultat simple.
!              = 6,7  => lit  ancien fichier de resultat complet (sans w).
!              = 8,9  => lit +ancien(<01/96) fich. res.  complet (avec w).
!----- nnt(Sortie) = nnt(Entree)
!  modif : 02/07/98

!! START_OF_USE_SECTION

      use const_mod,  only: zero, rho0

!nb & fl      use para0_mod,  only: imax, jmax, kmax, ixjmax, nsmax, ijkmax
      use para0_mod,  only: imax, jmax, kmax, ixjmax, nsmax_TS, ijkmax
      use para_mod,   only: nvmax
      use moment_mod, only: sxa, sxc0, sxc1, sxc2, sxg, sxn, sxst, sxxa
     &  , sxxc0, sxxc1, sxxc2, sxxg, sxxn, sxxst, sxya, sxyc0, sxyc1
     &  , sxyc2, sxyg, sxyn, sxyst, sya, syc0, syc1, syc2, syg, syn
     &  , syst, syya, syyc0, syyc1, syyc2, syyg, syyn, syyst, vicmom
     &  , copy_to_vicmon, copy_from_vicmon

      use vareq_mod,  only: uslpfx, vslpfx

      use bloc0_mod,  only: avsdz, avudz, b, bvf, eta, fqajc, scal
     &  , tpstot, u, ub, v, vb, w, fss, q2turb, tmu, ks2, dz
      use bloc_mod,   only: umoy, vmoy, q, refexp, numit, unsvol, aire
      use datadc_mod, only: nvrajc, nvral, nvras, nvrau, nvrb, nvreac
     &  , nvret, nvrfc, nvrfq, nvrfs, nvrfw, nvrhac, nvrhg, nvrhn
     &  , nvrmom, nvrn2, nvrqs, nvrs, nvrt, nvrtbq, nvrtgx, nvrtgy
     &  , nvrtke, nvrts, nvru, nvrub, nvrug, nvrum, nvrusl, nvrv, nvrvb
     &  , nvrvg, nvrvm, nvrvsl, nvrw, nvrxzo
      use ice_mod,    only: albq, fsbbq, hgbq, hnbq, qstobq, tbq, tenagx
     &  , tenagy, ts, xzo
      use dynami_mod, only: ug, vg
      use varno_mod,  only: titcv, nvrl, ltyp, krlm

#if ( BATHY >= 1 )
      use update_clio_bathy_tools, only: mean_neighbours_with_mask
     & , mise_a_zero, tms_prev, salglobal, aire_prev, dz_prev
      use bloc0_mod,   only: tms
      use bloc_mod,    only: q2tmin
#endif

      use global_constants_mod, only: dblp=>dp, ip

      use newunit_clio_mod, only: clio3_out_id
      !use newunit_clio_mod, only: clio3_out_id, mouchard_id

!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "varno.com"
!dmr [NOEQUI]   #include "vareq.com"
! [SCRPTCOM] #include "ice.com"
! [SCRPTCOM] #include "dynami.com"

!dmr [NOEQUI] moment.com provides vicmom that maps out to 35 different variables
! [SCRPTCOM] #include "moment.com"

!! END_OF_INCLUDE_SECTION

      integer(ip),  parameter :: nszmax = ijkmax

!- dummy variables :
       character*(*) ccfile
!- variables locales equivalentes :
       real(dblp), dimension(nszmax) :: vloc = 0.0_dblp !dmr on dirait que vloc est intent(out) ...

!dmr [NOEQUI]      equivalence ( q(1,1,1), vloc(1) ) ! q(:,:,:) is not used in this subroutine

!dmr [NOEQUI]!- pour ancien fichier :
!dmr [NOEQUI]      dimension wloc(imax,jmax,kmax)
!dmr [NOEQUI]      equivalence ( w(1,1,1) ,  wloc(1,1,1) )

!- local variables :
       character*3 cc3
       character*6 cc6exp
       character*120 ccline

       integer(ip) :: nnt, nn99,  i, j, k, kkerr, kkr, kksize, llost, n
     & , nbvar, nnc, nnc1, nnct, nnr, nnvrd, ns, nv, nb_var

       real(dblp) :: egajc, hmajc

!--instructions "data" :
! [SCRPTCOM] #include "datadc.com"

      logical :: logic_elmt

#if ( BATHY >= 1 )
!nb to use for global salinity check
      real(kind=dblp) :: volglobal
      real(kind=dblp) :: salglobal_time, saldiff, scaldiff
      real(kind=dblp) :: salglobal_time_needed
      real(kind=dblp) :: volglobal_time, voldeep
#endif

      integer(ip) :: ccfile3_id
      
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initialisation & Ouverture du fichier binaire a lire .          |
!-----------------------------------------------------------------------

!- Unite(fortran) du Fichier "ccfile" a lire :
      nnct = len(ccline)
      ccline(:6) = 'Rbin :'
      nnc = 6

      open(newunit=ccfile3_id,file=ccfile,status='old',form='UNFORMATTED')

      if (nnt.ge.4) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Lecture, fichier binaire "ccfile", Ancien format                |
!-----------------------------------------------------------------------

!--lecture de la premiere partie :
      read (ccfile3_id, end=281) numit
      read (ccfile3_id, end=281) tpstot
      read (ccfile3_id, end=281) eta
      read (ccfile3_id, end=281) ub
      read (ccfile3_id, end=281) vb
      read (ccfile3_id, end=281) u
      read (ccfile3_id, end=281) v
      read (ccfile3_id, end=281) scal

      read (ccfile3_id, end=281) umoy
      read (ccfile3_id, end=281) vmoy
      read (ccfile3_id, end=281) hgbq
      read (ccfile3_id, end=281) fsbbq
      read (ccfile3_id, end=281) qstobq
      read (ccfile3_id, end=281) albq
      read (ccfile3_id, end=281) hnbq
      read (ccfile3_id, end=281) ts
      read (ccfile3_id, end=281) ug
      read (ccfile3_id, end=281) vg
      read (ccfile3_id, end=281) tbq
      read (ccfile3_id, end=281) xzo
      read (ccfile3_id, end=281) tenagx
      read (ccfile3_id, end=281) tenagy
      read (ccfile3_id, end=281) sxg
      read (ccfile3_id, end=281) syg
      read (ccfile3_id, end=281) sxxg
      read (ccfile3_id, end=281) syyg
      read (ccfile3_id, end=281) sxyg
      read (ccfile3_id, end=281) sxn
      read (ccfile3_id, end=281) syn
      read (ccfile3_id, end=281) sxxn
      read (ccfile3_id, end=281) syyn
      read (ccfile3_id, end=281) sxyn
      read (ccfile3_id, end=281) sxa
      read (ccfile3_id, end=281) sya
      read (ccfile3_id, end=281) sxxa
      read (ccfile3_id, end=281) syya
      read (ccfile3_id, end=281) sxya
      read (ccfile3_id, end=281) sxc0
      read (ccfile3_id, end=281) syc0
      read (ccfile3_id, end=281) sxxc0
      read (ccfile3_id, end=281) syyc0
      read (ccfile3_id, end=281) sxyc0
      read (ccfile3_id, end=281) sxc1
      read (ccfile3_id, end=281) syc1
      read (ccfile3_id, end=281) sxxc1
      read (ccfile3_id, end=281) syyc1
      read (ccfile3_id, end=281) sxyc1
      read (ccfile3_id, end=281) sxc2
      read (ccfile3_id, end=281) syc2
      read (ccfile3_id, end=281) sxxc2
      read (ccfile3_id, end=281) syyc2
      read (ccfile3_id, end=281) sxyc2
      read (ccfile3_id, end=281) sxst
      read (ccfile3_id, end=281) syst
      read (ccfile3_id, end=281) sxxst
      read (ccfile3_id, end=281) syyst
      read (ccfile3_id, end=281) sxyst


!-----
      nvrl(nvret) = mod(nvrl(nvret),2) + 2
      nvrl(nvrub) = mod(nvrl(nvrub),2) + 2
      nvrl(nvrvb) = mod(nvrl(nvrvb),2) + 2
      nvrl(nvru)  = mod(nvrl(nvru),2)  + 2
      nvrl(nvrv)  = mod(nvrl(nvrv),2)  + 2
      nvrl(nvrt)  = mod(nvrl(nvrt),2)  + 2
      nvrl(nvrs)  = mod(nvrl(nvrs),2)  + 2

      nvrl(nvrum)  = mod(nvrl(nvrum),2)  + 2
      nvrl(nvrvm)  = mod(nvrl(nvrvm),2)  + 2
      nvrl(nvrhg)  = mod(nvrl(nvrhg),2)  + 2
      nvrl(nvrfq)  = mod(nvrl(nvrfq),2)  + 2
      nvrl(nvrqs)  = mod(nvrl(nvrqs),2)  + 2
      nvrl(nvrhn)  = mod(nvrl(nvrhn),2)  + 2
      nvrl(nvral)  = mod(nvrl(nvral),2)  + 2
      nvrl(nvrts)  = mod(nvrl(nvrts),2)  + 2
      nvrl(nvrug)  = mod(nvrl(nvrug),2)  + 2
      nvrl(nvrvg)  = mod(nvrl(nvrvg),2)  + 2
      nvrl(nvrtbq) = mod(nvrl(nvrtbq),2) + 2
      nvrl(nvrxzo) = mod(nvrl(nvrxzo),2) + 2
      nvrl(nvrtgx) = mod(nvrl(nvrtgx),2) + 2
      nvrl(nvrtgy) = mod(nvrl(nvrtgy),2) + 2
      nvrl(nvrmom) = mod(nvrl(nvrmom),2) + 2

      if (nnt.eq.6 .or.nnt.eq.7) then
!--lecture de la seconde partie :
        read (ccfile3_id, end=282) cc6exp
        read (ccfile3_id, end=282) b
        read (ccfile3_id, end=282) bvf
        read (ccfile3_id, end=282) avsdz
        read (ccfile3_id, end=282) avudz
        read (ccfile3_id, end=282) fqajc
        read (ccfile3_id, end=282)
     &     (((fss(i,j,ns),i=1,imax),j=1,jmax),ns=1,nsmax_TS)
!-----
        nvrl(nvrb)   = mod(nvrl(nvrb),2)   + 2
        nvrl(nvrn2)  = mod(nvrl(nvrn2),2)  + 2
        nvrl(nvras)  = mod(nvrl(nvras),2)  + 2
        nvrl(nvrau)  = mod(nvrl(nvrau),2)  + 2
        nvrl(nvrajc) = mod(nvrl(nvrajc),2) + 2
        nvrl(nvrfc)  = mod(nvrl(nvrfc),2)  + 2
        nvrl(nvrfs)  = mod(nvrl(nvrfs),2)  + 2
!--fin de la seconde partie .
        write(clio3_out_id,
     &   '(A,I11,3A)') 'Exp '//cc6exp//', Iter', numit,
     &      ' , File ', ccfile, ' : FULL RESULTS read - OK.'

      elseif (nnt.ge.8) then
!--lecture de la seconde partie, ancien fichier :
        read (ccfile3_id, end=282) cc6exp
        read (ccfile3_id, end=282) b
        read (ccfile3_id, end=282) bvf
        read (ccfile3_id, end=282) avsdz
        read (ccfile3_id, end=282) avudz
        read (ccfile3_id, end=282) w
!-----
        nvrl(nvrb)   = mod(nvrl(nvrb),2)   + 2
        nvrl(nvrn2)  = mod(nvrl(nvrn2),2)  + 2
        nvrl(nvras)  = mod(nvrl(nvras),2)  + 2
        nvrl(nvrau)  = mod(nvrl(nvrau),2)  + 2
        nvrl(nvrw)   = mod(nvrl(nvrw),2)   + 2
        read (ccfile3_id, end=286) egajc
        nvrl(nvreac) = mod(nvrl(nvreac),2) + 2
        read (ccfile3_id, end=286) hmajc
        nvrl(nvrhac) = mod(nvrl(nvrhac),2) + 2
        read (ccfile3_id, end=287)
     &     (((fss(i,j,ns),i=1,imax),j=1,jmax),ns=1,nsmax_TS)
        nvrl(nvrfc)  = mod(nvrl(nvrfc),2)  + 2
        nvrl(nvrfs)  = mod(nvrl(nvrfs),2)  + 2
!--fin de la seconde partie .
        write(clio3_out_id,
     &   '(A,I11,3A)') 'Exp '//cc6exp//', Iter', numit,
     &      ' , File ', ccfile, ' : OLD FILE read I & II - OK.'

      else
        read (ccfile3_id,end=285) cc6exp
        write(clio3_out_id,
     &   '(A,I11,3A)') 'Exp '//cc6exp//', Iter', numit,
     &      ' , File ', ccfile, ' : Restart (I)  read - OK.'
      endif

      goto 290

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Cas de lecture Incomplete :
 281  continue
        write(clio3_out_id,'(A)') 
     &   ' 1ere partie du fichier Incomplete ! =>  '
        write(clio3_out_id,'(3A)')
     &    'ARRET dans "redrunb", fichier : ', ccfile
        STOP

 282  continue
        write(clio3_out_id,'(A)') 
     &   ' 2nd  partie du fichier Incomplete ! =>  '
        write(clio3_out_id,'(3A)')
     &    'ARRET dans "redrunb", fichier : ', ccfile
        STOP

 285  continue
        write(clio3_out_id,'(3A,I11,A)')
     &    'Lecture de ', ccfile,
     &       ' terminee ( sans Ref., Iter', numit, ' ) .'
        goto 500

 286  continue
        write(clio3_out_id,'(A,I11,3A)')
     &    'Exp '//cc6exp//', Iter', numit,
     &      ' , File ', ccfile, ' : read I & II except E/Hajc & fss!'
        goto 290

 287  continue
        write(clio3_out_id,'(A,I11,3A)')
     &    'Exp '//cc6exp//', Iter', numit,
     &      ' , File ', ccfile, ' : read I & II except fss !'

 290  continue

!dmr [NOEQUI]
      call copy_to_vicmon()

      else
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Lecture, fichier binaire "ccfile", Nouveau format               |
!-----------------------------------------------------------------------

      write(*,*) "PING! Reading restart file", ccfile

!--Lecture de l'en-tete
      read(ccfile3_id) numit, tpstot, cc6exp
      read(ccfile3_id) nbvar

      if (nbvar.gt.nvmax) then
        write(clio3_out_id,'(2(A,I4))')
     &    'ARRET, redrunb : Nb.var=', nbvar,
     &     ' > nvmax=', nvmax
        goto 950
      endif

      llost = 0
      nnvrd = 0
      do n=1,nbvar
!--Lecture variable par variable .

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--3.1 Lecture du No,taille & nom de la variable & verifications :
!-------------
      read(ccfile3_id,end=910,err=920) nv, kkr, cc3

!--Verification :
      if (nv.lt.1 .or. nv.gt.nvmax) then
        write(clio3_out_id,'(2(A,I3))')
     &    'ARRET, redrunb : Var. nv= ', nv,
     &     ' out of range 1,nvmax=', nvmax
        goto 950
      elseif (titcv(nv).ne.cc3) then
        write(clio3_out_id,'(A,I3,4A)')
     &    'ARRET, redrunb : nv=', nv,
     &        ', Var.= ',titcv(nv), ' <-> on File= ', cc3
        goto 950
      endif
!- taille (supposee) du tableau (=Nb d'elements) :
      kksize = krlm(nv)
      if (ltyp(nv).ge.0) kksize = kksize * ixjmax
      if (kksize.gt.kkr) then
        write(clio3_out_id,'(A,I3,2A,2(A,I8))')
     &    'WARNING, redrunb : nv=', nv,
     &        ', Var.= ', titcv(nv), ', Nb.Elm. to read', kksize,
     &        ' > on file', kkr
        kksize = kkr
      elseif (kksize.eq.0) then
        write(clio3_out_id,'(3A,I3,A)')
     &    'WARNING, redrunb : Lecture de ',
     &         titcv(nv), ' nv=', nv, ' Pas Prevue (Tab.Vide) !'
        llost = llost + 1
        kksize = 1
        call redtab(vloc(1),kksize,ccfile3_id,kkerr)
        vloc(1) = 0.
        nnr = mod(nvrl(nv),2)
        goto 320
      elseif (kksize.lt.kkr) then
        write(clio3_out_id,'(A,I3,2A,2(A,I8))')
     &    'WARNING, redrunb : nv=', nv,
     &        ', Var.= ', titcv(nv), ', Nb.Elm.lu=', kksize,
     &        ' < on file=', kkr
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--3.2 Recherche du vrai tableau ; lecture par appel a "redtab" :
!-------------
      nnr = mod(nvrl(nv),2) + 2
!nb & fl      if (nv.le.nsmax) then
      if (nv.le.nsmax_TS) then
           call redtab(scal(1,1,1,nv),kksize,ccfile3_id,kkerr)
      elseif (nv.ge.nvrfw .and. nv.le.nvrfw+nsmax_TS) then ! nb & fl changed nsmax
           ns = nv-nvrfw
           call redtab(fss(1,1,ns),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvret) then
           call redtab(eta(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrub) then
           call redtab(ub(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrvb) then
           call redtab(vb(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrajc) then
           call redtab(fqajc(1,1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrb) then
           call redtab(b(1,1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvru) then
           call redtab(u(1,1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrv) then
           call redtab(v(1,1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrtke) then
           call redtab(q2turb(1,1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrw) then
           call redtab(w(1,1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrn2) then
           call redtab(bvf(1,1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvras) then
           call redtab(avsdz(1,1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrau) then
           call redtab(avudz(1,1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrusl) then
           call redtab(uslpfx(1,1,-1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrvsl) then
           call redtab(vslpfx(1,1,-1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrum) then
           call redtab(umoy(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrvm) then
           call redtab(vmoy(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrhg) then
           call redtab(hgbq(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrfq) then
           call redtab(fsbbq(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrqs) then
           call redtab(qstobq(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvral) then
           call redtab(albq(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrhn) then
           call redtab(hnbq(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrts) then
           call redtab(ts(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrug) then
           call redtab(ug(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrvg) then
           call redtab(vg(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrtbq) then
           call redtab(tbq(1,1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrxzo) then
           call redtab(xzo(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrtgx) then
           call redtab(tenagx(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrtgy) then
           call redtab(tenagy(1,1),kksize,ccfile3_id,kkerr)
      elseif (nv.eq.nvrmom) then
           call redtab(vicmom(1,1,1),kksize,ccfile3_id,kkerr)

!dmr [NOEQUI]
           call copy_from_vicmon()
!dmr [NOEQUI]
!-----
      else
!- cas d'une variable dont la lecture n'est pas prevue :
        kksize = 1
        call redtab(vloc(1),kksize,ccfile3_id,kkerr)
!- reinitialise a zero "vloc" :
        vloc(1) = 0.
        write(clio3_out_id,'(3A,I3,A)')
     &    'WARNING, redrunb : Lecture de ',
     &         titcv(nv), ' nv=', nv, ' Pas Prevue !'
        nnr = nnr - 2
        llost = llost + 1
!-----
      endif

#if ( BATHY >= 1 )
!nb Here compute salinity globale before
      if (nv.eq.2) then !if salinity
      volglobal=0.0
      salglobal=0.
      do i=1,imax
        do j=1,jmax
          do k=1,kmax
          volglobal=volglobal+aire_prev(i,j)*dz_prev(k)*tms_prev(i,j,k)
          salglobal=salglobal+aire_prev(i,j)*dz_prev(k)*tms_prev(i,j,k)
     &              *scal(i,j,k,2)*rho0 !m2*m*g/kg*kg/m3 -> g
          enddo
        enddo
      enddo
      write(*,*) 'volume ocean et salinite globale', volglobal,
     &           salglobal, salglobal/volglobal/rho0
      endif
!#endif

!#if ( BATHY >= 1 )
      if (nv.le.nsmax_TS) then
           logic_elmt = mean_neighbours_with_mask(scal(:,:,:,nv),1) ! last argument 1 = tms2D
      elseif (nv.ge.nvrfw .and. nv.le.nvrfw+nsmax_TS) then ! nb & fl
changed nsmax
           ns = nv-nvrfw
           logic_elmt = mean_neighbours_with_mask(fss(:,:,ns),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvret) then
           logic_elmt = mean_neighbours_with_mask(eta(:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvrub) then
           logic_elmt = mise_a_zero(ub(:,:),2) ! last argument 1 = tms2D
           ub(:,:) = 0.0d0
      elseif (nv.eq.nvrvb) then
           logic_elmt = mise_a_zero(vb(:,:),2) ! last argument 1 = tms2D
           vb(:,:) = 0.0d0
      elseif (nv.eq.nvrajc) then
           logic_elmt = mise_a_zero(fqajc(:,:,:),2) ! last argument 1 = tms2D
      elseif (nv.eq.nvrb) then
           logic_elmt = mean_neighbours_with_mask(b(:,:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvru) then
           logic_elmt = mise_a_zero(u(:,:,:),2) ! last argument 1 = tms2D
           u(:,:,:) = 0.0d0
      elseif (nv.eq.nvrv) then
           logic_elmt = mise_a_zero(v(:,:,:),2) ! last argument 1 = tms2D
           v(:,:,:) = 0.0d0
      elseif (nv.eq.nvrtke) then
           ! logic_elmt = mise_a_zero(q2turb(:,:,:),1) ! last argument 1 = tms2D !!!! HUUUUUUUMMMMMM PAS SUR DU TOUT !!!!!
           q2turb = q2tmin
      elseif (nv.eq.nvrw) then
           logic_elmt = mise_a_zero(w(:,:,:),2) ! last argument 1 = tms2D
           w(:,:,:) = 0.0d0
      elseif (nv.eq.nvrn2) then
           logic_elmt = mean_neighbours_with_mask(bvf(:,:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvras) then
           logic_elmt = mean_neighbours_with_mask(avsdz(:,:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvrau) then
           logic_elmt = mean_neighbours_with_mask(avudz(:,:,:),2) ! last argument 1 = tms2D
      elseif (nv.eq.nvrusl) then
           logic_elmt = mise_a_zero(uslpfx(:,:,:),2) ! last argument 1 = tms2D
      elseif (nv.eq.nvrvsl) then
           logic_elmt = mise_a_zero(vslpfx(:,:,:),2) ! last argument 1 = tms2D
      elseif (nv.eq.nvrum) then
           logic_elmt = mise_a_zero(umoy(:,:),2) ! last argument 1 = tms2D
           umoy(:,:) = 0.0d0
      elseif (nv.eq.nvrvm) then
           logic_elmt = mise_a_zero(vmoy(:,:),2) ! last argument 1 = tms2D
           vmoy(:,:) = 0.0d0
      elseif (nv.eq.nvrhg) then
           logic_elmt = mean_neighbours_with_mask(hgbq(:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvrfq) then
           logic_elmt = mise_a_zero(fsbbq(:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvrqs) then
           logic_elmt = mise_a_zero(qstobq(:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvral) then
           logic_elmt = mean_neighbours_with_mask(albq(:,:), 1) ! last argument 1 = tms2D
      elseif (nv.eq.nvrhn) then
           logic_elmt = mean_neighbours_with_mask(hnbq(:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvrts) then
           logic_elmt = mean_neighbours_with_mask(ts(:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvrug) then
           logic_elmt = mise_a_zero(ug(:,:),2) ! last argument 1 = tms2D
           ug(:,:) = 0.0d0
      elseif (nv.eq.nvrvg) then
           logic_elmt = mise_a_zero(vg(:,:),2) ! last argument 1 = tms2D
           vg(:,:) = 0.0d0
      elseif (nv.eq.nvrtbq) then
           logic_elmt = mean_neighbours_with_mask(tbq(:,:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvrxzo) then
           logic_elmt = mean_neighbours_with_mask(xzo(:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvrtgx) then
           logic_elmt = mean_neighbours_with_mask(tenagx(:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvrtgy) then
           logic_elmt = mean_neighbours_with_mask(tenagy(:,:),1) ! last argument 1 = tms2D
      elseif (nv.eq.nvrmom) then
           do nb_var=1,UBOUND(vicmom,DIM=3) 
c           logic_elmt = mean_neighbours_with_mask(vicmom(:,:,:),1) ! last argument 1 = tms2D
             logic_elmt = mean_neighbours_with_mask(vicmom(:,:,nb_var),1) ! last argument 1 = tms2D
           enddo
      endif
!#endif

!#if ( BATHY >= 1 )
!nb&aq check global salinity
      if (nv.eq.2) then !if salinity
      volglobal_time=0.
      salglobal_time=0.
      salglobal_time_needed=0.
      voldeep=0.
      do i=1,imax
        do j=1,jmax
          do k=1,kmax
          volglobal_time=volglobal_time+aire(i,j)*dz(k)*tms(i,j,k)
          salglobal_time=salglobal_time+aire(i,j)*dz(k)*
     &     tms(i,j,k)*scal(i,j,k,2)*rho0 !m2*m*g/kg*kg/m3 -> g
           if (k.le.7) then ! only deep ocean (k=1=bottom, k=7-> 1000m)
             voldeep=voldeep+aire(i,j)*dz(k)*tms(i,j,k)
           endif
          enddo
        enddo
      enddo
!      write(*,*) 'volume ocean, salinite globale et salinity
!     &            global now thersf, mean salinity before and now', volglobal, 
!     &            salglobal, salglobal_time, salglobal/volglobal,
!     &            salglobal_time/volglobal

!      write(*,*) 'salinity before and now, mean salinity',
!     & salglobal, salglobal_time, salglobal_time/volglobal/rho0

#if ( F_PALAEO_FWF == 0 || APPLY_UNCORFWF == 0 )
!nb if no flux from ice sheet elevation change 
!or a corrected flux
!we modify the salinity
!and keep total salt mass constant
      !saldiff=0.
      saldiff=(salglobal_time-salglobal) !in g

#else
!nb if there is a fresh water flux (uncorrected) coming from ice sheets we want to
!keep global mean salinity the same
! so we want salinity_now= mass_salt_now/Vnow 
!                        = salinity_before = mass_salt_before/Vbefore
       !salglobal_time_needed= salglobal/volglobal*volglobal_time
       salglobal_time_needed= salglobal*(volglobal_time/volglobal)
       saldiff=(salglobal_time-salglobal_time_needed) !in g 
#endif

      !write(*,*) 'saldiff ', saldiff
      write(*,*) 'ratio volume now / volume before ',
     &             volglobal_time/volglobal
      write(*,*) 'ratio volume before / volume now ', 
     &             volglobal/volglobal_time
      write(*,*) 'ratio sal', salglobal_time/salglobal
      write(*,*) 'salglobal now, before', salglobal_time, salglobal
    

      if (abs(saldiff).gt.0) then
!      write(*,*) 'saldiff ', saldiff

      scaldiff=saldiff/voldeep/rho0


       do i=1,imax
        do j=1,jmax
!          do k=1,kmax
          do k=1,7
           if ((tms(i,j,k).ne.0).and.(aire(i,j).gt.0)) then
!          write(*,*) 'test ', scal(i,j,k,2), 
!     &                saldiff*aire(i,j)*dz(k)*tms(i,j,k)/volglobal,
!     &                saldiff*aire(i,j)*dz(k)*tms(i,j,k)/volglobal
!     &                  /(aire(i,j)*dz(k)*tms(i,j,k)*rho0)
!     &                  aire(i,j),dz(k),tms(i,j,k),i,j

!          scal(i,j,k,2)=scal(i,j,k,2)-(saldiff
!     &                  *aire(i,j)*dz(k)*tms(i,j,k)/volglobal !g
!     &                  /(aire(i,j)*dz(k)*tms(i,j,k)*rho0))!/(m2*m*kg/m3) ->g/kg

          scal(i,j,k,2)=scal(i,j,k,2)-scaldiff
!     &                  *(aire(i,j)*dz(k)*tms(i,j,k))/voldeep !g
!     &                  /(aire(i,j)*dz(k)*tms(i,j,k)*rho0))!/(m2*m*kg/m3) ->g/kg

           endif
          enddo
        enddo
       enddo
      endif



! verif
      salglobal_time=0.
       do i=1,imax
        do j=1,jmax
          do k=1,kmax
          salglobal_time=salglobal_time+aire(i,j)*dz(k)*
     &     tms(i,j,k)*scal(i,j,k,2)*rho0 !m2*m*g/kg*kg/m3 -> g
          enddo
        enddo
      enddo

      write(*,*) 'verif correction salinity salglobal new meansalinity',
     & salglobal, salglobal_time, salglobal_time/volglobal_time/rho0

      write(*,*) 'ratio sal', salglobal_time/salglobal
      write(*,*) 'salglobal now modif, before', salglobal_time, salglobal
       
      endif !if salinity

#endif



 320  continue

      if (kkerr.eq.-1) then
        write(clio3_out_id,'(3A,I3,A)')
     &    'ARRET, redrunb : Error reading ',
     &         titcv(nv), ' nv=', nv
        goto 950
      elseif (kkerr.eq.0) then
        write(clio3_out_id,'(3A,I3,A)')
     &    'ARRET, redrunb : End_of_File reading ',
     &         titcv(nv), ' nv=', nv
        goto 950
      endif
!- fin de la lecture de la "n"ieme variable .
      if (ltyp(nv).ne.99) nvrl(nv) = nnr

      if (nnr.ge.2) then
!--Consigne les noms des variables luees sur "ccline" :
        nnvrd = nnvrd + 1
        nnc1 = 1 + nnc
        nnc = min(3+nnc1,nnct)
        if(nnc1.le.nnc) ccline(nnc1:nnc) = ' '//titcv(nv)
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      enddo
!--fin de la lecture des resultats.

      write(clio3_out_id,'(A,I11,3A,I3,A)')
     &  'Exp '//cc6exp//', Iter', numit,
     &    ' , File ', ccfile, ' :',  nnvrd, ' variables read - OK.'
      if (nn99.eq.2) then
        write(clio3_out_id,'(A)') ccline(:nnc)
        write(clio3_out_id,
     &   '(A,I11,3A,I3,A)') 'Exp '//cc6exp//', Iter', numit,
     &    ' , File ', ccfile, ' :',  nnvrd, ' variables read - OK.'
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif

      if (refexp.eq.'      ') refexp = cc6exp

 500  continue
!--Fin de la lecture.
      close(ccfile3_id)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  7 ) Traitement des cas d'un changement de stockage / bathymetrie    |
!-----------------------------------------------------------------------

!     if(nnt.lt.0) then
!       do 710 j=1,jmax
!        do 710 i=1,imax
!         ub(i,j) = ub(i,j) * hu(i,j)
!         vb(i,j) = vb(i,j) * hu(i,j)
!710    continue
!     endif

!--Tester si la bathymetrie a ete lue :
      if (unsvol.gt.zero) then
!--precaution : vitesse nulle en dehors du domaine .
        do k=1,kmax
         do j=1,jmax
          do i=1,imax
            u(i,j,k) = u(i,j,k) * tmu(i,j,k)
            v(i,j,k) = v(i,j,k) * tmu(i,j,k)
          enddo
         enddo
        enddo
        do j=1,jmax
         do i=1,imax
           ub(i,j) = ub(i,j) * tmu(i,j,ks2)
           vb(i,j) = vb(i,j) * tmu(i,j,ks2)
         enddo
        enddo
      endif

!dmr [NOEQUI]
      q(:,:,:) = RESHAPE(vloc, (/ imax, jmax, kmax /))
!dmr [NOEQUI]
      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  9 ) Traitement des cas "problematiques" .                           |
!-----------------------------------------------------------------------

 910  continue
      write(clio3_out_id,'(A)') 'ARRET, redrunb : End of File !'
      goto 950

 920  continue
      write(clio3_out_id,'(A)') 'ARRET, redrunb : Read Error  !'

 950  continue
      write(clio3_out_id,'(A,I11,3A,I3,A)')
     &  'Exp '//cc6exp//', Iter', numit,
     &    ' , File ', ccfile, ' :',  nnvrd, ' variables read & STOP'
      write(clio3_out_id,'(A)') ccline(:nnc)

      if (nn99.eq.2) close(99)
      close(ccfile3_id)

      stop
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine redrunb -
      end
