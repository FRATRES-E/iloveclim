!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:51 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:51 CET 2009

      SUBROUTINE redforc(nn99)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! mise en place des tableaux servant au forcage du modele (ocean uniquement)
! cas nn99=2 : impression de controle sur fichier "mouchard" ;
! Preparation du rappel en surface ou en profondeur.
! kforc -> quels fichiers de donnees a lire.
!  modif : 30/09/99

!! START_OF_USE_SECTION

      use const_mod, only: yeaday, zero

      use para0_mod, only: ijkmax, imax, ixjmax, jmax, kmax, nsmax
      use para_mod,  only:
      use bloc0_mod, only: js1, js2, ks1, ks2, tms, dts, spvr, scalr
     &  , phiss, rappes, unsdz, rappel, rapint
      use bloc_mod , only: kforc, mdforc, ijrap, nrap
      use reper_mod, only: yforc, rapp1, rapp0, kbath1, kbath2, filcor, unstyr
      use global_constants_mod, only: dblp=>dp, ip

      use newunit_clio_mod, only: clio3_out_id, mouchard_id

!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "reper.com"

!! END_OF_INCLUDE_SECTION

!--variables locales :

      character*50 ccfile
      character*70 line
      real(dblp), dimension(imax,jmax) :: zrap

      real(dblp)  :: ccmult, cctsr, cctsrk
      integer(ip) :: i, j, k, k2yf, kk1, kk2, kkerr, kknb, kkyf, n, n2yf
     &             , nnline, ns, nsrp, nn99
      real(dblp)  :: spvr2, xxmsq, xxns

      integer(ip) :: tsobs_om_id, ccfile_id, ccfile2_id, filcor_id
      
      line = ' '
      nnline = 1

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Forcage Constant (ex: moyenne annuelle) .                       |
!-----------------------------------------------------------------------

!---------------
! Options : kforc(for all Scalar) ;
! kk=mod(abs(kforc),100) -> Constant Field ; kforc/100 -> Time Dependant
!     kforc >= 0   : read file "tsobs.om" (spvr+T+S for all levels )
! 1+mod(kk-1,nsmax)  : Nb of Scalar to read ; and for each of theses :
! 1+(kk-1/nsmax) = k : Nb of level, from the surface to the bottom
!  File_name : "ts1up.om" (k=1), "ts2up.om" (k=2), ...
!-
! concerning each scalar (N), yforc(N) :
!   non.zero : read file "flux"N".om", and multiply surface flux by yflux.
!---------------


!--2.1 Lecture du forcage Dynamique :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--tension de vent en surface (Wind Stress) :
!ic0  open(unit=22,file='wsxy.om',status='old',form='UNFORMATTED')
!ic0  read (22, end=900) phisu
!ic0  read (22, end=900) phisv
!ic0  close(22)
!ic0  line = line(:nnline)//'wsxy.om '
!ic0  nnline = nnline + 8

!--2.2 Lecture des temperatures et salinites observees :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!- spvr (=Special Value) doit verifier, pour tout scalaire ns :
!  spvr < mini { Inf(Scal)ns - kmax , 2*Inf(Scal)ns - Sup(Scal)ns }
      spvr = -100.

      kkyf = mod(abs(kforc),100)
      if (kforc.ge.0) then
!- Lecture de spvr & kmax*2 niveaux :
        kknb = 2 * ijkmax
        open(newunit=tsobs_om_id,file='inputdata/clio/tsobs.om',status='old'
     & ,form='UNFORMATTED')
        read (tsobs_om_id, end=900) spvr
        call redtab( scalr(1,1,1,1), kknb, tsobs_om_id , kkerr)
        if (kkerr.ne.1) goto 900
        close(tsobs_om_id)
        line = line(:nnline)//'tsobs.om '
        nnline = nnline + 9
!- precaution : valeur speciale si rappel non nul :
        do ns=3,nsmax
         do k=1,kmax
          if (rapp1(k).ne.0.) then
            do j=1,jmax
             do i=1,imax
               scalr(i,j,k,ns) = spvr
             enddo
             enddo
            endif
        enddo
       enddo
      endif

      if (kkyf.eq.0 .or. kkyf.gt.kmax*nsmax) then
        k2yf = 0
        n2yf = 0
        kk1 = kmax
      else
!- n2yf = Nb de scalaire a lire ;
!- k2yf = Nb de niveaux(pour chaque Scal) a partir de la surface ;
!  et donc  kkyf = (k2yf-1)*nsmax + n2yf
        k2yf = 1 + (kkyf-1) / nsmax
        n2yf = 1 + mod(kkyf-1,nsmax)

!- Lecture de spvr & (k2yf*n2yf) niveaux :
        kk1 = kmax - k2yf + 1
        write(ccfile,'(A,I1,A)') 'ts', k2yf, 'up.om'

        open(newunit=ccfile_id, file=ccfile, status='old', form='UNFORMATTED')
        read (ccfile_id) spvr2
        if (kforc.ge.0 .and. spvr.ne.spvr2) then
          write(clio3_out_id,*) 'ARRET dans routine Redforc :'
          write(clio3_out_id,*) 'SPVR Differentes entre tsobs.om et '//ccfile
          close(ccfile_id)
          stop
        endif
        spvr = spvr2
        do ns=1,n2yf
         do k=kk1,kmax
          call redtab( scalr(1,1,k,ns), ixjmax, ccfile_id , kkerr)
          if (kkerr.ne.1) goto 900
         enddo
        enddo
        close(ccfile_id)
        line = line(:nnline)//ccfile
        nnline = nnline + 9
      endif

      if (kforc.ge.0) then
        kk1 = 1
      else
        kk2 = kmax - k2yf
        do ns=1,nsmax
         if (ns.gt.n2yf) kk2 = kmax
         do k=1,kk2
          do j=1,jmax
           do i=1,imax
             scalr(i,j,k,ns) = spvr
           enddo
          enddo
         enddo
        enddo
      endif

! correction celcius-kelvin
      do k=ks1,ks2
       do j=1,jmax
        do i=1,imax
          if (scalr(i,j,k,1).ne.spvr) scalr(i,j,k,1) =
     &                                scalr(i,j,k,1) + 273.15d0
        enddo
       enddo
      enddo

!--2.3 Mise en place des tableaux "rappes & rappel" (A Delta-t*unstyr pres) :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--Rappel explicite differencie T / S : Rappes .
        do ns=0,nsmax
         do j=1,jmax
          do i=1,imax
           rappes(i,j,ns) = rapp0(ns)
          enddo
         enddo
        enddo

!--Rappel implicite commun a T & S : Rappel .
      do k=kk1,kmax
       do j=1,jmax
        do i=1,imax
          if(j.eq.js1.and.kbath1(i).ge.(ks2+ks1-k)) then
!- Rappel sur le mur S :
            rappel(i,j,k) = max(rapp1(k), rapp0(nsmax+2))
          elseif(j.eq.js2.and.kbath2(i).ge.(ks2+ks1-k)) then
!- Rappel sur le mur N :
            rappel(i,j,k) = max(rapp1(k), rapp0(nsmax+1))
          else
!- Rappel a l'interieur du bassin :
            rappel(i,j,k) = rapp1(k)
          endif
        enddo
       enddo
      enddo

!--2.4 Flux en surface pour chaque scalaire :
!------------------------------------------------------------------

!- Convention Flux : + vers le Haut ; unites : Scal() x L(m) / Temps(s)
!-  dans modele : phiss = (Flux->haut) x (Delta_T) / Dz(1er_niv)

      do ns=1,nsmax
        if (yforc(ns).ne.zero) then
!- lecture fichier flux :
          write(ccfile,'(A,I1,A)') 'flux', ns, '.om'
          open(newunit=ccfile2_id,file=ccfile,status='old',form='UNFORMATTED')
          call redtab( phiss(1,1,ns), ixjmax, ccfile2_id , kkerr)
          close(ccfile2_id)
          line = line(:nnline)//ccfile
          nnline = nnline + 9
        endif
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Forcage Saisonnier / Mensuel .                                  |
!-----------------------------------------------------------------------

!--LECTURE des donnes saisoniers de T, S et tension du vent
!--tension de vent mensuelle en surface (Wind Stress) :
!sai  open(unit=23,file='wsxy.mens.3x3.om',
!sai &    status='old',form='UNFORMATTED')
!sai   read(23) txmens
!sai   read(23) tymens
!sai  close(23)

!--temperature mensuelle et salinite saisoniere observees :
!sai  open(unit=25,file='tlev.mens.3x3.om',
!sai &    status='old',form='UNFORMATTED')
!sai   read(25) tmens
!sai  close(25)

!sai  open(unit=26,file='slev.sais.3x3.om',
!sai &    status='old',form='UNFORMATTED')
!sai   read(26) smens
!sai  close(26)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  7 ) Modification du Forcage, Conversion pour utilisation directe .  |
!-----------------------------------------------------------------------

!- modification localisee par appel a "correct" :
      if (mdforc.eq.1) call correct(nn99,filcor)

!- modification du rappel explicite par la lecture de filcor
      if (mdforc.eq.2) then
        write(clio3_out_id,*) 'Restoring not applied in some regions'
        write(clio3_out_id,*) 'file = ', filcor
        open(newunit = filcor_id, file = filcor, status = 'old')
        read(filcor_id,*) zrap
        do ns=0,nsmax
         do j=1,jmax
          do i=1,imax
           rappes(i,j,ns)=rappes(i,j,ns)*zrap(i,j)
          enddo
         enddo
        enddo
        write(clio3_out_id,*) 'restoring 10,10 = ',
     &    rappes(10,10,0),rappes(10,10,1) ,rappes(10,10,2)
       endif

!- prepare le calcul du Flux derive de l'Advec.H.Obs :
!adh  nflag = 1
!adh  call initflx(nflag, nn99)

!--Mise en place definitive des tableaux "rappel & rappes" (= Delta-t / tau) :
      cctsr = unstyr / (yeaday * 86400.)
      do k=1,kmax
        cctsrk = dts(k) * cctsr
        do j=1,jmax
         do i=1,imax
          rappel(i,j,k) = tms(i,j,k) * cctsrk * rappel(i,j,k)
         enddo
        enddo
      enddo
      cctsrk = dts(ks2) * cctsr
      do ns=0,nsmax
       do j=1,jmax
        do i=1,imax
         rappes(i,j,ns) = cctsrk * rappes(i,j,ns)
        enddo
       enddo
      enddo

!--Mise en place du tableau "rapint" (=Delta-t/tau) ; si pas d'Obs -> rapint=0
      do ns=1,nsmax
       do k=1,kmax
        cctsrk = dts(k) * cctsr
        do n=1,nrap(k,ns)
         i = 1 + mod(ijrap(n,k,ns),imax)
         j = 1 + ijrap(n,k,ns)/imax
         rapint(n,k,ns) = cctsrk * rapint(n,k,ns)
     &                  * min(tms(i,j,k), (scalr(i,j,k,ns)-spvr))
 
        enddo
       enddo
      enddo

!--Terre ou Pas d'Obs -> Rappes = 0
      if (nn99.eq.2) then
       write(mouchard_id,'(A)') 'Points sans Observations en surface :'
       do j=1,jmax
        do i=1,imax
          xxmsq = nsmax * tms(i,j,ks2)
          do ns=1,nsmax
            xxns = min(tms(i,j,ks2), (scalr(i,j,ks2,ns)-spvr))
            rappes(i,j,ns) = rappes(i,j,ns) * xxns
            xxmsq = xxmsq - xxns
          enddo
          xxns = min(tms(i,j,ks2), (scalr(i,j,ks2,2)-spvr))
          rappes(i,j,0) = rappes(i,j,0) * xxns
!dmr debug          if (xxmsq.gt.epsil) write(99,'(2I4,3x,1P10E13.5)')
!dmr debug     &        i,j,(scalr(i,j,ks2,ns),ns=1,nsmax)
        enddo
       enddo
      else
        do ns=0,nsmax
         nsrp = max(ns,2-ns)
         do j=1,jmax
          do i=1,imax
           rappes(i,j,ns) = rappes(i,j,ns)
     &              * min(tms(i,j,ks2), (scalr(i,j,ks2,nsrp)-spvr))
          enddo
         enddo
        enddo
      endif
!--Pas d'Obs -> Rappel = 0
      do ns=1,nsmax
       do k=1,kmax
        do j=1,jmax
         do i=1,imax
          rappel(i,j,k) = min(rappel(i,j,k), (scalr(i,j,k,ns)-spvr))
         enddo
        enddo
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--conversion & changements d'unites (mais pas de signe !), limiteurs de flux :
!- unitfx(0) = Year : Fx. en m/y ; ATTENTION : phiss(0) en m et phimnx en m/s.
!- unitfx(1) = rho.Cp : Fx. en W/m2 ; unitfx(2) = Year / 34.7 g/l : Fx. en m/y
!- conversion tableau phiss <- utilise directement ds equation des scalaires.
      do ns=0,nsmax
        if (yforc(ns).ne.zero) then
          ccmult = yforc(ns) * dts(ks2) * unsdz(ks2)
          if (ns.eq.0) ccmult = yforc(ns) * dts(ks2)
          do j=1,jmax
           do i=1,imax
            phiss(i,j,ns) = ccmult * phiss(i,j,ns)
           enddo
          enddo
        endif
!ic0      ccmult = dts(ks2) * unsdz(ks2) / abs(unitfx(ns))
!ic0      if (ns.eq.0) ccmult = 1.0 / abs(unitfx(ns))
!ic0      do j=1,jmax
!ic0       do i=1,imax
!ic0        phimnx(i,j,0,ns) = ccmult * phimnx(i,j,0,ns)
!ic0        phimnx(i,j,1,ns) = ccmult * phimnx(i,j,1,ns)
!ic0       enddo
!ic0      enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  8 ) Impression de controle sur le fichier "mouchard" .              |
!-----------------------------------------------------------------------

      if (nn99 == 2) then
!--ecriture de controle :
        write(99,*) 'DeltaT(Surf)/tau , Coeff rappel Expl. FW, T, S,'
     &              //' rappel N, S :'
        write(mouchard_id,*) cctsrk, (rapp0(k),k=0,nsmax+2)
        write(mouchard_id,*) 'coeff rappel de 1 a kmax :'
        write(mouchard_id,*) rapp1
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  9 ) Sortie de la routine .
!-----------------------------------------------------------------------

      write(clio3_out_id,'(A)') 'Files read :'//line(:nnline)

      return

 900  continue                      ! Error handling for READ errors from various files.

      write(clio3_out_id,*) 'Arret routine "redforc" :'
      write(clio3_out_id,*) 'Probleme de lecture apres les fichiers suivant :'
      write(clio3_out_id,*) line(:nnline)

      stop
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine redforc -
      end
