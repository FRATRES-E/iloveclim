!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:58 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:58 CET 2009

      SUBROUTINE wrtbox(nbox, nboxmx, nnvusl, nnvvsl,
     &  sc0ref, flxbox, scabox, wmerid, scamdz, slpbox, vinbox, volbox,
     &  sfxbox, fxfact, ttfx, titbox, titref, titexp)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! appele par "meridflu" : Ecriture sur fichiers des bilans de boite.
!  modif : 23/03/99

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      
!      use newunit_clio_mod, only: clio3_out_id
      
!! END_OF_USE_SECTION

!- parametres locaux : <- en arguments
!     parameter ( nboxmx = 5 )

!--dummy variables :
      
      integer(kind=ip) :: nboxmx, ncount
      real(kind=dblp), dimension(0:nsmax)                   ::sc0ref
      real(kind=dblp), dimension(0:14,nboxmx,0:kmax,0:nsmax)::flxbox
      real(kind=dblp), dimension(nboxmx,0:kmax,nsmax)       ::scabox
     &               ,wmerid
      real(kind=dblp), dimension(nboxmx,kmax,nsmax)         ::scamdz
      real(kind=dblp), dimension(0:14,nboxmx,0:nsmax)       ::slpbox
      real(kind=dblp), dimension(nboxmx,0:kmax,0:nsmax)     ::vinbox
      real(kind=dblp), dimension(nboxmx,0:kmax)             ::volbox
      real(kind=dblp), dimension(nboxmx,0:nsmax)            ::sfxbox
      real(kind=dblp), dimension(0:nsmax)                   ::fxfact
      
c~       dimension sc0ref(0:nsmax)
c~       dimension flxbox(0:14,nboxmx,0:kmax,0:nsmax)
c~       dimension scabox(nboxmx,0:kmax,nsmax)
c~       dimension wmerid(nboxmx,0:kmax,nsmax)
c~       dimension scamdz(nboxmx,kmax,nsmax)
c~       dimension slpbox(0:14,nboxmx,0:nsmax)
c~       dimension vinbox(nboxmx,0:kmax,0:nsmax), volbox(nboxmx,0:kmax)
c~       dimension sfxbox(nboxmx,0:nsmax), fxfact(0:nsmax)
      
      character(len=8) :: ttfx(0:nsmax)
      character(len=20) :: titbox(nboxmx)
      character(len=27) :: titref(0:nsmax)
      character*(*) titexp

!- variables locales equivalentes :
!- variables locales :
      character(len=2) :: cc2, titm2c(0:4)
      character(len=10) :: titmur(0:4)
      character(len=20) :: racfil
      
      
!--- More locales
      real(kind=dblp) :: ccff
      integer(kind=ip):: imiddl, k, m, mm, mm1, mm2, n, nbil, nbox
     &                 , nfich, nm, nncc, nnvusl, nnvvsl, ns
      real(kind=dblp) :: salent, salsor, sorm1, sorm2, sormoy, sormoy1
     &                 , sormoy2, tement, temsor      

!-------Tableaux pour le calcul des Bilans sur les boites :
!-  Notes : Modification des noms de variables :
!-var_name < 6cc : ibond1/2 <- ibound1/2, scabox <- wscalm, scamdz <- wvscalm
!-  flxbox :
! 1er indice m : mod(m,3)= 1 => Entrant ; 2 => Sortant ; 0 => Entrant - Sortant
! 2eme indice nbil : numero de la boite
! 3eme indice k    : k > 0 = > niveau  ; = 0 => Somme sur tous les Niv.
! 4eme indice ns   : 0 <-> Masse ; entre 1 et nsmax <-> scalaire(ns)
!-- correspondance flxbox <-> w(v)nor/sud/ouest/est/tot :
!  avec : flxbox( 0/1/2 , nbil,k,-) <->   wtot(k,3/1/2,nbil)
!         flxbox( 3/4/5 , nbil,k,-) <->   wsud(k,3/1/2,nbil)
!         flxbox( 6/7/8 , nbil,k,-) <->   wnor(k,3/1/2,nbil)
!         flxbox( 9/10/11,nbil,k,-) <-> wouest(k,3/1/2,nbil)
!         flxbox(12/13/14,nbil,k,-) <->   west(k,3/1/2,nbil)
!   et    flxbox( 0/1/2 , nbil,0,-) <->  wvtot(3/1/2,nbil)
!         flxbox( 3/4/5 , nbil,0,-) <->  wvsud(3/1/2,nbil)
!         flxbox( 6/7/8 , nbil,0,-) <->  wvnor(3/1/2,nbil)
!         flxbox( 9/10/11,nbil,0,-) <-> wvwest(3/1/2,nbil)    <- wvouest
!         flxbox(12/13/14,nbil,0,-) <->  wvest(3/1/2,nbil)
! et aussi        wmerid(nbil,0,ns) <->  wvmeri(nbil,ns)
!-----

      data ncount / 0 /

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--initialisation :
      ncount = ncount + 1

!--definition de la racine des noms de fichiers :
      nncc = 0
      do n=1,len(refexp)
        if (refexp(n:n).ne.' ') nncc = n
      enddo
      write(racfil,'(A,I3,A)') refexp(:nncc)//'_', ncount, '.'
      if (ncount.lt.100) racfil(nncc+2:nncc+2) = '0'
      if (ncount.lt.10 ) racfil(nncc+3:nncc+3) = '0'
      nncc = nncc + 5

!- Ouverture et ecriture de la 1ere ligne :
      do nbil=1,nbox
         write(cc2,'(A1,I1)') 'b', nbil
!        open(60+nbil, file='box'//cc2//'.out', status='unknown')
         open(60+nbil, file=racfil(:nncc)//cc2//'.box',
     &        status='unknown')
         write(60+nbil,'(2A)') 'Fichier Resultat : ', titexp
         write(60+nbil,'(A,I3,2A)') ' Bilan sur la boite No=', nbil,
     &                              '  : ', titbox(nbil)
      enddo

!- titre pour chaque type de Flux : ttfx(0 -> nsmax)

!- titre pour chaque Mur :
      titmur(0) = 'Global ,  '
      titmur(1) = 'mur Sud,  '
      titmur(2) = 'mur Nord, '
      titmur(3) = 'mur West, '
      do nm=1,4
        titm2c(nm) = titmur(nm)(5:5)//':'
      enddo
        titm2c(0)  = titmur(0)(:1)//':'
!-

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  8 ) Ecriture des bilans de boites (Flux Masse & Scal.) sur fichiers |
!-----------------------------------------------------------------------

!-----
!--Ecriture sur fichiers des flux (Masse & Scal.) par Niveau.

      do nbil=1,nbox
       nfich=60+nbil
       write(nfich,'(A,1P9E14.6)')' scal_ref =',(sc0ref(ns),ns=1,nsmax)
!      write(clio3_out_id,*) "nbil,nbox,nfich",nbil,nbox,nfich
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--Flux de Masse (ns=0)
       ns = 0
       do nm=0,4
!---
        mm = nm*3
        write(nfich,*)
        write(nfich,'(5A)') titmur(nm), 'Flux ',
     &          ttfx(ns), ' : Net, Entree, Sortie :'
        do k=ks1,ks2
         write(nfich,1019) k, (flxbox(m,nbil,k,ns),m=mm,2+mm)
        enddo
         write(nfich,1018) (flxbox(m,nbil,0,ns),m=mm,2+mm)
!---
       enddo
!--fin de la boucle sur nm = "No_Mur".

!--Flux pour chaque scalaire :
       do ns=1,nsmax
        do nm=0,4
!---
         mm = nm*3
         mm1 = mm + 1
         mm2 = mm + 2
         write(nfich,*)
         write(nfich,'(5A)') titmur(nm), 'Flux ',
     &           ttfx(ns), ' : Net, Entree, Sortie', titref(ns)
         do k=ks1,ks2
          write(nfich,1019) k, (flxbox(m,nbil,k,ns),m=mm,mm2),
     &     (flxbox(m,nbil,k,ns)/max(flxbox(m,nbil,k,0),epsil),m=mm1,mm2)
         enddo
         write(nfich,1018) (flxbox(m,nbil,0,ns),m=mm,mm2),
     &    (flxbox(m,nbil,0,ns)/max(flxbox(m,nbil,0,0),epsil),m=mm1,mm2)
!---
        enddo
       enddo
!--fin des boucles sur nm = "No_Mur" et sur "ns".

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!---
       write(nfich,*)
       write(nfich,*) "Global : T_moy , dz*T_moy , Debit*T_moy"
       do k=ks1,ks2
         write(nfich,1019) k, scabox(nbil,k,1),
     &                        scamdz(nbil,k,1),wmerid(nbil,k,1)
       enddo
       write(nfich,'(E12.5)') wmerid(nbil,0,1)
!---
       write(nfich,*)
       write(nfich,*) "Global : S_moy , dz*S_moy , Debit*S_moy"
       do k=ks1,ks2
         write(nfich,1019) k,scabox(nbil,k,2),
     &                       scamdz(nbil,k,2),wmerid(nbil,k,2)
       enddo
       write(nfich,'(E12.5)') wmerid(nbil,0,2)
!---
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

 1018 format(2X,2(2X,3(1PE12.5,1X)))
 1019 format(I2,2(2X,3(1PE12.5,1X)))
 1105 format(A,1P5E13.5)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Ecriture sur fichiers, bilan sur les boites, flux Downsloping :
!-----

!     if (nbox.ge.1 .and.nvrl(nvrusl).ge.2 .and.nvrl(nvrvsl).ge.2) then
      if (nbox.ge.1 .and.nnvusl.ge.2 .and.nnvvsl.ge.2) then
!-----

      do nbil=1,nbox
!--Ecriture sur fichier, boite par boite :
        nfich = 60 + nbil
!- Masse :
          ns = 0
          write(nfich,*)
          write(nfich,'(5A)') 'Bilan Dwnslp.,  Flux ',
     &                    ttfx(ns), ' : Net, Entree, Sortie :'
        do nm=0,4
          mm = nm*3
          write(nfich,1105) titm2c(nm), (slpbox(m,nbil,ns),m=mm,mm+2)
        enddo
        do ns=1,nsmax
!- Scalaire "ns" :
          write(nfich,*)
          write(nfich,'(5A)') 'Bilan Dwnslp.,  Flux ', ttfx(ns),
     &     ' : Net, In, Out, Delta', titref(ns)(3:5),'(=Flx/Debit) :'
         do nm=0,4
          mm = nm*3
          write(nfich,1105) titm2c(nm), (slpbox(m,nbil,ns),m=mm,mm+2),
     &     (slpbox(m,nbil,ns)/max(slpbox(m,nbil,0),epsil),m=mm+1,mm+2)
         enddo
        enddo
!--fin de l'ecriture.
      enddo

!-----Fin des bilans sur les boites des flux de Downsloping .
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  9 ) Autres bilans de boites .
!-----------------------------------------------------------------------

      do nbil=1,nbox
!- Ecriture :
      nfich=60+nbil

!-----
!- Moyenne de la densite + des scalaires a l'interieur de la boite :
      write(nfich,*)
      write(nfich,'(A)')
     &     ' Volume (m^3)  & Moyenne de Rho, scal(ns=1,nsmax) :'
      ccff = dx * dy
      do k=ks1,ks2
        write(nfich,1021) k, volbox(nbil,k)*ccff,
     &                   (vinbox(nbil,k,ns),ns=0,nsmax)
      enddo
        write(nfich,1020) volbox(nbil,0)*ccff,
     &                   (vinbox(nbil,0,ns),ns=0,nsmax)
!- Moyenne des Flux de surface (Air/Mer) :
      write(nfich,'(2A)') 'Surf.Sup (m^2)',
     &  ' & Moyenne des Flux de Surf. [Unit=Flx(ns) <- class.par] :'
        write(nfich,1020) volbox(nbil,ks2)*unsdz(ks2)*ccff,
     &                (sfxbox(nbil,ns)*fxfact(ns),ns=0,nsmax)
 1020 format(2X,1PE12.5,1X,1P6E13.5)
 1021 format(I2,1PE12.5,1X,1P6E13.5)
!-----

      write(nfich,*)
      write(nfich,*) 'Echange total Entrant(=Sortant) (Sv) :',
     &                flxbox(1,nbil,0,0)*1.d-6
      sormoy = flxbox(0,nbil,0,1)*4.d6
      write(nfich,'(A,1PE12.5)')' Flux total de chaleur(W) : ',sormoy
      write(nfich,'(A,1PE12.5)')' Flux total de Sel        : ',
     &                          flxbox(0,nbil,0,2)
      sormoy1 = flxbox(1,nbil,0,1) / max(flxbox(1,nbil,0,0),epsil)
      sormoy2 = flxbox(2,nbil,0,1) / max(flxbox(2,nbil,0,0),epsil)
!--conversion K -> deg.C :
      imiddl = imax / 2
!     write(clio3_out_id,*) 'imiddl, jeq, ks2, scal :',
!    &            imiddl, jeq, ks2, scal(imiddl,jeq,ks2,1)
      if (scal(imiddl,jeq,ks2,1).gt.100.0) then
        sormoy1 = sormoy1 - 273.15
        sormoy2 = sormoy2 - 273.15
      endif
      write(nfich,*) 'Temp.(=Flx/Debit) In/Out :',sormoy1,sormoy2
      sormoy1 = flxbox(1,nbil,0,2) / max(flxbox(1,nbil,0,0),epsil)
      sormoy2 = flxbox(2,nbil,0,2) / max(flxbox(2,nbil,0,0),epsil)
      write(nfich,*) 'Sali.(=Flx/Debit) In/Out :',sormoy1,sormoy2

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Transport et Transferts(T,S) lies a la circulation meridienne :
      sorm1=0.0
      sorm2=0.0
      do k=ks1,ks2
         sorm1 = sorm1 + max(zero, flxbox(0,nbil,k,0))
         sorm2 = sorm2 + max(zero,-flxbox(0,nbil,k,0))
      enddo
      tement=0.0
      salent=0.0
      temsor=0.0
      salsor=0.0
      do k=ks1,ks2
         if (abs(flxbox(0,nbil,k,0)).gt.1e-2) then
           if (flxbox(0,nbil,k,0).gt.zero) then
            tement=tement+wmerid(nbil,k,1)/sorm1
            salent=salent+wmerid(nbil,k,2)/sorm1
           else
            temsor=temsor-wmerid(nbil,k,1)/sorm2
            salsor=salsor-wmerid(nbil,k,2)/sorm2
           endif
         endif
      enddo
!--conversion K -> deg.C :
      if (scal(imiddl,jeq,ks2,1).gt.100.0) then
        tement = tement - 273.15
        temsor = temsor - 273.15
      endif
      sorm1 = sorm1 * 1.d-6
      sorm2 = sorm2 * 1.d-6
      write(nfich,*)
      write(nfich,*) 'echange meridien In/Out (Sv) :',sorm1,sorm2
      write(nfich,'(A,1PE12.5)')' Flux(<= circ.Merid) de chaleur(W) : '
     &                         , wmerid(nbil,0,1)*4.d6
      write(nfich,'(A,1PE12.5)')' Flux(<= circ.Merid) de Sel        : '
     &                         , wmerid(nbil,0,2)
      write(nfich,*) 'temp.(=Flx/Debit) In/Out :',tement,temsor
      write(nfich,*) 'sali.(=Flx/Debit) In/Out :',salent,salsor

      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!- fermeture des fichiers :
      do nbil=1,nbox
        nfich = 60 + nbil
        close(nfich)
      enddo

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine wrtbox -
      end
