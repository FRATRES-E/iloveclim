!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:45 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:45 CET 2009

      SUBROUTINE etat(mixage, nn99)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--CALCUL DE LA POUSSEE ET CALCUL DE N2 (=bvf) -- Ajust.Conv. selon "mixage"
!   b = gpes * (rho - rhoref) / rho0 ; rhoref = rho(Tref,Sref,z)
!     rho0 = 1030 kg/m3 ; gpes = 9.80 (common static)
!     rho est donne par la formule de C.Eckart (unite : kg/l)
!        (Amerrican Journal of Science, Vol256, April 58, pp225-240)
! Teste de Comparaison : utilise pour g & rho0 & P(z) :
! rho0 = 1.03 kg/l | gpes = 9.8 | P(atm) = z * gpes * rho0 / 101325 .
! -- Calcul de N2 (=Br.Vais.Freq)(bvf) suivant la formule numerique :
! BVF(k-1/2) = [b(T(k-1),S(k-1),z(k-1/2)) - b(T(k),S(k),z(k-1/2))] / dz(k-1/2)
! -- Integration de rho.g.dz par boite pesante.
!  effectue les raccords de grille pour les tableaux w, scal, q et bvf
! -- mixage :  dirige Ajust.Conv. : 0 : rien ;
!  1 : Suppr. TOUTES les instab. ; 2 : melange par paires avec lstab passages
!  3 : Permute et melange partiel (=ajcmix)
! Tient compte de l'Accelerateur de Convergence dans l'Ajust.Convectif.
!  modif : 22/07/97  --  Version Type 7  / Celcius .  --

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip      

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use densit_mod
      use newunit_clio_mod, only: clio3_out_id, mouchard_id

!! END_OF_USE_SECTION
      integer(kind=ip)  :: mixage
      integer(kind=ip) :: nn99

!--variables locales :
      real(kind=dblp), dimension(imax)  :: qintf
      real(kind=dblp), dimension(kmax)  :: rhoref
      real(kind=dblp), dimension(imax,kmax)  :: bup, bdw, cfb1, cfb2

!- pour Ajust. Conv. (tous, 1, 2, 3) :
      real(kind=dblp), dimension(nsmax) :: scalmi
      real(kind=dblp), dimension(kmax)  :: bbajc, mixtab, ffiajc
      real(kind=dblp), dimension(imax,kmax)  :: ffajc

!--Constantes de l'equation d'etat :
c~       dimension eckart(0:10)
      real(kind=dblp), dimension(0:10)  :: eckart
      real(kind=dblp)  :: tkc, tref, sref, tavr, savr, cksi
      
      data eckart /
     &             0.698d0  , 5890.0d0 ,   38.0d0 , 0.375d0 ,  3.0d0  ,
     &             1779.5d0 ,  11.25d0 , 0.0745d0 ,   3.8d0 , 0.01d0  ,
     &           101325.0d0 /

      data tkc / 273.15d0 /
!ic0  data tkc / 0.d0 /

      data tref, sref / 4.d0 , 34.0d0 /

!--Pour l'ecriture en Temp.Potentielle :
      data tavr, savr, cksi / 2.d0 , 34.7d0 , 0.9d-4 /
!- cksi = (T - teta) / |z|, evalue pour S= 35, T=2 , |z|=5km, formule du Gill.



!--even more locales ...
      integer(kind=ip):: n, k, k2, k1, kcrois, j, i, kloc, kins, kdwmix
     &                    , kupmix, ktopmi, ns, kbotmi, nstab, kk2, kk1
     &                    , kk, kkd, nk, kdwm
     
      real(kind=dblp) :: gamma, ccpotp, unscpp, ccpotr, ccpz, ccb1, ccb2
     &                    , ccajc, zmix, tloc, sloc, unshmi, bdwmix
     &                    , bbsav, dqd2, tmix, smix, ajc0mx, bdwins
     &                    , ajchmx, sscdwm, ddssc


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Mise en place des Coeffs. a la 1ere Iter :                      |
!-----------------------------------------------------------------------

      if (numit.le.nstart) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!       write(clio3_out_id,*) 'etat : version 7  K , I-Boite '
        write(clio3_out_id,*) 'etat : version 7  C , I-Boite '
     &          //', (lstab,mixage) = ', lstab, mixage

!---------
!--1.1 Equation d'etat en Temperature Potentielle :
      do n=0,10
        cstrho(n) = eckart(n)
      enddo
      gamma = cksi * eckart(10) / ( gpes * rho0 )
      ccpotp = ( eckart(2) - 2.0*eckart(3)*tavr ) * gamma
      unscpp = 1.0 / (1.0 + ccpotp)
      ccpotr = ( eckart(6) - 2.0*eckart(7)*tavr - eckart(9)*savr )
     &         * gamma * unscpp
      cstrho(0)  = cstrho(0) + ccpotr
      cstrho(5)  = cstrho(5) - ccpotr * cstrho(1)
      cstrho(6)  = cstrho(6) - ccpotr * cstrho(2)
      cstrho(7)  = cstrho(7) - ccpotr * cstrho(3)
      cstrho(8)  = cstrho(8) + ccpotr * cstrho(4)
      cstrho(10) = cstrho(10) * unscpp

!---------
!--1.2 Mise en place des coefficients dependant de la profondeur.

!- pour l'Equ. d'etat (cfb*k*) :
      gravit = gpes * 1000. / rho0
      ccpz = gpes * rho0 / cstrho(10)
      cstrho(1) = cstrho(1) + 1.
      do k=1,kmax
        cfb1z(k)  = cstrho(1) - ccpz * z(k)
      enddo
      do k=1,kmax+1
        cfb1z4(k) = cstrho(1) - ccpz * zw(k)
      enddo

!---------
!--mise en place de la "densite" de reference ; coeff bilan Ajust.Conv. :
          ccb1 = cstrho(4)*sref
     &         + (cstrho(2)-cstrho(3)*tref)*tref
          ccb2 = cstrho(5)
     &         + (cstrho(6) - cstrho(7)*tref)*tref
     &         - (cstrho(8) + cstrho(9)*tref)*sref
!     ccajc = rho0 / ( DFLOAT(abs(ninfo)) * dts(kmax) )
      ccajc = rho0 / dts(kmax)
      zwtau(1+kmax) = zw(1+kmax)
      do k=kmax,1,-1
        rhoref(k) = 1./(cstrho(0)+ ccb2/(ccb1+cfb1z(k)) )
        bref(k) = gravit * rhoref(k)
        dztau(k) = dz(k) * dts(ks2) / dts(k)
        zwtau(k) = zwtau(k+1) - dztau(k)
        rho0dz(k) = ccajc * dztau(k)
      enddo

!- pour le mixage par paire (m2) :
      do k=2,kmax
        cfm2up(k) = dztau(k)   / (dztau(k-1) + dztau(k))
        cfm2dw(k) = dztau(k-1) / (dztau(k-1) + dztau(k))
        zmix = 0.5 * (zw(k-1) + zw(k+1))
        rhozdz(0,k) = (z(k)   - zmix) * rho0dz(k)
        rhozdz(1,k) = (z(k-1) - zmix) * rho0dz(k-1)
      enddo

!- pour l'Ajust.Conv. par permutation :
      do k2=1,kmax
       do k1=1,kmax
!       dzsdz(k1,k2) = (dz(k1) * unsdz(k2)) * (dts(k2) / dts(k1))
        dzsdz(k1,k2) = dztau(k1) / dztau(k2)
       enddo
      enddo
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      if (nn99.eq.2) then
!--Ecriture de controle sur "mouchard" :
        write(mouchard_id,'(A,3F18.14)') 'tavr,savr, cksi :', tavr, savr, cksi
        write(mouchard_id,'(A,3F18.14)') 'ccpotp,ccpotr :', ccpotp, ccpotr
        write(mouchard_id,*) 'Coeff Eq.d''Etat : Eckart | en Temp.Pot :'
        write(mouchard_id,'(1P2E17.10)') (eckart(k),cstrho(k),k=0,10)
        write(mouchard_id,'(A,3F18.14)') 'tref,sref, rhoref(tref,sref,k) :',
     &                           tref, sref
        write(mouchard_id,'(1P5E16.9)') (rhoref(k),k=1,kmax)
        write(mouchard_id,*)
        write(mouchard_id,*) 'Coeff dztau(k)=dz(k)*dt(surf)/dt(k) :'
        write(mouchard_id,'(5F16.9)') (dztau(k),k=1,kmax)
        write(mouchard_id,*)
      endif

      kcrois = 1
      do k=2,kmax
        if (dztau(k).gt.dztau(k-1)) kcrois = k
      enddo
      if (kcrois.ne.1.and.lstab.eq.-3) then
        write(clio3_out_id,
     &   '(A,I3,A)') 'Ajusclio/sources/t.Conv(-3=lstab) & Acc.Conv(k=',
     &                      kcrois, ' ) INCOMPATIBLE !'
        write(clio3_out_id,'(A)') ' => ARRET , routine "etat" '
        stop
      endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      endif
clio/sources/
!--Initialisation :
      ninstb = 0

!--Traitement differencie selon le type d'ajustement convectif :

      if (mixage.eq.1) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Suppression de TOUTES les instabilites (lstab = -1) :           |
!-----------------------------------------------------------------------

!--convention : mixtab(k) : plus bas niveau melange avec la boite "k"
!               kdwmix : plus bas  niveau a melanger avec la boite courante.
!               kupmix : plus haut niveau a melanger avec la boite courante.
!               ktopmi : plus haut niveau melange de toute la colonne

!--Debut de la boucle externe sur l'indice de latitude j :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,k,kins,ns,kloc,tloc,sloc,ccb1,ccb2,bup,bdw)
!$DIR SPP LOOP_PRIVATE(mixtab,kdwmix,kupmix,ktopmi,zmix,unshmi)
!$DIR SPP LOOP_PRIVATE(scalmi,bdwmix,bbsav,bbajc,kbotmi,qintf,dqd2)
      jloop: do j=js1,js2
!-----

!--2.1 Calcul de b(i,j,k) dans tout le domaine :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      do k=ks1,ks2
        do i=is1(j),is2(j)
!- conversion Degr.Celcius et tranfert -> variables locales (tloc,sloc) :
          kloc = max(k,kfs(i,j))
          tloc = scal(i,j,kloc,1) - tkc
          sloc = scal(i,j,kloc,2)
          ccb1 = cstrho(4)*sloc
     &         + (cstrho(2)-cstrho(3)*tloc)*tloc
          ccb2 = cstrho(5)
     &         + (cstrho(6) - cstrho(7)*tloc)*tloc
     &         - (cstrho(8) + cstrho(9)*tloc)*sloc
          b(i,j,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) ) - bref(k)
          bup(i,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z4(k+1)) )
          bdw(i,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z4( k )) )
       enddo
      enddo

      iloop: do i=is1(j),is2(j)
!- boucle NON VECTORISABLE sur l'indice de longitude "i" :

!--2.2 Detecte la 1ere instabilite :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        do kins=kfs(i,j)+1,ks2
!- teste la stabilite de kins / kins-1 :
          if (bdw(i,kins).gt.bup(i,kins-1)) then
!- initialisation :
            do k=ks1,ks2
              mixtab(k) = k
            enddo
            kdwmix = kins - 1
            mixtab(kins) = kdwmix
            kupmix = kins
            ktopmi = kins
            fqajc(i,j,kins) = fqajc(i,j,kins) + 1.
            goto 230    ! exit loop on kins, but continue loop in i (iloop).
                        ! i.e., skip the cycle command right after the loop.
          endif
        enddo
        cycle iloop

 230    continue        ! required as a landing point from exit from previous loop
                        ! and also from a loop deeper down

!--2.3 mixing des boites depuis kdwmix (=mixtab(kins)) --> kupmix
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        zmix   = 0.5 * ( zw(kupmix+1) + zw(kdwmix) )
        unshmi = 1.0 / ( zwtau(kupmix+1) - zwtau(kdwmix) )
        do ns=1,nsmax
          scalmi(ns) = 0.
          do k=kdwmix,kupmix
            scalmi(ns) = scalmi(ns) + dztau(k)*scal(i,j,k,ns)
          enddo
          scalmi(ns) = scalmi(ns) * unshmi
          do k=kdwmix,kupmix
            scal(i,j,k,ns) = scalmi(ns)
          enddo
        enddo
!- fin du mixing .

!--2.4 calcule les nouveaux "b" et fait le bilan de l'ajustement convectif :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
          tloc = scal(i,j,kins,1) - tkc
          sloc = scal(i,j,kins,2)
          ccb1 = cstrho(4)*sloc
     &         + (cstrho(2)-cstrho(3)*tloc)*tloc
          ccb2 = cstrho(5)
     &         + (cstrho(6) - cstrho(7)*tloc)*tloc
     &         - (cstrho(8) + cstrho(9)*tloc)*sloc
        bdwmix = gravit / (cstrho(0) + ccb2/(ccb1+cfb1z4(kdwmix)) )
        do k=kdwmix,kupmix
          bbsav = b(i,j,k)
          bup(i,k) = bdwmix
          bdw(i,k) = bdwmix
          b(i,j,k) = gravit /
     &               (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) ) - bref(k)
          bbajc(k) = (bbsav-b(i,j,k)) * (z(k)-zmix) * rho0dz(k)
          fqajc(i,j,1) = fqajc(i,j,1) + bbajc(k)
        enddo
        bup(i,kupmix) = gravit /
     &                  (cstrho(0) + ccb2/(ccb1+cfb1z4(kupmix+1)) )
        if ( kdwmix.eq.kfs(i,j) .and. kdwmix.gt.ks1 ) then
          do k=ks1,kdwmix-1
            b(i,j,k) = gravit /
     &               (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) ) - bref(k)
          enddo
        endif
!       if (kupmix.ge.(ks2-1)) hmajc(i,j) = max(hmajc(i,j),-z(kdwmix))
!       ninstb = ninstb + kupmix - kdwmix

!--2.5 Test de stab. - determine kdwmix & kupmix :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        kbotmi = max(kfs(i,j)+1,kdwmix)
        do kins=kbotmi,ks2
!- teste la stabilite de kins / kins-1 :
          if (mixtab(kins).eq.kins .and.
     &        bdw(i,kins).gt.bup(i,kins-1)) then
            kdwmix = mixtab(kins-1)
            mixtab(kins) = kdwmix
            kupmix = kins
            fqajc(i,j,kins) = fqajc(i,j,kins) + 1.
!- Jonction avec une zone (initialement instable) deja melangee :
            do k=kins+1,ktopmi
              if (mixtab(k).eq.kins) then
                kupmix = k
                mixtab(k) = kdwmix
              endif
            enddo
            ktopmi = kupmix
            goto 230    ! exit loop on kins, and return above before stage 2.3
          endif
       enddo

!--fin de la boucle sur "i" .
      enddo iloop

!--Fin du bloc 'suppression de TOUTES les instabilites'.

!--2.6 Calcul de la pression reduite q :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--calcul de la pression (q) (transfere depuis 'uve')
!-   (integration de rho.g par la Methode des Trapezes) :
!
!     do i=is1(j),is2(j)
!       q(i,j,ks2)= 0.5*dz(ks2)*b(i,j,ks2)
!     enddo
!
!     do k=ks2-1,ks1,-1
!      do i=is1(j),is2(j)
!       q(i,j,k)=q(i,j,k+1)+0.5*dz4(k+1)*(b(i,j,k+1)+b(i,j,k))
!      enddo
!     enddo
!
!---------------------
!--calcul de la pression (q) , integration de rho.g "par boite pesante" :

      do i=1,imax
        qintf(i) = 0.0
      enddo

      do k=ks2,ks1,-1
       do i=is1(j),is2(j)
         dqd2 = 0.5*dz(k)*b(i,j,k)
         q(i,j,k) = qintf(i) + dqd2
         qintf(i) = q(i,j,k) + dqd2
       enddo
      enddo


!--2.7 calcul de N2 (transfere depuis 'stab') :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      do k=ks1+1,ks2
       do i=is1(j),is2(j)
        bvf(i,j,k) = tms(i,j,k-1) * unsdzw(k) * (bup(i,k-1)-bdw(i,k))
       enddo
      enddo

!--Fin de la boucle externe sur l'indice de latitude j .
      enddo jloop

      elseif(mixage.eq.2) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Suppression des instabilites : melange par paires (lstab > 0) . |
!-----------------------------------------------------------------------

!--Debut de la boucle externe sur l'indice de latitude j :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,k,kloc,tloc,sloc,ccb1,ccb2,bup,bdw)
!$DIR SPP LOOP_PRIVATE(ns,nstab,kk,nk,kk1,kk2,kkd)
!$DIR SPP LOOP_PRIVATE(scalmi,tmix,smix,bbajc,ffajc,qintf,dqd2)
      sloop: do j=js1,js2
!-----

!--3.1 Calcul de b(i,j,k) dans tout le domaine :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      do k=ks1,ks2
        do i=is1(j),is2(j)
!- conversion Degr.Celcius et tranfert -> variables locales (tloc,sloc) :
          kloc = max(k,kfs(i,j))
          tloc = scal(i,j,kloc,1) - tkc
          sloc = scal(i,j,kloc,2)
          ccb1 = cstrho(4)*sloc
     &              + (cstrho(2)-cstrho(3)*tloc)*tloc
          ccb2 = cstrho(5)
     &              + (cstrho(6) - cstrho(7)*tloc)*tloc
     &              - (cstrho(8) + cstrho(9)*tloc)*sloc
          b(i,j,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) ) - bref(k)
          bup(i,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z4(k+1)) )
          bdw(i,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z4( k )) )
          ffajc(i,k) = 1.
       enddo
      enddo

!--3.2 Debut du traitement de l'Ajust.Convectif : boites a melanger ?
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!- le nombre de passages = lstab
      nloop: do nstab=1,lstab
       kk2 = ks2 - mod(numit+nstab,2)
!--debut boucle NON VECTORISABLE sur l'indice de longitude "i" :
       longi: do i=is1(j),is2(j)

        kk1 = kfs(i,j) + 1
        do kk=kk2,kk1,-2
!--3.2 Detecte les boites a melanger : teste la stabilite de kk / kk-1 :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
          if (bdw(i,kk).gt.bup(i,kk-1)) then
!--3.3 mixing des boites   kk / kk-1 :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
            kkd = kk - 1
            do ns=1,nsmax
              scalmi(ns) = cfm2up(kk) * scal(i,j,kk,ns)
     &                   + cfm2dw(kk) * scal(i,j,kkd,ns)
            enddo
!--3.4 Mise en place des nouveaux scalaires et calcule des nouveaux "b" :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
            tmix = scalmi(1) - tkc
            smix = scalmi(2)
            ccb1 = cstrho(4)*smix
     &           + (cstrho(2)-cstrho(3)*tmix)*tmix
            ccb2 = cstrho(5)
     &           + (cstrho(6) - cstrho(7)*tmix)*tmix
     &           - (cstrho(8) + cstrho(9)*tmix)*smix
            do nk=0,1
             k = kk - nk
             do ns=1,nsmax
               scal(i,j,k,ns) = scalmi(ns)
             enddo
             bbajc(k) = b(i,j,k)
             bdw(i,k) = gravit / (cstrho(0) + ccb2/(ccb1+cfb1z4(k)))
             b(i,j,k) = gravit / (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) )
     &                 - bref(k)
             bbajc(k) = rhozdz(nk,kk) * (bbajc(k) - b(i,j,k))
            enddo
            bup(i,kkd) =  bdw(i,kk)
            bup(i,kk) = gravit / (cstrho(0) + ccb2/(ccb1+cfb1z4(kk+1)))
            if ( kkd.eq.kfs(i,j) .and. kkd.gt.ks1) then
              do k=ks1,kkd-1
               b(i,j,k) = gravit / (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) )
     &                  - bref(k)
              enddo
            endif
!--3.5 bilan de l'ajustement convectif :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
            fqajc(i,j,1) = fqajc(i,j,1) + bbajc(kkd) + bbajc(kk)
            fqajc(i,j,kk) = fqajc(i,j,kk) + ffajc(i,kk)
            ffajc(i,kk) = 0.
!           hmajc(i,j) = max(hmajc(i,j),-z(kk))
!           if ( hmajc(i,j).ge.-zw(kk+1) )
!    &           hmajc(i,j) = max(hmajc(i,j),-z(kkd))
!           ninstb = ninstb + 1
          endif
!- fin du mixing .
        enddo

       enddo longi

      enddo nloop


!--Fin du bloc 'suppression des instabilites par paires'.

!--3.6 Calcul de la pression reduite q :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--calcul de la pression (q) (transfere depuis 'uve')
!-   (integration de rho.g par la Methode des Trapezes) :
!
!     do i=is1(j),is2(j)
!       q(i,j,ks2)= 0.5*dz(ks2)*b(i,j,ks2)
!     enddo
!
!     do k=ks2-1,ks1,-1
!      do i=is1(j),is2(j)
!       q(i,j,k)=q(i,j,k+1)+0.5*dz4(k+1)*(b(i,j,k+1)+b(i,j,k))
!      enddo
!     enddo
!
!---------------------
!--calcul de la pression (q) , integration de rho.g "par boite pesante" :

      do i=1,imax
        qintf(i) = 0.0
      enddo

      do k=ks2,ks1,-1
       do i=is1(j),is2(j)
        dqd2 = 0.5*dz(k)*b(i,j,k)
        q(i,j,k) = qintf(i) + dqd2
        qintf(i) = q(i,j,k) + dqd2
       enddo
      enddo


!--3.7 calcul de N2 (transfere depuis 'stab') :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      do k=ks1+1,ks2
       do i=is1(j),is2(j)
        bvf(i,j,k) = tms(i,j,k-1) * unsdzw(k) * (bup(i,k-1)-bdw(i,k))
       enddo
      enddo

!--Fin de la boucle externe sur l'indice de latitude j .
      enddo sloop

      elseif (mixage.eq.3) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Suppression des instabilites par permutation (lstab = -3) :     |
!-----------------------------------------------------------------------

!--convention : kins : niveau (= interface) instable.
!               kdwm : plus bas atteint par l'Ajust. Conv.
!-----
      ajc0mx = 1.0 - ajcmix

!--Debut de la boucle externe sur l'indice de latitude j :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,k,kins,ns,kloc,tloc,sloc,cfb1,cfb2,bup,bdw)
!$DIR SPP LOOP_PRIVATE(kdwm,bdwins,ajchmx,sscdwm,ddssc,scalmi)
!$DIR SPP LOOP_PRIVATE(zmix,bbsav,bbajc,ffiajc,qintf,dqd2)
      jloop2: do j=js1,js2
!-----

!--4.1 Calcul de b(i,j,k) dans tout le domaine :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      do k=ks1,ks2
        do i=is1(j),is2(j)
!- conversion Degr.Celcius et tranfert -> variables locales (tloc,sloc) :
          kloc = max(k,kfs(i,j))
          tloc = scal(i,j,kloc,1) - tkc
          sloc = scal(i,j,kloc,2)
          cfb1(i,k) = cstrho(4)*sloc
     &              + (cstrho(2)-cstrho(3)*tloc)*tloc
          cfb2(i,k) = cstrho(5)
     &              + (cstrho(6) - cstrho(7)*tloc)*tloc
     &              - (cstrho(8) + cstrho(9)*tloc)*sloc
          b(i,j,k) = gravit /
     &      (cstrho(0)+ cfb2(i,k)/(cfb1(i,k)+cfb1z(k)) ) - bref(k)
          bup(i,k) = gravit /
     &      (cstrho(0)+ cfb2(i,k)/(cfb1(i,k)+cfb1z4(k+1)) )
          bdw(i,k) = gravit /
     &      (cstrho(0)+ cfb2(i,k)/(cfb1(i,k)+cfb1z4( k )) )
        enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

       longi2: do i=is1(j),is2(j)
!--debut boucle NON VECTORISABLE sur l'indice de longitude "i" :

        kinsloop: do kins=kfs(i,j)+1,ks2
          ffiajc(kins) = 1.
          if (bdw(i,kins).gt.bup(i,kins-1))  then
!- debut du traitemnent de l'instabilite de kins / kins-1 :

!--bilan de l'Ajust.Conv : Interface (k/k-1) concernees par le melange :
          fqajc(i,j,kins) = fqajc(i,j,kins) + 1.
          ffiajc(kins) = 0.

!--4.2 repere les boites a permuter et a melanger partiellement  :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Recherche de "kdwm" tel que bup(kdwm-1) > b(Tkins,Skins,zw(kdwm))
         kdwm = kfs(i,j)
         do k=kins-1,kfs(i,j)+1,-1
          bdwins = gravit /
     &       (cstrho(0)+ cfb2(i,kins)/(cfb1(i,kins)+cfb1z4(k)) )
          if (bdwins.le.bup(i,k-1) ) then
            kdwm = k
            exit
          endif
          fqajc(i,j,k) = fqajc(i,j,k) + ffiajc(k)
          ffiajc(k) = 0.
         enddo

!--4.3 Permutation circulaire de H=dztau(kins) + melange kdwm --> kins :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        ajchmx = ajcmix / ( zwtau(kins+1) - zwtau(kdwm) )
        do ns=1,nsmax
!- Permute & calcule Moy. & New Value :
          sscdwm = scal(i,j,kins,ns)
          scalmi(ns) = dztau(kins) * sscdwm
          scal(i,j,kins,ns) = scal(i,j,kins-1,ns)
          do k=kdwm,kins-1
            ddssc = dzsdz(kins,k) * (sscdwm - scal(i,j,k,ns))
            sscdwm = scal(i,j,k,ns)
            scalmi(ns) = scalmi(ns) + dztau(k) * sscdwm
            scal(i,j,k,ns) = scal(i,j,k,ns) + ddssc
          enddo
          scalmi(ns) = scalmi(ns) * ajchmx
!- ajoute le melange partiel(=ajcmix) depuis kdwm -> kins :
          do k=kdwm,kins
            scal(i,j,k,ns) = ajc0mx * scal(i,j,k,ns) + scalmi(ns)
          enddo
       enddo
!- fin du mixing .

!--4.4 calcule les nouveaux "b" et fait le bilan de l'ajustement convectif :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        zmix   = 0.5 * ( zw(kins+1) + zw(kdwm) )
        do k=kdwm,kins
          tloc = scal(i,j,k,1) - tkc
          sloc = scal(i,j,k,2)
          bbajc(k) = b(i,j,k)
          cfb1(i,k) = cstrho(4)*sloc
     &              + (cstrho(2)-cstrho(3)*tloc)*tloc
          cfb2(i,k) = cstrho(5)
     &              + (cstrho(6) - cstrho(7)*tloc)*tloc
     &              - (cstrho(8) + cstrho(9)*tloc)*sloc
          b(i,j,k) = gravit /
     &      (cstrho(0)+ cfb2(i,k)/(cfb1(i,k)+cfb1z(k)) ) - bref(k)
          bup(i,k) = gravit /
     &      (cstrho(0)+ cfb2(i,k)/(cfb1(i,k)+cfb1z4(k+1)) )
          bdw(i,k) = gravit /
     &      (cstrho(0)+ cfb2(i,k)/(cfb1(i,k)+cfb1z4( k )) )
          bbajc(k) = (bbajc(k)-b(i,j,k)) * (z(k)-zmix) * rho0dz(k)
          fqajc(i,j,1) = fqajc(i,j,1) + bbajc(k)
       enddo
        if ( kdwm.eq.kfs(i,j) .and. kdwm.gt.ks1 ) then
          do k=ks1,kdwm-1
           b(i,j,k) = gravit /
     &     (cstrho(0)+ cfb2(i,kdwm)/(cfb1(i,kdwm)+cfb1z(k))) - bref(k)
          enddo
        endif
!       if (kins.ge.(ks2-1)) hmajc(i,j) = max(hmajc(i,j),-z(kdwm))
!       ninstb = ninstb + kins - kdwm

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


       endif
        enddo kinsloop


       enddo longi2

!--Fin du bloc 'suppression des instabilites par permutation'.

!--4.6 Calcul de la pression reduite q :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--calcul de la pression (q) (transfere depuis 'uve')
!-   (integration de rho.g par la Methode des Trapezes) :
!
!     do i=is1(j),is2(j)
!       q(i,j,ks2)= 0.5*dz(ks2)*b(i,j,ks2)
!     enddo
!
!     do k=ks2-1,ks1,-1
!      do i=is1(j),is2(j)
!       q(i,j,k)=q(i,j,k+1)+0.5*dz4(k+1)*(b(i,j,k+1)+b(i,j,k))
!      enddo
!     enddo
!
!---------------------
!--calcul de la pression (q) , integration de rho.g "par boite pesante" :

      do i=1,imax
        qintf(i) = 0.0
      enddo
      do k=ks2,ks1,-1
        do i=is1(j),is2(j)
         dqd2 = 0.5*dz(k)*b(i,j,k)
         q(i,j,k) = qintf(i) + dqd2
         qintf(i) = q(i,j,k) + dqd2
        enddo
       enddo


!--4.7 calcul de N2 (transfere depuis 'stab') :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

       do k=ks1+1,ks2
        do i=is1(j),is2(j)
         bvf(i,j,k) = tms(i,j,k-1) * unsdzw(k) * (bup(i,k-1)-bdw(i,k))
        enddo
       enddo

!--Fin de la boucle externe sur l'indice de latitude j .
      enddo jloop2

      else
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Sans Ajustement Convectif : mixage= 0 (lstab = 0).              |
!-----------------------------------------------------------------------

!--Debut de la boucle externe sur l'indice de latitude j :
!$DIR SPP LOOP_PARALLEL
!$DIR SPP LOOP_PRIVATE(i,k,kloc,tloc,sloc,ccb1,ccb2)
!$DIR SPP LOOP_PRIVATE(bup,bdw,qintf,dqd2)
      latj: do j=js1,js2
!-----

!--5.1 Calcul de b(i,j,k) dans tout le domaine :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      do k=ks1,ks2
        do i=is1(j),is2(j)
!- conversion Degr.Celcius et tranfert -> variables locales (tloc,sloc) :
          kloc = max(k,kfs(i,j))
          tloc = scal(i,j,kloc,1) - tkc
          sloc = scal(i,j,kloc,2)
          ccb1 = cstrho(4)*sloc
     &              + (cstrho(2)-cstrho(3)*tloc)*tloc
          ccb2 = cstrho(5)
     &              + (cstrho(6) - cstrho(7)*tloc)*tloc
     &              - (cstrho(8) + cstrho(9)*tloc)*sloc
          b(i,j,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z(k)) ) - bref(k)
          bup(i,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z4(k+1)) )
          bdw(i,k) = gravit /
     &      (cstrho(0)+ ccb2/(ccb1+cfb1z4( k )) )
         enddo
       enddo

!--5.6 Calcul de la pression reduite q :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--calcul de la pression (q) (transfere depuis 'uve')
!-   (integration de rho.g par la Methode des Trapezes) :
!
!       do i=is1(j),is2(j)
!        q(i,j,ks2)= 0.5*dz(ks2)*b(i,j,ks2)
!       enddo
!
!       do k=ks2-1,ks1,-1
!        do i=is1(j),is2(j)
!         q(i,j,k)=q(i,j,k+1)+0.5*dz4(k+1)*(b(i,j,k+1)+b(i,j,k))
!        enddo
!       enddo
!
!---------------------
!--calcul de la pression (q) , integration de rho.g "par boite pesante" :

        do i=1,imax
         qintf(i) = 0.0
        enddo

        do k=ks2,ks1,-1
         do i=is1(j),is2(j)
          dqd2 = 0.5*dz(k)*b(i,j,k)
          q(i,j,k) = qintf(i) + dqd2
          qintf(i) = q(i,j,k) + dqd2
         enddo
        enddo


!--5.7 calcul de N2 (transfere depuis 'stab') :
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      do k=ks1+1,ks2
       do i=is1(j),is2(j)
        bvf(i,j,k) = tms(i,j,k-1) * unsdzw(k) * (bup(i,k-1)-bdw(i,k))
       enddo
      enddo


      enddo latj

!----------------------

!--Fin du traitement differencie selon "mixage" .
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  7 ) raccord cyclique + raccord de grille pour w, q, bvf, scal .     |
!-----------------------------------------------------------------------

      call raccord(w(1,1,1), zero, kmax, 4)
      call raccord(q(1,1,1), zero, kmax, 0)
      call raccord(b(1,1,1), zero, kmax, 0)
      call raccord(bvf(1,1,1), zero, kmax, 4)
      call raccord(scal(1,1,1,1), zero, kmax*nsmax, 0)

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine etat -
      end
