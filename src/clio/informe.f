!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009

!dmr --- 2015-10-23
!dmr --- Added the computations of freshwater for the Arctic by Frazer Davies
!dmr --- Following flag set to one activates that in the output
!dmr --- FRAZER_ARCTIC -> moved to choixcomposantes
!dmr ---

      SUBROUTINE informe(nn99)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!    prepare(titres .en tete...) et ouvre le fichier "evolu".
!    The computation are performed in inforun
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  modif : 16/11/99

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod, only: epsil

      use para0_mod, only: imax, jmax, kmax, nsmax
      use para_mod,  only: nchinf, nbsmax, nchsep, ninfmx
      use bloc0_mod, only: dx, dy, js1, js2, ks1, ks2, nitrun, spvr
     &                   , is1, iu1, dz, dzw, z, zw, scalr, is2, iu2
     &                   , daeta, fqajc, fss
      use bloc_mod,  only: ninfo, nstart, ntmoy, refexp, zsurf, zvolo
     &                   , zvols, zvolv, zvolw, ctmi, zsurfs, zsurfo, zsurfv
      use ice_mod,   only: vwx
      use dynami_mod,only:
      use reper_mod, only: dlat, fmtw, jmnor, jmsud, knor, ksud, nferme
     &                   , nocean, ylat1, zmdeta, nvinfo, nvhsf, titvar
     &                   , ktsum, scalwr, vinfom

      use newunit_clio_mod, only: clio3_out_id, mouchard_id, evolu_id
     &                   , testevolu_id

!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "reper.com"
! [SCRPTCOM] #include "ice.com"
! [SCRPTCOM] #include "dynami.com"

!! END_OF_INCLUDE_SECTION

!--variables locales :

      character*(nchinf) titinf
      character(len=8) :: fmtinf, cc8
      character(len=30) :: fmtr, fmtitr

      integer(kind=ip) :: nn99, i, j, k, kz, n, nfrinf, nn, nn0
     &                  , nninfo, nnntot, ns, nv, nv0

      real(kind=dblp), dimension(kmax) :: surfs, surfo, surfv

      real(kind=dblp)  :: ddxy, volo, vols, volv, volw, xx1, xxx
     &                  , yylatn, yylats, zlim1

      integer(kind=ip) :: run_param_id

!--variables locales a conserver d'un appel a l'autre -> dans "reper .com" .

!-----

 1000 format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
 1111 format(3(F7.1,1X,F7.3,1X),I3,A)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initialisation, ouverture et ecriture de l'entete du fichier.   |
!-----------------------------------------------------------------------

!--selon ntmoy (<-run.param) , definit "ktsum(nv)" pour "inforun" :
!  ntmoy=0 pas de moy. ; =1 moy. de qq. var. ; =2 moy de ttes les var (ss numit)
!  ktsum(nv) : =1 -> moyenne sur ninfo iter., =0 -> output ponctuel

!- initialise ktsum(nv) :
      do nv=1,ninfmx
        ktsum(nv) = 0
      enddo

!--facteur d'echelle pour "D.eta" (unite = 10e-6 m/s) :
        zmdeta = 1.D+6

!--lecture des 3 1eres lignes du fichier "run.param" :
        open(newunit=run_param_id,file='run.param',status='old')
        read(run_param_id,'(A6)') refexp
        read(run_param_id,*)
        read(run_param_id,'(A8,1X,I3)') fmtinf, nfrinf
        close(run_param_id)

!--definition des titres "titvar":
        nv = 1
        titvar(nv) = 'NoIt'    !1
        nv = nv + 1
        titvar(nv) = 'T yr'    !2
        if (ntmoy.ge.1) ktsum(nv) = 1
        nv = nv + 1
        titvar(nv) = 'EgAjC'   !3
        ktsum(nv) = 1
        nv = nv + 1
        titvar(nv) = 'V_AjC'   !4
        ktsum(nv) = 1
        nv = nv + 1
        titvar(nv) = 'D.Eta'   !5
        ktsum(nv) = 1
        nv = nv + 1
!       titvar(nv) = 'EgEta'
        titvar(nv) = 'M.Eta'   !6
        if (ntmoy.eq.2) ktsum(nv) = 1
!-------
        nv0 = nv
        nv = nv + 1
        titvar(nv) = 'DrHSF'   !7
        nv = nv + 1
        titvar(nv) = 'InHSF'   !8
        nv = nv + 1
        titvar(nv) = 'BeHSF'   !9
        nv = nv0 + nvhsf + 1
!L30    titvar(nv) = 'FLOS'
!L30    nv = nv + 1
        titvar(nv) = 'DANS'    !10
        nv = nv + 1
        titvar(nv) = 'DANN'    !11
        nv = nv + 1
        titvar(nv) = 'ISCS'    !12
        nv = nv + 1
        titvar(nv) = 'ISCN'    !13
        nv = nv + 1
        titvar(nv) = 'FRAS'    !14
        nv = nv + 1
        titvar(nv) = 'FRAN'    !15
        nv = nv + 1
        titvar(nv) = 'PNWS'
        nv = nv + 1
        titvar(nv) = 'PNWN'
        nv = nv + 1
#if ( FRAZER_ARCTIC == 1 )
        titvar(nv) = 'BASS'
        nv = nv + 1
        titvar(nv) = 'BASN'
        nv = nv + 1
#endif                         /* ! if FRAZER_ARCTIC +2 */
        titvar(nv) = 'ADGIN'
        nv = nv + 1
        titvar(nv) = 'ADPro'
        nv = nv + 1
        titvar(nv) = 'ADOut'    !20
        nv = nv + 1
        titvar(nv) = 'AABpr'
        nv = nv + 1
        titvar(nv) = 'AABex'
        nv = nv + 1
        titvar(nv) = 'AABat'
        nv = nv + 1
        titvar(nv) = 'Fc30A'
        nv = nv + 1
        titvar(nv) = 'Fs30A'
        nv = nv + 1
        titvar(nv) = 'Fsa30A'    !26
        nv = nv + 1
        titvar(nv) = 'Fsd30A'
        nv = nv + 1
        titvar(nv) = 'FsaNA'
        nv = nv + 1
        titvar(nv) = 'FsdNA'
        nv = nv + 1
        titvar(nv) = 'dtsaltA'
        nv = nv + 1
        titvar(nv) = 'Mov30A'
        nv = nv + 1
        titvar(nv) = 'Fsber'    !32
#if ( FRAZER_ARCTIC == 1 )
        nv = nv + 1
        titvar(nv) = 'BSFF'
        nv = nv + 1
        titvar(nv) = 'BSFP'
        nv = nv + 1
        titvar(nv) = 'BSFN'
        nv = nv + 1
        titvar(nv) = 'FSFF'     !30
        nv = nv + 1
        titvar(nv) = 'FSFP'
        nv = nv + 1
        titvar(nv) = 'FSFN'
        nv = nv + 1
        titvar(nv) = 'BAFF'
        nv = nv + 1
        titvar(nv) = 'BAFP'
        nv = nv + 1
        titvar(nv) = 'BAFN'     !35
        nv = nv + 1
        titvar(nv) = 'ADVE'
        nv = nv + 1
        titvar(nv) = 'DIFF'
        nv = nv + 1
        titvar(nv) = 'AOW'
        nv = nv + 1
        titvar(nv) = 'AOE'
        nv = nv + 1
        titvar(nv) = 'AOFC'     !40
        nv = nv + 1
        titvar(nv) = 'CAFF'
        nv = nv + 1
        titvar(nv) = 'CAFP'
        nv = nv + 1
        titvar(nv) = 'CAFN'
        nv = nv + 1
        titvar(nv) = 'KSFF'
        nv = nv + 1
        titvar(nv) = 'KSFP'     !45
        nv = nv + 1
        titvar(nv) = 'KSFN'
#endif                         /*  !46+2 if FRAZER_ARCTIC activated */
        if (ntmoy.ge.1) then
          do nn=1+nv0,nv
            ktsum(nn) = 1
          enddo
        endif
!-------
        nv0 = nv
        nv = nv + 1
!ic0    titvar(nv) = 'T'
        titvar(nv) = 'T-c'      !27
        scalwr(1) = -273.15
        nv = nv + 1
        titvar(nv) = 'T1-o'
        nv = nv + 1
        titvar(nv) = '|T-o|'
        nv = nv + 1
        titvar(nv) = 'S-30'     !30
        scalwr(2)  =   -30.
        nv = nv + 1
        titvar(nv) = 'S1-o'
        nv = nv + 1
        titvar(nv) = '|S-o|'   !32
        nv = nv + 1
        titvar(nv) = 'A'
        nv = nv + 1
        titvar(nv) = 'A1-o'
        nv = nv + 1
        titvar(nv) = '|A-o|'
        nv = nv + 1
        titvar(nv) = 'B'
        nv = nv + 1
        titvar(nv) = 'B1-o'
        nv = nv + 1
        titvar(nv) = '|B-o|'   !38
!---
        nv = nv0+3*nsmax + 1 ! 26+3*2 + 1= 33
#if ( OCYCC == 1 )
!dmr&pb --- Number of tracers is inconsistent between versions
!dmr&pb ---  quick fix : we do NOT provide titvars for the 11
!dmr&pb --- additionnal tracers of the oceanic carbon cycle
        nv = nv0 + 3*2 + 1
#endif
#if ( ISM != 0 && ISOOCN != 0)
!mab: test to fix problem with evolu file, nsmax can be everything
!mab: when applying GRISLI and isotopes
        nv = nv0 + 3*4 + 1
#endif

        titvar(nv) = '|w|'
        if (ntmoy.eq.2) ktsum(nv) = 1
        nv = nv + 1
        titvar(nv) = '|u|'
        if (ntmoy.eq.2) ktsum(nv) = 1
        nv = nv + 1
        titvar(nv) = '|v|'
        if (ntmoy.eq.2) ktsum(nv) = 1
        nv = nv + 1
        titvar(nv) = 'K.E'  !36
        if (ntmoy.eq.2) ktsum(nv) = 1
!       do k=ks1,ks2
!         nv = nv + 1
!         write(titvar(nv),'(A3,I2)') 'ke.',k
!       enddo
        if (ninfo.ge.0) then
!- Moyenne des valeurs absolues des differences :
          do k=ks1,ks2
            nv = nv + 1
            write(titvar(nv),'(A3,I2)') 'T-o',k
          enddo
          do k=ks1,ks2
            nv = nv + 1
            write(titvar(nv),'(A3,I2)') 'S-o',k
          enddo
        else
            nv = nv + 1
            write(titvar(nv),'(A3,I2)') 'T-o',1
!- Moyenne des differences :
          do k=1+ks1,ks2
            nv = nv + 1
            write(titvar(nv),'(A1,I2,A2)') 'T',k,'-o'
          enddo
            nv = nv + 1
            write(titvar(nv),'(A3,I2)') 'S-o',1
          do k=1+ks1,ks2
            nv = nv + 1
            write(titvar(nv),'(A1,I2,A2)') 'S',k,'-o'
          enddo
        endif
        if (ntmoy.eq.2) then
          do nn=1+nv0,nv
            ktsum(nn) = 1
          enddo
        endif
!-------

        nocean = nv  ! 76 ?
        nv = nv + 1
        titvar(nv) = 'AIEFN'
        nv = nv + 1
        titvar(nv) = 'AIEFS'
        nv = nv + 1
        titvar(nv) = 'A15N'
        nv = nv + 1
        titvar(nv) = 'A15S'
        nv = nv + 1
        titvar(nv) = 'A85N'
        nv = nv + 1
        titvar(nv) = 'A85S'
        nv = nv + 1
        titvar(nv) = 'ALEN'
        nv = nv + 1
        titvar(nv) = 'ALES'
        nv = nv + 1
        titvar(nv) = 'VOLN'
        nv = nv + 1
        titvar(nv) = 'VOLS'
        nv = nv + 1
        titvar(nv) = 'VONN'
        nv = nv + 1
        titvar(nv) = 'VONS'
        nv = nv + 1
        titvar(nv) = 'ECGN'
        nv = nv + 1
        titvar(nv) = 'ECGS'
        nv = nv + 1
        titvar(nv) = 'FRAG'
        nv = nv + 1
#if ( FRAZER_ARCTIC == 1 )
        titvar(nv) = 'FSIF'
        nv = nv + 1
#endif
        titvar(nv) = 'SPNG'
        nv = nv + 1
#if ( FRAZER_ARCTIC == 1 )
        titvar(nv) = 'KSIF'
        nv = nv + 1
#endif
        titvar(nv) = 'BERG'
        nv = nv + 1
#if ( FRAZER_ARCTIC == 1 )
        titvar(nv) = 'BSIF'
        nv = nv + 1
        titvar(nv) = 'CAIF'
        nv = nv + 1
        titvar(nv) = 'BAIF'
        nv = nv + 1
#endif
        titvar(nv) = 'ThEx'
        nv = nv + 1
        titvar(nv) = 'ISMM'
        nv = nv + 1
        titvar(nv) = 'IcbN'
        nv = nv + 1
        titvar(nv) = 'IcbS'
        if (ntmoy.ge.1) then
          do nn=1+nocean,nv
            ktsum(nn) = 1
          enddo
        endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        nvinfo = nv
        if(nvinfo.gt.ninfmx) then
          write(clio3_out_id,*) 'Arret ! Depassement Nombre Max de variables'
     &      //' traitees par la routine "informe"'
          stop
        endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Definition et Ecriture de l'entete :
!- nombre d'enregistrements :
        nninfo = abs(ninfo)
        nn0 = (nstart - 1) / nninfo
        if(nstart.eq.1) nn0 = -1
        nferme = (nstart - 1 + nitrun) / nninfo
        nnntot = nferme - nn0
        nn0 = nninfo * (1 + nn0)
        nferme = nninfo * nferme

!- definition des formats :
        write(fmtw,'(A,I3,A,I1,A)')
     &    '(',nfrinf,'(A',nchsep,','//fmtinf//'))'
        write(fmtr,'(A,I3,A,I1,A)')
     &    '(',nfrinf,'(',nchsep,'X,'//fmtinf//'))'
        write(fmtitr,'(A,I3,A,I1,A)') '(',ninfmx,'A',nchinf,')'

!--Ouverture du fichier "evolu" :
!ray    irecl = nvinfo*nchinf + nchinf
!ray    open(90,file='evolu',status='unknown',RECL=irecl)
        open(newunit=evolu_id,file='outputdata/ocean/evolu',status='unknown')
#if ( ISM >= 2 )
        open(newunit=testevolu_id,file='outputdata/ocean/testevolu'
     &  ,form='formatted',status='unknown')
#endif
!- ecriture de 2 lignes d entete :
        write(evolu_id,1000) fmtr, spvr, nvinfo, nnntot, 0, nfrinf
#if ( ISM >= 2 )
        write(testevolu_id,1000) fmtr, spvr, nvinfo, nnntot, 0, nfrinf
#endif
        xxx = 0.001 * DFLOAT(nninfo)
        xx1 = 0.001 * DFLOAT(nn0)
        write(cc8,'(A,I1)') ' ntmoy=', ntmoy
        write(evolu_id,1111) DFLOAT(nchinf), 0., xx1, xxx, 0., 0., 0, cc8
#if ( ISM >= 2 )
        write(testevolu_id,1111) DFLOAT(nchinf), 0., xx1, xxx, 0., 0., 0, cc8
#endif

!- ecriture de 2 lignes de titre :
        write(evolu_id,'(A,I8,A,I8,A,I5)')
     &     'Evolution chronologique - Experience '//refexp
     &   //'   de', nn0, ' a', nferme, ' pas', nninfo
        write(evolu_id,fmtitr) (titvar(nv),nv=1,nvinfo)
#if ( ISM >= 2 )
        write(testevolu_id,'(A,I8,A,I8,A,I5)')
     &     'Evolution chronologique - Experience '//refexp
     &   //'   de', nn0, ' a', nferme, ' pas', nninfo
        write(testevolu_id,fmtitr) (titvar(nv),nv=1,nvinfo)
#endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--preparation de "titvar" pour l ecriture parmi les valeurs numeriques :
        do nv=2,nvinfo
          titinf = titvar(nv)(:nchinf)
          titvar(nv) = '  '//titinf
        enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

! initialisation of vwx
        do j=1,57
         do k=1,kmax+1
          do n=0,nbsmax
           vwx(j,k,n)=0.0
          enddo
         enddo
        enddo
        write(mouchard_id,*) 'vwx',vwx(43,1,3)
!--calcul des differents "volumes" intervenant dans les moyennes :
        ddxy = dx * dy
        vols = 0.
        volo = 0.
        volv = 0.
        volw = 0.
        do k=ks1,ks2
          surfs(k) = 0.
          surfo(k) = 0.
          surfv(k) = 0.
          do j=js1,js2
            do i=is1(j),is2(j)
              surfs(k) = surfs(k) + ctmi(i,j,k,0)
              surfo(k) = surfo(k) +
     &          min( ctmi(i,j,k,0), (scalr(i,j,k,1)-spvr) )
            enddo
            do i=iu1(j),iu2(j)
              surfv(k) = surfv(k) + ctmi(i,j,k,1)
            enddo
          enddo
          vols = vols + surfs(k) * dz(k)
          volo = volo + surfo(k) * dz(k)
          volv = volv + surfv(k) * dz(k)
          if(k.lt.ks2) volw = volw + surfs(k) * dzw(k+1)
        enddo
        zvols = 0.
        zvolo = 0.
        zvolv = 0.
        zvolw = 0.
        if(vols.gt.epsil) zvols = 1. / vols
        if(volo.gt.epsil) zvolo = 1. / volo
        if(volv.gt.epsil) zvolv = 1. / volv
        if(volw.gt.epsil) zvolw = 1. / volw
        vols = vols * ddxy
        volo = volo * ddxy
        volv = volv * ddxy
        volw = volw * ddxy
        do k=ks1,ks2
          if(surfs(k).gt.epsil) then
            zsurfs(k) = 1. / surfs(k)
          else
            zsurfs(k) = 0.
          endif
          if(surfo(k).gt.epsil) then
            zsurfo(k) = 1. / surfo(k)
          else
            zsurfo(k) = 0.
          endif
          if(surfv(k).gt.epsil) then
            zsurfv(k) = 1. / surfv(k)
          else
            zsurfv(k) = 0.
          endif
          surfs(k)  = surfs(k) * ddxy
          surfo(k)  = surfo(k) * ddxy
          surfv(k)  = surfv(k) * ddxy
        enddo
        zsurf = zsurfs(ks2)
        zmdeta = zmdeta * zsurf

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--initialisation des variables utilisees uniquement pour "informe":
        do j=1,jmax
         do i=1,imax
          daeta(i,j) = 0.
         enddo
        enddo

!--Remise a zero de fss, Energ. Aj.Conv & Freq.Aj.Conv (avant 1ere it) :
        do k=1,kmax
         do j=1,jmax
          do i=1,imax
            fqajc(i,j,k) = 0.
          enddo
         enddo
        enddo
        do ns=0,nsmax
         do j=1,jmax
          do i=1,imax
           fss(i,j,ns) = 0.
          enddo
         enddo
        enddo

!--Initialisation of the arrays for the accumulation
        do nv=1,nvinfo
          vinfom(nv)=0.
        enddo

!--Definition of the limits for the downsloping budget
      ksud = ks2
      zlim1 = -1000.0
      do kz=ks2,ks1,-1
        if (z(kz).ge.zlim1) ksud = kz
      enddo
      knor = ks2
      zlim1 = -450.0
      do kz=ks2,ks1,-1
        if (z(kz).ge.zlim1) knor = kz
      enddo
      yylats = -55.
      jmsud = nint( (yylats-ylat1)/dlat ) + 1
      yylatn = 70.
      jmnor = nint( (yylatn-ylat1)/dlat ) + 1

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        if (nn99.eq.2) then
!--------
          write(mouchard_id,*) 'nvinfo, nocean, ntmoy :', nvinfo, nocean, ntmoy
          write(mouchard_id,*) 'zvols , vols :', zvols, vols
          write(mouchard_id,*) 'zvolo , volo :', zvolo, volo
          write(mouchard_id,*) 'zvolv , volv :', zvolv, volv
          write(mouchard_id,*) 'zvolw , volw :', zvolw, volw
          write(mouchard_id,*) 'zsurfs, surfs (ks2) :', zsurfs(ks2), surfs(ks2)
          write(mouchard_id,'(1P6E12.5)') (zsurfs(k),k=ks1,ks2)
          write(mouchard_id,'(1P6E12.5)') ( surfs(k),k=ks1,ks2)
          write(mouchard_id,*) 'zsurfo, surfo (ks2) :', zsurfo(ks2), surfo(ks2)
          write(mouchard_id,'(1P6E12.5)') (zsurfo(k),k=ks1,ks2)
          write(mouchard_id,'(1P6E12.5)') ( surfo(k),k=ks1,ks2)
          write(mouchard_id,*) 'zsurfv, surfv (ks2) :', zsurfv(ks2), surfv(ks2)
          write(mouchard_id,'(1P6E12.5)') (zsurfv(k),k=ks1,ks2)
          write(mouchard_id,'(1P6E12.5)') ( surfv(k),k=ks1,ks2)
          write(mouchard_id,'(A,2I4,2F9.2)')' downsloping jmnor,knor,lat_N,zw :'
     &                   , jmnor, knor, ylat1+(jmnor-1)*dlat, zw(knor)
          write(mouchard_id,'(A,2I4,2F9.2)')' downsloping jmsud,ksud,lat_S,zw :'
     &                   , jmsud, ksud, ylat1+(jmsud-1)*dlat, zw(ksud)
          write(mouchard_id,*)
!--------
          write(mouchard_id,'(A,I4,2A,I4,A)') 'File "evolu" :',nnntot,' time ',
     &    'records of',nvinfo,' (=nvinfo) var. ; nv,title,ktsum :'
          write(mouchard_id,'(4(A,I3,1X,A,I2,A1))')
     &      (' no=',n, titvar(n), ktsum(n), ',', n=1,nvinfo)
          write(mouchard_id,*)
!--------
        endif
      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine informe -
      end subroutine informe
