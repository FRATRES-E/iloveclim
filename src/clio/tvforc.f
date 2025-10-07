!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:56 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:56 CET 2009

      SUBROUTINE tvforc(nn99)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  called by "clio" : set-up time-dependant, zonally uniform forcing
!  1ere itt : Lit le fichier "chronos" et initialise ;
!- Convention :
!   kvar  : 1 Wst ; 2 Phis(ns) ; 3 Surf_Obs ; 4 Rappes(ns)
!   nsvar : No_Scal concerne (si Wst, 1=wstx,2=wsty).
!   ktyp(<- ktyp/3) :  = 0 remplace , = 1 additione , = 2 multiplie
!   kspv(=mod(ktyp,3)) =0 spv sans effet, =1 modif si spv, =2 modif sauf si spv
!  modification locale ([is,ie]x[js,je] -> 2nd ligne) definie par valeur unique
!   ou globale par valeur correspondante lue sur fichier(-> 2nd ligne).
!--------
!  modif : 14/05/98


!! START_OF_USE_SECTION

      use const_mod, only: epsil, one, radian, yeaday, zero

      use para0_mod, only: imax, jmax, nsmax, ltest
      use para_mod,  only:
      use bloc0_mod, only: scalr, ks2, nsav, tms, tpstot
      use bloc_mod,  only: ninfo, nstart, numit
      use reper_mod, only: dxwi, dyai, dywj, icheck, xwi1, xwpoln, ywj1
     &                   , jsep
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

      integer(ip), parameter :: ntmax = 1000 , jzmax = 9

!- dummy variables :

!--variables a conserver d'un appel a l'autre :

      integer(ip), dimension(imax,jmax), save :: jztab
      integer(ip), save :: ntabc, kvar, nsvar, ktyp, kspv, jzm, ntm

      real(dblp),                         save :: ytm1, dytm
      real(dblp), dimension(jzmax),       save :: ysepz, t0zvar
      real(dblp), dimension(0:jzmax),     save :: t1zvar
      real(dblp), dimension(jzmax,ntmax), save :: tzvar

!- variables locales :
      integer(ip), dimension(2)        :: kloc
      real(dblp), dimension(imax,jmax) :: varloc
      real(dblp), dimension(jmax)      :: ccxam
      character*30 fmtvfc
      character*70 line

      integer(ip) :: nn99, i, j, jj, jjz, n, nnsmn, nnsmx, nnvmn, nnvmx
     &             , ns, nt, ntj
      real(dblp)  :: ccya, ccyvr, degre, spv, xw, xxi1, ya, yrtps, ytj
     &             , ytt, yw, yyalim, yyj1, yyv

      integer(ip) :: chronos_id

      if (numit.eq.nstart) then
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Ininitialisation Lecture du fichier "chronos" .
!-----------------------------------------------------------------------

      write(clio3_out_id,*)
     &  'tvforc : Time dependant Forcing ; read file=chronos'
      if (nn99.eq.2) write(99,'(2A)')
     &  'tvforc : Time dependant Forcing ; read file=chronos :'

      open(newunit=chronos_id, file='chronos', status = 'old')
      read(chronos_id,*)
      read(chronos_id,*)
      read(chronos_id,*)
      read(chronos_id,*)

      read(chronos_id,*) ntabc
      ntabc = 1
      do nt=1,ntabc
        read(29,'(A)') line
        if (nn99.eq.2) write(mouchard_id,'(A)') line

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Definition du type de modification :
      read(chronos_id,*) kvar, nsvar, ktyp, spv
      read(chronos_id,*) jzm, ntm, ytm1, dytm
      kspv = mod(ktyp,3)
      ktyp = ktyp / 3

!- Verification & Exclusion :
      nnvmn = 1
      nnvmx = 4
      nnsmn = 1
      nnsmx = nsmax
      if (kvar.eq.1) nnsmx = 2
      if (jzm.gt.jzmax .or. ntm.gt.ntmax) then
        write(clio3_out_id,*) 'tvforc : ARRET , depassement d''indice !'
        write(clio3_out_id,*) 'Max : jzmax,ntmax =', jzmax,ntmax
        write(clio3_out_id,*) 'real:  jzm , ntm  =', jzm , ntm
        close(chronos_id)
        stop
      elseif (kvar.lt.nnvmn .or. kvar.gt.nnvmx .or.
     &        nsvar.lt.nnsmn .or. nsvar.gt.nnsmx) then
        write(clio3_out_id,*)
     &   'tvforc : ARRET , Mauvais choix de variable :'
        write(clio3_out_id,*) 'min : nnvmn,nnsmn =', nnvmn, nnsmn
        write(clio3_out_id,*) 'Max : nnvmx,nnsmx =', nnvmx, nnsmx
        write(clio3_out_id,*) 'real: kvar, nsvar =', kvar, nsvar
        close(chronos_id)
        stop
      elseif (kvar.ne.3) then
!- Cas pas encore traites :
        write(clio3_out_id,'(A,2I4,A)') 'tvforc : ARRET ,
     &                      Cas kvar,kspv=', kvar,kspv, ' Inactif !'
        close(chronos_id)
        stop
      endif
      read(chronos_id,*) (ysepz(j),j=1,jzm-1)
      read(chronos_id,'(A)') fmtvfc

!- Lecture de la serie chronologique :
      do n=1,ntm
        read(chronos_id,fmtvfc) (tzvar(j,n),j=1,jzm)
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Mise en place (1ere itt) du tab. (jztab) definissant les zones :
!-----------------------------------------------------------------------

!- initialisation :
      if (ktyp.eq.2) then
        t1zvar(0) = 1.
      elseif (ktyp.eq.1) then
        t1zvar(0) = 0.
      else
        t1zvar(0) = spv
      endif
      do j=1,jzm
        t0zvar(j) = 0.
        if (ktyp.eq.3) t0zvar(j) = 1.
      enddo

!--Mise en place de de "jztab" - Differenciation Centre / Coin des mailles :
      xxi1 = xwi1
      yyj1 = ywj1
      if (kvar.eq.1) then
!--Coins des mailles :
        xxi1 = xwi1 - 0.5 * dxwi
        yyj1 = ywj1 - 0.5 * dywj
      endif
!-----
      degre = 1.0 / radian
      do j=1,jmax
        yw = yyj1 + dywj * DFLOAT(j-1)
        ccxam(j)= sin(yw * radian)
        jjz = 1
        do jj=1,jzm-1
          if (yw.ge.ysepz(jj)) jjz = jjz + 1
        enddo
        do i=1,imax
!         yvrlat(i,j) = yw
          jztab(i,j) = jjz
        enddo
      enddo
      if (ltest.gt.2) then
!-----
      yyalim = 90.0 - abs(dyai)
      do i=1,imax
       xw = xxi1 + dxwi * DFLOAT(i-1)
       ya = 90.0 + xwpoln - xw
       if (abs(ya).ge.yyalim) cycle
       ccya = cos(ya * radian)
       do j=jsep(i),jmax
         ccyvr = ccya * ccxam(j)
         yyv = degre * asin(ccyvr)
!        yvrlat(i,j) = yyv
         jztab(i,j) = 1
         do jj=1,jzm-1
          if (yyv.ge.ysepz(jj)) jztab(i,j) = jztab(i,j) + 1
         enddo
       enddo
      enddo
!-----
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--Traitement des valeurs speciales :
      kloc(1) = 1
      kloc(2) = 0
      if (kspv.eq.2) then
        kloc(1) = 0
        kloc(2) = 1
      endif
      if (kvar.eq.3 .and. kspv.ne.0) then
        do j=1,jmax
         do i=1,imax
          if (scalr(i,j,ks2,nsvar).eq.spv.or.tms(i,j,ks2).eq.0) then
            jztab(i,j) = jztab(i,j) * kloc(1)
          else
            jztab(i,j) = jztab(i,j) * kloc(2)
          endif
         enddo
        enddo
      elseif (kspv.ne.0) then
        write(clio3_out_id,*) 'pas programme !'
      endif
!-----

      enddo
      close(chronos_id)

      if (nn99.eq.2) then
        write(mouchard_id,*) 'tvforc, tableau jztab :'
        do j=jmax,1,-1
          write(mouchard_id,'(122Z1,I3)') (jztab(i,j),i=1,imax),j
        enddo
        write(mouchard_id,*)
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin du traitement specifique a la 1ere itt.
      endif

!     do 500 nt=1,ntabc
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Definition du nouveau tableau :
!-----------------------------------------------------------------------

!- positionnement (ntj,ytj) dans la serie temp. (unite = yr)
      yrtps = tpstot / (86400. * yeaday)
      ytj = 1. + (yrtps - ytm1) / dytm
      ntj = max( 1, min(ntm-1,nint(ytj-0.5)) )
      ytj = ytj - ntj
      ytj = max( zero, min(one,ytj) )

!- mise en place des valeurs sur les differentes zones :
      do j=1,jzm
        ytt = (1.-ytj)*tzvar(j,ntj) + ytj*tzvar(j,ntj+1)
        t1zvar(j) = ytt
        if (ktyp.eq.2 .and. abs(t0zvar(j)).gt.epsil) then
          t1zvar(j) = ytt / t0zvar(j)
        elseif (ktyp.eq.1) then
          t1zvar(j) = ytt - t0zvar(j)
        endif
        t0zvar(j) = ytt
      enddo

      if (nn99.eq.2 .and. mod(numit,ninfo).eq.0) then
        write(mouchard_id,'(2A,I9,F12.3,I8,F10.6)') 'tvforc:',
     &   ' it,tps(yr), ntj,ytj / val(lat)=', numit, yrtps, ntj,ytj
        write(mouchard_id,'(1P6E12.4)') (t0zvar(j),j=1,jzm)
!       write(mouchard_id,'(1P6E12.4)') (tzvar(j,ntj),j=1,jzm)
!       write(mouchard_id,'(1P6E12.4)') (tzvar(j,ntj+1),j=1,jzm)
!       write(mouchard_id,'(1P6E12.4)') (t1zvar(j),j=1,jzm)
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Traitements des differents cas de tableau modifie :
!-----------------------------------------------------------------------

!     if (kvar.lt.2) then
!       write(clio3_out_id,*) 'pas programme !'
!     elseif (kvar.lt.3) then
!--Traitement de scalr(*,*,ks2,ns) :
        ns = nsvar
        if (ktyp.eq.3) then
!- Multiplie :
          do j=1,jmax
           do i=1,imax
            scalr(i,j,ks2,ns) = scalr(i,j,ks2,ns) * t1zvar(jztab(i,j))
           enddo
          enddo
        elseif (ktyp.eq.1) then
!- Ajoute :
          do j=1,jmax
           do i=1,imax
            scalr(i,j,ks2,ns) = scalr(i,j,ks2,ns) + t1zvar(jztab(i,j))
           enddo
          enddo
        else
!- Remplace :
          do j=1,jmax
           do i=1,imax
            scalr(i,j,ks2,ns) = t1zvar(jztab(i,j))
           enddo
          enddo
        endif
        if (nn99.eq.2 .and. icheck.ge.1 .and. mod(numit,nsav).eq.0)
     &   write(mouchard_id,'(A,I4,A,/,(1P10E11.3))')
     &   ' scalr(i=', icheck, ') :', (scalr(icheck,j,ks2,ns),j=1,jmax)
!     else
!       write(clio3_out_id,*) 'pas programme !'
!     endif

!--fin du traitement d'un tableau 2 D.
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--fin du traitement de la modification.
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! 500  continue

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine tvforc -
      end
