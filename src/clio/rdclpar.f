!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:51 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:51 CET 2009

      SUBROUTINE rdclpar(nn99, nresp, ncoast, ntrsm, kkwdiv,
     &              cmultr,clistr, titv,titavr,titrsm, suffix,fmt,fmtx)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Appele par "class" et "statis".
!  Preparation (titres, format ...) pour sortie sur fichier ASCII, Fmt standard
!   par lecture du fichier "class.par".
!-----
!  modif : 28/09/99

!! START_OF_USE_SECTION

      use const_mod, only:

      use para0_mod, only: imax, nsmax
      use para_mod , only: nvmax
      use vareq_mod, only: fmt1, cmult, nabs, irl1, irl2, jrl1, jrl2, spv
     &                   , cadit, cliss

      use bloc0_mod , only: ims1, ims2, js1, js2
      use bloc_mod  , only: nslpdw
      use datadc_mod, only: nvrb, nvret, nvrs, nvrt, nvrw

      use varno_mod,  only: ltyp, titcv, nvrl
      use global_constants_mod, only: dblp=>dp, ip

      use newunit_clio_mod, only: mouchard_id, clio3_out_id
!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "varno.com"
!dmr [NOEQUI]   #include "vareq.com"

!! END_OF_INCLUDE_SECTION


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

!--dummy variables :
      integer(ip), dimension(0:nvmax) :: ncoast
      integer(ip), dimension(0:*)     :: ntrsm
      real(dblp),  dimension(0:*)     :: cmultr, clistr

      character*4 suffix(-2:nvmax)
      character*30 fmt(-1:nvmax)
      character*34 titv(-2:nvmax), titavr(0:19), titrsm(0:3*nsmax+1)
      character*80 fmtx

!--variables locales :
      character*3 cct, cc3
      character*8 cc8
      character*34 cc34, notitr
      character*120 line

      integer(ip) :: nn99, nresp, kkwdiv, ll0, ll1, llm, ltp, n, na3
     &             , nflg, nfrc, nn, nn3, nnbis, no, nrespp, nv

      real(dblp)  :: speval, xx1, xx2, xx3, xx4

      integer(ip) :: class_par_id

!--instructions "data" :
! [SCRPTCOM] #include "datadc.com"

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initialisation & Definition prealables                          |
!-----------------------------------------------------------------------

      nrespp = mod(iabs(nresp),10)

!--definition des tailles des tableaux a imprimer :
      do 110 ltp=0,3
        irl1(ltp) = ims1 - 1 + mod(ltp,2)
        irl2(ltp) = ims2 + 1
        if (ltp.le.1) then
          jrl1(ltp) = js1 - 1
        else
          jrl1(ltp) = js1
        endif
        jrl2(ltp) = js2 + 1
 110  continue
!--Mise en place des tableaux tmgl et csgl ; Def. taille tableaux en sortie :
!     call defgl(irl1,irl2,jrl1,jrl2,nunion)

      if (nn99.eq.2) then
        write(mouchard_id,*) 'Taille des tableaux en sortie : '
        write(mouchard_id,*) ' irl1, irl2, jrl1, jrl2 :'
        write(mouchard_id,'(20I4)') irl1, irl2, jrl1, jrl2
!       write(mouchard_id,*) ' imdl1(jdl1->jdl2) : '
!       write(mouchard_id,'(20I4)')  (imdl1(j),j=jdl1,jdl2)
!       write(mouchard_id,*) ' imdl2, nunion : '
!       write(mouchard_id,*) imdl2, nunion
      endif

!--initialisation :
      nfrc = 0
      nslpdw = 0
      notitr = ' --- title not defined ---        '
      do 210 nv=-2,nvmax
        write(suffix(nv),'(A1,I3)') '.', 200+nv
        fmt1(nv) = 'E12.4'
        titv(nv) = notitr
 210  continue
        fmt1(-2) = 'E13.6'
      do 215 nv=0,3*nsmax+1
        cmultr(nv) = 1.
        clistr(nv) = 0.
        titrsm(nv) = notitr
        ntrsm(nv) = 0
 215  continue
      do 220 nv=0,nvmax
        nvrl(nv) = 0
        spv(nv) = 1.E9
        cmult(nv) = 1.
        cadit(nv) = 0.
        cliss(nv) = 0.
        nabs(nv) = 0
 220  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Lecture sur "class.par" des elements specifiques pour ecriture  |
!-----------------------------------------------------------------------

!--lecture sur fichier "class.par", pour chaque variable, du titre reduit,
!- du numero correspondant, de la Valeur SPeciale, des Coeff Multiplicatif
!- et Additif, du format elementaire utilise et du titre asssocie.

      open(newunit=class_par_id,file='class.par',status='old')
 1120 format(A3,1X,I3,1X,3(E11.4,1X),F8.3,1X,I3,1X,A8,1X,A3,1X,A34,1X)
      nflg = 0
 240  continue
      read(class_par_id,'(A)') line
      cct = line(:3)
      if(cct.eq.'fin') goto 250
      if(nflg.eq.1) then
        if (cct.eq.'nfr') then
          read(line,1120) cct, nn3, speval
          nfrc = nn3
          if(nn3.eq.0) nfrc = imax
          goto 240
        endif
        read(line,1120) cct, nn3, xx1, xx2, xx3, xx4,
     &                  na3, cc8, cc3, cc34
        nv = -nn3 - 10
        if (nv.ge.0 .and. nv.lt.30 .and.
     &      mod(nv,10).le.nsmax .and. nv.ne.10) then
           nv = mod(nv,10) + nsmax * (nv/10) + nv/20
           titrsm(nv) = cc34
           cmultr(nv) = xx2
           if(nresp.ne.-1) clistr(nv) = xx4
           ntrsm(nv) = na3
        elseif(nn3.ge.0.and.nn3.le.nvmax) then
          nv = nn3
          if (ltyp(nv).eq.99) goto 240
          spv(nv)   = xx1
          cmult(nv) = xx2
          cadit(nv) = xx3
          cliss(nv) = xx4
          nabs(nv) = na3
        elseif (nn3.ne.-1 .and. nn3.ne.-2 .and. nn99.eq.2) then
          write(mouchard_id,'(A,I4,3(1X,A))') 'WARNING(rdclpar), bad Var_Nb :',
     &          nn3, cct, cc3, cc34
        endif
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
        if(nn3.ge.-2.and.nn3.le.nvmax) then
          no = nn3
          if (titcv(no).ne.cct) then
            write(clio3_out_id,'(A)') 'class : Difference de Titre => Arret !'
            write(clio3_out_id,'(A,I3,4A)') ' variable No =', no,
     &          ' : ', titcv(no), ' & from "class.par" : ', cct
            stop
          endif
          titv(no) = cc34
          fmt1(no)(:8) = cc8
!- Verifie que le suffix "cc3" n'est pas deja utilise :
          nnbis = -3
          do 245 n=-2,nvmax
           if (suffix(n).eq.'.'//cc3) nnbis = n
 245      continue
          if (nnbis.eq.-3) then
            write(suffix(no),'(A1,A3)') '.', cc3
          elseif (no.ge.1 .and. no.eq.nnbis+1) then
!--2 variables 2D consecutives sur un meme fichier :
            ll0 = mod(ltyp(nnbis),12)
            ll1 = mod(ltyp(no),   12)
            if (ll0.ge.8.and.ll1.ge.8) then
              write(suffix(no),'(A1,A3)') '.', cc3
            else
             write(clio3_out_id,'(2A,2(I3,A),4A)')
     &        'WARNING : same suffix for Var',
     &       '(no=',nnbis,',',no,') ',titcv(nnbis),' & ',titcv(no),' !'
            endif
          else
             write(clio3_out_id,'(2A,2(I3,A),4A)')
     &        'WARNING : same suffix for Var',
     &       '(no=',nnbis,',',no,') ',titcv(nnbis),' & ',titcv(no),' !'
          endif
        endif
!-----
      endif
      if(cct.eq.'deb') nflg = 1
      goto 240
 250  continue
      close(class_par_id)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Modification & Definition des formats & titres annexes          |
!-----------------------------------------------------------------------

!--modification (suivant le cas traite) :

      if (nresp.eq.3) then
!--modif des titres :
        titv(nvrt)(14:17) = 'Obs.'
        titv(nvrs)(14:17) = 'Obs.'
        titv(nvrb)(10:13) = 'Obs.'
!       titv(nvrub) = 'Zonal. Wind Stess wstx  (mm2/s2)  '
!       titv(nvrvb) = 'Merid. Wind Stess wsty  (mm2/s2)  '
        titv(nvret) = ' Vert. Ekman Veloc.'//titv(nvrw)(20:34)
!       cmult(nvrub) = -1.d+06
!       cmult(nvrvb) = -1.d+06
        cmult(nvret) = cmult(nvrw)
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      do 260 nv=1,nvmax
        if (ltyp(nv).eq.99) goto 260
!- ncoast = 0 / 1 / -1 : pour tmi(llm) et le niveau du fond kniv(+/-llm).
        ltp = mod(abs(ltyp(nv)),4)
        llm = min(ltp,1)
        ncoast(nv) = sign( llm, nabs(nv) )
        if (nv.eq.nvrw) kkwdiv = nabs(nv)
        nabs(nv) = abs(nabs(nv))
        if (nresp.ge.0) nabs(nv) =  4 * ( (3+nabs(nv)) / 4 )
 260  continue
      if (nrespp.eq.5.or.nresp.eq.-1) then
        do 265 nv=0,nvmax
          if (ltyp(nv).eq.99) goto 265
          spv(nv) = 9.E9
          cadit(nv) = 0.
          cliss(nv) = 0.
          nabs(nv) = 4
          fmt1(nv) = '1PE12.4'
 265    continue
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--mise en place des formats complets :
      fmtx = '(2Hk=,I2,2(A2,1H(,I3,1H,,I3,2H):,'//fmt1(-2)
     &      //'),2(A4,'//fmt1(-2)//'))'

      do 280 nv=-1,nvmax
        write(fmt(nv),'(A,I3,A)') '(', nfrc, '('//fmt1(nv)(:8)//'))'
 280  continue

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--titre secondaire pour Moyenne Zonale :
      do 290 n=0,9
        nn = 10 + n
        titavr(n)  = titv(0)
        titavr(nn) = titv(-1)
 290  continue
      do 295 n=0,4
        nn = 10 + n
        titavr(n) = ' Grid'//titv(0)(3:14)
        titavr(nn) = 'Grid Zon.Mean '
 295  continue
        titavr(1)  = 'Abs.V G.Zon.'//titv(0)(10:14)
        titavr(5)  = 'Abs.V'//titv(0)(3:14)
        titavr(11) = 'AbsV G.Zon.Mn '
        titavr(15) = 'AbsV Zon.Mean '

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine rdclpar -
      end
