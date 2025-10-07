!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:47 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:47 CET 2009

      SUBROUTINE foroutp(irn,jrn)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! accompagne l'insertion du bloc "varno .com"
! initialisation de krl1(nv), krlm(nv), ltyp(nv)
!-----
! Ctke [Ctk0] => ligne specifique a la version avec [sans] TKE .
! Cice [Cic0] => ligne specifique a la version avec [sans] glace marine .
! Cncd [Cnc0] => ligne specifique a la version avec [sans] NetCDF .
!-----
!  modif : 27/09/99

!! START_OF_USE_SECTION

      use const_mod,  only:
      use para0_mod,  only: imax, jmax, kmax, nsmax, ltest
      use para_mod,   only: krlajc, krlal, krlas, krlau, krlavi, krlavs
     &  , krlaxt, krlayt, krlb, krlet, krlfc, krlfq, krlfs, krlfs3
     &  , krlfs4, krlfw, krlhac, krlhat, krlhdt, krlhg, krlhn, krlmin
     &  , krln2, krlps, krlpsx, krlpsy, krlqs, krls, krls3, krls4
     &  , krlslx, krlsly, krlt, krltke, krlts, krlu, krlub, krlug
     &  , krlugm, krlusl, krlv, krlvaf, krlvat, krlvb, krlvdf, krlvdt
     &  , krlvg, krlvgm, krlvsl, krlw, krlwgm, kv2dmx, nkb0, nv2dmx
     &  , nvmax
      use bloc0_mod,  only:
      use bloc_mod,   only:
      use datadc_mod, only: nvrajc, nvral, nvras, nvrau, nvravi, nvravs
     &  , nvraxt, nvrayt, nvrb, nvreac, nvret, nvrfc, nvrfq, nvrfs
     &  , nvrfs3, nvrfs4, nvrfw, nvrhac, nvrhat, nvrhdt, nvrhg, nvrhn
     &  , nvrmom, nvrn2, nvrps, nvrpsx, nvrpsy, nvrqs, nvrs, nvrs3
     &  , nvrs4, nvrslx, nvrsly, nvrt, nvrtbq, nvrtgx, nvrtgy, nvrtke
     &  , nvrts, nvru, nvrub, nvrug, nvrugm, nvrum, nvrusl, nvrv, nvrvaf
     &  , nvrvat, nvrvb, nvrvdf, nvrvdt, nvrvg, nvrvgm, nvrvm, nvrvsl
     &  , nvrw, nvrwgm, nvrxzo

      use varno_mod,  only: ltyp, krlm, krl1, titcv, nvrl, kvar2d

      use global_constants_mod, only: dblp=>dp, ip

      use newunit_clio_mod, only: clio3_out_id      

!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

c~ #include "type.com"
! [SCRPTCOM] #include "para.com"
! [SCRPTCOM] #include "bloc.com"
! [SCRPTCOM] #include "varno.com"

!! END_OF_INCLUDE_SECTION

!--dummy variables :
      integer(ip), dimension(imax,8) :: irn
      integer(ip), dimension(jmax,8) :: jrn

!! --- local variables
      integer(ip) :: i,j,kk2d,lln,nn,nn2d,nv

!--instructions "data" :
! [SCRPTCOM] #include "datadc.com"

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initilisation
!-----------------------------------------------------------------------

      do 110 nv=-2,nvmax
        titcv(nv) = '---'
 110  continue
!- default_value(krl1) must remains < krlmin
      do 130 nv=0,nvmax
        ltyp(nv) = 99
        krlm(nv) = 0
        krl1(nv) = -99
        nvrl(nv) = 0
        kvar2d(nv) = 0
 130  continue

!--Definition du numero de chaque variable "nvr..."
!-   -> instruction data dans "varno .com"

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Definition du titre (reduit a 3c) de chaque variable :          |
!-----------------------------------------------------------------------

      titcv(-2) = 'Bil'
      titcv(-1) = 'TrM'
      titcv( 0) = ' Av'
      titcv(nvrt  ) = ' Av'
      titcv(nvrt  ) = ' T '
      titcv(nvrs  ) = ' S '
      titcv(nvrs3 ) = ' sA'
      titcv(nvrs4 ) = ' sB'
      titcv(nvrw  ) = ' W '
      titcv(nvru  ) = ' U '
      titcv(nvrv  ) = ' V '
      titcv(nvrtke) = 'q2t'
!tk0  titcv(nvrtke) = ' Ke'
      titcv(nvrub ) = ' Ub'
      titcv(nvrvb ) = ' Vb'
      titcv(nvret ) = 'eta'
      titcv(nvrps ) = 'PsH'
      titcv(nvrfw ) = 'FxW'
      titcv(nvrfc ) = 'FxT'
      titcv(nvrfs ) = 'FxS'
      titcv(nvrfs3) = 'FxA'
      titcv(nvrfs4) = 'FxB'
      titcv(nvreac) = 'Eac'
      titcv(nvrajc) = 'AjC'
      titcv(nvrau ) = 'Avu'
      titcv(nvras ) = 'Avs'
      titcv(nvrb  ) = ' b '
      titcv(nvrn2 ) = ' N2'

      titcv(nvrusl) = 'Usl'
      titcv(nvrvsl) = 'Vsl'
      titcv(nvrhac) = 'Hac'

      titcv(nvrhg ) = ' hg'
      titcv(nvrfq ) = ' fq'
      titcv(nvrqs ) = ' qs'
      titcv(nvral ) = ' al'
      titcv(nvrhn ) = ' hn'
      titcv(nvrts ) = ' ts'
      titcv(nvrum ) = ' um'
      titcv(nvrvm ) = ' vm'
      titcv(nvrug ) = ' ug'
      titcv(nvrvg ) = ' vg'
      titcv(nvrtbq) = 'tbq'
      titcv(nvrxzo) = ' z0'
      titcv(nvrtgx) = 'tgx'
      titcv(nvrtgy) = 'tgy'
      titcv(nvrmom) = 'MOM'

      titcv(nvraxt) = 'axt'
      titcv(nvrayt) = 'ayt'
      titcv(nvrhat) = 'hat'
      titcv(nvrhdt) = 'hdt'
      titcv(nvrvat) = 'vat'
      titcv(nvrvdt) = 'vdt'
      titcv(nvrvaf) = 'vaf'
      titcv(nvrvdf) = 'vdf'

      titcv(nvravi) = 'Avi'
      titcv(nvravs) = 'AvT'
      titcv(nvrslx) = 'slx'
      titcv(nvrsly) = 'sly'
      titcv(nvrpsx) = 'Psx'
      titcv(nvrpsy) = 'Psy'
      titcv(nvrugm) = 'Ugm'
      titcv(nvrvgm) = 'Vgm'
      titcv(nvrwgm) = 'Wgm'

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Definition des "localisation" de la vaviable sur la maille :    |
!-----------------------------------------------------------------------

      ltyp(0) = -1
      do 310 nn=0,nsmax-1
        ltyp(nn+nvrt) =  0
 310  continue
      ltyp(nvrw  ) =  4
      ltyp(nvru  ) =  3 + 12
      ltyp(nvrv  ) =  3 + 24
      ltyp(nvrtke) =  4
!tk0  ltyp(nvrtke) =  3
      ltyp(nvrub ) = 11 + 12
      ltyp(nvrvb ) = 11 + 24
      ltyp(nvret ) =  8
      ltyp(nvrps ) = 11
      do 330 nn=0,nsmax
        ltyp(nn+nvrfw) =  8
 330  continue
      ltyp(nvreac) =  8
      ltyp(nvrajc) =  4
      ltyp(nvrau ) =  7
      ltyp(nvras ) =  4
      ltyp(nvrb  ) =  0
      ltyp(nvrn2 ) =  4

      ltyp(nvrusl) =  9
      ltyp(nvrvsl) = 10
      ltyp(nvrhac) =  8

      ltyp(nvravi) =  4
      ltyp(nvravs) =  4
      ltyp(nvrslx) =  4
      ltyp(nvrsly) =  4
      ltyp(nvrpsx) =  0
      ltyp(nvrpsy) =  0
      ltyp(nvrugm) =  0 + 12
      ltyp(nvrvgm) =  0 + 24
      ltyp(nvrwgm) =  0

      ltyp(nvrhg ) =  8
      ltyp(nvrfq ) =  8
      ltyp(nvrqs ) =  8
      ltyp(nvral ) =  8
      ltyp(nvrhn ) =  8
      ltyp(nvrts ) =  8
      ltyp(nvrum ) = 11 + 12
      ltyp(nvrvm ) = 11 + 24
      ltyp(nvrug ) = 11 + 12
      ltyp(nvrvg ) = 11 + 24
      ltyp(nvrtbq) =  8
      ltyp(nvrxzo) =  8
      ltyp(nvrtgx) = 11 + 12
      ltyp(nvrtgy) = 11 + 24
      ltyp(nvrmom) =  8

      ltyp(nvraxt) =  1
      ltyp(nvrayt) =  2
      ltyp(nvrhat) =  0
      ltyp(nvrhdt) =  0
      ltyp(nvrvat) =  0
      ltyp(nvrvdt) =  0
      ltyp(nvrvaf) =  4
      ltyp(nvrvdf) =  4

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) correspondance varible"nv" <-> niveaux"krl" associes :
!-----------------------------------------------------------------------

!- Definit la correspondance par les tableaux krl1(nv), krlm(nv)
!  N.B.: krl1(nv) >= krlmin uniquement si place reservee dans "vrl()" pour "nv"
      krl1(nvrt  ) = krlt
      krl1(nvrs  ) = krls
      krl1(nvrs3 ) = krls3
      krl1(nvrs4 ) = krls4
      krl1(nvrw  ) = krlw
      krl1(nvru  ) = krlu
      krl1(nvrv  ) = krlv
      krl1(nvrtke) = krltke
      krl1(nvrub ) = krlub
      krl1(nvrvb ) = krlvb
      krl1(nvret ) = krlet
      krl1(nvrps ) = krlps
      krl1(nvrfw ) = krlfw
      krl1(nvrfc ) = krlfc
      krl1(nvrfs ) = krlfs
      krl1(nvrfs3) = krlfs3
      krl1(nvrfs4) = krlfs4
      krl1(nvreac) = krlajc
      krl1(nvrajc) = krlajc
      krl1(nvrau ) = krlau
      krl1(nvras ) = krlas
      krl1(nvrb  ) = krlb
      krl1(nvrn2 ) = krln2

      krl1(nvrusl) = krlusl
      krl1(nvrvsl) = krlvsl
      krl1(nvrhac) = krlhac

      krl1(nvrhg ) = krlhg
      krl1(nvrfq ) = krlfq
      krl1(nvrqs ) = krlqs
      krl1(nvral ) = krlal
      krl1(nvrhn ) = krlhn
      krl1(nvrts ) = krlts

      krl1(nvrug ) = krlug
      krl1(nvrvg ) = krlvg

      krl1(nvraxt) = krlaxt
      krl1(nvrayt) = krlayt
      krl1(nvrhat) = krlhat
      krl1(nvrhdt) = krlhdt
      krl1(nvrvat) = krlvat
      krl1(nvrvdt) = krlvdt
      krl1(nvrvaf) = krlvaf
      krl1(nvrvdf) = krlvdf

      krl1(nvravi) = krlavi
      krl1(nvravs) = krlavs
      krl1(nvrslx) = krlslx
      krl1(nvrsly) = krlsly
      krl1(nvrpsx) = krlpsx
      krl1(nvrpsy) = krlpsy
      krl1(nvrugm) = krlugm
      krl1(nvrvgm) = krlvgm
      krl1(nvrwgm) = krlwgm


!--Mise en place de "krlm" = nv de niveaux de la variable "nv" :
      do 450 nv=1,nvmax
        if (ltyp(nv).eq.99) goto 450
        krlm(nv) = kmax
        lln = mod(ltyp(nv),12) / 4
        if(lln.eq.2) krlm(nv) = 1
 450  continue
      krlm(nvrusl) = 2 + nsmax
      krlm(nvrvsl) = 2 + nsmax
      krlm(nvrtbq) = nkb0
      krlm(nvrmom) = 35

!--Mise en place de "kvar2d" = indice reperant les variables non 3D(i,j,k) :
!     write( 6,'(A)') ' titcv  nv   ltyp  krlm '
!    &              //' krl1 kvar2d kk2d  nn2d  nvrl'
      nn2d = 0
      kk2d = 0
      do 470 nv=1,nvmax
        if (ltyp(nv).eq.99) goto 470
        if (krl1(nv).lt.krlmin) goto 470
        if (krlm(nv).ne.kmax) then
          kvar2d(nv) = 1 + min(kk2d,kv2dmx-krlm(nv))
          nn2d = nn2d + 1
          kk2d = kk2d + krlm(nv)
        endif
!       write( 6,'(A,8I6)') titcv(nv), nv, ltyp(nv), krlm(nv),
!    &       krl1(nv), kvar2d(nv), kk2d, nn2d, nvrl(nv)
 470  continue
      if (kk2d.gt.kv2dmx) then
        write(clio3_out_id,'(2(A,I6))')'ARRET, foroutp : Nb.Var.2D =', kk2d,
     &                   ' > kv2dmx =', kv2dmx
        stop
      elseif (nn2d.gt.nv2dmx) then
        write(clio3_out_id,'(2(A,I6))')'ARRET, foroutp : Nb.Var.2D =', nn2d,
     &                   ' > nv2dmx =', nv2dmx
        stop
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Definition des indices entourant le point (i,j) (pour le Lissage) :
!-----------------------------------------------------------------------

!     \  | i-1 |  i  | i+1 |
!    j+1 |  8  |  4  |  6  |
!     j  |  1  |  0  |  2  |
!    j-1 |  5  |  3  |  7  |
!---------------------------

      do 530 j=1,jmax
        jrn(j,3) = j - 1
        jrn(j,5) = j - 1
        jrn(j,7) = j - 1
        jrn(j,1) = j
        jrn(j,2) = j
        jrn(j,4) = j + 1
        jrn(j,6) = j + 1
        jrn(j,8) = j + 1
 530  continue
        jrn(1,3) = 1
        jrn(1,5) = 1
        jrn(1,7) = 1
        jrn(jmax,4) = jmax
        jrn(jmax,6) = jmax
        jrn(jmax,8) = jmax

      do 550 i=1,imax
        irn(i,1) = i - 1
        irn(i,5) = i - 1
        irn(i,8) = i - 1
        irn(i,3) = i
        irn(i,4) = i
        irn(i,2) = i + 1
        irn(i,6) = i + 1
        irn(i,7) = i + 1
 550  continue
      if (ltest.ge.1) then
        irn(1,1) = imax - 2
        irn(1,5) = imax - 2
        irn(1,8) = imax - 2
        irn(imax,2) = 3
        irn(imax,6) = 3
        irn(imax,7) = 3
      else
        irn(1,1) = 1
        irn(1,5) = 1
        irn(1,8) = 1
        irn(imax,2) = imax
        irn(imax,6) = imax
        irn(imax,7) = imax
      endif

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine foroutp -
      end
