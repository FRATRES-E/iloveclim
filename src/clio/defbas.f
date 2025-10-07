!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:44 CET 2009

      SUBROUTINE defbas(kdefb,nn99,titbas)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Mise en place des limites de bassin (3 Oceans)  E (=isbas) et W (=iebas)
!  [ pour le calcul des Transports Meridien de masse et de scalaire ]
! a partir des limites de zones (is/ezon) et de la bathym (tms & tmu).
!- Mise en place des Coeff. pour Moy.Zonale(Vraie_Latitude)
! -- appelee par "class" apres "defgrid".
!  kdefb < 0 => iebas(k)=iebas(ks2) independant de k ; autre => depend de k
! 1er bassin : kdefb = 0 => Indien seul ; |kdefb| = 1 => Indien+Pacifique
!-----
!  modif : 29/03/99

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use reper_mod
      use newunit_clio_mod, only: clio3_out_id, mouchard_id
      
!! END_OF_USE_SECTION

!- dummy variables :
      character*(*) titbas(0:nbsmax)

!- variables locales :
      real(kind=dblp), dimension(jmax) :: ccxam, ssxa
c~       dimension ccxam(jmax), ssxa(jmax)
!     dimension vvv(imax,jmax)
      character(len=30) :: fmtl
      character(len=70) :: line


      integer(kind=ip):: kdefb, nn99, icycl, j, jj, i, nb, k, jjs, iis
     >                , n, nflag, iml, jml, kml, nfrcl, jj1, jj2, jj3
     >                , ii, nz, isav, llm, jj0

      real(kind=dblp) :: yw, yyalim, xw, degre, yy, spvl, ya, ccya, ssya
     >                , xa, ccyvr, yvr, aageo, bbgeo, unscyv, dyvr, yvrf
     >                , yyj, yyjf, ywu1, xwu1

      integer(ip) :: jmaxl_dat_id


!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  0 ) Initialisation .                                                |
!-----------------------------------------------------------------------

      degre = 1.0 / radian
      icycl = imax - 2
      jmvlat = 0

      do j=1,jmax
       yy = ylat1 + dlat * DFLOAT(j-1)
       jj = min(j,js2)
       do i=1,imax
        jmaxl(i,j) = j
        yvrlat(i,j) = yy
        jgeogr(i,j,0) = jj
        jgeogr(i,j,1) = jj
        rgeogr(i,j,0) = 1.0
        rgeogr(i,j,1) = 1.0
        ageogr(i,j) = 1.0
        bgeogr(i,j) = 0.0
!       vvv(i,j) = -9.0
       enddo
      enddo

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Definition des limites Est & West de chaque bassin .            |
!-----------------------------------------------------------------------

!--Mise en place des limites de bassin (3 Oceans)  E (=isbas) et W (=iebas)
      do nb=0,nbsmax
       titbas(nb) = titzon(nb)
       do j=1,jmax
        if (j.lt.ju1-1 .or. j.gt.ju2+1) iebas(j,0,nb) = isbas(j,nb)-1
        do k=1,kmax
          iebas(j,k,nb) = iebas(j,0,nb)
        enddo
       enddo
      enddo

      if_nbsmax_eq_3: if (nbsmax.eq.3) then

      if (abs(kdefb).eq.1) then
!-----
!--Reunion (Indien + Pacifique) = bassin 1  :
      do j=1,jmax
        isbas(j,1) = min(isbas(j,1), isbas(j,2))
        iebas(j,0,1) = max(iebas(j,0,1),iebas(j,0,2))
        do k=1,kmax
          iebas(j,k,1) = iebas(j,0,1)
        enddo
      enddo
      titbas(1) = '(Indo.Pac)'

      elseif (abs(kdefb).eq.2) then
!-----
!--Reunion (Indien + Atlantique) = bassin 1  :
      do j=1,jmax
        if (isbas(j,3).le.iebas(j,0,3)) then
          isbas(j,1) = isbas(j,3)-icycl
          if (isbas(j,1).gt.iebas(j,0,1))
     &      iebas(j,0,1) = iebas(j,0,3)-icycl
        endif
        do k=1,kmax
          iebas(j,k,1) = iebas(j,0,1)
        enddo
      enddo
      titbas(1) = '(Ind.+Atl)'
!-----
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Definition des limites Nord et Sud de chaque bassin .           |
!-----------------------------------------------------------------------

! Calcul Transp.Merid.Scal. sur chaque portion de bassin avec bord E&W etanche
!- ksud = 0 => Conserve toutes les portions / bord E & W impermeables
!- ksud = 1 => Conserve que la portion comprenant yysuez = 29.N
      ksud = 1

      if (ksud.eq.1) then
!--Definition des limites Sud de separation des bassins :
!-  a partir des frontieres E de bassin (iebas) et masque Terre / Mer .
      do nb=1,nbsmax
        jjs = js1
        iis = iebas(jsbas(nb),0,nb)
        do j=jsbas(nb),js1,-1
         if (iebas(j,0,nb).lt.iis) then
          do i=iebas(j,0,nb),iis-1
           if (tmu(icl(i+1),j,ks2).eq.one .and. jjs.eq.js1) jjs = j+1
          enddo
         elseif (iebas(j,0,nb).gt.iis) then
          do i=iis+1,iebas(j,0,nb)
           if (tmu(icl(i+1),j,ks2).eq.one .and. jjs.eq.js1) jjs = j+1
          enddo
         else
           if (tmu(icl(iis+1),j,ks2).eq.one.and.jjs.eq.js1) jjs = j+1
         endif
         iis = iebas(j,0,nb)
        enddo
        jsbas(nb) = jjs
      enddo
      jsbas(0) = jsbas(nbsmax)

!--A partir des 2 fins de frontiere -> limite Sud de bassin :
      if (abs(kdefb).eq.1) then
!- Indien <- Indo.Pac
        jsbas(3) = max(jsbas(2),jsbas(3))
        jsbas(2) = max(jsbas(1),jsbas(2))
        jsbas(1) = jsbas(3)
      elseif (abs(kdefb).eq.2) then
!- Indien <- Indo.Atl
        jsbas(3) = max(jsbas(2),jsbas(3))
        jsbas(2) = max(jsbas(1),jsbas(2))
        jsbas(1) = jsbas(2)
      else
        jsbas(3) = max(jsbas(2),jsbas(3))
        jsbas(2) = max(jsbas(1),jsbas(2))
        jsbas(1) = max(jsbas(0),jsbas(1))
      endif

      else
!--Definition des limites Sud de separation des bassins :
!- = indice du point precedant le 1er flux.N-S non nul
      loop_250: do nb=1,nbsmax
        loop_245: do j=ju1,ju2
          if (isbas(j,nb).gt.iebas(j,0,nb)) cycle loop_245
          loop_240: do i=isbas(j,nb),1+iebas(j,0,nb)
!         loop_240: do i=isbas(j,nb),iebas(j,0,nb)
!           if (tms(i,j,ks2).eq.one.and.tms(i,j-1,ks2).eq.one) exit loop_245
            if (tmu(icl(i),j,ks2).eq.one) exit loop_245
          enddo loop_240
        enddo loop_245
        jsbas(nb) = j - 1
      enddo loop_250
!-----
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Limite Sud pour bassin Global :
!- = indice du point precedant le 1er flux.N-S non nul
      loop_260outer: do j=ju1,ju2
       loop_260inner: do i=iu1(j),iu2(j)
!      loop_260inner: do 260 i=isf1(j),isf2(j)
!       if (tms(i,j,ks2).eq.one .and. tms(i,j-1,ks2).eq.one) exit loop_260outer
        if (tmu(icl(i),j,ks2).eq.one) exit loop_260outer
       enddo loop_260inner
      enddo loop_260outer
      jsbas(0) = j - 1

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--mise en place des limites Nord des bassins :
!- = indice du point suivant le 1er flux.N-S non nul
      jebas(0) = js1
      loop_280: do nb=1,nbsmax
        loop_275: do j=ju2,ju1,-1
          if (isbas(j,nb).gt.iebas(j,0,nb)) cycle loop_275
          loop_270: do i=isbas(j,nb),1+iebas(j,0,nb)
!         loop_270: do i=isbas(j,nb),iebas(j,0,nb)
!           if (tms(i,j,ks2).eq.one.and.tms(i,j-1,ks2).eq.one) exit loop_275
            if (tmu(icl(i),j,ks2).eq.one) exit loop_275
          enddo loop_270
        enddo loop_275
        jebas(nb) = min(j+1,jmax)
        jebas(0) = max(jebas(0),jebas(nb))
      enddo loop_280

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Supression des zones sans frontieres E/W impermeables :         |
!-----------------------------------------------------------------------

      if (kdefb.lt.0) then
!--Cas "kdefb < 0" : selon le masque en surface,
      k = ks2
      do nb=1,nbsmax
        do j=js1,js2
!- supression des segments de bassin "ouverts" :
         if (tmu(  icl(isbas(j,nb))  , j ,k).eq.one .or.
     &       tmu(  icl(isbas(j,nb))  ,j+1,k).eq.one .or.
     &       tmu(icl(iebas(j,0,nb)+1), j ,k).eq.one .or.
     &       tmu(icl(iebas(j,0,nb)+1),j+1,k).eq.one)
     &       iebas(j,k,nb) = isbas(j,nb) - 1
!        if((tms( icl(isbas(j,nb)-1) ,j,k).eq.one .and.
!    &       tms( icl(isbas(j,nb))   ,j,k).eq.one ).or.
!    &      (tms(icl(iebas(j,0,nb))  ,j,k).eq.one .and.
!    &       tms(icl(iebas(j,0,nb)+1),j,k).eq.one ))
!    &       iebas(j,k,nb) = isbas(j,nb) - 1
        enddo
      enddo
      do nb=1,nbsmax
       do k=ks1,ks2-1
        do j=js1,js2
         iebas(j,k,nb) = iebas(j,ks2,nb)
        enddo
       enddo
      enddo

      else
!- cas "kdefb >= 0" : selon le masque du niveau "k" correspondant,
      do nb=1,nbsmax
       do k=ks1,ks2
        do j=js1,js2
!- supression des segments de bassin "ouverts" :
         if (tmu(  icl(isbas(j,nb))  , j ,k).eq.one .or.
     &       tmu(  icl(isbas(j,nb))  ,j+1,k).eq.one .or.
     &       tmu(icl(iebas(j,0,nb)+1), j ,k).eq.one .or.
     &       tmu(icl(iebas(j,0,nb)+1),j+1,k).eq.one)
     &       iebas(j,k,nb) = isbas(j,nb) - 1
!        if((tms( icl(isbas(j,nb)-1) ,j,k).eq.one .and.
!    &       tms( icl(isbas(j,nb))   ,j,k).eq.one ).or.
!    &      (tms(icl(iebas(j,0,nb))  ,j,k).eq.one .and.
!    &       tms(icl(iebas(j,0,nb)+1),j,k).eq.one ))
!    &       iebas(j,k,nb) = isbas(j,nb) - 1
        enddo
       enddo
      enddo

!--Fin du traitement des zones "ouvertes" .
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4.a) Impression de control sur fichier "mouchard" .
!-----------------------------------------------------------------------

      if (nn99.eq.2) then
        write(mouchard_id,'(A)') 'Bassins : Limite Sud et Nord '
     &                //'(+ Lim.Nord Zones) jsbas,jebas & jezon :'
        do nb=0,nbsmax
          write(mouchard_id,'(A,I2,2(1X,2A,2I4))') '  nb=', nb, titbas(nb),
     &                   ' jsbas,jebas =', jsbas(nb), jebas(nb),
     &                 ' | '//titzon(nb), ' jezon =', jezon(nb)
        enddo
        write(mouchard_id,'(2A)') ' Limites E/W de bassins :',
     &   ' j, isbas(j,nb),(iebas(j,k,nb),k=0,ks2) avec nb=1,2,3'
        write(mouchard_id,'(6X,9(5X,A,I3,3X))') ('nb =', nb, nb=1,nbsmax)
        do j=jmax,1,-1
          write(mouchard_id,'(A,I4,9(A,I4,A,2I4))') 'j=', j,
     &    (' |',isbas(j,n),',',iebas(j,0,n),iebas(j,ks2,n),n=1,nbsmax)
        enddo
        write(mouchard_id,*)
        write(mouchard_id,'(A)') 'Details des sections "permeables" (iebas) :'
        do nb=1,nbsmax
         nflag = 0
         do j=jmax,1,-1
          if (iebas(j,0,nb).ne.iebas(j,ks2,nb)) then
           if (nflag.eq.0) then
             nflag = 1
             write(mouchard_id,'(A,I2,A)')
     &         '-> j,isbas, iebas(0 -> ks2) [ bassin nb=', nb, ' ]'
           endif
           write(mouchard_id,'(2(I4,A1),40I4)') j, ',', isbas(j,nb), ',',
     &                            (iebas(j,k,nb),k=0,ks2)
          endif
         enddo
        enddo
        write(mouchard_id,*)
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4.b) Lecture de "jmaxl" qui definit ligne_brisee_isolatitude .
!-----------------------------------------------------------------------

!--Lecture du fichier jmaxl.dat, indices "jmaxl" pour PSM & TRM Fct(Lat.vraie):
      open(newunit=jmaxl_dat_id, file='jmaxl.dat', status='old', err=480)
      read(jmaxl_dat_id,'(2A)',err=475) fmtl, line
      read(line,*) spvl, iml, jml, kml, nfrcl
      read(jmaxl_dat_id,*)
      read(jmaxl_dat_id,*)
      read(jmaxl_dat_id,'(A)',err=475) ttvlat
      jj1 = max(1,-jml)
      jj2 = max(jml,1)
      jj3 = sign(1,jml)
      do j=jj1,jj2,jj3
        read(jmaxl_dat_id,fmtl,err=475) (jmaxl(i,j),i=1,iml)
      enddo
      jmvlat = abs(jml)
 475  continue                      ! After error while opening file 'jmaxl.dat'
      close(jmaxl_dat_id)
 480  continue                      ! After error while reading from file 'jmaxl.dat'
      if (jmvlat.eq.0) ttvlat = ' '
      if (nn99.eq.2 .and. jmvlat.ge.1) write(mouchard_id,'(3A,I4)')
     &  ' - read IsoLat_lines on "jmaxl.dat": ttvlat=', ttvlat,
     &  ' ; jmvlat=', jmvlat

      endif if_nbsmax_eq_3

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  5 ) Coeff. pour Moy.Zonale(Vraie_Latitude), Centre des mailles (lpt=0)
!-----------------------------------------------------------------------

      if (ltest.le.2) return

!- calcul des COS & SIN une seule fois :
       do j=jeq,jmax
         yw = ywj1 + dywj * DFLOAT(j-1)
         ccxam(j)= sin(yw * radian)
         ssxa(j) = cos(yw * radian)
!        xa = 90.0 + yw
!        ccxa(j) = cos(xa * radian)
!        ssxa(j) = sin(xa * radian)
       enddo

!- Modif de jgeogr & rgeogr sur la partie de grille AA uniquement :
      yyalim = 90.0 - abs(dyai)
      loop_590: do i=1,imax
        xw = xwi1 + dxwi * DFLOAT(i-1)
        ya = 90.0 + xwpoln - xw
        if (abs(ya).ge.yyalim) cycle loop_590
        ccya = cos(ya * radian)
        ssya = sin(ya * radian)
        loop_580: do j=jsep(i),jmax
          yw = ywj1 + dywj * DFLOAT(j-1)
          xa = 90.0 + yw
!- Vraie Lat :
          ccyvr = ccya * ccxam(j)
          yvr = degre * asin(ccyvr)
          yvrlat(i,j) = yvr

!- matrice de rotation a appliquer aux vecteurs : angle(ugeo,ugrid)
          ccyvr = 1.0 - ccyvr*ccyvr
          if (ccyvr.le.epsil) then
!- point singulier (ex. Pole Nord) => pas de vitesse :
            aageo = 1.0
            bbgeo = 0.0
          else
!- point regulier :
            unscyv = 1. / sqrt(ccyvr)
            aageo = ssxa(j) * unscyv
            bbgeo = ssya * ccxam(j) * unscyv
          endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Evaluation de la taille de la maille ds la direction meridienne
          dyvr = dywj*max(abs(aageo)*cmy(i,j,0),abs(bbgeo)*cmx(i,j,0))
          yvrf = yvr - 0.5*dyvr
          if (dyvr.le.epsil) then
            write(clio3_out_id,'(A,2I3,1PE10.3)') 'ARRET defbas :'
     &          //' i,j,dyvr=',i,j,dyvr
            write(clio3_out_id,'(A,1P4E10.3)') 'cmx,y aa/bbgeo :',
     &        cmx(i,j,0), cmy(i,j,0), aageo, bbgeo
            stop
          endif
!- maille entre 2 lat.j : jgeogr = indice Inferieur ou egal correspondant :
          jj0  = nint( (yvr - ywj1)/dywj )
          yyj  = ywj1 + dywj * DFLOAT(jj0)
          yyjf = yyj - (0.5 + epsil) * dywj
          if (jj0.eq.js2 .or. yvrf.lt.yyjf) then
            jgeogr(i,j,0) = jj0
            yyj = ywj1 + dywj * DFLOAT(jj0 - 1)
          else
            jgeogr(i,j,0) = 1 + jj0
          endif
!- maille entre 2 lat.j : rgeogr = proportion de la maille sur "jgeogr" :
          rgeogr(i,j,0) = (yyj + 0.5*dywj - yvrf)/dyvr
          rgeogr(i,j,0) = min(one,rgeogr(i,j,0))

!--verification :
          if (rgeogr(i,j,0).lt.zero) then
             llm = 0
             write(clio3_out_id,'(A,6I4)') 'ARRET Defbas : Pb. Coeff < 0 !'
     &          //' i,j,llm,jgeogr=', i, j, llm, jgeogr(i,j,llm)
             write(clio3_out_id,'(A,1P3E13.6)') 'rgeogr , yvr, yvrf :',
     &                                rgeogr(i,j,llm),yvr,yvrf
             stop
!         elseif (nn99.eq.2 .and. i.eq.icheck .and. j.eq.jcheck) then
!           write(mouchard_id,'(A,3I4)') 'i,j,jgeogr =', i,j,jgeogr(i,j,0)
!           write(mouchard_id,'(A,1P4E10.3)') ' cmx,cmy aageo,bbgeo :',
!    &       cmx(i,j,0), cmy(i,j,0), aageo, bbgeo
!           write(mouchard_id,'(A,1P4E10.3)') ' dyvr,yvr, yyj, ddyy :',
!    &       dyvr, yvr, yyj, (yyj + 0.5*dywj - yvrf)/dyvr
          endif

        enddo loop_580
      enddo loop_590

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  6 ) Coeff. pour Moy.Zonale(Vraie_Latitude), Coin des mailles (lpt=3)
!-----------------------------------------------------------------------

!- calcul des COS & SIN une seule fois :
       ywu1 = ywj1 - 0.5 * dywj
       xwu1 = xwi1 - 0.5 * dxwi
       do j=jeq,jmax
         yw = ywu1 + dywj * DFLOAT(j-1)
         xa = 90.0 + yw
         ccxam(j)= sin(yw * radian)
         ssxa(j) = cos(yw * radian)
!        ccxa(j) = cos(xa * radian)
!        ssxa(j) = sin(xa * radian)
       enddo

!- Modif de j,r,a,b/geogr sur la partie de grille AA uniquement :
      yyalim = 90.0 - abs(dyai)
      loop_690: do i=1,imax
        xw = xwu1 + dxwi * DFLOAT(i-1)
        ya = 90.0 + xwpoln - xw
        if (abs(ya).ge.yyalim) cycle loop_690
        ccya = cos(ya * radian)
        ssya = sin(ya * radian)
        loop_680: do j=jsep(i),jmax
          yw = ywu1 + dywj * DFLOAT(j-1)
          xa = 90.0 + yw
!- Vraie Lat :
          yvr = degre * asin(ccya * ccxam(j))

!- matrice de rotation a appliquer aux vecteurs : angle(ugeo,ugridG)
          ccyvr = ccya * ccxam(j)
          ccyvr = 1.0 - ccyvr*ccyvr
          if (ccyvr.le.epsil) then
!- point singulier (ex. Pole Nord) => pas de vitesse :
            ageogr(i,j) = 0.0
            bgeogr(i,j) = 0.0
            dyvr = cmy(i,j,3)
          else
!- point regulier :
            unscyv = 1. / sqrt(ccyvr)
            ageogr(i,j) = ssxa(j) * unscyv
            bgeogr(i,j) = ssya * ccxam(j) * unscyv
            dyvr = abs(ageogr(i,j))*cmy(i,j,3)
          endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Evaluation de la taille de la maille ds la direction meridienne
          dyvr = dywj * max(dyvr, abs(bgeogr(i,j))*cmx(i,j,3) )
          yvrf = yvr - 0.5*dyvr
          if (dyvr.le.epsil) then
            write(clio3_out_id,'(A,2I3,1PE10.3)') 'ARRET defbas :'
     &          //' i,j,dyvr=',i,j,dyvr
            write(clio3_out_id,'(A,1P4E10.3)') 'cmx,y a/bgeogr :',
     &        cmx(i,j,3), cmy(i,j,3), ageogr(i,j), bgeogr(i,j)
            stop
          endif
!- maille entre 2 lat.j : jgeogr = indice Inferieur ou egal correspondant :
          jj0 = nint( (yvr - ywu1)/dywj )
          yyj = ywu1 + dywj * DFLOAT(jj0)
          yyjf = yyj - (0.5 + epsil) * dywj
          if (jj0.eq.ju2 .or. yvrf.lt.yyjf) then
            jgeogr(i,j,1) = jj0
            yyj = ywu1 + dywj * DFLOAT(jj0 - 1)
          else
            jgeogr(i,j,1) = 1 + jj0
          endif
!- maille entre 2 lat.j : rgeogr = proportion de la maille sur "jgeogr" :
          rgeogr(i,j,1) = (yyj + 0.5*dywj - yvrf)/dyvr
          rgeogr(i,j,1) = min(one,rgeogr(i,j,1))

!--verification :
          if (rgeogr(i,j,1).lt.zero) then
             llm = 1
             write(clio3_out_id,'(A,6I4)') 'ARRET Defbas : Pb. Coeff < 0 !'
     &          //' i,j,llm,jgeogr=', i, j, llm, jgeogr(i,j,llm)
             write(clio3_out_id,'(A,1P3E13.6)') 'rgeogr , yvr, yvrf :',
     &                                rgeogr(i,j,llm),yvr,yvrf
             stop
!         elseif (nn99.eq.2 .and. i.eq.icheck .and. j.eq.jcheck) then
!           write(mouchard_id,'(A,3I4)') 'i,j,jgeogr =', i,j,jgeogr(i,j,0)
!           write(mouchard_id,'(A,1P4E10.3)') ' cmx,cmy ageog,bgeog :',
!    &       cmx(i,j,0), cmy(i,j,0), ageogr(i,j), bgeogr(i,j)
!           write(mouchard_id,'(A,1P4E10.3)') ' dyvr,yvr, yyj, ddyy :',
!    &       dyvr, yvr, yyj, (yyj + 0.5*dywj - yvrf)/dyvr
          endif

        enddo loop_680
      enddo loop_690

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  7 ) Impression de control sur fichier "mouchard" .
!-----------------------------------------------------------------------

      ii = icheck
      jj = jcheck
      if (nn99.eq.2 .and. ii.ge.1 .and. ii.le.imax
     &              .and. jj.ge.1 .and. jj.le.jmax) then
        write(mouchard_id,'(2A,2I4)') 'Coeff. Moy.Zonale(Vraie_Latitude) :',
     &                  ' icheck,jcheck =', ii, jj
        write(mouchard_id,'(A,I4,3F10.6)') ' centre : jgeogr(0),rgeogr(0) =',
     &                jgeogr(ii,jj,0), rgeogr(ii,jj,0)
        write(mouchard_id,'(A,I4,3F10.6)') '   coin : jgeogr(1),rgeogr(1),'
     &  //'a/bgeogr', jgeogr(ii,jj,1), rgeogr(ii,jj,1)
     &              , ageogr(ii,jj), bgeogr(ii,jj)

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- Cherche tous les points contribuant a Z.Aver.(jcheck) :
      loop_750: do llm=0,1
       loop_740: do nz=0,nbsmax
         if (ii.lt.iszon(jj,nz) .or. ii.gt.iezon(jj,nz) ) cycle loop_740
         if (llm.eq.0) then
           write(mouchard_id,'(A,I3,A,I2,A)') 'Contribution a Z.Aver(jj=',
     &        jj, ', nb=', nz,' ), (centre) ; i,j,jgeogr,rgeogr :'
         else
           write(mouchard_id,'(A,I3,A,I2,A)') 'Contribution a Z.Aver(jj=',
     &        jj, ', nb=', nz,' ),  (coin)  ; i,j,jgeogr,rgeogr :'
         endif
         do j=js1,jezon(nz)
          isav = 0
          do i=iszon(j,nz),iezon(j,nz)
           if (jgeogr(i,j,llm).eq.jj.and.rgeogr(i,j,llm).eq.one) then
              if (isav.eq.0) isav = i
           else
            if (isav.ne.0) then
!- write isav -> i-1
              if (isav.eq.i-1) then
                write(mouchard_id,1710) i-1, j, jj, one
              else
                write(mouchard_id,1720) isav,' ->',i-1,j,jj,one
              endif
              isav = 0
            endif
            if (jgeogr(i,j,llm).eq.jj.and.rgeogr(i,j,llm).gt.epsil)
     &         write(mouchard_id,1710) i,j,jgeogr(i,j,llm), rgeogr(i,j,llm)
            if ( 1+jgeogr(i,j,llm).eq.jj .and.
     &           1.-rgeogr(i,j,llm).gt.epsil )
     &         write(mouchard_id,1710) i,j,jgeogr(i,j,llm), 1.-rgeogr(i,j,llm)
           endif
          enddo
          if (isav.ne.0 .and. isav.eq.i-1) then
            write(mouchard_id,1710) i-1, j, jj, one
          elseif (isav.ne.0) then
            write(mouchard_id,1720) isav, ' ->', i-1, j, jj, one
          endif
         enddo
       enddo loop_740
      enddo loop_750
!-----
!     llm=0
!     nz=0
!     do j=js1,jezon(nz)
!       do i=iszon(j,nz),iezon(j,nz)
!        if (jgeogr(i,j,llm).eq.jj.and.rgeogr(i,j,llm).gt.epsil)
!    &      vvv(i,j) = rgeogr(i,j,llm)
!        if (1+jgeogr(i,j,llm).eq.jj .and.
!    &      1.-rgeogr(i,j,llm).gt.epsil ) vvv(i,j)=1.-rgeogr(i,j,llm)
!       enddo
!     enddo
!     open(newunit=check_defbas_id,file='check',status='unknown')
!     read(check_defbas_id,*)
!     read(check_defbas_id,*)
!     read(check_defbas_id,*)
!     write(check_defbas_id,*) 'jj=', jj
!     do j=jmax,1,-1
!       write(check_defbas_id,'(125F7.3)') (vvv(i,j),i=1,imax)
!     enddo
!     write(check_defbas_id,*)
!     close(check_defbas_id)
!-----
      endif

      return
 1710 format(2X,I4,2X,2I4,F10.6)
 1720 format(I3,A,2I3,I4,F10.6)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine defbas -
      end
