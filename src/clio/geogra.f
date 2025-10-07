!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:47 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:47 CET 2009

      SUBROUTINE geogra(js1, js2, jeq, jdl1, jdl2, ijsdl, ijudl,
     &                  iberp, ibera)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!   Version 3x3 deg -- called by the routines  : defgrid, feedom, unibin.
!   initialisation of the indexes defining the geographical sectors:
!   zones (is/ezon), basins (is/ebas,js/ebas), straits(is/ehsf,js/ehsf)
! + separation between the 2 grids (jsep) & correspondance [jdl1, jdl2, ijsdl, ijudl]
! Input : js1,js2.
!-----
!  modif : 24/09/99

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use reper_mod

      use newunit_clio_mod, only: clio3_out_id

!! END_OF_USE_SECTION

!--- locales declaration

      integer(kind=ip):: js1, js2, jeq, jdl1, jdl2, ijsdl, ijudl, iberp
     >                 , ibera, ii, icycl, i, iic0, j, jj86, jj80
     >                 , jjberi, nz, nb, jjsuez, iis, iie, nvd, iidxx
     >                 , jj

      real(kind=dblp) :: xx, yy, xu0, yu0,  xxw, yyw, yyberi, yysuez
     >                 , xxi1, yyj1



!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Contants defining the  resolution and the axes .              |
!-----------------------------------------------------------------------

!--distance in longitude and latitude between 2 pts of the grid (in degree) :
      dlat  = 3.0
      dlong = 3.0
#if ( HRCLIO == 1 )
      dlat  = 1.5
      dlong = 1.5
#endif
!dmr [NOEQUI]
      dxwi = dlong
      dywj = dlat
!dmr [NOEQUI]

!--location of  the 1st scalar point (1,1) lat,long in degree :
      xlon1 =  25.5
      ylat1 = -79.5
#if ( HRCLIO == 1 )
      xlon1 =  24.75
      ylat1 = -78.75
#endif

!dmr [NOEQUI]
      xwi1 = xlon1
      ywj1 = ylat1
!dmr [NOEQUI]

!- index j (pt Velocity) corresponding to the geographical equator :
      jeq = nint(0.5 - ylat1 / dlat) + 1

!--2nd grille : grid North Atlantic + Arctic (AA)
!- geographical location (long E) of the North pole of the second grid
      xwpoln = untour - 111.0

!--In case of 2 connected grids :
!-    (xa  = longitudeAA = index i / ya  = latitudeAA = index j )
!- location of the 1st scalar point (1,jsepar) lat,long in degree (axes AA):
      xalon1 = 90. - 0.5 * dlat
      yalat1 = -46.5

!--In case of 1 grid composed of 2 grids:
!-    (xa  = longitudeAA = index j / ya  = latitudeAA = index i )
      dxaj =  dywj
      dyai = -dxwi
!- equator <=> longAA=90 E : xaj1 + (jeq-1)*dxaj - 0.5*dxaj = 90.
      xaj1 = 90.0 - dxaj * DFLOAT(jeq-1) + 0.5 * dxaj
!- xwpoln <=> latAA=90 N : yai1 + (ipoln-1)*dyai - 0.5*dyai = 90.
      ii = nint( (xwpoln - xwi1)/dxwi + 0.5  ) + 1
      yai1 = 90.0 - dyai * DFLOAT(ii-1) + 0.5 * dyai

!--verticaly (for the outputs = reverese order compred to the du model), in levels :
      dniv  = -1.0
      zniv1 = -1.0

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) Differents types of  basin .                      |
!-----------------------------------------------------------------------

!--Initialisation of "icl" (cyclic along "i")(jcl1 & jcl2 <- in "defgrid"):
      if (ltest.lt.1) then
!--bassin type rectangular :
        icycl = 1
        do i=-imax,imax+imax
          icl(i) = min(max(1,i),imax)
        enddo
      else
!--bassin with cyclic condition , (periode = icycl = imax - 2) :
        icycl = imax - 2
        iic0 = 2*icycl - 2
        do i=-imax,imax+imax
          icl(i) = 2 + mod(i+iic0,icycl)
        enddo
      endif

      if (ltest.lt.2) then
!--bassin one one  grid [ MODIFIE jeq !! ] :
        jeq  = jsepar - 1
        jdl1 = js2
        jdl2 = js1
      else
!--bassin on 2 grids, definition of the limits and  sectors that
!-- corresponds to each other
!- [ iw = ijsdl - ja et xw - xlpoln = 90. - ya ] thus :
        ijsdl = nint( (xwpoln+90.-xlon1-yalat1) / dlong ) + 1 + jsepar
        ijudl = ijsdl + 1
!- part of the equator (en long) located in the Atlantic Ocean : [ 60 W, 20 E ]
        jdl1 = nint( (untour+20.-xlon1) / dlong ) + 1
        jdl2 = nint( (untour-60.-xlon1) / dlong ) + 1
        jdl1 = ijsdl - jdl1
        jdl2 = ijsdl - jdl2
      endif

!--Separation between between the 2 grilds (if different metrics ) :
      if (jsepar.eq.jmax) then
        do i=1,imax
          jsep(i) = jsepar
        enddo
      else
!- limit (included in AA) between the 2 grids : x = 296.E - y.N, jeq<j<jsepar
        do i=1,imax
          xx = xlon1 + dlong * DFLOAT(i-1)
          yy = 296. - xx
          j = nint( (yy - ylat1) / dlat  ) + 1
          jsep(i) = min(max(j,jeq),jsepar+1)
        enddo
      endif

!--Separation Arctic + GIN Sea : velocity = jnorth ; scal. = jnorth - 1
      do i=1,imax
        jnorth(i) = jmax
      enddo
      if (ltest.ge.3) then
        xu0 = xlon1 - 1.5 * dlong
        yu0 = ylat1 - 1.5 * dlat
        jj86 = nint( (86.0 - yu0) / dlat)
        jj80 = nint( (80.0 - yu0) / dlat)
        jj80 = nint( (80.0 - yu0) / dlat)
        do i=1,imax
          xx = xu0 +  dlong * DFLOAT(i)
          if (xx.le.329.0) then
            yy = 86.0
            jnorth(i) = nint( (yy - yu0) / dlat )
          elseif (xx.le.335.0) then
            yy = 84.0 - xx + 330.0
            jnorth(i) = nint( (yy - yu0) / dlat )
          elseif (xx.le.341.0) then
            yy = 69.0 - xx + 336.0
            jnorth(i) = nint( (yy - yu0) / dlat )
          else
            yy = 65.0
            jnorth(i) = nint( (yy - yu0) / dlat )
          endif
        enddo
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


!- Bering Strait :
      if (ltest.eq.3 .and. iberp.ne.0) then
!- location velocity point : 170 W, 66 N  <-> grid AA : (200.9, 12.1)
        xxw = -170.
        yyw =  66.
!       yyrw = yyw * radian
!       xxrw = mod((untour + xxw - xwpoln), untour) * radian
!       yya = asin(cos(yyrw) * cos(xxrw)) / radian
!- (jberp=jsepar, jbera=jmax)
        iberp = nint( (untour + xxw - xwi1) / dxwi + 1.5 )
        ibera = nint( (12.1d0 - yai1) / dyai + 1.5 )
!- Modification of the location of the exit to the Arctic (because of B-grid) :
        ibera = ibera + 1
#if ( HRCLIO == 1 )
        iberp=113
        ibera=204
#endif
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) Definitions of the "Zones" .                                       |
!-----------------------------------------------------------------------

!--initialisation of thelimits of the de "zones" (=3 oceans) for the means
!  iszon(j,nz) : index i of the first box of longitude j of the zone nz
!  iezon(j,nz) : index i of the last box of longitude j of the zone nz

!- initialisation :
      do j=1,jmax
        iszon(j,0) = 2
        iezon(j,0) = imax - 1
        iezon(j,nbsmax) = imax - 1
        do nz=1,(nbsmax-1)
          iezon(j,nz) = 1
        enddo
      enddo
      do nz=0,nbsmax
        titzon(nz) = '(Global)  '
      enddo

      if(nbsmax.eq.3) then
        titzon(1) = '(Indian)  '
        titzon(2) = '(Pacific) '
        titzon(3) = '(Atlantic)'
!- South of Bering : 3 zones , 1 : Indian, 2 : Pacific, 3 Atlantic
      yyberi = 67.0
      jjberi = nint( (yyberi - ylat1) / dlat ) + 1
      do j=js1,jjberi
        yy = ylat1 + dlat * DFLOAT(j-1)
!- boundary Indian/Pacific : xx = f(yy) ; then iezon(-,1)
        if (yy.le.-31.0) then
          xx = 143.5
        elseif (yy.le.-6.0) then
          xx = 112.5 - yy
        elseif (yy.le.9.0) then
          xx = 103.5 - 0.25 * yy
        else
          xx = 105.0
        endif
        iezon(j,1) = nint( (xx-xlon1) / dlong - 0.5 ) + 1
!- boundary Pacific/Atlantic : iezon(-,2)
        if (yy.le.-60.0) then
          xx = 297.0
        elseif (yy.le.-54.0) then
          xx = 237.0 - yy
        elseif (yy.le.-42.0) then
          xx = 291.0
        elseif (yy.le.-36.0) then
          xx = 333.0 + yy
        elseif (yy.le.3.0) then
          xx = 297.0
        elseif (yy.le.22.5) then
          xx = 303.0 - yy - yy
        else
          xx = 258.0
        endif
        iezon(j,2) = nint( (xx-xlon1) / dlong - 0.5 ) + 1
       enddo
      endif

!--Initialisation of the beginning and end of zones from iezon:
      do j=1,jmax
        iszon(j,1)  = 2
        do nz=2,nbsmax
          iszon(j,nz) = iezon(j,nz-1) + 1
        enddo
      enddo

!--Initialisation of isbas(-,-) & iebas(-,0,-) [=limits E & W of the basins]
!-  form the limits of the zones (iszon,iezon) :
      do nb=0,nbsmax
       do j=1,jmax
         isbas(j,nb) = iszon(j,nb)
         iebas(j,0,nb) = iezon(j,nb)
       enddo
      enddo

!--Initialisation of the limits North & South for the separation between basins
!-  and the end of the zone (jezon) - latitude without communications between basins : Suez
      yysuez = 29.0
      jjsuez = nint( (yysuez - ylat1) / dlat ) + 1
      do nz=0,nbsmax
        jsbas(nz) = jjsuez
        jebas(nz) = js2
        jezon(nz) = jmax
      enddo
      jsbas(0) = js1

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      if (ltest.ge.2) then
!--Partial average North Atl.(put Northward of the Indian) :
!-  in order to conserve the separation :
      jezon(1) = jjsuez
      do j=jjsuez+1,js2
        iszon(j,1) = iszon(j,nbsmax)
        iezon(j,1) = iezon(j,nbsmax)

        yy = ylat1 + dlat * DFLOAT(j-1)
!- suppression of the Hudson Bay and Baffin bay
        if (yy.gt.68.and.yy.le.86.0) then
          xx = 335.0
        elseif (yy.le.50.0) then
!- suppression of the Golf of Mexico :
          xx = 193.0 + yy + yy
        else
          xx = 293.0
        endif
        iis = int( (xx-xlon1) / dlong ) + 1
        iszon(j,1) = max(iis, iszon(j,1))

!- suppression of the Mediterrannean :
        if (yy.le.43.0) then
          xx = 355.0
        else
          xx = 312.0 + yy
        endif
        iie = int( (xx-xlon1) / dlong ) + 1
        iezon(j,1) = min(iie, iezon(j,1))

       enddo
      endif

!--end fo the bloc defining the limits between zones .

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  4 ) Definitions of the Straits .                                      |
!-----------------------------------------------------------------------

      if (ltest.lt.1) then
!--sqaured Basin :
        nvhsf = 0
        ndhsf = 0
      else
!--

!--indexes ishsf,iehsf,jshsf,jehsf defining the straits (for the computation of Flux) :
      ndhsf = 0
      xxi1 = xlon1 - 0.5 * dlong
      yyj1 = ylat1 - 0.5 * dlat
!- Drake passage : long = 70 W , lat de 72 S a 50 S
      ndhsf = ndhsf + 1
      nvd = min(nhsfmx,ndhsf)
      tithsf(nvd) = ' Drake     E+W,E,W :'
      ishsf(nvd) = nint( (290. - xxi1 )/dlong ) + 1
      iehsf(nvd) = 0
      jshsf(nvd) = nint( (-72. - ylat1)/dlat  ) + 1
      jehsf(nvd) = nint( (-50. - ylat1)/dlat  ) + 1
!- Indonesian passage : long = 114 E , lat de 24 S a 0 N (et -> 24 N)
      ndhsf = ndhsf + 1
      nvd = min(nhsfmx,ndhsf)
      tithsf(nvd) = ' Indonesia E+W,E,W :'
      ishsf(nvd) = nint( (114. - xxi1 )/dlong ) + 1
      iehsf(nvd) = 0
      jshsf(nvd) = nint( (-24. - ylat1)/dlat  ) + 1
      jehsf(nvd) = nint( ( 24. - ylat1)/dlat  ) + 1
!- Bering Strait :
      ndhsf = ndhsf + 1
      nvd = min(nhsfmx,ndhsf)
      tithsf(nvd) = ' Bering    N+S,N,S :'
      ishsf(nvd) = max(iberp - 1, 1)
      iehsf(nvd) = min(iberp, imax)
      jshsf(nvd) = jsepar
      jehsf(nvd) = 0
!- Dubble Passage Greenland / Island / Scotland :
      yy = 66.0
!  location of Islande , width of the passages (grid GG) :
      xx = 340.0
      iidxx = nint(6.0 / dlong)
      ii = 1 + nint( (xx - xlon1) / dlong )
      jj = 1 + nint( (yy - yyj1 ) / dlat  )
!-
      ndhsf = ndhsf + 1
      nvd = min(nhsfmx,ndhsf)
      tithsf(nvd) = ' Dan.Str.  N+S,N,S :'
      ishsf(nvd) = ii - iidxx
      iehsf(nvd) = ii
      jshsf(nvd) = jj
      jehsf(nvd) = 0
!-
      ndhsf = ndhsf + 1
      nvd = min(nhsfmx,ndhsf)
      tithsf(nvd) = ' Scot/Isl  N+S,N,S :'
      ishsf(nvd) = ii + 1
      iehsf(nvd) = ii + 2 * iidxx
      jshsf(nvd) = jj
      jehsf(nvd) = 0
!- Dubble Passage Greenland / Spitzberg / Norvege :
!   modif Lat(GG)(Spitzberg) old=81.N -> new=83.N <= consistant with "evolu"
!- Fram Strait : lat(GG) = 83.N , long(GG) de 339.E a 346.E
      ndhsf = ndhsf + 1
      nvd = min(nhsfmx,ndhsf)
      tithsf(nvd) = ' Fram Str. N+S,N,S :'
      ishsf(nvd) = nint( (339. - xxi1) / dlong ) + 1
      iehsf(nvd) = nint( (346. - xxi1) / dlong ) + 1
      jshsf(nvd) = nint( ( 83. - yyj1) / dlat  ) + 1
      jehsf(nvd) = 0
!- North-West Passage (Arctic - Baffin): lat(GG)=84.N, long(GG) 324.E -> 330.E
      ndhsf = ndhsf + 1
      nvd = min(nhsfmx,ndhsf)
      tithsf(nvd) = ' Arct-Baff N+S,N,S :'
      ishsf(nvd) = nint( (324. - xxi1) / dlong ) + 1
      iehsf(nvd) = nint( (330. - xxi1) / dlong ) + 1
      jshsf(nvd) = nint( ( 84. - yyj1) / dlat  ) + 1
      jehsf(nvd) = 0
!- Spitzberg/Norway (Norw.Sea - Barents): lat(GG)=83.N, long(GG) 347.E -> 360.E
      ndhsf = ndhsf + 1
      nvd = min(nhsfmx,ndhsf)
      tithsf(nvd) = ' Spitz/Nor N+S,N,S :'
      ishsf(nvd) = nint( (347. - xxi1) / dlong ) + 1
      iehsf(nvd) = nint( (360. - xxi1) / dlong ) + 1
      jshsf(nvd) = nint( ( 83. - yyj1) / dlat  ) + 1
      jehsf(nvd) = 0
!- Gibraltar Strait  : long(GG) = 355.0 E , lat(GG) from 35 N a 43 N
      ndhsf = ndhsf + 1
      nvd = min(nhsfmx,ndhsf)
      tithsf(nvd) = ' Gibralt.  E+W,E,W :'
      ishsf(nvd) = nint( (355. - xxi1 )/dlong ) + 1
      iehsf(nvd) = 0
      jshsf(nvd) = nint( ( 35. - ylat1)/dlat  ) + 1
      jehsf(nvd) = nint( ( 43. - ylat1)/dlat  ) + 1
!- Passage Australia / New-Zeeland : long de 145 E a 176 E, lat = 38 S
      ndhsf = ndhsf + 1
      nvd = min(nhsfmx,ndhsf)
      tithsf(nvd) = ' Austr/N-Z N+S,N,S :'
      ishsf(nvd) = nint( (145. - xlon1) / dlong ) + 1
      iehsf(nvd) = nint( (176. - xlon1) / dlong ) + 1
      jshsf(nvd) = nint( (-38. - yyj1 ) / dlat  ) + 1
      jehsf(nvd) = 0
!- Nb. particular passages (-> computation of  Horiz.Stream-Funct. => evolu , ".x")
      nvhsf = 3
      if (imax.lt.120) ndhsf = nvhsf
      if (ndhsf.gt.nhsfmx) then
        write(clio3_out_id,*)
     &   'STOP, geogra, Pb.under-Dimensionned : ndhsf,nhsfmx =',
     &            ndhsf, nhsfmx
        stop
      endif

!-----
      endif

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- end of the routine geogra -
      end
