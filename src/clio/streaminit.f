!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:55 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:55 CET 2009

      SUBROUTINE streaminit
!
! A)
! i1(2)zon arrays determining basins in model grid coordinates
! South of Bering Strait: 3 zones
! nb=1: Indian, nb=2: Pacific, nb=3: Atlantic, nb=0: Global Ocean
! North of Bering Strait included in Atlantic part
! B)
! reading data file: jmaxl.dat ==> declared in 'reper .com'
! connects the model j-coordinates to true lattitudes and is used in computing
! contributions to streamFUNCTIONs for North Atlantic and Arctic Ocean:
! v.dx_(model j) and/or u.dy_(model j) ==> v_eff.dy_(true latitude)
! this part is independant of land-sea masks !
!

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use reper_mod
      
      use stream_mod, only: i1zon, i2zon
      
      use newunit_clio_mod, only: mouchard_id      
!! END_OF_USE_SECTION

c~       integer i1zon(jmax,0:nbsmax),i2zon(jmax,0:nbsmax)
c~       common/streamzon/ i1zon,i2zon


      character*30 fmtl
      character*70 line

      integer(kind=ip):: i, iml, j, jbering, jj1, jj2, jj3, jml, kml
     &                 , nb, nfrcl
      real(kind=dblp) :: spvl, xx, ybering, yy

      integer(kind=ip):: jmaxldat_id

! A) --------------------- A)
! code taken from geogra.f for defining borders of basins i1zon/i2zon
! may be different from arrays iszon and iezon of geogra.f !
      dlat  = 3.0
      dlong = 3.0
      xlon1 =  25.5
      ylat1 = -79.5
#if (HRCLIO == 1 )
      dlat  = 1.5
      dlong = 1.5
      xlon1 =  24.75
      ylat1 = -78.75
#endif
      do j=1,jmax
        i1zon(j,0) = 2
        i2zon(j,0) = imax - 1
        i2zon(j,nbsmax) = imax - 1
        do nb=1,nbsmax-1
          i2zon(j,nb) = 1
        enddo
      enddo
      ybering = 67.0
      jbering = nint( (ybering - ylat1) / dlat ) + 1
      do j=js1,jbering
        yy = ylat1 + dlat * DFLOAT(j-1)
!- border Indian/Pacific: xx = f(yy) ; then i2zon(-,1)
        if (yy.le.-31.0) then
          xx = 143.5
        elseif (yy.le.-6.0) then
          xx = 112.5 - yy
        elseif (yy.le.9.0) then
          xx = 103.5 - 0.25 * yy
        else
          xx = 105.0
        endif
        i2zon(j,1) = nint( (xx-xlon1) / dlong - 0.5 ) + 1
!- border Pacific/Atlantic: i2zon(-,2)
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
          xx = 291.0
          xx = 258.0
        endif
        i2zon(j,2) = nint( (xx-xlon1) / dlong - 0.5 ) + 1
      enddo
      do j=1,jmax
        i1zon(j,1)  = 2
        do nb=2,nbsmax
          i1zon(j,nb) = i2zon(j,nb-1) + 1
        enddo
      enddo
! B) --------------------- B)
! read file : jmaxl.dat
      do j=1,jmax
       do i=1,imax
        jmaxl(i,j)=j
       enddo
      enddo
      jmvlat = 0



      open(newunit=jmaxldat_id, file='inputdata/clio/jmaxl.dat'
     >   , status='old', err=480)
      read(jmaxldat_id,'(2A)',err=475) fmtl, line
      read(line,*) spvl, iml, jml, kml, nfrcl
      read(jmaxldat_id,*)
      read(jmaxldat_id,*)
      read(jmaxldat_id,'(A)',err=475) ttvlat
!      write(mouchard_id,*)'File jmaxl.dat'
!      write(mouchard_id,*) spvl, iml, jml, kml, nfrcl
      jj1 = max(1,-jml)
      jj2 = max(jml,1)
      jj3 = sign(1,jml)
      do j=jj1,jj2,jj3
        read(jmaxldat_id,fmtl,err=475) (jmaxl(i,j),i=1,iml)
!        write(mouchard_id,fmtl,err=475) (jmaxl(i,j),i=1,iml)
      enddo
      close(jmaxldat_id)
      jmvlat = abs(jml)
      write(mouchard_id,*)'File jmaxl.dat; jmvlat =', jmvlat
      return
 475  continue
      STOP ' File jmaxl.dat not properly read '
 480  continue
      STOP ' File jmaxl.dat not properly opened '
      return
      end
!------------------------------------------------------------
      SUBROUTINE streaminitout(nstreamout)



      use comemic_mod, only: nyears
      use para0_mod,   only: kmax
      use para_mod,    only: nbsmax
      use stream_mod, only: uuu, vvv, jmtt, kmtt, streamdat_id
     > , heatsaltfluxdat_id


#if ( COMATM == 0 )
      integer :: nyears,irunlabel,iatm,ilan,iice,iobtrop,
     &                   iobclin,nwrskip

      common /ec_timepar/ nyears,irunlabel,iatm,ilan,iice,iobtrop,
     &                   iobclin,nwrskip
#endif

      integer, intent(in) :: nstreamout


      integer :: streamctl_id, heatsaltfluxctl_id

!
c~       integer :: jmtt, kmtt
c~       parameter (jmtt=57,kmtt=kmax+1)
c~       real*4 uuu(jmtt,0:kmtt,0:nbsmax)
c~       real*4 vvv(jmtt,0:nbsmax+4)
c~       common/streamflu/ uuu,vvv



! A) ------------------------ c write str*.ctl for GrADS



!) for meridional overturning + (stream.ctl; stream.dat)
      open(newunit=streamctl_id,file='outputdata/ocean/stream.ctl')
      write(streamctl_id,fmt="('dset   ^stream.dat')")
      write(streamctl_id,fmt="('options big_endian')")
      write(streamctl_id,fmt="('undef ',1p,e12.4)") -1.0e20
      write(streamctl_id
     > ,fmt="('title Clio streamFUNCTIONs for merid. overturning')")
      write(streamctl_id
     > ,fmt="('xdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
#if ( HRCLIO == 0 )
      write(streamctl_id
     > ,fmt="('ydef ',i3,' linear ',2f7.2)") 57,-81.0,3.0
#elif ( HRCLIO == 1 )
      write(streamctl_id
     > ,fmt="('ydef ',i3,' linear ',2f7.2)") 57,-79.5,1.5
#endif
      write(streamctl_id,fmt="('zdef ',i3,' levels')") kmtt+1
#if ( HRCLIO == 0 )
      write(streamctl_id
     > ,fmt="(' 5500     5126.18 4385.22 3661.11 2963.25 2307.36')")
      write(streamctl_id
     > ,fmt="(' 1717.90  1225.11  850.19  588.88  415.07  299.29')")
      write(streamctl_id
     > ,fmt="('  219.86   163.28  121.52   89.75   64.96   45.20')")
      write(streamctl_id,fmt="('   29.17    15.98    5.00    0.00')")
#elif ( HRCLIO == 1 )
      write(streamctl_id
     > ,fmt="('5500  5154.36 4809.36 4465.10 4121.78')")
      write(streamctl_id
     > ,fmt="('3779.64  3439.02 3100.38 2764.42 2432.10')")
      write(streamctl_id
     > ,fmt="('2104.96  1785.36 1477.08 1186.20 921.90')")
      write(streamctl_id
     > ,fmt="('695.80   517.28  386.54  294.50  229.48')")
      write(streamctl_id
     > ,fmt="('182.18   146.40  118.30  95.44  76.30')")
      write(streamctl_id,fmt="('59.64 45.40 32.48 20.76 10.00')")
#endif
      if (nstreamout.eq.1) then
       write(streamctl_id
     > ,fmt="('tdef ',i5,' linear 1jan0001  30dy')") nyears*12
      else
       write(streamctl_id
     > ,fmt="('tdef ',i5,' linear 1jan0001  1YR')") nyears
      endif
      write(streamctl_id,fmt="('vars 4')")
      write(streamctl_id
     > ,fmt="('psig       22  99 StreamFUNCTION Global   Sv')")
      write(streamctl_id
     > ,fmt="('psii       22  99 StreamFUNCTION Indian   Sv')")
      write(streamctl_id
     > ,fmt="('psip       22  99 StreamFUNCTION Pacific  Sv')")
      write(streamctl_id
     > ,fmt="('psia       22  99 StreamFUNCTION Atlantic Sv')")
      write(streamctl_id,fmt="('endvars')")
      close(streamctl_id)
!) for heat flux and salt flux  + (heatsaltflux.ctl; heatsaltflux.dat)
      open(newunit=heatsaltfluxctl_id
     >   , file='outputdata/ocean/heatsaltflux.ctl')
      write(heatsaltfluxctl_id,fmt="('dset   ^heatsaltflux.dat')")
      write(heatsaltfluxctl_id,fmt="('options big_endian')")
      write(heatsaltfluxctl_id,fmt="('undef ',1p,e12.4)") -1.0e20
      write(heatsaltfluxctl_id
     > ,fmt="('title Clio streamFUNCTIONs for heat flux')")
      write(heatsaltfluxctl_id
     > ,fmt="('xdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
#if ( HRCLIO == 0 )
      write(heatsaltfluxctl_id
     > ,fmt="('ydef ',i3,' linear ',2f7.2)") 57,-81.0,3.0
#elif ( HRCLIO == 1 )
      write(heatsaltfluxctl_id
     > ,fmt="('ydef ',i3,' linear ',2f7.2)") 57,-81.0,3.0
#endif
      write(heatsaltfluxctl_id
     > ,fmt="('zdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
      if (nstreamout.eq.1) then
       write(heatsaltfluxctl_id
     > ,fmt="('tdef ',i5,' linear 1jan0001  30dy')") nyears*12
      else
       write(heatsaltfluxctl_id
     > ,fmt="('tdef ',i5,' linear 1jan0001  1YR')") nyears
      endif
      write(heatsaltfluxctl_id,fmt="('vars 8')")
      write(heatsaltfluxctl_id
     > ,fmt="('hg       1  99 Heat Flux Global   PW')")
      write(heatsaltfluxctl_id
     > ,fmt="('hi       1  99 Heat Flux Indian   PW')")
      write(heatsaltfluxctl_id
     > ,fmt="('hp       1  99 Heat Flux Pacific  PW')")
      write(heatsaltfluxctl_id
     > ,fmt="('ha       1  99 Heat Flux Atlantic PW')")
      write(heatsaltfluxctl_id
     > ,fmt="('sg       1  99 Salt Flux Global   psuSv')")
      write(heatsaltfluxctl_id
     > ,fmt="('si       1  99 Salt Flux Indian   psuSv')")
      write(heatsaltfluxctl_id
     > ,fmt="('sp       1  99 Salt Flux Pacific  psuSv')")
      write(heatsaltfluxctl_id
     > ,fmt="('sa       1  99 Salt Flux Atlantic psuSv')")
      write(heatsaltfluxctl_id,fmt="('endvars')")
      close(heatsaltfluxctl_id)


! B) ------------------------ c opening outputfiles



      open(newunit=streamdat_id,file='outputdata/ocean/stream.dat',
     &         form='unformatted',access='direct',
     &         recl=Size(uuu)*Kind(uuu(1,0,0)))
      open(newunit=heatsaltfluxdat_id,file='outputdata/ocean/heatsaltflux.dat',
     &         form='unformatted',access='direct',
     &         recl=Size(vvv)*Kind(vvv(1,0)))


      return
      end
!------------------------------------------------------------
