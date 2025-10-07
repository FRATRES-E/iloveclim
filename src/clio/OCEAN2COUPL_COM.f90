!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Copyright for initial code Hugues Goosse, Thierry Fichefet, https://www.elic.ucl.ac.be/modx/index.php?id=289
!      Branched from version 1.2

!   Further development:
!   Copyright 2023 Didier M. Roche (a.k.a. dmr)

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#include "choixcomposantes.h"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      MODULE OCEAN2COUPL_COM

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       USE global_constants_mod, ONLY: dblp=>dp, ip
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      IMPLICIT NONE
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   History
! dmr           Change from 0.0.0: Created a modular version from the initial no-module legacy code
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      CHARACTER(LEN=5), PARAMETER :: version_mod ="0.1.0"
      

      PRIVATE
      PUBLIC :: ec_oc2co, initseaalb
      
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Tentative list of variables communicated from ocean to atm:
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      
      
      CONTAINS
      
      
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_oc2co(ist)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** communicate oceanic data to the coupler
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!! START_OF_USE_SECTION

        use comemic_mod, only: fracto
        use comcoup_mod, only: couptcc
        use comatm,      only: nlat, nlon
        use comsurf_mod, only: fractn, nse, noc, tsurfn, albesn, albocef

#if ( ISOOCN >= 2 )
      USE isoatm_mod, ONLY: ratio_oceanatm
      USE iso_param_mod, ONLY: ieau17, ieaud
      USE para0_mod, only: owisostrt, owisostop, oc2atwisoindx
#endif

#if ( DOWNSTS == 1 )
      use vertDownsc_mod,    only: tsurfn_d
#endif

! afq -- this is to force an ice sheet albedo
!        if the ocean is covered by an ice sheet
#if ( IMSK == 1 )
      USE input_icemask, only:icemaskalb
      use comland_mod, only: albsnow
#endif

      use para0_mod, only: imax, jmax
      use para_mod,  only:
      use bloc0_mod, only: ks2, scal
      use bloc_mod , only:
      use ice_mod  , only: albq, ts, hgbq, hnbq
!! END_OF_USE_SECTION

      implicit none

!! START_OF_INCLUDE_SECTION

!! END_OF_INCLUDE_SECTION

      integer ix,jy,ist

#if ( ISOOCN >= 2 )
      INTEGER :: iz, iatmwiso
#endif

      real*8 zfld(imax,jmax),zalb,zalbp
      real*8 sst(nlat,nlon),albq1(nlat,nlon),ssi(nlat,nlon)
      real*8 seaalb(nlat,nlon),hic(nlat,nlon),hsn(nlat,nlon)
      real*8 tzero

      tzero=273.15d0

! *** SST
      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = scal(ix,jy,ks2,1)
        enddo
      enddo
      call oc2at(zfld,sst)

#if ( ISOOCN >= 2 )

! [NOTA] Loop here should be on the ocean part, hence owisostrt -> owisostop

      do iz = owisostrt, owisostop

! *** 17O, 18O, DH
      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = scal(ix,jy,ks2,iz)
        enddo
      enddo
      ! [NOTA] in the after-following call, the index should be oceanic
      !        need to convert the index from oceanic to atmospheric
      !        (only 17 -> 2H)
      call oc2atwisoindx(iz,iatmwiso)

      call oc2at(zfld,ratio_oceanatm(:,:,iatmwiso))

      enddo ! on iz, wisos ...

#endif
! *** ICE fraction
      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = 1.0-albq(ix,jy)
        enddo
      enddo
      call oc2at(zfld,albq1)
! *** STI
      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = ts(ix,jy)
        enddo
      enddo
      call oc2at(zfld,ssi)
! *** hic
      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = hgbq(ix,jy)
        enddo
      enddo
      call oc2at(zfld,hic)
! *** hsn
      do ix = 1, imax
        do jy = 1, jmax
          zfld(ix,jy) = hnbq(ix,jy)
        enddo
      enddo
      call oc2at(zfld,hsn)

      call detseaalb(seaalb)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      do ix = 1, nlon
        do jy = 1, nlat
          fractn(jy,ix,noc)=(1.0d0-albq1(jy,ix))*fracto(jy,ix)
          fractn(jy,ix,nse)=albq1(jy,ix)*fracto(jy,ix)

#if ( DOWNSTS == 1 )
          tsurfn_d(jy,ix,nse,:)=min(tzero,ssi(jy,ix))
          tsurfn_d(jy,ix,noc,:)=max(tzero-1.8d0,sst(jy,ix))
#endif
          tsurfn(jy,ix,nse)=min(tzero,ssi(jy,ix))
          tsurfn(jy,ix,noc)=max(tzero-1.8d0,sst(jy,ix))

          call ec_shine(tzero-0.15,tzero-0.25,tsurfn(jy,ix,nse),hic(jy,ix),hsn(jy,ix),zalb,zalbp)
          albesn(jy,ix,nse)=(1.0-couptcc(jy,ix))*zalbp + couptcc(jy,ix)*zalb
          albesn(jy,ix,noc)=albocef*seaalb(jy,ix)

#if ( IMSK == 1 )
          if (icemaskalb(jy,ix).gt.0.9) then
! afq -- we have an ice sheet over the ocean, force a snow albedo:
             albesn(jy,ix,nse)=albsnow(jy)
             albesn(jy,ix,noc)=albsnow(jy)
          endif
#endif

        enddo
      enddo


      return

      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE detseaalb(seaalb)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** calculates albedos as a FUNCTION of the time of the year, linearly
! *** interpolating between seasonal mean values
! *** albsea is the albedo of the open sea
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


        use comemic_mod, only: iday, imonth, albsea
        use comatm,      only: nlat, nlon

      implicit none


      integer i,j,id1,is1,is2
      real*8  albseaz(nlat),seaalb(nlat,nlon)
      real*8  sfrac

! *** interpolate between seasonal means

      id1=(imonth-1)*30+iday-14
      if (id1.lt.1) id1=id1+360

      is1=(id1+89)/90
      is2=is1+1
      if (is2.eq.5) is2=1

      sfrac=(id1-((is1-1)*90.+1.))/90.

      do j=1,nlat
        albseaz(j)=albsea(j,is1)+(albsea(j,is2)-albsea(j,is1))*sfrac
      enddo

      do j=1,nlon
        do i=1,nlat
          seaalb(i,j)=albseaz(i)
        enddo
      enddo

      return
      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE initseaalb
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comatm,      only: nlat
      use comemic_mod, only: albsea

      use global_constants_mod, only: ip

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      implicit none


      integer i,is
      
      integer(kind=ip) :: albedo_dat_id
      
! *** read climatological zonal mean albedos for each season

      open (newunit=albedo_dat_id,file='inputdata/clio/albedo.dat')

      read(albedo_dat_id,*)
      do i=1,nlat
        read(albedo_dat_id,45)(albsea(i,is),is=1,4)
      enddo
45    format(4(2x,f7.4))
      close(albedo_dat_id)
      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE oc2at(fin,fout)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use comcoup_mod
      use comatm, only: nlat, nlon

#if ( BATHY >= 2 )
!cnb
      use update_clio_bathy_tools, only: la_date
#endif

      implicit none

      integer ji,jk,i,j
      real*8  fin(ijocn),fout(nlat,nlon),zsum

      ji = 0
      do i=1,nlat
        do j=1,nlon
          zsum = 0.
          ji=ji+1
          do jk = 1, kamax
            zsum = zsum + wo2a(ji,jk) * fin(indo2a(ji,jk))
          enddo
          fout(i,j) = zsum
        enddo
      enddo

      end

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_shine(tfsn,tfsg,ts,hgbq,hnbq,zalb,zalbp)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** This SUBROUTINE computes albedo of snow-sea ice following SHINE &
! *** HENDERSSON-SELLERS [1985]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!  tfsn   : melting point temperature of snow (273.15 in CLIO)
!  tfsg   : melting point temperature of ice (273.05 in CLIO)
!  ts     : surface temperature
!  hgbq   : ice thickness
!  hnbq   : snow thickness
!  zalb   : ice/snow albedo for overcast sky
!  zalbp  : ice/snow albedo for clear sky



      USE comatm,  only:
      USE comphys, only: albice, alphd, alphdi, alphs, cgren

      implicit none

      real*8 tfsn,tfsg,ts,hgbq,hnbq,zalb,zalbp
      real*8 al
!     real*8 albin,albis,albice,alphd,alphdi,alphs,cgren



!driess    albin = 0.45
!driess    albis = 0.45
!driess    albice = 0.45
!     albice = 0.53
!     alphd  = 0.80
!     alphdi = 0.72
!     alphs  = 0.65
!driess    alphd  = 0.72
!driess    alphdi = 0.64
!driess    alphs  = 0.55
!     albin = 0.43
!     albis = 0.43
!     albice = 0.43
!     alphd  = 0.70
!     alphdi = 0.62
!     alphs  = 0.53
!     cgren = 0.04

!  albin: Albedo of melting ice in the arctic.
!  albis: Albedo of melting ice in the antarctic (SHINE
!         & HENDERSSON-SELLERS, 1985).
!  albice: Albedo of melting ice.
!  alphd : Albedo of snow (thickness > 0.05m)
!  alphdi: Albedo of thick bare ice
!  alphs : Albedo of melting snow
!  cgren: Correction of the snow or ice albedo to take into account
!         effects of cloudiness (GRENFELL & PEROVICH, 1984)

      if (hnbq.gt.0.0) then

! ***  Case of ice covered by snow.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (ts.lt.tfsn) then

! ***    Freezing snow.

          if (hnbq.gt.0.05) then
            zalbp = alphd
          else
            if (hgbq.gt.1.5) then
              zalbp = alphdi+(hnbq*(alphd-alphdi)/0.05)
            else if (hgbq.gt.1.0.and.hgbq.le.1.5) then
                   al = 0.472+2.0*(alphdi-0.472)*(hgbq-1.0)
            else if (hgbq.gt.0.05.and.hgbq.le.1.0) then
                   al = 0.2467+(0.7049*hgbq)-(0.8608*(hgbq*hgbq))+(0.3812*(hgbq*hgbq*hgbq))
            else
              al = 0.1+3.6*hgbq
            endif
            if (hgbq.le.1.5) zalbp=al+(hnbq*(alphd-al)/0.05)
          endif
        else
!
! ***    Melting snow.
!
          if (hnbq.ge.0.1) then
            zalbp = alphs
          else
            zalbp = albice+((alphs-albice)/0.1)*hnbq
          endif
        endif
      else
!
! *** Case of ice free of snow.
!
        if (ts.lt.tfsg) then
!
! ***    Freezing ice.
!
          if (hgbq.gt.1.5) then
            zalbp = alphdi
          else if (hgbq.gt.1..and.hgbq.le.1.5) then
                 zalbp = 0.472+2.*(alphdi-0.472)*(hgbq-1.)
          else if (hgbq.gt.0.05.and.hgbq.le.1.) then
                 zalbp = 0.2467+(0.7049*hgbq)-(0.8608*(hgbq*hgbq))+(0.3812*(hgbq*hgbq*hgbq))
          else
            zalbp = 0.1+3.6*hgbq
          endif
        else
!
! *** Melting ice.
!
          if (hgbq.gt.1.5) then
            zalbp = albice
          else if (hgbq.gt.1..and.hgbq.le.1.5)  then
                 zalbp = 0.472+(2.*(albice-0.472)*(hgbq-1.))
          else if (hgbq.gt.0.05.and.hgbq.le.1.) then
                 zalbp = 0.2467+0.7049*hgbq-(0.8608*(hgbq*hgbq))+(0.3812*(hgbq*hgbq*hgbq))
          else
            zalbp = 0.1+3.6*hgbq
          endif
        endif
      endif
      zalb=zalbp+cgren

!cnb test increase albedo sea ice
!      zalb=zalb*2
!      zalbp=zalbp*2
!      if (zalb.gt.1) then
!        zalb=1.0
!      endif
!      if (zalbp.gt.1) then
!        zalbp=1.0
!      endif

      return
      end
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      END MODULE OCEAN2COUPL_COM
