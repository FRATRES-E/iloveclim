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

      MODULE COUPL2OCEAN_COM

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
      PUBLIC :: ec_co2oc

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Tentative list of variables communicated from atm to ocean:
!
!       winstua        --  Zonal      wind stress
!       winstva        --  Meridional wind stress
!       sumoswr        --  Shortwave radiation over the ocean
!       sumofwf        --  freshwater sum for the ocean: Evap - Rain - Snow - Runoff
!       flux_gris2ecb
!       wf_ice_sheet
!       sumiceb
!       sumohfx        --  net heat flux over ocean without solar flux
!       sumisno        --  net snowflux over sea-ice
!       sumiswr        --  short wave radiation over sea-ice
!       sumihfx        --  net heat flux over sea-ice without solar flux
!       sumicof        --  coefficients sensible-heat flux for implicit ice-temperature calculation
!       couptcc        --  proportion of cloud cover
!       sumuv10        --
!       sumpress       --  surface atmospheric pressure
!       sumu10         --  total (geostr+div) zonal wind vector at 10 meter heigth
!       sumv10         --  total (geostr+div) meridional wind vector at 10 meter heigth
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      CONTAINS


      SUBROUTINE ec_co2oc(ist)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** communicate data from coupler to ocean
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MODULE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use global_constants_mod, only: dblp=>dp

      use comcoup_mod, only: sumicof, sumihfx, sumisno, sumiswr, sumofwf, sumohfx, sumohsn, sumohss, sumoswr, sumpress, sumu10,   &
                             sumuv10, sumv10, winstua, winstva, couptcc, sumro
#if ( ICEBERG == 2 && ISM != 2  )
     use comcoup_mod,  only: sumiceb
#endif      

#if ( ISOOCN >= 1 )
      USE iso_param_mod, ONLY : ieau, neauiso, ieau18, delta, rsmow
#endif

#if ( CONSEAU == 1 || F_PALAEO_FWF == 1 )
      use varsCONSEAU_mod, only:  flux_gris2ecb, fwatconseau
#endif
#if ( F_PALAEO_FWF == 2 )
      use bloc0_mod, only:  fwatconseau
      use comland_mod, only: wf_ice_sheet
#endif
      use const_mod,   only: pi, zero

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use para0_mod,   only: imax, jmax
      use para_mod,    only:
      use bloc0_mod,   only: ks2, ihyster, tms, xfreshsv, temp_vent, temp_press
      use bloc_mod,    only: xang1, xang2
      use ice_mod,     only: fwat, ficebergn, ficebergs, tairox, tairoy, tenagy, tenagx, fsolcn, ratbqo, hnplbq, fsolg, ratbqg,   &
                             vabq, fevabq, fcscn, flecn, cloud, fwruno

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( ICEBERG  == 1 )
      use iceberg_mod, only: wind10_u, wind10_v
#elif ( ICEBERG == 2 && ISM != 2  )
      use iceberg_mod, only: iceflux_in_snow, wind10_u, wind10_v
#elif ( ICEBERG == 2 && ISM == 2 )
      use iceberg_mod, only: wind10_u, wind10_v
#endif      
     
      use dynami_mod,  only: area

      use coord_mod, only: zlatt, zlont

#if ( forced_winds == 1 )
      use offline_wind_forcing, only: get_winds
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       BY REFERENCE VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        INTEGER, INTENT(in) :: ist

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       LOCAL VARIABLES
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       REAL(kind=dblp), DIMENSION(imax,jmax) :: zfld ! dmr generic placeholder variable

#if ( ICEBERG > 0 )
       REAL(kind=dblp) :: w10u_old
#endif
       INTEGER(kind=ip):: ix,jy

#if ( ISOOCN >= 1 )
! variable to store the water and its isotopes from the atm to the ocn
! this variable is on the atmospheric water isotope count
       REAL(kind=dblp), DIMENSION(neauiso) :: variso
       CHARACTER(LEN=10) :: varisonm
       INTEGER(kind=ip) :: i,j,k, iz
#endif


!dmr ---  Former implicit variables
        REAL(KIND=dblp)   :: x1, x2, y1, y2, areat, sv
        INTEGER(KIND=ip)  :: irunlabel, iyear, iday2
!dmr ---
        
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       MAIN BODY OF THE ROUTINE
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


#if ( forced_winds == 0 )      
! *** zonal wind stress
      call at2oc(winstua,zfld)
      do ix = 1, imax
        do jy = 1, jmax
           tairox(ix,jy) = zfld(ix,jy)
           tenagx(ix,jy) = zfld(ix,jy)
        enddo
      enddo
! *** meridional wind stress
      call at2oc(winstva,zfld)
      do ix = 1, imax
        do jy = 1, jmax
           tairoy(ix,jy) = zfld(ix,jy)
           tenagy(ix,jy) = zfld(ix,jy)
        enddo
      enddo
 
#else
      iday2=mod(ist,360)
      if (iday2.eq.0) iday2 = 360
      tairox(:,:) = get_winds(0,iday2)
      tenagx(:,:) = tairox(:,:)
      tairoy(:,:) = get_winds(1,iday2)
      tenagy(:,:) =  tairoy(:,:)
#endif      
      
! *** short wave radiation over ocean
      call at2oc(sumoswr,zfld)
      do ix = 1, imax
        do jy = 1, jmax
          fsolcn(ix,jy) = zfld(ix,jy)
        enddo
      enddo
! *** net fresh water flux over ocean/sea-ice
#if ( ISOOCN >= 1 )
! here the loop is:
!     1 = water
!     2 = 1H_216O
!     3 = 1H_217O
!     4 = 1H_218O
!     5 = 2H1H16O

      DO iz=ieau, neauiso
      call at2oc(sumofwf(:,:,iz),zfld)
           areat=0.
           x1=280*pi/180
           x2=360*pi/180
           y1=20*pi/180
           y2=50*pi/180

      do ix = 1, imax
        do jy = 1, jmax
        ! [NOTA] in fwat, the last dimension is owatert
        !        namely:
        !        1 = water
        !        2 = salinity [never used outside this loop, water flux]
        !        3 -> 5 = isotopes
        !        So it is strangely coherent with atmosphere for the
        !         wrong reason ... ;-)
          fwat(ix,jy,iz) = -zfld(ix,jy)*86400.*1000.

!#if ( ISM == 1 )
!        if (flgism) then
!          fwat(ix,jy)=fwat(ix,jy)+(frwism(ix,jy)*86400.*1000.)
!        endif
!#endif
        enddo
      enddo

      enddo ! iso_loop
#else

      call at2oc(sumofwf,zfld)
           areat=0.
           x1=280*pi/180
           x2=360*pi/180
           y1=20*pi/180
           y2=50*pi/180
      do ix = 1, imax
        do jy = 1, jmax
          fwat(ix,jy) = -zfld(ix,jy)*86400.*1000.
#if ( ISM == 1 )
        if (flgism) then
          fwat(ix,jy)=fwat(ix,jy)+(frwism(ix,jy)*86400.*1000.)
        endif
#endif
        enddo
      enddo

#endif

! dmr Adding the runoff into the coupler, to be used as a diagnostic variable only!

      call at2oc(sumro,zfld)

      do ix = 1, imax
        do jy = 1, jmax
          fwruno(ix,jy) = -zfld(ix,jy)*86400.*1000.
        enddo
      enddo




!nb also for F_PALAEO_FWF
#if ( CONSEAU == 1 || F_PALAEO_FWF == 1 )
      call at2oc(flux_gris2ecb,zfld)
      do ix = 1, imax
        do jy = 1, jmax
!--- dmr&afq [ fwatconseau ] = m.s-1
!--- afq           fwatconseau(ix,jy) = -zfld(ix,jy)
           fwatconseau(ix,jy) = zfld(ix,jy)

!           if (fwatconseau(ix,jy) .ne. 0) then
!         write(*,*) 'fwatconseau dans ec_co2oc',ix,jy,fwatconseau(ix,jy)
!           endif
        enddo
      enddo

#endif

#if ( F_PALAEO_FWF == 2 )
      call at2oc(wf_ice_sheet,zfld)
      do ix = 1, imax
        do jy = 1, jmax
           fwatconseau(ix,jy) = zfld(ix,jy)
        enddo
      enddo
#endif


!For iceberg flux from excess snow summer over the year since Jan. 1st
#if ( ICEBERG == 2 && ISM != 2 )
      call at2oc(sumiceb,zfld)
      do ix = 1, imax
        do jy = 1, jmax
           iceflux_in_snow(ix,jy) = zfld(ix,jy)
        enddo
      enddo

#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr -- 04.10.11: Cleanup of unused loop when ihyster == 0
!dmr ---  Code below for use in hysteresis loops
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ecrespin:mettre la valeur du flux eau douce dans Sv
      if (ihyster.ne.0) then
         sv=xfreshsv(irunlabel+iyear-ihyster+1)
!  write(mouchard_id,*) 'sv xfreshsv',sv, irunlabel+iyear-ihyster+1

         do ix = 1, imax
            do jy = 1, jmax
               if ((zlont(ix,jy).le.x2).and.(zlont(ix,jy).ge.x1).and.(zlatt(ix,jy).le.y2).and.(zlatt(ix,jy).ge.y1)) then
                  areat=areat+(area(ix,jy)*tms(ix,jy,ks2))
               endif
            enddo
         enddo

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!     write(*,*)'total',areat
         do ix = 1, imax
            do jy = 1, jmax
               if ((zlont(ix,jy).le.x2).and.(zlont(ix,jy).ge.x1).and.(zlatt(ix,jy).le.y2).and.(zlatt(ix,jy).ge.y1)) then
#if ( ISOOCN >= 1 )
                  fwat(ix,jy,ieau)=fwat(ix,jy,ieau)+((sv*1.E+6/areat)*86400.*1000.)
                  WRITE(*,*) "Warning ! No isotopes in hysteresis !! "
#else
                  fwat(ix,jy)=fwat(ix,jy)+((sv*1.E+6/areat)*86400.*1000.)
#endif
               endif
            enddo
         enddo

      endif
!-----|--1--------2---------3---------4---------5---------6---------7-|

! *** net heat flux over ocean without solar flux
      call at2oc(sumohfx,zfld)
      do ix = 1, imax
        do jy = 1, jmax
          ratbqo(ix,jy) = zfld(ix,jy)
        enddo
      enddo
! *** net snowflux over sea-ice
#if ( ISOOCN >= 1 )

      ! [NOTA] water isotopes dimensionality is the same here as previously
      DO iz=ieau, neauiso
      call at2oc(sumisno(:,:,iz),zfld)
! --- Water isotopes conservation is problematic here!
! --- I do not see a proper way to do this except by suppressing the max
! ---  requirement below

! dmr --- UPDATED on 2015-02-11: need to do it anyhow if not water is not conserved at ALL !!

      do ix = 1, imax
        do jy = 1, jmax
          hnplbq(ix,jy,iz) = max (zero,zfld(ix,jy)*86400.*1000.)
        enddo
      enddo
! --- Water isotopes conservation?
      ENDDO

#else
      call at2oc(sumisno,zfld)
      do ix = 1, imax
        do jy = 1, jmax
          hnplbq(ix,jy) = max (zero,zfld(ix,jy)*86400.*1000.)
        enddo
      enddo
#endif
! *** short wave radiation over sea-ice
      call at2oc(sumiswr,zfld)
      do ix = 1, imax
        do jy = 1, jmax
          fsolg(ix,jy) = zfld(ix,jy)
        enddo
      enddo
! *** net heat flux over sea-ice without solar flux
      call at2oc(sumihfx,zfld)
      do ix = 1, imax
        do jy = 1, jmax
          ratbqg(ix,jy) = zfld(ix,jy)
        enddo
      enddo
! *** coefficients sensible-heat flux for implicit ice-temperature calculation
      call at2oc(sumicof,zfld)
      do ix = 1, imax
        do jy = 1, jmax
          vabq(ix,jy) = zfld(ix,jy)
        enddo
      enddo
! [ADD] cloud cover in ice model was zero all the time ...
      call at2oc(couptcc,zfld)
      do ix = 1, imax
        do jy = 1, jmax
          cloud(ix,jy) = zfld(ix,jy)
        enddo
      enddo

! *** 10 meter height Wind Magnitude (m/s) for iceberg model
!dmr ??? fixed 10-meters winds ? ?
      call at2oc(sumuv10,zfld)
      do ix = 1, imax
        do jy = 1, jmax
         temp_vent(ix,jy) = zfld(ix,jy)
        enddo
      enddo

! *** Surface pressure
      call at2oc(sumpress,zfld)
      do ix = 1, imax
        do jy = 1, jmax
         temp_press(ix,jy) = zfld(ix,jy)
        enddo
      enddo
! *** heat due to iceberg melting
      ficebergn=sumohsn
      ficebergs=sumohss
!     write(mouchard_id,*) 'iceberg',ficebergn,ficebergs

#if ( ICEBERG > 0 )
!dmr --- The two lines before can be used as to calculate the icebergs volume
!dmr --- and thus the heat from the excess snow stuff

! *** total (geostr+div) wind vector at 10 meter heigth
      call at2oc(sumu10,zfld)
      do ix = 1, imax
        do jy = 1, jmax
          wind10_u(ix,jy) = zfld(ix,jy)
        enddo
      enddo

      call at2oc(sumv10,zfld)
      do ix = 1, imax
        do jy = 1, jmax
          wind10_v(ix,jy) = zfld(ix,jy)
        enddo
      enddo

! *** Rotation of the wind vector
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      do ix = 1, imax
        do jy = 1, jmax
         w10u_old = wind10_u(ix,jy)
         wind10_u(ix,jy)  = w10u_old*xang2(ix,jy)-wind10_v(ix,jy)*xang1(ix,jy)
         wind10_v(ix,jy)  = wind10_v(ix,jy)*xang2(ix,jy)+w10u_old*xang1(ix,jy)
        enddo
      enddo


#endif /* on ICEBERG > 0 */


! Rotation of the wind stress

      do ix = 1, imax
        do jy = 1, jmax
         tairox(ix,jy)  = tenagx(ix,jy)*xang2(ix,jy)-tenagy(ix,jy)*xang1(ix,jy)
         tairoy(ix,jy)  = tenagy(ix,jy)*xang2(ix,jy)+tenagx(ix,jy)*xang1(ix,jy)
        enddo
      enddo
      do ix = 1, imax
        do jy = 1, jmax
         tenagx(ix,jy)  = tairox(ix,jy)*1.4
         tenagy(ix,jy)  = tairoy(ix,jy)*1.4
        enddo
      enddo


      do ix = 1, imax
        do jy = 1, jmax
          fevabq(ix,jy) = 0.0
          fcscn(ix,jy) = 0.0
          flecn(ix,jy) = 0.0

        enddo
      enddo

!      zsumz1=0.0
!      do j=1,jmax
!        do i=is1(j),is2(j)
!          zsumz1=zsumz1+area(i,j)*tms(i,j,ks2)
!        enddo
!      enddo
!      zsumz=0.0
!      do j=1,jmax
!        do i=is1(j),is2(j)
!          zsumz=zsumz+area(i,j)*tms(i,j,ks2)*fwat(i,j)
!        enddo
!      enddo
!      write(97,*) zsumz/zsumz1,zsumz1

      return

      END SUBROUTINE ec_co2oc

! --- dmr Interpolation routine from atmosphere to ocean
!         Uses "wa2o" = weights from atmosphere to ocean

      SUBROUTINE at2oc(fin,fout)

      use comcoup_mod, only: ijocn, ijatm, komax, inda2o, wa2o
      use comatm, only: nlat, nlon

      implicit none

      ! --- dmr input variable on atmospheric grid      
      real(kind=dblp), dimension(nlat,nlon), intent(in) :: fin(nlat,nlon)
      ! --- dmr output variable on oceanic grid
      real(kind=dblp), dimension(ijocn),     intent(out):: fout(ijocn)
      
      real(kind=dblp), dimension(ijatm) :: finf
      real(kind=dblp)                   :: zsum

      integer(kind=ip)                  :: ji,jk,i,j
      
      ji=0
      do i=1,nlat
        do j=1,nlon
          ji=ji+1
          finf(ji)=fin(i,j)
        enddo
      enddo

      do ji = 1, ijocn
      zsum = 0.
        do jk = 1, komax
          zsum = zsum + wa2o(ji,jk) * finf(inda2o(ji,jk))
        enddo
        fout(ji) = zsum
      enddo

      end subroutine at2oc

      END MODULE COUPL2OCEAN_COM

