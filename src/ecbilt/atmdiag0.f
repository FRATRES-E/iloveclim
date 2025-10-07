!dmr -- Added optional components choice - Fri Oct 28 11:02:20 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Fri Oct 28 11:02:20 CEST 2011
!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_atmout(istep)


      USE comatm, only: 
      use comemic_mod, only: iatm

#if ( ISOATM >= 1 )
      USE iso_param_mod, ONLY : ieau
#endif
      implicit none

      integer istep

      if ( mod(istep,iatm) .eq. 0) then
        call ec_selectout(istep)
      endif

      call ec_outamean(istep)
      call ec_outiocht(istep)

      return
      end subroutine ec_atmout

!23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine ec_selectout(istep)
!-----------------------------------------------------------------------
! *** this routine selects which kind of outputs should be written to output
! *** itel counts the number of days in the output interval
! *** meantype = 1: output monthly mean fields.
! *** meantype = 2: output seasonal mean fields.
! *** meantot = 1: computes whole period monthly or seasonal mean fields.
! *** meanyl  = 1: computes yearly monthly or seasonal mean fields.
! *** ioutdaily = 1 output instantanous fields.
! *** instcount: counter for output instantaneous fields.
! *** ixout: frequency for output instantanous fields in days.
! *** written by xueli wang.
!-----------------------------------------------------------------------
      use Atmosphere_Output

      USE comatm,  only: radius, om, nlat, nlon, iwater
      USE comdyn,  only: psi, pp, psig, vtot, qgpv, qprime, utot
      use comphys, only: cosnow, dysnow, tosnow, corain, dyrain, torain
     &                 , temp0g, temp2g, temp4g, dragane, roair, uv10rws
     &                 , q0, uvw10, cdragw
      use comdiag, only: newtotvar, minterv, meantype, meantot, ixout
     &                 , instcount, ioutyearly, itel, nvl, meanyl
     &                 , ioutdaily
      use comemic_mod, only: iatm, iseason, iday, imonth
      use comoutlocal_mod
      use comsurf_mod

#if ( ISOATM >= 1 )
      USE iso_param_mod, ONLY : ieau, ieau18, ieau17, ieaud
#endif
#if ( CLM_INDICES >= 1 )
      USE CLIMATE_INDICES_MOD, only: DAILYSTEP_FOR_CLIM_INDICES
#endif
      implicit none

      integer i,j,k,l,istep
      integer idmean,idstd
      real*8  facstr,costt,sintt,pfac,psifac,qpfac

      ivlevel(1) = 200
      ivlevel(2) = 500
      ivlevel(3) = 800
      itlevel(0) = 100
      itlevel(1) = 350
      itlevel(2) = 650
      itlevel(3) = 1000

!** some computation and unit transformation

      call ec_sptogg(psi(1,1),psig(1,1,1),pp)
      call ec_sptogg(psi(1,2),psig(1,1,2),pp)
      call ec_sptogg(psi(1,3),psig(1,1,3),pp)
      call ec_sptogg(qprime(1,1),qgpv(1,1,1),pp)
      call ec_sptogg(qprime(1,2),qgpv(1,1,2),pp)
      call ec_sptogg(qprime(1,3),qgpv(1,1,3),pp)

!  *** compute the precipitation, evaporation and runoffs in cm/year
      facstr = roair*uv10rws
      pfac=100.*3600.*24.*360.
      psifac=radius*radius*om
      qpfac=om
      do i=1,nlat
        costt=cos(dragane(i))
        sintt=sin(dragane(i))
        do j=1,nlon
          do k=1,3
            qgpv(i,j,k) = qgpv(i,j,k)*qpfac
            psig(i,j,k) = psig(i,j,k)*psifac
          enddo

#if ( ISOATM >= 1 )
c~           dyrain1(i,j) = (dyrain(i,j,ieau)+dysnow(i,j,ieau))*pfac
c~           corain1(i,j) = (corain(i,j,ieau)+cosnow(i,j,ieau))*pfac
c~           torain1(i,j) = (torain(i,j,ieau)+tosnow(i,j,ieau))*pfac
          iso18torain(i,j) = (torain(i,j,ieau18)+tosnow(i,j,ieau18))*pfac
          iso17torain(i,j) = (torain(i,j,ieau17)+tosnow(i,j,ieau17))*pfac
          isodtorain(i,j) = (torain(i,j,ieaud)+tosnow(i,j,ieaud))*pfac
c~           snow1(i,j)   = tosnow(i,j,ieau)*pfac
          iso18tosnow(i,j) = tosnow(i,j,ieau18)*pfac
          iso17tosnow(i,j) = tosnow(i,j,ieau17)*pfac
          isodtosnow(i,j) = tosnow(i,j,ieaud)*pfac
c~           evap1(i,j)   = evap(i,j,ieau)*pfac
c~ #else
#endif
          dyrain1(i,j) = (dyrain(i,j,iwater)+dysnow(i,j,iwater))*pfac
          corain1(i,j) = (corain(i,j,iwater)+cosnow(i,j,iwater))*pfac
          torain1(i,j) = (torain(i,j,iwater)+tosnow(i,j,iwater))*pfac
          snow1(i,j)   = tosnow(i,j ,iwater)*pfac
          evap1(i,j)   = evap(i,j,iwater)*pfac

          runofl1(i,j) = arunofl(i,j)*pfac
          runofo1(i,j) = arunofo(i,j)*pfac
          eminp1(i,j)  = evap1(i,j)-torain1(i,j)
          hesw(i,j)    = hesw0(i,j)+hesw1(i,j)+hesw2(i,j)+hesws(i,j)
          nlrads(i,j)  = ulrads(i,j)-dlrads(i,j)
          swrs1(i,j)   = hesws(i,j) !mohr
          lwrs1(i,j)   = dlrads(i,j) !mohr
          if (q0(i).gt.0d0) then
            albep(i,j)   = 1d0 - hesw(i,j)/q0(i)
          else
            albep(i,j)   = 1d0
          endif
          winstu1(i,j)= cdragw(i,j)*facstr*uvw10(i,j)*
     &                  (utot(i,j,3)*costt-vtot(i,j,3)*sintt)
          winstv1(i,j)= cdragw(i,j)*facstr*uvw10(i,j)*
     &                  (utot(i,j,3)*sintt+vtot(i,j,3)*costt)
        enddo
      enddo
!  *** compute temperature in K (?)
      do i=1,nlat
        do j=1,nlon
          tsurf1(i,j) = tsurf(i,j)+newtotvar(Surface_Temperature,6)
          temp4g1(i,j)=temp4g(i,j)+newtotvar(Surface_Temperature,6)
          temp2g1(i,j)=temp2g(i,j)+newtotvar(Surface_Temperature,6)
          tempsg1(i,j)=tempsg(i,j)+newtotvar(Surface_Temperature,6)
          temp0g1(i,j)=temp0g(i,j)+newtotvar(Surface_Temperature,6)
        enddo
      enddo

      if(meantype.eq.2) then
        if(istep.gt.(11*30*iatm)) then
          iseason = iseason + 1
          if (iseason.gt.90) iseason = 1
          itel = itel + 1
        endif
      else
        itel = itel + 1
      endif

      instcount = instcount + 1

!  *** write the instant data
c~       if (ioutdaily .eq. 1.and.instcount.eq.ixout) then
      if (ioutdaily .eq. 1) then
        call ec_outputinst
      endif

#if ( CLM_INDICES >= 1 )
      call DAILYSTEP_FOR_CLIM_INDICES()
#endif

      if (meantype .eq. 1) then
        minterv=30
        if (meantot .eq. 1) then
          call ec_outputmtl
        endif
        if (meanyl .eq. 1) then
          call ec_outputmyl
        endif
      endif

!      if (meantype .eq. 2 .and. istep .gt. 11*30*iatm) then
!        minterv=90
!        if (meantot .eq. 1) then
!          call ec_outputmtl
!        endif
!        if (meanyl .eq. 1) then
!          call ec_outputmyl
!        endif
!      endif

      if (ioutyearly .eq. 1) then
        call ec_outputyrl
      endif

      if (itel.eq.minterv) itel = 0
      if (instcount.eq.ixout) instcount = 0

      return
      end subroutine ec_selectout

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-|
      SUBROUTINE ec_outputinst
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-|
! *** this routine output instantanous fields.
! *** complete re-write Didier M. Roche, 2020-01-24
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-|

      USE Atmosphere_Output, only: open, close, Instantaneous_Data
     &  , output, write
     &  , Surface_Temperature         ! tsurf1
     &  , Temperature
     &  , Stratospheric_Temperature   ! temp0g1
     &  , Two_Meter_Temperature       ! tempsg1
     &  , Wind_U                      ! u200, u500, u800
     &  , Wind_V                      ! v200, v500, v800
     &  , Surface_Pressure            ! pground
     &  , Vertical_Pressure_Wind      ! omegg
     &  , U_Stress                    ! winstu1
     &  , V_Stress                    ! winstv1
     &  , Wind_at_10_Meter            ! uv10
     &  , Ageostrophic_Wind_U         ! udivg
     &  , Ageostrophic_Wind_V         ! vdivg
     &  , Stream_Function             ! psig
     &  , Velocity_Potential          ! chig
     &  , QG_Potential_Vorticity      ! qgpv
     &  , Geopotential_Height         ! geopg
     &  , Heating_Force               ! vhforg1, vhforg2
     &  , Potential_Vorticity_Forcing ! vforg1, vforg2, vforg3
     &  , Large_Scale_Precipitation   ! dyrain1
     &  , Convective_Precipitation    ! corain1
     &  , Total_Precipitation         ! torain1
     &  , Total_Snow_Fall             ! snow1
     &  , Surface_Evaporation         ! evap1
     &  , Evap_Minus_Precip           ! eminp1
     &  , Surface_Sensible_Heat_Flux  ! hflux
     &  , Surface_Latent_Heat_Flux    ! eflux
     &  , Surface_Solar_Radiation     ! hesws
     &  , Surface_Shortwave_Radiation ! swrs1
     &  , Surface_Longwave_Radiation  ! lwrs1
     &  , Top_Solar_Radiation         ! hesw
     &  , Top_Thermal_Radiation       ! ulrad0
     &  , Surface_Thermal_Radiation   ! nlrads
     &  , Surface_Albedo              ! alb2es
     &  , Planetary_Albedo            ! albep
     &  , Bottom_Moisture             ! abmoisg
     &  , Land_Snow_Depth             ! adsnow
     &  , Sea_Ice_Thickness           ! ahic
     &  , Land_Surface_Runoff         ! runofl1
     &  , Ocean_Surface_Runoff        ! runofo1
     &  , Specific_Humidity           ! rmoisg
     &  , Relative_Humidity           ! relhum
     &  , Drag_Coefficient_W          ! cdragw
     &  , Drag_Coefficient_V          ! cdragv
     &  , Total_Cloud_Cover           ! tccd
     &  , User_Assigned_T1            ! dumt1
     &  , User_Assigned_T2            ! dumt2
     &  , User_Assigned_U1            ! dumu1
     &  , User_Assigned_U2            ! dumu2
#if ( ISOATM >= 1 )
     &  , Total_Precip18              ! iso18torain
     &  , Total_Precip17              ! iso17torain
     &  , Total_Precipd               ! isodtorain
     &  , Total_Snow18                ! iso18tosnow
     &  , Total_Snow17                ! iso17tosnow
     &  , Total_Snowd                 ! isodtosnow
     &  , Spec_Hum_d18                ! rmoisg, ieau18
     &  , Spec_Hum_d17                ! rmoisg, ieau17
     &  , Spec_Hum_dd                 ! rmoisg, ieaud
#endif
#if ( ISOATM >= 2 )
     &  , Land_Snow_D18               ! adsnow, ieau18
     &  , Land_Snow_D17               ! adsnow, ieau17
     &  , Land_Snow_Dd                ! adsnow, ieaud
     &  , Bottom_Mois18               ! abmoisg, ieau18
     &  , Bottom_Mois17               ! abmoisg, ieau17
     &  , Bottom_Moisd                ! abmoisg, ieaud
#endif



      use comoutlocal_mod, only: tsurf1, temp0g1, temp2g1, temp4g1, tempsg1
     &  , tempsg1, winstu1, winstv1, dyrain1, corain1, torain1, snow1, evap1
     &  , eminp1, swrs1, lwrs1, hesw, nlrads, albep, runofl1, runofo1

      use comdyn,          only: u200, u500, u800, v200, v500, v800, omegg
     &  , udivg, vdivg, psig, chig, qgpv, geopg

      use comphys,         only: uv10, vhforg1, vhforg2, vforg1, vforg2, vforg3
     &  , rmoisg, relhum, cdragw, tccd, dumt1, dumt2, dumu1, dumu2

      use comsurf_mod,     only: pground, hflux, eflux, hesws, ulrad0, alb2es
     &  , abmoisg, adsnow, ahic, cdragv


      use comemic_mod, only: iday, imonth, iyear
      use comdiag,     only: newtotvar

      use comatm, only: nlat, nlon, iwater

#if ( ISOATM >= 1 )
      USE iso_param_mod, ONLY : ieau, ieau18, ieau17, ieaud
      USE comoutlocal_mod, only: iso18torain, iso18tosnow, iso17torain
     &                   , iso17tosnow, isodtorain, isodtosnow
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-|

      IMPLICIT NONE

      INTEGER, PARAMETER :: IOFlag = 1

      LOGICAL :: need_to_write=.TRUE.
      REAL(kind=8), dimension(1:nlat,1:nlon,1:4) :: mean, stddev
      INTEGER :: realyear,realmonth,realday

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8-|


      CALL open(Instantaneous_Data)

      IF (output(newtotvar(Surface_Temperature,IOFlag))) THEN
            mean(:,:,1) = tsurf1(1:nlat,1:nlon)
            CALL write(Surface_Temperature,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Temperature,IOFlag))) THEN
         mean(1:nlat,1:nlon,1) = temp0g1(1:nlat,1:nlon)
         mean(1:nlat,1:nlon,2) = temp2g1(1:nlat,1:nlon)
         mean(1:nlat,1:nlon,3) = temp4g1(1:nlat,1:nlon)
         mean(1:nlat,1:nlon,4) = tempsg1(1:nlat,1:nlon)
         CALL write(Temperature,mean(1:nlat,1:nlon,1:4))
      END IF

      IF (output(newtotvar(Stratospheric_Temperature,IOFlag))) THEN
         mean(1:nlat,1:nlon,1) = temp0g1(1:nlat,1:nlon)
         CALL write(Stratospheric_Temperature,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Two_Meter_Temperature,IOFlag))) THEN
         mean(1:nlat,1:nlon,1) = tempsg1(1:nlat,1:nlon)
         CALL write(Two_Meter_Temperature,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Wind_U,IOFlag))) THEN
         mean(1:nlat,1:nlon,1) = u200(1:nlat,1:nlon)
         mean(1:nlat,1:nlon,2) = u500(1:nlat,1:nlon)
         mean(1:nlat,1:nlon,3) = u800(1:nlat,1:nlon)
         CALL write(Wind_U,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(Wind_V,IOFlag))) THEN
         mean(1:nlat,1:nlon,1) = v200(1:nlat,1:nlon)
         mean(1:nlat,1:nlon,2) = v500(1:nlat,1:nlon)
         mean(1:nlat,1:nlon,3) = v800(1:nlat,1:nlon)
         CALL write(Wind_V,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(Surface_Pressure,IOFlag))) THEN
         mean(1:nlat,1:nlon,1) = pground(1:nlat,1:nlon)
         CALL write(Surface_Pressure,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Vertical_Pressure_Wind,IOFlag))) THEN
         mean(1:nlat,1:nlon,1:3)=omegg(1:nlat,1:nlon,1:3)
         CALL write(Vertical_Pressure_Wind,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(U_Stress,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=winstu1(1:nlat,1:nlon)
         CALL write(U_Stress,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(V_Stress,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=winstv1(1:nlat,1:nlon)
         CALL write(V_Stress,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Wind_at_10_Meter,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=uv10(1:nlat,1:nlon)
         CALL write(Wind_at_10_Meter,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_U,IOFlag))) THEN
         mean(1:nlat,1:nlon,1:3)=udivg(1:nlat,1:nlon,1:3)
         CALL write(Ageostrophic_Wind_U,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_V,IOFlag))) THEN
         mean(1:nlat,1:nlon,1:3)=vdivg(1:nlat,1:nlon,1:3)
         CALL write(Ageostrophic_Wind_V,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(Stream_Function,IOFlag))) THEN
         mean(1:nlat,1:nlon,1:3)=psig(1:nlat,1:nlon,1:3)
         CALL write(Stream_Function,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(Velocity_Potential,IOFlag))) THEN
         mean(1:nlat,1:nlon,1:3)=chig(1:nlat,1:nlon,1:3)
         CALL write(Velocity_Potential,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(QG_Potential_Vorticity,IOFlag))) THEN
         mean(1:nlat,1:nlon,1:3)=qgpv(1:nlat,1:nlon,1:3)
         CALL write(QG_Potential_Vorticity,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(Geopotential_Height,IOFlag))) THEN
         mean(1:nlat,1:nlon,1:3)=geopg(1:nlat,1:nlon,1:3)
         CALL write(Geopotential_Height,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(Heating_Force,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=vhforg1(1:nlat,1:nlon)
         mean(1:nlat,1:nlon,2)=vhforg2(1:nlat,1:nlon)
         CALL write(Heating_Force,mean(1:nlat,1:nlon,1:2))
      END IF

      IF (output(newtotvar(Potential_Vorticity_Forcing,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=vforg1(1:nlat,1:nlon)
         mean(1:nlat,1:nlon,2)=vforg2(1:nlat,1:nlon)
         mean(1:nlat,1:nlon,3)=vforg3(1:nlat,1:nlon)
         CALL write(Potential_Vorticity_Forcing,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(Large_Scale_Precipitation,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=dyrain1(1:nlat,1:nlon)
         CALL write(Large_Scale_Precipitation,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Convective_Precipitation,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=corain1(1:nlat,1:nlon)
         CALL write(Convective_Precipitation,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Total_Precipitation,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=torain1(1:nlat,1:nlon)
         CALL write(Total_Precipitation,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Total_Snow_Fall,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=snow1(1:nlat,1:nlon)
         CALL write(Total_Snow_Fall,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Surface_Evaporation,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=evap1(1:nlat,1:nlon)
         CALL write(Surface_Evaporation,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Evap_Minus_Precip,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=eminp1(1:nlat,1:nlon)
         CALL write(Evap_Minus_Precip,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Surface_Sensible_Heat_Flux,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=hflux(1:nlat,1:nlon)
         CALL write(Surface_Sensible_Heat_Flux,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Surface_Latent_Heat_Flux,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=eflux(1:nlat,1:nlon)
         CALL write(Surface_Latent_Heat_Flux,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Surface_Solar_Radiation,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=hesws(1:nlat,1:nlon)
         CALL write(Surface_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Surface_Shortwave_Radiation,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=swrs1(1:nlat,1:nlon)
         CALL write(Surface_Shortwave_Radiation,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Surface_Longwave_Radiation,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=lwrs1(1:nlat,1:nlon)
         CALL write(Surface_Longwave_Radiation,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Top_Solar_Radiation,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=hesw(1:nlat,1:nlon)
         CALL write(Top_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
      END IF


      IF (output(newtotvar(Top_Thermal_Radiation,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=ulrad0(1:nlat,1:nlon)
         CALL write(Top_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Surface_Thermal_Radiation,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=nlrads(1:nlat,1:nlon)
         CALL write(Surface_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Surface_Albedo,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=alb2es(1:nlat,1:nlon)
         CALL write(Surface_Albedo,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Planetary_Albedo,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=albep(1:nlat,1:nlon)
         CALL write(Planetary_Albedo,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Bottom_Moisture,IOFlag))) THEN
c~ #if ( ISOATM >= 2 )
c~          mean(1:nlat,1:nlon,1)=abmoisg(1:nlat,1:nlon,ieau)
c~ #else
         mean(1:nlat,1:nlon,1)=abmoisg(1:nlat,1:nlon,iwater)
c~ #endif
         CALL write(Bottom_Moisture,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Land_Snow_Depth,IOFlag))) THEN
c~ #if ( ISOATM >= 2 )
c~          mean(1:nlat,1:nlon,1)=adsnow(1:nlat,1:nlon,ieau)
c~ #else
         mean(1:nlat,1:nlon,1)=adsnow(1:nlat,1:nlon,iwater)
c~ #endif
         CALL write(Land_Snow_Depth,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Sea_Ice_Thickness,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=ahic(1:nlat,1:nlon)
         CALL write(Sea_Ice_Thickness,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Land_Surface_Runoff,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=runofl1(1:nlat,1:nlon)
         CALL write(Land_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Ocean_Surface_Runoff,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=runofo1(1:nlat,1:nlon)
         CALL write(Ocean_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Specific_Humidity,IOFlag))) THEN
c~ #if (ISOATM >= 1 )
c~          mean(1:nlat,1:nlon,1)=rmoisg(1:nlat,1:nlon,ieau)
c~ #else
         mean(1:nlat,1:nlon,1)=rmoisg(1:nlat,1:nlon, iwater)
c~ #endif
         CALL write(Specific_Humidity,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Relative_Humidity,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=relhum(1:nlat,1:nlon)
         CALL write(Relative_Humidity,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Drag_Coefficient_W,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=cdragw(1:nlat,1:nlon)
         CALL write(Drag_Coefficient_W,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Drag_Coefficient_V,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=cdragv(1:nlat,1:nlon)
         CALL write(Drag_Coefficient_V,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Total_Cloud_Cover,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=tccd(1:nlat,1:nlon)
         CALL write(Total_Cloud_Cover,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(User_Assigned_T1,IOFlag))) THEN
         mean(1:nlat,1:nlon,1:3)=dumt1(1:nlat,1:nlon,1:3)
         CALL write(User_Assigned_T1,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(User_Assigned_T2,IOFlag))) THEN
         mean(1:nlat,1:nlon,1:3)=dumt2(1:nlat,1:nlon,1:3)
         CALL write(User_Assigned_T2,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(User_Assigned_U1,IOFlag))) THEN
         mean(1:nlat,1:nlon,1:3)=dumu1(1:nlat,1:nlon,1:3)
         CALL write(User_Assigned_U1,mean(1:nlat,1:nlon,1:3))
      END IF

      IF (output(newtotvar(User_Assigned_U2,IOFlag))) THEN
         mean(1:nlat,1:nlon,1:3)=dumu2(1:nlat,1:nlon,1:3)
         CALL write(User_Assigned_U2,mean(1:nlat,1:nlon,1:3))
      END IF

#if ( ISOATM >= 1 )

      IF (output(newtotvar(Total_Precip18,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=iso18torain(1:nlat,1:nlon)
         CALL write(Total_Precip18,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Total_Precip17,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=iso17torain(1:nlat,1:nlon)
         CALL write(Total_Precip17,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Total_Precipd,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=isodtorain(1:nlat,1:nlon)
         CALL write(Total_Precipd,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Total_Snow18,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=iso18tosnow(1:nlat,1:nlon)
         CALL write(Total_Snow18,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Total_Snow17,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=iso17tosnow(1:nlat,1:nlon)
         CALL write(Total_Snow17,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Total_Snowd,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=isodtosnow(1:nlat,1:nlon)
         CALL write(Total_Snowd,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Spec_Hum_d18,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=rmoisg(1:nlat,1:nlon,ieau18)
         CALL write(Spec_Hum_d18,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Spec_Hum_d17,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=rmoisg(1:nlat,1:nlon,ieau17)
         CALL write(Spec_Hum_d17,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Spec_Hum_dd,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=rmoisg(1:nlat,1:nlon,ieaud)
         CALL write(Spec_Hum_dd,mean(1:nlat,1:nlon,1:1))
      END IF

#endif


#if ( ISOATM >= 2 )

      IF (output(newtotvar(Land_Snow_D18,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=adsnow(1:nlat,1:nlon,ieau18)
         CALL write(Land_Snow_D18,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Land_Snow_D17,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=adsnow(1:nlat,1:nlon,ieau17)
         CALL write(Land_Snow_D17,mean(1:nlat,1:nlon,1:1))
      END IF

       IF (output(newtotvar(Land_Snow_Dd,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=adsnow(1:nlat,1:nlon,ieaud)
         CALL write(Land_Snow_Dd,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Bottom_Mois18,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=abmoisg(1:nlat,1:nlon,ieau18)
         CALL write(Land_Snow_D18,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Bottom_Mois17,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=abmoisg(1:nlat,1:nlon,ieau17)
         CALL write(Land_Snow_D17,mean(1:nlat,1:nlon,1:1))
      END IF

      IF (output(newtotvar(Bottom_Moisd,IOFlag))) THEN
         mean(1:nlat,1:nlon,1)=abmoisg(1:nlat,1:nlon,ieaud)
         CALL write(Land_Snow_Dd,mean(1:nlat,1:nlon,1:1))
      END IF

#endif

      CALL close

      RETURN
      END SUBROUTINE ec_outputinst


!23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE ec_outputmtl
!-----------------------------------------------------------------------
! *** this routine calls another routine which computes the  whole period
! *** monthly or seasonal mean fields.
! *** written by xueli wang and nanne weber, april 1995.
!-----------------------------------------------------------------------
      USE Atmosphere_Output


      USE comatm,  only: iwater
      USE comdyn,  only: u200, u500, u800, v200, v500, v800, omegg
     &                 , udivg, vdivg, psig, chig, qgpv, geopg
      use comphys, only: uv10, vhforg1, vhforg2, vforg1, vforg2, vforg3
     &                 , rmoisg, relhum, cdragw, tccd, dumt1, dumt2
     &                 , dumu1, dumu2

      use comdiag
c~ Something like ...
c~       use comdiag, only: sxtsurf, sxtstrat, sxtemp2g, sxtemp4g, sxtempsg
c~      &                 , sxt2m, sxu200, sxu500, sxu800, sxv200, sxv500
c~      &                 , sxv800, sxpground, sxomeg1, sxomeg2, sxomeg3
c~      &                 , sxwinstu1, sxwinstv1, sxuv10, sxudivg1, sxudivg2
c~      &                 , sxudivg3, sxvdivg1, sxvdivg2, sxvdivg3, sxgrpsi1
c~      &                 , sxgrpsi2, sxgrpsi3, sxchi1, sxchi2, sxchi3
c~      &                 , sxqgpv1, sxqgpv2, sxqgpv3, sxgh1, sxgh2, sxgh3,
c~      &                 , sxvhforg1, sxvhforg2, sxvhforg3

      use comemic_mod, only: iday, imonth, iseason, iyear, nyears
      use comcoup_mod
      use comoutlocal_mod
      use comsurf_mod

#if (ISOATM >= 1 )
       USE iso_param_mod, ONLY : ieau, ieau18, ieau17, ieaud
#endif

      IMPLICIT NONE

      INTERFACE
         SUBROUTINE ec_totstat(xx,sumxx,sumxxsq,xxm,xxdev,compute)
         IMPLICIT NONE
         REAL*8, DIMENSION(:,:),   INTENT(in)    :: xx
         REAL*8, DIMENSION(:,:,:), INTENT(inout) :: sumxx, sumxxsq
         REAL*8, DIMENSION(:,:),   INTENT(out)   :: xxm, xxdev
         LOGICAL,                  INTENT(in)    :: compute
         END SUBROUTINE ec_totstat
      END INTERFACE

      INTEGER, PARAMETER :: IOFlag = 3

      LOGICAL :: need_to_write
      REAL(kind=8), dimension(1:nlat,1:nlon,1:4) :: mean, stddev
      INTEGER :: realyear,realmonth

      realyear=iyear
      realmonth=imonth
!mab: wrong!!!
!      if ( irunlabeld/=360 ) realyear=iyear+1
!      if ( irunlabeld/=360 ) realmonth=(irunlabeld/30)+1+(ndays/30)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!

      IF (itel == 0) RETURN

      IF (meantype == 1) THEN
         need_to_write = ( realyear == nyears .AND. imonth .GE. realmonth .AND.
     &                     iday == 30 )
         IF (.NOT. need_to_write) need_to_write = ( iyear == nyears .AND.
     &                     imonth .LT. realmonth .AND. iday == 30 )
         IF (need_to_write) CALL open(Total_Monthly_Means)
      ELSE
         need_to_write = ( realyear == nyears .AND. imonth .GE. realmonth .AND.
     &                     iseason == 90 )
         IF (.NOT. need_to_write) need_to_write = ( iyear == nyears .AND.
     &                     imonth .LT. realmonth .AND. iseason == 90 )
         IF (need_to_write) CALL open(Total_Seasonal_Means)
      END IF

      IF (output(newtotvar(Surface_Temperature,IOFlag))) THEN
         CALL ec_totstat(tsurf1(1:nlat,1:nlon),
     &        sxtsurf(1:nlat,1:nlon,1:12), sytsurf(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,1),    stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Temperature,IOFlag))) THEN
         CALL ec_totstat(temp0g1(1:nlat,1:nlon),
     &        sxtstrat(1:nlat,1:nlon,1:12),sytstrat(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(temp2g1(1:nlat,1:nlon),
     &        sxtemp2g(1:nlat,1:nlon,1:12),sytemp2g(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,2),   stddev  (1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(temp4g1(1:nlat,1:nlon),
     &        sxtemp4g(1:nlat,1:nlon,1:12),sytemp4g(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,3),   stddev  (1:nlat,1:nlon,3),
     &        need_to_write)
         CALL ec_totstat(tempsg1(1:nlat,1:nlon),
     &        sxtempsg(1:nlat,1:nlon,1:12),sytempsg(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,4),   stddev  (1:nlat,1:nlon,4),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Temperature,mean(1:nlat,1:nlon,1:4))
         END IF
      END IF

      IF (output(newtotvar(Stratospheric_Temperature,IOFlag))) THEN
         CALL ec_totstat(temp0g1(1:nlat,1:nlon),
     &        sxtstrat(1:nlat,1:nlon,1:12), sytstrat(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),    stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
           CALL write(Stratospheric_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Two_Meter_Temperature,IOFlag))) THEN
         CALL ec_totstat(tempsg1(1:nlat,1:nlon),
     &        sxt2m(1:nlat,1:nlon,1:12),syt2m (1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Two_Meter_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Wind_U,IOFlag))) THEN
         CALL ec_totstat(u200(1:nlat,1:nlon),
     &        sxu200(1:nlat,1:nlon,1:12),syu200(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(u500(1:nlat,1:nlon),
     &        sxu500(1:nlat,1:nlon,1:12),syu500(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(u800(1:nlat,1:nlon),
     &        sxu800(1:nlat,1:nlon,1:12),syu800(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Wind_U,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Wind_V,IOFlag))) THEN
         CALL ec_totstat(v200(1:nlat,1:nlon),
     &        sxv200(1:nlat,1:nlon,1:12),syv200(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(v500(1:nlat,1:nlon),
     &        sxv500(1:nlat,1:nlon,1:12),syv500(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(v800(1:nlat,1:nlon),
     &        sxv800(1:nlat,1:nlon,1:12),syv800(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Wind_V,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Surface_Pressure,IOFlag))) THEN
         CALL ec_totstat(pground(1:nlat,1:nlon),
     &        sxpground(1:nlat,1:nlon,1:12),sypground(1:nlat,1:nlon,1:12),
     &        mean     (1:nlat,1:nlon,1),   stddev   (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Pressure,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Vertical_Pressure_Wind,IOFlag))) THEN
         CALL ec_totstat(omegg(1:nlat,1:nlon,1),
     &        sxomeg1(1:nlat,1:nlon,1:12),syomeg1(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(omegg(1:nlat,1:nlon,2),
     &        sxomeg2(1:nlat,1:nlon,1:12),syomeg2(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,2),   stddev (1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(omegg(1:nlat,1:nlon,3),
     &        sxomeg3(1:nlat,1:nlon,1:12),syomeg3(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,3),   stddev (1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Vertical_Pressure_Wind,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(U_Stress,IOFlag))) THEN
         CALL ec_totstat(winstu1(1:nlat,1:nlon),
     &        sxwinstu1(1:nlat,1:nlon,1:12),sywinstu1(1:nlat,1:nlon,1:12),
     &        mean     (1:nlat,1:nlon,1),   stddev   (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(U_Stress,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(V_Stress,IOFlag))) THEN
         CALL ec_totstat(winstv1(1:nlat,1:nlon),
     &        sxwinstv1(1:nlat,1:nlon,1:12),sywinstv1(1:nlat,1:nlon,1:12),
     &        mean     (1:nlat,1:nlon,1),   stddev   (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(V_Stress,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Wind_at_10_Meter,IOFlag))) THEN
         CALL ec_totstat(uv10(1:nlat,1:nlon),
     &        sxuv10(1:nlat,1:nlon,1:12),syuv10(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Wind_at_10_Meter,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_U,IOFlag))) THEN
         CALL ec_totstat(udivg(1:nlat,1:nlon,1),
     &        sxudivg1(1:nlat,1:nlon,1:12),syudivg1(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(udivg(1:nlat,1:nlon,2),
     &        sxudivg2(1:nlat,1:nlon,1:12),syudivg2(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,2),   stddev  (1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(udivg(1:nlat,1:nlon,3),
     &        sxudivg3(1:nlat,1:nlon,1:12),syudivg3(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,3),   stddev  (1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Ageostrophic_Wind_U,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_V,IOFlag))) THEN
         CALL ec_totstat(vdivg(1:nlat,1:nlon,1),
     &        sxvdivg1(1:nlat,1:nlon,1:12),syvdivg1(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(vdivg(1:nlat,1:nlon,2),
     &        sxvdivg2(1:nlat,1:nlon,1:12),syvdivg2(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,2),   stddev  (1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(vdivg(1:nlat,1:nlon,3),
     &        sxvdivg3(1:nlat,1:nlon,1:12),syvdivg3(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,3),   stddev  (1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Ageostrophic_Wind_V,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Stream_Function,IOFlag))) THEN
         CALL ec_totstat(psig(1:nlat,1:nlon,1),
     &        sxgrpsi1(1:nlat,1:nlon,1:12),sygrpsi1(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(psig(1:nlat,1:nlon,2),
     &        sxgrpsi2(1:nlat,1:nlon,1:12),sygrpsi2(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,2),   stddev  (1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(psig(1:nlat,1:nlon,3),
     &        sxgrpsi3(1:nlat,1:nlon,1:12),sygrpsi3(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,3),   stddev  (1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Stream_Function,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Velocity_Potential,IOFlag))) THEN
         CALL ec_totstat(chig(1:nlat,1:nlon,1),
     &        sxchi1(1:nlat,1:nlon,1:12),sychi1(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(chig(1:nlat,1:nlon,1),
     &        sxchi2(1:nlat,1:nlon,1:12),sychi2(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(chig(1:nlat,1:nlon,1),
     &        sxchi3(1:nlat,1:nlon,1:12),sychi3(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Velocity_Potential,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(QG_Potential_Vorticity,IOFlag))) THEN
         CALL ec_totstat(qgpv(1:nlat,1:nlon,1),
     &        sxqgpv1(1:nlat,1:nlon,1:12),syqgpv1(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(qgpv(1:nlat,1:nlon,2),
     &        sxqgpv2(1:nlat,1:nlon,1:12),syqgpv2(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,2),   stddev (1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(qgpv(1:nlat,1:nlon,3),
     &        sxqgpv3(1:nlat,1:nlon,1:12),syqgpv3(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,3),   stddev (1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(QG_Potential_Vorticity,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Geopotential_Height,IOFlag))) THEN
         CALL ec_totstat(geopg(1:nlat,1:nlon,1),
     &        sxgh1(1:nlat,1:nlon,1:12),sygh1 (1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(geopg(1:nlat,1:nlon,2),
     &        sxgh2(1:nlat,1:nlon,1:12),sygh2 (1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(geopg(1:nlat,1:nlon,3),
     &        sxgh3(1:nlat,1:nlon,1:12),sygh3 (1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Geopotential_Height,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Heating_Force,IOFlag))) THEN
         CALL ec_totstat(vhforg1(1:nlat,1:nlon),
     &        sxvhforg1(1:nlat,1:nlon,1:12),syvhforg1(1:nlat,1:nlon,1:12),
     &        mean     (1:nlat,1:nlon,1),   stddev   (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(vhforg2(1:nlat,1:nlon),
     &        sxvhforg2(1:nlat,1:nlon,1:12),syvhforg2(1:nlat,1:nlon,1:12),
     &        mean     (1:nlat,1:nlon,2),   stddev   (1:nlat,1:nlon,2),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Heating_Force,mean(1:nlat,1:nlon,1:2))
         END IF
      END IF

      IF (output(newtotvar(Potential_Vorticity_Forcing,IOFlag))) THEN
         CALL ec_totstat(vforg1(1:nlat,1:nlon),
     &        sxvforg1(1:nlat,1:nlon,1:12),syvforg1(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(vforg2(1:nlat,1:nlon),
     &        sxvforg2(1:nlat,1:nlon,1:12),syvforg2(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,2),   stddev  (1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(vforg3(1:nlat,1:nlon),
     &        sxvforg3(1:nlat,1:nlon,1:12),syvforg3(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,3),   stddev  (1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
         CALL write(Potential_Vorticity_Forcing,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Large_Scale_Precipitation,IOFlag))) THEN
         CALL ec_totstat(dyrain1(1:nlat,1:nlon),
     &        sxdyrain(1:nlat,1:nlon,1:12),sydyrain(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
           CALL write(Large_Scale_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Convective_Precipitation,IOFlag))) THEN
         CALL ec_totstat(corain1(1:nlat,1:nlon),
     &        sxcorain(1:nlat,1:nlon,1:12),sycorain(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Convective_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Precipitation,IOFlag))) THEN
         CALL ec_totstat(torain1(1:nlat,1:nlon),
     &        sxtorain(1:nlat,1:nlon,1:12),sytorain(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snow_Fall,IOFlag))) THEN
         CALL ec_totstat(snow1(1:nlat,1:nlon),
     &        sxsnow(1:nlat,1:nlon,1:12),sysnow(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Snow_Fall,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Evaporation,IOFlag))) THEN
         CALL ec_totstat(evap1(1:nlat,1:nlon),
     &        sxevap(1:nlat,1:nlon,1:12),syevap(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Evaporation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Evap_Minus_Precip,IOFlag))) THEN
         CALL ec_totstat(eminp1(1:nlat,1:nlon),
     &        sxeminp(1:nlat,1:nlon,1:12),syeminp(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Evap_Minus_Precip,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Sensible_Heat_Flux,IOFlag))) THEN
         CALL ec_totstat(hflux(1:nlat,1:nlon),
     &        sxhflux(1:nlat,1:nlon,1:12),syhflux(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
          CALL write(Surface_Sensible_Heat_Flux,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Latent_Heat_Flux,IOFlag))) THEN
         CALL ec_totstat(eflux(1:nlat,1:nlon),
     &        sxeflux(1:nlat,1:nlon,1:12),syeflux(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Latent_Heat_Flux,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Solar_Radiation,IOFlag))) THEN
         CALL ec_totstat(hesws(1:nlat,1:nlon),
     &        sxhesws(1:nlat,1:nlon,1:12),syhesws(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Top_Solar_Radiation,IOFlag))) THEN
         CALL ec_totstat(hesw(1:nlat,1:nlon),
     &        sxhesw(1:nlat,1:nlon,1:12),syhesw(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Top_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

! mohr
! be sure it is going as output for swrs flag... !mohr
      IF (output(newtotvar(Surface_Shortwave_Radiation,IOFlag))) THEN
         CALL ec_totstat(swrs1(1:nlat,1:nlon),
     &        sxswrs(1:nlat,1:nlon,1:12),syswrs(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Shortwave_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Longwave_Radiation,IOFlag))) THEN
         CALL ec_totstat(lwrs1(1:nlat,1:nlon),
     &        sxlwrs(1:nlat,1:nlon,1:12),sylwrs(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Longwave_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF


      IF (output(newtotvar(Surface_Albedo,IOFlag))) THEN
         CALL ec_totstat(alb2es(1:nlat,1:nlon),
     &        sxalbes(1:nlat,1:nlon,1:12),syalbes(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Albedo,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Planetary_Albedo,IOFlag))) THEN
         CALL ec_totstat(albep(1:nlat,1:nlon),
     &        sxalbep(1:nlat,1:nlon,1:12),syalbep(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Planetary_Albedo,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Thermal_Radiation,IOFlag))) THEN
         CALL ec_totstat(nlrads(1:nlat,1:nlon),
     &        sxnlrads(1:nlat,1:nlon,1:12),synlrads(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
           CALL write(Surface_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Top_Thermal_Radiation,IOFlag))) THEN
         CALL ec_totstat(ulrad0(1:nlat,1:nlon),
     &        sxulrad1(1:nlat,1:nlon,1:12),syulrad1(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Top_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Bottom_Moisture,IOFlag))) THEN
c~ #if ( ISOATM >= 2 )
c~          CALL ec_totstat(abmoisg(1:nlat,1:nlon,ieau),
c~      &        sxbmoisg(1:nlat,1:nlon,1:12),sybmoisg(1:nlat,1:nlon,1:12),
c~      &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
c~      &        need_to_write)
c~ #else
         CALL ec_totstat(abmoisg(1:nlat,1:nlon,iwater),
     &        sxbmoisg(1:nlat,1:nlon,1:12),sybmoisg(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
c~ #endif
         IF (need_to_write) THEN
            CALL write(Bottom_Moisture,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Land_Snow_Depth,IOFlag))) THEN
c~ #if ( ISOATM >= 2 )
c~          CALL ec_totstat(adsnow(1:nlat,1:nlon,ieau),
c~      &        sxdsnow(1:nlat,1:nlon,1:12),sydsnow(1:nlat,1:nlon,1:12),
c~      &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1),
c~      &        need_to_write)
c~ #else
         CALL ec_totstat(adsnow(1:nlat,1:nlon,iwater),
     &        sxdsnow(1:nlat,1:nlon,1:12),sydsnow(1:nlat,1:nlon,1:12),
     &        mean   (1:nlat,1:nlon,1),   stddev (1:nlat,1:nlon,1),
     &        need_to_write)
c~ #endif
         IF (need_to_write) THEN
            CALL write(Land_Snow_Depth,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Sea_Ice_Thickness,IOFlag))) THEN
         CALL ec_totstat(ahic(1:nlat,1:nlon),
     &        sxhic(1:nlat,1:nlon,1:12),syhic (1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Sea_Ice_Thickness,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Land_Surface_Runoff,IOFlag))) THEN
         CALL ec_totstat(runofl1(1:nlat,1:nlon),
     &        sxrunofl(1:nlat,1:nlon,1:12),syrunofl(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Ocean_Surface_Runoff,IOFlag))) THEN
         CALL ec_totstat(runofo1(1:nlat,1:nlon),
     &        sxrunofo(1:nlat,1:nlon,1:12),syrunofo(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
             CALL write(Ocean_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Specific_Humidity,IOFlag))) THEN
c~ #if (ISOATM >= 1 )
c~          CALL ec_totstat(rmoisg(1:nlat,1:nlon,ieau),
c~      &        sxrmoisgw3(1:nlat,1:nlon,1:12),syrmoisgw3(1:nlat,1:nlon,1:12),
c~      &        mean      (1:nlat,1:nlon,1),   stddev    (1:nlat,1:nlon,1),
c~      &        need_to_write)
c~ #else
         CALL ec_totstat(rmoisg(1:nlat,1:nlon,iwater),
     &        sxrmoisgw3(1:nlat,1:nlon,1:12),syrmoisgw3(1:nlat,1:nlon,1:12),
     &        mean      (1:nlat,1:nlon,1),   stddev    (1:nlat,1:nlon,1),
     &        need_to_write)
c~ #endif
         IF (need_to_write) THEN
            CALL write(Specific_Humidity,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Relative_Humidity,IOFlag))) THEN
         CALL ec_totstat(relhum(1:nlat,1:nlon),
     &        sxrelhum(1:nlat,1:nlon,1:12),syrelhum(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Relative_Humidity,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Drag_Coefficient_W,IOFlag))) THEN
         CALL ec_totstat(cdragw(1:nlat,1:nlon),
     &        sxcdragw(1:nlat,1:nlon,1:12),sycdragw(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Drag_Coefficient_W,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Drag_Coefficient_V,IOFlag))) THEN
         CALL ec_totstat(cdragv(1:nlat,1:nlon),
     &        sxcdragv(1:nlat,1:nlon,1:12),sycdragv(1:nlat,1:nlon,1:12),
     &        mean    (1:nlat,1:nlon,1),   stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Drag_Coefficient_V,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Cloud_Cover,IOFlag))) THEN
         CALL ec_totstat(tccd(1:nlat,1:nlon),
     &        sxtcc(1:nlat,1:nlon,1:12),sytcc (1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Cloud_Cover,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_T1,IOFlag))) THEN
         CALL ec_totstat(dumt1(1:nlat,1:nlon,1),
     &        sxdt11(1:nlat,1:nlon,1:12),sydt11(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(dumt1(1:nlat,1:nlon,2),
     &        sxdt12(1:nlat,1:nlon,1:12),sydt12(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(dumt1(1:nlat,1:nlon,3),
     &        sxdt13(1:nlat,1:nlon,1:12),sydt13(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_T1,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_T2,IOFlag))) THEN
         CALL ec_totstat(dumt2(1:nlat,1:nlon,1),
     &        sxdt21(1:nlat,1:nlon,1:12),sydt21(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(dumt2(1:nlat,1:nlon,2),
     &        sxdt22(1:nlat,1:nlon,1:12),sydt22(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(dumt2(1:nlat,1:nlon,3),
     &        sxdt23(1:nlat,1:nlon,1:12),sydt23(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_T2,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_U1,IOFlag))) THEN
         CALL ec_totstat(dumu1(1:nlat,1:nlon,1),
     &        sxdu11(1:nlat,1:nlon,1:12),sydu11(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(dumu1(1:nlat,1:nlon,2),
     &        sxdu12(1:nlat,1:nlon,1:12),sydu12(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(dumu1(1:nlat,1:nlon,3),
     &        sxdu13(1:nlat,1:nlon,1:12),sydu13(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_U1,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_U2,IOFlag))) THEN
         CALL ec_totstat(dumu2(1:nlat,1:nlon,1),
     &        sxdu21(1:nlat,1:nlon,1:12),sydu21(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_totstat(dumu2(1:nlat,1:nlon,2),
     &        sxdu22(1:nlat,1:nlon,1:12),sydu22(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,2),   stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_totstat(dumu2(1:nlat,1:nlon,3),
     &        sxdu23(1:nlat,1:nlon,1:12),sydu23(1:nlat,1:nlon,1:12),
     &        mean  (1:nlat,1:nlon,3),   stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_U2,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF


#if ( ISOATM >= 1 )
      IF (output(newtotvar(Total_Precip18,IOFlag))) THEN
         CALL ec_totstat(iso18torain(1:nlat,1:nlon),
     &        sxtorain18(1:nlat,1:nlon,1:12),sytorain18(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Precip18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Precip17,IOFlag))) THEN
         CALL ec_totstat(iso17torain(1:nlat,1:nlon),
     &        sxtorain17(1:nlat,1:nlon,1:12),sytorain17(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Precip17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Precipd,IOFlag))) THEN
         CALL ec_totstat(isodtorain(1:nlat,1:nlon),
     &        sxtoraind(1:nlat,1:nlon,1:12),sytoraind(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Precipd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snow18,IOFlag))) THEN
         CALL ec_totstat(iso18tosnow(1:nlat,1:nlon),
     &        sxtosnow18(1:nlat,1:nlon,1:12),sytosnow18(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Snow18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snow17,IOFlag))) THEN
         CALL ec_totstat(iso17tosnow(1:nlat,1:nlon),
     &        sxtosnow17(1:nlat,1:nlon,1:12),sytosnow17(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Snow17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snowd,IOFlag))) THEN
         CALL ec_totstat(isodtosnow(1:nlat,1:nlon),
     &        sxtosnowd(1:nlat,1:nlon,1:12),sytosnowd(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Snowd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Spec_Hum_d18,IOFlag))) THEN
         CALL ec_totstat(rmoisg(1:nlat,1:nlon,ieau18),
     &        sxrmoisg18(1:nlat,1:nlon,1:12),syrmoisg18(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Spec_Hum_d18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Spec_Hum_d17,IOFlag))) THEN
         CALL ec_totstat(rmoisg(1:nlat,1:nlon,ieau17),
     &        sxrmoisg17(1:nlat,1:nlon,1:12),syrmoisg17(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Spec_Hum_d17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Spec_Hum_dd,IOFlag))) THEN
         CALL ec_totstat(rmoisg(1:nlat,1:nlon,ieaud),
     &        sxrmoisgd(1:nlat,1:nlon,1:12),syrmoisgd(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Spec_Hum_dd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

#endif
#if ( ISOATM >= 2 )
      IF (output(newtotvar(Land_Snow_D18,IOFlag))) THEN
         CALL ec_totstat(adsnow(1:nlat,1:nlon,ieau18),
     &        sxdsnow18(1:nlat,1:nlon,1:12),sydsnow18(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_D18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Land_Snow_D17,IOFlag))) THEN
         CALL ec_totstat(adsnow(1:nlat,1:nlon,ieau17),
     &        sxdsnow17(1:nlat,1:nlon,1:12),sydsnow17(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_D17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Land_Snow_Dd,IOFlag))) THEN
         CALL ec_totstat(adsnow(1:nlat,1:nlon,ieaud),
     &        sxdsnowd(1:nlat,1:nlon,1:12),sydsnowd(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_Dd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Bottom_Mois18,IOFlag))) THEN
         CALL ec_totstat(abmoisg(1:nlat,1:nlon,ieau18),
     &        sxbmoisg18(1:nlat,1:nlon,1:12),sybmoisg18(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_D18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Bottom_Mois17,IOFlag))) THEN
         CALL ec_totstat(abmoisg(1:nlat,1:nlon,ieau17),
     &        sxbmoisg17(1:nlat,1:nlon,1:12),sybmoisg17(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_D17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Bottom_Moisd,IOFlag))) THEN
         CALL ec_totstat(abmoisg(1:nlat,1:nlon,ieaud),
     &        sxbmoisgd(1:nlat,1:nlon,1:12),sybmoisgd(1:nlat,1:nlon,1:12),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_Dd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
#endif

      IF (need_to_write) CALL close

      RETURN
      END subroutine ec_outputmtl

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_totstat(x,sumx1,sumy1,xmean,xstd,compute)
! *** --------------------------------------------------------------------
! *** This routine computes whole period monthly or seasonal mean and standard
! *** deviation around it.
! *** all arrays are assumed to have the same shape.
! *** ------------------------------------------------------------------------


      USE comatm
      USE comdiag
      use comemic_mod, only: nyears, imonth, irunlabeld, iyear, iday
      use comcoup_mod


      IMPLICIT NONE
      REAL*8, DIMENSION(:,:),   INTENT(in)    :: x
      REAL*8, DIMENSION(:,:,:), INTENT(inout) :: sumx1, sumy1
      REAL*8, DIMENSION(:,:),   INTENT(out)   :: xmean, xstd
      LOGICAL,                  INTENT(in)    :: compute

      INTEGER :: k,ncase,nk
      INTEGER :: realyear,realmonth

      realyear=iyear
      realmonth=1
      if ( irunlabeld/=360 ) realyear=iyear+1
      if ( irunlabeld/=360 ) realmonth=(irunlabeld/30)+1
      IF (meantype == 1) THEN
         k     = imonth
         ncase = nyears*30
         nk    = 12
         IF (realyear == 1 .AND. imonth == realmonth .AND. iday == 1) THEN
            sumx1(:,:,1:nk)   = 0.0
            sumy1(:,:,1:nk)   = 0.0
         END IF
      ELSE                      ! IF (meantype == 2) THEN
         k     = Mod(imonth/3,4) + 1
         ncase = (nyears-1)*90
         nk    = 4
         IF (realyear == 1 .AND. imonth == 12 .AND. iday == 1) then
            sumx1(:,:,1:nk)   = 0.0
            sumy1(:,:,1:nk)   = 0.0
         END IF
      END IF

      sumx1(:,:,k) = sumx1(:,:,k) + x(:,:)
      sumy1(:,:,k) = sumy1(:,:,k) + x(:,:)**2

      IF (compute) THEN
         !IF ( iyear == nyears .AND. imonth .GE. realmonth .AND. iday == 30 ) ncase=ncase+30
         xmean(:,:) = sumx1(:,:,k)/ncase
         xstd (:,:) = sumy1(:,:,k) - sumx1(:,:,k)**2/ncase
         WHERE (xstd(:,:) <= 0.0)
            xstd(:,:) = 0.0
         ELSE WHERE
            xstd(:,:) = Sqrt(xstd(:,:)/(ncase - 1.0))
         END WHERE
      END IF

      RETURN
      END SUBROUTINE ec_totstat

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_outputmyl
!-----------------------------------------------------------------------
! *** this routine calls another routine which computes the yearly monthly
! *** mean fields.
! *** written by xueli wang, september 1995.
!-----------------------------------------------------------------------
      USE Atmosphere_Output


      USE comatm
      USE comdyn
      use comphys
      use comdiag
      use comemic_mod, only: iday
      use comcoup_mod
      use comoutlocal_mod
      use comsurf_mod

#if (ISOATM >= 1 )
       USE iso_param_mod, ONLY : ieau, ieau17, ieau18, ieaud
#endif

      IMPLICIT NONE

      INTERFACE
         SUBROUTINE ec_stat(xx,sumxx,sumxxsq,xxm,xxdev,compute)
         IMPLICIT NONE
         REAL*8, DIMENSION(:,:), INTENT(in)    :: xx
         REAL*8, DIMENSION(:,:), INTENT(inout) :: sumxx, sumxxsq
         REAL*8, DIMENSION(:,:), INTENT(out)   :: xxm, xxdev
         LOGICAL,                INTENT(in)    :: compute
         END SUBROUTINE ec_stat
      END INTERFACE

      INTEGER, PARAMETER :: IOFlag = 2

      LOGICAL :: need_to_write
      real(kind=8), dimension(1:nlat,1:nlon,1:4) :: mean, stddev

      IF (itel == 0) RETURN

!mab      need_to_write = ( itel == minterv )
      need_to_write = ( (mod(iday,30).eq.0) )

      IF (need_to_write) THEN
         IF (meantype == 1) THEN
            CALL open(Monthly_Means)
         ELSE
            CALL open(Seasonal_Means)
         END IF
      END IF


      IF (output(newtotvar(Surface_Temperature,IOFlag))) THEN
         CALL ec_stat(tsurf1(1:nlat,1:nlon),
     &        s1tsurf(1:nlat,1:nlon),  s2tsurf(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Temperature,IOFlag))) THEN
         CALL ec_stat(temp0g1(1:nlat,1:nlon),
     &        s1tstrat(1:nlat,1:nlon),  s2tstrat(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(temp2g1(1:nlat,1:nlon),
     &        s1temp2g(1:nlat,1:nlon),  s2temp2g(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(temp4g1(1:nlat,1:nlon),
     &        s1temp4g(1:nlat,1:nlon),  s2temp4g(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3),
     &        need_to_write)
         CALL ec_stat(tempsg1(1:nlat,1:nlon),
     &        s1tempsg(1:nlat,1:nlon),  s2tempsg(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,4),stddev  (1:nlat,1:nlon,4),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Temperature,mean(1:nlat,1:nlon,1:4))
         END IF
      END IF

      IF (output(newtotvar(Stratospheric_Temperature,IOFlag))) THEN
         CALL ec_stat(temp0g1(1:nlat,1:nlon),
     &        s1tstrat(1:nlat,1:nlon),  s2tstrat(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
           CALL write(Stratospheric_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Two_Meter_Temperature,IOFlag))) THEN
         CALL ec_stat(tempsg1(1:nlat,1:nlon),
     &        s1t2m(1:nlat,1:nlon),  s2t2m (1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Two_Meter_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Wind_U,IOFlag))) THEN
         CALL ec_stat(u200(1:nlat,1:nlon),
     &        s1u200(1:nlat,1:nlon),  s2u200(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(u500(1:nlat,1:nlon),
     &        s1u500(1:nlat,1:nlon),  s2u500(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(u800(1:nlat,1:nlon),
     &        s1u800(1:nlat,1:nlon),  s2u800(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Wind_U,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Wind_V,IOFlag))) THEN
         CALL ec_stat(v200(1:nlat,1:nlon),
     &        s1v200(1:nlat,1:nlon),  s2v200(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(v500(1:nlat,1:nlon),
     &        s1v500(1:nlat,1:nlon),  s2v500(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(v800(1:nlat,1:nlon),
     &        s1v800(1:nlat,1:nlon),  s2v800(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Wind_V,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Surface_Pressure,IOFlag))) THEN
         CALL ec_stat(pground(1:nlat,1:nlon),
     &        s1pground(1:nlat,1:nlon),  s2pground(1:nlat,1:nlon),
     &        mean     (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Pressure,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Vertical_Pressure_Wind,IOFlag))) THEN
         CALL ec_stat(omegg(1:nlat,1:nlon,1),
     &        s1omeg(1:nlat,1:nlon,1),s2omeg(1:nlat,1:nlon,1),
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(omegg(1:nlat,1:nlon,2),
     &        s1omeg(1:nlat,1:nlon,2),s2omeg(1:nlat,1:nlon,2),
     &        mean  (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(omegg(1:nlat,1:nlon,3),
     &        s1omeg(1:nlat,1:nlon,3),s2omeg(1:nlat,1:nlon,3),
     &        mean  (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Vertical_Pressure_Wind,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(U_Stress,IOFlag))) THEN
         CALL ec_stat(winstu1(1:nlat,1:nlon),
     &        s1winstu1(1:nlat,1:nlon),  s2winstu1(1:nlat,1:nlon),
     &        mean     (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(U_Stress,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(V_Stress,IOFlag))) THEN
         CALL ec_stat(winstv1(1:nlat,1:nlon),
     &        s1winstv1(1:nlat,1:nlon),  s2winstv1(1:nlat,1:nlon),
     &        mean     (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(V_Stress,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Wind_at_10_Meter,IOFlag))) THEN
         CALL ec_stat(uv10(1:nlat,1:nlon),
     &        s1uv10(1:nlat,1:nlon),  s2uv10(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Wind_at_10_Meter,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_U,IOFlag))) THEN
         CALL ec_stat(udivg(1:nlat,1:nlon,1),
     &        s1udivg(1:nlat,1:nlon,1),s2udivg(1:nlat,1:nlon,1),
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(udivg(1:nlat,1:nlon,2),
     &        s1udivg(1:nlat,1:nlon,2),s2udivg(1:nlat,1:nlon,2),
     &        mean   (1:nlat,1:nlon,2),stddev (1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(udivg(1:nlat,1:nlon,3),
     &        s1udivg(1:nlat,1:nlon,3),s2udivg(1:nlat,1:nlon,3),
     &        mean   (1:nlat,1:nlon,3),stddev (1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Ageostrophic_Wind_U,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_V,IOFlag))) THEN
         CALL ec_stat(vdivg(1:nlat,1:nlon,1),
     &        s1vdivg(1:nlat,1:nlon,1),s2vdivg(1:nlat,1:nlon,1),
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(vdivg(1:nlat,1:nlon,2),
     &        s1vdivg(1:nlat,1:nlon,2),s2vdivg(1:nlat,1:nlon,2),
     &        mean   (1:nlat,1:nlon,2),stddev (1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(vdivg(1:nlat,1:nlon,3),
     &        s1vdivg(1:nlat,1:nlon,3),s2vdivg(1:nlat,1:nlon,3),
     &        mean   (1:nlat,1:nlon,3),stddev (1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Ageostrophic_Wind_V,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Stream_Function,IOFlag))) THEN
         CALL ec_stat(psig(1:nlat,1:nlon,1),
     &        s1psi(1:nlat,1:nlon,1),s2psi (1:nlat,1:nlon,1),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(psig(1:nlat,1:nlon,2),
     &        s1psi(1:nlat,1:nlon,2),s2psi (1:nlat,1:nlon,2),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(psig(1:nlat,1:nlon,3),
     &        s1psi(1:nlat,1:nlon,3),s2psi (1:nlat,1:nlon,3),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Stream_Function,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Velocity_Potential,IOFlag))) THEN
         CALL ec_stat(chig(1:nlat,1:nlon,1),
     &        s1chi(1:nlat,1:nlon,1),s2chi (1:nlat,1:nlon,1),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(chig(1:nlat,1:nlon,2),
     &        s1chi(1:nlat,1:nlon,2),s2chi (1:nlat,1:nlon,2),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(chig(1:nlat,1:nlon,3),
     &        s1chi(1:nlat,1:nlon,3),s2chi (1:nlat,1:nlon,3),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Velocity_Potential,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(QG_Potential_Vorticity,IOFlag))) THEN
         CALL ec_stat(qgpv(1:nlat,1:nlon,1),
     &        s1qgpv(1:nlat,1:nlon,1),s2qgpv(1:nlat,1:nlon,1),
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(qgpv(1:nlat,1:nlon,2),
     &        s1qgpv(1:nlat,1:nlon,2),s2qgpv(1:nlat,1:nlon,2),
     &        mean  (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(qgpv(1:nlat,1:nlon,3),
     &        s1qgpv(1:nlat,1:nlon,3),s2qgpv(1:nlat,1:nlon,3),
     &        mean  (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(QG_Potential_Vorticity,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Geopotential_Height,IOFlag))) THEN
         CALL ec_stat(geopg(1:nlat,1:nlon,1),
     &        s1gh(1:nlat,1:nlon,1),s2gh  (1:nlat,1:nlon,1),
     &        mean(1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(geopg(1:nlat,1:nlon,2),
     &        s1gh(1:nlat,1:nlon,2),s2gh  (1:nlat,1:nlon,2),
     &        mean(1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(geopg(1:nlat,1:nlon,3),
     &        s1gh(1:nlat,1:nlon,3),s2gh  (1:nlat,1:nlon,3),
     &        mean(1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Geopotential_Height,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Heating_Force,IOFlag))) THEN
         CALL ec_stat(vhforg1(1:nlat,1:nlon),
     &        s1vhforg1(1:nlat,1:nlon),  s2vhforg1(1:nlat,1:nlon),
     &        mean     (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(vhforg2(1:nlat,1:nlon),
     &        s1vhforg2(1:nlat,1:nlon),  s2vhforg2(1:nlat,1:nlon),
     &        mean     (1:nlat,1:nlon,2),stddev   (1:nlat,1:nlon,2),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Heating_Force,mean(1:nlat,1:nlon,1:2))
         END IF
      END IF

      IF (output(newtotvar(Potential_Vorticity_Forcing,IOFlag))) THEN
         CALL ec_stat(vforg1(1:nlat,1:nlon),
     &        s1vforg1(1:nlat,1:nlon),  s2vforg1(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(vforg2(1:nlat,1:nlon),
     &        s1vforg2(1:nlat,1:nlon),  s2vforg2(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(vforg3(1:nlat,1:nlon),
     &        s1vforg3(1:nlat,1:nlon),  s2vforg3(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
         CALL write(Potential_Vorticity_Forcing,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Large_Scale_Precipitation,IOFlag))) THEN
         CALL ec_stat(dyrain1(1:nlat,1:nlon),
     &        s1dyrain(1:nlat,1:nlon),  s2dyrain(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
           CALL write(Large_Scale_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Convective_Precipitation,IOFlag))) THEN
         CALL ec_stat(corain1(1:nlat,1:nlon),
     &        s1corain(1:nlat,1:nlon),  s2corain(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Convective_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Precipitation,IOFlag))) THEN
         CALL ec_stat(torain1(1:nlat,1:nlon),
     &        s1torain(1:nlat,1:nlon),  s2torain(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snow_Fall,IOFlag))) THEN
         CALL ec_stat(snow1(1:nlat,1:nlon),
     &        s1snow(1:nlat,1:nlon),  s2snow(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
             CALL write(Total_Snow_Fall,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Evaporation,IOFlag))) THEN
         CALL ec_stat(evap1(1:nlat,1:nlon),
     &        s1evap(1:nlat,1:nlon),  s2evap(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Evaporation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Evap_Minus_Precip,IOFlag))) THEN
         CALL ec_stat(eminp1(1:nlat,1:nlon),
     &        s1eminp(1:nlat,1:nlon),  s2eminp(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Evap_Minus_Precip,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Sensible_Heat_Flux,IOFlag))) THEN
         CALL ec_stat(hflux(1:nlat,1:nlon),
     &        s1hflux(1:nlat,1:nlon),  s2hflux(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
          CALL write(Surface_Sensible_Heat_Flux,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Latent_Heat_Flux,IOFlag))) THEN
         CALL ec_stat(eflux(1:nlat,1:nlon),
     &        s1eflux(1:nlat,1:nlon),  s2eflux(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Latent_Heat_Flux,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Solar_Radiation,IOFlag))) THEN
         CALL ec_stat(hesws(1:nlat,1:nlon),
     &        s1hesws(1:nlat,1:nlon),  s2hesws(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

!mohr
      IF (output(newtotvar(Surface_Shortwave_Radiation,IOFlag))) THEN
         CALL ec_stat(swrs1(1:nlat,1:nlon),
     &        s1swrs(1:nlat,1:nlon),  s2swrs(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Shortwave_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Longwave_Radiation,IOFlag))) THEN
         CALL ec_stat(lwrs1(1:nlat,1:nlon),
     &        s1lwrs(1:nlat,1:nlon),  s2lwrs(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Longwave_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

!! I have to calculate the proper variables for s1swrs, s2swrs and
!! s1lwrs, s2lwrs, by the moment I keep this
!! to be sure that I get an output.
!!mohr

      IF (output(newtotvar(Top_Solar_Radiation,IOFlag))) THEN
         CALL ec_stat(hesw(1:nlat,1:nlon),
     &        s1hesw(1:nlat,1:nlon),  s2hesw(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Top_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Top_Thermal_Radiation,IOFlag))) THEN
         CALL ec_stat(ulrad0(1:nlat,1:nlon),
     &        s1ulrad1(1:nlat,1:nlon),  s2ulrad1(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Top_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Thermal_Radiation,IOFlag))) THEN
         CALL ec_stat(nlrads(1:nlat,1:nlon),
     &        s1nlrads(1:nlat,1:nlon),  s2nlrads(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
           CALL write(Surface_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Albedo,IOFlag))) THEN
         CALL ec_stat(alb2es(1:nlat,1:nlon),
     &        s1albes(1:nlat,1:nlon),  s2albes(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Surface_Albedo,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Planetary_Albedo,IOFlag))) THEN
         CALL ec_stat(albep(1:nlat,1:nlon),
     &        s1albep(1:nlat,1:nlon),  s2albep(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Planetary_Albedo,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Bottom_Moisture,IOFlag))) THEN
c~ #if ( ISOATM >= 2 )
c~          CALL ec_stat(abmoisg(1:nlat,1:nlon,ieau),
c~      &        s1bmoisg(1:nlat,1:nlon),  s2bmoisg(1:nlat,1:nlon),
c~      &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
c~      &        need_to_write)
c~ #else
         CALL ec_stat(abmoisg(1:nlat,1:nlon,iwater),
     &        s1bmoisg(1:nlat,1:nlon),  s2bmoisg(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
c~ #endif
         IF (need_to_write) THEN
            CALL write(Bottom_Moisture,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Land_Snow_Depth,IOFlag))) THEN
c~ #if ( ISOATM >= 2 )
c~          CALL ec_stat(adsnow(1:nlat,1:nlon,ieau),
c~      &        s1dsnow(1:nlat,1:nlon),  s2dsnow(1:nlat,1:nlon),
c~      &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
c~      &        need_to_write)
c~ #else
         CALL ec_stat(adsnow(1:nlat,1:nlon,iwater),
     &        s1dsnow(1:nlat,1:nlon),  s2dsnow(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write)
c~ #endif
         IF (need_to_write) THEN
            CALL write(Land_Snow_Depth,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Sea_Ice_Thickness,IOFlag))) THEN
         CALL ec_stat(ahic(1:nlat,1:nlon),
     &        s1hic(1:nlat,1:nlon),  s2hic (1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Sea_Ice_Thickness,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Land_Surface_Runoff,IOFlag))) THEN
         CALL ec_stat(runofl1(1:nlat,1:nlon),
     &        s1runofl(1:nlat,1:nlon),  s2runofl(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Ocean_Surface_Runoff,IOFlag))) THEN
         CALL ec_stat(runofo1(1:nlat,1:nlon),
     &        s1runofo(1:nlat,1:nlon),  s2runofo(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Ocean_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Specific_Humidity,IOFlag))) THEN
c~ #if (ISOATM >= 1 )
c~          CALL ec_stat(rmoisg(1:nlat,1:nlon,ieau),
c~      &        s1rmoisgw3(1:nlat,1:nlon),  s2rmoisgw3(1:nlat,1:nlon),
c~      &        mean      (1:nlat,1:nlon,1),stddev    (1:nlat,1:nlon,1),
c~      &        need_to_write)
c~ #else
         CALL ec_stat(rmoisg(1:nlat,1:nlon,iwater),
     &        s1rmoisgw3(1:nlat,1:nlon),  s2rmoisgw3(1:nlat,1:nlon),
     &        mean      (1:nlat,1:nlon,1),stddev    (1:nlat,1:nlon,1),
     &        need_to_write)
c~ #endif
         IF (need_to_write) THEN
            CALL write(Specific_Humidity,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Relative_Humidity,IOFlag))) THEN
         CALL ec_stat(relhum(1:nlat,1:nlon),
     &        s1relhum(1:nlat,1:nlon),  s2relhum(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Relative_Humidity,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Drag_Coefficient_W,IOFlag))) THEN
         CALL ec_stat(cdragw(1:nlat,1:nlon),
     &        s1cdragw(1:nlat,1:nlon),  s2cdragw(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Drag_Coefficient_W,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Drag_Coefficient_V,IOFlag))) THEN
         CALL ec_stat(cdragv(1:nlat,1:nlon),
     &        s1cdragv(1:nlat,1:nlon),  s2cdragv(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Drag_Coefficient_V,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Cloud_Cover,IOFlag))) THEN
         CALL ec_stat(tccd(1:nlat,1:nlon),
     &        s1tcc(1:nlat,1:nlon),  s2tcc (1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Cloud_Cover,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_T1,IOFlag))) THEN
         CALL ec_stat(dumt1(1:nlat,1:nlon,1),
     &        s1dt1(1:nlat,1:nlon,1),s2dt1 (1:nlat,1:nlon,1),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(dumt1(1:nlat,1:nlon,2),
     &        s1dt1(1:nlat,1:nlon,2),s2dt1 (1:nlat,1:nlon,2),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(dumt1(1:nlat,1:nlon,3),
     &        s1dt1(1:nlat,1:nlon,3),s2dt1 (1:nlat,1:nlon,3),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_T1,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_T2,IOFlag))) THEN
         CALL ec_stat(dumt2(1:nlat,1:nlon,1),
     &        s1dt2(1:nlat,1:nlon,1),s2dt2 (1:nlat,1:nlon,1),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(dumt2(1:nlat,1:nlon,2),
     &        s1dt2(1:nlat,1:nlon,2),s2dt2 (1:nlat,1:nlon,2),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(dumt2(1:nlat,1:nlon,3),
     &        s1dt2(1:nlat,1:nlon,3),s2dt2 (1:nlat,1:nlon,3),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_T2,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_U1,IOFlag))) THEN
         CALL ec_stat(dumu1(1:nlat,1:nlon,1),
     &        s1du1(1:nlat,1:nlon,1),s2du1 (1:nlat,1:nlon,1),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(dumu1(1:nlat,1:nlon,2),
     &        s1du1(1:nlat,1:nlon,2),s2du1 (1:nlat,1:nlon,2),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(dumu1(1:nlat,1:nlon,3),
     &        s1du1(1:nlat,1:nlon,3),s2du1 (1:nlat,1:nlon,3),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_U1,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_U2,IOFlag))) THEN
         CALL ec_stat(dumu2(1:nlat,1:nlon,1),
     &        s1du2(1:nlat,1:nlon,1),s2du2 (1:nlat,1:nlon,1),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         CALL ec_stat(dumu2(1:nlat,1:nlon,2),
     &        s1du2(1:nlat,1:nlon,2),s2du2 (1:nlat,1:nlon,2),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write)
         CALL ec_stat(dumu2(1:nlat,1:nlon,3),
     &        s1du2(1:nlat,1:nlon,3),s2du2 (1:nlat,1:nlon,3),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(User_Assigned_U2,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

#if ( ISOATM >= 1 )
      IF (output(newtotvar(Total_Precip18,IOFlag))) THEN
         CALL ec_stat(iso18torain(1:nlat,1:nlon),
     &      s1torain18(1:nlat,1:nlon),  s2torain18(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Precip18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Precip17,IOFlag))) THEN
         CALL ec_stat(iso17torain(1:nlat,1:nlon),
     &      s1torain17(1:nlat,1:nlon),s2torain17(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Precip17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Precipd,IOFlag))) THEN
         CALL ec_stat(isodtorain(1:nlat,1:nlon),
     &      s1toraind(1:nlat,1:nlon),s2toraind(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Precipd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snow18,IOFlag))) THEN
         CALL ec_stat(iso18tosnow(1:nlat,1:nlon),
     &      s1tosnow18(1:nlat,1:nlon),s2tosnow18(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Snow18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snow17,IOFlag))) THEN
         CALL ec_stat(iso17tosnow(1:nlat,1:nlon),
     &      s1tosnow17(1:nlat,1:nlon),s2tosnow17(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Snow17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snowd,IOFlag))) THEN
         CALL ec_stat(isodtosnow(1:nlat,1:nlon),
     &      s1tosnowd(1:nlat,1:nlon),s2tosnowd(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Total_Snowd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Spec_Hum_d18,IOFlag))) THEN
         CALL ec_stat(rmoisg(1:nlat,1:nlon,ieau18),
     &        s1rmoisg18(1:nlat,1:nlon),s2rmoisg18(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Spec_Hum_d18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Spec_Hum_d17,IOFlag))) THEN
         CALL ec_stat(rmoisg(1:nlat,1:nlon,ieau17),
     &        s1rmoisg17(1:nlat,1:nlon),s2rmoisg17(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Spec_Hum_d17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Spec_Hum_dd,IOFlag))) THEN
         CALL ec_stat(rmoisg(1:nlat,1:nlon,ieaud),
     &        s1rmoisgd(1:nlat,1:nlon),s2rmoisgd(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Spec_Hum_dd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

#endif
#if ( ISOATM >= 2 )
      IF (output(newtotvar(Land_Snow_D18,IOFlag))) THEN
         CALL ec_stat(adsnow(1:nlat,1:nlon,ieau18),
     &        s1dsnow18(1:nlat,1:nlon),s2dsnow18(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_D18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Land_Snow_D17,IOFlag))) THEN
         CALL ec_stat(adsnow(1:nlat,1:nlon,ieau17),
     &        s1dsnow17(1:nlat,1:nlon),s2dsnow17(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_D17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Land_Snow_Dd,IOFlag))) THEN
         CALL ec_stat(adsnow(1:nlat,1:nlon,ieaud),
     &        s1dsnowd(1:nlat,1:nlon),s2dsnowd(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_Dd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Bottom_Mois18,IOFlag))) THEN
         CALL ec_stat(abmoisg(1:nlat,1:nlon,ieau18),
     &        s1bmoisg18(1:nlat,1:nlon),s2bmoisg18(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_D18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Bottom_Mois17,IOFlag))) THEN
         CALL ec_stat(abmoisg(1:nlat,1:nlon,ieau17),
     &        s1bmoisg17(1:nlat,1:nlon),s2bmoisg17(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_D17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Bottom_Moisd,IOFlag))) THEN
         CALL ec_stat(abmoisg(1:nlat,1:nlon,ieaud),
     &        s1bmoisgd(1:nlat,1:nlon),s2bmoisgd(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write)
         IF (need_to_write) THEN
            CALL write(Land_Snow_Dd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
#endif

      IF (need_to_write) CALL close

      RETURN
      END SUBROUTINE ec_outputmyl

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_stat(xx,sumxx,sumxxsq,xxm,xxdev,compute)
!-----------------------------------------------------------------------
! *** this function sums xx and its square and computes the means and
! *** standard deviation at the end of a summation interval. The return
! *** value of this function indicates whether the output arrays xxm and
! *** xxdev contain valid data.
! *** all arrays are assumed to have the same shape.
!-----------------------------------------------------------------------

      USE comatm
      USE comdiag
      use comemic_mod, only:

      IMPLICIT NONE
      REAL*8, DIMENSION(:,:), INTENT(in)    :: xx
      REAL*8, DIMENSION(:,:), INTENT(inout) :: sumxx, sumxxsq
      REAL*8, DIMENSION(:,:), INTENT(out)   :: xxm, xxdev
      LOGICAL,                INTENT(in)    :: compute

      IF (itel /= 0) THEN
         IF (itel == 1) THEN
            sumxx  (:,:) = xx(:,:)
            sumxxsq(:,:) = xx(:,:)**2
         ELSE
            sumxx  (:,:) = sumxx  (:,:) + xx(:,:)
            sumxxsq(:,:) = sumxxsq(:,:) + xx(:,:)**2
         END IF

         IF (compute) THEN
            xxm  (:,:) = sumxx  (:,:)/minterv
            xxdev(:,:) = sumxxsq(:,:) - sumxx(:,:)**2/minterv
            WHERE (xxdev(:,:) < 0.0)
               xxdev(:,:) = 0.0
            ELSE WHERE
               xxdev (:,:) = Sqrt(xxdev(:,:)/(minterv - 1.0))
            END WHERE
         END IF
      END IF

      RETURN
      END SUBROUTINE ec_stat


!23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine ec_outamean(istep)
!-----------------------------------------------------------------------

! *** computes mean values of ocean basins: m/yr
!-----------------------------------------------------------------------


      USE comatm
      USE comdyn
      use comphys
      use comdiag
      use comemic_mod, only: iatm, nstpyear
      use comcoup_mod
      use comsurf_mod
      use comunit
      use newunit_mod, only: ocbasin_id

#if ( ISOATM >= 1 )
      USE iso_param_mod, only: ieau
#endif

      implicit none

      integer i,j,k,l,istep,ibas
      real*8  precmean(nbasa),evapmean(nbasa),runomean(nbasa)
      real*8 sum1(nlat,nbasa),sum2(nlat,nbasa),sum3(nlat,nbasa)
      real*8 ttest,unit

      unit=86400.

      if (istep.eq.1) then
        do ibas=1,nbasa
          precan(ibas)=0.0
          evapan(ibas)=0.0
          runoan(ibas)=0.0
        enddo
      endif

      do l=1,nbasa
         precmean(l)=0.
         evapmean(l)=0.
         runomean(l)=0.
      enddo

! *** calculate mean values of all basins
      do ibas=1,nbasa
        do i=1,nlat
          sum1(i,ibas)=0d0
          sum2(i,ibas)=0d0
          sum3(i,ibas)=0d0
        enddo
      enddo

      do j=1,nlon
       do i=1,nlat
         if (iocbasa(i,j).gt.0) then
c~ #if ( ISOATM >= 1 )
c~           sum1(i,iocbasa(i,j))=sum1(i,iocbasa(i,j))+fractoc(i,j)
c~      &         *(torain(i,j,ieau)+tosnow(i,j,ieau))
c~           sum2(i,iocbasa(i,j))=sum2(i,iocbasa(i,j))+fractoc(i,j)
c~      &         *evap(i,j,ieau)
c~ #else

c~ #endif
          sum1(i,iocbasa(i,j))=sum1(i,iocbasa(i,j))+fractoc(i,j)
     &         *(torain(i,j,iwater)+tosnow(i,j,iwater))
          sum2(i,iocbasa(i,j))=sum2(i,iocbasa(i,j))+fractoc(i,j)
     &         *evap(i,j,iwater)
          sum3(i,iocbasa(i,j))=sum3(i,iocbasa(i,j))+fractoc(i,j)
     &         *arunofo(i,j)
         endif

c~ #if ( ISOATM >= 1 )
c~          sum1(i,nbasa)=sum1(i,nbasa)+fractoc(i,j)
c~      &                  *(torain(i,j,ieau)+tosnow(i,j,ieau))
c~          sum2(i,nbasa)=sum2(i,nbasa)+fractoc(i,j)*evap(i,j,ieau)
c~ #else
c~ #endif
         sum1(i,nbasa)=sum1(i,nbasa)+fractoc(i,j)*(torain(i,j,iwater)
     &                +tosnow(i,j,iwater))
         sum2(i,nbasa)=sum2(i,nbasa)+fractoc(i,j)*evap(i,j,iwater)     
         sum3(i,nbasa)=sum3(i,nbasa)+fractoc(i,j)*arunofo(i,j)
        enddo
      enddo

      do i=1,nlat
       do ibas=1,nbasa
        precmean(ibas)=precmean(ibas)+sum1(i,ibas)*darea(i)
        evapmean(ibas)=evapmean(ibas)+sum2(i,ibas)*darea(i)
        runomean(ibas)=runomean(ibas)+sum3(i,ibas)*darea(i)
       enddo
      enddo
      do ibas=1,nbasa
        precan(ibas)=precan(ibas) + precmean(ibas)
        evapan(ibas)=evapan(ibas) + evapmean(ibas)
        runoan(ibas)=runoan(ibas) + runomean(ibas)
      enddo

      ttest=float(istep)/float(nstpyear)

      if (ttest.eq.int(ttest)) then
        do ibas=1,nbasa
          precan(ibas)=precan(ibas)*unit/(arocbasa(ibas)*float(iatm))
          evapan(ibas)=evapan(ibas)*unit/(arocbasa(ibas)*float(iatm))
          runoan(ibas)=runoan(ibas)*unit/(arocbasa(ibas)*float(iatm))
        enddo
        write(ocbasin_id,*) 'precip'
        write(ocbasin_id,*)  precan
        write(ocbasin_id,*) 'evap'
        write(ocbasin_id,*)  evapan
        write(ocbasin_id,*) 'runoff'
        write(ocbasin_id,*)  runoan
        write(ocbasin_id,*)
        call flush(ocbasin_id)
        do ibas=1,nbasa
          precan(ibas)=0.0
          evapan(ibas)=0.0
          runoan(ibas)=0.0
        enddo
      endif
      end subroutine ec_outamean

!23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine ec_outiocht(istep)
!-----------------------------------------------------------------------

! *** computes annual mean implied ocean heat transport: PW (10^15 W)
!-----------------------------------------------------------------------


      USE comatm
      USE comdyn
      use comphys
      use comdiag
      use comemic_mod, only: day, iatm, nstpyear
      use comcoup_mod
      use comsurf_mod
      use comunit
      use newunit_mod, only: ocheattr_id

      implicit none

      integer i,j,k,l,istep,ibas
      real*8  hfmean(nlat,nbasa)
      real*8  hatl(nlat),hpac(nlat),hind(nlat),hall(nlat)
      real*8  ttest,unit,dum

      unit=1e-15
      if (istep.eq.1) then
        do ibas=1,nbasa
          do i=1,nlat
            hfmeanan(i,ibas)=0d0
          enddo
        enddo
      endif
      do ibas=1,nbasa
        do i=1,nlat
          hfmean(i,ibas)=0d0
        enddo
      enddo

! *** calculate zonal mean net surface heatflux of all basins

      do j=1,nlon
        do i=1,nlat
          dum=fractoc(i,j)*(heswsn(i,j,noc)+dlradsn(i,j,noc)-
     &        ulradsn(i,j,noc)-efluxn(i,j,noc)-hfluxn(i,j,noc))
          if (iocbasa(i,j).gt.0) then
            hfmean(i,iocbasa(i,j))=hfmean(i,iocbasa(i,j))+dum
          endif
          hfmean(i,nbasa)=hfmean(i,nbasa)+dum
        enddo
      enddo

      do ibas=1,nbasa
        do i=1,nlat
          hfmeanan(i,ibas)=hfmeanan(i,ibas) + hfmean(i,ibas)*darea(i)
        enddo
      enddo

      if(mod(nint(day*real(iatm))+1,nstpyear).eq.0) then
!if (mod(istep,nstpyear).eq.0) then
        do ibas=1,nbasa
          do i=1,nlat
            hfmeanan(i,ibas)=hfmeanan(i,ibas)/real(nstpyear)
          enddo
        enddo
        do i=1,nlat
          hatl(i)=0d0
          hpac(i)=0d0
          hind(i)=0d0
          hall(i)=0d0
        enddo
        hatl(nlat)=-hfmeanan(nlat,1)
        hpac(nlat)=-hfmeanan(nlat,2)
        hind(nlat)=-hfmeanan(nlat,4)
        hall(nlat)=-hfmeanan(nlat,nbasa)
        do i=nlat-1,1,-1
          hatl(i)=hatl(i+1)-hfmeanan(i,1)-hfmeanan(i,5)
          hpac(i)=hpac(i+1)-hfmeanan(i,2)-hfmeanan(i,6)
          hind(i)=hind(i+1)-hfmeanan(i,4)
          hall(i)=hall(i+1)-hfmeanan(i,nbasa)
        enddo
        write(ocheattr_id) (real(hatl(i)*unit),i=1,nlat)
        write(ocheattr_id) (real(hpac(i)*unit),i=1,nlat)
        write(ocheattr_id) (real(hind(i)*unit),i=1,nlat)
        write(ocheattr_id) (real(hall(i)*unit),i=1,nlat)
        call flush(ocheattr_id)
        do ibas=1,nbasa
          do i=1,nlat
            hfmeanan(i,ibas)=0d0
          enddo
        enddo
      endif
      return
      end subroutine ec_outiocht

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_outputyrl
!-----------------------------------------------------------------------
! *** this routine calls another routine which computes the yearly mean
! *** written by camiel severijns, march 2001.
!-----------------------------------------------------------------------
      use Atmosphere_Output

      USE comatm
      USE comdyn
      use comphys
      use comdiag
      use comemic_mod, only: iday, imonth
      use comcoup_mod
      use comoutlocal_mod
      use comsurf_mod

! --- BdB 05-2019: added variables for rewriting files
      use comemic_mod, only: iyear, nwrskip, nyears
     &                     , new_year_atm, current_int_atm


#if (ISOATM >= 1 )
       USE iso_param_mod, ONLY : ieau,ieau17,ieau18,ieaud
#endif

      IMPLICIT NONE

      REAL*8, DIMENSION(1:nlat,1:nlon), SAVE ::
     &     tsurfs, tsurfs2, tstrats, tstrats2,
     &     temp2gs, temp2gs2, temp4gs, temp4gs2,
     &     tempsgs, tempsgs2, tstrat2s, tstrat2s2,
     &     t2ms, t2ms2, u200s, u200s2, u500s, u500s2,
     &     u800s, u800s2, v200s, v200s2, v500s, v500s2,
     &     v800s, v800s2, pgrounds, pgrounds2, omeg1s,
     &     omeg1s2, omeg2s, omeg2s2, omeg3s, omeg3s2,
     &     winstu1s, winstu1s2, winstv1s, winstv1s2,
     &     uv10s, uv10s2, udivg1s, udivg1s2, udivg2s,
     &     udivg2s2, udivg3s, udivg3s2, vdivg1s, vdivg1s2,
     &     vdivg2s, vdivg2s2, vdivg3s, vdivg3s2, grpsi1s,
     &     grpsi1s2, grpsi2s, grpsi2s2, grpsi3s, grpsi3s2,
     &     chi1s, chi1s2, chi2s, chi2s2, chi3s, chi3s2,
     &     qgpv1s, qgpv1s2, qgpv2s, qgpv2s2, qgpv3s,
     &     qgpv3s2, gh1s, gh1s2, gh2s, gh2s2, gh3s, gh3s2,
     &     vhforg1s, vhforg1s2, vhforg2s, vhforg2s2
      REAL*8, DIMENSION(1:nlat,1:nlon), SAVE ::
     &     vforg1s, vforg1s2, vforg2s, vforg2s2, vforg3s,
     &     vforg3s2, dyrains, dyrains2, corains, corains2,
     &     torains, torains2, snows, snows2, evaps, evaps2,
     &     eminps, eminps2, hfluxs, hfluxs2, efluxs,
     &     efluxs2, heswss, heswss2, hesw_s, hesw_s2,
     &     swrss, swrss2, !mohr
     &     lwrss, lwrss2, !mohr
     &     albess, albess2, albeps, albeps2, nlradss,
     &     nlradss2, ulrad1s, ulrad1s2, bmoisgs, bmoisgs2,
     &     dsnows, dsnows2, hics, hics2, runofls, runofls2,
     &     runofos, runofos2, rmoisgw3s, rmoisgw3s2,
     &     relhums, relhums2, cdragws, cdragws2, cdragvs,
     &     cdragvs2, tccs, tccs2, dt11s, dt11s2, dt12s,
     &     dt12s2, dt13s, dt13s2, dt21s, dt21s2, dt22s,
     &     dt22s2, dt23s, dt23s2, du11s, du11s2, du12s,
     &     du12s2, du13s, du13s2, du21s, du21s2, du22s,
     &     du22s2, du23s, du23s2
#if ( ISOATM >= 1 )
     &    ,torains18, torains218
     &    ,torains17, torains217
     &    ,torainsd, torains2d
     &    ,tosnows18, tosnows218
     &    ,tosnows17, tosnows217
     &    ,tosnowsd, tosnows2d
     &    ,rmoisg18s, rmoisg18s2
     &    ,rmoisg17s, rmoisg17s2
     &    ,rmoisgds,  rmoisgds2
#endif
#if ( ISOATM >= 1 )
     &    ,dsnows18, dsnows218
     &    ,dsnows17, dsnows217
     &    ,dsnowsd, dsnows2d
     &    ,bmoisgs18, bmoisgs218
     &    ,bmoisgs17, bmoisgs217
     &    ,bmoisgsd, bmoisgs2d
#endif

      INTERFACE
         SUBROUTINE ec_yrlstat(xx,sumxx,sumxxsq,xxm,xxdev,compute,
     &     samples)
         IMPLICIT NONE
         REAL*8, DIMENSION(:,:), INTENT(in)    :: xx
         REAL*8, DIMENSION(:,:), INTENT(inout) :: sumxx, sumxxsq
         REAL*8, DIMENSION(:,:), INTENT(out)   :: xxm, xxdev
         LOGICAL,                INTENT(in)    :: compute
         INTEGER,                INTENT(in)    :: samples
         END SUBROUTINE ec_yrlstat
      END INTERFACE

      INTEGER, PARAMETER :: IOFlag = 4

      LOGICAL :: need_to_write
      REAL*8, DIMENSION(1:nlat,1:nlon,1:4) :: mean, stddev
      INTEGER, SAVE :: samples = 0
      INTEGER :: realmonth

      realmonth=12
!mab: wrong!!!!!!!!!!!!!!!!
!       if ( irunlabeld/=360 ) realmonth=(irunlabeld/30)
!      if (iday.eq.1) print*,'realmonth: ',realmonth,irunlabeld
!     &,imonth,iyear,irunlabel

      samples = samples + 1 ! Increment sampel counter

      need_to_write = ( imonth == realmonth .AND. iday == 30 )

! --- BdB 05-2019: calculate start time
      IF (need_to_write) THEN
        new_year_atm = new_year_atm + 1
        CALL open(Yearly_Means)
      END IF

      IF (output(newtotvar(Surface_Temperature,IOFlag))) THEN
         CALL ec_yrlstat(tsurf1(1:nlat,1:nlon),
     &        tsurfs(1:nlat,1:nlon),   tsurfs2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1), stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Temperature,IOFlag))) THEN
         CALL ec_yrlstat(temp0g1(1:nlat,1:nlon),
     &        tstrats(1:nlat,1:nlon),  tstrats2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(temp2g1(1:nlat,1:nlon),
     &        temp2gs(1:nlat,1:nlon),  temp2gs2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(temp4g1(1:nlat,1:nlon),
     &        temp4gs(1:nlat,1:nlon),  temp4gs2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3),
     &        need_to_write,samples)
         CALL ec_yrlstat(tempsg1(1:nlat,1:nlon),
     &        tempsgs(1:nlat,1:nlon),  tempsgs2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,4),stddev  (1:nlat,1:nlon,4),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Temperature,mean(1:nlat,1:nlon,1:4))
         END IF
      END IF

      IF (output(newtotvar(Stratospheric_Temperature,IOFlag))) THEN
         CALL ec_yrlstat(temp0g1(1:nlat,1:nlon),
     &        tstrat2s(1:nlat,1:nlon),  tstrat2s2(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
           CALL write(Stratospheric_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Two_Meter_Temperature,IOFlag))) THEN
         CALL ec_yrlstat(tempsg1(1:nlat,1:nlon),
     &        t2ms(1:nlat,1:nlon),  t2ms2 (1:nlat,1:nlon),
     &        mean(1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Two_Meter_Temperature,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Wind_U,IOFlag))) THEN
         CALL ec_yrlstat(u200(1:nlat,1:nlon),
     &        u200s(1:nlat,1:nlon),  u200s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(u500(1:nlat,1:nlon),
     &        u500s(1:nlat,1:nlon),  u500s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(u800(1:nlat,1:nlon),
     &        u800s(1:nlat,1:nlon),  u800s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Wind_U,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Wind_V,IOFlag))) THEN
         CALL ec_yrlstat(v200(1:nlat,1:nlon),
     &        v200s(1:nlat,1:nlon),  v200s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(v500(1:nlat,1:nlon),
     &        v500s(1:nlat,1:nlon),  v500s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(v800(1:nlat,1:nlon),
     &        v800s(1:nlat,1:nlon),  v800s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Wind_V,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Surface_Pressure,IOFlag))) THEN
         CALL ec_yrlstat(pground(1:nlat,1:nlon),
     &        pgrounds(1:nlat,1:nlon),  pgrounds2(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Pressure,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Vertical_Pressure_Wind,IOFlag))) THEN
         CALL ec_yrlstat(omegg(1:nlat,1:nlon,1),
     &        omeg1s(1:nlat,1:nlon),  omeg1s2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(omegg(1:nlat,1:nlon,2),
     &        omeg2s(1:nlat,1:nlon),  omeg2s2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,2),stddev (1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(omegg(1:nlat,1:nlon,3),
     &        omeg3s(1:nlat,1:nlon),  omeg3s2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,3),stddev (1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Vertical_Pressure_Wind,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(U_Stress,IOFlag))) THEN
         CALL ec_yrlstat(winstu1(1:nlat,1:nlon),
     &        winstu1s(1:nlat,1:nlon),  winstu1s2(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(U_Stress,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(V_Stress,IOFlag))) THEN
         CALL ec_yrlstat(winstv1(1:nlat,1:nlon),
     &        winstv1s(1:nlat,1:nlon),  winstv1s2(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(V_Stress,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Wind_at_10_Meter,IOFlag))) THEN
         CALL ec_yrlstat(uv10(1:nlat,1:nlon),
     &        uv10s(1:nlat,1:nlon),  uv10s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Wind_at_10_Meter,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_U,IOFlag))) THEN
         CALL ec_yrlstat(udivg(1:nlat,1:nlon,1),
     &        udivg1s(1:nlat,1:nlon),  udivg1s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(udivg(1:nlat,1:nlon,2),
     &        udivg2s(1:nlat,1:nlon),  udivg2s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(udivg(1:nlat,1:nlon,3),
     &        udivg3s(1:nlat,1:nlon),  udivg3s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Ageostrophic_Wind_U,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Ageostrophic_Wind_V,IOFlag))) THEN
         CALL ec_yrlstat(vdivg(1:nlat,1:nlon,1),
     &        vdivg1s(1:nlat,1:nlon),  vdivg1s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(vdivg(1:nlat,1:nlon,2),
     &        vdivg2s(1:nlat,1:nlon),  vdivg2s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(vdivg(1:nlat,1:nlon,3),
     &        vdivg3s(1:nlat,1:nlon),  vdivg3s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Ageostrophic_Wind_V,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Stream_Function,IOFlag))) THEN
         CALL ec_yrlstat(psig(1:nlat,1:nlon,1),
     &        grpsi1s(1:nlat,1:nlon),  grpsi1s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(psig(1:nlat,1:nlon,2),
     &        grpsi2s(1:nlat,1:nlon),  grpsi2s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(psig(1:nlat,1:nlon,3),
     &        grpsi3s(1:nlat,1:nlon),  grpsi3s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Stream_Function,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Velocity_Potential,IOFlag))) THEN
         CALL ec_yrlstat(chig(1:nlat,1:nlon,1),
     &        chi1s(1:nlat,1:nlon),  chi1s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(chig(1:nlat,1:nlon,1),
     &        chi2s(1:nlat,1:nlon),  chi2s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(chig(1:nlat,1:nlon,1),
     &        chi3s(1:nlat,1:nlon),  chi3s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Velocity_Potential,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(QG_Potential_Vorticity,IOFlag))) THEN
         CALL ec_yrlstat(qgpv(1:nlat,1:nlon,1),
     &        qgpv1s(1:nlat,1:nlon),  qgpv1s2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(qgpv(1:nlat,1:nlon,2),
     &        qgpv2s(1:nlat,1:nlon),  qgpv2s2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,2),stddev (1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(qgpv(1:nlat,1:nlon,3),
     &        qgpv3s(1:nlat,1:nlon),  qgpv3s2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,3),stddev (1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(QG_Potential_Vorticity,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Geopotential_Height,IOFlag))) THEN
         CALL ec_yrlstat(geopg(1:nlat,1:nlon,1),
     &        gh1s(1:nlat,1:nlon),  gh1s2 (1:nlat,1:nlon),
     &        mean(1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(geopg(1:nlat,1:nlon,2),
     &        gh2s(1:nlat,1:nlon),  gh2s2 (1:nlat,1:nlon),
     &        mean(1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(geopg(1:nlat,1:nlon,3),
     &        gh3s(1:nlat,1:nlon),  gh3s2 (1:nlat,1:nlon),
     &        mean(1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Geopotential_Height,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Heating_Force,IOFlag))) THEN
         CALL ec_yrlstat(vhforg1(1:nlat,1:nlon),
     &        vhforg1s(1:nlat,1:nlon),  vhforg1s2(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,1),stddev   (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(vhforg2(1:nlat,1:nlon),
     &        vhforg2s(1:nlat,1:nlon),  vhforg2s2(1:nlat,1:nlon),
     &        mean    (1:nlat,1:nlon,2),stddev   (1:nlat,1:nlon,2),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Heating_Force,mean(1:nlat,1:nlon,1:2))
         END IF
      END IF

      IF (output(newtotvar(Potential_Vorticity_Forcing,IOFlag))) THEN
         CALL ec_yrlstat(vforg1(1:nlat,1:nlon),
     &        vforg1s(1:nlat,1:nlon),  vforg1s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(vforg2(1:nlat,1:nlon),
     &        vforg2s(1:nlat,1:nlon),  vforg2s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,2),stddev  (1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(vforg3(1:nlat,1:nlon),
     &        vforg3s(1:nlat,1:nlon),  vforg3s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,3),stddev  (1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
         CALL write(Potential_Vorticity_Forcing,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(Large_Scale_Precipitation,IOFlag))) THEN
         CALL ec_yrlstat(dyrain1(1:nlat,1:nlon),
     &        dyrains(1:nlat,1:nlon),  dyrains2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Large_Scale_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Convective_Precipitation,IOFlag))) THEN
         CALL ec_yrlstat(corain1(1:nlat,1:nlon),
     &        corains(1:nlat,1:nlon),  corains2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Convective_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Precipitation,IOFlag))) THEN
         CALL ec_yrlstat(torain1(1:nlat,1:nlon),
     &        torains(1:nlat,1:nlon),   torains2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1), stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Total_Precipitation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snow_Fall,IOFlag))) THEN
         CALL ec_yrlstat(snow1(1:nlat,1:nlon),
     &        snows(1:nlat,1:nlon),  snows2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Total_Snow_Fall,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Evaporation,IOFlag))) THEN
         CALL ec_yrlstat(evap1(1:nlat,1:nlon),
     &        evaps(1:nlat,1:nlon),  evaps2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Evaporation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Evap_Minus_Precip,IOFlag))) THEN
         CALL ec_yrlstat(eminp1(1:nlat,1:nlon),
     &        eminps(1:nlat,1:nlon),  eminps2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Evap_Minus_Precip,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Sensible_Heat_Flux,IOFlag))) THEN
         CALL ec_yrlstat(hflux(1:nlat,1:nlon),
     &        hfluxs(1:nlat,1:nlon),  hfluxs2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Sensible_Heat_Flux,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Latent_Heat_Flux,IOFlag))) THEN
         CALL ec_yrlstat(eflux(1:nlat,1:nlon),
     &        efluxs(1:nlat,1:nlon),  efluxs2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Latent_Heat_Flux,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Solar_Radiation,IOFlag))) THEN
         CALL ec_yrlstat(hesws(1:nlat,1:nlon),
     &        heswss(1:nlat,1:nlon),  heswss2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

!mohr
      IF (output(newtotvar(Surface_Shortwave_Radiation,IOFlag))) THEN
         CALL ec_yrlstat(swrs1(1:nlat,1:nlon),
     &        swrss(1:nlat,1:nlon),  swrss2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Shortwave_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Longwave_Radiation,IOFlag))) THEN
         CALL ec_yrlstat(lwrs1(1:nlat,1:nlon),
     &        lwrss(1:nlat,1:nlon),  lwrss2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Longwave_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

!!mohr

      IF (output(newtotvar(Top_Solar_Radiation,IOFlag))) THEN
         CALL ec_yrlstat(hesw(1:nlat,1:nlon),
     &        hesw_s(1:nlat,1:nlon),  hesw_s2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Top_Solar_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Albedo,IOFlag))) THEN
         CALL ec_yrlstat(alb2es(1:nlat,1:nlon),
     &        albess(1:nlat,1:nlon),  albess2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Albedo,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Planetary_Albedo,IOFlag))) THEN
         CALL ec_yrlstat(albep(1:nlat,1:nlon),
     &        albeps(1:nlat,1:nlon),  albeps2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Planetary_Albedo,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Surface_Thermal_Radiation,IOFlag))) THEN
         CALL ec_yrlstat(nlrads(1:nlat,1:nlon),
     &        nlradss(1:nlat,1:nlon),  nlradss2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Surface_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Top_Thermal_Radiation,IOFlag))) THEN
         CALL ec_yrlstat(ulrad0(1:nlat,1:nlon),
     &        ulrad1s(1:nlat,1:nlon),  ulrad1s2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Top_Thermal_Radiation,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Bottom_Moisture,IOFlag))) THEN
c~ #if ( ISOATM >= 2 )
c~          CALL ec_yrlstat(abmoisg(1:nlat,1:nlon,ieau),
c~      &        bmoisgs(1:nlat,1:nlon),  bmoisgs2(1:nlat,1:nlon),
c~      &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
c~      &        need_to_write,samples)
c~ #else
         CALL ec_yrlstat(abmoisg(1:nlat,1:nlon,iwater),
     &        bmoisgs(1:nlat,1:nlon),  bmoisgs2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
c~ #endif
         IF (need_to_write) THEN
            CALL write(Bottom_Moisture,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Land_Snow_Depth,IOFlag))) THEN
c~ #if ( ISOATM >= 2 )
c~          CALL ec_yrlstat(adsnow(1:nlat,1:nlon,ieau),
c~      &        dsnows(1:nlat,1:nlon),  dsnows2(1:nlat,1:nlon),
c~      &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
c~      &        need_to_write,samples)
c~ #else
         CALL ec_yrlstat(adsnow(1:nlat,1:nlon,iwater),
     &        dsnows(1:nlat,1:nlon),  dsnows2(1:nlat,1:nlon),
     &        mean  (1:nlat,1:nlon,1),stddev (1:nlat,1:nlon,1),
     &        need_to_write,samples)
c~ #endif
         IF (need_to_write) THEN
            CALL write(Land_Snow_Depth,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Sea_Ice_Thickness,IOFlag))) THEN
         CALL ec_yrlstat(ahic(1:nlat,1:nlon),
     &        hics(1:nlat,1:nlon),  hics2 (1:nlat,1:nlon),
     &        mean(1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Sea_Ice_Thickness,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Land_Surface_Runoff,IOFlag))) THEN
         CALL ec_yrlstat(runofl1(1:nlat,1:nlon),
     &        runofls(1:nlat,1:nlon),  runofls2(1:nlat,1:nlon),
     &        mean(1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Land_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Ocean_Surface_Runoff,IOFlag))) THEN
         CALL ec_yrlstat(runofo1(1:nlat,1:nlon),
     &        runofos(1:nlat,1:nlon),  runofos2(1:nlat,1:nlon),
     &        mean(1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Ocean_Surface_Runoff,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Specific_Humidity,IOFlag))) THEN
c~ #if (ISOATM >= 1 )
c~          CALL ec_yrlstat(rmoisg(1:nlat,1:nlon,ieau),
c~      &        rmoisgw3s(1:nlat,1:nlon),  rmoisgw3s2(1:nlat,1:nlon),
c~      &        mean     (1:nlat,1:nlon,1),stddev    (1:nlat,1:nlon,1),
c~      &        need_to_write,samples)
c~ #else
         CALL ec_yrlstat(rmoisg(1:nlat,1:nlon,iwater),
     &        rmoisgw3s(1:nlat,1:nlon),  rmoisgw3s2(1:nlat,1:nlon),
     &        mean(1:nlat,1:nlon,1),stddev    (1:nlat,1:nlon,1),
     &        need_to_write,samples)
c~ #endif
         IF (need_to_write) THEN
            CALL write(Specific_Humidity,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Relative_Humidity,IOFlag))) THEN
         CALL ec_yrlstat(relhum(1:nlat,1:nlon),
     &        relhums(1:nlat,1:nlon),  relhums2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Relative_Humidity,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Drag_Coefficient_W,IOFlag))) THEN
         CALL ec_yrlstat(cdragw(1:nlat,1:nlon),
     &        cdragws(1:nlat,1:nlon),  cdragws2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Drag_Coefficient_W,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Drag_Coefficient_V,IOFlag))) THEN
         CALL ec_yrlstat(cdragv(1:nlat,1:nlon),
     &        cdragvs(1:nlat,1:nlon),  cdragvs2(1:nlat,1:nlon),
     &        mean   (1:nlat,1:nlon,1),stddev  (1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Drag_Coefficient_V,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Cloud_Cover,IOFlag))) THEN
         CALL ec_yrlstat(tccd(1:nlat,1:nlon),
     &        tccs(1:nlat,1:nlon),  tccs2 (1:nlat,1:nlon),
     &        mean(1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Total_Cloud_Cover,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_T1,IOFlag))) THEN
         CALL ec_yrlstat(dumt1(1:nlat,1:nlon,1),
     &        dt11s(1:nlat,1:nlon),  dt11s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(dumt1(1:nlat,1:nlon,2),
     &        dt12s(1:nlat,1:nlon),  dt12s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(dumt1(1:nlat,1:nlon,3),
     &        dt13s(1:nlat,1:nlon),  dt13s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(User_Assigned_T1,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_T2,IOFlag))) THEN
         CALL ec_yrlstat(dumt2(1:nlat,1:nlon,1),
     &        dt21s(1:nlat,1:nlon),  dt21s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(dumt2(1:nlat,1:nlon,2),
     &        dt22s(1:nlat,1:nlon),  dt22s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(dumt2(1:nlat,1:nlon,3),
     &        dt23s(1:nlat,1:nlon),  dt23s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(User_Assigned_T2,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_U1,IOFlag))) THEN
         CALL ec_yrlstat(dumu1(1:nlat,1:nlon,1),
     &        du11s(1:nlat,1:nlon),  du11s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(dumu1(1:nlat,1:nlon,2),
     &        du12s(1:nlat,1:nlon),  du12s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(dumu1(1:nlat,1:nlon,3),
     &        du13s(1:nlat,1:nlon),  du13s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(User_Assigned_U1,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

      IF (output(newtotvar(User_Assigned_U2,IOFlag))) THEN
         CALL ec_yrlstat(dumu2(1:nlat,1:nlon,1),
     &        du21s(1:nlat,1:nlon),  du21s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         CALL ec_yrlstat(dumu2(1:nlat,1:nlon,2),
     &        du22s(1:nlat,1:nlon),  du22s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,2),stddev(1:nlat,1:nlon,2),
     &        need_to_write,samples)
         CALL ec_yrlstat(dumu2(1:nlat,1:nlon,3),
     &        du23s(1:nlat,1:nlon),  du23s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,3),stddev(1:nlat,1:nlon,3),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(User_Assigned_U2,mean(1:nlat,1:nlon,1:3))
         END IF
      END IF

#if ( ISOATM >= 1 )
      IF (output(newtotvar(Total_Precip18,IOFlag))) THEN
         CALL ec_yrlstat(iso18torain(1:nlat,1:nlon),
     &      torains18(1:nlat,1:nlon),torains218(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Total_Precip18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Precip17,IOFlag))) THEN
         CALL ec_yrlstat(iso17torain(1:nlat,1:nlon),
     &      torains17(1:nlat,1:nlon),torains217(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Total_Precip17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Precipd,IOFlag))) THEN
         CALL ec_yrlstat(isodtorain(1:nlat,1:nlon),
     &      torainsd(1:nlat,1:nlon),torains2d(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Total_Precipd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snow18,IOFlag))) THEN
         CALL ec_yrlstat(iso18tosnow(1:nlat,1:nlon),
     &      tosnows18(1:nlat,1:nlon),tosnows218(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Total_Snow18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snow17,IOFlag))) THEN
         CALL ec_yrlstat(iso17tosnow(1:nlat,1:nlon),
     &      tosnows17(1:nlat,1:nlon),tosnows217(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Total_Snow17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Total_Snowd,IOFlag))) THEN
         CALL ec_yrlstat(isodtosnow(1:nlat,1:nlon),
     &      tosnowsd(1:nlat,1:nlon),tosnows2d(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Total_Snowd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Spec_Hum_d18,IOFlag))) THEN
         CALL ec_yrlstat(rmoisg(1:nlat,1:nlon,ieau18),
     &        rmoisg18s(1:nlat,1:nlon),rmoisg18s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Spec_Hum_d18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Spec_Hum_d17,IOFlag))) THEN
         CALL ec_yrlstat(rmoisg(1:nlat,1:nlon,ieau17),
     &        rmoisg17s(1:nlat,1:nlon),rmoisg17s2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Spec_Hum_d17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

      IF (output(newtotvar(Spec_Hum_dd,IOFlag))) THEN
         CALL ec_yrlstat(rmoisg(1:nlat,1:nlon,ieaud),
     &        rmoisgds(1:nlat,1:nlon),rmoisgds2(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Spec_Hum_dd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF

#endif
#if ( ISOATM >= 2 )
      IF (output(newtotvar(Land_Snow_D18,IOFlag))) THEN
         CALL ec_yrlstat(adsnow(1:nlat,1:nlon,ieau18),
     &        dsnows18(1:nlat,1:nlon),dsnows218(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Land_Snow_D18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Land_Snow_D17,IOFlag))) THEN
         CALL ec_yrlstat(adsnow(1:nlat,1:nlon,ieau17),
     &        dsnows17(1:nlat,1:nlon),dsnows217(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Land_Snow_D17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Land_Snow_Dd,IOFlag))) THEN
         CALL ec_yrlstat(adsnow(1:nlat,1:nlon,ieaud),
     &        dsnowsd(1:nlat,1:nlon),dsnows2d(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Land_Snow_Dd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Bottom_Mois18,IOFlag))) THEN
         CALL ec_yrlstat(abmoisg(1:nlat,1:nlon,ieau18),
     &        bmoisgs18(1:nlat,1:nlon),bmoisgs218(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Land_Snow_D18,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Bottom_Mois17,IOFlag))) THEN
         CALL ec_yrlstat(abmoisg(1:nlat,1:nlon,ieau17),
     &        bmoisgs17(1:nlat,1:nlon),bmoisgs217(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Land_Snow_D17,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
      IF (output(newtotvar(Bottom_Moisd,IOFlag))) THEN
         CALL ec_yrlstat(abmoisg(1:nlat,1:nlon,ieaud),
     &        bmoisgsd(1:nlat,1:nlon),bmoisgs2d(1:nlat,1:nlon),
     &        mean (1:nlat,1:nlon,1),   stddev(1:nlat,1:nlon,1),
     &        need_to_write,samples)
         IF (need_to_write) THEN
            CALL write(Land_Snow_Dd,mean(1:nlat,1:nlon,1:1))
         END IF
      END IF
#endif
      IF (need_to_write) THEN
         CALL close
         samples = 0 ! Restart counting samples

! ---    BdB 05-2019: reset year counter
         if (mod(iyear,nwrskip).eq.0 .or. iyear.eq.nyears) then
          new_year_atm = 0
          current_int_atm = current_int_atm + nwrskip
         end if
      END IF

      RETURN
      END SUBROUTINE ec_outputyrl

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_yrlstat(x,sumx1,sumy1,xmean,xstd,compute,samples)
! *** --------------------------------------------------------------------
! *** This routine computes yearly mean and standard  deviation around it.
! *** all arrays are assumed to have the same shape.
! *** ------------------------------------------------------------------------


      USE comatm
      USE comdiag
      use comemic_mod, only: iday, imonth
      use comcoup_mod

      IMPLICIT NONE
      REAL*8, DIMENSION(:,:), INTENT(in)    :: x
      REAL*8, DIMENSION(:,:), INTENT(inout) :: sumx1, sumy1
      REAL*8, DIMENSION(:,:), INTENT(out)   :: xmean, xstd
      LOGICAL,                INTENT(in)    :: compute
      INTEGER,                INTENT(in)    :: samples


      IF (imonth == 1 .AND. iday == 1) THEN
         sumx1(:,:) = 0.0
         sumy1(:,:) = 0.0
      END IF

      sumx1(:,:) = sumx1(:,:) + x(:,:)
      sumy1(:,:) = sumy1(:,:) + x(:,:)**2

      IF (compute) THEN
         xmean(:,:) = sumx1(:,:)/samples
         xstd (:,:) = sumy1(:,:) - sumx1(:,:)**2/samples
         WHERE (xstd(:,:) <= 0.0)
            xstd(:,:) = 0.0
         ELSE WHERE
            xstd(:,:) = Sqrt(xstd(:,:)/(samples - 1.0))
         END WHERE
      END IF

      RETURN
      END SUBROUTINE ec_yrlstat
