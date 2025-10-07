!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of [insert sub-component name here, in following Foobar]
!!      Foobar is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
!!      as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!!
!!      Foobar is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!!      of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!!
!!      You should have received a copy of the GNU General Public License along with Foobar.
!!      If not, see <http://www.gnu.org/licenses/>.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"
!dmr -- Adding the choice of components through the pre-processing options

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: atmrad_mod
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module atmrad_mod is collecting the radiation related routines initially located in atmphy0.f
!
!>     @date Creation date: Thu Nov 21 08:15:37 CET 2019
!>     @date Last modification: Fri Nov 22 08:35:30 CET 2019
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module atmrad_mod

       implicit none
       private

       public :: ec_solar, celest

      ! NOTE_avoid_public_variables_if_possible

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ec_solar
!
!>     @brief Calculate the incoming solar radiation
!
!      DESCRIPTION:
!!        Calculates incoming solar radiation as a FUNCTION of latitude for each day of the year, given the orbital parameters
!!          One year has 360 days.  Reference: A. Berger, JAS, 35, 2362-2367,1978
!>
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ec_solar(istep) result(returnValue)

        use comatm,               only: nlat, tanfi, sinfi

#if ( HOURLY_RAD > 0 )
        use comatm,               only: cosfi
        use comemic_mod,          only: iatm
        use comphys,              only: solflux_hourly, n_steps_day
#endif
        use comphys,              only: ecc, ecc2, eccf, facttsi, iscencel, isceno3, iscentsi, iscenvol, iso3strt, istsistrt,   &
                    isvolstrt, obl, oblf, omweb, omwebf, perh, so, solarc, solardref, solarm, solarcl, solarvol, solarcl,       &
                    solartsi, ghg, kosz, solarf, dayfr, q0
        use comemic_mod,          only: iday, imonth, iyear, nstpyear
        use comrunlabel_mod,      only: irunlabelf
        use comunit,              only: iuo

        use global_constants_mod, only: dblp=>dp, ip, deg2rad=>deg_to_rad, pi=>pi_dp, sp
        use newunit_mod, only: info_id

!! #define DINSOL

#ifdef DINSOL

        use dinsol_mod, only: dinsol_init, get_orbit_param, orb_param

        use comphys,              only: irunlabeloffset, iscenyear
        use comemic_mod,          only: iyear
        use comrunlabel_mod,      only: irunlabelf
    
#if ( F_PALAEO >= 1 )
       use palaeo_timer_mod,     only: celest_year
#endif


#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in] istep
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer, intent(in) :: istep

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        logical                          :: returnValue, success
        integer(kind=ip)                 :: i,j,l,NVE,indxvol,indxtsi
        real(kind=dblp)                  :: beta,alam,alam0,ala,ala0, fac1,fac2,fac3,ro,roref,deltal, sindl, cosdl, tandl, rkosz &
                         , rkosha1, ha1, solard, tsi, ksw
#if ( HOURLY_RAD > 0 )
! dmr hourly
        integer(kind=ip)                 :: hi
        real(kind=dblp)                  :: heure, h, cz, solflux
#endif
        real(kind=dblp), dimension(nlat) :: solarcf
        real(kind=dblp), dimension(2)    :: alpho3sw

#ifdef DINSOL
        logical :: valid_request
        type(orb_param) :: OB
        real(kind=sp) :: annee
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! Present-day orbital parameters: eccentricity ecc, obliquity obl and
! angle om between Vernal Equinox and Perihelion (angles all given
! in degrees and converted to radians). Solarc is the solar constant.
! NVE is day of the Vernal Equinox, set at 21 MARCH
! Implementatie van Nanne (Weber)
!
! [NOTA/DELETE] solarc is a candidate for deletion as a global variable
!               could easily be moved to a local variable of the present
!               routine. Solar Constant value is transmitted through solarm

      if ( iscencel.eq.1) then
! --- Transient Celest Routine
        if (mod(istep,nstpyear).eq.1) then
            success = celest(.false.)

#ifdef DINSOL
            call dinsol_init(valid_request)
#endif

        endif
        
        ecc=ecc2
        obl=asin(so)
        omweb=(perh+180.00)*deg2rad

#ifdef DINSOL

#if ( F_PALAEO >= 1 )
      annee = -1.0*celest_year
#else
      annee=(iscenyear+(irunlabelf-irunlabeloffset)+iyear)
#endif

        OB = get_orbit_param(annee)
!        write(*,*) "OB%ecc = ", OB%ecc
!        write(*,*) "OB%obl = ", OB%oblq
!        write(*,*) "OB%omweb = ", OB%prcs

! dmr replace DINSOL values in the ECBilt values
       ecc=OB%ecc
       obl=OB%oblq*deg2rad
       omweb=OB%prcs*deg2rad
        
#endif

!~         write(*,*) "INSOL OUTPUT =="
        
!~         write(*,*) "ecc = ", ecc
!~         write(*,*) "obl = ", obl/deg2rad
!~         write(*,*) "omweb = ", omweb/deg2rad
        
!      elseif (iscencel.eq.2) then
!! --- Transient Bretagnon Routine
!        if (mod(istep,nstpyear).eq.1) call bretagnon
!        ecc=ecc2
!        obl=asin(so)
!        omweb=(perh+180.00)*deg2rad
      else
! --- Fixed Values from namelist
        ecc=eccf
        obl=oblf*deg2rad
        omweb=(omwebf+180.00)*deg2rad
      endif

! --- Write Orbital Parameters every year to info
      if (mod(istep,nstpyear).eq.1) then
        write(info_id,'(A,I6,3(X,F12.6))') 'Orbital Parameters: ', irunlabelf+iyear,ecc,obl,omweb
      endif

! --- day of the Vernal Equinox, set at 21 MARCH
      NVE=30+30+21

! First compute alam0 (the starting point). Then convert days to
! the true longitude ala, using Berger s formulae.
! Longitude ala loops from 0 (Vernal Equinox) to 359, ro is earth-sun
! distance relative to the major axis, del is declination.
!
      ala0=0.
      beta=(1.-ecc**2.)**0.5
      fac1=(0.5*ecc+0.125*ecc**3.)*(1.+beta)*sin(ala0-omweb)
      fac2=0.25*ecc**2.*(0.5+beta)*sin(2.*(ala0-omweb))
      fac3=0.125*ecc**3.*(1./3.+beta)*sin(3.*(ala0-omweb))
      alam0=ala0-2.*(fac1-fac2+fac3)


      l=(imonth-1)*30+iday-NVE
      if (l.lt.0) l=l+360
      alam=alam0+l*1.0*pi/180.


      fac1=(2.*ecc-0.25*ecc**3.)*sin(alam-omweb)
      fac2=1.25*ecc**2.*sin(2.*(alam-omweb))
      fac3=(13./12.)*ecc**3.*sin(3.*(alam-omweb))
      ala=alam+fac1+fac2+fac3
      ro=(1.-ecc**2.)/(1.+ecc*cos(ala-omweb))
      deltal=asin(sin(obl)*sin(ala))


      sindl=sin(deltal)
      cosdl=cos(deltal)
      tandl=tan(deltal)


! factor voor variable afstand Aarde-Zon (Berger, p.2367; Velds, p. 99)
      solard=1./ro**2.
!     solardref=1./roref**2.
      solardref=solard
      ksw=0.042
      alpho3sw(1)=0.247
      alpho3sw(2)=0.324
      solarcl(:)=solarm
      indxtsi=1
!     indxtsi=istsistrt-842
      indxvol=1
        if (iscentsi.eq.0.)solartsi(:)=0.
        if ((iscentsi.eq.1).AND.((irunlabelf+iyear).ge.istsistrt) )then
          indxtsi=irunlabelf+iyear-0
          if (indxtsi.gt.2001) indxtsi=2001
          if (iyear.eq.0) indxtsi=indxtsi+1
          solarc=solarm+solartsi(indxtsi)*facttsi
          solarcl(:)=solarc
        if (mod(istep,nstpyear).eq.1) then
         write(info_id,11) 'sol forcing ',iyear,indxtsi,solarcl(15), solartsi(indxtsi)*facttsi
        endif
      endif

      if ((iscenvol.eq.1).AND.((irunlabelf+iyear).ge.isvolstrt) ) then
          indxvol=irunlabelf+iyear-0
          if (indxvol.gt.2001) indxvol=2001
          if (iyear.eq.0) indxvol=1
        do i=1,10
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,1)
        enddo
        do i=11,16
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,2)
        enddo
        do i=17,22
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,3)
        enddo
        do i=23,32
          solarcl(i)=solarcl(i)+solarvol(indxvol,imonth,4)
        enddo
      endif
      if ((isceno3.eq.1).AND.((irunlabelf+iyear).ge.iso3strt))then

        do i=1,16
        solarcl(i)=solarcl(i)+(4*alpho3sw(1)*ksw*(ghg(20)-25.))
        enddo
        do i=17,32
        solarcl(i)=solarcl(i)+(4*alpho3sw(2)*ksw*(ghg(20)-25.))
        enddo
      endif

      if (mod(istep,nstpyear).eq.1) then
       write(info_id,12) 'vol forcing ',iyear,indxvol,solarcl(15), solarvol(indxvol,imonth,2),solarvol(indxvol,imonth,4)
      endif

11      format(A12,2i6,2f12.3)
12      format(A12,2i6,3f12.3)
      do i=1,nlat
! zonneconstante is 1370 in sw parameterisatie
       solarcf(i)=solarcl(i)/1370.d0
! beide effecten samen
       solarf(i)=solarcf(i)*solard
      enddo


      do i=1,nlat
         rkosha1=-tanfi(i)*tandl
         rkosha1=sign(min(abs(rkosha1),1.d0),rkosha1)
         ha1=acos(rkosha1)
         rkosz=sinfi(i)*sindl*(ha1 - tan(ha1))/pi
         if (ha1 .ne. 0.d0) then
            kosz(i) = rkosz*pi/ha1
         else
            kosz(i) = rkosz
         endif
         dayfr(i) = ha1/pi
         q0(i)=rkosz*solarcl(i)*solard

#if ( HOURLY_RAD > 0 )
!=======================================================================
!   loops over the hours
!   currently, use of a 2 hours timestep, 12 steps per day
!   as an input for CARAIB
!   In the future, could be modified to have a more generic approach
!=======================================================================

         if (mod(istep,iatm).eq.0) then

         do hi = 1, n_steps_day


           heure = (float(hi) - 0.5) / (12./24.)
           h = pi * (heure - 12.)/12.

           cz = sinfi(i)*sindl+cosfi(i)*cosdl*dcos(h)

           if(dabs(h) .ge. ha1) then
              solflux = 0.0
           else
              solflux = cz*solarcl(i)*solard
           endif


           solflux_hourly(hi,i) = solflux

!~            if (i.eq.28) then
!~              WRITE(*,'(A,I7, 2F25.5)') "heure", day_of_year, heure, solflux_hourly(hi,i)
!~            endif

        enddo

!~         if (i.eq.28) then
!~            WRITE(*,'(A)') "heure"
!~         endif

        endif
#endif

      enddo

        returnValue = .true.

        return

      end function ec_solar

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: celest
!
!>     @brief Computation of orbital parameters according to Berger, 1978
!
!      DESCRIPTION:
!
!>     ! --- Input:
!! --- ^^^^^^
!! ---   - time elapsed since 1950 (>0 after 1950) (annee, units : kyr)
!! ---
!! --- Output:
!! --- ^^^^^^
!! ---   - longitude of the perihelion measured from the moving vernal
!! ---     point (perh, units : deg)
!! ---   - eccentricity (ecc)
!! ---   - sine of the obliquity (so)
!! ---
!! --- Read:
!! --- ^^^^^
!! ---   File: mbcs2_cor
!! ---   ~~~~~
!! ---     - amplitude for dvlpt of (e-pi) system (cfr. theoretical part)
!! ---       (ae)
!! ---     - mean rate for dvlpt of (e-pi) system (cfr. theoretical part)
!! ---       (y -> be, units : rad.yr-1)
!! ---     - phase for dvlpt of (e-pi) system (cfr. theoretical part)
!! ---       (z -> ce, units : rad)
!! ---
!! ---     - amplitude for dvlpt of obliquity (aob, units : arcsecond)
!! ---     - mean rate for dvlpt of obliquity
!! ---       (y -> bob, units : rad.yr-1)
!! ---     - phase for dvlpt of obliquity
!! ---       (z -> cob, units : rad)
!! ---
!! ---     - amplitude for dvlpt of general precession in longitude
!! ---       (aop, units : arcsecond)
!! ---     - mean rate for dvlpt of general precession in longitude
!! ---       (y -> bop, units : rad.yr-1)
!! ---     - phase for dvlpt of general precession in longitude
!! ---       (z -> cop, units : rad)
!! ---
!! --- Write:
!! --- ^^^^^^
!! ---   File: mout.d (iwr1)
!! ---   ~~~~~
!! ---     - number of terms for dvlpt of (e-pi) system
!! ---       (cfr. theoretical part) (nef)
!! ---     - number of terms for dvlpt of obliquity (nob)
!! ---     - number of terms for dvlpt of general precession in
!! ---       longitude (nop)
!! ---     - climatic precession (pre)
!! ---     - absolute value of year (or number of page-AB)
!! ---       (ipage, units : kyr)
!! ---     - year (ikyr, units : kyr)
!! ---     - eccentricity (ecc)
!! ---     - longitude of the perihelion measured from the moving vernal
!! ---       point (perh, units : deg)
!! ---     - obliquity (xob, units : deg)
!! ---  ?  - caloric equator (xeq, units : deg)
!! ---
!! ---   File: bcmlog (iwr6)
!! ---   ~~~~~
!! ---     - year (ikyr, units : kyr)
!! ---     - eccentricity (ecc)
!! ---     - longitude of the perihelion measured from the moving vernal
!! ---       point (perh, units : deg)
!! ---     - obliquity (xob, units : deg)
!! ---  ?  - caloric equator (xeq, units : deg)
!! ---
!! --- Computation:
!! --- ^^^^^^^^^^^^
!! ---   - pi (pi, units : rad)
!! ---   - conversion factor deg --> rad (pir, units : rad.deg-1)
!! ---   - conversion factor arcsecond --> rad (pirr, units : rad.arcsecond-1)
!! ---   - independent term for dvlpt of obliquity (xod, units : deg)
!! ---   - independent term for dvlpt of general precession in longitude
!! ---     (xop, units : deg)
!! ---   - linear term (precessional constant) for dvlpt of general
!! ---     precession in longitude (prm, units : arcsecond.yr-1)
!! ---   - eccentricity * sine (longitude of perihelion in a fixed
!! ---     reference frame) (xes)
!! ---   - eccentricity * cosine (longitude of perihelion in a fixed
!! ---     reference frame) (xec)
!! ---   - argument of the trigonometrical FUNCTION in the expansions of
!! ---     (e-pi) system or general precession in longitude or obliquity
!! ---     (arg, units : rad)
!! ---   - longitude of perihelion in a fixed reference frame
!! ---     (rp, units : rad; drp, units : deg)
!! ---   - general precession in longitude (prg, units : arcsecond;
!! ---     dprg, units : deg)
!! ---
!! --- Bibliography:
!! --- ^^^^^^^^^^^^^
!! ---   - Berger, A., Long-term variations of daily insolation and
!! ---     Quaternary climatic changes
!! ---     J. Atmos. Sc., 35 (12), 2362-2367, 1978.
!! ---
!! --- Last changes made (by,when):
!! --- ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!! ---   DUTRIEUX Alexis and LOUTRE Marie-France, 05/11/1992
!! ---
!! --- ===================================================================
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function celest(init_flg) result(returnValue)

!~ ###       use comatm,               only:
       use comphys,              only: ecc2, irunlabeloffset, iscenyear, perh, so

       use comemic_mod,          only: iyear
       use comrunlabel_mod,      only: irunlabelf
       use comunit,              only: iuo

#if ( F_PALAEO >= 1 )
       use palaeo_timer_mod,     only: celest_year
#endif

       use global_constants_mod, only: dblp=>dp, ip
       use newunit_mod, only: mbcs2_cor_id, info_id

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in] void
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical, intent(in) :: init_flg

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical                                   :: returnValue
       integer(kind=ip), save                    :: iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
!~        common/wr0/iwr1,iwr2,iwr3,iwr4,iwr5,iwr6,iwr7,iwr8,iwr9,iwr0
       real(kind=dblp), dimension(-16:+16), save :: errtst
!~        common/tst/errtst

       real(kind=dblp), dimension(19) ,save      :: ae, be, ce
       real(kind=dblp), dimension(104),save      :: aob, bob, cob
       real(kind=dblp), dimension(177),save      :: aop, bop, cop
       real(kind=dblp),                save      :: pi1, pir, pirr, annee, y, z, xod, xop, prm, t, xes, xec, arg, rp    &
                      , drp, prg, dprg, pre, xob, xeq
       integer(kind=ip),               save      :: nef, nob, nop, ipage, ikyr, i


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


! --- Future (>0)
! --- Past (<0)

#if ( F_PALAEO >= 1 )
! --- Inconsistency in different files : BP is POSITIVE in palaeo_timer,
! ---  negative here !!!!
      annee = -1.0*celest_year/1000.0
#else

      annee=(iscenyear+(irunlabelf-irunlabeloffset)+iyear)/1000.0d0
!      annee=(-1950.0)/1000.0d0
!      annee=(iscenyear)/1000.0d0
#endif
      print*, "insol", annee, iyear
      write(info_id,'(A,F15.5,A,F15.5,A)')"Celest Year ", 1950 + annee*1000.0," AD == ", sign(annee,real(iscenyear))*(-1.0)*1000.0 &
                                        ," BP"

!~ --- [TODO] Following lines looks pretty much like an init ... do we really need to re-read iuo+38 every single year?
      if (init_flg) then
!*****
! --- 0. Computation of errtst (= needed to resolve error test)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=-16,16,1
         errtst(i)=10.d0**(i)
      enddo
! ----

      pi1 = dacos(-1.0d0)
      pir = pi1 / 180.0d0
      pirr = pir / 3600.0d0
!
! --- 1.earth orbital elements :
! --- ^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      rewind (mbcs2_cor_id)
!
! --- Eccentricity:
! --- ~~~~~~~~~~~~~
!
      nef = 19
      do i=1,nef
        read(mbcs2_cor_id,9000) ae(i),y,z
        be(i) = y * pirr
        ce(i) = z * pir
      enddo
!
! --- Obliquity:
! --- ~~~~~~~~~~
!
      xod = 23.320556d0
      nob = 104
      do i=1,nob
        read(mbcs2_cor_id,9005) aob(i),y,z
        bob(i) = y * pirr
        cob(i) = z * pir
      enddo
!
! --- General precession in longitude:
! --- ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      xop = 3.392506d0
      prm = 50.439273d0
      nop = 177
      do i=1,177
        read(mbcs2_cor_id,9005) aop(i),y,z
        bop(i) = y * pirr
        cop(i) = z * pir
     enddo


     endif
!~ --- [TODO] Previous lines need to be moved to an init function or only call once ...

!
! --- 2.numerical value for ecc pre xob:
! --- ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
      t = annee * 1000.0d0
      xes = 0.0d0
      xec = 0.0d0
      do i=1,nef
        arg = be(i) * t + ce(i)
        xes = xes + ae(i) * dsin(arg)
        xec = xec + ae(i) * dcos(arg)
      enddo
      ecc2 = dsqrt(xes * xes + xec * xec)
      if ((dabs(xec).lt.errtst(-8)).and.(dabs(xes).lt.errtst(-8))) then
        rp = 0.0d0
      else
        rp = datan2(xes,xec)
        if (rp.lt.0.0d0) then
          rp = rp + 2.0d0 * pi1
        endif
      endif
      drp = rp / pir
!
      prg = prm * t
      do i=1,nop
        arg = bop(i) * t + cop(i)
        prg = prg + aop(i) * dsin(arg)
      enddo
      dprg = prg / 3600.0d0 + xop
      perh = drp + dprg
      perh = dmod(perh,360.0d0)
      if (perh.lt.0.0d0) then
        perh = perh + 360.0d0
      endif
!
      pre = ecc2 * dsin(perh * pir)
!
      xob = xod
      do i=1,nob
        arg = bob(i) * t + cob(i)
        xob = xob + aob(i) / 3600.0d0 * dcos(arg)
      enddo
!
      so = dsin(xob * pir)
      xeq = (datan(4.0d0 * pre / (pi1 * so))) / pir
!
      ipage = dabs(t / 1000.)
!
      ikyr = t / 1000.
!
!     *** change by Hans Renssen, June 13 2001
!      write(info_id,*) 'celest: annee, irunlabelf, iyear, ecc2, so, perh'
!      write(info_id,*) annee,irunlabelf,iyear,ecc2,so,perh
        call flush(info_id)
!     *** end change HR
!
 9000 format(13x,f11.8,f20.7,f20.6)
 9005 format(7x,f13.7,2x,f10.6,2x,f10.4)
 9010 format(//,1x,'long term daily insolation',/,1x,'number of terms in:',1x,'eccentricity',i5,2x,'obliquity',i5,2x, &
                'general precession',i5,//)
 9015 format('precession climatique= ',f8.5,5x,'page = ',i4)
 9020 format(1x,'date = ',i6,3x,'eccen = ',f9.6,3x,'long per = ', f7.2,3x,'obliq = ',f7.3,3x,'cal eq = ',f5.2,/)
!

        returnValue = .true.

        return

      end function celest

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module atmrad_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
