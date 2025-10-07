!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of iLOVECLIM/VECODE
!!      iLOVECLIM/VECODE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
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
!      MODULE: veget_sub_mod
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module veget_sub_mod is a fortran90 module port of initial veget_sub.f
!
!>     @date Creation date: October, 17th, 2019
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module veget_sub_mod

       use global_constants_mod, only: dblp=>dp, ip, str_len, stdout

       implicit none
       private

       public :: ccstat, ccstatR, initcpar, ccparam, climpar, ccdyn, ccdynR

      ! NOTE_avoid_public_variables_if_possible

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: initcpar
!
!>     @brief This function is a fortran90 port of the initial 77 subroutine
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function INITCPAR() result(returnValue)

#if ( CYCC == 2 )
        USE veget_iso, ONLY: c13frac,c13frac4
        USE C_res_mod, ONLY: c13init
#endif

        use comatm, only: nlat, nlon
        use veget_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in] void
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical         :: returnValue

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! initialisation of variables
         ades=0.0011
         acr=28
         a=7000.
         bet=0.002
         gamm=0.00017
!        gamm2=0.00025
         gdd0_min=1000.
         gdd0_max=1800.
!jmc0    gdd0_min=1000.
!jmc0    gdd0_max=1700.
         fmax=0.9
         nppmax=1.3
         v1=0.000664
         v2=0.119
         v3=3.73
         c1t=0.046
         c2t=0.58
         c3t=1.6
         c1g=0.069
         c2g=0.38
         c3g=1.6
         d1t=0.22
         d2t=7.19
         d3t=5.5
         d1g=0.6
         d2g=0.41
         d3g=6.0
         e1t=17.9
         e2t=167.3
         e3t=15.
         e1g=0.67
         e2g=50.5
         e3g=100.
         f1t=0.43
         f2t=24.3
         f3t=13.
         f1g=0.34
         f2g=17.8
         f3g=50.
         k2t=1.
         k3t=0.025
         k0t=0.6
         k0g=0.2
         k2g=0.55
         k4g=0.025
         k3g=0.025
         t3g=1.
         t1tn=4
         t1td=1
         deng=20
         dentd=20
         dentn=6
         ps5=0.04
         soilt=5
!dmr --- The following parameters are not declared in VECODE -- CLIMBER
         acwd=100
         acwt=100
         acwg=100
         acwn=100
         zrd=1.
         zrt=1.
         zrg=0.6
         zrn=0.6
         rsd=0.
         rst=300.
         rsg=130.
         rsn=160.
!dmr ---
#if ( CYCC == 2 )
!dmr&nb --- Adding the carbon isotopes initialization ...

!         c14init=100. /* Moved to ccstat_isotopes */
!         c14tdec=1./8240. /* Moved to c14dec in carbone_co2 */
!        the same decay rate of c14 is determined in ocn_bio.f: C14dec

         c13frac=1-18./1000.
         c13frac4=1-5./1000.
         c13init=c13frac
#endif

!dmr --- following was done in CLIMBER, but seems inconsitent with iLOVECLIM
!cnb - initialisation des reservoirs de carbon
!      if (KTVM.eq.2) then
!
!      do i=1,IT
!        do k=1,NS
!
!        b1(i,k)=b1t(i,k)*st(i,k)+b1g(i,k)*sg(i,k)
!        b2(i,k)=b2t(i,k)*st(i,k)+b2g(i,k)*sg(i,k)
!        b3(i,k)=b3t(i,k)*st(i,k)+b3g(i,k)*sg(i,k)
!        b4(i,k)=b4t(i,k)*st(i,k)+b4g(i,k)*sg(i,k)
!
!        enddo
!      enddo
!
!      endif
!cnb ---

!dmr&nb ---

        returnValue = .true.

        return

      end function initcpar

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ccstat
!
!>     @brief This function is a fortran90 port of the initial 77 subroutine
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function CCSTAT(fracgr,darea) result(returnValue)

       use comemic_mod, only: iyear
       use comatm, only: nlat, nlon
       use veget_mod
       use comrunlabel_mod, only: irunlabelf

#if ( FROG_EXP > 0 )
       use veget_mod, only : Fv, Fv_t, Fv_g
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  fracgr   fraction of grass-type plants
!>    @param[in] darea    area of the cells
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real(kind=dblp), dimension(nlat, nlon), intent(in) :: fracgr
       real(kind=dblp), dimension(nlat),       intent(in) :: darea

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical                                            :: returnValue
       real(kind=dblp)                                    :: tempor1
       integer(kind=ip)                                   :: indxv

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!...1) calculation of initial carbon cycle parameters

        returnValue =  CCPARAM()

! calculation of equilibrium storages


!~ [TODO] lat, lon indexes are coming from a global module [UGLY]
!~        change to have it by reference for the function


! leaves biomass
! b1t is leaves phytomass for trees, b1g - for grass (kg C/m2)
! t1t is residence time of carbon in trees, t1g - in grass (years)

        b1t(lat,lon)=k1t*t1t*nppt
        b1g(lat,lon)=k1g*t1g*nppg

!   stems and roots biomass

        b2t(lat,lon)=(1-k1t)*t2t*nppt
        b2g(lat,lon)=(1-k1g)*t2g*nppg

!   litter

        b3t(lat,lon)=(k0t*b1t(lat,lon)/t1t+k2t/t2t*b2t(lat,lon))*t3t
        b3g(lat,lon)=(k0g*b1g(lat,lon)/t1g+k2g/t2g*b2g(lat,lon))*t3g

! mortmass and soil organic matter

        b4t(lat,lon)=(k3t/t3t*b3t(lat,lon))*t4t
        b4g(lat,lon)=(k4g/t2g*b2g(lat,lon)+k3g/t3g*b3g(lat,lon))*t4g
#if ( FROG_EXP > 0 )
        Fv_t(lat,lon)=(k3t/t3t*b3t(lat,lon))*t4t
        Fv_g(lat,lon)=(k4g/t2g*b2g(lat,lon)+k3g/t3g*b3g(lat,lon))*t4g
#endif

! initialization of fraction dynamic variables

!
! until year 1992: distribution min & max is prescribed:
! after year 1992: vegetation is allowed to grow in the grid cell where no deforestation
! takes place in 1992. Vegetation can decline everywhere
!
        indxv=992
        if ((iscendef.eq.1).AND.((irunlabelf+iyear).ge.ivegstrt)) then
         if (irunlabelf+iyear.eq.ivegstrt) st_const(lat,lon)=st(lat,lon)
         tempor1=9.0E+19
         if((irunlabelf+iyear).le.1992) then
          indxv=(irunlabelf+iyear)-1000
          if(indxv.eq.0) indxv=1
          if (farea(lat,lon,indxv).lt.9.0E+19) tempor1=st_const(lat,lon)-farea(lat,lon,indxv)
         else
          if((farea(lat,lon,indxv).lt.9.0E+19).AND.(farea(lat,lon,indxv).gt.0.D0)) tempor1=st_const(lat,lon)-farea(lat,lon,indxv)
         endif
         st(lat,lon)=min(forshare_st,tempor1)
        else
         st(lat,lon)=forshare_st
        endif
        if(st(lat,lon).lt.0.) st(lat,lon)=0.

        sd(lat,lon)=desshare_st
!~ [DEPRECATED]
!~ #if ( ISM == 1 )
!~         if (flgism) then
!~          if (sd(lat,lon).lt.soiltype(lat,lon,ca)) then
!~           sd(lat,lon)=soiltype(lat,lon,ca)
!~           if ((sd(lat,lon)+st(lat,lon)).gt.1.) st(lat,lon)=1-sd(lat,lon)
!~          endif
!~         endif
!~ #endif
!~ [DEPRECATED]
        snlt(lat,lon)=nlshare_st
        if(sd(lat,lon).lt.0.) sd(lat,lon)=0.
        sg(lat,lon)=1.-st(lat,lon)-sd(lat,lon)
        if(sg(lat,lon).lt.0.) sg(lat,lon)=0.

        returnValue = CLIMPAR(fracgr,darea)

        return

      end function ccstat

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ccparam
!
!>     @brief This function is a fortran90 port of the initial 77 subroutine
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function CCPARAM() result(returnValue)

#if ( CYCC == 2 )
        USE veget_iso
#endif

#if ( COMATM == 1 )
        use comatm, only: nlat, nlon
        use veget_mod
#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in] void
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical         :: returnValue
       real(kind=dblp) :: npp1,npp2,avefor,differ,pcr
       real(kind=dblp) :: db1,db2,db3

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


! calculation of current cycle parameters

! potential trees share

       avefor=ave_pr*ave_pr*ave_pr*ave_pr
       differ=gdd0-gdd0_min
       db1=-bet*differ
       db2=gamm*differ
       db3=differ*differ

       if(differ.lt.0) then
          forshare_st=0
       else
          forshare_st=(1-exp(db1))*avefor/(avefor+a*db3*exp(db2))
       endif
       if (forshare_st.gt.fmax) forshare_st=fmax

! potential desert share - desshare_st

       desshare_st=0

! northern deserts

       if(gdd0.lt.100) desshare_st=1

       if(gdd0.ge.100.and.gdd0.lt.gdd0_min) desshare_st=(gdd0_min-gdd0)/(gdd0_min-100.)

! southern deserts

          if (gdd0.ge.gdd0_max) then

            pcr=acr*exp(gamm2/2.*differ)

            if (ave_pr05.le.pcr) then
                desshare_st=1
                forshare_st=0
            else
                db2=(ave_pr05-pcr)/exp(gamm2*differ)
                desshare_st=1.03/(1+ades*db2*db2)-0.03
                if (desshare_st.lt.0) desshare_st=0
            endif

         endif

! calculation of NPP, Lieth's formula

        db1=-v1*ave_pr
!        if (gdd0.ge.gdd0_max) db1=-v1*ave_pr05
        db2=-v2*ave_t
        npp1=(1.-exp(db1))
        npp2=1./(1.+v3*exp(db2))
        if(npp1.lt.npp2) then
                npp=nppmax*npp1
            else
                npp=nppmax*npp2
        endif

! CO2 enrichment factor
!
!-AM (may 2007)
! Modified in order to account for different enhancement factors for tree and grass
!  -> betat & betag
! also division by log(2) is now included in beta_x (see veget.f)
!       npp=npp*(1.0+((0.25/LOG(2.))*LOG(co2ghg/280.)))
        nppt=npp*(1.0+(betat*LOG(co2ghg/280.)))
        nppg=npp*(1.0+(betag*LOG(co2ghg/280.)))

! allocation factors and residence time of leaves biomass

        k1t=c1t+c2t/(1+c3t*nppt)
        k1g=c1g+c2g/(1+c3g*nppg)

        t1t=d1t+d2t/(1+d3t*nppt)
        t1g=d1g+d2g/(1+d3g*nppg)

!   residence time of stems and roots biomass

        t2t=e1t+e2t/(1+e3t*nppt)
        t2g=e1g+e2g/(1+e3g*nppg)

!   residence time of fast carbon pool

        t3t=16.*exp(-ps5*(ave_t-soilt))
        t3g=40.*exp(-ps5*(ave_t-soilt))

! residence time of slow soil organic matter

        t4t=900.*exp(-ps5*(ave_t-soilt))
        t4g=t4t

#if ( PF_CC > 0 )

!dmr --- should be used only on permafrost areas

!kc **tuneable residence times**
!kc new residence time for permafrost fast soil
!kc boost South edges soil carbon

        t6t=t3t*(60*fr_ndx+80)
        t6g=t3g*(60*fr_ndx+80)

!kc new residence time for permafrost slow soil

!        t5t=t4t
        t5t=t4t*(0.1*fr_ndx+0.1)
        t5g=t5t

!print *,'t3t',t3t,'t3g',t3g,'t4t',t4t
#endif

!calculation of potential nedleleaves trees ratio

        nlshare_st=(t1t-t1td)/(t1tn-t1td)
        if (nlshare_st.gt.1) nlshare_st=1
        if (nlshare_st.lt.0) nlshare_st=0

#if ( CYCC == 2 )
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul des C4

!dmr&nb A_REPRENDRE

! calculation of c4 grass fraction
        !! [BUG] TATMSMIN NEVER INITIALIZED !!!
        TATMSMIN(lat,lon) = 0.0d0

!        if(TATMSMIN(lat,lon).lt.14) then
        if(TATMSMIN(lat,lon).lt.12) then
            g4share_st=0
        else
!          if(TATMSMIN(lat,lon).lt.17) then
!           g4share_st=1-(17-TATMSMIN(lat,lon))/(17.-14.)
          if(TATMSMIN(lat,lon).lt.15.5) then
           g4share_st=1-(15.5-TATMSMIN(lat,lon))/3.5
        else
           g4share_st=1
        endif
        endif

!dmr&nb A_REPRENDRE

!dmr&nb --- Fin ajout du calcul des C4
#endif

        returnValue = .true.

        return

      end function ccparam

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ccparam
!
!>     @brief This function is a fortran90 port of the initial 77 subroutine
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function climpar(fracgr,darea) result(returnValue)

#if ( CYCC == 2 )
        USE veget_iso
#endif

#if ( COMATM == 1 )
        use comatm, only: nlat,nlon
        use veget_mod
#endif

#if ( FROG_EXP > 0 )
       use veget_mod, only : Fv, Fv_t, Fv_g
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  fracgr   fraction of grass-type plants
!>    @param[in] darea    area of the cells
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real(kind=dblp), dimension(nlat, nlon), intent(in) :: fracgr
       real(kind=dblp), dimension(nlat),       intent(in) :: darea

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical         :: returnValue
       real(kind=dblp) :: tempor1

#if ( CYCC == 2 )
!dmr&nb --- 31 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14

!dmr&nb A_REPRENDRE
       real(kind=dblp) :: tempor2=0.0_dblp
#endif

       integer(ip)     :: KTVM

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         KTVM=2
! calculation of annual averaged LAI - lai

         laig=b1g(lat,lon)*deng
         lait=b1t(lat,lon)*(dentn*snlt(lat,lon)+dentd*(1-snlt(lat,lon)))

         BLAI(lat,lon,1)=lait
         BLAI(lat,lon,2)=laig

! calculation of annual carbon uptake

        if (KTVM.eq.2) then
           tempor1=b1(lat,lon)+b2(lat,lon)+b3(lat,lon)+b4(lat,lon)

#if ( CYCC ==2 )
!dmr&nb --- 31 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14

!dmr&nb A_REPRENDRE
          tempor2=bc13(lat,lon)
#endif
        else
           tempor1=0.D0
        endif

        b1(lat,lon)=b1t(lat,lon)*st(lat,lon)+b1g(lat,lon)*sg(lat,lon)
        b2(lat,lon)=b2t(lat,lon)*st(lat,lon)+b2g(lat,lon)*sg(lat,lon)
        b3(lat,lon)=b3t(lat,lon)*st(lat,lon)+b3g(lat,lon)*sg(lat,lon)
        b4(lat,lon)=b4t(lat,lon)*st(lat,lon)+b4g(lat,lon)*sg(lat,lon)
        b12(lat,lon)=b1(lat,lon)+b2(lat,lon)
        b34(lat,lon)=b3(lat,lon)+b4(lat,lon)

#if ( FROG_EXP > 0 )
        Fv(lat,lon)=Fv_t(lat,lon)*st(lat,lon)+Fv_g(lat,lon)*sg(lat,lon)
#endif

#if ( CYCC ==2 )
!dmr&nb --- 31 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14

!dmr&nb A_REPRENDRE

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        bc14(lat,lon)=(b1t14(lat,lon)+b2t14(lat,lon)+b3t14(lat,lon)+b4t14(lat,lon))*st(lat,lon)                                 &
                     +(b1g14(lat,lon)+b2g14(lat,lon)+b3g14(lat,lon)+b4g14(lat,lon))*sg(lat,lon)

!~        print *,lat,lon,b1t14(lat,lon),b1g14(lat,lon),b4t14(lat,lon)

        bc13(lat,lon)=(b1t13(lat,lon)+b2t13(lat,lon)+b3t13(lat,lon)+b4t13(lat,lon))*st(lat,lon)                                 &
                     +(b1g13(lat,lon)+b2g13(lat,lon)+b3g13(lat,lon)+b4g13(lat,lon))*sg(lat,lon)

!dmr&nb --- Fin ajout du calcul isotopes C13, C14
#endif

        anup(lat,lon)=(b1(lat,lon)+b2(lat,lon)+b3(lat,lon)+b4(lat,lon)-tempor1)
        stock(lat,lon)=b1(lat,lon)+b2(lat,lon)+b3(lat,lon)+b4(lat,lon)
#if ( CYCC == 1 )
        stockloch(lat,lon)=fracgr(lat,lon)*darea(lat)*stock(lat,lon)*1E-12
        anuploch(lat,lon)=fracgr(lat,lon)*darea(lat)*anup(lat,lon)*1E-12
        anuploch(lat,lon)=fco2veg*anuploch(lat,lon)
#elif ( CYCC ==2 )
!dmr&nb --- 31 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14

!dmr&nb A_REPRENDRE

        anup13(lat,lon)=bc13(lat,lon)-tempor2
!dmr&nb --- Fin ajout du calcul isotopes C13, C14
#endif

!...      NET PRIMARY PRODUCTION

!       pnpp(lat,lon)=npp*(1-sd(lat,lon))
        pnpp(lat,lon)=nppt*st(lat,lon)+nppg*sg(lat,lon)

        returnValue = .true.

        return

      end function climpar

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ccdyn
!
!>     @brief This function is a fortran90 port of the initial 77 subroutine
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ccdyn(fracgr,darea) result(returnValue)

#if ( CYCC == 2 )
       use veget_iso, ONLY : B4T14, B3T14, B4T13, B3T13, B2T14, B1T14, B2T13, B1T13, B4G14, B3G14, B4G13, B3G13, B2G14, B1G14  &
                    , B2G13, B1G13, BC13, BC14, G4SHARE_ST, SG4, C13FRAC, C13FRAC4
       use carbone_co2, ONLY: C14ATM, C14DEC, PA_C
       use mod_sync_time, ONLY: KENDY
       use C_res_mod, ONLY: c13atm
#endif

       use comemic_mod, only: iyear
       use veget_mod
       use comatm, only: nlat, nlon
       use comrunlabel_mod, only: irunlabelf

#if ( FROG_EXP > 0 )
       use veget_mod, only : Fv, Fv_t, Fv_g
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  fracgr   fraction of grass-type plants
!>    @param[in] darea    area of the cells
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real(kind=dblp), dimension(nlat, nlon), intent(in) :: fracgr
       real(kind=dblp), dimension(nlat),       intent(in) :: darea

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical         :: returnValue

! temporal var*/
       integer(kind=ip):: indxv
       real(kind=dblp) :: tempor1,tempor2,db2,fd,dst,dd,nld,dstime, dsd,temp_sg,temp_st

#if ( CYCC == 2 )
!dmr&nb 12 Novembre 2009
!dmr&nb --- Ajout du calcul des isotopes C13, C14
        real(kind=dblp):: tempor3,tempor4,tempor5,tempor6
!dmr&nb --- Fin ajout du calcul des isotopes C13, C14
        real(kind=dblp):: G4D
#endif

#if ( PF_CC > 0 )
        real(kind=dblp):: b4t_hold, b4g_hold, b3t_hold, b3g_hold
        real(kind=dblp):: w_frac, p_frac
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! calculation of current carbon cycle parameters

        returnValue = ccparam()


#if( PF_CC > 0 )

       p_frac = pf_percent(lat,lon)
       w_frac = 1.0 - pf_percent(lat,lon)
#endif

! calculation of fraction dynamic variables

        fd=forshare_st-st(lat,lon)
        dd=desshare_st-sd(lat,lon)
        nld=nlshare_st-snlt(lat,lon)
#if ( CYCC == 2 )
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul des C4

        g4d=g4share_st-sg4(lat,lon)
!dmr&nb --- Fin ajout du calcul des C4
#endif

        temp_st=st(lat,lon)
        temp_sg=sg(lat,lon)


! calculation of forest dynamics; exponential filtre
        if (abs(fd).lt.100.*tiny(fd)) then
          fd = 0.0d0
          dst = 0.0d0
        else
          dst=fd*(1.d0-exp(-1./t2t))
        endif

        st(lat,lon)=temp_st+dst
        if (st(lat,lon).lt.0.) st(lat,lon)=0.

! desert dynamics; exponential filtre
        dsd=dd*(1.-exp(-1./t2g))
! calculation of characteristic time of desert propagation
        tempor1=sd(lat,lon)+dsd+st(lat,lon)
        if (tempor1.gt.0.9) then
             dstime=(t2g*(1.-tempor1)+t2t*(tempor1-0.9))*10.
             dsd=dd*(1.-exp(-1./dstime))
        endif
        sd(lat,lon)=sd(lat,lon)+dsd
        if (sd(lat,lon).lt.0.) sd(lat,lon)=0.

#if ( ISM == 1 )
        if (flgism) then
         if (sd(lat,lon).lt.soiltype(lat,lon,ca)) then
           sd(lat,lon)=soiltype(lat,lon,ca)
           if ((sd(lat,lon)+st(lat,lon)).gt.1.) st(lat,lon)=1.-sd(lat,lon)
           if (st(lat,lon).lt.0.) st(lat,lon)=0.
         endif
        endif
#endif
!
!
        indxv=irunlabelf+iyear-i0dfor
        tempor1=0.
!
! ----------------------
! Constant vegetation:
! ----------------------
        if ((iscendef.eq.-1).AND.(indxv.ge.0)) then
! Defines ref. tree and desert fractions:
         if (indxv.eq.0) then
          st_const(lat,lon)=st(lat,lon)
          sd_const(lat,lon)=sd(lat,lon)
         else
          st(lat,lon)=st_const(lat,lon)
          sd(lat,lon)=sd_const(lat,lon)
         endif
        endif
!
! Available fraction for vegetation
        tempor2=1.-sd(lat,lon)
        if(tempor2.lt.0.) tempor2=0.
!
! ----------------------
! Deforestation scenario:
! ----------------------
! Forests/trees are replaced with grassland (cropland) in accordance with R&F 1700-1992 scenario.
!
! Once the end of the deforestation scenario file is reached the land use is kept at its state
!  as in the last year.
!
! Three versions: version A1 & A2 (AM & MFL, july 2008);
!                 version B (VB & AM, sept. 2008)
!                 version C (ED & AM for MILMO);
!
! SCENARIO version A ================= goes from here ====>>
!
!        if ((iscendef.eq.1).AND.(indxv.gt.0)) then
!         if(indxv.gt.ndfor) indxv=ndfor
! Version A1
! Scenario conc/efor
! actual fraction of trees is the minimum of
!     [potential fraction , fraction not occupied by cultures]
!         if(farea(lat,lon,indxv).lt.9.0E+19)
!     &       tempor1=tempor2-farea(lat,lon,indxv)
!         if(tempor1.lt.0.) tempor1=0.
!         st(lat,lon)=min(st(lat,lon),tempor1)
! Version A2
! Scenario Conc/Efor
! the tree fraction is reduced by the crop fraction
!         if(farea(lat,lon,indxv).lt.9.0E+19)
!     &    st(lat,lon)=st(lat,lon)-farea(lat,lon,indxv)
!         if(st(lat,lon).lt.0.) st(lat,lon)=0.
!        endif
!
! SCENARIO version A ======================== to there ====<<
!
! If replacing one version by the other be careful to comment and decomment
!  all instructions comprised between the tags "from here" and "to there".
!
! SCENARIO version B ================= goes from here ====>>
!
! Scenario as implemented by Victor Brovkin,
!  complies with the rules edicted in EMIC intercomparison project.
! Deforestation is computed with respect to a reference vegetation distribution;
! the reference distribution corresponds to that as in year ivgstrt.
! It goes like this:
!  - crop = 0 => sd=sd_ref; st=st_ref & sg=1-st-sd=sg_ref;
!  - crop > 0 => sd=sd_ref; st=max(0,st_ref-crop) ; sg=1-sd-st;
!
        if ((iscendef.eq.1).AND.(indxv.ge.0)) then
! Defines ref. tree and desert fractions:
         if (indxv.eq.0) then
          st_const(lat,lon)=st(lat,lon)
          sd_const(lat,lon)=sd(lat,lon)
         else
          if(indxv.gt.ndfor) indxv=ndfor
          sd(lat,lon)=sd_const(lat,lon)
          tempor2=1.-sd_const(lat,lon)
          st(lat,lon)=st_const(lat,lon)
          if(farea(lat,lon,indxv).lt.9.0E+19) then
           st(lat,lon)=st_const(lat,lon)-farea(lat,lon,indxv)
          endif
          if(st(lat,lon).lt.0.) st(lat,lon)=0.
         endif
        endif
!
! SCENARIO version B ======================== to there ====<<
!
! If replacing one version by the other be careful to comment and decomment
!  all instructions comprised between the tags "from here" and "to there".
!
! SCENARIO version C ================= goes from here ====>>
!
! Scenario as formerly written (MILMO)
! The tree fraction is the mimimum among:
!      - potential tree fraction
!      - reference tree-fraction (st_const) less the crop fraction (farea)
! This is not exactly the scenario preconised by V.B. for EMICs
!        (since veg. cover does not necessarily correspond to the reference cover: st smaller, sd evolves...)
! It may produce deforestation fluxes even if farea=0 (if st < st_const)
!
! this part NEEDs TO BE VERIFIED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        if ((iscendef.eq.1).AND.(indxv.ge.0)) then
!         if (indxv.eq.0) then
! Defines ref. tree fraction:
!          st_const(lat,lon)=st(lat,lon)
!         else
!          if(indxv.le.ndfor) then
!           if(farea(lat,lon,indxv).lt.9.0E+19)
!     &       tempor1=st_const(lat,lon)-farea(lat,lon,indxv)
!          else
!           indxv=ndfor
!           if((farea(lat,lon,indxv).lt.9.0E+19).and.(farea(lat,lon,indxv).gt.0.D0))
!     &       tempor1=st_const(lat,lon)-farea(lat,lon,indxv)
!          endif
!          if(tempor1.lt.0.) tempor1=0.
!          st(lat,lon)=min(st(lat,lon),tempor1)
!         endif
!
! SCENARIO version C ======================== to there ====<<
!
!
! do not forget to update dst (needed for corrections below)
        dst=st(lat,lon)-temp_st

        sg(lat,lon)=tempor2-st(lat,lon)
        if (sg(lat,lon).lt.0) sg(lat,lon)=0.
#if ( CYCC == 2 )
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul des C4

        sg4(lat,lon)=g4share_st-g4d*exp(-1./t2g)
!dmr&nb --- Fin ajout du calcul des C4
#endif

        snlt(lat,lon)=nlshare_st-nld*exp(-1./t2t)

! calculation of dynamics of storages

! calculation of changes of storages due to conservation law

! correction for trees

        tempor1=b4t(lat,lon)
        tempor2=b3t(lat,lon)
#if ( CYCC == 2 )
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14

        tempor3=b4t14(lat,lon)
        tempor4=b3t14(lat,lon)
        tempor5=b4t13(lat,lon)
        tempor6=b3t13(lat,lon)
!dmr&nb --- Fin ajout du calcul isotopes C13, C14
#endif

        if(st(lat,lon).gt.0) then
!dmr&nb --- Si les arbres avancent ...
         if(dst.gt.0) then
          if(dst.le.temp_sg) then !cnb si la progression est inferieur a espace d herbes
            b4t(lat,lon)=(b4t(lat,lon)*temp_st+b4g(lat,lon)*dst)/st(lat,lon)
            b3t(lat,lon)=(b3t(lat,lon)*temp_st+b3g(lat,lon)*dst)/st(lat,lon)
#if ( CYCC == 2 )
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14

            b4t14(lat,lon)=(b4t14(lat,lon)*temp_st+b4g14(lat,lon)*dst)/st(lat,lon)
            b3t14(lat,lon)=(b3t14(lat,lon)*temp_st+b3g14(lat,lon)*dst)/st(lat,lon)

            b4t13(lat,lon)=(b4t13(lat,lon)*temp_st+b4g13(lat,lon)*dst)/st(lat,lon)
            b3t13(lat,lon)=(b3t13(lat,lon)*temp_st+b3g13(lat,lon)*dst)/st(lat,lon)

!dmr&nb --- Fin ajout du calcul isotopes C13, C14
#endif
          else !cnb sinon
            b4t(lat,lon)=(b4t(lat,lon)*temp_st+b4g(lat,lon)*temp_sg)/st(lat,lon) !cnb
            b3t(lat,lon)=(b3t(lat,lon)*temp_st+b3g(lat,lon)*temp_sg)/st(lat,lon) !cnb
#if ( CYCC == 2 )
!dmr&nb --- 12 novembre 2009 !cnb
!dmr&nb --- Ajout du calcul isotopes C13, C14 !cnb

            b4t14(lat,lon)=(b4t14(lat,lon)*temp_st+b4g14(lat,lon)*temp_sg)/st(lat,lon) !cnb
            b3t14(lat,lon)=(b3t14(lat,lon)*temp_st+b3g14(lat,lon)*temp_sg)/st(lat,lon) !cnb

            b4t13(lat,lon)=(b4t13(lat,lon)*temp_st+b4g13(lat,lon)*temp_sg)/st(lat,lon) !cnb
            b3t13(lat,lon)=(b3t13(lat,lon)*temp_st+b3g13(lat,lon)*temp_sg)/st(lat,lon) !cnb
#endif
          endif !cnb
         endif
         b2t(lat,lon)=b2t(lat,lon)*temp_st/st(lat,lon)
         b1t(lat,lon)=b1t(lat,lon)*temp_st/st(lat,lon)
#if ( CYCC == 2 )
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14

         b2t14(lat,lon)=b2t14(lat,lon)*temp_st/st(lat,lon)
         b1t14(lat,lon)=b1t14(lat,lon)*temp_st/st(lat,lon)

         b2t13(lat,lon)=b2t13(lat,lon)*temp_st/st(lat,lon)
         b1t13(lat,lon)=b1t13(lat,lon)*temp_st/st(lat,lon)
!dmr&nb --- Fin ajout du calcul isotopes C13, C14
#endif
        endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! correction for grass
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       if(sg(lat,lon).gt.0) then
!dmr&nb --- Si les arbres avancent ...
               if(dst.gt.0) then
                if(dst.le.temp_sg) then !cnb si la progression des arbres est inferieur a l espace d herbes
                b4g(lat,lon)=b4g(lat,lon)*(temp_sg-dst)/sg(lat,lon)
                b3g(lat,lon)=b3g(lat,lon)*(temp_sg-dst)/sg(lat,lon)
!nb                if (b4g(lat,lon).lt.0) b4g(lat,lon)=0 !cnb
!nb                if (b3g(lat,lon).lt.0) b3g(lat,lon)=0 !cnb
#if ( CYCC == 2 )
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14

                 b4g14(lat,lon)=b4g14(lat,lon)*(temp_sg-dst)/sg(lat,lon)
                 b3g14(lat,lon)=b3g14(lat,lon)*(temp_sg-dst)/sg(lat,lon)

                 b4g13(lat,lon)=b4g13(lat,lon)*(temp_sg-dst)/sg(lat,lon)
                 b3g13(lat,lon)=b3g13(lat,lon)*(temp_sg-dst)/sg(lat,lon)

!dmr&nb --- Fin ajout du calcul isotopes C13, C14
#endif
               else !cnb si les arbres progressent au dela des herbes y a plus d herbes
                 b4g(lat,lon)=0.0
                 b3g(lat,lon)=0.0
#if ( CYCC == 2 )
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14
                 b4g14(lat,lon)=0.0
                 b3g14(lat,lon)=0.0

                 b4g13(lat,lon)=0.0
                 b3g13(lat,lon)=0.0
!dmr&nb --- Fin ajout du calcul isotopes C13, C14
#endif
                endif !cnb si progression arbres inferieur a espace d herbes

               else
!dmr&nb --- Si les arbres n'avancent pas (reculent)
          b4g(lat,lon)=(b4g(lat,lon)*temp_sg-tempor1*dst)/sg(lat,lon)
          b3g(lat,lon)=(b3g(lat,lon)*temp_sg-tempor2*dst)/sg(lat,lon)
#if ( CYCC == 2 )
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14

           b4g14(lat,lon)=(b4g14(lat,lon)*temp_sg-tempor3*dst)/sg(lat,lon)
           b3g14(lat,lon)=(b3g14(lat,lon)*temp_sg-tempor4*dst)/sg(lat,lon)

           b4g13(lat,lon)=(b4g13(lat,lon)*temp_sg-tempor5*dst)/sg(lat,lon)
           b3g13(lat,lon)=(b3g13(lat,lon)*temp_sg-tempor6*dst)/sg(lat,lon)

!dmr&nb --- Fin ajout du calcul isotopes C13, C14
#endif
               endif
        b2g(lat,lon)=b2g(lat,lon)*temp_sg/sg(lat,lon)
        b1g(lat,lon)=b1g(lat,lon)*temp_sg/sg(lat,lon)
#if ( CYCC == 2 )
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14

         b2g14(lat,lon)=b2g14(lat,lon)*temp_sg/sg(lat,lon)
         b1g14(lat,lon)=b1g14(lat,lon)*temp_sg/sg(lat,lon)

         b2g13(lat,lon)=b2g13(lat,lon)*temp_sg/sg(lat,lon)
         b1g13(lat,lon)=b1g13(lat,lon)*temp_sg/sg(lat,lon)

!dmr&nb --- Fin ajout du calcul isotopes C13, C14
#endif
       endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! slow soil organic matter
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( PF_CC > 0 )

!kc place holder for b4 terms

        b4t_hold=b4t(lat,lon)
        b4g_hold=b4g(lat,lon)
!        b4t14_hold=b4t14(lat,lon)
!        b4g14_hold=b4g14(lat,lon)
!        b4t13_hold=b4t13(lat,lon)
!        b4g13_hold=b4g13(lat,lon)

!kc for non permafrost effected fraction

        b4t(lat,lon)=w_frac*(b4t_hold+k3t/t3t*b3t(lat,lon)-b4t_hold/t4t)
        b4g(lat,lon)=w_frac*(b4g_hold+k4g/t2g*b2g(lat,lon)+k3g/t3g*b3g(lat,lon)-b4g_hold/t4g)

!kc for permafrost fraction of cell for b4 and c14 and c13

        IF ( p_frac.GT.0.0) THEN
        b5t(lat,lon)=p_frac*(b4t_hold+k3t/t3t*b3t(lat,lon)-b4t_hold/t5t)
        b5g(lat,lon)=p_frac*(b4g_hold+k4g/t2g*b2g(lat,lon)+k3g/t3g*b3g(lat,lon)-b4g_hold/t5g)
        ENDIF
#else


        b4t(lat,lon)=b4t(lat,lon)+k3t/t3t*b3t(lat,lon)-b4t(lat,lon)/t4t
        b4g(lat,lon)=b4g(lat,lon)+k4g/t2g*b2g(lat,lon)+k3g/t3g*b3g(lat,lon)-b4g(lat,lon)/t4g

#if ( FROG_EXP > 0 )
        Fv_t(lat,lon)=Fv_t(lat,lon)+k3t/t3t*b3t(lat,lon) !-b4t(lat,lon)/t4t
        Fv_g(lat,lon)=Fv_g(lat,lon)+k4g/t2g*b2g(lat,lon)+k3g/t3g*b3g(lat,lon) !-b4g(lat,lon)/t4g

#endif


#endif
#if ( CYCC == 2 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        b4t14(lat,lon)=b4t14(lat,lon)+k3t/t3t*b3t14(lat,lon)-b4t14(lat,lon)/t4t
        b4g14(lat,lon)=b4g14(lat,lon)+k4g/t2g*b2g14(lat,lon)+k3g/t3g*b3g14(lat,lon)-b4g14(lat,lon)/t4g

        b4t13(lat,lon)=b4t13(lat,lon)+k3t/t3t*b3t13(lat,lon)-b4t13(lat,lon)/t4t
        b4g13(lat,lon)=b4g13(lat,lon)+k4g/t2g*b2g13(lat,lon)+k3g/t3g*b3g13(lat,lon)-b4g13(lat,lon)/t4g

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr&nb --- Fin ajout du calcul isotopes C13, C14
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#endif

#if( PF_CC > 0 )

!kc add permafrost carbon back in to b4 components

        b4t(lat,lon)=b4t(lat,lon)+b5t(lat,lon)
        b4g(lat,lon)=b4g(lat,lon)+b5g(lat,lon)

#endif

!   fast soil organic matter

#if( PF_CC > 0 )

!kc place holder for fast soil terms
        b3t_hold=b3t(lat,lon)
        b3g_hold=b3g(lat,lon)
!        b3t14_hold=b3t14(lat,lon)
!        b3g14_hold=b3g14(lat,lon)
!        b3t13_hold=b3t13(lat,lon)
!        b3g13_hold=b3g13(lat,lon)

!kc non-permafrost fraction of cell, b3, permafrost effected b6

        b3t(lat,lon)=w_frac*(b3t_hold+b1t(lat,lon)/t1t*k0t+k2t/t2t*b2t(lat,lon)-b3t_hold/t3t)
        b3g(lat,lon)=w_frac*(b3g_hold+b1g(lat,lon)/t1g*k0g+k2g/t2g*b2g(lat,lon)-b3g_hold/t3g)


        IF ( p_frac.GT.0.0) THEN
        b6t(lat,lon)=p_frac*(b3t_hold+b1t(lat,lon)/t1t*k0t+k2t/t2t*b2t(lat,lon)-b3t_hold/t6t)
        b6g(lat,lon)=p_frac*(b3g_hold+b1g(lat,lon)/t1g*k0g+k2g/t2g*b2g(lat,lon)-b3g_hold/t6g)
        ENDIF


#else
        b3t(lat,lon)=b3t(lat,lon)+b1t(lat,lon)/t1t*k0t+k2t/t2t*b2t(lat,lon)-b3t(lat,lon)/t3t
        b3g(lat,lon)=b3g(lat,lon)+b1g(lat,lon)/t1g*k0g+k2g/t2g*b2g(lat,lon)-b3g(lat,lon)/t3g

#endif

#if ( CYCC == 2 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        b3t14(lat,lon)=b3t14(lat,lon)+b1t14(lat,lon)/t1t*k0t+k2t/t2t*b2t14(lat,lon)-b3t14(lat,lon)/t3t

        b3g14(lat,lon)=b3g14(lat,lon)+b1g14(lat,lon)/t1g*k0g+k2g/t2g*b2g14(lat,lon)-b3g14(lat,lon)/t3g

        b3t13(lat,lon)=b3t13(lat,lon)+b1t13(lat,lon)/t1t*k0t+k2t/t2t*b2t13(lat,lon)-b3t13(lat,lon)/t3t

        b3g13(lat,lon)=b3g13(lat,lon)+b1g13(lat,lon)/t1g*k0g+k2g/t2g*b2g13(lat,lon)-b3g13(lat,lon)/t3g

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr&nb --- Fin ajout du calcul isotopes C13, C14
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#endif

#if ( PF_CC > 0 )
!kc sum b3 components back together

        b3t(lat,lon)=b3t(lat,lon)+b6t(lat,lon)
        b3g(lat,lon)=b3g(lat,lon)+b6g(lat,lon)

!        IF ((lat.EQ.27).AND.(lon.EQ.11)) then
!        WRITE(*,*) "/ /", b4t(lat,lon)+b4g(lat,lon),
!     &  p_frac, lat, lon
!        ENDIF
!        PAUSE

!        b3t14(lat,lon)=b3t14(lat,lon)+b6t14(lat,lon)
!        b3g14(lat,lon)=b3g14(lat,lon)+b6g14(lat,lon)
!        b3t13(lat,lon)=b3t13(lat,lon)+b6t13(lat,lon)
!        b3g13(lat,lon)=b3g13(lat,lon)+b6g13(lat,lon)

!kc Soil respiration

!      rsoil_g(lat,lon)=w_frac*((-b3g_hold/t3g)+(-b4g_hold/t4g))
!     >                +p_frac*((-b3g_hold/t6g)+(-b4g_hold/t5g))
!      rsoil_t(lat,lon)=w_frac*((-b3g_hold/t3t)+(-b4g_hold/t4t))
!     >                +p_frac*((-b3g_hold/t6t)+(-b4g_hold/t5t))

!kc Veg respiration

!       rveg_g(lat,lon)=-b2g(lat,lon)/t2g
!       rveg_t(lat,lon)=(-b1t(lat,lon)/t1t)+(-b2t(lat,lon)/t2t)

#endif



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! leaves biomass
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        b1t(lat,lon)=b1t(lat,lon)+k1t*nppt-b1t(lat,lon)/t1t
        b1g(lat,lon)=k1g*nppg*t1g

#if ( CYCC == 2 )
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14

        b1t14(lat,lon)=b1t14(lat,lon)+k1t*nppt*(c14atm/PA_C)-b1t14(lat,lon)/t1t

          IF(KENDY.EQ.1) then
            IF(lat.eq.1.and.lon.eq.1) WRITE(*,*) 'ccdyn.f Flux al C14'
          ENDIF

        b1g14(lat,lon)=k1g*nppg*(c14atm/PA_C)*t1g


        b1t13(lat,lon)=b1t13(lat,lon)+k1t*nppt*c13atm*c13frac-b1t13(lat,lon)/t1t


        b1g13(lat,lon)=k1g*nppg*c13atm*(c13frac*(1-sg4(lat,lon))+c13frac4*sg4(lat,lon))*t1g

!dmr&nb --- Fin ajout du calcul isotopes C13, C14
#endif

!   stems and roots biomass

        b2t(lat,lon)=b2t(lat,lon)+(1-k1t)*nppt-b2t(lat,lon)/t2t
        b2g(lat,lon)=b2g(lat,lon)+(1-k1g)*nppg-b2g(lat,lon)/t2g

#if ( CYCC == 2 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        b2t14(lat,lon)=b2t14(lat,lon)+(1-k1t)*nppt*(c14atm/PA_C)-b2t14(lat,lon)/t2t

        b2g14(lat,lon)=b2g14(lat,lon)+(1-k1g)*nppg*(c14atm/PA_C)-b2g14(lat,lon)/t2g

        b2t13(lat,lon)=b2t13(lat,lon)+(1-k1t)*nppt*c13atm*c13frac-b2t13(lat,lon)/t2t

        b2g13(lat,lon)=b2g13(lat,lon)+(1-k1g)*nppg*c13atm*(c13frac*(1-sg4(lat,lon))+c13frac4*sg4(lat,lon))-b2g13(lat,lon)/t2g

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr&nb --- Fin ajout du calcul isotopes C13, C14
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!dmr&nb --- 12 novembre 2009
!dmr&nb --- Ajout du calcul isotopes C13, C14
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   c14 annual decay
        IF(KENDY.EQ.1) then
          b1t14(lat,lon)=b1t14(lat,lon)*(1-c14dec)
          b2t14(lat,lon)=b2t14(lat,lon)*(1-c14dec)
          b3t14(lat,lon)=b3t14(lat,lon)*(1-c14dec)
          b4t14(lat,lon)=b4t14(lat,lon)*(1-c14dec)
          b1g14(lat,lon)=b1g14(lat,lon)*(1-c14dec)
          b2g14(lat,lon)=b2g14(lat,lon)*(1-c14dec)
          b3g14(lat,lon)=b3g14(lat,lon)*(1-c14dec)
          b4g14(lat,lon)=b4g14(lat,lon)*(1-c14dec)
          IF(lat.eq.1.and.lon.eq.1)  WRITE(*,*) 'ccdyn.f lnd C14 decay'
        ENDIF

!dmr&nb --- Fin ajout du calcul isotopes C13, C14
#endif

        returnValue = climpar(fracgr,darea)

        returnValue = .true.

        return

      end function ccdyn

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ccstatR
!
!>     @brief This function is a fortran90 port of the initial 77 subroutine
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ccstatR(fracgr,darea) result(returnValue)

       use comemic_mod, only:
       use comatm, only: nlat, nlon
       use veget_mod
       use comrunlabel_mod, only: irunlabelf

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  fracgr   fraction of grass-type plants
!>    @param[in] darea    area of the cells
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real(kind=dblp), dimension(nlat, nlon), intent(in) :: fracgr
       real(kind=dblp), dimension(nlat),       intent(in) :: darea

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical         :: returnValue
       real(kind=dblp) :: tempor1

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        returnValue = CCPARAM()


        stR(lat,lon)=forshare_st
        sdR(lat,lon)=desshare_st
        snltR(lat,lon)=nlshare_st
        if(sdR(lat,lon).lt.0.) sdR(lat,lon)=0.
        sgR(lat,lon)=1.-stR(lat,lon)-sdR(lat,lon)
        if(sgR(lat,lon).lt.0.) sgR(lat,lon)=0.

        returnValue=CLIMPAR(fracgr,darea)

        returnValue = .true.

        return

      end function ccstatR

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ccdynR
!
!>     @brief This function is a fortran90 port of the initial 77 subroutine
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ccdynR(fracgr,darea) result(returnValue)

      use comemic_mod, only:
      use veget_mod
      use comatm, only: nlat, nlon
      use comrunlabel_mod, only: irunlabelf

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  fracgr   fraction of grass-type plants
!>    @param[in] darea    area of the cells
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real(kind=dblp), dimension(nlat, nlon), intent(in) :: fracgr
       real(kind=dblp), dimension(nlat),       intent(in) :: darea

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical          :: returnValue
! temporal var*/
        integer(kind=ip):: indxv
        real(kind=dblp) :: tempor1,tempor2,db2,fd,dst,dd,nld,dstime, dsd,temp_sg,temp_st

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! calculation of current carbon cycle parameters

        returnValue = ccparam()

! calculation of fraction dynamic variables

        fd=forshare_st-stR(lat,lon)
        dd=desshare_st-sdR(lat,lon)
        nld=nlshare_st-snltR(lat,lon)
        temp_st=stR(lat,lon)
        temp_sg=sgR(lat,lon)


        if (abs(fd).lt.100.*tiny(fd)) then
          fd = 0.0d0
          dst = 0.0d0
        else
          dst=fd*(1.d0-exp(-1./t2t))
        endif

! calculation of forest dynamics; exponential filtre
!-AM    dst=forshare_st-fd*exp(-1./t2t)-st(lat,lon)

!
!
         stR(lat,lon)=stR(lat,lon)+dst
!
        if (stR(lat,lon).lt.0.) stR(lat,lon)=0.
!
        snltR(lat,lon)=nlshare_st-nld*exp(-1./t2t)

! desert dynamics; exponential filtre
        dsd=desshare_st-dd*exp(-1./t2g)-sdR(lat,lon)
        tempor1=sdR(lat,lon)+dsd+stR(lat,lon)

! calculation of characteristic time of desert propagation
        if (tempor1.gt.0.9) then
             dstime=t2g*(1-tempor1)*10.+t2t*(tempor1-0.9)*10.
             dsd=desshare_st-dd*exp(-1./dstime)-sdR(lat,lon)
        endif

        sdR(lat,lon)=sdR(lat,lon)+dsd
        if (sdR(lat,lon).lt.0) sdR(lat,lon)=0.
!
!
        sgR(lat,lon)=1.-stR(lat,lon)-sdR(lat,lon)

        if (sgR(lat,lon).lt.0) sgR(lat,lon)=0.

        returnValue = .true.

        return

      end function ccdynR

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      end module veget_sub_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
