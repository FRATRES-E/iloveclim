!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of the downscaling procedure of the ECBilt-zoom function
!!      atmphys_d is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
!!      as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!!
!!      atmphys_d is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
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
!      MODULE: atmphys_d
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module atmphys-d contains the downscaling-enabled routines derived from atmphys0.f
!
!>     @date Creation date: January, 06th, 2016
!>     @date Last modification: SLastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module atmphys_d

#if ( DOWNSTS == 1 )

       implicit none
       private

       public :: ec_ptmoisgp_d
       public :: ec_ptmoisgp_dreduce

      ! NOTE_avoid_public_variables_if_possible

      contains


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ec_ptmoisgp_d
!
!>     @brief This function is a downscaling enabled port of the original ec_ptmoisgp function, with explicit topography reference
!
!      DESCRIPTION:
!
!>     Computation of ground pressure and temperature in order to calculate the maximum precipitable water content in latlon i,j
!>     qmount contains the topography for this purpose assuming temperature varies linearly with log of the pressure
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ec_ptmoisgp_d(temp2g_ij,temp4g_ij,geopg_ij2,tmount,qmax,dqmdt,qmount_ipoints) result(returnValue)

#if ( COMATM == 1 )
        use comatm,  only: nlat, nlon, nsh, nsh2, nvl, grav, nm, ntl, rgas, rlogtl12, tlevel, plevel
        use comdyn,  only: geopg
        use comphys, only: temp4g, hmoisr, gpm500, qmount
        use comunit
#endif

        use newunit_mod, only: error_id
        
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  temp2g_ij temp2g variable at the i,j location
!>    @param[in]  temp4g_ij temp4g variable at the i,j location
!>    @param[in]  geopg_ij2 geopg  variable at the i,j location and level 2
!>    @param[out] tmount
!>    @param[out] qmax
!>    @param[out] dqmdt
!>    @param[out] qmount_ipoints
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       double precision              , intent(in)    :: temp2g_ij, temp4g_ij, geopg_ij2
       double precision, dimension(:), intent(in)    :: qmount_ipoints
       double precision, dimension(:), intent(out)   :: tmount, qmax, dqmdt

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Through commons variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( COMATM == 0 )
#include "comatm.h"
! dmr comdyn.h provides:
! dmr          geopg
#include "comdyn.h"
! dmr comphys.h provides:
! dmr         temp4g
! dmr         hmoisr
! dmr         gpm500
! dmr         qmount
#include "comphys.h"
! dmr comunit.h provides:
! dmr         iuo
#include "comunit.h"
#endif

! dmr comrunlabel.h is not used??
! #include "comrunlabel.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical             :: returnValue
       integer             :: k

! dmr  z500        : is the mean altitude of the 500 hPa level (decameters?)
! dmr  hfac        : is ... do not know!
! dmr  hred        : is ... do not know!
! dmr  pfac        : is ... do not know!
! dmr  alpha       : is the local slope of temperature decrease in log(P)
! dmr  t500        : is the mean temperature at the 500 hPa level
! dmr  hmount      : is the local topography (meters)
! dmr  delta_qmount: is the step in the intermediate points for vertical sampling

       real(kind=8)                        :: z500, hfac, hred, pfac, alpha, t500, hmount
       real(kind=8)                        :: delta_qmount

       double precision, parameter         :: magic_coeff = 0.6d0 ! used to tune the amount of subgrid precipitation ...

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   1.0: compute constants independent of surface topography
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr  grav  :     is gravity acc. defined in comatm and initial0
! dmr  gpm500:     is from comphys.h and initial0: mean 500 hPa geopotential height [m]
       z500 = gpm500*grav
! dmr  rgas  :     Gas Constant, defined as 287. in initial0
       hfac = 2.0_8/rgas

! dmr  hmoisr:     reduction factor of mountain heights in order to tune the amount of water that is allowed to pass a topographic
!                  barier, value is currently 1.0d0 in initial0 (!)
       hred = hmoisr*grav*magic_coeff

! dmr  pressure and temperature at second level
       pfac = log(plevel(2)/tlevel(2))

! dmr calculate temperature at t500 assuming the temperature varies
! dmr linearly with log(p) : T = Tref + alpha * log (p/pref)

       alpha = (temp2g_ij - temp4g_ij)*rlogtl12
       t500  = temp4g_ij + alpha*pfac

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   2.0: starts part with topography influence
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr   2.2: loop over the points to populate the output

      do k=1, UBOUND(qmount_ipoints(:),DIM=1)

! dmr calculate reduced ground height in decameters
! dmr reduction occurs in order to tune the amount of moisture which
! dmr is allowed to pass a topographic barier

!~         hmount=qmount_ipoints(k)*hred
!~         if (hmount.lt.0d0) hmount=0d0
        hmount = max(qmount_ipoints(k)*hred,0.d0)

! dmr calculate the groundpressure assuming that the mean geopotential
! dmr height at 500 hPa is gpm500 decameter
! dmr calculate 10 mtr temperature in K


        tmount(k)=t500**2 - hfac*alpha*(hmount-geopg_ij2-z500)
! - dmr removing safeguarding

!~         if (tmount(k).lt.0) then
!~           write(error_id,*) 'in latlon ', "unknown, this is _d version ..."
!~           write(error_id,*) tmount(k),hmount,t500,geopg_ij2
!~           call ec_error(18)
!~         else
          tmount(k)=sqrt(tmount(k))
!~         endif

!        qmax=ec_detqmax_d(tmount,i,j,dqmdt)
        qmax(k)=ec_detqmax_d(tmount(k),dqmdt(k),temp2g_ij,temp4g_ij,geopg_ij2,hmount,alpha,t500,z500)

      enddo

      returnValue = .true.

      end function ec_ptmoisgp_d

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ec_ptmoisgp_dreduce
!
!>     @brief This function is a downscaling enabled port of the original ec_ptmoisgp function, with explicit topography reference
!
!      DESCRIPTION:
!
!>     Computation of ground pressure and temperature in order to calculate the maximum precipitable water content in latlon i,j
!>     qmount contains the topography for this purpose assuming temperature varies linearly with log of the pressure
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ec_ptmoisgp_dreduce(temp2g_ij,temp4g_ij,geopg_ij2,tmount,qmax,dqmdt,qmount_ipoints,k_min,k_max) result(returnValue)

#if ( COMATM == 1 )
        use comatm, only: nlat, nlon, nsh, nsh2, nvl, grav, nm, ntl, rgas, rlogtl12, tlevel, plevel
        use comdyn, only: geopg
        use comphys,only: gpm500, hmoisr
        use comunit,only: iuo
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  temp2g_ij temp2g variable at the i,j location
!>    @param[in]  temp4g_ij temp4g variable at the i,j location
!>    @param[in]  geopg_ij2 geopg  variable at the i,j location and level 2
!>    @param[out] tmount
!>    @param[out] qmax
!>    @param[out] dqmdt
!>    @param[out] qmount_ipoints
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       double precision              , intent(in)    :: temp2g_ij, temp4g_ij, geopg_ij2
       double precision, dimension(:), intent(in)    :: qmount_ipoints
       integer                       , intent(in)    :: k_min, k_max
       double precision, dimension(:), intent(out)   :: tmount, qmax, dqmdt

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Through commons variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( COMATM == 0 )
#include "comatm.h"
! dmr comdyn.h provides:
! dmr          geopg
#include "comdyn.h"
! dmr comphys.h provides:
! dmr         temp4g
! dmr         hmoisr
! dmr         gpm500
! dmr         qmount
#include "comphys.h"
! dmr comunit.h provides:
! dmr         iuo
#include "comunit.h"
#endif


! dmr comrunlabel.h is not used??
! #include "comrunlabel.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical             :: returnValue
       integer             :: k

! dmr  z500        : is the mean altitude of the 500 hPa level (decameters?)
! dmr  hfac        : is ... do not know!
! dmr  hred        : is ... do not know!
! dmr  pfac        : is ... do not know!
! dmr  alpha       : is the local slope of temperature decrease in log(P)
! dmr  t500        : is the mean temperature at the 500 hPa level
! dmr  hmount      : is the local topography (meters)
! dmr  delta_qmount: is the step in the intermediate points for vertical sampling

       real(kind=8)                        :: z500, hfac, hred, pfac, alpha, t500, hmount
       real(kind=8)                        :: delta_qmount

       double precision, parameter         :: magic_coeff = 0.6d0 ! used to tune the amount of subgrid precipitation ...

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   1.0: compute constants independent of surface topography
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr  grav  :     is gravity acc. defined in comatm and initial0
! dmr  gpm500:     is from comphys.h and initial0: mean 500 hPa geopotential height [m]
       z500 = gpm500*grav
! dmr  rgas  :     Gas Constant, defined as 287. in initial0
       hfac = 2.0_8/rgas

! dmr  hmoisr:     reduction factor of mountain heights in order to tune the amount of water that is allowed to pass a topographic
!                  barier, value is currently 1.0d0 in initial0 (!)
       hred = hmoisr*grav*magic_coeff

! dmr  pressure and temperature at second level
       pfac = log(plevel(2)/tlevel(2))

! dmr calculate temperature at t500 assuming the temperature varies
! dmr linearly with log(p) : T = Tref + alpha * log (p/pref)

       alpha = (temp2g_ij - temp4g_ij)*rlogtl12
       t500  = temp4g_ij + alpha*pfac

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   2.0: starts part with topography influence
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr   2.2: loop over the points to populate the output

!~       do k=1, UBOUND(qmount_ipoints(:),DIM=1)
      do k=k_min, k_max

! dmr calculate reduced ground height in decameters
! dmr reduction occurs in order to tune the amount of moisture which
! dmr is allowed to pass a topographic barier

!~         hmount=qmount_ipoints(k)*hred
!~         if (hmount.lt.0d0) hmount=0d0
        hmount = max(qmount_ipoints(k)*hred,0.d0)

! dmr calculate the groundpressure assuming that the mean geopotential
! dmr height at 500 hPa is gpm500 decameter
! dmr calculate 10 mtr temperature in K


        tmount(k)=t500**2 - hfac*alpha*(hmount-geopg_ij2-z500)

! - dmr removing safeguarding
!~         if (tmount(k).lt.0) then
!~           write(error_id,*) 'in latlon ', "unknown, this is _d version ..."
!~           write(error_id,*) tmount(k),hmount,t500,geopg_ij2
!~           call ec_error(18)
!~         else
          tmount(k)=sqrt(tmount(k))
!~         endif

!        qmax=ec_detqmax_d(tmount,i,j,dqmdt)
        qmax(k)=ec_detqmax_d(tmount(k),dqmdt(k),temp2g_ij,temp4g_ij,geopg_ij2,hmount,alpha,t500,z500)

      enddo

      returnValue = .true.

      end function ec_ptmoisgp_dreduce

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: ec_detqmax_d
!
!>     @brief This function is a downscaling enabled port of the original ec_detqmax_d function, with topo added
!
!      DESCRIPTION:
!
!>     determines the maximum water content in latlon point i,j for given ground- and 650 and 350 hPa temperature by linear
!>     interpolation in qmtabel
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function ec_detqmax_d(tmount,dqmdt,temp2g_ij,temp4g_ij,geopg_ij2_l,hmount_l,alpha_l,t500_l,z500_l) result(returnValue)

#if ( COMATM == 1 )
      use comatm, only: alogpl2tl2, nlat, nlon, nsh2, nvl, alogtl12, alogtl1pl2, grav, nsh, nm, ntl, rgas, rlogtl12
      use comphys,only: iqmtab, jqmtab, kqmtab, rdtqmi, rdtqmj, rdtqmk, tqmimin, tqmjmin, tqmkmin, tqmi, tqmj, tqmk, qmtabel
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real(kind=8), intent(in)    :: tmount, temp4g_ij, temp2g_ij, hmount_l, alpha_l, t500_l, z500_l, geopg_ij2_l
       real(kind=8), intent(out)   :: dqmdt

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Through commons variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#if ( COMATM == 0 )
! dmr comphys.h provides:
! dmr         tqmi
! dmr         tqmj
! dmr         tqmk
! dmr         iqmtab
! dmr         jqmtab
! dmr         kqmtab
! dmr         tqmimin
! dmr         tqmjmin
! dmr         tqmkmin
! dmr         qmtabel
! dmr         rdtqmi
! dmr         rdtqmj
! dmr         rdtqmk

#include "comphys.h"

#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real(kind=8)        :: returnValue
       real(kind=8)        :: ti,tj,tk,dqmdi,dqmdj,dqmdk,qmax
       real(kind=8)        :: dtgdt
       real(kind=8)        :: mountred ! afq precip threshold on qmax depends on local elevation

       integer             :: ii,jj,kk

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       ti=temp4g_ij
       tj=tmount-temp4g_ij
       tk=temp4g_ij-temp2g_ij


!~        if (ti.lt.tqmi(0)) then
!~          ti=tqmi(0)
!~ !         write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!~ !         call ec_error(121)
!~        endif
!~        if (ti.gt.tqmi(iqmtab)) then
!~          ti=tqmi(iqmtab)
!~ !         write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!~ !         call ec_error(121)
!~        endif

       ti = min(max(ti,tqmi(0)),tqmi(iqmtab))

!~        if (tj.lt.tqmj(0)) then
!~          tj=tqmj(0)
!~ !         write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!~ !         call ec_error(121)
!~        endif
!~        if (tj.gt.tqmj(jqmtab)) then
!~          tj=tqmj(jqmtab)
!~ !         write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!~ !         call ec_error(121)
!~        endif

       tj = min(max(tj,tqmj(0)),tqmj(jqmtab))

!~        if (tk.lt.tqmk(0)) then
!~          tk=tqmk(0)
!~ !         write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!~ !         call ec_error(121)
!~        endif
!~        if (tk.gt.tqmk(kqmtab)) then
!~          tk=tqmk(kqmtab)
!~ !         write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!~ !         call ec_error(121)
!~        endif

       tk = min(max(tk,tqmk(0)),tqmk(kqmtab))


       ii=min(iqmtab-1,int((ti-tqmimin)*rdtqmi))
       jj=min(jqmtab-1,int((tj-tqmjmin)*rdtqmj))
       kk=min(kqmtab-1,int((tk-tqmkmin)*rdtqmk))

       dqmdi=(qmtabel(ii+1,jj,kk)-qmtabel(ii,jj,kk))*rdtqmi
       dqmdj=(qmtabel(ii,jj+1,kk)-qmtabel(ii,jj,kk))*rdtqmj
       dqmdk=(qmtabel(ii,jj,kk+1)-qmtabel(ii,jj,kk))*rdtqmk

       qmax = qmtabel(ii,jj,kk) + (ti-tqmi(ii))*dqmdi + (tj-tqmj(jj))*dqmdj + (tk-tqmk(kk))*dqmdk

       qmax = min(max(qmax,0d0),0.2d0)

!~        if (qmax.lt.0d0) then
!~          qmax=0d0
!~        elseif (qmax.gt.0.2) then
!~ !dmr --- Tentative fix to try and stop the crash of model
!~          qmax = 0.2
!~ !dmr --- Tentative fix to try and stop the crash of model
!~ !        write(error_id,*) 'in latlon ',i,j,' qmax ',qmax, 'tentative fix'
!~ !        call ec_error(121)
!~        endif

!       alpha=(temp2g_ij-temp4g_ij)*rlogtl12
!       t500=temp4g_ij+alpha*alogpl2tl2
!       z500=gpm500*grav
!       hmount=qmount(i,j)*hmoisr*grav

       dtgdt=(rgas*t500_l*alogtl1pl2 + (hmount_l-geopg_ij2_l-z500_l))/(rgas*tmount*alogtl12)

       dqmdt=dqmdi + dqmdj * (dtgdt - 1d0) + dqmdk

       ! afq -- hmount_l is in decameters (!)
       ! here we assume a linear relation with altitude from 0 to 2500 m
       !mountred = min (0.8 + (hmount_l / (10.*2000.)) * (1. - 0.8), 1.0)
       !mountred = min (0.8 + 0.2 * hmount_l / 20000., 1.0)
       mountred = min (0.75 + 0.25 * hmount_l / 10000., 1.0)

       returnValue = mountred*1.00*qmax
!      ec_detqmax=0.85*qmax

      end function ec_detqmax_d

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#endif
end module atmphys_d

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
