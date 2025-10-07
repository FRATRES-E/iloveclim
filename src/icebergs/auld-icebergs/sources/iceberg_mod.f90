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
!      MODULE: iceberg_mod
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module iceberg_mod is the FORTRAN 90 transcription of the previous include
!
!>     @date Creation date: June, 19th, 2018
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

      module iceberg_mod

!~        use AnotherModule_mod, only: some_variable
!~        use AnotherModule_mod, only: some_otherfunction

       use global_constants_mod, only: dblp=>dp, ip, dip
       use para0_mod, only: imax, jmax, kmax

       implicit none

       private :: imax, jmax, kmax

!~        public :: someFunction
!~        public :: Foo

      ! NOTE_avoid_public_variables_if_possible




! mab: added flag_calv_coupl, - 28.7.2011
! mab: changed icbmax from 6000000 to 12 000 000 because as soon as GRISLI input
! mab: is used, linit will be 120 per year(10 size classes per month) so if a run
! mab: is for 100 000 years -> 12 000 000
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  fichier "iceberg.com" : incorpore par instruction 'include' dans les routines
!      ICEBST, deficeb, icebtraj, icebdyn
!  modif : 20/04/00

! JONO description of all icebergs vars here would be nice...

! JONO_icb_cpl_wind 	JONO_icb_fix_wind
! [windu, windv] in icebdyn.f will no longer be a function of (fixed) [uwind,vwind]
! instead now a function of ecbiltvar [utot(i,j,3),vtot3] -> [utot10(i,j),vtot10]
! -coupling0-> [sumu10(i,j),sumv10] -ec_co2oc-> [wind10_u(ii,jj),wind10_v] (iceberg var)
!
! wind10_u(,),wind10_v:	Sum of geostrophic and divergent wind at 800hPa (utot(i,j,3))
!			interpolated to 10m with uv10rwv=0.8 (default)	(utot10(i,j))
!			summed over atm. timesteps in coupling0		(sumu10(i,j))
!			translated and rotated in ec_co2oc		(wind10_u(ii,jj))
! JONO APW introduced productionvars xn0(1000)... and common flag_read_positb
!
! JA heat2layers 27-11-2006
! introducing dVol_icb(i,j,k) is volume (m3) change (melt;)
! per icb timestep(=1day) per water layer (ks2=surface, ksubm is icb bottom layer)
! will be plugged straight into thersf.f SO NOT RESISTANT TO TIMESTEP CHANGES

      integer(ip), parameter :: icbmax = 500000, nsrmax = 20000, nmoyenne = 1000, linitmax = 1000

!--variable (en common ?) :
!
! mab: tables with icbmax
!      xn(icbmax), yn(icbmax), hiceb(icbmax), wiceb(icbmax),
!      uiceb(icbmax), viceb(icbmax),pond_icb(icbmax),termBIGx(icbmax),
!      termBIGy(icbmax),kiceb(icbmax),siceb(icbmax), id_abs_table(icbmax)
!      startxn(icbmax),startyn(icbmax),flagstartwritten(icbmax)
!      starthiceb(icbmax),startwiceb(icbmax)


      integer(dip) :: icbwrout

      integer(dip), dimension(linitmax):: pond_icb0
      integer(dip), dimension(icbmax)  :: id_abs_table, pond_icb

      real(dblp)   :: remsim, roiceb, dticeb, repuls
      real(dblp), dimension(imax,jmax) :: urep, vrep
      real(dblp), dimension(icbmax)    :: xn, yn, xn_old, yn_old, hiceb, wiceb, startxn, startyn, starthiceb, startwiceb
      real(dblp), dimension(linitmax)  :: xn0, yn0, hiceb0, wiceb0
      real(dblp), dimension(icbmax)    :: uiceb, viceb
      real(dblp), dimension(nsrmax)    :: t_lach, Debk, berg_mk
      real(dblp), dimension(imax,jmax) :: fonte_icb, vol_icb, vol_icb_cum, fonte_icb_cum, fonte_icb_mois, vol_icb_as, vol_icb_fm  &
                                        , wind10_u, wind10_v, disticeb
      real(dblp), dimension(jmax)      :: fonte_zon, vol_zon
      real(dblp), dimension(imax,jmax,kmax):: dVol_icb = 0.0_dblp, dVol_calv = 0.0_dblp
      real(dblp), dimension(icbmax)    :: vol_orig, vol_melt, termBIGx,termBIGy
      real(dblp), dimension(20)        :: poids
      real(dblp), dimension(nsrmax)    :: srcup, srcdown, srcx, srcy, srcw, srch, pond
      real(dblp)                       :: d_cotes

      integer(dip)                     :: nbtot_write_iceb, nbtot_write_max
! mab: for opening new track and traj files so that they dont get too long
      integer(ip)                      :: filenumber, numtracki
      integer(ip), dimension(icbmax)   :: flagstartwritten
! mab: to count the number of times there would have been a warning in
!      fort.2058
      integer(dip)                     :: wind10u, wind10v, seaicevel

      integer(ip), dimension(nsrmax)   :: nsr_tmp, nsr_frq
      integer(ip)                      :: nfterm, nsource, nfintab, nftemp1, nftemp2, nftemp3, nftemp4, nbcase, nfmoy, nbrjour

      integer(ip), dimension(icbmax)   :: kiceb, siceb
      integer(ip)                      :: lmx, nbricb, nbricb_plus, nbricb_moins, nbricb_plus_n, nbricb_plus_s, nbricb_moins_n    &
                                        , nbricb_moins_s, linit, nfricb, nficeb, nbrmois, mois, necriture, numiticb, nitrunicb    &
                                        ,it_icb, it_icbl, nfrqicb, nbrmoins_cum_s, nbrplus_cum_s, nbrmoins_cum_n       &
                                        ,nbrplus_cum_n, nbricb_n, nbricb_s, nbrmoins_cum, nbrplus_cum

! mab:ndbug = iterations between debug writing etc
      character(len=30)                      :: fmticb, fmttrm
      character(len=30), dimension(nsrmax)   :: srcetat
      character(len=20)                      :: icbexp

! JONO_out march 2004 ndbug declared as common int4
      integer(dip)                           :: siceb_1, siceb_2, siceb_3, siceb_4, siceb_5, siceb_6, siceb_7, siceb_8, siceb_9   &
                                               ,siceb_10, ndbug, necriflag
      integer(dip)                           :: id_max_used
      logical                                :: flag_write_iceb, flag_read_positb, flag_calv_coupl, flag_output_icb


!~       contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [Routine or Function Name Here]
!
!>     @brief This subroutine / function is adding a wonderful knew unused functionality
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~       function someFunction(inParam, outParam, inOutParam) result(returnValue)

!~        use AnotherModule_mod, only: some_variable
!~        use AnotherModule_mod, only: some_otherfunction

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~        real, intent(in)    :: inParam
!~        real, intent(inout) :: inOutParam
!~        real, intent(out)   :: outParam

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~        real                :: returnValue
!~        real                :: someVariable


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!>    @bug Description of the stupid sticky bug that we know exist there but is not corrected yet!

!~       end function someFunction

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module iceberg_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
