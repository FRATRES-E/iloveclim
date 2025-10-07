!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of [iLOVECLIM/COUPLER]
!!      iLOVECLIM/COUPLER is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
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
!      MODULE: [ipcc_output_mod]
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module [ipcc_output_mod] is handling the creation of a particuliar type of fried noodles ...
!
!>     @date Creation date: October, 01st, 2019
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

      module ipcc_output_mod

       use global_constants_mod, only: dblp=>dp

       implicit none
       public

       public :: IPCC_output

      ! NOTE_avoid_public_variables_if_possible
      ! Following here are private in the end ...

       real(dblp)                       :: tmc_ipcc,ts_ipcc,fr_ipcc,cland_ipcc,coc_ipcc,ctot_ipcc, moc_ipcc,th_ipcc             &
                                        , oh_ipcc, co2_ipcc
       real(dblp)                       :: tmc0,tmc,tsurfmean,moc,thex,cland
       real(dblp), dimension(3)         :: CARSTOK = 0.0_dblp !!! [BUG] <- this variable should be linked to something, not left at zero ...


      contains


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [IPCC_output]
!
!>     @brief This function is dealing with adding IPCC--specific output for the whole iLOVECLIM system
!
!      DESCRIPTION:
!
!>     *** this routine writes the current state of ecbilt to datafiles
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      function IPCC_output(ist,jst) result(returnValue)
!-----------------------------------------------------------------------
!     *** this routine writes the current state of ecbilt to datafiles
!-----------------------------------------------------------------------

      use comdiag,               only: irad
      use comemic_mod,           only: iatm, iyear, nstpyear, nyears
      use comrunlabel_mod,       only: irunlabelf
      use comunit,               only: iuo

!~       use ipcc_output_mod, only: tmc_ipcc,ts_ipcc,fr_ipcc,cland_ipcc,coc_ipcc,ctot_ipcc,moc_ipcc,th_ipcc,tmc0,tmc,    &
!~                                  tsurfmean,moc,thex,cland,oh_ipcc, CARSTOK, co2_ipcc

      use commons_mod,           only: ulrad1nT,fswutoaGA,fswdtoaGA
      use atmos_composition_mod, only: get_PGA_CO2

      implicit none


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  ist Blablabla ...
!>    @param[in]  jst Blablabla ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer, intent(in)    :: ist, jst
       logical                :: returnValue

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer                 :: istep

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      istep=(ist-1)*iatm+jst

      if (mod(istep,nstpyear).eq.1) then
         ts_ipcc=0.
         th_ipcc=0.
         moc_ipcc=0.
         oh_ipcc=0.
         co2_ipcc=get_PGA_CO2()
      endif
      if (mod(istep,iatm).eq.0) then
         ts_ipcc=ts_ipcc+(tsurfmean/360.)
      endif

      if (iyear.ne.0) then
         if (mod(istep,(iatm*30)).eq.0) then
            th_ipcc=th_ipcc+(thex/12.)
            moc_ipcc=moc_ipcc+(moc/12.)
            oh_ipcc=oh_ipcc+((tmc-tmc0)*1.4E+18)
            tmc0=tmc
         endif
      endif

      if (mod(istep,nstpyear).eq.0) then
         if (irad.ne.1) ulrad1nT=0.
         cland_ipcc=Carstok(3)
         coc_ipcc=Carstok(2)
         ctot_ipcc=Carstok(1)+Carstok(2)+Carstok(3)
         write(iuo+48,'(I4,1X,F8.3,1X,F8.3,1X,E15.5,1X,F9.4,1X,F6.2,1X,F8.3,1X,3(E17.7,1X))')                         &
               iyear+irunlabelf,(ulrad1nT-fswutoaGA+fswdtoaGA),ts_ipcc,oh_ipcc,th_ipcc,moc_ipcc,co2_ipcc,ctot_ipcc    &
              ,cland_ipcc,coc_ipcc
         ulrad1nT=0.
         fswutoaGA=0.
         fswdtoaGA=0.

      endif

      if ((mod(istep,nstpyear).eq.0).and.(iyear.eq.nyears)) then
         close (iuo+48)
         close(36)
      endif

      returnValue = .FALSE.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      end function IPCC_output


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module ipcc_output_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
