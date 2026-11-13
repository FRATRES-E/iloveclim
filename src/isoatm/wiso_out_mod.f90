!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2026 Didier M. Roche (a.k.a. dmr) | iLOVECLIM / FRATRES coding group

!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at

!       http://www.apache.org/licenses/LICENSE-2.0

!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!   Style sheet: v1.0.0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [wiso_out_mod]
!
!>     @brief Water isotope outputs of ECBilt in the standard atmospheric netCDF files (ROADMAP 2e, wiso:2e:0003).
!
!      DESCRIPTION:
!
!>     The variables are declared in outp_wisoatm.param (installed by the isoatm component, read after
!!     outp_atmos.param, numbered after the standard variables; same format and output flags) and found here by their
!!     short name; a variable absent from that file is not written.
!!
!!     Quantities, for the water and for 17O, 18O and 2H (iwat17, iwat18, iwat2h):
!!       p  total precipitation (torain + tosnow)    r  rain (torain)    s  snow (tosnow)
!!       e  evaporation (evap)                       v  vapour column (rmoisg)
!!     summed at every time step, at the same point of the step as the other ECBilt outputs (ec_selectout), over the
!!     interval of each file: yearly-monthly means (atmmmyl), whole-run monthly means (atmmmwp, per calendar month),
!!     yearly means (atmym); instantaneous outputs (atminst) use the values of the step.
!!
!!     Variables (short names; X = p, r, s, e, v):
!!       amounts   pp18 pp17 ppd (p), sn18 sn17 snd (s), ev18 ev17 evd (e) in cm.yr-1; q18 q17 qd (v) in the units of q:
!!                 means over the interval, consistent with pp, snow, evap and q of outp_atmos.param (wiso:2e:0001), so
!!                 that ratios formed offline are amount-weighted;
!!       deltas    d18X, d17X, ddX (permil), dxX = ddX - 8 d18X (permil), e17X = ln(R17) - 0.528 ln(R18) (per meg),
!!                 from the sums: R = sum(isotopologue)/sum(water), delta = R - 1 (rsmow = 1), never an average of deltas;
!!                 set to fill_value where the water sum is below WO_TINY of its field maximum (no amount to weigh).
!!
!!     Contract
!!       wiso_out_step()              : ec_selectout, every step, before the output routines.
!!       wiso_out_write(file, write)  : ec_outputinst / ec_outputmtl / ec_outputmyl / ec_outputyrl, with their file
!!                                      open, before their close; writes when write is true.
!
!      REFERENCES: ROADMAP 2e, D40 (permil and per meg at the output boundary only).
!
!      AUTHORS: dmr&clo
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module wiso_out_mod

#if ( WISOATM == 1 )

        use global_constants_mod, only: ip, dblp=>dp
        use comatm,               only: nlat, nlon, iwater, iwat17, iwat18, iwat2h

        implicit none

        private

        public :: WO_INST, WO_MTL, WO_MYL, WO_YRL
        public :: wiso_out_step, wiso_out_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Constants and state
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip), parameter :: WO_INST = 1, WO_MTL = 2, WO_MYL = 3, WO_YRL = 4    !< output files

        integer(ip), parameter :: QP = 1, QR = 2, QS = 3, QE = 4, QV = 5, NQ = 5       !< quantities
        integer(ip), parameter :: SW = 1, S17 = 2, S18 = 3, S2H = 4, NSP = 4          !< water, 17O, 18O, 2H
        character(len=1), parameter :: QCHAR(NQ) = ['p', 'r', 's', 'e', 'v']

        real(dblp), parameter :: PFAC    = 100.0_dblp*3600.0_dblp*24.0_dblp*360.0_dblp  !< m/s -> cm/yr (as ec_selectout)
        real(dblp), parameter :: WO_TINY = 1.0e-10_dblp     !< water sum below this fraction of its maximum: no delta
        real(dblp), parameter :: THETA17 = 0.528_dblp       !< reference slope of the 17O-excess

        ! output variables: short name, quantity, kind (0 amount, 1 d18, 2 d17, 3 dD, 4 d-excess, 5 17O-excess), species
        integer(ip), parameter :: NAMOUNT = 12, NDELTA = 5*NQ, NVAR = NAMOUNT + NDELTA
        character(len=8), save :: vname(NVAR)
        integer(ip),      save :: vq(NVAR), vkind(NVAR), vsp(NVAR), vidx(NVAR)
        logical,          save :: ready = .false.

        real(dblp), allocatable, save :: smyl(:,:,:,:), syrl(:,:,:,:), smtl(:,:,:,:,:), sinst(:,:,:,:)
        real(dblp),              save :: nmyl = 0.0_dblp, nyrl = 0.0_dblp, nmtl(12) = 0.0_dblp
        logical,                 save :: myl_reset = .true.

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Every step: sums of the interval of each file
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wiso_out_step()

          use comemic_mod, only: imonth, iday
          use comdiag,     only: diag_stepinday

          if (.not. ready) call setup()
          call current(sinst)

          if (myl_reset) then
            smyl = 0.0_dblp
            nmyl = 0.0_dblp
            myl_reset = .false.
          endif
          smyl = smyl + sinst
          nmyl = nmyl + 1.0_dblp

          smtl(:,:,:,:,imonth) = smtl(:,:,:,:,imonth) + sinst
          nmtl(imonth) = nmtl(imonth) + 1.0_dblp

          if (imonth == 1 .and. iday == 1 .and. diag_stepinday == 1) then
            syrl = 0.0_dblp
            nyrl = 0.0_dblp
          endif
          syrl = syrl + sinst
          nyrl = nyrl + 1.0_dblp

        end subroutine wiso_out_step

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Writing, with the file of the calling output routine open
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine wiso_out_write(file, do_write)

          use comemic_mod, only: imonth

          integer(ip), intent(in) :: file
          logical,     intent(in) :: do_write

          if (.not. do_write) return
          if (.not. ready) call setup()
          select case (file)
          case (WO_INST)
            call current(sinst)
            call write_all(sinst, 1.0_dblp)
          case (WO_MYL)
            call write_all(smyl, nmyl)
            myl_reset = .true.
          case (WO_MTL)
            call write_all(smtl(:,:,:,:,imonth), nmtl(imonth))
          case (WO_YRL)
            call write_all(syrl, nyrl)
          end select

        end subroutine wiso_out_write

        subroutine write_all(sm, n)

          use atmoutp_mod, only: out_write => write
          use comdiag,     only: fill_value

          real(dblp), intent(in) :: sm(:,:,:,:)      !< (nlat,nlon,NQ,NSP) sums
          real(dblp), intent(in) :: n                !< number of steps summed

          real(dblp)  :: f(nlat, nlon, 1), r17(nlat, nlon), r18(nlat, nlon), r2h(nlat, nlon), w(nlat, nlon)
          logical     :: m(nlat, nlon)
          integer(ip) :: v, q

          if (n <= 0.0_dblp) return
          do v = 1, NVAR
            if (vidx(v) == 0) cycle
            q = vq(v)
            if (vkind(v) == 0) then
              f(:,:,1) = sm(:,:,q,vsp(v))/n
              if (q /= QV) f(:,:,1) = f(:,:,1)*PFAC
            else
              w = sm(:,:,q,SW)
              m = w > WO_TINY*maxval(w) .and. w > 0.0_dblp
              r17 = 1.0_dblp
              r18 = 1.0_dblp
              r2h = 1.0_dblp
              where (m)
                r17 = sm(:,:,q,S17)/w
                r18 = sm(:,:,q,S18)/w
                r2h = sm(:,:,q,S2H)/w
              end where
              select case (vkind(v))
              case (1)
                f(:,:,1) = (r18 - 1.0_dblp)*1.0e3_dblp
              case (2)
                f(:,:,1) = (r17 - 1.0_dblp)*1.0e3_dblp
              case (3)
                f(:,:,1) = (r2h - 1.0_dblp)*1.0e3_dblp
              case (4)
                f(:,:,1) = ((r2h - 1.0_dblp) - 8.0_dblp*(r18 - 1.0_dblp))*1.0e3_dblp
              case (5)
                f(:,:,1) = (log(r17) - THETA17*log(r18))*1.0e6_dblp
              end select
              where (.not. m) f(:,:,1) = fill_value
            endif
            call out_write(vidx(v), f)
          enddo

        end subroutine write_all

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Helpers
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        ! values of the step (m/s for the fluxes, column amount for the vapour)
        subroutine current(c)

          use comphys,     only: rmoisg, torain, tosnow
          use comsurf_mod, only: evap

          real(dblp), intent(out) :: c(:,:,:,:)

          integer(ip) :: s, k

          do s = 1, NSP
            k = species(s)
            c(:,:,QP,s) = torain(:,:,k) + tosnow(:,:,k)
            c(:,:,QR,s) = torain(:,:,k)
            c(:,:,QS,s) = tosnow(:,:,k)
            c(:,:,QE,s) = evap(:,:,k)
            c(:,:,QV,s) = rmoisg(:,:,k)
          enddo

        end subroutine current

        pure function species(s) result(k)
          integer(ip), intent(in) :: s
          integer(ip)             :: k
          select case (s)
          case (S17)
            k = iwat17
          case (S18)
            k = iwat18
          case (S2H)
            k = iwat2h
          case default
            k = iwater
          end select
        end function species

        ! variable table, and the output numbers of those declared in outp_wisoatm.param
        subroutine setup()

          use comdiag, only: nametotvar, numtotvar

          character(len=3), parameter :: ANAME(3) = ['pp ', 'sn ', 'ev ']
          integer(ip),      parameter :: AQ(3)    = [QP, QS, QE]
          character(len=3), parameter :: DNAME(5) = ['d18', 'd17', 'dd ', 'dx ', 'e17']
          character(len=2), parameter :: SNAME(2:4) = ['17', '18', 'd ']
          integer(ip) :: v, a, s, q, d, i

          v = 0
          do a = 1, 3                                    ! pp18 pp17 ppd, sn.., ev..
            do s = S17, S2H
              v = v + 1
              vname(v) = trim(ANAME(a))//trim(SNAME(s))
              vq(v) = AQ(a)
              vkind(v) = 0
              vsp(v) = s
            enddo
          enddo
          do s = S17, S2H                                ! q18 q17 qd
            v = v + 1
            vname(v) = 'q'//trim(SNAME(s))
            vq(v) = QV
            vkind(v) = 0
            vsp(v) = s
          enddo
          do q = 1, NQ                                   ! d18p d17p ddp dxp e17p, ... r, s, e, v
            do d = 1, 5
              v = v + 1
              vname(v) = trim(DNAME(d))//QCHAR(q)
              vq(v) = q
              vkind(v) = d
              vsp(v) = SW
            enddo
          enddo

          vidx = 0
          do v = 1, NVAR
            do i = 1, numtotvar
              if (trim(nametotvar(i,2)) == trim(vname(v))) vidx(v) = i
            enddo
          enddo

          allocate(smyl(nlat, nlon, NQ, NSP), syrl(nlat, nlon, NQ, NSP), smtl(nlat, nlon, NQ, NSP, 12))
          allocate(sinst(nlat, nlon, NQ, NSP))
          smyl = 0.0_dblp
          syrl = 0.0_dblp
          smtl = 0.0_dblp
          ready = .true.

        end subroutine setup

#endif

      end module wiso_out_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
