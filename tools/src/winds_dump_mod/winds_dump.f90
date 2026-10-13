! tools/src/winds_dump_mod/winds_dump.f90 --- A small module to dump winds in binary file for replay
      module wind_dump_mod
        use global_constants_mod, only: dblp=>dp
        implicit none
        private
        public :: wind_dump
        integer, save :: udump = -1
      contains
        subroutine wind_dump(istep)
          use comatm,  only: nlat, nlon, nsh2, dtime, om, radius, iwater
          use comdyn,  only: psi, chi, u800, v800, udivg, vdivg, divg
          use comphys, only: rmoisg, umoisr
          integer, intent(in) :: istep
          if (udump < 0) then
            open(newunit=udump, file='wind_dump.bin', form='unformatted', access='stream', status='replace')
            write(udump) int(nlat,4), int(nlon,4), int(nsh2,4), real(dtime,8), real(om,8), real(radius,8), &
                         real(umoisr,8)
          endif
          write(udump) int(istep,4), real(psi(:,3),8), real(chi(:,3),8),                   &
                       real(u800 + udivg(:,:,3),8), real(v800 + vdivg(:,:,3),8),           &
                       real(divg(:,:,3),8), real(rmoisg(:,:,iwater),8)
          flush(udump)
        end subroutine wind_dump
      end module wind_dump_mod
