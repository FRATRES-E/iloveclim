!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright 2026 iLOVECLIM / LUDUS coding group

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

#include "choixcomposantes.h"

#ifndef EVOLU_NETCDF
#define EVOLU_NETCDF 0
#endif

#if ( EVOLU_NETCDF == 1 )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [evolu_netcdf_mod]
!>     @brief Self-describing NetCDF output for the CLIO chronological diagnostics ("evolu").
!      DESCRIPTION:
!>     Companion to evolu_diag_mod. Same nvinfo global diagnostics, same cadence, same already-averaged values as the
!!     text "evolu", into a self-describing NetCDF-4 file. Global integrals -> no spatial axis: one unlimited 'time'
!!     dimension and nvinfo scalar var(time). CF-like metadata come from a tabular namelist mapping each INTERNAL
!!     registry name (titvar) to a CF-valid NetCDF name + attributes; level families expand from a '%' template.
!>     @note dmr&clo -- Direct netcdf binding (FROG pattern), not io_nc_mod. Coexists with the text evolu (guarded).
!>     @note dmr&clo -- DEV MODE: partial namelist allowed; uncovered diagnostics are listed and skipped unless
!!           strict_coverage=.true. Time: 'months since 1583-01-01', calendar '360_day', +1 day/record (reserve R10).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module evolu_netcdf_mod

       use global_constants_mod, only: dblp=>dp, ip
       use netcdf

       implicit none
       private
       public :: evolu_nc_init, evolu_nc_write, evolu_nc_close

       integer(ip), parameter :: NFIELD = 5
       integer(ip), parameter :: LEN_NAME = 64, LEN_LONG = 128, LEN_UNIT = 32

       type :: evolu_nc_var
         character(len=LEN_NAME) :: nc_name   = ''
         character(len=LEN_LONG) :: long_name = ''
         character(len=LEN_NAME) :: std_name  = ''
         character(len=LEN_UNIT) :: units     = ''
         integer(ip)             :: reg_index = 0
         integer(ip)             :: varid     = -1
         real(dblp)              :: depth     = -1._dblp
       end type evolu_nc_var

       type :: nml_row
         character(len=LEN_NAME) :: internal  = ''
         character(len=LEN_NAME) :: nc_name   = ''
         character(len=LEN_LONG) :: long_name = ''
         character(len=LEN_NAME) :: std_name  = ''
         character(len=LEN_UNIT) :: units     = ''
       end type nml_row

       integer(ip), save :: ncid       = -1
       integer(ip), save :: time_dimid = -1, time_varid = -1
       integer(ip), save :: nvar_nc    = 0
       integer(ip), save :: rec_count  = 0
       real(dblp),  save :: day_count  = 0._dblp
       type(evolu_nc_var), allocatable, save :: ncvars(:)

      contains

       subroutine evolu_nc_init(filename, nmlpath, nvinfo, titvar, level_depths, ks1, ks2, resolve_index, strict_coverage)
         use newunit_clio_mod, only: clio3_out_id
         character(len=*), intent(in) :: filename, nmlpath
         integer(ip),      intent(in) :: nvinfo, ks1, ks2
         character(len=*), intent(in) :: titvar(nvinfo)
         real(dblp),       intent(in) :: level_depths(ks1:ks2)
         logical,          intent(in) :: strict_coverage
         interface
           integer(ip) function resolve_index(name)
             import :: ip
             character(len=*), intent(in) :: name
           end function
         end interface
         type(nml_row), allocatable :: rows(:), exp_rows(:)
         integer(ip),   allocatable :: exp_lev(:)
         type(nml_row)              :: one_exp(ks2-ks1+1)
         integer(ip)                :: one_lev(ks2-ks1+1)
         integer(ip) :: nrows, i, e, nexp, nout, ridx, iv, nlev
         logical, allocatable :: covered(:)
         character(len=LEN_NAME) :: intname

         nlev = ks2 - ks1 + 1
         call read_namelist(nmlpath, rows, nrows, clio3_out_id)

         allocate(exp_rows(nvinfo + nlev), exp_lev(nvinfo + nlev))
         nexp = 0
         do i = 1, nrows
           call expand_row(rows(i), ks1, ks2, one_exp, nout, one_lev)
           do e = 1, nout
             nexp = nexp + 1
             if (nexp > size(exp_rows)) then
               write(clio3_out_id,*) 'STOP evolu_netcdf: expanded rows exceed bound'
               error stop 'evolu_netcdf: expansion overflow'
             endif
             exp_rows(nexp) = one_exp(e) ; exp_lev(nexp) = one_lev(e)
           enddo
         enddo

         allocate(ncvars(nexp), covered(nvinfo)) ; covered = .false. ; nvar_nc = 0
         do i = 1, nexp
           intname = exp_rows(i)%internal
           ridx    = resolve_index(trim(intname))
           if (ridx < 1 .or. ridx > nvinfo) then
             write(clio3_out_id,*) 'STOP evolu_netcdf: resolved index out of range for ', trim(intname)
             error stop 'evolu_netcdf: bad resolved index'
           endif
           if (covered(ridx)) then
             write(clio3_out_id,*) 'STOP evolu_netcdf: diagnostic covered twice: ', trim(intname)
             error stop 'evolu_netcdf: duplicate coverage'
           endif
           covered(ridx) = .true. ; nvar_nc = nvar_nc + 1
           ncvars(nvar_nc)%nc_name   = exp_rows(i)%nc_name
           ncvars(nvar_nc)%long_name = exp_rows(i)%long_name
           ncvars(nvar_nc)%std_name  = exp_rows(i)%std_name
           ncvars(nvar_nc)%units     = exp_rows(i)%units
           ncvars(nvar_nc)%reg_index = ridx
           if (exp_lev(i) > 0) ncvars(nvar_nc)%depth = level_depths(ks1 + exp_lev(i) - 1)
         enddo

         do i = 1, nvinfo
           if (.not. covered(i)) then
             if (strict_coverage) then
               write(clio3_out_id,*) 'STOP evolu_netcdf: diagnostic not covered: ', trim(titvar(i))
               error stop 'evolu_netcdf: incomplete coverage'
             else
               write(clio3_out_id,*) 'evolu_netcdf (dev): NOT in namelist, skipped: ', trim(titvar(i))
             endif
           endif
         enddo

         call nc_check( nf90_create(trim(filename), NF90_NETCDF4, ncid), 'create '//trim(filename) )
         call nc_check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dimid), 'def_dim time' )
         call nc_check( nf90_def_var(ncid, 'time', NF90_DOUBLE, time_dimid, time_varid), 'def_var time' )
         call nc_check( nf90_put_att(ncid, time_varid, 'units',    'months since 1583-01-01'), 'att time units' )
         call nc_check( nf90_put_att(ncid, time_varid, 'calendar', '360_day'),              'att time cal' )
         call nc_check( nf90_put_att(ncid, time_varid, 'standard_name', 'time'),            'att time std' )
         call nc_check( nf90_put_att(ncid, time_varid, 'axis', 'T'),                        'att time axis' )

         do iv = 1, nvar_nc
           call nc_check( nf90_def_var(ncid, trim(ncvars(iv)%nc_name), NF90_DOUBLE, time_dimid, ncvars(iv)%varid), &
                          'def_var '//trim(ncvars(iv)%nc_name) )
           call nc_check( nf90_put_att(ncid, ncvars(iv)%varid, 'long_name', trim(ncvars(iv)%long_name)), &
                          'att long '//trim(ncvars(iv)%nc_name) )
           if (trim(ncvars(iv)%std_name) /= '' .and. trim(ncvars(iv)%std_name) /= '-') &
             call nc_check( nf90_put_att(ncid, ncvars(iv)%varid, 'standard_name', trim(ncvars(iv)%std_name)), &
                            'att std '//trim(ncvars(iv)%nc_name) )
           if (trim(ncvars(iv)%units) /= '' .and. trim(ncvars(iv)%units) /= '-') &
             call nc_check( nf90_put_att(ncid, ncvars(iv)%varid, 'units', trim(ncvars(iv)%units)), &
                            'att units '//trim(ncvars(iv)%nc_name) )
           if (ncvars(iv)%depth >= 0._dblp) &
             call nc_check( nf90_put_att(ncid, ncvars(iv)%varid, 'depth', ncvars(iv)%depth), &
                            'att depth '//trim(ncvars(iv)%nc_name) )
         enddo

         call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'title', 'iLOVECLIM CLIO chronological diagnostics (evolu)'), 'gatt title' )
         call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'source', 'iLOVECLIM / CLIO'), 'gatt source' )
         call nc_check( nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'CF-1.8'), 'gatt conv' )

         call nc_check( nf90_enddef(ncid), 'enddef' )
         rec_count = 0 ; day_count = 0._dblp
         deallocate(rows, exp_rows, exp_lev, covered)
       end subroutine evolu_nc_init

       subroutine evolu_nc_write(vinfor, nvinfo)
         integer(ip), intent(in) :: nvinfo
         real(dblp),  intent(in) :: vinfor(nvinfo)
         integer(ip) :: iv
         if (ncid == -1) return
         rec_count = rec_count + 1
         day_count = day_count + 1._dblp
         call nc_check( nf90_put_var(ncid, time_varid, day_count, start=(/rec_count/)), 'put time' )
         do iv = 1, nvar_nc
           call nc_check( nf90_put_var(ncid, ncvars(iv)%varid, vinfor(ncvars(iv)%reg_index), start=(/rec_count/)), &
                          'put '//trim(ncvars(iv)%nc_name) )
         enddo
       end subroutine evolu_nc_write

       subroutine evolu_nc_close()
         if (ncid /= -1) then
           call nc_check( nf90_close(ncid), 'close' )
           ncid = -1
           if (allocated(ncvars)) deallocate(ncvars)
         endif
       end subroutine evolu_nc_close

       subroutine split_fields(line, fields, nf)
         character(len=*),                    intent(in)  :: line
         character(len=*), dimension(NFIELD), intent(out) :: fields
         integer(ip),                         intent(out) :: nf
         integer(ip) :: i, p, start
         character(len=len(line)) :: seg
         fields(:) = '' ; nf = 0 ; start = 1
         do i = 1, NFIELD
           if (start > len_trim(line) + 1) exit
           p = index(line(start:), ';')
           if (p == 0) then
             seg = line(start:) ; start = len(line) + 1
           else
             seg = line(start:start+p-2) ; start = start + p
           endif
           nf = nf + 1 ; fields(nf) = trim(adjustl(seg))
           if (p == 0) exit
         enddo
       end subroutine split_fields

       subroutine read_namelist(path, rows, nrows, unit_out)
         character(len=*),           intent(in)  :: path
         type(nml_row), allocatable, intent(out) :: rows(:)
         integer(ip),                intent(out) :: nrows
         integer(ip),                intent(in)  :: unit_out
         character(len=512) :: line
         character(len=LEN_LONG), dimension(NFIELD) :: f
         integer(ip) :: u, ios, nf, ncount, pass
         logical :: exists
         inquire(file=path, exist=exists)
         if (.not. exists) then
           write(unit_out,*) 'STOP evolu_netcdf: namelist not found: ', trim(path)
           error stop 'evolu_netcdf: namelist not found'
         endif
         do pass = 1, 2
           ncount = 0
           open(newunit=u, file=path, status='old', action='read')
           do
             read(u,'(A)', iostat=ios) line
             if (ios /= 0) exit
             if (len_trim(line) == 0) cycle
             if (line(1:1) == '#') cycle
             if (verify(line, ' ') == 0) cycle
             line = adjustl(line)
             if (line(1:1) == '#') cycle
             call split_fields(line, f, nf)
             if (nf /= NFIELD) then
               write(unit_out,*) 'STOP evolu_netcdf: malformed row (', nf, ' fields): ', trim(line)
               error stop 'evolu_netcdf: malformed namelist row'
             endif
             ncount = ncount + 1
             if (pass == 2) then
               rows(ncount)%internal  = f(1) ; rows(ncount)%nc_name  = f(2)
               rows(ncount)%long_name = f(3) ; rows(ncount)%std_name = f(4)
               rows(ncount)%units     = f(5)
             endif
           enddo
           close(u)
           if (pass == 1) then ; nrows = ncount ; allocate(rows(nrows)) ; endif
         enddo
       end subroutine read_namelist

       subroutine expand_row(row, ks1, ks2, out_rows, nout, lev_index)
         type(nml_row),               intent(in)  :: row
         integer(ip),                 intent(in)  :: ks1, ks2
         type(nml_row), dimension(:), intent(out) :: out_rows
         integer(ip),                 intent(out) :: nout
         integer(ip),   dimension(:), intent(out) :: lev_index
         integer(ip) :: k, m
         character(len=8) :: idx_disp, idx_pad, idx_plain
         if (index(row%internal, '%') == 0) then
           out_rows(1) = row ; lev_index(1) = 0 ; nout = 1 ; return
         endif
         nout = 0
         do k = ks1, ks2
           m = k - ks1 + 1
           write(idx_disp, '(I2)')   m
           write(idx_pad,  '(I2.2)') m
           write(idx_plain,'(I0)')   m
           nout = nout + 1
           out_rows(nout)%internal  = subst_keep(row%internal,  '%', idx_disp(1:2))
           out_rows(nout)%nc_name   = subst_trim(row%nc_name,   '%', trim(idx_pad))
           out_rows(nout)%long_name = subst_trim(row%long_name, '%', trim(adjustl(idx_plain)))
           out_rows(nout)%std_name  = row%std_name
           out_rows(nout)%units     = row%units
           lev_index(nout)          = m
         enddo
       end subroutine expand_row

       function subst_keep(s, c, repl) result(r)
         character(len=*), intent(in) :: s, repl
         character(len=1), intent(in) :: c
         character(len=LEN_LONG) :: r
         integer(ip) :: p
         p = index(s, c)
         if (p == 0) then ; r = s ; else ; r = s(1:p-1) // repl // s(p+1:) ; endif
       end function subst_keep

       function subst_trim(s, c, repl) result(r)
         character(len=*), intent(in) :: s, repl
         character(len=1), intent(in) :: c
         character(len=LEN_LONG) :: r
         integer(ip) :: p
         p = index(s, c)
         if (p == 0) then ; r = s ; else ; r = s(1:p-1) // repl // s(p+1:) ; endif
       end function subst_trim

       subroutine nc_check(status, wherep)
         use newunit_clio_mod, only: clio3_out_id
         integer(ip),      intent(in) :: status
         character(len=*), intent(in) :: wherep
         if (status /= nf90_noerr) then
           write(clio3_out_id,*) 'STOP evolu_netcdf: ', trim(wherep), ' : ', trim(nf90_strerror(status))
           error stop 'evolu_netcdf: NetCDF error'
         endif
       end subroutine nc_check

      end module evolu_netcdf_mod

#endif /* EVOLU_NETCDF */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
