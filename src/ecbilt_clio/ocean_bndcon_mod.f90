!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!   Copyright for initial code Hugues Goosse, Thierry Fichefet, https://www.elic.ucl.ac.be/modx/index.php?id=289
!      Branched from version 1.2 (detseaalb, initseaalb, oc2at, ec_shine and the apply part of ec_oc2co)

!   Further development:
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
!      MODULE: [ocean_bndcon_mod]
!
!>     @brief Ocean / sea-ice surface boundary conditions for the atmosphere, held on the coupler side (T21 grid).
!
!      DESCRIPTION:
!>     The ocean model (CLIO) fills ocn_bndcon on the atmospheric grid through its gather routine (ec_oc2co_gather in
!!     OCEAN2COUPL_COM). ocean_bndcon_apply then derives the surface-type state seen by ECBilt (fractn, tsurfn, albesn
!!     for noc and nse). ocn_bndcon is the only ocean -> atmosphere interface, hence what is recorded and replayed in
!!     ECBilt standalone mode (ROADMAP step 0a, driven by ocean_coupling_mod). Nothing in this module depends on CLIO.
!!
!!     Contract (public entry points):
!!       ocean_bndcon_apply()           : ocn_bndcon -> fractn, tsurfn, albesn over ocean and sea ice.
!!       ocean_bndcon_write(kday,path)  : append ocn_bndcon as the state of absolute day kday (record mode).
!!       ocean_bndcon_read(kday,path)   : set ocn_bndcon to the recorded state of absolute day kday (replay mode).
!!       ocean_bndcon_read_clim(kday,file) : set ocn_bndcon to the day of year of kday in a 360-record file (climatology).
!!       initseaalb()                   : read the seasonal zonal-mean open-sea albedos (albsea).
!!       oc2at(fin,fout)                : interpolate a field from the CLIO grid to the atmospheric grid.
!!
!>     Files: one NetCDF file per model year, <path>ocean_bndcon_<year>.nc (year as i6.6, "m" prefix if negative),
!!     record = day of year (1..360) of the ocean state, fields stored exactly as double, layout (lon, lat, time).
!!     kday is the number of days completed since model year 0: the state gathered at the start of day i of a run
!!     labelled irunlabel is kday = irunlabel*360 + i-1. A climatology file has the same layout with exactly 360 records;
!!     it is produced outside the model (e.g. NCO nces over complete recorded years) and its year is ignored.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module ocean_bndcon_mod

        use global_constants_mod, only: str_len, ip, dblp=>dp
        use comatm,               only: nlat, nlon

        implicit none

        private

        public :: ocn_bndcon
        public :: N_OCN_BNDCON, BNDCON_SST, BNDCON_SIC, BNDCON_TSI, BNDCON_HIC, BNDCON_HSN
        public :: ocean_bndcon_apply, ocean_bndcon_write, ocean_bndcon_read, ocean_bndcon_read_clim, initseaalb, oc2at

! dmr&clo   Fields of ocn_bndcon (third index).
        integer(ip), parameter :: N_OCN_BNDCON = 5   ! number of fields
        integer(ip), parameter :: BNDCON_SST   = 1   ! sea surface temperature              [K]
        integer(ip), parameter :: BNDCON_SIC   = 2   ! sea-ice fraction of the ocean part   [1]   (1 - CLIO albq)
        integer(ip), parameter :: BNDCON_TSI   = 3   ! sea-ice surface temperature          [K]
        integer(ip), parameter :: BNDCON_HIC   = 4   ! sea-ice thickness                    [m]
        integer(ip), parameter :: BNDCON_HSN   = 5   ! snow thickness on sea ice            [m]

! dmr&clo   NetCDF names, units and long names of the fields, in the order above.
        character(len=8),  dimension(N_OCN_BNDCON), parameter :: BNDCON_NAME  = [ character(len=8) ::                   &
                                                                  "sst", "sic", "tsi", "hic", "hsn" ]
        character(len=8),  dimension(N_OCN_BNDCON), parameter :: BNDCON_UNIT  = [ character(len=8) ::                   &
                                                                  "K", "1", "K", "m", "m" ]
        character(len=40), dimension(N_OCN_BNDCON), parameter :: BNDCON_LNAME = [ character(len=40) ::                  &
                                                                  "sea surface temperature",                            &
                                                                  "sea-ice fraction of the ocean part",                 &
                                                                  "sea-ice surface temperature",                        &
                                                                  "sea-ice thickness",                                  &
                                                                  "snow thickness on sea ice" ]

        integer(ip),       parameter :: DAYS_PER_YEAR = 360                     ! 360-day model calendar
        character(len=*),  parameter :: BNDCON_FILE_ROOT = "ocean_bndcon_"      ! file name: <path><root><year>.nc

! dmr&clo   The ocean / sea-ice state on the atmospheric grid, as gathered from the ocean model or replayed.
        real(dblp), dimension(nlat,nlon,N_OCN_BNDCON) :: ocn_bndcon

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Derive the atmospheric surface state over ocean and sea ice (formerly the second half of ec_oc2co).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine ocean_bndcon_apply()

          use comemic_mod,    only: fracto
          use comcoup_mod,    only: couptcc
          use comsurf_mod,    only: fractn, nse, noc, tsurfn, albesn, albocef
#if ( DOWNSTS == 1 )
          use vertDownsc_mod, only: tsurfn_d
#endif
#if ( IMSK == 1 )
          use input_icemask,  only: icemaskalb   ! afq: force an ice-sheet albedo where an ice sheet covers the ocean
          use comland_mod,    only: albsnow
#endif

          real(dblp), parameter :: TZERO = 273.15_dblp   ! [K]

          integer(ip) :: ix, jy
          real(dblp)  :: zalb, zalbp
          real(dblp)  :: seaalb(nlat,nlon)

          call detseaalb(seaalb)

          do ix = 1, nlon
            do jy = 1, nlat
              fractn(jy,ix,noc) = (1.0_dblp-ocn_bndcon(jy,ix,BNDCON_SIC))*fracto(jy,ix)
              fractn(jy,ix,nse) = ocn_bndcon(jy,ix,BNDCON_SIC)*fracto(jy,ix)
#if ( DOWNSTS == 1 )
              tsurfn_d(jy,ix,nse,:) = min(TZERO,ocn_bndcon(jy,ix,BNDCON_TSI))
              tsurfn_d(jy,ix,noc,:) = max(TZERO-1.8_dblp,ocn_bndcon(jy,ix,BNDCON_SST))
#endif
              tsurfn(jy,ix,nse) = min(TZERO,ocn_bndcon(jy,ix,BNDCON_TSI))
              tsurfn(jy,ix,noc) = max(TZERO-1.8_dblp,ocn_bndcon(jy,ix,BNDCON_SST))
              call ec_shine(TZERO-0.15_dblp, TZERO-0.25_dblp, tsurfn(jy,ix,nse), ocn_bndcon(jy,ix,BNDCON_HIC),           &
                            ocn_bndcon(jy,ix,BNDCON_HSN), zalb, zalbp)
              albesn(jy,ix,nse) = (1.0_dblp-couptcc(jy,ix))*zalbp + couptcc(jy,ix)*zalb
              albesn(jy,ix,noc) = albocef*seaalb(jy,ix)
#if ( IMSK == 1 )
              if (icemaskalb(jy,ix).gt.0.9_dblp) then
                ! afq: ice sheet over the ocean, force a snow albedo
                albesn(jy,ix,nse) = albsnow(jy)
                albesn(jy,ix,noc) = albsnow(jy)
              endif
#endif
            enddo
          enddo

        end subroutine ocean_bndcon_apply

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Record / replay of ocn_bndcon (see the module description for the file layout and kday).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine ocean_bndcon_write(kday, path)

          use io_nc_mod,            only: IO_NC_FILE, IO_NC_AXIS, IO_GRID_VAR
          use global_constants_mod, only: rad_to_deg
          use comatm,               only: phi

          integer(ip),      intent(in) :: kday   !< absolute day of the ocean state
          character(len=*), intent(in) :: path   !< directory, with trailing "/"

          type(IO_NC_FILE), target, save :: bndcon_file
          type(IO_NC_AXIS),         save :: ax_lon, ax_lat, ax_time
          type(IO_GRID_VAR),        save :: bndcon_var(N_OCN_BNDCON)
          integer(ip),              save :: cur_year = -huge(1_ip)

          integer(ip)            :: year, irec, k, j
          real(dblp)             :: lons(nlon), lats(nlat), fld(nlon,nlat)
          character(len=str_len) :: fname

          year = obnd_year(kday)
          irec = obnd_rec(kday)

          if (year /= cur_year) then
            fname = obnd_filename(path, year)
            call bndcon_file%init(FileName=fname,                                                                         &
                                  OPT_TitleFile="iLOVECLIM ocean/sea-ice boundary conditions for ECBilt (T21), record mode")
            do j = 1, nlon
              lons(j) = 360.0_dblp*real(j-1,dblp)/real(nlon,dblp)
            enddo
            lats(:) = phi(:)*rad_to_deg
            call ax_lon%init("lon", OPT_vals_to_write1D=lons, OPT_AxisUnit="degrees_east")
            call ax_lat%init("lat", OPT_vals_to_write1D=lats, OPT_AxisUnit="degrees_north")
            call ax_time%init("time", OPT_AxisUnit="days since start of year (state at start of day)",                    &
                              OPT_isTime=.true., OPT_calendar="360_day")
            call ax_lon%wrte(fname)
            call ax_lat%wrte(fname)
            call ax_time%wrte(fname)
            do k = 1, N_OCN_BNDCON
              call bndcon_var(k)%init(trim(BNDCON_NAME(k)), bndcon_file, "lon lat time",                                  &
                                      OPT_longname=trim(BNDCON_LNAME(k)), OPT_units=trim(BNDCON_UNIT(k)))
            enddo
            cur_year = year
          endif

          do k = 1, N_OCN_BNDCON
            fld(:,:) = transpose(ocn_bndcon(:,:,k))
            call bndcon_var(k)%wrte(fld, irec)
          enddo

        end subroutine ocean_bndcon_write

        subroutine ocean_bndcon_read(kday, path)

          integer(ip),      intent(in) :: kday   !< absolute day of the ocean state
          character(len=*), intent(in) :: path   !< directory, with trailing "/"

          call obnd_read_record(obnd_filename(path, obnd_year(kday)), obnd_rec(kday), .false.)

        end subroutine ocean_bndcon_read

        subroutine ocean_bndcon_read_clim(kday, fname)

          integer(ip),      intent(in) :: kday    !< absolute day; only its day of year is used
          character(len=*), intent(in) :: fname   !< climatology file (360 records)

          call obnd_read_record(fname, obnd_rec(kday), .true.)

        end subroutine ocean_bndcon_read_clim

        subroutine obnd_read_record(fname, irec, is_clim)

          use io_nc_mod, only: IO_NC_FILE, IO_GRID_VAR

          character(len=*), intent(in) :: fname     !< file to read from
          integer(ip),      intent(in) :: irec      !< record (day of year)
          logical,          intent(in) :: is_clim   !< .true.: the file must hold exactly DAYS_PER_YEAR records

          type(IO_NC_FILE), target, save :: bndcon_file
          type(IO_GRID_VAR),        save :: bndcon_var(N_OCN_BNDCON)
          character(len=str_len),   save :: cur_fname = ""
          integer(ip),              save :: cur_nrec  = 0

          integer(ip) :: k
          real(dblp)  :: fld(nlon,nlat)

          if (fname /= cur_fname) then
            call bndcon_file%open(fname)
            ! grid check: nrec reads the length of any named dimension
            if (bndcon_file%nrec("lon") /= nlon .or. bndcon_file%nrec("lat") /= nlat) then
              write(*,*) "ocean_bndcon_read: ", trim(fname), " is not on the ", nlon, " x ", nlat, " (lon, lat) grid"
              stop 1
            endif
            cur_nrec = bndcon_file%nrec()
            if (is_clim .and. cur_nrec /= DAYS_PER_YEAR) then
              write(*,*) "ocean_bndcon_read: climatology file ", trim(fname), " has ", cur_nrec, " records instead of ",   &
                         DAYS_PER_YEAR
              stop 1
            endif
            do k = 1, N_OCN_BNDCON
              call bndcon_var(k)%init(trim(BNDCON_NAME(k)), bndcon_file, "lon lat time")
            enddo
            cur_fname = fname
          endif

          if (irec > cur_nrec) then
            write(*,*) "ocean_bndcon_read: record ", irec, " not in ", trim(fname), " (", cur_nrec, " records).",        &
                       " Was the recording complete for this year?"
            stop 1
          endif

          do k = 1, N_OCN_BNDCON
            call bndcon_var(k)%read(fld, irec)
            ocn_bndcon(:,:,k) = transpose(fld(:,:))
          enddo

        end subroutine obnd_read_record

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Indexing helpers: model year and record (day of year) of an absolute day, file name of a year.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        function obnd_year(kday) result(year)
          integer(ip), intent(in) :: kday
          integer(ip)             :: year

          year = (kday - modulo(kday, DAYS_PER_YEAR)) / DAYS_PER_YEAR   ! floor(kday/360), also for kday < 0
        end function obnd_year

        function obnd_rec(kday) result(irec)
          integer(ip), intent(in) :: kday
          integer(ip)             :: irec

          irec = modulo(kday, DAYS_PER_YEAR) + 1_ip
        end function obnd_rec

        function obnd_filename(path, year) result(fname)
          character(len=*), intent(in) :: path
          integer(ip),      intent(in) :: year
          character(len=str_len)       :: fname

          character(len=8) :: cyear

          if (year >= 0) then
            write(cyear,'(i6.6)') year
          else
            write(cyear,'("m",i6.6)') -year
          endif
          fname = trim(path)//BNDCON_FILE_ROOT//trim(cyear)//'.nc'
        end function obnd_filename

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Open-sea albedo: seasonal climatology and its interpolation to the current day (legacy CLIO code, moved).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine detseaalb(seaalb)

          ! albedo of the open sea as a function of the time of the year, linearly interpolated between seasonal means

          use comemic_mod, only: iday, imonth, albsea

          real(dblp), intent(out) :: seaalb(nlat,nlon)

          integer(ip) :: i, j, id1, is1, is2
          real(dblp)  :: albseaz(nlat)
          real(dblp)  :: sfrac

          id1 = (imonth-1)*30+iday-14
          if (id1.lt.1) id1 = id1+360

          is1 = (id1+89)/90
          is2 = is1+1
          if (is2.eq.5) is2 = 1

          sfrac = (id1-((is1-1)*90.0_dblp+1.0_dblp))/90.0_dblp

          do j = 1, nlat
            albseaz(j) = albsea(j,is1)+(albsea(j,is2)-albsea(j,is1))*sfrac
          enddo

          do j = 1, nlon
            do i = 1, nlat
              seaalb(i,j) = albseaz(i)
            enddo
          enddo

        end subroutine detseaalb

        subroutine initseaalb()

          ! read climatological zonal mean albedos for each season

          use comemic_mod, only: albsea

          integer(ip) :: i, is
          integer(ip) :: albedo_dat_id

          open(newunit=albedo_dat_id, file='inputdata/clio/albedo.dat')

          read(albedo_dat_id,*)
          do i = 1, nlat
            read(albedo_dat_id,45) (albsea(i,is), is=1,4)
          enddo
45        format(4(2x,f7.4))
          close(albedo_dat_id)

        end subroutine initseaalb

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Ocean grid -> atmospheric grid interpolation (legacy CLIO code, moved).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine oc2at(fin, fout)

          use comcoup_mod, only: ijocn, kamax, wo2a, indo2a

          real(dblp), intent(in)  :: fin(ijocn)          !< field on the CLIO grid (passed by sequence association)
          real(dblp), intent(out) :: fout(nlat,nlon)     !< field on the atmospheric grid

          integer(ip) :: ji, jk, i, j
          real(dblp)  :: zsum

          ji = 0
          do i = 1, nlat
            do j = 1, nlon
              zsum = 0.0_dblp
              ji = ji+1
              do jk = 1, kamax
                zsum = zsum + wo2a(ji,jk) * fin(indo2a(ji,jk))
              enddo
              fout(i,j) = zsum
            enddo
          enddo

        end subroutine oc2at

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&clo   Snow / sea-ice albedo following Shine & Henderson-Sellers (1985) (legacy CLIO code, moved).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        subroutine ec_shine(tfsn, tfsg, ts, hgbq, hnbq, zalb, zalbp)

          ! albin : albedo of melting ice in the Arctic
          ! albis : albedo of melting ice in the Antarctic (Shine & Henderson-Sellers, 1985)
          ! albice: albedo of melting ice
          ! alphd : albedo of snow (thickness > 0.05 m)
          ! alphdi: albedo of thick bare ice
          ! alphs : albedo of melting snow
          ! cgren : correction of the snow or ice albedo for the effect of cloudiness (Grenfell & Perovich, 1984)

          use comphys, only: albice, alphd, alphdi, alphs, cgren

          real(dblp), intent(in)  :: tfsn    !< melting point temperature of snow (273.15 in CLIO)
          real(dblp), intent(in)  :: tfsg    !< melting point temperature of ice (273.05 in CLIO)
          real(dblp), intent(in)  :: ts      !< surface temperature
          real(dblp), intent(in)  :: hgbq    !< ice thickness
          real(dblp), intent(in)  :: hnbq    !< snow thickness
          real(dblp), intent(out) :: zalb    !< ice/snow albedo for overcast sky
          real(dblp), intent(out) :: zalbp   !< ice/snow albedo for clear sky

          real(dblp) :: al

          if (hnbq.gt.0.0_dblp) then
            ! case of ice covered by snow
            if (ts.lt.tfsn) then
              ! freezing snow
              if (hnbq.gt.0.05_dblp) then
                zalbp = alphd
              else
                if (hgbq.gt.1.5_dblp) then
                  zalbp = alphdi+(hnbq*(alphd-alphdi)/0.05_dblp)
                else if (hgbq.gt.1.0_dblp.and.hgbq.le.1.5_dblp) then
                  al = 0.472_dblp+2.0_dblp*(alphdi-0.472_dblp)*(hgbq-1.0_dblp)
                else if (hgbq.gt.0.05_dblp.and.hgbq.le.1.0_dblp) then
                  al = 0.2467_dblp+(0.7049_dblp*hgbq)-(0.8608_dblp*(hgbq*hgbq))+(0.3812_dblp*(hgbq*hgbq*hgbq))
                else
                  al = 0.1_dblp+3.6_dblp*hgbq
                endif
                if (hgbq.le.1.5_dblp) zalbp = al+(hnbq*(alphd-al)/0.05_dblp)
              endif
            else
              ! melting snow
              if (hnbq.ge.0.1_dblp) then
                zalbp = alphs
              else
                zalbp = albice+((alphs-albice)/0.1_dblp)*hnbq
              endif
            endif
          else
            ! case of ice free of snow
            if (ts.lt.tfsg) then
              ! freezing ice
              if (hgbq.gt.1.5_dblp) then
                zalbp = alphdi
              else if (hgbq.gt.1.0_dblp.and.hgbq.le.1.5_dblp) then
                zalbp = 0.472_dblp+2.0_dblp*(alphdi-0.472_dblp)*(hgbq-1.0_dblp)
              else if (hgbq.gt.0.05_dblp.and.hgbq.le.1.0_dblp) then
                zalbp = 0.2467_dblp+(0.7049_dblp*hgbq)-(0.8608_dblp*(hgbq*hgbq))+(0.3812_dblp*(hgbq*hgbq*hgbq))
              else
                zalbp = 0.1_dblp+3.6_dblp*hgbq
              endif
            else
              ! melting ice
              if (hgbq.gt.1.5_dblp) then
                zalbp = albice
              else if (hgbq.gt.1.0_dblp.and.hgbq.le.1.5_dblp) then
                zalbp = 0.472_dblp+(2.0_dblp*(albice-0.472_dblp)*(hgbq-1.0_dblp))
              else if (hgbq.gt.0.05_dblp.and.hgbq.le.1.0_dblp) then
                zalbp = 0.2467_dblp+0.7049_dblp*hgbq-(0.8608_dblp*(hgbq*hgbq))+(0.3812_dblp*(hgbq*hgbq*hgbq))
              else
                zalbp = 0.1_dblp+3.6_dblp*hgbq
              endif
            endif
          endif
          zalb = zalbp+cgren

        end subroutine ec_shine

      end module ocean_bndcon_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
