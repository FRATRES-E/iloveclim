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

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: [ec_iso2co]
!
!>     @author  Thomas Extier (tex)
!>     @author  Didier M. Roche (dmr)

!>     @brief This module communicates the variables needed to calculate the d18Oatm to the coupler
!
!>     @date Creation date: June, 12th 2018
!>     @date Last modification: June, 26th 2018
!>     @author Last modified by : tex&dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

module ec_iso2co

#if (OXYISO > 0)

  use global_constants_mod, only : dp, ip, nbdays => days_year365d_i
  use taillesGrilles,       only : nlat => iEcb, nlon => jEcb
  use ec_co2ca,             only : npft

  implicit none

  real(kind=dp), dimension(nlat,nlon,nbdays)      :: rhu_photo
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: Tnew
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: temp_photo
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: epsieq
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: d18Oleaf
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: Ci
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: GPPnew
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: gpp_C3
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: gpp_C4
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: GPPO2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: epsil
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: PHOo2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: Pmehler
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: O2prod
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: fphoto
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: f_darkleaves
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: a_darksoil
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: alpha18ter
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: d18Oter
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: diffO2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: m_18o2
  real(kind=dp), dimension(nlat,nlon,npft)        :: epsieq_PFT2
  real(kind=dp), dimension(nlat,nlon,nbdays,npft) :: d18Oleaf_PFT2

#if (WAXISO > 0)
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: epsiplus_2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: dDlw_2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: dDblw_2
  real(kind=dp), dimension(nlat,nlon,nbdays,npft) :: dDwax_2
#endif

#if (D17ISO > 0)
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: alpha17_2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: lambda_t2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: R18Olw2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: R17Olw2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: d17Oleaf2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: d17Oter_local2
  real(kind=dp), dimension(nlat,nlon,nbdays,npft) :: d17Oleaf_PFT2
  real(kind=dp), dimension(nlat,nlon,nbdays,npft) :: d17Otl_PFT2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: R18ter2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: R17ter2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: d17Oter
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: R18Ogw2
  real(kind=dp), dimension(nlat,nlon,nbdays)      :: R17Ogw2
#endif

#endif

  contains

#if (OXYISO > 0)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7-|
  subroutine write_d18O_iso
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7-|
! ***
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7-|

  use global_constants_mod, only  : str_len
  use taillesGrilles,       only  : nlat => iEcb, nlon => jEcb
  use ec_co2ca,             only  : frac_land, j_antarctique, npft
  use comatm,               only  : nsh2, nvl, ntl, nsh, nm
  use file_libs,            only  : fileDescriptor, open_f, close_f

#if ( COMATM == 1 )
  use comsurf_mod,          only  : nld, fractn
  use comemic_mod,          only  : iyear
#endif

  implicit none

#if ( COMATM == 0 )
#include "comsurf.h" /* nld, fractn */
#include "comemic.h" /* iyear */
#endif

  character(len=10)       :: cyear
  integer(ip)             :: i,j,k,l,ilat,ird, gauss_asc_id
  real(dp)                :: phi(nlat),rj,rw,dumwei,ari(nlat/2)
  real(dp), parameter     :: pi = 2*acos(0.)
  real(dp), parameter     :: radian_to_degree = 180_dp/pi
  real(dp), allocatable, dimension(:) :: templon,templat

  type(fileDescriptor)    :: rhu_photo_f, Tnew_f, temp_photo_f &
 & , epsieq_f, d18Oleaf_f, Ci_f, gppo2_f, epsil_f, phoo2_f, pmehler_f &
 & , o2prod_f, fphoto_f, fdarkl_f, a_darks_f, a18ter_f, d18Oter_f     &
 & , diffO2_f, m_18o2_f, gpp_C3_f, gpp_C4_f, epsieq_PFT2_f            &
 & , d18Oleaf_PFT2_f

  character(len=str_len)   ::    							          &
              rhu_photo_nm     = "outputdata/d18Oatm/rhu_photo",      &
              Tnew_nm          = "outputdata/d18Oatm/Tnew",           &
              temp_photo_nm    = "outputdata/d18Oatm/temp_photo",     &
              epsieq_nm        = "outputdata/d18Oatm/epsieq",         &
              d18Oleaf_nm      = "outputdata/d18Oatm/d18Oleaf",       &
              Ci_nm            = "outputdata/d18Oatm/Ci",             &
              gppo2_nm         = "outputdata/d18Oatm/gpp_oxy",        &
              epsil_nm         = "outputdata/d18Oatm/epsil",          &
              phoo2_nm         = "outputdata/d18Oatm/pho_oxy",        &
              pmehler_nm       = "outputdata/d18Oatm/Pmehler",        &
              o2prod_nm        = "outputdata/d18Oatm/O2prod",         &
              fphoto_nm        = "outputdata/d18Oatm/fphoto",         &
              fdarkl_nm        = "outputdata/d18Oatm/f_darkleaves",   &
              a_darks_nm       = "outputdata/d18Oatm/a_darksoil",     &
              a18ter_nm        = "outputdata/d18Oatm/alpha18ter",     &
              d18Oter_nm       = "outputdata/d18Oatm/d18Oter",        &
              diffO2_nm        = "outputdata/d18Oatm/diffO2",         &
              m_18o2_nm        = "outputdata/d18Oatm/m_18o2",         &
              gpp_C3_nm        = "outputdata/d18Oatm/gpp_C3",         &
              gpp_C4_nm        = "outputdata/d18Oatm/gpp_C4",         &
              epsieq_PFT2_nm   = "outputdata/d18Oatm/epsieq_PFT",     &
              d18Oleaf_PFT2_nm = "outputdata/d18Oatm/d18Oleaf_PFT"

  logical, parameter       ::                                         &
               rhu_photo_fm     = .true. ,  						  &
               Tnew_fm          = .true.,                             &
               temp_photo_fm    = .true.,                             &
               epsieq_fm        = .true.,                             &
               d18Oleaf_fm      = .true.,                             &
               Ci_fm            = .true.,                             &
               gppo2_fm         = .true.,                             &
               epsil_fm         = .true.,                             &
               phoo2_fm         = .true.,                             &
               pmehler_fm       = .true.,                             &
               o2prod_fm        = .true.,                             &
               fphoto_fm        = .true.,                             &
               fdarkl_fm        = .true.,                             &
               a_darks_fm       = .true.,                             &
               a18ter_fm        = .true.,                             &
               d18Oter_fm       = .true.,                             &
               diffO2_fm        = .true.,                             &
               m_18o2_fm        = .true.,                             &
               gpp_C3_fm        = .true.,                             &
               gpp_C4_fm        = .true.,                             &
               epsieq_PFT2_fm   = .true.,                             &
               d18Oleaf_PFT2_fm = .true.



   open(newunit=gauss_asc_id, file='inputdata/gauss.asc',status='old',form='formatted')
      rewind(gauss_asc_id)
      ilat=nlat/2
   10 continue
        read(gauss_asc_id,220,end=15) rj,rw
        ird=int(rj)
        if (ird.eq.ilat) then
          do k=1,ird
            read(gauss_asc_id,220) ari(k),dumwei
          enddo
          goto 20
        else
          goto 10
        endif
   15 continue
   20 continue

      do j=1,ilat
        phi(j)=-ari(ilat+1-j)
        phi(ilat+j)=ari(j)

      enddo

      do j=1,nlat
        phi(j)=asin(phi(j))
      enddo

 220  format(f18.10,f17.10)
    close(gauss_asc_id)

	allocate(templon(nlon))
      templon = (/ ((360.0d0*l)/nlon,l=0,nlon-1) /)
    allocate(templat(nlat))
	  templat = phi(1:nlat)*radian_to_degree



  write(cyear,'(i10)'),iyear

  rhu_photo_f%isFormatted = rhu_photo_fm
  call open_f(rhu_photo_f, trim(rhu_photo_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
!  Tnew_f%isFormatted = Tnew_fm
!  call open_f(Tnew_f, trim(Tnew_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  temp_photo_f%isFormatted = temp_photo_fm
  call open_f(temp_photo_f, trim(temp_photo_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
!  epsieq_f%isFormatted = epsieq_fm
!  call open_f(epsieq_f, trim(epsieq_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
!  epsieq_PFT2_f%isFormatted = epsieq_PFT2_fm
!  call open_f(epsieq_PFT2_f, trim(epsieq_PFT2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  d18Oleaf_f%isFormatted = d18Oleaf_fm
  call open_f(d18Oleaf_f, trim(d18Oleaf_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  d18Oleaf_PFT2_f%isFormatted = d18Oleaf_PFT2_fm
  call open_f(d18Oleaf_PFT2_f, trim(d18Oleaf_PFT2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  Ci_f%isFormatted = Ci_fm
  call open_f(Ci_f, trim(Ci_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  gppo2_f%isFormatted = gppo2_fm
  call open_f(gppo2_f, trim(gppo2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
!  epsil_f%isFormatted = epsil_fm
!  call open_f(epsil_f, trim(epsil_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
!  phoo2_f%isFormatted = phoo2_fm
!  call open_f(phoo2_f, trim(phoo2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
!  pmehler_f%isFormatted = pmehler_fm
!  call open_f(pmehler_f, trim(pmehler_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
!  o2prod_f%isFormatted = o2prod_fm
!  call open_f(o2prod_f, trim(o2prod_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  fphoto_f%isFormatted = fphoto_fm
  call open_f(fphoto_f, trim(fphoto_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  fdarkl_f%isFormatted = fdarkl_fm
  call open_f(fdarkl_f, trim(fdarkl_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  a_darks_f%isFormatted = a_darks_fm
  call open_f(a_darks_f, trim(a_darks_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  a18ter_f%isFormatted = a18ter_fm
  call open_f(a18ter_f, trim(a18ter_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
!  d18Oter_f%isFormatted = d18Oter_fm
!  call open_f(d18Oter_f, trim(d18Oter_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
!  diffO2_f%isFormatted = diffO2_fm
!  call open_f(diffO2_f, trim(diffO2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
!  m_18o2_f%isFormatted = m_18o2_fm
!  call open_f(m_18o2_f, trim(m_18o2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
!  gpp_C3_f%isFormatted = gpp_C3_fm
!  call open_f(gpp_C3_f, trim(gpp_C3_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
!  gpp_C4_f%isFormatted = gpp_C4_fm
!  call open_f(gpp_C4_f, trim(gpp_C4_nm)//trim(adjustl(cyear))//".dat", o_stat="new")


  do i=1,nlat
    do j=1,nlon

      if(fractn(i,j,nld).GE.frac_land.AND.i.GE.j_antarctique) then
		 write(rhu_photo_f%id,'(367f10.4)') templon(j),templat(i),rhu_photo(i,j,1:365)
!		 write(Tnew_f%id,'(367f10.4)') templon(j),templat(i),Tnew(i,j,1:365)
		 write(temp_photo_f%id,'(367f10.4)') templon(j),templat(i),temp_photo(i,j,1:365)
!		 write(epsieq_f%id,'(367f10.4)') templon(j),templat(i),epsieq(i,j,1:365)
!		 write(epsieq_PFT2_f%id,'(367f10.4)') templon(j),templat(i),epsieq_PFT2(i,j,1:26)
		 write(d18Oleaf_f%id,'(367f10.4)') templon(j),templat(i),d18Oleaf(i,j,1:365)
!		 write(d18Oleaf_PFT2_f%id,'(367f10.4)') templon(j),templat(i),d18Oleaf_PFT2(i,j,1:365,1:26)
		 write(Ci_f%id,'(367f10.4)') templon(j),templat(i),Ci(i,j,1:365)
		 write(gppo2_f%id,'(367f10.4)') templon(j),templat(i),GPPO2(i,j,1:365)
!		 write(epsil_f%id,'(367f10.4)') templon(j),templat(i),epsil(i,j,1:365)
!		 write(phoo2_f%id,'(367f10.4)') templon(j),templat(i),PHOo2(i,j,1:365)
!		 write(pmehler_f%id,'(367f10.4)') templon(j),templat(i),Pmehler(i,j,1:365)
!		 write(o2prod_f%id,'(367f10.4)') templon(j),templat(i),O2prod(i,j,1:365)
		 write(fphoto_f%id,'(367f10.4)') templon(j),templat(i),fphoto(i,j,1:365)
		 write(fdarkl_f%id,'(367f10.4)') templon(j),templat(i),f_darkleaves(i,j,1:365)
		 write(a_darks_f%id,'(367f10.4)') templon(j),templat(i),a_darksoil(i,j,1:365)
		 write(a18ter_f%id,'(367f10.4)') templon(j),templat(i),alpha18ter(i,j,1:365)
! 		 write(d18Oter_f%id,'(367f10.4)') templon(j),templat(i),d18Oter(i,j,1:365)
!		 write(diffO2_f%id,'(367f10.4)') templon(j),templat(i),diffO2(i,j,1:365)
! 		 write(m_18o2_f%id,'(367f10.4)') templon(j),templat(i),m_18o2(i,j,1:365)
!		 write(gpp_C3_f%id,'(367f10.4)') templon(j),templat(i),gpp_C3(i,j,1:365)
!		 write(gpp_C4_f%id,'(367f10.4)') templon(j),templat(i),gpp_C4(i,j,1:365)
	  endif

	enddo
  enddo

  call close_f(rhu_photo_f)
!  call close_f(Tnew_f)
  call close_f(temp_photo_f)
!  call close_f(epsieq_f)
!  call close_f(epsieq_PFT2_f)
  call close_f(d18Oleaf_f)
!  call close_f(d18Oleaf_PFT2_f)
  call close_f(Ci_f)
  call close_f(gppo2_f)
!  call close_f(epsil_f)
!  call close_f(phoo2_f)
!  call close_f(pmehler_f)
!  call close_f(o2prod_f)
  call close_f(fdarkl_f)
  call close_f(a_darks_f)
  call close_f(a18ter_f)
!  call close_f(d18Oter_f)
!  call close_f(diffO2_f)
!  call close_f(m_18o2_f)
!  call close_f(gpp_C3_f)
!  call close_f(gpp_C4_f)

  deallocate(templon)
  deallocate(templat)

  end subroutine write_d18O_iso

#endif

#if (WAXISO > 0)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7-|
  subroutine write_dDwax_iso
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7-|
! ***
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7-|

  use global_constants_mod, only  : str_len
  use taillesGrilles,       only  : nlat => iEcb, nlon => jEcb
  use ec_co2ca,             only  : frac_land, j_antarctique, npft
  use comatm,               only  : nsh2, nvl, ntl, nsh, nm
  use file_libs,            only  : fileDescriptor, open_f, close_f

#if ( COMATM == 1 )
  use comsurf_mod,          only  : nld, fractn
  use comemic_mod,          only  : iyear
#endif

  implicit none

#if ( COMATM == 0 )
#include "comsurf.h" /* nld, fractn */
#include "comemic.h" /* iyear */
#endif

  character(len=10)       :: cyear
  integer(ip)             :: i,j,k,l,ilat,ird, gauss_asc_id
  real(dp)                :: phi(nlat),rj,rw,dumwei,ari(nlat/2)
  real(dp), parameter     :: pi = 2*acos(0.)
  real(dp), parameter     :: radian_to_degree = 180_dp/pi
  real(dp), allocatable, dimension(:) :: templon,templat

  type(fileDescriptor)    :: epsiplus_2_f, dDlw_2_f, dDblw_2_f, dDwax_2_f

  character(len=str_len)  ::                                          &
              epsiplus_2_nm  = "outputdata/dDwax/epsiplus",           &
              dDlw_2_nm      = "outputdata/dDwax/dDlw",               &
              dDblw_2_nm     = "outputdata/dDwax/dDblw",              &
              dDwax_2_nm     = "outputdata/dDwax/dDwax"

   logical, parameter     ::                                          &
              epsiplus_2_fm  = .true.,                                &
              dDlw_2_fm      = .true.,                                &
              dDblw_2_fm     = .true.,                                &
              dDwax_2_fm     = .true.

   open(newunit=gauss_asc_id, file='inputdata/gauss.asc',status='old',form='formatted')
      rewind(gauss_asc_id)
      ilat=nlat/2
   10 continue
        read(gauss_asc_id,220,end=15) rj,rw
        ird=int(rj)
        if (ird.eq.ilat) then
          do k=1,ird
            read(gauss_asc_id,220) ari(k),dumwei
          enddo
          goto 20
        else
          goto 10
        endif
   15 continue
   20 continue

      do j=1,ilat
        phi(j)=-ari(ilat+1-j)
        phi(ilat+j)=ari(j)

      enddo

      do j=1,nlat
        phi(j)=asin(phi(j))
      enddo

 220  format(f18.10,f17.10)
    close(gauss_asc_id)

	allocate(templon(nlon))
      templon = (/ ((360.0d0*l)/nlon,l=0,nlon-1) /)
    allocate(templat(nlat))
	  templat = phi(1:nlat)*radian_to_degree



  write(cyear,'(i10)'),iyear

  epsiplus_2_f%isFormatted = epsiplus_2_fm
  call open_f(epsiplus_2_f, trim(epsiplus_2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  dDlw_2_f%isFormatted = dDlw_2_fm
  call open_f(dDlw_2_f, trim(dDlw_2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  dDblw_2_f%isFormatted = dDblw_2_fm
  call open_f(dDblw_2_f, trim(dDblw_2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  dDwax_2_f%isFormatted = dDwax_2_fm
  call open_f(dDwax_2_f, trim(dDwax_2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")

  do i=1,nlat
    do j=1,nlon

      if(fractn(i,j,nld).GE.frac_land.AND.i.GE.j_antarctique) then
         write(epsiplus_2_f%id,'(367f10.4)') templon(j),templat(i),epsiplus_2(i,j,1:365)
         write(dDlw_2_f%id,'(367f10.4)') templon(j),templat(i),dDlw_2(i,j,1:365)
         write(dDblw_2_f%id,'(367f10.4)') templon(j),templat(i),dDblw_2(i,j,1:365)
         write(dDwax_2_f%id,'(9600f10.4)') templon(j),templat(i),dDwax_2(i,j,1:365,1:26)
	  endif

	enddo
  enddo


  call close_f(epsiplus_2_f)
  call close_f(dDlw_2_f)
  call close_f(dDblw_2_f)
  call close_f(dDwax_2_f)

  deallocate(templon)
  deallocate(templat)

  end subroutine write_dDwax_iso
#endif

#if (D17ISO > 0)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7-|
  subroutine write_d17O_iso
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7-|
! ***
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7-|

  use global_constants_mod, only  : str_len
  use taillesGrilles,       only  : nlat => iEcb, nlon => jEcb
  use ec_co2ca,             only  : frac_land, j_antarctique, npft
  use comatm,               only  : nsh2, nvl, ntl, nsh, nm
  use file_libs,            only  : fileDescriptor, open_f, close_f

#if ( COMATM == 1 )
  use comsurf_mod,          only  : nld, fractn
  use comemic_mod,          only  : iyear
#endif

  implicit none

#if ( COMATM == 0 )
#include "comsurf.h" /* nld, fractn */
#include "comemic.h" /* iyear */
#endif

  character(len=10)       :: cyear
  integer(ip)             :: i,j,k,l,ilat,ird, gauss_asc_id
  real(dp)                :: phi(nlat),rj,rw,dumwei,ari(nlat/2)
  real(dp), parameter     :: pi = 2*acos(0.)
  real(dp), parameter     :: radian_to_degree = 180_dp/pi
  real(dp), allocatable, dimension(:) :: templon,templat

  type(fileDescriptor)    :: alpha17_2_f,lambda_t2_f     &
   & , R18Olw2_f,R17Olw2_f,d17Oleaf2_f,d17Oter_local2_f  &
   & ,d17Oleaf_PFT2_f,d17Otl_PFT2_f,R18ter2_f,R17ter2_f  &
   & ,d17Oter_f,R18Ogw2_f,R17Ogw2_f

  character(len=str_len)  ::                                               &
              alpha17_2_nm        = "outputdata/D17O/alpha17",             &
              lambda_t2_nm        = "outputdata/D17O/lambda_t",            &
              R18Olw2_nm          = "outputdata/D17O/R18Olw",              &
              R17Olw2_nm          = "outputdata/D17O/R17Olw",              &
              d17Oleaf2_nm        = "outputdata/D17O/d17Oleaf",            &
              d17Oter_local2_nm   = "outputdata/D17O/d17Oter_local",       &
              d17Oleaf_PFT2_nm    = "outputdata/D17O/d17Oleaf_PFT",        &
              d17Otl_PFT2_nm      = "outputdata/D17O/d17Otl_PFT",          &
              R18ter2_nm          = "outputdata/D17O/R18ter",              &
              R17ter2_nm          = "outputdata/D17O/R17ter",              &
              d17Oter_nm          = "outputdata/D17O/d17Oter",             &
              R18Ogw2_nm          = "outputdata/D17O/R18Ogw",              &
              R17Ogw2_nm          = "outputdata/D17O/R17Ogw"


   logical, parameter     ::                                          &
              alpha17_2_fm         = .true.,                          &
              lambda_t2_fm         = .true.,                          &
              R18Olw2_fm           = .true.,                          &
              R17Olw2_fm           = .true.,                          &
              d17Oleaf2_fm         = .true.,                          &
              d17Oter_local2_fm    = .true.,                          &
              d17Oleaf_PFT2_fm     = .true.,                          &
              d17Otl_PFT2_fm       = .true.,                          &
              R18ter2_fm           = .true.,                          &
              R17ter2_fm           = .true.,                          &
              d17Oter_fm           = .true.,                          &
              R18Ogw2_fm           = .true.,                          &
              R17Ogw2_fm           = .true.


   open(newunit=gauss_asc_id, file='inputdata/gauss.asc',status='old',form='formatted')
      rewind(gauss_asc_id)
      ilat=nlat/2
   10 continue
        read(gauss_asc_id,220,end=15) rj,rw
        ird=int(rj)
        if (ird.eq.ilat) then
          do k=1,ird
            read(gauss_asc_id,220) ari(k),dumwei
          enddo
          goto 20
        else
          goto 10
        endif
   15 continue
   20 continue

      do j=1,ilat
        phi(j)=-ari(ilat+1-j)
        phi(ilat+j)=ari(j)

      enddo

      do j=1,nlat
        phi(j)=asin(phi(j))
      enddo

 220  format(f18.10,f17.10)
    close(gauss_asc_id)

	allocate(templon(nlon))
      templon = (/ ((360.0d0*l)/nlon,l=0,nlon-1) /)
    allocate(templat(nlat))
	  templat = phi(1:nlat)*radian_to_degree



  write(cyear,'(i10)'),iyear

  alpha17_2_f%isFormatted = alpha17_2_fm
  call open_f(alpha17_2_f, trim(alpha17_2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  lambda_t2_f%isFormatted = lambda_t2_fm
  call open_f(lambda_t2_f, trim(lambda_t2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  R18Olw2_f%isFormatted = R18Olw2_fm
  call open_f(R18Olw2_f, trim(R18Olw2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  R17Olw2_f%isFormatted = R17Olw2_fm
  call open_f(R17Olw2_f, trim(R17Olw2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  d17Oleaf2_f%isFormatted = d17Oleaf2_fm
  call open_f(d17Oleaf2_f, trim(d17Oleaf2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  d17Oter_local2_f%isFormatted = d17Oter_local2_fm
  call open_f(d17Oter_local2_f, trim(d17Oter_local2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  d17Oleaf_PFT2_f%isFormatted = d17Oleaf_PFT2_fm
  call open_f(d17Oleaf_PFT2_f, trim(d17Oleaf_PFT2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  d17Otl_PFT2_f%isFormatted = d17Otl_PFT2_fm
  call open_f(d17Otl_PFT2_f, trim(d17Otl_PFT2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  R18ter2_f%isFormatted = R18ter2_fm
  call open_f(R18ter2_f, trim(R18ter2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  R17ter2_f%isFormatted = R17ter2_fm
  call open_f(R17ter2_f, trim(R17ter2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  d17Oter_f%isFormatted = d17Oter_fm
  call open_f(d17Oter_f, trim(d17Oter_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  R18Ogw2_f%isFormatted = R18Ogw2_fm
  call open_f(R18Ogw2_f, trim(R18Ogw2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")
  R17Ogw2_f%isFormatted = R17Ogw2_fm
  call open_f(R17Ogw2_f, trim(R17Ogw2_nm)//trim(adjustl(cyear))//".dat", o_stat="new")


  do i=1,nlat
    do j=1,nlon

      if(fractn(i,j,nld).GE.frac_land.AND.i.GE.j_antarctique) then
         write(alpha17_2_f%id,'(367f10.4)') templon(j),templat(i),alpha17_2(i,j,1:365)
         write(lambda_t2_f%id,'(367f10.4)') templon(j),templat(i),lambda_t2(i,j,1:365)
         write(R18Olw2_f%id,'(367f10.4)') templon(j),templat(i),R18Olw2(i,j,1:365)
         write(R17Olw2_f%id,'(367f10.4)') templon(j),templat(i),R17Olw2(i,j,1:365)
         write(d17Oleaf2_f%id,'(367f10.4)') templon(j),templat(i),d17Oleaf2(i,j,1:365)
         write(d17Oter_local2_f%id,'(367f10.4)') templon(j),templat(i),d17Oter_local2(i,j,1:365)
         write(d17Oleaf_PFT2_f%id,'(9600f16.4)') templon(j),templat(i),d17Oleaf_PFT2(i,j,1:365,1:26)
         write(d17Otl_PFT2_f%id,'(9600f16.4)') templon(j),templat(i),d17Otl_PFT2(i,j,1:365,1:26)
         write(R18ter2_f%id,'(367f10.4)') templon(j),templat(i),R18ter2(i,j,1:365)
         write(R17ter2_f%id,'(367f10.4)') templon(j),templat(i),R17ter2(i,j,1:365)
         write(d17Oter_f%id,'(367f10.4)') templon(j),templat(i),d17Oter(i,j,1:365)
         write(R18Ogw2_f%id,'(367f10.4)') templon(j),templat(i),R18Ogw2(i,j,1:365)
         write(R17Ogw2_f%id,'(367f10.4)') templon(j),templat(i),R17Ogw2(i,j,1:365)

      endif

    enddo
  enddo


  call close_f(alpha17_2_f)
  call close_f(lambda_t2_f)
  call close_f(R18Olw2_f)
  call close_f(R17Olw2_f)
  call close_f(d17Oleaf2_f)
  call close_f(d17Oter_local2_f)
  call close_f(d17Oleaf_PFT2_f)
  call close_f(d17Otl_PFT2_f)
  call close_f(R18ter2_f)
  call close_f(R17ter2_f)
  call close_f(d17Oter_f)
  call close_f(R18Ogw2_f)
  call close_f(R17Ogw2_f)


  deallocate(templon)
  deallocate(templat)

  end subroutine write_d17O_iso
#endif


end module ec_iso2co

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
