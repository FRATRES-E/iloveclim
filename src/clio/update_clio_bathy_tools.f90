!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of iLOVECLIM/LUDUS
!!      update_clio_bathy_tools is free software: you can redistribute it and/or modify it under the terms of the GNU General
!!      Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
!!      version.
!!
!!      update_clio_bathy_tools is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
!!      implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!!      See the GNU General Public License for more details.
!!
!!      You should have received a copy of the GNU General Public License along with update_clio_bathy_tools.
!!      If not, see <http://www.gnu.org/licenses/>.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"
!dmr -- Adding the choice of components through the pre-processing options

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: update_clio_bathy_tools
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module update_clio_bathy_tools is handling the update, writing and reading tools for all necessary aspects of
!!             regarding the update of the bathymetry in CLIO as a transient feature
!
!>     @date Creation date: January, 21st, 2019
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

      module update_clio_bathy_tools

       ! Size of maximum character length & double precision definition
       use global_constants_mod, only: str_len, dblp=>dp, ip

       ! CLIO dimension variables
       use para0_mod, only: imax, jmax, kmax

       !
       implicit none
       private

       public :: get_clio_masks_fnm
       public :: write_clio_masks, read_clio_masks
       public :: update_diff_masks
       public :: mean_neighbours_with_mask
#if ( OCYCC == 1 )
       public :: mean_neighbours_with_mask_CC
#endif
       public :: mise_a_zero
       public :: mise_a_valeur
       public :: update_prev

      ! NOTE_avoid_public_non_parameter_variables_if_possible

      ! MODULE VARIABLES

      ! tracer and speed masks at the previous update: could be before updating the bathy or at the previous restart point
      real(dblp), dimension(imax,jmax,kmax),public :: tms_prev, tmu_prev
      real(dblp), dimension(imax,jmax)      :: Stms_prev, Stmu_prev
      integer(ip), dimension(imax,jmax,kmax),public :: diff_tms, diff_tmu
      integer(ip), dimension(imax,jmax)     :: Sdiff_tms, Sdiff_tmu

      integer(ip), parameter :: masktms_2D = 1, masktmu_2D = 2, masktms_3D = 1, masktmu_3D = 2

      character(len=str_len), parameter :: genmask_resname="res_masks", gen_ext=".om"

      real(kind=dblp), public :: salglobal

      real(dblp), dimension(imax,jmax),public :: aire_prev
      real(dblp), dimension(kmax),public :: dz_prev

#if ( BATHY >= 2 )
      ! frequency of bathymetry update in years
      integer(ip), parameter :: update_time_bathy = 50000
      integer(ip), parameter, public :: reload_bathy = update_time_bathy * 360

      !Bering
      integer(ip),parameter, public ::  bering_nb = 43
      real(dblp), dimension(bering_nb), public :: bering_date
      real(dblp), dimension(bering_nb), public :: bering_value 

      !kamax
      integer(ip),parameter, public ::  kamax_nb = 43
      real(dblp), dimension(kamax_nb), public :: kamax_date
      real(dblp), dimension(kamax_nb), public :: kamax_value 

#endif

#if ( BATHY >= 2 || NC_BERG == 2 )
      integer,public :: la_date
#endif

!#if ( F_PALAEO_FWF == 1 && APPLY_UNCORFWF == 0 ) 
#if ( F_PALAEO_FWF >= 1 || CONSEAU == 1 ) 
      real(dblp), dimension(imax,jmax), public :: sum_flux_out! flux in Sv for output
#endif


      ! do we need to update fields?
      logical :: need_update = .false.

      interface mean_neighbours_with_mask
          module procedure mean_neighbours_with_mask2D, mean_neighbours_with_mask3D
      end interface mean_neighbours_with_mask

#if ( OCYCC == 1 )
      interface mean_neighbours_with_mask_CC
          module procedure mean_neighbours_with_mask2D_CC, mean_neighbours_with_mask3D_CC
      end interface mean_neighbours_with_mask_CC
#endif

      interface mise_a_zero
          module procedure mise_a_zero_2D, mise_a_zero_3D
      end interface mise_a_zero

      interface mise_a_valeur
          module procedure mise_a_valeur_2D, mise_a_valeur_3D
      end interface mise_a_valeur



      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: get_clio_masks_fnm()
!
!>     @brief This function returns the name of the restart mask file based on the irunlabel given
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function get_clio_masks_fnm(in_runlabel) result(masks_nm)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  in_runlabel: the irunlabel clio used for the current timestep filename
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer, intent(in) :: in_runlabel
       character(str_len)  :: masks_nm

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       character(str_len) :: int_to_str

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       write(int_to_str,'(I6.6)') in_runlabel

       masks_nm = trim(genmask_resname)//trim(int_to_str)//trim(gen_ext)

      end function get_clio_masks_fnm

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: write_clio_masks
!
!>     @brief This subroutine let CLIO write tracer and velocity masks into a separated restart file in the view of having it
!>              updated upon restart if needed
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine write_clio_masks(resfileName)

       use file_libs, only: fileDescriptor, open_f, close_f

       !dmr --- tracer & speed masks to be written out
       use bloc0_mod, only: tms, tmu, dz
       use bloc_mod, only: aire

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  resfileName  The name of the restart file to be created that should contain the masks
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       character(*), intent(in)    :: resfileName

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       type(fileDescriptor)   :: rescliomasks_f

! dmr --- File name variable (nm)
       character(len=str_len) :: rescliomasks_nm

! dmr --- Formatted flag: .false. is unformatted file
       logical, parameter     :: rescliomasks_fm   = .false.

       logical :: success

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       rescliomasks_nm=resfileName
!       rescliomasks_nm='startdata/'//resfileName
       rescliomasks_f%isFormatted=rescliomasks_fm

       call open_f(rescliomasks_f,rescliomasks_nm)

       write(*,*) "File: "//rescliomasks_f%fileName//" is opened at: ",rescliomasks_f%id

       write(rescliomasks_f%id) tms
       write(rescliomasks_f%id) tmu

       write(rescliomasks_f%id) aire
       write(rescliomasks_f%id) dz

       call close_f(rescliomasks_f)

      end subroutine write_clio_masks

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: read_clio_masks
!
!>     @brief This subroutine let CLIO read the previous tracer and velocity masks into a separated restart file and load it in
!>              dedicated module variable (the latter are private)
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine read_clio_masks(resfileName)

       use file_libs, only: fileDescriptor, open_f, close_f

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  resfileName  The name of the restart file to be created that should contain the masks
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       character(*), intent(in)    :: resfileName

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       type(fileDescriptor)   :: rescliomasks_f

! dmr --- File name variable (nm)
       character(len=str_len) :: rescliomasks_nm

! dmr --- Formatted flag: .false. is unformatted file
       logical, parameter     :: rescliomasks_fm   = .false.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       rescliomasks_nm=resfileName
!       rescliomasks_nm='startdata/'//resfileName
       rescliomasks_f%isFormatted=rescliomasks_fm

       call open_f(rescliomasks_f,rescliomasks_nm)

       write(*,*) "File: "//rescliomasks_f%fileName//" is opened at: ",rescliomasks_f%id

       read(rescliomasks_f%id) tms_prev
       read(rescliomasks_f%id) tmu_prev

       read(rescliomasks_f%id) aire_prev
       read(rescliomasks_f%id) dz_prev

!       write(*,*) 'tms_prev', tms_prev(:,:,kmax)

       call close_f(rescliomasks_f)

       Stms_prev(:,:) = tms_prev(:,:,kmax)
       Stmu_prev(:,:) = tmu_prev(:,:,kmax)


      end subroutine read_clio_masks

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: update_diff_masks
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

      function update_diff_masks(tms_new, tmu_new) result(state)

        use ncio, only: nc_create, nc_write_dim, nc_write

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  tms_new, tmu_new to update the difference with from the tms_prev, tmu_prev
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real, dimension(imax,jmax,kmax), intent(in)    :: tms_new, tmu_new
       logical                                        :: state

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       character(len=str_len) :: filename

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       diff_tms(:,:,:) = int(tms_new(:,:,:) - tms_prev(:,:,:))
       diff_tmu(:,:,:) = int(tmu_new(:,:,:) - tmu_prev(:,:,:))

       Sdiff_tms(:,:) = diff_tms(:,:,kmax)
       Sdiff_tmu(:,:) = diff_tmu(:,:,kmax)


!       write(*,*) "diff_tms pathologique", diff_tms(40,44,:)
!       write(*,*) "diff_tmu pathologique", diff_tmu(40,44,:)


!       write(*,*) "Sprev", tms_prev(39,44,15), tms_prev(40,44,15), tms_prev(41,44,15)
!       write(*,*) "Snnew", tms_new(39,44,15), tms_new(40,44,15), tms_new(41,44,15)
!       write(*,*) "Sprev", tms_prev(40,43,15), tms_prev(40,44,15), tms_prev(40,45,15)
!       write(*,*) "Snnew", tms_new(40,43,15), tms_new(40,44,15), tms_new(40,45,15)


!       write(*,*) "******"
!       write(*,*) "Uprev", tmu_prev(39,44,15), tmu_prev(40,44,15), tmu_prev(41,44,15)
!       write(*,*) "Unnew", tmu_new(39,44,15), tmu_new(40,44,15), tmu_new(41,44,15)
!       write(*,*) "Uprev", tmu_prev(40,43,15), tmu_prev(40,44,15), tmu_prev(40,45,15)
!       write(*,*) "Unnew", tmu_new(40,43,15), tmu_new(40,44,15), tmu_new(40,45,15)
!       read(*,*)

       need_update = sum(abs(diff_tms)).ge.1

       write(*,*)
       write(*,*) "Updated CLIO masks successfully"
       write(*,*) "Need_update= ", need_update, sum(abs(diff_tms))
       write(*,*)

       state = need_update

#if (1)
       filename="test_masks_diff_CLIO.nc"

       call nc_create(filename,overwrite=.TRUE.,netcdf4=.TRUE.)
       call nc_write_dim(filename,"x",x=1,dx=1,nx=imax,units="index")
       call nc_write_dim(filename,"y",x=1,dx=1,nx=jmax,units="index")
       call nc_write_dim(filename,"z",x=1,dx=1,nx=kmax,units="index")
       call nc_write(filename,"diff_tms",diff_tms,dim1="x",dim2="y",dim3="z")
       call nc_write(filename,"diff_tmu",diff_tmu,dim1="x",dim2="y",dim3="z")
       call nc_write(filename,"Sdiff_tms",Sdiff_tms,dim1="x",dim2="y")
       call nc_write(filename,"Sdiff_tmu",Sdiff_tmu,dim1="x",dim2="y")
#endif

!cnb update tms_prev et tmu_prev
!      tmu_prev(:,:,:)=tmu_new(:,:,:)
!      tms_prev(:,:,:)=tms_new(:,:,:)

      end function update_diff_masks

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: update_prev
!
!>     @brief This function 
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function update_prev(tms_new, tmu_new) result(state)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       real, dimension(imax,jmax,kmax), intent(in)  :: tms_new, tmu_new
       logical                                   :: state

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|




!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       tmu_prev(:,:,:)=tmu_new(:,:,:)
       tms_prev(:,:,:)=tms_new(:,:,:)


       Stms_prev(:,:) = tms_prev(:,:,kmax)
       Stmu_prev(:,:) = tmu_prev(:,:,kmax)

       write(*,*)
       write(*,*) "Updated CLIO prev  masks"
       write(*,*)


       state=.true.

      end function update_prev

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: mise_a_zero_2D
!
!>     @brief This function allows setting a 2D field to ... zero!
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function mise_a_zero_2D(var, typ_mask) result(state)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real,    dimension(:,:), intent(inout)    :: var
       integer,                 intent(in)       :: typ_mask
       logical                                   :: state

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer            :: i,j, i_min, i_max, j_min, j_max
       integer, parameter :: nb_neigh = 1
       integer, parameter :: new_ocean = 1

       integer(ip), dimension(imax,jmax) :: where_update
       real(dblp),  dimension(imax,jmax) :: mask_prev

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       select case (typ_mask)
         case(masktms_2D)
             where_update(:,:) = Sdiff_tms
             mask_prev(:,:) = Stms_prev(:,:)
         case(masktmu_2D)
             where_update(:,:) = Sdiff_tmu
             mask_prev(:,:) = Stmu_prev(:,:)
       end select

       i_min = LBOUND(var, dim=1)
       i_max = UBOUND(var, dim=1)

       j_min = LBOUND(var, dim=2)
       j_max = UBOUND(var, dim=2)

       do i=i_min, i_max
         do j=j_min, j_max

            if (where_update(i,j).eq.new_ocean) then

               var(i,j) = 0.0d0

            endif ! on where_update
         enddo
       enddo

      end function mise_a_zero_2D

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: mise_a_zero_3D
!
!>     @brief This function allows setting a 3D field to ... zero!
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function mise_a_zero_3D(var, typ_mask) result(state)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real,    dimension(:,:,:), intent(inout)  :: var
       integer,                   intent(in)     :: typ_mask
       logical                                   :: state

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer            :: i,j, i_min, i_max, j_min, j_max, k, k_min, k_max
       integer, parameter :: new_ocean = 1

       integer(ip), dimension(imax,jmax,kmax) :: where_update
       real(dblp),  dimension(imax,jmax,kmax) :: mask_prev

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       select case (typ_mask)
         case(masktms_3D)
             where_update(:,:,:) = diff_tms(:,:,:)
             mask_prev(:,:,:) = tms_prev(:,:,:)
         case(masktmu_3D)
             where_update(:,:,:) = diff_tmu(:,:,:)
             mask_prev(:,:,:) = tmu_prev(:,:,:)
       end select

       i_min = LBOUND(var, dim=1)
       i_max = UBOUND(var, dim=1)

       j_min = LBOUND(var, dim=2)
       j_max = UBOUND(var, dim=2)

       k_min = LBOUND(var, dim=3)
       k_max = UBOUND(var, dim=3)

       do k=k_min, k_max

       do i=i_min, i_max
         do j=j_min, j_max

            if (where_update(i,j,k).eq.new_ocean) then

                var(i,j,k) = 0.0d0

            endif ! on where_update
         enddo
       enddo
       enddo

      end function mise_a_zero_3D

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: mise_a_zero_2D
!
!>     @brief This function allows setting a 2D field to ... zero!
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function mise_a_valeur_2D(var, typ_mask, valeur) result(state)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real,    dimension(:,:), intent(inout)    :: var
       integer,                 intent(in)       :: typ_mask
       logical                                   :: state
       real                                      :: valeur

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer            :: i,j, i_min, i_max, j_min, j_max
       integer, parameter :: new_ocean = 1


       integer(ip), dimension(imax,jmax) :: where_update
       real(dblp),  dimension(imax,jmax) :: mask_prev

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       select case (typ_mask)
         case(masktms_2D)
             where_update(:,:) = Sdiff_tms
             mask_prev(:,:) = Stms_prev(:,:)
         case(masktmu_2D)
             where_update(:,:) = Sdiff_tmu
             mask_prev(:,:) = Stmu_prev(:,:)
       end select

       i_min = LBOUND(var, dim=1)
       i_max = UBOUND(var, dim=1)

       j_min = LBOUND(var, dim=2)
       j_max = UBOUND(var, dim=2)

       do i=i_min, i_max
         do j=j_min, j_max

            if (where_update(i,j).eq.new_ocean) then

               var(i,j) = valeur

            endif ! on where_update
         enddo
       enddo

      end function mise_a_valeur_2D

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function mise_a_valeur_3D(var, typ_mask, valeur) result(state)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real,    dimension(:,:,:), intent(inout)  :: var
       integer,                   intent(in)     :: typ_mask
       logical                                   :: state
       real                                      :: valeur

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer            :: i,j, i_min, i_max, j_min, j_max, k, k_min, k_max
       integer, parameter :: new_ocean = 1

       integer(ip), dimension(imax,jmax,kmax) :: where_update
       real(dblp),  dimension(imax,jmax,kmax) :: mask_prev

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       select case (typ_mask)
         case(masktms_3D)
             where_update(:,:,:) = diff_tms(:,:,:)
             mask_prev(:,:,:) = tms_prev(:,:,:)
         case(masktmu_3D)
             where_update(:,:,:) = diff_tmu(:,:,:)
             mask_prev(:,:,:) = tmu_prev(:,:,:)
       end select

       i_min = LBOUND(var, dim=1)
       i_max = UBOUND(var, dim=1)

       j_min = LBOUND(var, dim=2)
       j_max = UBOUND(var, dim=2)

       k_min = LBOUND(var, dim=3)
       k_max = UBOUND(var, dim=3)

       do k=k_min, k_max

       do i=i_min, i_max
         do j=j_min, j_max

            if (where_update(i,j,k).eq.new_ocean) then

                var(i,j,k) = valeur

            endif ! on where_update
         enddo
       enddo
       enddo

      end function mise_a_valeur_3D

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: mean_neighbours_with_mask2D
!
!>     @brief This function allows computing the mean value of a field using neighbours and a mask
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function mean_neighbours_with_mask2D(var, typ_mask) result(state)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real,    dimension(:,:), intent(inout)    :: var
       integer,                 intent(in)       :: typ_mask
       logical                                   :: state

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer            :: i,j, i_min, i_max, j_min, j_max, w, z, i_new, j_new
       !integer, parameter :: nb_neigh = 1
       integer            :: nb_neigh
       integer, parameter :: new_ocean = 1
       integer            :: nb_mean
       real(dblp)         :: value_mean
       logical            :: good_one, good_two

       integer(ip), dimension(imax,jmax) :: where_update
       real(dblp),  dimension(imax,jmax) :: mask_prev

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       select case (typ_mask)
         case(masktms_2D)
             where_update(:,:) = Sdiff_tms
             mask_prev(:,:) = Stms_prev(:,:)
         case(masktmu_2D)
             where_update(:,:) = Sdiff_tmu
             mask_prev(:,:) = Stmu_prev(:,:)
       end select

       i_min = LBOUND(var, dim=1)
       i_max = UBOUND(var, dim=1)

       j_min = LBOUND(var, dim=2)
       j_max = UBOUND(var, dim=2)

       do i=i_min, i_max
         do j=j_min, j_max

            if (where_update(i,j).eq.new_ocean) then
               value_mean = 0.0
               nb_mean = 0

               !on cherche autour des voisins avec des valeurs
               !oceaniques
               nb_neigh=0
               !write(*,*) "nb_neigh, nb_mean", nb_neigh, nb_mean
               do while ((nb_neigh.le.10) .and. (nb_mean.eq.0))
                nb_neigh=nb_neigh+1
                !write(*,*) "nb_neigh, nb_mean", nb_neigh, nb_mean

               do w=nb_neigh*(-1), nb_neigh
                 do z=nb_neigh*(-1), nb_neigh
                    i_new = i+w
                    j_new = j+z
                    good_one = in_range(i_new,i_min,i_max,circular=.FALSE.)
                    good_two = in_range(j_new,j_min,j_max,circular=.TRUE.)
                    if ((good_one.and.good_two).and.(mask_prev(i_new,j_new).eq.1.0)) then
                       value_mean = value_mean + var(i_new,j_new)
                       nb_mean = nb_mean + 1
                    endif

                 enddo
               enddo

               end do !while

             if (nb_mean.gt.0) then
!                write(*,*) "ping 2D, updated one cell", var(i,j), value_mean / nb_mean, nb_mean, i, j, imax, jmax
               var(i,j) = value_mean / nb_mean
!nb
             else
                write(*,*) "ping 2D, nb_mean=0"
!               var(i,j) = 0.0
             endif

            endif ! on where_update
         enddo
       enddo

      end function mean_neighbours_with_mask2D

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: mean_neighbours_with_mask3D
!
!>     @brief This function allows computing the mean value of a field using neighbours and a mask
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function mean_neighbours_with_mask3D(var, typ_mask) result(state)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real,    dimension(:,:,:), intent(inout)  :: var
       integer,                   intent(in)     :: typ_mask
       logical                                   :: state

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer            :: i,j, i_min, i_max, j_min, j_max, w, z, i_new, j_new, k, k_min, k_max
       !integer, parameter :: nb_neigh = 1
       integer            :: nb_neigh
       integer, parameter :: new_ocean = 1
       integer            :: nb_mean
       real(dblp)         :: value_mean
       logical            :: good_one, good_two

       integer(ip), dimension(imax,jmax,kmax) :: where_update
       real(dblp),  dimension(imax,jmax,kmax) :: mask_prev

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       select case (typ_mask)
         case(masktms_3D)
             where_update(:,:,:) = diff_tms(:,:,:)
             mask_prev(:,:,:) = tms_prev(:,:,:)
         case(masktmu_3D)
             where_update(:,:,:) = diff_tmu(:,:,:)
             mask_prev(:,:,:) = tmu_prev(:,:,:)
       end select

       i_min = LBOUND(var, dim=1)

       j_min = LBOUND(var, dim=2)
       j_max = UBOUND(var, dim=2)

       k_min = LBOUND(var, dim=3)
       k_max = UBOUND(var, dim=3)

       do k=k_min, k_max

       do i=i_min, i_max
         do j=j_min, j_max

            if (where_update(i,j,k).eq.new_ocean) then
               value_mean = 0.0
               nb_mean = 0

               !on cherche autour des voisins avec des valeurs
               !oceaniques
               nb_neigh=0
               !write(*,*) "nb_neigh, nb_mean", nb_neigh, nb_mean
               do while ((nb_neigh.le.10) .and. (nb_mean.eq.0))
                nb_neigh=nb_neigh+1
                !write(*,*) "nb_neigh, nb_mean", nb_neigh, nb_mean

                do w=nb_neigh*(-1), nb_neigh
                 do z=nb_neigh*(-1), nb_neigh
                    i_new = i+w
                    j_new = j+z
                    good_one = in_range(i_new,i_min,i_max,circular=.FALSE.)
                    good_two = in_range(j_new,j_min,j_max,circular=.TRUE.)
                    if ((good_one.and.good_two)) then
                      if ((mask_prev(i_new,j_new,k).eq.1.0)) then
                       value_mean = value_mean + var(i_new,j_new,k)
                       nb_mean = nb_mean + 1
                      endif
                    endif

                 enddo
                enddo

                !write(*,*) "nb_neigh, nb_mean", nb_neigh, nb_mean
               end do !while

             if (nb_mean.gt.0) then
!               write(*,*) "ping 3D, updated one cell", var(i,j,k), value_mean / nb_mean, nb_mean, i, j, imax, jmax, k
               var(i,j,k) = value_mean / nb_mean
!nb
             else
               write(*,*) "ping 3D, nb_mean=0,", i,j,k
!               var(i,j,k) = 0.0
             endif

            endif ! on where_update
         enddo
       enddo
       enddo

      end function mean_neighbours_with_mask3D
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( OCYCC == 1 )

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: mean_neighbours_with_mask2D_CC
!
!>     @brief This function allows computing the mean value of a field using neighbours and a mask
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function mean_neighbours_with_mask2D_CC(var, typ_mask) result(state)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       use declars_mod
       use loveclim_transfer_mod, only : diff_MGT, MGT_prev

       real,    dimension(:,:), intent(inout)    :: var
       integer,                 intent(in)       :: typ_mask
       logical                                   :: state

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer            :: i,j, i_min, i_max, j_min, j_max, w, z, i_new, j_new
       !integer, parameter :: nb_neigh = 1
       integer            :: nb_neigh
       integer, parameter :: new_ocean = 1
       integer            :: nb_mean
       real(dblp)         :: value_mean
       logical            :: good_one, good_two

       integer(ip), dimension(LT,NOC_CBR) :: where_update
       real(dblp),  dimension(LT,NOC_CBR) :: mask_prev

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!       select case (typ_mask)
!         case(masktms_2D)
             where_update(:,:) = diff_MGT(:,1,:) ! j=1 for surface in OCYCC
             mask_prev(:,:) = MGT_prev(:,1,:)
!             where_update(:,:) = Sdiff_tms
!             mask_prev(:,:) = Stms_prev(:,:)
!         case(masktmu_2D)
!             where_update(:,:) = Sdiff_tmu
!             mask_prev(:,:) = Stmu_prev(:,:)
!       end select

       i_min = LBOUND(var, dim=1)
       i_max = UBOUND(var, dim=1)

       j_min = LBOUND(var, dim=2)
       j_max = UBOUND(var, dim=2)

       do i=i_min, i_max
         do j=j_min, j_max

            if (where_update(i,j).eq.new_ocean) then
               value_mean = 0.0
               nb_mean = 0

               !on cherche autour des voisins avec des valeurs
               !oceaniques
               nb_neigh=0
               !write(*,*) "nb_neigh, nb_mean", nb_neigh, nb_mean
               do while ((nb_neigh.le.10) .and. (nb_mean.eq.0))
                nb_neigh=nb_neigh+1
                !write(*,*) "nb_neigh, nb_mean", nb_neigh, nb_mean

               do w=nb_neigh*(-1), nb_neigh
                 do z=nb_neigh*(-1), nb_neigh
                    i_new = i+w
                    j_new = j+z
                    good_one = in_range(i_new,i_min,i_max,circular=.FALSE.)
                    good_two = in_range(j_new,j_min,j_max,circular=.TRUE.)
                    if ((good_one.and.good_two).and.(mask_prev(i_new,j_new).eq.1.0)) then
                       value_mean = value_mean + var(i_new,j_new)
                       nb_mean = nb_mean + 1
                    endif

                 enddo
               enddo

               end do !while

             if (nb_mean.gt.0) then
!                write(*,*) "ping 2D, updated one cell", var(i,j), value_mean / nb_mean, nb_mean, i, j, imax, jmax
               var(i,j) = value_mean / nb_mean
!nb
             else
                write(*,*) "ping 2D, nb_mean=0"
!               var(i,j) = 0.0
             endif

            endif ! on where_update
         enddo
       enddo

      end function mean_neighbours_with_mask2D_CC

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: mean_neighbours_with_mask3D_CC
!
!>     @brief This function allows computing the mean value of a field using neighbours and a mask
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function mean_neighbours_with_mask3D_CC(var, typ_mask) result(state)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       use loveclim_transfer_mod, only : diff_MGT, MGT_prev
       use declars_mod

       real,    dimension(:,:,:), intent(inout)  :: var
       integer,                   intent(in)     :: typ_mask
       logical                                   :: state

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer            :: i,j, i_min, i_max, j_min, j_max, w, z, i_new, j_new, k, k_min, k_max
       !integer, parameter :: nb_neigh = 1
       integer            :: nb_neigh
       integer, parameter :: new_ocean = 1
       integer            :: nb_mean
       real(dblp)         :: value_mean
       logical            :: good_one, good_two

!nb index are different in OCYCC
!       integer(ip), dimension(jmax,kmax,imax) :: where_update
!       real(dblp),  dimension(jmax,kmax,imax) :: mask_prev
       integer(ip), dimension(LT,JT,NOC_CBR) :: where_update
       real(dblp),  dimension(LT,JT,NOC_CBR) :: mask_prev

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!       select case (typ_mask)
!         case(masktms_3D)
             where_update(:,:,:) = diff_MGT(:,:,:)
             mask_prev(:,:,:) = MGT_prev(:,:,:)
!         case(masktmu_3D)
!             where_update(:,:,:) = diff_tmu(:,:,:)
!             mask_prev(:,:,:) = tmu_prev(:,:,:)
!       end select

       i_min = LBOUND(var, dim=1)
       i_max = UBOUND(var, dim=1)

       j_min = LBOUND(var, dim=2)
       j_max = UBOUND(var, dim=2)

       k_min = LBOUND(var, dim=3)
       k_max = UBOUND(var, dim=3)

       do k=k_min, k_max

       do i=i_min, i_max
         do j=j_min, j_max

            if (where_update(i,j,k).eq.new_ocean) then
               value_mean = 0.0
               nb_mean = 0

               !on cherche autour des voisins avec des valeurs
               !oceaniques
               nb_neigh=0
               !write(*,*) "nb_neigh, nb_mean", nb_neigh, nb_mean
               do while ((nb_neigh.le.10) .and. (nb_mean.eq.0))
                nb_neigh=nb_neigh+1
                !write(*,*) "nb_neigh, nb_mean", nb_neigh, nb_mean

               do w=nb_neigh*(-1), nb_neigh
                 do z=nb_neigh*(-1), nb_neigh
                    i_new = i+w
                    j_new = j+z
                    good_one = in_range(i_new,i_min,i_max,circular=.FALSE.)
                    good_two = in_range(j_new,j_min,j_max,circular=.TRUE.)
                    if ((good_one.and.good_two).and.(mask_prev(i_new,j_new,k).eq.1.0)) then
                       value_mean = value_mean + var(i_new,j_new,k)
                       nb_mean = nb_mean + 1
                    endif

                 enddo
               enddo

               end do !while

             if (nb_mean.gt.0) then
!               write(*,*) "ping 3D CC, updated one cell", var(i,j,k), value_mean / nb_mean, nb_mean, i, j, imax, jmax, k
               var(i,j,k) = value_mean / nb_mean
!nb
             else
               write(*,*) "ping 3D CC,  nb_mean=0"
!               var(i,j,k) = 0.0
             endif

            endif ! on where_update
         enddo
       enddo
       enddo

      end function mean_neighbours_with_mask3D_CC

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
#endif



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: in_range
!
!>     @brief This function returns the given integer within a certain range to ensure appropriate array indexing
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function in_range(indx_val,indx_min,indx_max,circular) result(returnStatus)

!~        use AnotherModule_mod, only: some_variable
!~        use AnotherModule_mod, only: some_otherfunction

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer(ip), intent(inout)           :: indx_val
       integer(ip), intent(in)              :: indx_min, indx_max
       logical, optional, intent(in)        :: circular
       logical                              :: returnStatus
       integer, parameter                   :: big = 9999

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical :: circularly

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       if (.not.present(circular)) then
          circularly = .false.
       else
          circularly = circular
       endif

       returnStatus = .true.

       if (indx_val.gt.indx_max) then
          if (circularly) then
             indx_val = indx_min + indx_val-indx_max-1 ! indx_val - indx_max is positive => e.g. 65-64 = 1 - 1 = 0 => indx_min
          else
             indx_val = big
             returnStatus = .false.
          endif
       elseif (indx_val.lt.indx_min) then
          if (circularly) then
             indx_val = indx_max - (indx_min - indx_val) + 1
          else
             indx_val = big
             returnStatus = .false.
          endif
       endif
      end function in_range

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


















!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: [Routine name here]
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

      function someFunction(inParam, outParam, inOutParam) result(returnValue)

!~        use AnotherModule_mod, only: some_variable
!~        use AnotherModule_mod, only: some_otherfunction

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa   By reference variables ...
!>    @param[in]  inParam  The Beaufitul Parameter that does all the input !!
!>    @param[out] outParam The Beaufitul Parameter that does all the output!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real, intent(in)    :: inParam
       real, intent(inout) :: inOutParam
       real, intent(out)   :: outParam

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! tfa  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real                :: returnValue
       real                :: someVariable


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!>    @bug Description of the stupid sticky bug that we know exist there but is not corrected yet!

      end function someFunction

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|






!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module update_clio_bathy_tools

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
