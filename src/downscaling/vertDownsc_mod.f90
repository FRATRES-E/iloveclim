#include "choixcomposantes.h"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module contient les variables pour le downscaling vertical sur grille ECBilt pour akkumulation avant passage
!       à une sous-grille. Ce module ne doit dépendre en aucun cas de variables sous-grilles.
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche
!      Date   : 27 Mai 2009
!      Derniere modification : 27 Jan 2016
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       module vertDownsc_mod

       use comatm,            only: nlat, nlon
       use ecbilt_topography, only: nb_levls
#if (DOWNSCALING == 2)
       use input_subgrid2L, only: max_nb_points
#endif

       implicit none

! dmr 2016-01-22: moved here the last variables from input_downsts.f90

        REAL, DIMENSION(nlat,nlon) :: tland_max, tland_min, nethfxland_max, nethfxland_min

! dmr  [XXX] what is the reason of this strange 3:3 definition? Should disappear in new version ...
! dmr        my guess is that it is a "ntyps" variable that is shrunk to one element (land surface type)
! dmr        New version should include all indices for portability ...

        REAL, DIMENSION(nlat,nlon,3:3) :: efluxn_max, efluxn_min, hfluxn_max&
                                    ,hfluxn_min, qsurfn_max, qsurfn_min,&
                                    q10n_max, q10n_min


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Started addition of the new variables that depend on the vertical sampling
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

                                       ! Ugly trick: since I do not have access to the "ntyps" number
                                       ! of surface type (I do not want to include comsurf.h here)
                                       ! I thus declare a new local one temporarily until comsurf.h is
                                       ! rewritten as a FORTRAN90 module.
       integer, parameter :: ntyps_l = 3

       double precision, dimension(nlat,nlon,nb_levls)                :: q10_d, relhum_d, pground_d, tempsg_d, meltheat_d
       double precision, dimension(nlat,nlon,ntyps_l:ntyps_l,nb_levls):: efluxn_d, hfluxn_d, q10n_d, qsurfn_d
       double precision, dimension(nlat,nlon,nb_levls)                :: tland_d, nethfxland_d, rmount_d, rmount_ps_d, qmount_d
       double precision, dimension(nlat,nlon,nb_levls)                :: temp_profile, tsurf_d
       double precision, dimension(nlat,nlon,nb_levls)                :: dyrain_d, dysnow_d, corain_d, cosnow_d, torain_d, tosnow_d
       double precision, dimension(nlat,nlon,ntyps_l,nb_levls)        :: tsurfn_d
#if ( DOWNSCALING == 2 )
       !double precision, dimension(nlat,nlon,nb_levls,max_nb_points) :: tsurf_d_interp_sg
       double precision, dimension(nb_levls,nlat,nlon,max_nb_points) :: tsurf_d_interp_sg
#endif

       contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: create_surftemp_var_d
!
!>     @brief This subroutine create a complete surface temperature variable from tland and tocean with nb_levls vertical virtual
!!             surfaces for the vertical downscaling and the SMB computations
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine create_surftemp_var_d()

       use comland_mod, only:
#if ( COMATM == 1 )
       use comsurf_mod, only: nld, fractn
#endif

#if ( COMATM == 0 )
#include "comsurf.h"
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!       <void>

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer :: ind,lat,lon

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the subroutine starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       do lat = LBOUND(tempsg_d,dim=1), UBOUND(tempsg_d,dim=1)
        do lon = LBOUND(tempsg_d,dim=2), UBOUND(tempsg_d,dim=2)
         if (fractn(lat,lon,nld).LT.0.05) then                      ! mainly ocean
           do ind = LBOUND(tempsg_d,dim=3), UBOUND(tempsg_d,dim=3)
             temp_profile(lat,lon,ind) = tempsg_d(lat,lon,ind)
           enddo
         else                                                      ! mainly land
           do ind = LBOUND(tempsg_d,dim=3), UBOUND(tempsg_d,dim=3)
             temp_profile(lat,lon,ind) = tland_d(lat,lon,ind)
           enddo
         endif
        enddo
       enddo

      end subroutine create_surftemp_var_d

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       end module vertDownsc_mod
