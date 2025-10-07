!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module a ete developpe comme interface des variables
!       sous-maille vers ECBilt
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche, Aurelien F. Quiquet
!      Date   : 4 Feb 2016
!      Derniere modification :
!-----|--1--------2---------3---------4---------5---------6---------7-|
!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"
!dmr -- Adding the choice of components through the pre-processing options
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Module de developpement en vue de l'extension du downscaling
!-----|--1--------2---------3---------4---------5---------6---------7-|

module input_subgrid2L

use taillesGrilles, only: iecb, jecb, sgd, sgnx, sgny, sgnxm, sgnym
#if ( ISM >= 2 )
use taillesGrilles, only: CNX,CNY
#endif

implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!       Ci-dessous on trouve un ensemble de declarations communes
!        pour les differentes sous-routines de l'interpolation
!        sous grille ==> LOVECLIM
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

!       Taille de la sous-grille:
!integer, parameter :: sgnx = 241, sgny = 241

!      Longitudes et latitudes ECBilt
real*8, dimension (iEcb) :: latEcb
real*8, dimension (jEcb) :: lonEcb
!      Tableaux des lat, lon de la sous grille consideree
real, dimension(sgnxm,sgnym,sgd) :: latSG, lonSG
!      Tableaux des distances au pole sur la sous grille
real, dimension(sgnxm,sgnym,sgd) :: distlat, distlon
!      Tableaux des lat, lon d'ECBilt sur la sous grille consideree
real, dimension (sgnxm,sgnym,sgd)  :: iEg, jEg !AFQ: vient de input_L2GRISLI, normalement en target... a voir
!!real, dimension (sgnx,sgny)  :: iEgr, jEgr !AFQ: vient de input_L2GRISLI, normalement en target... a voir

!       Tableau du nombre de cases sous grille sur la grille ECBilt
integer,          dimension(iecb,jEcb)          :: nbpointssg
double precision, dimension(iEcb,jEcb)          :: sub_grid_notflat
double precision, dimension(iEcb,jEcb)          :: area_sg_onecb

!       Tableaux de la topo sur la sous-grille
real (kind=8), dimension(sgnxm,sgnym,sgd)      :: topoSG
!       Tableaux de l'epaisseur sur la sous-grille
real (kind=8), dimension(sgnxm,sgnym,sgd)      :: epaisSG
!       Tableaux du masque glace dans la sous-grille
real (kind=8), dimension(sgnxm,sgnym,sgd)      :: masqueSGice
!       Tableaux de la surface des pixels de la sous-grille
real (kind=8), dimension(sgnxm,sgnym,sgd)      :: areaSG
#if ( SMB_TYP == 2 )
!       Tableaux du biais annuel en temperature sur la sous-grille
real (kind=8), dimension(sgnxm,sgnym,sgd)      :: annbiasSG
#endif
!      Definition des tableaux de variables climatiques interpolees
!        en moyenne annuelle
real (kind=4), dimension(sgnxm,sgnym,sgd)          :: tfyear
real (kind=4), dimension(sgnxm,sgnym,sgd)          :: pfyear
#if ( SMB_TYP >= 1 )
real (kind=4), dimension(sgnxm,sgnym,sgd)          :: SMB_iLCM
#endif
#if ( ISM >= 2 )
real (kind=4), dimension(:,:), pointer :: tfyearnord
real (kind=4), dimension(:,:), pointer :: pfyearnord
real (kind=4), pointer                 :: nivo_mer
real (kind=4), dimension(:,:), pointer :: SMB_iLCMnord
#else
real (kind=4), parameter                     :: nivo_mer = 0.
#endif
real (kind=4), dimension(sgnxm,sgnym,sgd)          :: relhumyear
#if ( ISM >= 2 )
real (kind=8), dimension(:,:), pointer :: topoGRIS_sansnivomernord
real (kind=8), dimension(sgnxm,sgnym,sgd)      :: topoGRIS
real (kind=8), dimension(sgnxm,sgnym,sgd)      :: masqueGRIS

real (kind=8), dimension(:,:), pointer :: epaisglaceGRISnord


!      Definition des tableaux contenant les latitudes, longitudes sur
!      des points de la grille CLIO
real, dimension(CNX, CNY) :: latClio, lonClio
!      Definition des tableaux contenant les indices de correspondance
!      de la grille CLIO sur la grille ISM
integer, dimension(sgnxm, sgnym,sgd) :: iCliog, jCliog
#if ( SHELFMELT == 1 )
! basal melting under ice shelves
real (kind=8), dimension(:,:,:), pointer :: bmshelf_clionord
#endif

#endif

! Declaration des variables necessaires a l'interpolation dgc-waffles:
integer, parameter :: nneigh = 4 ! nb of neighbours for interpolation: 4 9 25 49 81
integer, parameter :: nexted = 2  ! nb of points for domain extension

integer     , dimension (2,0:nneigh,sgnxm,sgnym,sgd) :: index_interpL2G
real(kind=8), dimension (    nneigh,sgnxm,sgnym,sgd) :: weights_interpL2G
real(kind=8), dimension (           sgnxm,sgnym,sgd) :: sumweights_interpL2G


#if ( DOWNSCALING == 2 || ISM > 1 )
integer, parameter                          :: max_nb_points = 250 ! should be 217 for GRISLI, Hemin40 with ECBilt @ T21, afq: or not?
integer, dimension(iEcb,jEcb,max_nb_points) :: nx_list_points_sg, ny_list_points_sg, ng_list_points_sg
double precision, dimension(sgnxm,sgnym,sgd)      :: not_flat_subgrid

! Le tableau suivant donne les coordonnées du point subgrid i,j dans la grille ECBilt
!   sous la forme (nlat,nlon,nbpointssg)
integer,          dimension(3,sgnxm,sgnym,sgd)                 :: coords_subgrid_in_ECB
integer,          dimension(2,max_nb_points,iEcb,jECb)   :: coords_ECB_in_subgrid

integer,dimension(sgd) :: e_lat_low, e_lat_hig, e_lon_low, e_lon_hig

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr&aq  Subgrid downscaling variables ...
! dmr     dynamic rain & snow on the subgrid ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

double precision, dimension(iEcb,jEcb,max_nb_points) :: dyrain_sg, dysnow_sg, corain_sg, cosnow_sg, torain_sg, tosnow_sg, tsurf_sg&
                                                      , weights_low_sg, area_sg, relhum_sg, topo_sg
#if ( DOWN_T2M == 1 )
double precision, dimension(iEcb,jEcb,max_nb_points) :: tempsg_sg
#endif

double precision, dimension(iEcb,jEcb,max_nb_points) :: difftopo_sg, epaisglace_sg, snow_age_sg
#if ( SMB_TYP == 2 )
double precision, dimension(iEcb,jEcb,max_nb_points) :: annbias_sg
#endif
integer, dimension(iEcb,jEcb,max_nb_points)          :: index_low_sg

#if ( CARAIB > 0 || CARAIB_FORC_W > 0 )
double precision, dimension(iEcb,jEcb,max_nb_points) :: lon_sg, lat_sg
#endif

#endif

#if ( F_PALAEO == 3 )

integer, parameter :: nb_steps_fism    = 161 ! number of snapshots in the prescribed topo
integer, parameter :: update_time_fism = 250 ! in yrs
integer(kind=4)    :: indx_fISM
! Ganopolski files on NH40 contains 40 kyears of data with a 250 yrs step.
! files start at 40 kyrs BP and stops at 0 yrs BP (161 snapshots)

double precision,dimension(sgnx,sgny,nb_steps_fism) :: forcedISM !PB_PALAEO!!!!!! TODOOOOOOOOOOOOOOOOOOOOOOOOOOOO
double precision,dimension(sgnx,sgny,nb_steps_fism) :: forcedMSK
#if ( F_PALAEO_FWF == 1 )
double precision,dimension(sgnx,sgny,nb_steps_fism) :: forcedTHI  ! in meters
double precision,dimension(sgnx,sgny)               :: fluxFWF_SG ! in m^3.s^-1
double precision,dimension(sgnx,sgny)               :: fluxFWF_SG_route ! in m^3.s^-1
double precision,dimension(sgnx,sgny)               :: surface_FWF ! in m^2
#endif

#endif

contains

#if ( DOWNSCALING == 2 )
       subroutine cluster_subgrid_in_ecbilt()

!$ USE OMP_LIB

        integer i,j,n

        if (MAXVAL(nbpointssg(:,:)) > max_nb_points ) then
          write(*,*) "need to raise the number of points parameter to ", MAXVAL(nbpointssg(:,:))
          stop
        endif



        do n=1,sgd
!$OMP PARALLEL
!$OMP DO PRIVATE (i,j)
           do j=1,sgny(n)
              do i=1,sgnx(n)

               nx_list_points_sg(coords_subgrid_in_ECB(1,i,j,n),coords_subgrid_in_ECB(2,i,j,n),coords_subgrid_in_ECB(3,i,j,n)) = i
               ny_list_points_sg(coords_subgrid_in_ECB(1,i,j,n),coords_subgrid_in_ECB(2,i,j,n),coords_subgrid_in_ECB(3,i,j,n)) = j
               ng_list_points_sg(coords_subgrid_in_ECB(1,i,j,n),coords_subgrid_in_ECB(2,i,j,n),coords_subgrid_in_ECB(3,i,j,n)) = n

            enddo
         enddo
!$OMP END DO
!$OMP END PARALLEL
      enddo


        return
       end subroutine cluster_subgrid_in_ecbilt

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine cluster_var_subgrid_in_ecbilt_dble(var_subgrid,var_ecbilt_sg)

!$ USE OMP_LIB

        double precision, dimension(sgnxm,sgnym,sgd), intent(in) :: var_subgrid
        double precision, dimension(iEcb,jEcb,max_nb_points), intent(out) :: var_ecbilt_sg

        integer i,j,n

        do n=1,sgd
!$OMP PARALLEL
!$OMP DO PRIVATE (i,j)
           do j=1,sgny(n)
              do i=1,sgnx(n)
                 var_ecbilt_sg(coords_subgrid_in_ECB(1,i,j,n),coords_subgrid_in_ECB(2,i,j,n),coords_subgrid_in_ECB(3,i,j,n)) &
                               =var_subgrid(i,j,n)
              enddo
           enddo
!$OMP END DO
!$OMP END PARALLEL
        enddo

        return
       end subroutine cluster_var_subgrid_in_ecbilt_dble

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       subroutine cluster_nlev_subgrid_in_ecbilt_dble(var_subgrid,var_ecbilt_sg, nlevs)

!$ USE OMP_LIB

        integer, intent(in) :: nlevs
        double precision, dimension(nlevs,sgnxm,sgnym,sgd), intent(in) :: var_subgrid
        double precision, dimension(nlevs,iEcb,jEcb,max_nb_points), intent(out) :: var_ecbilt_sg

        integer i,j,n

        do n=1,sgd
!$OMP PARALLEL
!$OMP DO PRIVATE (i,j)
           do j=1,sgny(n)
              do i=1,sgnx(n)
                 var_ecbilt_sg(:,coords_subgrid_in_ECB(1,i,j,n),coords_subgrid_in_ECB(2,i,j,n),coords_subgrid_in_ECB(3,i,j,n)) &
                         =var_subgrid(:,i,j,n)
              enddo
           enddo
!$OMP END DO
!$OMP END PARALLEL
        enddo

        return
       end subroutine cluster_nlev_subgrid_in_ecbilt_dble


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       subroutine cluster_var_subgrid_in_ecbilt_dbleT(var_subgrid,var_ecbilt_sg)

        !$ USE OMP_LIB

        double precision, dimension(sgnxm,sgnym,sgd), intent(in) :: var_subgrid
        double precision, dimension(max_nb_points,iEcb,jEcb), intent(out) :: var_ecbilt_sg

        integer i,j,n

        do n=1,sgd
!$OMP PARALLEL
!$OMP DO COLLAPSE(2)
           do j=1,sgny(n)
              do i=1,sgnx(n)
                 var_ecbilt_sg(coords_subgrid_in_ECB(3,i,j,n),coords_subgrid_in_ECB(1,i,j,n),coords_subgrid_in_ECB(2,i,j,n)) &
                         =var_subgrid(i,j,n)
              enddo
           enddo
!$OMP END DO
!$OMP END PARALLEL
        enddo
     
        return
       end subroutine cluster_var_subgrid_in_ecbilt_dbleT

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


       subroutine cluster_var_subgrid_in_ecbilt_intg(var_subgrid,var_ecbilt_sg)

!$ USE OMP_LIB

        integer, dimension(sgnxm,sgnym,sgd), intent(in) :: var_subgrid
        integer, dimension(iEcb,jEcb,max_nb_points), intent(out) :: var_ecbilt_sg

        integer i,j,n

        do n=1,sgd
!$OMP PARALLEL
!$OMP DO PRIVATE (i,j)
           do j=1,sgny(n)
              do i=1,sgnx(n)
                 var_ecbilt_sg(coords_subgrid_in_ECB(1,i,j,n),coords_subgrid_in_ECB(2,i,j,n),coords_subgrid_in_ECB(3,i,j,n)) &
                         =var_subgrid(i,j,n)
              enddo
           enddo
!$OMP END DO
!$OMP END PARALLEL
        enddo
        
        return
       end subroutine cluster_var_subgrid_in_ecbilt_intg

       subroutine cast_ecb_subgrid_in_2d(var_ecbsg_in,var_subgrid_out)

!$ USE OMP_LIB

        double precision, dimension(iEcb,jEcb,max_nb_points), intent(in) :: var_ecbsg_in
        double precision, dimension(sgnxm,sgnym,sgd), intent(out) :: var_subgrid_out

        integer i,j,k,i_gris,j_gris,n_gris

!$OMP PARALLEL
!$OMP DO PRIVATE (i_gris,j_gris,n_gris)
        do i=1,iEcb
          do j=1,jEcb

            if (nbpointssg(i,j).gt.0) then !afq
              do k = 1, nbpointssg(i,j)

                i_gris = nx_list_points_sg(i,j,k)
                j_gris = ny_list_points_sg(i,j,k)
                n_gris = ng_list_points_sg(i,j,k)

                var_subgrid_out(i_gris,j_gris,n_gris) = var_ecbsg_in(i,j,k)

              enddo
            endif

          enddo
        enddo
!$OMP END DO
!$OMP END PARALLEL

        return
       end subroutine cast_ecb_subgrid_in_2d

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: reset_rain_snow_sub_grid
!
!>     @brief This subroutine resets rain / snow variables in the downscaling for sub-grid
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subroutine reset_rain_snow_sub_grid()

         relhum_sg(:,:,:) = 0.0d0
         dyrain_sg(:,:,:) = 0.0d0
         dysnow_sg(:,:,:) = 0.0d0
         corain_sg(:,:,:) = 0.0d0
         cosnow_sg(:,:,:) = 0.0d0
         torain_sg(:,:,:) = 0.0d0
         tosnow_sg(:,:,:) = 0.0d0

        return
       end subroutine reset_rain_snow_sub_grid
#endif
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module input_subgrid2L
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
