!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce sous-programme initialise les variables de la sous-grille
!       dans la version du calcul d'interpolation LOVECLIM - sous-grille
!       provient essentiellement de Init_L2GRISLI (rev 953)
!
!      Auteur : Didier M. Roche, Aurelien F. Quiquet
!      Date   : 5 Feb 2019
!      Derniere modification : -
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!afq -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"
!afq -- Adding the choice of components through the pre-processing options

       subroutine Init_subgrid

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Variables d'entree  :
!       Variables de sortie :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       use taillesGrilles,  only: sgd,sgnx,sgny,iEcb,jEcb
       use input_subgrid2L, only: latEcb,lonEcb,latSG,lonSG,iEg,jEg,distlat,distlon,nbpointssg,topoSG,epaisSG,areaSG
#if ( DOWNSCALING == 2 || ISM >= 2 )
       use input_subgrid2L, only: max_nb_points,                                                                                  &
            coords_subgrid_in_ECB, coords_ECB_in_subgrid, e_lat_low, e_lat_hig, e_lon_low, e_lon_hig,                             &
            nneigh, nexted, index_interpL2G,weights_interpL2G, sumweights_interpL2G
#endif
       use input_fichs_subgrid, only: file_lonlat_subgrid, file_topo_subgrid, file_epais_subgrid, file_area_subgrid,              &
            file_bias_subgrid
       use comdyn, only: rmount
       use output_ecbilt, only: topoECB

#if ( ISM >= 2 )
       use input_subgrid2L, only: CNX,CNY,latClio,lonClio,icliog,jcliog
#if ( SHELFMELT == 1 )
       use input_subgrid2L, only: bmshelf_CLIOnord
#endif
#endif

#if ( DOWNSCALING == 2 )
       use input_subgrid2L, only: cluster_subgrid_in_ecbilt,sub_grid_notflat
#if ( SMB_TYP == 2)
       use input_subgrid2L, only: annbiasSG
#endif
#endif
!afq -- outdated 11/2022       use interpol_ecb_subgrid_mod, only: init_interpol 
       use interpolate_mod, only: interpolate_init
       use transfer_ecb_subgrid_mod, only: init_subgrid_writing
#if ( CARAIB > 0 || CARAIB_FORC_W > 0 )
       use input_subgrid2L, only: cluster_var_subgrid_in_ecbilt_dble
       use input_subgrid2L, only: lat_sg, lon_sg, max_nb_points
#endif

       implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Variables locales
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical :: success

       integer i,j,nbd,ios
       integer :: Topo_id,IceThickness_id,Area_id,Bias_id,CoordT21_id
       character(len=80) :: filin

#if ( DOWNSCALING == 2 )
       integer,dimension(iEcb,jEcb) :: nbpointssg_loc
#endif
#if ( CARAIB > 0 || CARAIB_FORC_W > 0 )
       double precision, dimension(iEcb,jEcb,max_nb_points) :: lon_sgloc,lat_sgloc
#endif

!--- #if ( DOWNSTS == 1 )
!---       character(len=256) :: file_nm
!--- #endif

#if ( ISM >= 2 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Associe les pointeurs des variables L2G aux variables cibles
!        dans GRISLI : a terme le module input_GRISLI doit etre
!        remplace par un autre coherent avec le vrai modele GRISLI
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       call faitPointer
#endif

! afq -- we need to init the topo outside the subgrid:
       topoECB(:,:) = rmount(:,:)
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Lecture du fichier de longitudes, latitudes de la subgrid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       do nbd=1,sgd ! loop on the different subgrids

! afq TODO pour l'instant on laisse tout en ASCII avec les trucs de GRISLI, a mettre en netcdf un jour
! note : distlat, distlon ne servent à rien
          call readGrilleGRISLI(sgnx(nbd),sgny(nbd),distlat(1:sgnx(nbd),1:sgny(nbd),nbd),distlon(1:sgnx(nbd),1:sgny(nbd),nbd), &
               latSG(1:sgnx(nbd),1:sgny(nbd),nbd),lonSG(1:sgnx(nbd),1:sgny(nbd),nbd),file_lonlat_subgrid(nbd),10)

! lecture du fichier topo           
          filin=file_topo_subgrid(nbd)
          open(newunit=Topo_id,file=filin,iostat=ios)
          do j=1,sgny(nbd)
             do i=1, sgnx(nbd)
                read(Topo_id,*) topoSG(i,j,nbd)
             enddo
          enddo
          close(Topo_id)
          
! lecture du fichier epaisseur glace
          filin=file_epais_subgrid(nbd)
          open(newunit=IceThickness_id,file=filin,iostat=ios)
          do j=1,sgny(nbd)
             do i=1, sgnx(nbd)
                read(IceThickness_id,*) epaisSG(i,j,nbd)
             enddo
          enddo
          close(IceThickness_id)
              
! lecture du fichier aire des pixels
#if ( ISM >= 2 )
          areaSG(:,:,nbd) = 1.
#else
          filin=file_area_subgrid(nbd)
          open(newunit=Area_id,file=filin,iostat=ios)
          do j=1,sgny(nbd)
             do i=1, sgnx(nbd)
                read(Area_id,*) areaSG(i,j,nbd)
             enddo
          enddo
          close(Area_id)
#endif
          
! lecture du fichier de biais

#if ( SMB_TYP == 2 )
          filin=file_bias_subgrid(nbd)
          open(newunit=Bias_id,file=filin,iostat=ios)
          do j=1,sgny(nbd)
             do i=1, sgnx(nbd)
                read(Bias_id,*) annbiasSG(i,j,nbd)
             enddo
          enddo
          close(Bias_id)
#endif
       
       enddo

       topoSG(:,:,:) = max(topoSG(:,:,:),0.)


#if ( ISM >= 2 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Lecture du fichier de longitudes, latitudes de CLIO
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       call readLatLonClio(CNX,CNY,latClio,lonClio,                     &
           "inputdata/lcm2ism/lat.dat", "inputdata/lcm2ism/lon.dat")
#endif

       nbpointssg(:,:) = 0

#if ( CARAIB > 0 || CARAIB_FORC_W > 0 )
       lon_sgloc(:,:,:) = 0.
       lat_sgloc(:,:,:) = 0.
#endif
       
       do nbd=1,sgd
          
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Calcul des tableaux des longitudes, latitudes d'ECBilt sur
!        la sous grille consideree
!       Attention cette routine suppose que la grille d'ECBilt est en
!        5.625 x 5.625 (ecrit en dur dans la routine)
!       Puis calcul du nombre de points sous grille sur chaque case ECBilt
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          call calcijEcbg(sgnx(nbd),sgny(nbd),                                                                             &
                          latSG(1:sgnx(nbd),1:sgny(nbd),nbd),lonSG(1:sgnx(nbd),1:sgny(nbd),nbd),                           &
                          iEg(1:sgnx(nbd),1:sgny(nbd),nbd),jEg(1:sgnx(nbd),1:sgny(nbd),nbd))
       
#if ( DOWNSCALING == 2 || ISM >= 2 ) 
          call nbPointISM_Ecb(sgnx(nbd),sgny(nbd),iEcb,jEcb,                                                                  &
                              jEg(1:sgnx(nbd),1:sgny(nbd),nbd),iEg(1:sgnx(nbd),1:sgny(nbd),nbd),nbpointssg_loc,max_nb_points, &
                              coords_subgrid_in_ECB(:,1:sgnx(nbd),1:sgny(nbd),nbd),                                           &
                              coords_ECB_in_subgrid,                                                                          &
                              e_lat_low(nbd), e_lat_hig(nbd), e_lon_low(nbd), e_lon_hig(nbd))
          nbpointssg(:,:)=nbpointssg(:,:)+nbpointssg_loc(:,:)
#endif

       enddo !on the number of subgrids

!      Ci-dessous pour verification du calcul du nombre de points d .
!      CALL output_grd2D(nbpointssg,iEcb,jEcb,'checknbPoints.dat',1)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Init des poids pour les calculs de transferts de grille
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       filin='inputdata/coord-T21.dat'
! lecture des coordonnees d'ECBilt dans un fichier externe
       open(newunit=CoordT21_id,file=filin,iostat=ios)
       read(CoordT21_id,*) (latEcb(i),i=iEcb,1,-1)
       read(CoordT21_id,*) (lonEcb(j),j=1,jEcb, 1)
       close(CoordT21_id)

!afq -- outdated 11/2022        success = init_interpol()
#if ( DOWNSCALING == 2 )
       do nbd=1,sgd

          success = interpolate_init(sgnx(nbd),sgny(nbd),iEcb,jEcb,latEcb,lonEcb,                                     &
               latSG(1:sgnx(nbd),1:sgny(nbd),nbd),lonSG(1:sgnx(nbd),1:sgny(nbd),nbd),3,nneigh,nexted,                 &
               index_interpL2G(:,:,1:sgnx(nbd),1:sgny(nbd),nbd),                                                      &
               weights_interpL2G(:,1:sgnx(nbd),1:sgny(nbd),nbd), sumweights_interpL2G(1:sgnx(nbd),1:sgny(nbd),nbd),   &
               sub_grid_notflat)

       enddo
#endif

#if ( DOWNSCALING == 2 )
       call cluster_subgrid_in_ecbilt()
#endif

#if ( ISM >= 2 )
        
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Calcul des tableaux des longitudes, latitudes de CLIO sur
!        la grille de l'ISM considere
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       do nbd=1,sgd
          call calcijCliog(sgnx(nbd),sgny(nbd),latSG(:,:,nbd),lonSG(:,:,nbd),latClio,lonClio,iCliog(:,:,nbd), &
            jCliog(:,:,nbd))
       enddo
       
#endif /* On ISM >= 2 */

!--- #if ( DOWNSTS == 1 )
!---        file_nm = "outputdata/downscaling/sgout_profil.nc"
!---        call write_nc2d_subgrid_init(file_nm)
!--- #if ( DOWNSCALING == 2 )
!---        file_nm = "outputdata/downscaling/sgout_subgrid.nc"
!---        call write_nc2d_subgrid_init(file_nm)
!--- #endif
!--- #endif

!afq -- may 2020, we do no longer store the "profile" outputs, only sg:
#if ( DOWNSCALING == 2 )
       success = init_subgrid_writing()
#endif

!afq -- for CARAIB outputs, we need the lon/lat on the clustered grid:
#if ( CARAIB > 0 || CARAIB_FORC_W > 0 )
       do nbd=1,sgd
          call cluster_var_subgrid_in_ecbilt_dble(lonSG(:,:,nbd),lon_sgloc)
          call cluster_var_subgrid_in_ecbilt_dble(latSG(:,:,nbd),lat_sgloc)
          lon_sg(:,:,:)=lon_sg(:,:,:)+lon_sgloc(:,:,:) ! the grids are disconnected
          lat_sg(:,:,:)=lat_sg(:,:,:)+lat_sgloc(:,:,:) ! the grids are disconnected
      enddo
#endif

#if ( SHELFMELT == 1 )
       bmshelf_clionord(:,:,:) = 0.
#endif
       end subroutine Init_subgrid
