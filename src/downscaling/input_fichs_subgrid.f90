!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module est une construction pour le couplage de GRISLI
!       initialement developpe pour definir les variables d'entrees 
!       du modele GRISLI Nord en provenance et vers le modele LOVECLIM 
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche 
!      Date   : 15 Juin 2008
!      Derniere modification : 21 Juillet 2008
!-----|--1--------2---------3---------4---------5---------6---------7-|
!afq -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Module de declaration externe en vue d'un couplage 
!            LOVECLIM --> GRISLI Nord
!-----|--1--------2---------3---------4---------5---------6---------7-|
module input_fichs_subgrid

  use taillesGrilles,  only: sgd

  implicit none


#if ( GRIDS_DOWN == 1)
  
!     Fichier externe contenant les lat, lon global
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_lonlat_subgrid=(/"inputdata/downscaling/coord_global.dat"/)
! lb unique ncfile for topo, area, IS mask, bias
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_downscal=(/"inputdata/downscaling/downs_global.nc"/)

#elif (GRIDS_DOWN == 2)

!     Fichier externe contenant les lat, lon HN for GRISLI Nord
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_lonlat_subgrid=(/"inputdata/downscaling/coord_NH.dat"/)
! lb unique ncfile for topo, area, IS mask, bias
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_downscal=(/"inputdata/downscaling/downs_NH.nc"/)

#elif (GRIDS_DOWN == 3)

!     Fichier externe contenant les lat, lon HN for GRISLI South
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_lonlat_subgrid=(/"inputdata/downscaling/coord_SH.dat"/)
! lb unique ncfile for topo, area, IS mask, bias
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_downscal=(/"inputdata/downscaling/downs_SH.nc"/)

#elif (GRIDS_DOWN == 4)

!     Fichier externe contenant les lat, lon HN for Europe
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_lonlat_subgrid=(/"inputdata/downscaling/coord_europe.dat"/)
! lb unique ncfile for topo, area, IS mask, bias
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_downscal=(/"inputdata/downscaling/downs_europe.nc"/)

! Uncomment to add a specific area
!#elif (GRIDS_DOWN == 5)
 
!     Fichier externe contenant les lat, lon HN for where you want
!  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
!       file_lonlat_subgrid=(/"inputdata/downscaling/??.dat"/)
! lb unique ncfile for topo, area, IS mask, bias
!  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
!       file_downscal=(/"inputdata/downscaling/??.nc"/)



#else
  
  write(*,*) "more than 4 subgrids not yet possible..."
  stop

#endif

end module input_fichs_subgrid
