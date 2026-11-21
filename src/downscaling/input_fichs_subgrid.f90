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


!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Module de declaration externe en vue d'un couplage 
!            LOVECLIM --> GRISLI Nord
!-----|--1--------2---------3---------4---------5---------6---------7-|
module input_fichs_subgrid

  use taillesGrilles, only: sgd

  implicit none

#define NB_SUBGRIDS 1

#if ( NB_SUBGRIDS == 1)
  
!     Fichier externe contenant les lat, lon de GRISLI Nord
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_lonlat_subgrid=(/"inputdata/downscaling/coord-nord-40km.dat"/)
! afq, pour l'instant pour la topo je laisse un truc en ascii, on mettre du ncdf plus tard
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_topo_subgrid=(/"inputdata/downscaling/S-nord-40km.dat"/)
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_epais_subgrid=(/"inputdata/downscaling/H-nord-40km.dat"/)
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_area_subgrid=(/"inputdata/downscaling/area-nord-40km.dat"/)
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_bias_subgrid=(/"inputdata/downscaling/bias-nord-40km.dat"/)

#elif ( NB_SUBGRIDS == 2 )
  
!     Fichier externe contenant les lat, lon de GRISLI Nord
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_lonlat_subgrid=(/"inputdata/downscaling/coord-nord-40km.dat","inputdata/downscaling/coord-sud-40km.dat"/)
! afq, pour l'instant pour la topo je laisse un truc en ascii, on mettre du ncdf plus tard
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_topo_subgrid=(/"inputdata/downscaling/S-nord-40km.dat","inputdata/downscaling/S-sud-40km.dat"/)
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_epais_subgrid=(/"inputdata/downscaling/H-nord-40km.dat","inputdata/downscaling/H-sud-40km.dat"/)
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_area_subgrid=(/"inputdata/downscaling/area-nord-40km.dat","inputdata/downscaling/area-sud-40km.dat"/)
  CHARACTER(LEN = 80), PARAMETER, DIMENSION (sgd) ::               &
       file_bias_subgrid=(/"inputdata/downscaling/bias-nord-40km.dat","inputdata/downscaling/bias-sud-40km.dat"/)

#else
  
  write(*,*) "more than 2 subgrids not yet possible..."
  stop

#endif

end module input_fichs_subgrid
