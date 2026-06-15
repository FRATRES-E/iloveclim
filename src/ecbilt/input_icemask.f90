!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module contient les FLAGS en lien avec GRISLI
!       couplage LOVECLIM --> GRISLI
!       (dans l'environnement logiciel LUDUS)
!   
!      Auteur : Didier M. Roche 
!      Date   : 18/Juin/2007 <-- Mask of Precip, Lev Tarasov forcing !
!      Derniere modification : 18 Novembre 2008, Conversion en module
!-----|--1--------2---------3---------4---------5---------6---------7-|

      MODULE input_icemask
      USE comatm, only: nlat, nlon
! ### dmr ### ci-dessous : sale, mais problemes de compatibilite sinon
      REAL, DIMENSION(nlat,nlon), TARGET :: icemask
!###      REAL :: totalprecipcases
! afq -- we duplicate the icemask, for albedo purposes:
!        this is used by ec_co2oc to force a snow albedo is there is
!        an ice sheet over the ocean
      real, dimension(nlat,nlon) :: icemaskalb

      END MODULE input_icemask
