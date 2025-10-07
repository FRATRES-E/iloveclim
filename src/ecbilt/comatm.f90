!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Ce module est la conversion en module FORTRAN90/95 du fichier
!       comatm.h
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche (pour la conversion)
!      Date   : 18 Novembre 2008
!      Derniere modification : 14 Mai 2018, 17 juin 2022
!-----|--1--------2---------3---------4---------5---------6---------7-|
!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"
!dmr -- Adding the choice of components through the pre-processing options
!-----------------------------------------------------------------------
! *** File:     comatm.h
! *** Contents: General Parameter and Common declarations for ECbilt
!-----------------------------------------------------------------------
! *** PARAMETERS
!     nm  :   the truncation is of type T(riangular) nm.
!     nlon:   number of longitude points of the Gaussian grid
!     nlat:   number of latitude  points of the Gaussian grid
!     nvl :   number of vorticity levels in the vertical
!             (should be set to 3)
!     ntl :   number of temperature levels in the vertical
!             (equal to nvl-1)
!     nsh :   half of nsh2
!     nsh2:   number of coefficients needed to define one level of the
!             T nm model
!     ngp:    number of grid points of the Gaussian grid
!
! *** COMMON  /rvari/ pi,dp,om,rgas,grav,radius
!     pi :     value of pi
!     fzero:   value of f at 45 degrees
!     dp :     layer thicknesses [Pa]
!     om :     angular velocity of earth [rad/s]
!     rgas :   gas constant
!     grav :   gravity acceleration [m/s^2]
!     radius:  radius of earth [m]
!
! *** COMMON  /ipar/  iadyn,iaphys,initdate,initfield
!     iadyn :  with (1) atmospheric dynamics or without (0)
!     iaphys:  with (1) atmospheric physics or without (0)
!     iaphys:  with (1) atmospheric physics or without (0)
!     initfield:  initiqlisation form rest(0) or from file (1)
!     initdate: date of the beginning of the WHOLE experiment)
!               (actually it is the difference between the starting point
!                of the experiment and initdate)
!
! *** COMMON  /ggrid/ phi,dphi,dlab,darea,tarea,cosfi,cosfid,sinfi,
!                     sinfid,tanfi
!     phi :    Gauss points in radians
!     darea:   area [m**2] of gridboxes as a function of latitude
!     tarea:   total area of earth [m**2]
!     cosfi:   cosine of phi
!     sinfi:   sine of phi
!     tanfi:   tangens of phi
!
! *** COMMON  /ctstep/ dt,dtime,dtt
!     dt :     timestep in fraction of one day
!     dtime :  timestep in seconds
!     dtt :    timestep in non-dimensional units
!
! *** COMMON  /plev/  plevel,tlevel,rlogtl12
!     plevel:   contains value of the pressure at the nvl levels [Pa]
!     tlevel:   contains value of the pressure at the ntl levels [Pa]
!     rlogtl12: contains log(tlevel(1)/tlevel(2))
!-----------------------------------------------------------------------

      MODULE comatm

!-----|--1--------2---------3---------4---------5---------6---------7-|

      USE global_constants_mod, ONLY: dblp=>dp, ip

      IMPLICIT NONE

      INTEGER(ip), PARAMETER      :: nm=21, nlon=64, nlat=32, nvl=3, ntl=nvl-1
      INTEGER(ip), PARAMETER      :: nsh=((nm+1)*(nm+2))/2, nsh2=2*nsh
      INTEGER(ip), PARAMETER      :: ngp=nlon*nlat

      INTEGER(ip)                 :: iadyn,iaphys,initdate,initfield

      REAL(dblp)                  :: pi,fzero,dp,om,rgas,grav,radius
      REAL(dblp), DIMENSION(nlat) :: phi
      REAL(dblp), DIMENSION(nlat) :: cosfi,sinfi,tanfi
      REAL(dblp)                  ::  dt,dtime,dtt,rdtime
      REAL(dblp), DIMENSION(nvl)  :: plevel
      REAL(dblp), DIMENSION(ntl)  :: tlevel
      REAL(dblp)                  :: rlogtl12,p0,tarea,tareas
      REAL(dblp)                  :: alogtl12,alogtl1pl2,alogpl2tl2
      REAL(dblp), DIMENSION(nlat) :: darea,dareas
      
#if ( WISOATM == 0 )
      INTEGER(ip), parameter      :: nwisos = 1
      INTEGER(ip), parameter      :: iwater = 1
#elif ( WISOATM == 1 )
      INTEGER(ip), parameter      :: nwisos = 5
      INTEGER(ip), parameter      :: iwater = 1
      INTEGER(ip), parameter      :: iwat16 = 2
      INTEGER(ip), parameter      :: iwat17 = 3
      INTEGER(ip), parameter      :: iwat18 = 4
      INTEGER(ip), parameter      :: iwat2h = 5
#endif
      
      

!-----|--1--------2---------3---------4---------5---------6---------7-|

      END MODULE comatm
!      integer   nm,nsh,nlon,nlat,nsh2,ngp,nvl,ntl

!      common /ec_rvari/ pi,fzero,dp,om,rgas,grav,radius
!      common /ec_ipar/  iadyn,iaphys,initdate,initfield
!      common /ec_ggrid/ phi,cosfi,sinfi,tanfi,dareas,darea,tarea,tareas
!      common /ec_ctstep/ dt,dtime,dtt,rdtime
!      common /ec_plev/  plevel,tlevel,rlogtl12,p0,alogtl12,alogtl1pl2,
!     *               alogpl2tl2
