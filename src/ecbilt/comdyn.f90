!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module est la conversion en module FORTRAN90/95 du fichier
!       comatm.h
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche (pour la conversion)
!      Date   : 14 Mai 2018
!      Derniere modification : 14 Mai 2018, 21 Juin 2022
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** File:     comdyn.f90
! *** Contents: Common declarations for dynamical part of atmospheric model of ECbilt
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! *** COMMON /ec_intpar/ nshm, ll
!     nshm:   contains numbers 22 down to 1 for index 0 to 21
!     ll:     contains total wavenumber n of each spherical harmonic of
!             the corresponding index
!
! *** COMMON /ec_linopr/ rm, rinhel, diss
!     rm:     contains zonal wavenumber m of each spherical harmonic of
!             the corresponding index for zonal derivative operator
!     rinhel: Laplace and Helmholtz operator for Q-PSI inversion
!     diss:   dissipation coefficients for each spherical harmonic
!             diss(k,1) : hyperviscosity at the three levels
!                         (tdif sets timescale)
!             diss(k,2) : Ekman friction at lower level
!                         (tdis sets timescale)
!
! *** COMMON /ec_logpar/ lgdiss,infcomdyn.h
!     lgdiss: if .true. then orography and land-sea mask dependent
!             friction at the lower level plus Ekman friction,
!             else only Ekman friction
!     inf:    if .true. then artificial PV forcing read from file
!
! *** COMMON /ec_metras/
!     pp:     Legendre polynomials defined at Gausian latitudes
!     pd:     mu derivative of Legendre polynomials
!     pw:     weights for Legendre integrals
!
! *** COMMON /ec_phypar/ rdiss, ddisdx, ddisdy
!     rdiss:  landsea-mask/orography dependent friction
!     ddisdx: landsea-mask/orography dependent friction
!     ddisdy: landsea-mask/orography dependent friction
!
! *** COMMON /ec_modpar/ rrdef1,rrdef2,rl1,rl2,relt1,relt2 ,
!                      tdis,trel,tdif,addisl,addish,h0,idif
!     rrdef1: Rossby radius of deformation of 200-500 thickness
!     rrdef2: Rossby radius of deformation of 500-800 thickness
!     rl1:    one over Rossby rad. of def. squared of 200-500 thickness
!     rl2:    one over Rossby rad. of def. squared of 500-800 thickness
!     relt1:  nondimensional relaxation coefficient of 200-500 thickness
!     relt2:  nondimensional relaxation coefficient of 500-800 thickness
!     tdis:   Ekman dissipation timescale in days at lower level
!     trel:   relaxation time scale in days of the temperature
!     tdif:   dissipation timescale of scale-selective diffusion in
!             days for wavenumber nm
!     addisl: parameter used in the computation of the dissipation
!             timescale at the lower level over land
!     addish: parameter used in the computation of the dissipation
!             timescale at the lower level as a function of topography
!     h0:     scale factor for the topographically induced upward motion
!             at the lower level
!     idif:   determines scale-selectivity of hyperviscosity; power of
!             laplace operator
!
! *** COMMON /ec_sfield/ psi, psit, psito, qprime, dqprdt, for, ws
!     psi:    stream function at the nvl levels
!     psit:   thickness at the ntl levels
!     psito:  thickness at the ntl levels at the previous timestep
!     qprime: potential vorticity
!     dqprdt: time derivative of qprime
!     for:    constant potential vorticity forcing at the nvl levels
!     ws:     only used as portable workspace
!
! *** COMMON /ec_corog/ orog, rmount, rmountn
!     orog:   orography in m. divided by h0
!     rmount: mean orography in m.
!     rmountn orography in m for each surface type
!     dorodl: derivative of orog wrt lambda
!     dorodm: derivative of orag wrt sin(fi)
!
! *** COMMON /ec_zotras/ trigd, trigi, wgg
!             arrays used by the nag version of the fft
!
! *** COMMON /ec_cwind/ uv10,utot,vtot,u800,v800,u500,v500,u200,v200
!     uv10:   wind strength at 10 m, extrapolated from 800 hPa at grid
!             with imposed lower bound of uv10m m/s
!     uvw10:  wind strength at 10 m, extrapolated from 800 hPa at grid
!     uv10r:  reduction factor between 800 hPa and surface winds
!     uv10m:  minimum value of 10 m wind
!     utot :   total zonal wind
!     vtot :   total meridional wind
!     u800:   800 hPa geostrophic zonal wind
!     v800:   800 hPa geostrophic meridional wind
!     u500:   500 hPa geostrophic zonal wind
!     v500:   500 hPa geostrophic meridional wind
!     u200:   200 hPa geostrophic zonal wind
!     v200:   200 hPa geostrophic meridional wind
!
! *** COMMON /ec_cgpsi/  grpsi1,grpsi2,grpsi3
!     grpsi1: streamfunction on Gaussian grid at 200 mb
!     grpsi2: streamfunction on Gaussian grid at 500 mb
!     grpsi3: streamfunction on Gaussian grid at 800 mb
!
! *** COMMON /ec_cdwin/  udiv,vdiv,udivg,vdivg
!     udiv:   divergent wind east-west direction in spectral form
!     vdiv:   divergent wind north-south direction in spectral form
!     udivg:  divergent wind east-west direction on Gaussian grid
!     vdivg:  divergent wind north-south direction on Gaussian grid
!
! *** COMMON /ec_cdiv/   divs,divg
!     divs:   divergence in spectral form
!     divg:   divergence on Gaussian grid
!
! *** COMMON /ec_forcx/  iartif,ipvf1,ipvf2,ipvf3,ipvf4,ipvf5
! forcing parameter:
! iartif:     with (1) or without (0) artificial forcing
! ipvf1 :     with (1) or without (0) diabatic heating
! ipvf2 :     with (1) or without (0) advection of f by divergent wind
! ipvf3 :     with (1) or without (0) stretching term
! ipvf4 :     with (1) or without (0) advection of zeta by divergent
!             wind
! ipvf5 :     with (1) or without (0) advection of temperature by
!                                  divergent wind
! *** COMMON /ec_clfor/  forcgg1,forcggw1,forcggs1
!     forcgg1: daily variying artificial PV forcing on Gaussian grid
!              at 200 hPa only calculated from forcggw1 and forcggs1
!     forcggw1: artificial PV forcing winter season on Gaussian grid
!               at 200 hPa
!     forcggs1: artificial PV forcing summer season on Gaussian grid
!               at 200 hPa
!
! *** COMMON /ec_cvert/  omegs,omegg
!     omegs: vertical velocity omega in spectral form
!     omegg: vertical velocity omega on Gaussian grid
!
! *** COMMON /ec_cdfor/  dfor1,dfor2
!     dfor1: PV forcing of 200-500 thickness due to diabatical terms
!     dfor2: PV forcing of 500-800 thickness due to diabatical terms
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      MODULE comdyn

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use global_constants_mod, only: dblp=>dp, ip
      use comatm, only: nlat,nlon,nm,nsh,nsh2,nvl,ntl
      
      IMPLICIT NONE

      ! [NOTA] the following should be private, but are taken by other
      ! modules already. Keeping without for now 2022-06-21

      ! PRIVATE :: dblp, ip, nlat, nlon, nm, nsh, nsh2, nvl, ntl


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer(ip), dimension(0:nm)     :: nshm
      integer(ip), dimension(nsh)      :: ll

      integer(ip) :: ipert
      integer(ip) :: idif

      real(dblp),    dimension(nsh)      :: rm
      real(dblp),    dimension(nsh2,0:5) :: rinhel
      real(dblp),    dimension(nsh2,2)   :: diss
      real(dblp),    dimension(nlat,nsh) :: pp,pd,pw
      real(dblp),    dimension(nlat,nlon):: rdiss,ddisdx,ddisdy

      real(dblp) ::  rrdef1,rrdef2,rl1,rl2,relt1,relt2,tdis,trel,tdif
      real(dblp) ::  addisl,addish,h0

      real(dblp),    dimension(nsh2,nvl)     :: psi,qprime,dqprdt,for
      real(dblp),    dimension(nsh2,ntl)     :: psit,psito
      real(dblp),    dimension(nsh2)         :: ws, orog
      real(dblp),    dimension(nlat,nlon)    :: rmount, dorodl, dorodm, wgg
      real(dblp),    dimension(nlon,2)       :: trigd, trigi
      real(dblp),    dimension(nlat,nlon,3)  :: utot, vtot
      real(dblp),    dimension(nlat,nlon)    :: u800,v800,u500,v500,u200,v200
      real(dblp),    dimension(nsh2,nvl)     :: udiv,vdiv,chi
      real(dblp),    dimension(nlat,nlon,nvl):: udivg,vdivg,chig,psig,qgpv,geopg
      real(dblp),    dimension(nsh2,nvl)     :: divs,omegs
      real(dblp),    dimension(nlat,nlon,nvl):: divg,omegg

      integer(ip) :: iartif,ipvf1,ipvf2,ipvf3,ipvf4,ipvf5

      real(dblp),    dimension(nlat,nlon)    :: forcgg1,forcggw1,forcggs1,gekdis
      real(dblp),    dimension(nsh2)         :: dfor1 = 0.0_dblp, dfor2=0.0_dblp

      logical     :: lgdiss,inf


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      END MODULE comdyn
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~       common /ec_intpar/ nshm, ll, ipert
!~       common /ec_linopr/ rm, rinhel, diss
!~       common /ec_logpar/ lgdiss,inf
!~       common /ec_metras/ pp, pd, pw
!~       common /ec_phypar/ rdiss, ddisdx, ddisdy
!~       common /ec_modpar/ rrdef1,rrdef2,rl1,rl2,relt1,relt2 ,                                                                      &
!~      &                tdis,trel,tdif,addisl,addish,h0,idif
!~       common /ec_sfield/ psi, psit, psito, qprime, dqprdt, for, ws
!~       common /ec_corog/  orog, rmount , dorodl , dorodm
!~       common /ec_zotras/ trigd, trigi, wgg
!~       common /ec_cwind/  utot, vtot, u800, v800, u500,v500,u200,v200
!~       common /ec_cgpsi/  psig,qgpv,geopg
!~       common /ec_cdwin/  udiv,vdiv,udivg,vdivg
!~       common /ec_cdiv/   divs,divg,chig,chi
!~       common /ec_forcx/  iartif,ipvf1,ipvf2,ipvf3,ipvf4,ipvf5
!~       common /ec_clfor/  forcgg1,forcggw1,forcggs1
!~       common /ec_cvert/  omegs,omegg
!~       common /ec_cdfor/  dfor1,dfor2
!~       common /ec_ekman/  gekdis
