!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module est un portage en FORTRAN 90 du module CLIO initial
!       cree en FORTRAN 77.
!      (dans l'environnement logiciel LUDUS)
!
!      Auteur : ??, Didier M. Roche
!      Date   : ??
!      Derniere modification : 19 Aout 2014
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module dynami_mod

      use para0_mod, only: imax, jmax

      implicit none

      private :: imax, jmax

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! originally written with IMPLICIT INTEGER (I-N), REAL (A-H, O-Z)
! transferred as such in fortran90
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!               COMMONS FOR ICE DYNAMICS.
!               =========================
!
!
! nlminn   Minimum index for computation of ice drift (NH)
! nlmaxn   Maxiumum index for computation of ice drift (NH)
! nlmins   Minimum index for computation of ice drift (SH)
! nlmaxs   Maxiumum index for computation of ice drift (SH)
! zepsd1   First tolerance parameter
! zepsd2   Second tolerance parameter
! usdt     Inverse of the time step
! alpha    Coefficient for semi-implicit coriolis
! bound    Type of boundary conditions
! dm       Diffusion constant for dynamics
! om       Relaxation constant
! resl     Maximum value for the residual of relaxation
! nbitdf   Number of iterations for free drift
! nbiter   Number of sub-time steps for relaxation
! nbitdr   Max. number of iterations for relaxation
! cw       Drag coefficient for oceanic stress
! rhoco    rho0*cw
! rhoco2   rhoco*rhoco
! angvg    Turning angle for oceanic stress
! sangvg   sin(angvg)
! cangvg   cos(angvg)
! pstarh   First bulk-rheology parameter/2
! c        Second bulk-rhelogy parameter
! zetamn   Minimun value for viscosity
! creepl   Creep limit
! usecc2   1.0/(ecc*ecc)
! sber     Test if transport of ice at Bering or not
! iameth   Method for advection
! zfn      coriolis * 2
! wght     weight of the 4 neighbours to compute averages
! akappa   first group of metric coefficients
! alambd   second group of metric coefficients
! bkappa   third group of metric coefficients
! npo1i0   number of points where there is an ocenic velocity but
!          not an ice velocity
! ipo1i0   index i of these points
! jpo1i0   index j of these points
! ug       Ice velocity (x)
! vg       Ice velocity (y)
! uo       Ocean velosity used in ice dynamics (x)
! vo       Ocean velosity used in ice dynamics (y)
! uost     Fixed ocean velocity (x)
! vost     Fixed ocean velocity (y)
! hnm      Mean snow thickness
! hgm      Mean ice thickness
! uvdif    Diffusion velocity for scalars
! ren      Reynolds number for the grid
! gridsz   Grid size for diffusion constant
! dxs1     Lenght of the grid squares sides (x)
! dxs2     Lenght of the grid squares sides (y)
! dxc1     Lenght of the grid squares centres (x)
! dxc2     Lenght of the grid squares centres (y)
! zindfa   Mask for scalars
! area     Surface of a grid square
! dfhu     Modified diffusivity coefficient (x)
! dfhv     Modified diffusivity coefficient (y)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "latitd"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer      :: nlminn, nlmaxn, nlmins, nlmaxs

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "ctbqd1"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8) :: zepsd1, zepsd2, usdt, alpha, bound, dm, om, resl
      integer      :: nbitdf, nbiter, nbitdr

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "ctbqd2"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

     real(kind=8)  :: cw, rhoco, rhoco2, angvg, sangvg, cangvg, pstarh, c, zetamn, creepl, usecc2, sber, iameth


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "comgeo"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax)        :: zfn
      real(kind=8), dimension(imax,jmax,2,2)    :: wght, akappa, bkappa
      real(kind=8), dimension(imax,jmax,2,2,2,2):: alambd

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "o1i0"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer                                   :: npo1i0
      integer, dimension(5)                     :: ipo1i0, jpo1i0

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "combqd"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax)        :: ug=0.0d0, vg = 0.0d0, uo, vo, uost, vost, hnm, hgm

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "transf"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8)                              :: uvdif,ren,gridsz

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "comadv"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax)        :: dxs1, dxs2, dxc1, dxc2, zindfa, area

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ---  Block "comdff"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax)        :: dfhu, dfhv
!
      end module dynami_mod
