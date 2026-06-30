#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: veget_phys_params_mod
!
!>     @author  Didier M. Roche (dmr), refactoring by <author>
!
!>     @brief   Physical and biogeochemical parameters for the VECODE vegetation model.
!
!>     @details All scalars in this module are dimensionless or carry physical units as
!>              noted in the inline comments. They are initialised once by initcpar()
!>              (veget_sub_mod) and thereafter treated as read-only by every other
!>              routine.  No grid-level (nlat x nlon) data lives here.
!
!>     @date    Extracted from veget_mod : <date>
!>     @date    Last modification        : $LastChangedDate$
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module veget_phys_params_mod

        use global_constants_mod, only: dblp=>dp

        implicit none
        public   ! all parameters are intentionally public (read-only by convention)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Forest fraction potential — shape parameters
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: a          ! scaling factor in forshare formula            [°C^2·yr^4 mm^-4]
        real(dblp) :: bet        ! sensitivity of forest share to GDD            [°C^-1]
        real(dblp) :: gamm       ! sensitivity exponent numerator                [°C^-1]
        real(dblp) :: gamm2      ! sensitivity exponent for southern desert      [°C^-1]
        real(dblp) :: fmax       ! maximum allowed forest fraction               [–]
        real(dblp) :: gdd0_min   ! GDD0 lower threshold for forest existence     [°C·day·yr^-1]
        real(dblp) :: gdd0_max   ! GDD0 upper threshold (southern desert onset)  [°C·day·yr^-1]
        real(dblp) :: ades       ! southern desert shape parameter               [mm^-2]
        real(dblp) :: acr        ! critical precipitation scaling factor         [mm·yr^-1]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  NPP — Lieth formula parameters
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: nppmax     ! maximum NPP                                  [kgC m^-2 yr^-1]
        real(dblp) :: v1         ! Lieth precipitation sensitivity               [mm^-1]
        real(dblp) :: v2         ! Lieth temperature sensitivity                 [°C^-1]
        real(dblp) :: v3         ! Lieth temperature offset factor               [–]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  CO2 enrichment
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: betat      ! CO2 enrichment factor for trees (incl. 1/log2) [–]
        real(dblp) :: betag      ! CO2 enrichment factor for grass (incl. 1/log2) [–]
        real(dblp) :: co2ghg     ! atmospheric CO2 concentration                  [ppm]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Carbon allocation — trees
!    k1t : leaf allocation fraction     k0t : leaf turnover fraction to litter
!    k2t : stem turnover fraction       k3t : litter mineralisation fraction
!    c*t : k1t response to NPP         d*t : t1t response to NPP
!    e*t : t2t response to NPP         t1t–t4t : residence times (yr)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: c1t, c2t, c3t    ! leaf allocation k1t(NPP) coefficients  [–]
        real(dblp) :: d1t, d2t, d3t    ! leaf residence t1t(NPP) coefficients   [yr]
        real(dblp) :: e1t, e2t, e3t    ! stem residence t2t(NPP) coefficients   [yr]
        real(dblp) :: k0t               ! leaf turnover fraction to litter        [–]
        real(dblp) :: k1t               ! leaf allocation fraction                [–]
        real(dblp) :: k2t               ! stem turnover fraction                  [–]
        real(dblp) :: k3t               ! litter mineralisation fraction          [–]
        real(dblp) :: t1t               ! leaf residence time                     [yr]
        real(dblp) :: t2t               ! stem+root residence time                [yr]
        real(dblp) :: t3t               ! fast litter residence time              [yr]
        real(dblp) :: t4t               ! slow SOM residence time                 [yr]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Carbon allocation — grasses
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: c1g, c2g, c3g    ! leaf allocation k1g(NPP) coefficients  [–]
        real(dblp) :: d1g, d2g, d3g    ! leaf residence t1g(NPP) coefficients   [yr]
        real(dblp) :: e1g, e2g, e3g    ! stem residence t2g(NPP) coefficients   [yr]
        real(dblp) :: k0g               ! leaf turnover fraction to litter        [–]
        real(dblp) :: k1g               ! leaf allocation fraction                [–]
        real(dblp) :: k2g               ! stem turnover fraction                  [–]
        real(dblp) :: k3g               ! litter mineralisation fraction          [–]
        real(dblp) :: k4g               ! additional root-to-SOM transfer factor  [–]
        real(dblp) :: t1g               ! leaf residence time                     [yr]
        real(dblp) :: t2g               ! stem+root residence time                [yr]
        real(dblp) :: t3g               ! fast litter residence time              [yr]
        real(dblp) :: t4g               ! slow SOM residence time                 [yr]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Soil temperature effect on decomposition
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: ps5        ! Q10-like soil temperature sensitivity factor  [°C^-1]
        real(dblp) :: soilt      ! reference soil temperature                    [°C]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  LAI density parameters
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: deng       ! LAI per unit grass leaf biomass               [m^2 kgC^-1]
        real(dblp) :: dentd      ! LAI per unit deciduous tree leaf biomass      [m^2 kgC^-1]
        real(dblp) :: dentn      ! LAI per unit needleleaf tree leaf biomass     [m^2 kgC^-1]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Needleleaf fraction thresholds (function of leaf residence time)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: t1tn       ! t1t threshold for fully needleleaf            [yr]
        real(dblp) :: t1td       ! t1t threshold for fully broadleaf             [yr]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Desert propagation — density weighting factors
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: deng_prop  ! desert propagation weight for grass           [–]
        real(dblp) :: dentd_prop ! desert propagation weight for deciduous trees [–]
        real(dblp) :: dentn_prop ! desert propagation weight for needleleaf      [–]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Soil water / stomatal resistance parameters
!    acw* : available soil water (mm)   zr* : rooting depth factor (–)
!    rs*  : stomatal resistance (s/m)   rsd : desert stomatal resistance
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: acwd, acwt, acwg, acwn   ! available soil water           [mm]
        real(dblp) :: zrd,  zrt,  zrg,  zrn    ! rooting depth factors          [–]
        real(dblp) :: rsd,  rst,  rsg,  rsn    ! stomatal resistance             [s m^-1]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Dry-season minimum precipitation threshold
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: ps5_prc    ! minimum daily precip for vegetation (m)       [m day^-1]

#if ( CYCC == 2 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Carbon isotope fractionation factors (initialised in initcpar)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        real(dblp) :: c13frac    ! 13C fractionation factor for C3 plants        [–]
        real(dblp) :: c13frac4   ! 13C fractionation factor for C4 plants        [–]
#endif

      end module veget_phys_params_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
