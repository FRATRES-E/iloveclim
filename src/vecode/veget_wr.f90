#include "choixcomposantes.h"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      SUBROUTINE: veget_wr
!
!>     @author  dmr, refactoring by <author>
!
!>     @brief   Write vegetation diagnostics to ASCII files veget.zav and veget.outp,
!>              and trigger NetCDF output via the Vegetation_Output module.
!
!>     @details Called once per vegetation year (every nyvegt years) from veget.f90.
!>              On the first call (iyear == 0) it opens the output files and writes
!>              their headers.  On subsequent calls it computes zonal and global means,
!>              accumulates 2-D maps, and on every nwrveg-th year writes them out.
!
!>     @param[in] kveget    vegetation mode flag
!>     @param[in] nyvegt    call frequency (yr)
!>     @param[in] nwrveg    2-D map write frequency (yr)
!>     @param[in] iyr0vg    year offset for output labels
!>     @param[in] iyear     current model year
!>     @param[in] nyears    total run length (yr)
!>     @param[in] temveg    annual mean temperature accumulator    [°C]
!>     @param[in] gd0veg    annual GDD0 accumulator               [°C day yr^-1]
!>     @param[in] prcveg    precipitation accumulators            [mm yr^-1]
!>     @param[in] prcmin    minimum daily precipitation threshold  [m day^-1]
!>     @param[in] epss      land-fraction epsilon                  [–]
!>     @param[in] fracgr    land fraction per cell                 [–]
!>     @param[in] darea     latitude band area                    [m^2]
!>     @param[in] gd5veg    annual GDD5 accumulator               [°C day yr^-1]
!>     @param[in] titveg    run title string
!
!>     @date  Original : 01/10/00
!>     @date  Refactored to F90 : <date>
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine veget_wr(kveget, nyvegt, nwrveg, iyr0vg, iyear, nyears, &
                          temveg, gd0veg, prcveg,                         &
                          prcmin, epss, fracgr, darea, gd5veg,            &
                          titveg)

        use global_constants_mod, only: dblp=>dp, ip
        use comatm,               only: nlat, nlon
        use veget_mod,            only: numvegvar, st, sg, sd, snlt, blai, &
                                        anup, stock, b1, b2, b3, b4, pnpp, &
                                        newvegvar, st_moy
        use Vegetation_Output,    only: Yearly_Means, output,              &
                                        open  => open,                      &
                                        close => close,                     &
                                        write => write,                     &
                                        vegetzav_id, vegetoutp_id
        use comland_mod,          only: albland, forestfr
        use comemic_mod,          only: nwrskip, new_year_veg, current_int_veg
        use ipcc_output_mod,      only: cland
#if ( IMSK == 1 )
        use input_icemask
#endif
#if ( CLM_INDICES >= 2 )
        use global_constants_mod, only: sip
        use CLIMATE_INDICES_MOD,  only: SET_VEG_VARS
#endif

        implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Arguments
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip),       intent(in) :: kveget, nyvegt, nwrveg
        integer(ip),       intent(in) :: iyr0vg, iyear, nyears
        real(dblp),        intent(in) :: temveg(nlat,nlon)
        real(dblp),        intent(in) :: gd0veg(nlat,nlon)
        real(dblp),        intent(in) :: prcveg(nlat,nlon,2)
        real(dblp),        intent(in) :: gd5veg(nlat,nlon)
        real(dblp),        intent(in) :: prcmin, epss
        real(dblp),        intent(in) :: fracgr(nlat,nlon)
        real(dblp),        intent(in) :: darea(nlat)
        character(len=15), intent(in) :: titveg

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Saved state (replaces COMMON /cm0wvg/ and /cm1wvg/)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip), save :: nnctr
        real(dblp),  save :: vegsum  = 0.0_dblp
        real(dblp),  save :: vegsum2 = 0.0_dblp

        integer(ip), parameter :: nvgmax = 20
        real(dblp),  save      :: vegmap(nlat,nlon,nvgmax)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Local temporaries
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        integer(ip) :: i, j, k, kk, nbwr1, nbwr2
        real(dblp)  :: zavsum, totsum, vegsum2_loc
        real(dblp)  :: vegzav(0:nlat,nvgmax)
        real(dblp)  :: soiltype2(nlat,nlon)
        real(dblp)  :: ttscun, tts1

        ! output formatting — special value for non-land points
        real(dblp),        parameter :: vegspv = 999.0_dblp
        character(len=30), parameter :: fmtveg = '(65E18.5E3)                   '
        character(len=30), parameter :: fmtzav = '(65E18.5E3)                   '

        ! unit conversion factors per variable column
        ! Declared save (not parameter) so that cfmzav, titzav, titmap can be
        ! modified at runtime if needed (e.g. titmap(1) in the kveget==0 branch).
        ! Initialised via data statements — the only portable F90 way to supply
        ! compile-time repetition factors in array initialisers.
        real(dblp), save :: cfmzav(nvgmax)
        character(len=8),  save :: titzav(nvgmax)
        character(len=27), save :: titmap(nvgmax)

        !  1   2    3    4      5    6     7     8
        !  9   10   11   12–18  19   20
#if ( CYCC == 2 )
        data cfmzav / &
          100.0_dblp, 100.0_dblp, 100.0_dblp, 100.0_dblp, & ! 1-4  : fractions → %
          1.0_dblp,   1.0_dblp,                             & ! 5-6  : LAI
          25.0_dblp,                                        & ! 7    : albedo
          1.0_dblp,                                         & ! 8    : temperature
          0.001_dblp, 0.001_dblp, 0.001_dblp,              & ! 9-11 : GDD, precip
          1.0e-12_dblp, 1.0e-12_dblp, 1.0e-12_dblp,       & ! 12-14: carbon pools
          1.0e-12_dblp, 1.0e-12_dblp, 1.0e-12_dblp,       & ! 15-17: carbon pools
          1.0e-12_dblp,                                     & ! 18   : NPP
          100.0_dblp,                                       & ! 19   : ice fraction → %
          1.0_dblp                                          / ! 20   : isotope diagnostic
        data titzav / &
          'tree o/o', 'grass   ', 'desert  ', 't.needle', &
          't_lai   ', 'g_lai   ', 'albe o/o', 'temp oC ', &
          'gdd0 e-3', 'prc m/y ', 'prc_veg ', &
          'uptake G', 'stock Gt', 'b1 GtC  ', ' b2 GtC ', &
          'b3 GtC  ', ' b4 GtC ', 'NPP GtC/', &
          'icesheet', 'B1T13 ??' /
        data titmap / &
          'Fraction of Tree   (o/o) ; ', 'Fraction of Grass  (o/o) ; ', &
          'Fraction of Desert (o/o) ; ', 'Fract. of Needle L.(o/o) ; ', &
          ' L.A.I. for Tree    (1)  ; ', ' L.A.I. for Grass   (1)  ; ', &
          'An. Albedo, snow=0 (o/o) ; ',                                 &
          ' Annual Temperature (oC) ; ', ' Annual G.D.D.0  (/1000) ; ', &
          ' Annual Precipit.  (m/y) ; ', 'Reduced Prc. (Min&Dry_S) ; ', &
          ' Ann. C Uptake(kgC/y/m2) ; ',                                 &
          'Carbon Stock    (kgC/m2) ; ',                                 &
          ' B1 (kgC/m2) ;             ', ' B2 (kgC/m2) ;             ', &
          'B3 (kgC/m2) ;              ', 'B4 (kgC/m2) ;              ', &
          ' NPP (kgC/y/m2) ;          ',                                 &
          'icesheet fraction o/o ;    ',                                 &
          'B1T13 o/oo ?           ;   '                                  /
#else
        data cfmzav / &
          100.0_dblp, 100.0_dblp, 100.0_dblp, 100.0_dblp, & ! 1-4  : fractions → %
          1.0_dblp,   1.0_dblp,                             & ! 5-6  : LAI
          25.0_dblp,                                        & ! 7    : albedo
          1.0_dblp,                                         & ! 8    : temperature
          0.001_dblp, 0.001_dblp, 0.001_dblp,              & ! 9-11 : GDD, precip
          1.0e-12_dblp, 1.0e-12_dblp, 1.0e-12_dblp,       & ! 12-14: carbon pools
          1.0e-12_dblp, 1.0e-12_dblp, 1.0e-12_dblp,       & ! 15-17: carbon pools
          1.0e-12_dblp,                                     & ! 18   : NPP
          100.0_dblp,                                       & ! 19   : ice fraction → %
          0.001_dblp                                        / ! 20   : GDD5
        data titzav / &
          'tree o/o', 'grass   ', 'desert  ', 't.needle', &
          't_lai   ', 'g_lai   ', 'albe o/o', 'temp oC ', &
          'gdd0 e-3', 'prc m/y ', 'prc_veg ', &
          'uptake G', 'stock Gt', 'b1 GtC  ', ' b2 GtC ', &
          'b3 GtC  ', ' b4 GtC ', 'NPP GtC/', &
          'icesheet', 'gdd5 e-3' /
        data titmap / &
          'Fraction of Tree   (o/o) ; ', 'Fraction of Grass  (o/o) ; ', &
          'Fraction of Desert (o/o) ; ', 'Fract. of Needle L.(o/o) ; ', &
          ' L.A.I. for Tree    (1)  ; ', ' L.A.I. for Grass   (1)  ; ', &
          'An. Albedo, snow=0 (o/o) ; ',                                 &
          ' Annual Temperature (oC) ; ', ' Annual G.D.D.0  (/1000) ; ', &
          ' Annual Precipit.  (m/y) ; ', 'Reduced Prc. (Min&Dry_S) ; ', &
          ' Ann. C Uptake(kgC/y/m2) ; ',                                 &
          'Carbon Stock    (kgC/m2) ; ',                                 &
          ' B1 (kgC/m2) ;             ', ' B2 (kgC/m2) ;             ', &
          'B3 (kgC/m2) ;              ', 'B4 (kgC/m2) ;              ', &
          ' NPP (kgC/y/m2) ;          ',                                 &
          'icesheet fraction o/o ;    ',                                 &
          'Annual G.D.D.5  (/1000) ;  '                                  /
#endif

1000    format(A30,1x,E13.6,1x,I6,1x,I6,1x,I6,1x,I6)
1111    format(3(F7.1,1X,F7.3,1X),I3,A)
1112    format(3(F7.2,1X,F7.3,1X),I3,A)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Section 1 — First call: open output files and write headers
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (iyear == 0) then

          nnctr = 0
          do k = 1, len(titveg)
            if (titveg(k:k) /= ' ') nnctr = k
          end do
          vegsum  = 0.0_dblp
          vegsum2 = 0.0_dblp
          vegmap(:,:,:) = 0.0_dblp
          st_moy(:,:)   = 0.0_dblp

          ttscun = 1.0_dblp / 1000.0_dblp     ! time scale: kyr
          tts1   = real(nyvegt + iyr0vg, dblp) * ttscun
          nbwr1  = nyears / nyvegt
          nbwr2  = nyears / nwrveg

          if (kveget == 0) then
            titmap(1) = 'EcBilt Origin. Tree Frac.; '
            nbwr2 = (5 - nvgmax) * nbwr2
          else if (kveget < 0) then
            tts1  = real(iyr0vg, dblp) * ttscun
            nbwr1 = nbwr1 + 1
            nbwr2 = -nvgmax * (1 + nbwr2)
          else
            nbwr2 = -nvgmax * nbwr2
          end if

#if ( FAST_OUTPUT == 0 )
          open(newunit=vegetzav_id, &
               file='outputdata/vegetation/veget.zav', status='unknown')
          write(vegetzav_id, 1000) fmtzav, vegspv, 1+nlat, nvgmax, nbwr1, 65
          write(vegetzav_id, 1112) -87.19_dblp-5.625_dblp, 5.625_dblp, &
                                    8.0_dblp, 0.0_dblp,                  &
                                    tts1, real(nyvegt,dblp)*ttscun, -5

          open(newunit=vegetoutp_id, &
               file='outputdata/vegetation/veget.outp', status='unknown')
          write(vegetoutp_id, 1000) fmtveg, vegspv, nlon, nlat, nbwr2, 65
          write(vegetoutp_id, 1111) 0.0_dblp, 5.625_dblp, &
                                    -87.19_dblp, 5.625_dblp, &
                                    1.0_dblp, 1.0_dblp, 0
#endif
          if (kveget >= 0) return
          ! kveget < 0: fall through to write equilibrium state
        end if

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  No-vegetation mode: fill st from atmospheric tree fraction
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (kveget == 0) then
          do j = 1, nlon
            do i = 1, nlat
              if (fracgr(i,j) > epss) st(i,j) = forestfr(i,j)
            end do
          end do
        end if

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Section 2 — Compute zonal and global means
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        totsum      = 0.0_dblp
        vegzav(:,:) = 0.0_dblp

        do i = 1, nlat
          zavsum = 0.0_dblp

          do j = 1, nlon
            if (fracgr(i,j) <= epss) cycle
            zavsum = zavsum + fracgr(i,j)

            ! ice-sheet fraction (soiltype2)
            soiltype2(i,j) = 0.0_dblp
#if ( IMSK == 1 )
            if (icemask(i,j) > 0.9_dblp) soiltype2(i,j) = 1.0_dblp
#endif

            vegzav(i,1)  = vegzav(i,1)  + fracgr(i,j) * st(i,j)
            vegzav(i,2)  = vegzav(i,2)  + fracgr(i,j) * sg(i,j)
            vegzav(i,3)  = vegzav(i,3)  + fracgr(i,j) * sd(i,j)
            vegzav(i,4)  = vegzav(i,4)  + fracgr(i,j) * st(i,j) * snlt(i,j)
            vegzav(i,5)  = vegzav(i,5)  + fracgr(i,j) * st(i,j) * blai(i,j,1)
            vegzav(i,6)  = vegzav(i,6)  + fracgr(i,j) * sg(i,j) * blai(i,j,2)
            vegzav(i,7)  = vegzav(i,7)  + fracgr(i,j) * &
                             sum(albland(i,j,1:4))
            vegzav(i,8)  = vegzav(i,8)  + fracgr(i,j) * temveg(i,j)
            vegzav(i,9)  = vegzav(i,9)  + fracgr(i,j) * gd0veg(i,j)
            vegzav(i,10) = vegzav(i,10) + fracgr(i,j) * prcveg(i,j,1)
            vegzav(i,11) = vegzav(i,11) + fracgr(i,j) * prcveg(i,j,2)
            vegzav(i,12) = vegzav(i,12) + fracgr(i,j) * anup(i,j)  * darea(i)
            vegzav(i,13) = vegzav(i,13) + fracgr(i,j) * stock(i,j) * darea(i)
            vegzav(i,14) = vegzav(i,14) + fracgr(i,j) * b1(i,j)   * darea(i)
            vegzav(i,15) = vegzav(i,15) + fracgr(i,j) * b2(i,j)   * darea(i)
            vegzav(i,16) = vegzav(i,16) + fracgr(i,j) * b3(i,j)   * darea(i)
            vegzav(i,17) = vegzav(i,17) + fracgr(i,j) * b4(i,j)   * darea(i)
            vegzav(i,18) = vegzav(i,18) + fracgr(i,j) * pnpp(i,j) * darea(i)
            vegzav(i,19) = vegzav(i,19) + fracgr(i,j) * soiltype2(i,j)
            vegzav(i,20) = vegzav(i,20) + fracgr(i,j) * gd5veg(i,j)
          end do

          ! zonal mean normalisation
          totsum = totsum + zavsum * darea(i)
          do k = nvgmax, 1, -1
#if ( CYCC == 2 )
            if (((k >= 12 .and. k <= 18)) .or. k == 20) then
#else
            if (k >= 12 .and. k <= 18) then
#endif
              vegzav(i,k) = vegzav(i,k) * cfmzav(k)
              vegzav(0,k) = vegzav(0,k) + vegzav(i,k)
              cycle
            end if
            if (k >= 4 .and. k <= 6) then
              kk = (k-2) / 2
              if (vegzav(i,kk) > epss) then
                vegzav(0,k) = vegzav(0,k) + vegzav(i,k) * darea(i)
                vegzav(i,k) = cfmzav(k) * vegzav(i,k) / vegzav(i,kk)
              else
                vegzav(i,k) = vegspv
              end if
            else if (zavsum > epss) then
              vegzav(0,k) = vegzav(0,k) + vegzav(i,k) * darea(i)
              vegzav(i,k) = cfmzav(k) * vegzav(i,k) / zavsum
            else
              vegzav(i,k) = vegspv
            end if
          end do
        end do

        ! global mean
        do k = nvgmax, 1, -1
#if ( CYCC == 2 )
          if (((k >= 12 .and. k <= 18)) .or. k == 20) cycle
#else
          if (k >= 12 .and. k <= 18) cycle
#endif
          if (k >= 4 .and. k <= 6) then
            kk = (k-2) / 2
            if (vegzav(0,kk) > epss) then
              vegzav(0,k) = cfmzav(k) * vegzav(0,k) / vegzav(0,kk)
            else
              vegzav(0,k) = vegspv
            end if
          else if (totsum > epss) then
            vegzav(0,k) = cfmzav(k) * vegzav(0,k) / totsum
          else
            vegzav(0,k) = vegspv
          end if
        end do

#if ( FAST_OUTPUT == 0 )
        write(vegetzav_id,'(3A,I6)') 'Veget. Zonal_mean Output ; ', &
              titveg(:nnctr), ' ; year=', iyear + iyr0vg
        write(vegetzav_id,'(99A)') (titzav(k), k=1,nvgmax)
        do k = 1, nvgmax
          write(vegetzav_id, fmtzav) (vegzav(i,k), i=0,nlat)
        end do
        write(vegetzav_id,*)
#endif

        ! pass global carbon stock to IPCC diagnostics
        cland = vegzav(0,13)

#if ( FAST_OUTPUT == 0 )
        if (iyear + nyvegt > nyears) then
          close(vegetzav_id)
        else if (mod(iyear, nwrveg) == 0) then
          call flush(vegetzav_id)
        end if
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Section 3 — Accumulate 2-D vegetation map
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        vegsum = vegsum + 1.0_dblp
        do j = 1, nlon
          do i = 1, nlat
            if (fracgr(i,j) > epss) then
              vegmap(i,j,1)  = vegmap(i,j,1)  + st(i,j)
              vegmap(i,j,2)  = vegmap(i,j,2)  + sg(i,j)
              vegmap(i,j,3)  = vegmap(i,j,3)  + sd(i,j)
              vegmap(i,j,4)  = vegmap(i,j,4)  + snlt(i,j)
              vegmap(i,j,5)  = vegmap(i,j,5)  + blai(i,j,1)
              vegmap(i,j,6)  = vegmap(i,j,6)  + blai(i,j,2)
              vegmap(i,j,7)  = vegmap(i,j,7)  + sum(albland(i,j,1:4))
              vegmap(i,j,8)  = vegmap(i,j,8)  + temveg(i,j)
              vegmap(i,j,9)  = vegmap(i,j,9)  + gd0veg(i,j)
              vegmap(i,j,10) = vegmap(i,j,10) + prcveg(i,j,1)
              vegmap(i,j,11) = vegmap(i,j,11) + prcveg(i,j,2)
              vegmap(i,j,12) = vegmap(i,j,12) + anup(i,j)  * 1.0e12_dblp
              vegmap(i,j,13) = vegmap(i,j,13) + stock(i,j) * 1.0e12_dblp
              vegmap(i,j,14) = vegmap(i,j,14) + b1(i,j)    * 1.0e12_dblp
              vegmap(i,j,15) = vegmap(i,j,15) + b2(i,j)    * 1.0e12_dblp
              vegmap(i,j,16) = vegmap(i,j,16) + b3(i,j)    * 1.0e12_dblp
              vegmap(i,j,17) = vegmap(i,j,17) + b4(i,j)    * 1.0e12_dblp
              vegmap(i,j,18) = vegmap(i,j,18) + pnpp(i,j)  * 1.0e12_dblp
              vegmap(i,j,19) = vegmap(i,j,19) + soiltype2(i,j)
              vegmap(i,j,20) = vegmap(i,j,20) + gd5veg(i,j)
            else
              ! mark non-land cells in gd0veg for output sentinel
            end if
          end do
        end do

#if ( CLM_INDICES >= 2 )
        call SET_VEG_VARS(st, snlt, sg, sd, int(iyear, kind=sip))
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Section 4 — Write 2-D map every nwrveg years
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (mod(iyear, nwrveg) == 0) then

          vegsum2_loc = 1.0_dblp / vegsum
          st_moy(:,:) = st_moy(:,:) + vegmap(:,:,1) * vegsum2_loc
          vegsum2 = vegsum2 + 1.0_dblp

          if (iyear + nwrveg > nyears) then
            st_moy(:,:) = st_moy(:,:) / vegsum2
          end if

          do k = 1, nvgmax
            do j = 1, nlon
              do i = 1, nlat
                if (fracgr(i,j) > epss) then
                  vegmap(i,j,k) = cfmzav(k) * vegmap(i,j,k) * vegsum2_loc
                else
                  st_moy(i,j)   = 0.0_dblp
                  vegmap(i,j,k) = vegspv
                end if
              end do
            end do
          end do

          ! NetCDF output
          new_year_veg = new_year_veg + 1
          call open(Yearly_Means)
          do k = 1, numvegvar
            if (output(newvegvar(k,4))) call write(k, vegmap(1:nlat,1:nlon,k))
          end do

#if ( FAST_OUTPUT == 0 )
          do k = 1, nvgmax
            if (kveget == 0 .and. k /= 1 .and. k <= 6) cycle
            if (kveget == 0 .and. k == 19) cycle
            if (k == 11) then
              write(vegetoutp_id,'(3A,F6.3,A)') titmap(k), titveg(:nnctr), &
                    ' (>', prcmin*1000.0_dblp, ' mm/day)'
            else
              write(vegetoutp_id,'(2A)') titmap(k), titveg(:nnctr)
            end if
            write(vegetoutp_id,*) 'year=', iyear + iyr0vg
            do i = 1, nlat
              write(vegetoutp_id, fmtveg) (vegmap(i,j,k), j=1,nlon)
            end do
            write(vegetoutp_id,*)
          end do

          if (iyear + nwrveg > nyears) then
            close(vegetoutp_id)
          else
            call flush(vegetoutp_id)
          end if
#endif

          call close()

          ! reset NetCDF year counter on restart boundaries
          if (mod(iyear, nwrskip) == 0 .or. iyear == nyears) then
            new_year_veg    = 0
            current_int_veg = current_int_veg + nwrskip
          end if

          ! reset 2-D map accumulators
          vegsum        = 0.0_dblp
          vegmap(:,:,:) = 0.0_dblp

        end if  ! mod(iyear, nwrveg)

      end subroutine veget_wr

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
