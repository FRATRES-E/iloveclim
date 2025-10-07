! ################    CLIMATE DEFINITIONS   ###############

!
! 0
! 1
!
#define VEGGIE 0

#define FAST_OUTPUT 0

#define CLM_INDICES 0
#define BIOM_GEN 0
#define FROG_EXP 0

!
! 0 LOVECLIM standard version
! 1 include CARAIB vegetation model, but do not run it
! 2 include CARAIB vegetation model and run it
!
#define CARAIB 0

!
! 0 not activated
! 1 if CARAIB > 0 , then write the necessary files for offline CARAIB forcing
!
#define CARAIB_FORC_W 0

!
! 0 not activated
! 1 diagnostic hourly radiation computation (e.g. for CARAIB)
!
#define HOURLY_RAD 0

!
! 0 LOVECLIM standard version
! 1 include d18Oatm of O2 (requires CARAIB > 0 and ISOATM>=2)
!
#define OXYISO 0

!
! 0 LOVECLIM standard version
! 1 include d17O of O2 (requires CARAIB > 0 and OXYISO > 0)
!
#define D17ISO 0

!
! 0 LOVECLIM standard version
! 1 include dDwax isotopes (requires CARAIB > 0, OXYISO > 0 and ISOATM>=2)
!
#define WAXISO 0

!
! 0 no effect, keep the vecode function as normal
! 1 impose the 850 A.D. CMIP LUH2 reference vegetation distribution
!
#define VEG_LUH 0

!
! 0 do not use the icemask (LOVECLIM basic version LLN)
! 1 use the provided icemask (icemask.dat)
!
#define IMSK 1

!
! 0 use the original comatm.h (LLN compatibility)
! 1 use the comatm FORTRAN90 module (GRISLI compatibility)
!
#define COMATM 1

!
! 0 use the water routing provided by "labas.dat"
! 1 use the surface water routing as computed in "routageEAU"
!
#define ROUTEAU 1

!
! 0 LLN's version of evapotranspiration (BUG included)
! 1 Do not use evapotranspiration (no BUG fix evident)
!
#define EVAPTRS 1

!
! 0 LOVECLIM standard version
! 1 Do not allow oceanic evaporation through sea-ice
!
#define EVAPSI 1

!
! 0 No additional radiative forcing from dust
! 1 Additional radiative forcing from dust (after Claquin et al., 2003)
!
#define CLAQUIN 0

!
! 0 no call to palaeo_timer (default)
! 1 call to palaeo_timer: forced orbit, GHG and ice sheet
! 2 call to palaeo_timer: forced orbit and GHG. Be sure ISM==2||3!
! 3 call to palaeo_timer: forced orbit, GHG and ice sheet at ISM resolution (ISM==3)
!
#define F_PALAEO 0

! 0 do not include iceberg code
! 1 include iceberg code but do not use iceberg fluxes from land or icesheet model  (impose Armada of icebergs)
! 2 include iceberg code and feed with excess snow (if ISM!=2) or an ice flux from the icesheet model (if ISM==1)

#define ICEBERG 0

! ################  WATER ISOTOPES DEFINITIONS  ###############

!
! 0 do not compute water isotopes in the atmosphere
! 1 water isotopes are advected in the atmosphere without fractionnation
! 2 water isotopes computed following MJ79's methodology and definitions
! 3 water isotopes computed according to the full definition of Roche (2011)
!
#define ISOATM 0
! New iso switches should be 0/1
#define WISOATM 0
#define WISOATM_RESTART 0
!
! 0 use equilibrium solid vapor fractionnation
! 1 use Jouzel & Merlivat (1984) formulation for effective kinetic SV fractionnation
!   requires ISOATM == 2
#define FRAC_KINETIK 0

!
! 0 do not link surface isotopic content of land to atmosphere (fixed values for rmoisg)
! 1 water isotopes are computed in surface water content (bmoisg prognostic)
!
#define ISOLBM 0
! New iso switches should be 0/1
#define WISOLND 0
#define WISOLND_RESTART 0

!
! 0 do not link surface oceanic isotopic content to atmosphere (fixed values for ratio_oceanatm)
! 1 water isotopes are computed in ocean but not updated in the surface ocean for atmosphere (semi-coupled mode) Requires ISOATM >= 2
! 2 water isotopes are computed in ocean ("ratio_ocean" prognostic) Requires ISOATM >= 2
!
#define ISOOCN 0
! New iso switches should be 0/1
#define WISOOCN 0

!0 do not link ice-sheet isotopic content to icebergs
!1 do link ice-sheet isotopic content to icebergs
!2 prescribe isotopic content of icebergs
!
#define ISOBERG 0
! ################    ICE-SHEET DEFINITIONS    ###############

!
! 0 do not include any code for ice-sheets components
! 1 include agISM ice-sheet source code in compilation (do not start)
! 2 include GRISLI ice-sheet code and start it (Northern Hemisphere for now)
! 3 include GRISLI, start inits but do not run it (fixed ice sheet)
!
#define ISM 0

!
! 0 use PDD computed in GRISLI for the Mass Balance (auld method)
! 1 use the ITM computed in iLCM, using net mass balance for coupling
! 2 same as 1, with an additional geographical tuning of melt param. Requires DOWNSCALING==2
!
#define SMB_TYP 0

!
! 0  if ISM is not 2 this flag should be 0 as well
! 1  couple GRISLI and ECBilt using the snow computed in ECBilt
! 2  couple GRISLI and ECBilt using the total precipitation (snow and rain)
!
#define CPLTYP 0

!
! 0 Refreezing of ice melt and/or rain is discarded for SMB computation
! 1 Refreezing is estimated from annual rain,snow,melt and temperature (Janssens and Huybrechts 2000)
!   /!\ This has to be used with the downscaling (DOWNSCALING == 2 )
!
#define REFREEZING 0

!
! 0  transmit the climate fields from ECBilt with an interpolation only
! 1  use vertical downscaling for temperatures when coupling GRISLI and ECBilt
!
#define DOWNSTS 0

#define DOWNSCALING 0
!
! 0 use GRISLI basal melting rate (depends on module_choix in GRISLI)
! 1 use CLIO basal melting rate (make sure GRISLI is compiled with bmelt_clio_coupl)
!
#define SHELFMELT 0


! 0 calved ice of GRISLI is not further regarded
! 1 calved ice of GRISLI is given to iceberg module in the form of a flux [m3/s]
! 2 calving flux is put into the ocean directly as freswhwaterflux [m3/s]
#define CALVFLUX 0

!
! !!!! only used in combination with CALVFLUX == 2 !!!!
! 0 use take-up of latent heat due to iceberg meltflux as parameterized in CLIO
!(homogenously around Greenland
! 1 latent heat is taken-up at the calving site according to the amount of freshwater
! entering the ocean
#define HEATFWF 0

!
! 0 no watercycle between ECBilt and ISM is considerd
! 1 water conservation for runoff/calving between ECBilt and ISM with all flux as liquid water
! 2 water conservation for runoff/calving between ECBilt and ISM with calving as icebergs
!
#define CONSEAU 0
!
! 0 the calving flux has no latent heat effect on the ocean
! 1 the calving flux takes up latent heat from a thickness of kheatlimit
#define HEATCALV 0

! ################    CARBON DEFINITIONS   ###############

!
! 0 no carbon code included in the compilation
! 1 add the LOCH carbon code as in the sdandard version (no start of LOCH)
! 2 add the carbon cycle of iLOVECLIM (ocean, vegetation and/or atmosphere)
!
#define CYCC 0

!
! 0 no oceanic carbon cycle with CYCC = 2 (consistency not verified yet)
! 1 oceanic carbon cycle from Six & Maier-Raimer, 1996 with CYCC = 2
!
#define OCYCC 0

!
! 0 no interactive carbon cycle (CO2 for radiative code and CO2 for carbon cycle different)
! 1 interactive carbon cycle (CO2 for radiative code = CO2 from carbon cycle)
! 2 prescribed carbon cycle (CO2 from carbon cycle = CO2 for radiative code)
!
#define INTERACT_CYCC 0

!
! 0 remove the old 14C code from LOVECLIM (LLN code) <= mandatory option when cycc=1
! 1 include the old 14C code from LOVECLIM in compilation (no start)
!
#define OLDC14 0

!
! 0 do not include the 14C code from Veronique Mariotti <= mandatory option when cycc=1
! 1 include the 14C code from Veronique Mariotti taken from CLIMBER
!
#define KC14 0

!
! 0 production for 14C is constant
! 1 production for 14C is taken from a text file (may vary with time)
!
#define KC14P 0

!
! 0 no computation of nitrous oxide in the ocean
! 1 add nitrous oxide computation in the ocean
!
#define OXNITREUX 0

!
! 0 gaz coefficient exchange rate for the carbon independent of wind speed
! 1 gaz coefficient exchange rate for the carbon is function of the wind from ECBilt
! Beware : option number 1 is not functionnal yet!
!
#define WINDINCC 0


!
! 0 Wind averaged at the ocean surface 
! 1 ERA5 wind at the ocean surface
!
#define WINDS_ERA5 0

!
! 0 O2 is a diagnostic variable in the ocean only
! 1 O2 is computed as an exchange between the ocean and the atmosphere
! Beware : option number 1 not fully functionnal yet!
!
#define O2ATM 0

!
! 0 N2O is a diagnostic variable in the ocean only
! 1 N2O is computed as an exchange between the ocean and the atmosphere
! Beware : option number 1 not fully functionnal yet!
!
#define N2OATM 0

!
! 0 MEDUSA is not activated
! 1 MEDUSA is activated
!
#define MEDUSA 0

!
! 0 sediment_loopback is not activated : the riverine input is fixed, currently read in restart file
! 1 sediment_loopback is activated : the riverine input is computed from flux to sediments
!
#define sediment_loopback 0

!
! 0 PATH is not activated
! 1 PATH is activated with fixed particles BEWARE lim reads the particles netcdf file default NEMO concentrations at modern, change file... if necessary
! 2 PATH is activated with particles calculated in ocycc BEWARE opal is not treated in iloveclim !!!
! 3 PATH is activated with particles calculated in ocycc test BEWARE opal is not treated in iloveclim !!!
!
#define PATH 0

! I am the Neodymium
#define NEOD 0

!
! O BRINES are not activated
! 1 BRINES are activated
! 2 BRINES are activated and 2D frac values read from frac.nc
! 3 BRINES are activated and 2D frac values read from frac_XXX.nc depending on the date (only with F_PALAEO)
! 4 BRINES are activated and 1D frac value read from Frac_scenario.txt
!
#define BRINES 0

!
! O O2 isotopes in ocean are not activated
! 1 O2 isotopes in ocean are activated necessitates ocean carbon cycle and O2 isotopes in atmosphere
!
#define OOISO 0

!
! O Rayleigh equation to calculate O2 isotopes in ocean are not activated
! 1 Rayleigh equation to calculate O2 isotopes in ocean 
!
#define RAYLEIGH 0 

!
! O No iron limitation in the ocean
! 1 Iron limitation : estimation according to LFe in PISCES 
!
#define IRON_LIMITATION 0

!
! O coral module is not activated
! 1 coral module is activated
!
#define CORAL 0

!
! O coastal module is not activated
! 1 coastal module is activated
!
#define COASTAL 0

!
! O without carbon emissions
! 1 wih carbon emissions
!
#define CEMIS 0

!
! O with standard remineralisation profile for POC
! 1 wih 3D remineralisation profile
!
#define REMIN 0

!
! O with standard remineralisation profile for PIC (CaCO3)
! 1 wih 3D remineralisation profile
!
#define REMIN_CACO3 0

!
! O without aragonite (only calcite for CaCO3)
! 1 wih aragonite in addition to calcite
!
#define ARAG 0

! ################    CLIO OCEAN DEFINITIONS   ###############

!
! Bathy update
! 0 do not update the bathymetry during the run
! 1 update bathymetry at the beginning of the run only, update names: fractoc.dat, mozaic.w, bath_new.txt
! 2 update the bathymetry at each restart, update names with date in them
! 3 update the bathymetry during the run, update names with date in them, NOT WORKING!!!
#define BATHY 0

! 0 no geothermal heating along the seafloor
! 1 use the geothermal heat flux mapped by Lucazeau (2019) as a bottom boundary condition in the ocean
!   /!\ This option requires file 'geothermal_heating.nc' in inputdata/clio.
!
#define GEOTHERMAL 0

! ################    UN-CLASSIFIED DEFINITIONS   ###############

!
! 0 no ABEL code
! 1 ABEL code
!
#define ABEL 0

!
! 0 no Progress bar
! 1 Progress bar
!
#define PROGRESS 0

!
! 0 no uncorrected freshwater flux in thersf to the ocean
! 1 Laurentide Ice sheet melt, 0.09 Sv (change in thersf.f)
!   WAIS definition present
!
#define UNCORFLUX 0

!
! 0 whatever uncorrected fwf is computed in switches do NOT apply it to the ocean
! 1 apply the uncorrected fwf to the ocean, from switches UNCORFLUX, UNCORRUNOFF, WATGRISCONS
!
#define APPLY_UNCORFWF 0

!
! 0 no FWF defined scenario
! 1 F_PALAEO compliant FWF (use WITH F_PALAEO == 3 and APPLY_UNCORFWF == 1)
! 2 F_PALAEO compliant FWF with berg.nc (use WITH F_PALAEO == 2 and APPLY_UNCORFWF == 1 and BC switch 2)
!
#define F_PALAEO_FWF 0

!
! 0 do not apply specific correction for LGM (cf. code in ocycc and clio)
! 1 apply LGM specific corrections. EXPERIMENTAL, use with care!
!
#define LGMSWITCH 0

!
! 0 normal code for the evolu output
! 1 add the diagnostics of Frazer Davies
!
#define FRAZER_ARCTIC 0

!
! Wrapper construction
! Do not turn to 1 if unsure !
!
#define WRAP_EVOL 0

! ################    THIS IS THE END OF GLOBAL SWITCHES !!   ###############

#include "clio_switches.h"
#include "BC_switches.h"
#include "additional_flags.h"

! ################    THIS IS THE END !!   ###############
