! Module of carbonate production by coral reef
! Written by N. Bouttes, G. Munhoven, V. Brovkin, M. Berger, L. Kwiatkowski
! 16/05/2023
!----------------------------------------------------------------------
#include "choixcomposantes.h"

#if ( CORAL == 1 )
      module coral_mod

! use section
      use declars_mod
      use loveclim_transfer_mod, ONLY: ZX, TYER, TDAY,                  &
            SQRO2, MGT, NYR, mid_level, joursemaine, KLSR
      use C_res_mod, ONLY: coral_res_fich, c_nino_fich
      use global_constants_mod, only: dblp=>dp, ip
      use omega_mod, only : calc_omega_ar, omega_arag3D

! declaration section
       implicit none

! input output for corals
      REAL area
      REAL temp
      REAL sal
      REAL light_surf
      REAL arag_sat
      REAL depth 
      REAL :: mass_carb_new
      REAL :: net_carb
      !REAL :: mass_carb
      REAL d_sl

! local variables
      INTEGER, PARAMETER :: kmax_hypso = 300 !128
      REAL bati
      INTEGER mask !0=ocean, 1=land
      INTEGER t,i,j,k
!      REAL P_carb_an
!      dimension P_carb_an(LT, JT, NOC_CBR)
      REAL, dimension(LT, JT, NOC_CBR) :: P_carb_an
      !REAL mass_carb_old
      INTEGER compteur              !day of year (1 to 360)

      !REAL, dimension(LT, JT, NOC_CBR) :: hypso
      !REAL, dimension(LT, JT, NOC_CBR-2) :: hypso_read
      ! on subgrid depth
      REAL(kind=dblp), dimension(LT, kmax_hypso, NOC_CBR) :: hypso
!tbd      REAL(kind=dblp), dimension(LT, kmax_hypso, NOC_CBR-2) :: hypso_read
      REAL, dimension(LT, kmax_hypso, NOC_CBR) :: topoff
!tbd      REAL, dimension(LT, kmax_hypso, NOC_CBR-2) :: topoff_read
      REAL, dimension(LT, NOC_CBR) :: kd_490
!tbd      REAL, dimension(LT, NOC_CBR-2) :: kd_490_read
      REAL, dimension(kmax_hypso) :: level_hypso_read
      REAL, dimension(kmax_hypso) :: level_hypso
      REAL, dimension(kmax_hypso+1) :: level_bounds_hypso_read
      REAL, dimension(kmax_hypso+1) :: level_bounds_hypso
!      REAL level, latitude, longitude
      REAL, dimension(JT) ::  level
      REAL, dimension(NOC_CBR,LT) :: latitude
      REAL, dimension(NOC_CBR, LT) :: longitude
      character*256 testchar
      REAL, dimension(LT, JT, NOC_CBR) ::  mass_carb_prev
      REAL aragonite
      INTEGER d_hypso
      REAL level_hypso_inf, level_hypso_sup

!a supprimer      INTEGER year
      
      ! character*256 FILE_NAME 
      ! parameter(FILE_NAME = "toto.nc")
      ! integer status_nc, ncid
 

! constants
      REAL tau !constante de temps de wash out des coraux
      REAL dt                       ! time step
      REAL caco3_molar_mass

      INTEGER i_sl
      REAL dum
      INTEGER, PARAMETER :: t_sl=15
      REAL, dimension(t_sl) :: sea_level

      !REAL, dimension(JT) :: level_inf
      !REAL, dimension(JT) ::  level_sup
      !REAL depth, depth_inf, depth_sup
      REAL, dimension(LT, JT, NOC_CBR) ::  area_3D
      INTEGER kk

      ! a supprimer REAl par
      REAL PAR_m
      REAL dens_ocn
      REAL water_z_pp

      !REAL weathering_odic
      !REAL weathering_oalk
      ! weathering
      REAL C_sed
      REAL C_riv
      REAL C_car_a
      REAL A_riv

      REAL surface_ocean

      ! for bleaching
      REAL tau_bleach_moderat
      REAL tau_bleach_strong
      REAL dhw_thresh_moderat
      REAL dhw_thresh_strong

      REAL, dimension(LT,NOC_CBR) :: temp_too_low_all

      ! dissolution
       REAL :: lambda_diss

! output
      REAL total_mass_carb          ! in Pg CaCO3/yr
      REAL total_prod_carb          ! in Pg CaCO3/yr, production Gglobal
      REAL total_area_coral_an
      REAL total_area_coral_an_40m  ! same but limited to depth <= 40 m
      REAL total_prod_coral_an
      REAL total_prod_coral_an_40m  ! same but limited to depth <= 40 m
      REAL total_mass_coral_an
      REAL, dimension(LT,JT,NOC_CBR) :: coral_prod
      REAL, dimension(LT,JT,NOC_CBR) :: coral_area
      REAL, dimension(LT,JT,NOC_CBR) :: coral_cum_mass
      REAL, dimension(LT,kmax_hypso,NOC_CBR) :: coral_mass_subgrid
      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: coral_prod_out
      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: coral_mass_out
      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: coral_area_out
!tbd      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: omega_arag3D
      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: tau_bleach_out
      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: DHW_out
!tbd      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: PH
      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: TM_fix
!nb comprend pas pourquoi ca n est pas de la dimension de scal...
      REAL, dimension(LT,JT,NOC_CBR) :: corproda, corareaa
      REAL, dimension(LT,JT,NOC_CBR) :: omegaa, oco3a
      REAL, dimension(LT,JT,NOC_CBR) :: DHWa

      REAL coral_CO2
      character*256 outfilename
      INTEGER js
      REAL test_sum

! for bleaching
      REAL, dimension(LT,JT,NOC_CBR) :: tau_bleach
      REAL, dimension(LT,JT,NOC_CBR) :: timebleach
      REAL, dimension(LT,JT,NOC_CBR) :: DHW         ! degree heating week
      REAL, dimension(LT,JT,NOC_CBR) :: DHW_nb      ! nb of weeks when bleaching is activated (moderate)
      REAL, dimension(LT,JT,NOC_CBR) :: MMMclim     ! monthly mean max climatological temperature
!tbd      REAL, dimension(LT,JT,NOC_CBR) :: xsHS 
      REAL, dimension(LT,JT,NOC_CBR) :: temp_mois
      INTEGER, PARAMETER :: window_MMM=30           !30 years: window length for running mean of monthly temperature
      INTEGER, PARAMETER :: nb_mois=window_MMM*12
      REAL, dimension(LT,JT,NOC_CBR,nb_mois) :: temp_mois_all
      INTEGER indice_mois
!tbd      INTEGER, PARAMETER :: window_week=17 !number of weeks for DHW
      INTEGER, PARAMETER :: nb_hs=84 !window_week*joursemaine !17 weeks of 5 days = 85 days
      REAL, dimension(LT,JT,NOC_CBR,nb_hs) :: xsHS_all
      INTEGER indice_hs
      ! Temperature with variability for bleaching
      REAL, dimension(LT, JT, NOC_CBR) :: TM_pluswkvar
      INTEGER, SAVE :: k_mbiota_rand
      REAL :: nino3, nino3_var

      ! t0 for dissolution when no production
      REAL, dimension (LT, kmax_hypso, NOC_CBR) :: t0_diss_all
      REAL, dimension (LT, kmax_hypso, NOC_CBR) :: prod_before_all

      contains

! functions and subroutines
!----------------------------------------------------------------------
       subroutine mbiota_corals(i,j,n)
       use const_mod, only: gpes, rho0
       use carbonate_speciation_mod, only: calc_co3sat_arag, incche

       use marine_bio_mod, only: OCO3, OALK, ODIC, OC13, jprod, OPO4
       use loveclim_transfer_mod, only: TM, SM, mid_level, sabst_o, zx, &
            kendy, dvol, tday, tyer, tstoc, sqro2
       use mbiota_mod, only: scale_m, scanu


       integer(kind=ip), intent(in) :: i,j,n

! nb corals
       real(kind=dblp) :: area_coral, temp, sal, phos, light_surf,      &
       omega_arag, depth, kd, topof, i_pp, k1p, k2p, kbp, csatp, co2_pp &
       , hco3_pp, co3_pp, temp_coral_CO2, coral_13CO2, modif, p_bar,    &
       CO3sat_ar, sCO2, xpCO2, xCO2, xHCO3, xCO3                        &
       ,temp_too_low 

       integer(kind=ip):: js

!      REAL alk_pp
!      REAL tco2_pp

!nb corals
!nb called at the end of the marine bio day, after maphot

      !write(*,*) 'j', j

      ! depth
      !depth=mid_level(j+1) !-20 !zw is negative from 0 to -5500

      !temperature
      temp=TM(i,j,n) !25
      !temp=TM_fix(i,j,n) ! test with temperature fixed at end of first
                          !year
      temp_too_low=temp_too_low_all(i,n)
      !temp=25 !test

      ! salinity
      sal=SM(i,j,n) !34.7
      !sal=34.7 !test

      ! phosphate
      !print*, 'phosphate', OPO4(i,j,n)
      phos=OPO4(i,j,n) !mumol/L

      ! light at surface
      light_surf= SABST_O(i,n)


      ! light attenuation coefficient
      kd=kd_490(i,n)
      !kd=0.1 !test
      if (kd .le.0) then
          !write(*,*) 'kd', i,n, kd
          kd=0.1
          kd_490(i,n)=0.1
      endif


      !loop on subgrid depth
      !write(*,*) ' '
      !write(*,*) 'ZX ', (-1)*ZX(j), (-1)*ZX(j+1)
      do js=1,kmax_hypso

        !depth
        !depth=mid_level(j+1) !-20 !zw is negative from 0 to -5500

        !all depth in negative except ZX hence the (-1)*ZX
!        write(*,*) 'test depth ', level_bounds_hypso(js),
!     >             (-1)*ZX(j), (-1)*ZX(j+1)
        if ((level_bounds_hypso(js+1).le.(-1)*ZX(j)) .and.              &
           (level_bounds_hypso(js+1).gt.(-1)*ZX(j+1))) then
!        write(*,*) 'test depth ', level_bounds_hypso(js+1),
!     >             (-1)*ZX(j), (-1)*ZX(j+1)
          depth=level_hypso(js)
         ! write(*,*) 'depth', depth ! depth is negative

         !area: surface of continental bottom in the grid cell in m2
         !area_coral=SQRO2(i,n)
         !print*, 'area coral (1e3 km2)', i,n, area_coral*1e-9
         area_coral=hypso(i,js,n)
         !area_coral=hypso(i,js,n)*0.1 !10percent of the area for corals accounting for inter-reefal area
         !print*, 'area_coral in 1e3 km2',i,js,n, area_coral*1e-9

         ! corals can exist only if some surface area
         if (area_coral .gt. 1e-12) then
          !print*, 'area_coral in 1e3 km2',i,js,n, area_coral*1e-9

          !topographic factor
          topof=topoff(i,js,n)
          !topof=1 !test
          !if ((topof .le.1) .or. (topof .ge.0)) then
          !!  !print*, 'ok', i,js,n
          !else
          !  print*, 'topof',i,js,n, topof
          !endif

          !aragonite saturation
!         alk_pp     = OALK(i,j,n)*dens_ocn ! convert from eq/kg to eq/l
!         tco2_pp    = ODIC(i,j,n)*dens_ocn ! convert from mol/kg to
!         mol/l

          !i_pp=1


!         CALL calc_k(temp,sal,water_z_pp,k1p,k2p,kbp,csatp,
!     &                i_pp)


!         CALL calc_buff(alk_pp,tco2_pp,sal,k1p,k2p,kbp,
!     &                co2_pp,hco3_pp,co3_pp,
!     &                i_pp)


!nb just below to be removed: done in calc_omega
!nb       computes CO32- concentration at saturation
          ! needs temperature in Kelvin
          ! pressure in bar=rho * g * z = Kg/m3 * N/kg * m = N/m2 = Pa
          ! 1 bar = 1e5 Pa
          ! rho = b * rho0 /g ?
          ! en attendant use rho0
          ! beware mid_level is negative...
          !write(*,*) 'computation of p_bar', rho0, gpes,mid_level(j+1)
          !p_bar= rho0 * gpes * (-1) * mid_level(j+1)*1e-5
!          p_bar= rho0 * gpes * (-1) * depth*1e-5
!          CO3sat_ar=calc_co3sat_arag(temp+273.15, sal, p_bar)
          !CO3sat_ar=CO3sat_ar ! in mol/kg
!            omega_arag=OCO3(i,j,n)/CO3sat_ar
!nb end of to be removed


          omega_arag=omega_arag3D(i,j,n)

          !omega_arag=3 !test
          !omega_arag3D(i,j,n)=omega_arag
          !if (OCO3(i,j,n).gt.0) then
          ! print*, 'omega_arag ', OCO3(i,j,n)*1e6, CO3sat_ar*1e6,
          ! omega_arag
          !endif

!          print*,'coral ', 90-i*2.5, j, n, co3_pp/csatp/1.5
!     &                , temp_p, sal_p, water_z_p, k1, k2, kb
!     &                , co3_p*1.E6, csat*1.E6, co3_p/csat/1.5

!          OCO3(i,j,n)    = co3_pp*1.E6  ! convert to µmol/l

!          OCO3sat(i,j,n) = csatp*1.E6   ! convert to µmol/l

!          omega_arag     = co3_pp/csatp/1.5
           !omega_arag=2

      !for bleaching: computes tau_bleach and timebleach later used in coral to limit production
      !only after the first 30 years when the reference is computed for
      !MMMclim
      !now from the beginning using value from restart for MMMclim-> or
      !not
      ! if (NYR.gt.30) then
       if (NYR.gt.window_MMM) then
         call calc_bleach(i,j,n)
       endif


!      if (temp.gt.20) then
!      print*, 'coral input ', area_coral, temp, sal, light_surf,
!     >          omega_arag, depth, kd, topof
!      endif

      !call subroutine computing CaCO3 production by corals
      !call  corals(area,temp,sal,light_surf,omega_arag,depth
      !&
      !      ,mass_carb_new,kd,topof)
      call corals(area_coral,temp,sal,phos,light_surf,omega_arag,depth, &
         kd,topof,tau_bleach(i,j,n),timebleach(i,j,n),temp_too_low,     &
         coral_mass_subgrid(i,js,n), t0_diss_all(i,js,n),               &
         prod_before_all(i,js,n))

      !output of corals is net_carb in Pmol/day
      !if (net_carb.ne.0) then
      !   write(*,*) 'net_carb mbiota ', i,j,n, net_carb
      !endif
      coral_prod(i,j,n)=coral_prod(i,j,n)+net_carb !annual net production in Pmol/year
      !g_prod(i,j,n)=g_prod(i,j,n)+mass_carb_new !mass_carb_new
      !equivalent de prom
!      if (coral_prod(i,j,n).ne.0) then
!         write(*,*) 'coral_prod(i,j,n) mbiota ',i,j,n,coral_prod(i,j,n)
!      endif


      !New mass in g/day (net_carb in Pmol/day, MCaCO3 in g/mol))
      mass_carb_new=net_carb*caco3_molar_mass*1.0e15
      !if (mass_carb_new.ne.0) then
      !   write(*,*) 'mass_carb_new', i,j,n,mass_carb_new
      !endif


      !net mass
      !mass_carb_old=mass_carb
      !mass_carb=mass_carb_old+mass_carb_new


      !print*, 'corals ', mass_carb
      !if (net_carb.gt.0) then
      total_prod_coral_an=total_prod_coral_an+net_carb !Production in Pmol/day -> Pmol/year
      if ( depth .ge. -40) then
        total_prod_coral_an_40m=total_prod_coral_an_40m+net_carb !Production in Pmol/day -> Pmol/year
      endif
      total_mass_coral_an=total_mass_coral_an+mass_carb_new*1e-15!mass_carb en g *1e-15 g->Pg -> Pg/year
      !total_area_coral_an=total_area_coral_an+area_coral
      !!(area_coral*topof*1e-6) !in m2->km2 *1e-6
      coral_cum_mass(i,j,n)=coral_cum_mass(i,j,n)+mass_carb_new !g
      coral_mass_subgrid(i,js,n)=coral_mass_subgrid(i,js,n)             &
                                 +mass_carb_new
      !endif

      !if ((coral_cum_mass(i,j,n).gt.0).and.(KENDY.eq.1)) then ! if last
      !day of year and area with corals
      if ((coral_mass_subgrid(i,js,n).gt.0).and.(KENDY.eq.1)) then ! if last day of year and area with corals
      !if ((coral_mass_subgrid(i,js,n).gt.0)) then ! if last day of year
      !and area with corals
         total_area_coral_an=total_area_coral_an+(area_coral*topof*1e-6) !in m2 *1e-6 -> in km2
         if (depth .ge. -40) then
           total_area_coral_an_40m=total_area_coral_an_40m              &
                        +(area_coral*topof*1e-6) !in m2 *1e-6 -> in km2
         endif
         coral_area(i,j,n)=coral_area(i,j,n)+area_coral*topof !in m2
        ! write(*,*) 'coral area mbiota in 1e3 km2 ',
        ! coral_area(i,j,n)*1e-9
      endif

         endif ! area_coral > 0

        endif !within clio levels
      enddo !js


!      if (KENDY.eq.1) then
!      write(*,*) 'total_area_coral_an (km2) ', total_area_coral_an
!      write(*,*) 'total_area_coral_an (km2)  <= 40m', total_area_coral_an_40m
!      endif


!!!!! A supprimer ensuite, pour test
!      coral_prod_out(i,j,n)=coral_prod(i,j,n)!
!      if (coral_prod_out(i,j,n).ne.0) then
!         write(*,*) 'coral_prod_out(i,j,n) mbiota '
!     >      ,i,j,n,coral_prod_out(i,j,n)
!      endif

      !nb impact on ocean variables
      ! coral_CO2 is a daily carbonate production in Pmol/day in grid box
                                    ! g_prod has units of Pmol/day
                                    ! converted her to umol/kg/timestep
                                    ! SCALE_M = 10^-18
                                    ! SCANU = 10^-6
      !note TDAY=86400s et TSTOC=86400s -> TDAY/TSTOC=1
      !temp_coral_CO2 = g_prod(i,j,n)
      temp_coral_CO2 = net_carb !in Pmol/day

     ! SCALE_M=1e-18 to convert from 1e3 mumol/m3 to GtC
     ! SCANU =1e-6 to convert from mumpl/kg to mol/kg

     ! OALK in (eq/kg)
         modif=temp_coral_CO2*1./SCALE_M*SCANU/DVOL(i,j,n)/(TDAY/TSTOC)
         !modif=0
         !if (modif.ne.0) then
         ! print*, 'OALK avant ',i,j,n
         ! print*, OALK(i,j,n)
         ! print*, 'added modif', modif
          OALK(i,j,n)    = OALK(i,j,n) - 2.0 * modif
         ! print*, 'OALK apres', OALK(i,j,n)
         !endif

          ODIC(i,j,n)    = ODIC(i,j,n) - modif

          OC13(i,j,n)    = OC13(i,j,n) - 1.5 * modif

          coral_CO2      = coral_CO2 + temp_coral_CO2/(TDAY/TSTOC)*12.

          coral_13CO2    = coral_CO2*1.5


! dissolution of weathered bicarbonate in the surface water
! nb for now: homogenous value everywhere at ocean surface
! to be changed to river input
!C_riv and A_riv carbon and alkalinity brought by rivers to ocean

      if (j.eq.1) then !at the surface

          OALK(i,j,n) = OALK(i,j,n) + A_riv*SQRO2(i,n)/surface_ocean    &
            /SCALE_M*SCANU/DVOL(i,j,n)/(TYER/TSTOC)

!     >      weathering_oalk*SQRO2(i,n)/surface_ocean

          ODIC(i,j,n) = ODIC(i,j,n) + C_riv*SQRO2(i,n)/surface_ocean    &
            /SCALE_M*SCANU/DVOL(i,j,n)/(TYER/TSTOC)

!     >      weathering_odic*SQRO2(i,n)/surface_ocean
!          OC13(i,j,n) = OC13(i,j,n) + weathering_oc13(i,n)/1000.
!     >      /SCALE_M*SCANU/DVOL(i,j,n)/(TYER/TSTOC)
      endif !j

       end subroutine mbiota_corals
!-------------------------------------------------------------------------------------------------------------------------------------



!----------------------------------------------------------------------
! initialisation
      subroutine ini_coral

      use ncio,      only: nc_read
      use loveclim_transfer_mod, ONLY: joursemaine
      use loveclim_transfer_mod, only: TM

! local
      integer :: n

!constantes
      tau=4000 !in years
      dt=1 !pas de temps en jours
      PAR_m=0.4
      dens_ocn   = 1.03         ! kg/l
      water_z_pp = 25
      ! Molar mass of CaCO3 = 100.0869 g/mol
      caco3_molar_mass=100.0869
      !weathering constant in Pmol/day (1e15) (on veut environ 7.5-9 Tmol/year (1e12)).
      !TYER/TDAY=360 day/year
      ! all in form of hco3-
      !C_sed=enfouissement de CaCO3 dans les sediments, inluant les coraux
!      C_sed=7.5*1e-3/(TYER/TDAY) ! C_sed=G_coral a l equilibre
      !C_riv= carbon from carbonate dissolution brought by rivers to the
      !ocean all in HCO3- form
      C_riv=3.9*1e-5  !2*C_sed !Pmol/day
      ! C_sil_a= 0 ! consommation CO2 par alteration silicates
      ! C_vol+C_hyd =0 ! Carbon from volcanism and hydrothermals going
      ! to atmsophere
      ! C_car_a= consommation of CO2 by carbonate alteration
      C_car_a=C_riv /2. 
      !Alkalinity
      A_riv=C_riv
      !for bleaching
      !bleaching effect, time scales (yrs)
      tau_bleach_moderat =  5. !20. 
      tau_bleach_strong  = 20. !100.
      !DHW index bleaching thresholds, using weeks of 5 days
      !(joursemaine=5)
      dhw_thresh_moderat = 4.!  *7./joursemaine
      dhw_thresh_strong  = 8.!  *7./joursemaine


      !constante for dissolution
      lambda_diss=1/10.

! initialisation
      total_area_coral_an=0
      total_area_coral_an_40m=0
      total_prod_coral_an=0
      total_prod_coral_an_40m=0
      total_mass_coral_an=0

      coral_CO2=0

      coral_prod(:,:,:)=0.
      coral_area(:,:,:)=0.
      coral_cum_mass(:,:,:)=0.
      coral_mass_subgrid(:,:,:)=0.
      coral_prod_out(:,:,:)=0.
      coral_area_out(:,:,:)=0.
      coral_mass_out(:,:,:)=0.
      omega_arag3D(:,:,:)=0.
      tau_bleach_out(:,:,:)=0.
      DHW_nb(:,:,:)=0. ! nb of weeks with bleaching
      DHW_out(:,:,:)=0.
!tbd      PH(:,:,:)=0.
      TM_fix(:,:,:)=0.
      !dissolution
      t0_diss_all(:,:,:)=0.0
      prod_before_all(:,:,:)=0.0

! for bleaching
!      if (KLSR.eq.0) then ! only if fresh start, otherwise read in restart
        tau_bleach(:,:,:)=0.0
        timebleach(:,:,:)=0.0
        MMMclim(:,:,:)=0.0 !max of monthly temperature
!      endif
      !indice_mois=0
      indice_mois=nb_mois
      temp_mois(:,:,:)=0.
      DHW(:,:,:)=0.0
      MMMclim(:,:,:)=0.0 !max of monthly temperature
!tbd      xsHS(:,:,:)=0.
      indice_hs=1
!tbd  indice_hs=nb_hs
      xsHS_all(:,:,:,:)=0.
      nino3=0.0
      nino3_var=0.0

! limite temperature low
      temp_too_low_all(:,:)=0


! Read files

!area (surface and hypsometry) from gebco
      print*, 'read hypsometry'
      !!call nc_read_attr("GEBCO/hypsometry.nc", "title", testchar)
      !!write(*,*) "Title: ", trim(testchar)
!tbd      call nc_read("inputdata/hypsometry_clio.nc","hypso",hypso_read) !(:,:,:)
      call nc_read("inputdata/hypsometry_clio.nc","hypso",hypso) !(:,:,:)
      call nc_read("inputdata/hypsometry_clio.nc","level",              &
                   level_hypso_read) !(:,:,:)
      call nc_read("inputdata/hypsometry_clio.nc","level_bounds",       &
                   level_bounds_hypso_read) !(:,:,:)
!t      call nc_read("../../input/GEBCO/hypsometry_60_120.nc","level",level_hypso)
!t      call nc_read("../../input/GEBCO/hypsometry_60_120.nc","latitude",latitude)
!t      call nc_read("../../input/GEBCO/hypsometry_60_120.nc","longitude",longitude)
!t      d_hypso=level_hypso(2)-level_hypso(1)
!t      print*,'resolution vertical hypsometry ', d_hypso

      do js=1,kmax_hypso
        level_hypso(js)=level_hypso_read(kmax_hypso-js+1)
        level_bounds_hypso(js)=level_bounds_hypso_read(kmax_hypso+1-js+1)
      enddo
      level_bounds_hypso(kmax_hypso+1)=level_bounds_hypso_read(1)

!      print*, 'hypso size', ubound(hypso_read)
!      print*, 'level_hypso ', level_hypso
!      print*, 'level_bounds_hypso', level_bounds_hypso
!      print*, 'mid_level', mid_level
!      print*, 'ZX', ZX

     ! print*, 

      !test global value
      j=1
      test_sum=0
      do n=1,NOC_CBR
        do i=1,LT
          !if (MGT(i,j,n).eq.1) then
            !do js=51,153 !103m depth
            do js=148,250 !103m depth
            !do js=1,kmax_hypso !103m depth
            !do js=kmax_hypso-153,kmax_hypso-51 !103m depth
               !write(*,*) 'js', js, level_hypso(js)
!tbd               test_sum=test_sum+hypso_read(i,js,n)
               test_sum=test_sum+hypso(i,js,n)
            enddo
          !endif
        enddo
      enddo
      write(*,*) 'surface for levels between '                           &
                 , level_hypso(148), level_hypso(250)
      write(*,*) 'surface totale (1e3 km2)=',test_sum*1e-3*1e-6

!below tbd
      !nb tableau a bord replie :on recopie les colonnes croisees 
      ! valeur(1)=valeur(nmax-1) et
      ! valeur(nmax)=valeur(2)
!      do n=2,NOC_CBR-1
!        hypso(:,:,n)=hypso_read(:,:,n-1)
!      enddo
!      hypso(:,:,1)=hypso(:,:,NOC_CBR-1)
!      hypso(:,:,NOC_CBR)=hypso(:,:,2)

! topographic factor topof
      print*, 'read topof'
!tbd      call nc_read("inputdata/topof_clio.nc","topof",topoff_read) !(:,:,:)
      call nc_read("inputdata/topof_clio.nc","topof",topoff) !(:,:,:)
      !nb tableau a bord replie :on recopie les colonnes croisees 
      ! valeur(1)=valeur(nmax-1) et
      ! valeur(nmax)=valeur(2)
!tbd      do n=2,NOC_CBR-1
!tbd        topoff(:,:,n)=topoff_read(:,:,n-1)
!tbd      enddo
!tbd      topoff(:,:,1)=topoff(:,:,NOC_CBR-1)
!tbd      topoff(:,:,NOC_CBR)=topoff(:,:,2)



!Kd_490 (attenuation of insolation)
      print*, 'read Kd_490'
!tbd      call nc_read("inputdata/Kd_490_clio.nc","kd_490",kd_490_read)
      call nc_read("inputdata/Kd_490_clio.nc","kd_490",kd_490)
      !nb tableau a bord replie :on recopie les colonnes croisees
      ! valeur(1)=valeur(nmax-1) et
      ! valeur(nmax)=valeur(2)
!tbd      do n=2,NOC_CBR-1
!tbd        kd_490(:,n)=kd_490_read(:,n-1)
!tbd      enddo
!tbd      kd_490(:,1)=kd_490(:,NOC_CBR-1)
!tbd      kd_490(:,NOC_CBR)=kd_490(:,2)



      !For weathering flux in each ocean surface grid cell
      !total surface of ocean
      surface_ocean=0.
      j=1 !at the surface of ocean
      do n=1,NOC_CBR
        do i=1,LT
          if (MGT(i,j,n).eq.1) then
            surface_ocean=surface_ocean+SQRO2(i,n)
          endif
        enddo
      enddo 


!initialise TM_fix
      do n=1,NOC_CBR
       do j=1,JT
        do i=1,LT
          if (MGT(i,j,n).eq.1) then
            TM_fix(i,j,n)=TM(i,j,n)
          endif
        enddo
       enddo
      enddo

! call restart in read mode
!      call restart_coral(1) 
 

      end subroutine ini_coral

!----------------------------------------------------------------------

!----------------------------------------------------------------------
! coral production
      subroutine corals(area,temp,sal,phos,light_surf,omega_arag,depth, &
            kd,topof,tau_bleach_l, timebleach_l, temp_too_low, mass_carb&
            , t0_diss, prod_before)

!input output

      REAL area
      REAL temp
!      REAL tomin
!      REAL tomax
      REAL sal
      REAL light_surf
      REAL omega_arag !arag_sat
      REAL depth
      REAL mass_carb !, mass_carb_new
      !REAL d_sl
      REAL kd
      REAL topof
      REAL phos
!      REAL nitr
      REAL tau_bleach_l
      REAL timebleach_l
      REAL temp_too_low
      REAL t0_diss
      REAL prod_before


!local
!     to be moved to init and module
      REAL gmax_coral
      REAL zmax
      REAL Ik !RAD_IK
      REAL Imin
      REAL coral_dens
      REAL pk490
      REAL tmin, tmax
      REAL smin, smax
      REAL nmax, pmax
      REAL coef_ftemp
      REAL arag_k
      REAL dens_ocn
      REAL water_z_pp
      integer i_pp
      REAL tfmin, tfmax
      REAL a_temp, b_temp


!to be deleted      REAL area_prod
      REAL Iz
      !REAL mass_carb_old
      REAL P_carb, D_carb
      REAL g_coral
      REAL RAD_m
      REAL temp_factor


! maximum vertical accumulation rate for corals in m/yr -> per day for us (-> /360) -> m/day
!      gmax_coral=1.04/360 !mm/day
!      gmax_coral=1.0/1000./360 !m/day
!      gmax_coral=(4.0*1E-3)/360 !in m/day (from mm/year)
      gmax_coral=(3.0*1E-3)/360. !in m/day (from mm/year)

! Coefficient for Temperature dependency
      coef_ftemp = 0.2      
! target - caco3 flux of 0.105 Gtc/yr
! saturation light intensity Ik, in W/m2  (conversion between W/m2 and muE m-2 s-1 using 4.6 from Kirk, 1994)
      Ik=350./4.6 !RAD_IK ! from Guypour CLIMBER
      !Ik=50 !200 !50 !?? A voir ?? !pour moi : mE/m2/s
      ! Ik devrait etre en W/m2 ?

! minimum light intensity necessary for reef growth
      Imin=300./4.6 !in W/m2  (conversion between W/m2 and muE m-2 s-1 using 4.6 from Kirk, 1994)

! CaCO3 density, kg/m3 (density of 2.89g/cm3 and porosity 50%)
      coral_dens=1.445*1.0E3 ! in kg/m3
! light extinction coefficient  !!!
! to be deleted      pk490=0.1       
! minimum temperature for coral growth (grad. C)
      tmin= 18.1   
! maximum temp. for coral growth     
      tmax= 31.5
! minimum salinity for coral growth (--)
      smin= 30.0
! maximum salinity for coral growth     
      smax= 39.0
!maximum phosphate value for coral growth
      pmax=0.2 !micromol/L
!cc sea level (m)
!cc      isea_lev=250
!c supersaturation parameter
      arag_k=2.86

      dens_ocn   = 1.03         ! kg/l
!temperature min and max for linear production function (as a function
!of temperature)
      tfmin=18
      tfmax=31
      b_temp=1./(tfmax-tfmin)
      a_temp=1-(tfmax*b_temp)

!to be deleted      water_z_pp = 25

!to be deleted      area_prod=0

!!!      sum_prod_coral=0 


! temperature limitation
      if ((temp.lt.tmin).or.(temp.gt.tmax)                              &
          .or.(sal.lt.smin).or.(sal.gt.smax)                            &
          !.or.(phos.gt.pmax).or.(nitr.gt.nmax)) then
          !.or.(phos.gt.pmax)) then
          .or.(phos.gt.pmax).or.(temp_too_low.gt.0)) then

            P_carb = 0

            if (temp.gt.tmax) THEN
                                    ! Temperature excess always leads
                                    ! to strong bleaching
              tau_bleach_l = tau_bleach_strong
             timebleach_l = NYR

            endif

      else

! surace incoming shortwave radiation light_surf in W/m2, multiplied by PAR_m=0.4
! to have PAR stil in W/m2
              RAD_m=light_surf*PAR_m ! coupled version: RAD_m=SABST_O(i,n)*PAR_m ! en W/m2?
              pk490=Kd
              !Iz en W/m2 
              Iz=RAD_m*exp(-pk490*(-1*depth)) !depth passe en positif ! luminosite at depth
              g_coral = gmax_coral*tanh(Iz/Ik) !in m/day
!              print*, 'light limitation ',  gmax_coral, tanh(Iz/Ik),    &
!                       g_coral, Ik, Iz
              !g_coral = gmax_coral !test
              if (Iz.le.Imin) then
                 g_coral=0.0
              endif

!also limit in depth if less than 150m no more corals
!              print*, 'depth', depth
              if (depth.le.-150) then
                 g_coral=0.0
              endif
              


! limitation by temperature?
! my old version
!        epaisseur=0.004 !0.012 !0.004
!        x_opti=30 !25
!        y_opti=1
!        temp_factor=-epaisseur*(temp-x_opti)**2+y_opti
        !print*, 'temp_factor ', temp_factor
!        !temp_factor=1.0
!        if (temp_factor.lt.0) temp_factor=0.0

! Temperature dependency
!Manon s version
!              if ((temp.lt.24))then
!                g_coral = g_coral*coef_ftemp*(temp-18)
!              else if ((temp.lt.27).and.(temp.gt.24)) then
!                g_coral = g_coral*coef_ftemp*(24-18)
!              else if (temp.gt.27) then
!                g_coral = g_coral*(coef_ftemp*(temp-26)+coef_ftemp*(24-18))
!              endif

!new simpler version: linear
#if ( 1 )
                if (temp .lt. tfmin) then
                   temp_factor=0
                else if (temp .ge. tfmin .and. temp.le.tfmax) then
                  temp_factor=a_temp+b_temp*temp
                else
                  temp_factor=1
                endif
                g_coral = g_coral*temp_factor
#endif


! Limitation by recovering from bleaching 
#if ( 1 )
              IF (timebleach_l .NE. 0.) THEN
                g_coral = g_coral*(1.-EXP(-(NYR-timebleach_l)           &
                              /tau_bleach_l))
              ENDIF
#endif


! limitation by supersaturation
! Langdon & Atkinson, JGR, 2005
#if ( 1 ) 
              if (omega_arag.gt.1) then 
                g_coral = g_coral*(omega_arag-1.)/arag_k
              else 
                g_coral = 0.
              endif
              !print*, 'omega arag ', omega_arag, g_coral
#endif

! production in Pmol/day/grid cell 
! nb units: g_coral in m/day, area in m2 (*1e-6 to be in km2), coral_dens in kg/m3
! caco3_molar_mass in g/mol
! hence P_carb in Pmol/day
              P_carb = g_coral*topof                                    &
                            *area*1.0e-6                                &
                            *1.0e-6*coral_dens/caco3_molar_mass
!              print*, '*topof*area*coral_dens', P_carb
      endif

!Carbonate dissolution if omega < 1, sinon dissolution toute petite or no dissolution
!      arag_sat=5 ! test a commenter
!      if (arag_sat.lt.1) then
!           D_carb=mass_carb*(1.0-(arag_sat)) !ajout non linearité **2 ? dissolution complete si omega =0 ou pas (1.2 ou 1)?
!      else
!           D_carb=0.0 !0.001*mass_carb !0.0 ! 0 ou une toute petite valeur constante (genre 1percent) ?
!      endif

!      D_carb=0.0 !test

!nb    dissolution of all existing coral if no production
!       if (P_carb .eq. 0) then
!         D_carb=mass_carb/(caco3_molar_mass*1.0e15)
!       else
!         D_carb=0.0
!       endif

!nb    dissolution of half of all existing coral in 10 years if no production
!       if (P_carb .eq. 0) then ! si pas de production
!         if (prod_before .eq. 1) then ! si production avant
!            prod_before=0 ! plus de production
!            t0_diss=NYR-1 ! init time constante
!         endif
!         D_carb=mass_carb/(caco3_molar_mass*1.0e15)                     &
!               * exp(-lambda_diss*(NYR-t0_diss))
!       else
!         D_carb=0.0
!         prod_before=1
!       endif

! test : tout est dissous
!       D_carb=mass_carb/(caco3_molar_mass*1.0e15)

! For now no dissolution
      D_carb=0.0


! Net production in Pmol/day
      net_carb=P_carb-D_carb

!New mass in g/day
!      mass_carb_new=net_carb*caco3_molar_mass*1.0e15

!net mass
!      mass_carb_old=mass_carb
!      mass_carb=mass_carb_old+mass_carb_new


!      if (mass_carb_new .gt.0) then
!          print*, ' '
!          print*, 'in coral: prod, mass ',                             &
!          mass_carb_new, mass_carb
!          print*, 'light factor: ', tanh(Iz/Ik)
!          print*, 'omega arag factor: ', (omega_arag-1.)/arag_k
!      endif

      end subroutine corals
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      subroutine out_coral_global
! writing of global annual outputs

!nb write in Coral_output.txt
!      write (coral_res_fich,'(i6,5f14.5)')                              &
      write (coral_res_fich,'(i6,3f14.5,1f15.2)')                       &
         NYR, total_area_coral_an, total_prod_coral_an,                 &
         total_mass_coral_an
!         total_mass_coral_an, total_area_coral_an_40m,                  &
!         total_prod_coral_an_40m

      end subroutine
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!tbd moved in omega_mod
!       subroutine calc_omega(i,j,n)
!
!       use const_mod, only: gpes, rho0
!       use carbonate_speciation_mod, only: calc_co3sat_arag, incche
!       use loveclim_transfer_mod, only: TM, SM, mid_level
!       use marine_bio_mod, only: OCO3, OALK, ODIC
!       !use coral_mod, only: omega_arag3D
!
!       integer(kind=ip), intent(in) :: i,j,n
!
!       real(kind=dblp) :: p_bar
!       real(kind=dblp) :: sCO2, xpCO2, xCO2, xHCO3, xCO3, z_h
!       real(kind=dblp) :: omega_arag, CO3sat_ar
!
!!nb       computes CO32- concentration at saturation
!          ! needs temperature in Kelvin
!          ! pressure in bar=rho * g * z = Kg/m3 * N/kg * m = N/m2 = Pa
!          ! 1 bar = 1e5 Pa
!          ! rho = b * rho0 /g ?
!          ! en attendant use rho0
!          ! beware mid_level is negative...
!          !write(*,*) 'computation of p_bar', rho0, gpes,mid_level(j+1)
!          p_bar= rho0 * gpes * (-1) * mid_level(j+1)*1e-5
!          !p_bar=2.0
!          !write(*,*) 'p_bar', p_bar
!          CO3sat_ar=calc_co3sat_arag(TM(i,j,n)+273.15, SM(i,j,n), p_bar)
!          !CO3sat_ar=calc_co3sat_arag(30+273.15, 35.0, p_bar)
!          !CO3sat_ar=CO3sat_ar ! in mol/kg
!
!          !write(*,*) 'ODIC et OALK', ODIC(i,j,n), OALK(,j,n)
!          call incche(TM(i,j,n)+273.15,SM(i,j,n),p_bar,ODIC(i,j,n),     &
!          !call incche(30+273.15,35.0 ,p_bar,ODIC(i,j,n),                &
!            OALK(i,j,n),sCO2,xpCO2,xCO2,xHCO3,xCO3,z_h)
!
!          OCO3(i,j,n)=xCO3
!
!          omega_arag=OCO3(i,j,n)/CO3sat_ar
!          omega_arag3D(i,j,n)=omega_arag
!          !write(*,*),'j',j, OCO3(i,j,n)
!          !if (OCO3(i,j,n).gt.0) then
!          ! print*, 'omega_arag ', OCO3(i,j,n)*1e6, CO3sat_ar*1e6,
!          ! omega_arag
!          !endif
!
!          !nb keep ph=-log[h]
!          PH(i,j,n)=-log10(z_h)
!
!       end subroutine calc_omega
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      subroutine calc_monthly_temp(i,j,n)
! stores montly mean temperature
! over 30 years
! and update max value for grid i,j,n

!      use loveclim_transfer_mod, only: TM, KMON
      use loveclim_transfer_mod, only: KMON ! use temperature with added variability

      integer, intent(in) :: i,j,n
      REAL, dimension(nb_mois) :: temp_mois_temp

!      temp_mois(i,j,n)=temp_mois(i,j,n)+TM(i,j,n)
      temp_mois(i,j,n)=temp_mois(i,j,n)+TM_pluswkvar(i,j,n)
      if (KMON.eq.1) then ! si dernier jour du mois (jour 30)
          temp_mois(i,j,n)=temp_mois(i,j,n)/30. !monthly mean temperature
          ! shift all previous months and fill last one
          temp_mois_temp=temp_mois_all(i,j,n,:)
          temp_mois_all(i,j,n,1:indice_mois-1)=                         &
                    temp_mois_temp(2:indice_mois)
          temp_mois_all(i,j,n,indice_mois)=temp_mois(i,j,n)
          temp_mois(i,j,n)=0.0
          ! max of climatological monthly mean temperature
!          if (NYR.gt.window_MMM) then 
           MMMclim(i,j,n)=MAXVAL(temp_mois_all(i,j,n,:))
!          endif
          !write(*,*) 'MMMclim: ', MMMclim(i,j,n)
      endif


      end subroutine calc_monthly_temp
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      subroutine calc_DHW(i,j,n)
      ! computes degree heating weeks DHW (degree/week) following NOAA

!      use loveclim_transfer_mod, only: TM, KWEEK
!      use loveclim_transfer_mod, only: KWEEK

      integer, intent(in) :: i,j,n
      real HS
      real,dimension(nb_hs) :: xsHS_temp


      !first computes Hot Spot every day
!      if (TM(i,j,n).ge.MMMclim(i,j,n)) then
!        HS=TM(i,j,n)-MMMclim(i,j,n)
      if (TM_pluswkvar(i,j,n).ge.MMMclim(i,j,n)) then
        HS=TM_pluswkvar(i,j,n)-MMMclim(i,j,n)
      else
        HS=0.0
      endif
      !then computes excess Hot Spot if >1 degree and keep value in
      !matrix
      if (HS.ge.1) then
      !!!old version: tbd
        !xsHS(i,j,n)=xsHS(i,j,n)+HS
      !else
      !  xsHS(i,j,n)=0.0
      !!!end old version
        !first fill all values, then shift and replace last value
        if (indice_hs.lt.nb_hs) then !nb_hs=84 days
          xsHS_all(i,j,n,indice_HS)=HS
        else
          xsHS_temp(:)=xsHS_all(i,j,n,:)
          xsHS_all(i,j,n,1:indice_hs-1)=xsHS_temp(2:indice_hs)
          xsHS_all(i,j,n,indice_hs)=HS
        endif
      endif

!!! old version : tbd
!      if (KWEEK.eq.1) then ! if end of week
!        !average over the week to have the weekly excess hot spot
!        xsHS(i,j,n)=xsHS(i,j,n)/5 ! week of 5 days
!
!        !shifts and replaces last index value
!        xsHS_temp(:)=xsHS_all(i,j,n,:)
!        xsHS_all(i,j,n,1:indice_hs-1)=xsHS_temp(2:indice_hs)
!        xsHS_all(i,j,n,indice_hs)=xsHS(i,j,n)
!        !then computes Degree Heating Week=sum of excess
!        !temperature over 17 weeks of 5 days=85 days, equivalent to 12 weeks of 7 days=84days
!        DHW(i,j,n)=SUM(xsHS_all(i,j,n,:))
!        if (DHW(i,j,n).ge.dhw_thresh_moderat) then ! si plus que 4 degree for moderate bleaching
!          DHW_nb(i,j,n)=DHW(i,j,n)+1 ! compte le nombre de weeks qui declenche bleaching
!        endif
!        xsHS(i,j,n)=0.0 ! set back to 0
!        !write(*,*) 'DHW: ', DHW(i,j,n)
!      endif
!!! end old version


        !then computes Degree Heating Week=sum of excess
        !temperature over 84 days / 7 (in degree per week)
        DHW(i,j,n)=SUM(xsHS_all(i,j,n,:)/7.)
        if (DHW(i,j,n).ge.dhw_thresh_moderat) then ! si plus que 4 degree/week for moderate bleaching
          DHW_nb(i,j,n)=DHW_nb(i,j,n)+1 ! compte le nombre de weeks qui declenche bleaching
        endif
!tbd        xsHS(i,j,n)=0.0 ! set back to 0
      end subroutine calc_DHW
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      subroutine calc_bleach(i,j,n)

      !For bleaching, from Guy Munhoven based on NOAA
      !computes timebleach and tau_bleach used in corals

      integer, intent(in) :: i,j,n
     
                                    ! Degree-Heating Weeks control
                                    ! ----------------------------
          IF (DHW(i,j,n) .GE. dhw_thresh_strong) THEN
                                    ! Strong bleaching event detected!
            tau_bleach(i,j,n) = tau_bleach_strong
            timebleach(i,j,n) = NYR

          ELSEIF (DHW(i,j,n) .GE. dhw_thresh_moderat) THEN
                                    ! Moderate bleaching event detected
            IF (tau_bleach(i,j,n) .LT. tau_bleach_moderat) THEN
                                    ! Corals in healthy state (i.e., not
                                    ! recovering from a previous
                                    ! bleaching event)
                                    ! simply set the bleaching
                                    ! parameters
              tau_bleach(i,j,n) = tau_bleach_moderat
              timebleach(i,j,n) = NYR

            ELSEIF (tau_bleach(i,j,n) .LT. tau_bleach_strong) THEN
                                    ! Corals are recovering from a
                                    ! previous
                                    ! moderate event.
              IF ((NYR - timebleach(i,j,n)) .LE. 2.) THEN
                                    ! If this event was less than 2
                                    ! years ago,
                       
                tau_bleach(i,j,n) = tau_bleach_strong
                timebleach(i,j,n) = NYR

              ELSE
                                    ! This event was more than 2 years
                                    ! ago,
                                    ! leave in moderate state, but reset
                                    ! bleaching date.
                timebleach(i,j,n) = NYR

              ENDIF

            ELSE
                                    ! Corals are recovering from a
                                    ! previous
                                    ! strong event. Leave in recovery
                                    ! from
                                    ! strong but reset the bleaching
                                    ! date.
              timebleach(i,j,n) = NYR

            ENDIF

          ELSE
                                    ! No bleaching event triggered by
                                    ! DHW
                                    ! thresholds. Check if bleaching
                                    ! control
                                    ! can possibly be reset (i.e.,
                                    ! previous
                                    ! bleaching events were sufficiently
                                    ! long ago -- typically more than
                                    ! 4*tau).
            IF (timebleach(i,j,n) .NE. 0.) THEN
                                    ! Coral is currently in a recovery
                                    ! phase:
              IF ((NYR-timebleach(i,j,n)) .GT. 4.*tau_bleach(i,j,n))    &
                 THEN                  ! It was triggered more than 4*tau ago:
                timebleach(i,j,n) = 0.      ! Reset bleaching time
                tau_bleach(i,j,n) = 0.      ! Reset recovery time-scale
              ENDIF

            ENDIF

          ENDIF



                                    ! Report moderate and strong hot
                                    ! events,
                                    ! but limit print-out of hot events
                                    ! to
                                    ! temperatures above 25 degC to
                                    ! ignore
                                    ! numerous cold "hot events"
!          IF (TMMM(i,n).gt.25.) THEN
!            IF (DHW_curr(i,n).ge.dhw_thresh_strong) THEN
!              write(1000,*) 'Strong Hot event ', i, n, NYR, NJUL,
!     @                    DHW_curr(i,n)
!            ELSEIF (DHW_curr(i,n).ge.dhw_thresh_moderat) THEN
!              write(1000,*) 'Moderate Hot event ', i, n, NYR, NJUL,
!     @                    DHW_curr(i,n)
!            ELSE
!              CONTINUE
!            ENDIF
!          ENDIF



      end subroutine calc_bleach
!----------------------------------------------------------------------


!----------------------------------------------------------------------
      subroutine calc_temp_variability(i,j,n)
!from Guy Munhoven (in CLIMBER)
! Add parameterized weekly variability to the ocean surface temperature

      use loveclim_transfer_mod, only: TM

      integer, intent(in) :: i,j,n

      INTEGER, ALLOCATABLE, DIMENSION(:) :: iseed
      INTEGER iseed_size, i4seed
      INTEGER dt_values(8)
      DOUBLE PRECISION, SAVE :: gaus1, gaus2, gaus
      DOUBLE PRECISION       :: rand1, rand2
      DOUBLE PRECISION :: DT_rand

            IF (k_mbiota_rand .EQ. 0) THEN

                                    ! If the random number generator has
                                    ! never been called before: carry out
                                    ! standard conforming initialisation
                                    ! of the Fortran 90 random number
                                    ! generator.

                                    ! 1. Retrieve the size of the seed
              CALL RANDOM_SEED(SIZE=iseed_size)

                                    ! 2. Allocate memory for the seed
              ALLOCATE(iseed(iseed_size))

                                    ! 3. Generate a seed based upon the summary
                                    !    from DATE_AND_TIME

              DO i4seed = 1, iseed_size
                CALL DATE_AND_TIME(VALUES=dt_values)
                iseed(i4seed) = SUM(dt_values)
              ENDDO

                                    ! 4. Initialise the random number generator
                                    !    with the new seed

              CALL RANDOM_SEED(PUT=iseed)

                                    ! 5. clean up and set "ready-for-use flag"

              DEALLOCATE(iseed)

              k_mbiota_rand = 1

                                    ! 6. Set DT_rand to zero to start the AR(1)
                                    !    process below
              DT_rand = 0.


            ENDIF


                                    ! Generates a random number from a
                                    ! gaussian distribution with mean=0
                                    ! and variance=1.
                                    ! [GM???] corrected, needs *two* uniformly
                                    ! [GM???] distributed random variates

            IF (k_mbiota_rand .eq. 1) THEN
                                    ! If no unused random variate is in
                                    ! storage, generate two new ones with
                                    ! the Box-Muller transform.

              CALL RANDOM_NUMBER(rand1)
              CALL RANDOM_NUMBER(rand2)

              gaus1 = SQRT(-2.0D+00*LOG(rand1))                         &
                      *COS(2.0D+00*3.141592653589793D+00*rand2)

              gaus2 = SQRT(-2.0D+00*LOG(rand1))                         &
                      *SIN(2.0D+00*3.141592653589793D+00*rand2)

                                    ! Attribute the first one to 'gaus'
                                    ! and store the unused second one
              gaus = gaus1
              k_mbiota_rand = 2     ! Flag that there is one random number stored


            ELSE
                                    ! If there is an unused random number
              gaus = gaus2          ! stored, assign that one to 'gaus'
              k_mbiota_rand = 1     ! and flag that none is stored anymore.

            ENDIF
                                    ! Produce an AR(1) variate with
                                    ! - an auto-correlation parameter of 0.6
                                    ! - standard deviation of 0.3
                                    ! to add weekly variability to the otherwise
                                    ! smooth temperature evolution in the model.

!cnb ori from climber version            DT_rand=DT_rand*0.6+0.3*gaus
            DT_rand=DT_rand*0.8964+0.2758*gaus ! from recalculation of Guy

!            print *,'random',rand1,rand2,gaus,DT_rand


            TM_pluswkvar(i,j,n) = TM(i,j,n) + DT_rand

      end subroutine calc_temp_variability
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      subroutine calc_nino3()

      use loveclim_transfer_mod, only: TM, KMON, SQRO2

      integer i,j,n
      real nino3_temp
      real nino3_var_temp
      real sum_area

      j=1 ! surface

! nino3 = montly mean surface temperature averaged over 5S-5N
! and 150W - 90 W (soit 210E - 270E)
      do n=62,82 !NOC_CBR n=62-> 211.5E, n=82 ->271.5E
        do i=26,29 !LT  i=26->-4.5, i=29-> 4.5N
          if (MGT(i,j,n).eq.1) then
            nino3_temp=nino3_temp+TM(i,j,n)*SQRO2(i,n)
            sum_area=sum_area+SQRO2(i,n)
            nino3_var_temp=nino3_var_temp+TM_pluswkvar(i,j,n)*SQRO2(i,n)
          endif
        enddo
      enddo

      !area mean
      nino3=nino3+nino3_temp/sum_area
      nino3_var=nino3_var+nino3_var_temp/sum_area
      sum_area=0

      !write(*,*) 'nino3 day', nino3, sum_area

      if (KMON.eq.1) then ! si dernier jour du mois (jour 30)
          nino3=nino3/30. !monthly mean temperature
          nino3_var=nino3_var/30. !monthly mean temperature

          ! write output
          !write(*,*) 'nino3', nino3, nino3_var
          write (c_nino_fich,'(i6, 2F8.3)')                             &
          NYR, nino3, nino3_var

          nino3=0.0
          nino3_var=0.0
          sum_area=0.0
      endif

      end subroutine calc_nino3
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine restart_coral(choix)

      !use loveclim_transfer_mod, only: TM, KMON, SQRO2

       INTEGER :: choix

       INTEGER :: fich_num, nrecl
       CHARACTER*14, PARAMETER :: fich_res_name="rest_coral.dat"
       CHARACTER*24, PARAMETER ::                                       &
                         fich_res_name_old="startdata/rest_coral.dat"
       LOGICAL :: existe

! Find free number
       fich_num=298
       print*, "init_cor_fich", fich_num
       existe=.TRUE.

       DO WHILE (existe)
        INQUIRE(fich_num,OPENED=existe)
        print*, "init_cor_fich", fich_num, existe
        fich_num=fich_num+1 
       END DO 

! If write restart
!-----------------
       IF (choix.EQ.0) THEN
!        Size of written variables
         nrecl = 1*SIZE(coral_mass_subgrid)
         nrecl = nrecl * KIND(coral_mass_subgrid)

!        write restart
         OPEN(UNIT=fich_num, FILE=fich_res_name, STATUS='unknown',      &
              ACCESS='direct', RECL=nrecl, ACTION='write')

         WRITE(UNIT=fich_num, REC=1) coral_mass_subgrid

         CLOSE(UNIT=fich_num)

! If read restart
!-----------------
       ELSE IF (choix.EQ.1) THEN
!        Size of written variables
         nrecl = 1*SIZE(coral_mass_subgrid)
         nrecl = nrecl * KIND(coral_mass_subgrid)

!       read restart
        WRITE(*,*) "Lecture a partir du fichier:"
        WRITE(*,*) fich_res_name_old

        OPEN(UNIT=fich_num, FILE=fich_res_name_old, STATUS='old',       &
              ACCESS='direct', RECL=nrecl, ACTION='read')

        READ(UNIT=fich_num,REC=1) coral_mass_subgrid  

        CLOSE(UNIT=fich_num)

       !test
       write(*,*) 'coral_mass_subgrid', coral_mass_subgrid(:,:,:)
       coral_mass_subgrid(:,:,:)=coral_mass_subgrid(:,:,:)*10000
       write(*,*) 'coral_mass_subgrid', coral_mass_subgrid(:,:,:)

       ENDIF

      end subroutine restart_coral
!----------------------------------------------------------------------

      end module coral_mod
#endif



