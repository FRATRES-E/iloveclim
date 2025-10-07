! Coastal module
!
! Written by C. Contoux, N. Bouttes
!
!----------------------------------------------------------------------
#include "choixcomposantes.h"

#if (COASTAL == 1 )
      module coastal_mod

! use section
      use declars_mod
      use loveclim_transfer_mod, ONLY: ZX, TYER, TDAY,                  &
            SQRO2, MGT, NYR
      use C_res_mod, ONLY: coastal_res_fich

! declaration section
       implicit none

!input output coastal
!      INTEGER level
!      REAL hypso_coast

!local pas local
      INTEGER, PARAMETER :: kmax_hypso = 299
!      REAL bati
!      INTEGER mask !0=ocean, 1=land
!      INTEGER t,i,j,k
!      REAL, dimension(LT, JT, NOC_CBR) :: P_carb_an
      !REAL mass_carb_old
!      INTEGER compteur !jour de l anne (de 1 a 365)
      REAL, dimension(LT, kmax_hypso, NOC_CBR) :: hypso_coast
!nb not needed more      REAL, dimension(LT, kmax_hypso, NOC_CBR-2) :: hypso_coast_read
      REAL, dimension(kmax_hypso) :: level_hypso
!      REAL, dimension(JT) :: level
      REAL, dimension(LT) :: latitude
      REAL, dimension(NOC_CBR) ::  longitude
      character*256 testchar
      REAL, dimension(LT, JT, NOC_CBR) ::  surface_coast
!      INTEGER d_hypso
!      REAL level_hypso_inf, level_hypso_sup
!      REAL, dimension(NOC_CBR,LT) :: topof

!      INTEGER year
!      INTEGER i_sl
!      REAL dum
!      INTEGER, PARAMETER :: t_sl=15
!      REAL, dimension(t_sl) :: sea_level

!      !REAL depth, depth_inf, depth_sup
!      REAL, dimension(LT, JT, NOC_CBR) ::  area_3D
!      INTEGER kk

!     REAl par
!      !REAL, dimension() par(jmax,imax)
!      !REAL Kd_490
!      REAL, dimension(NOC_CBR, LT) ::  Kd_490 !(jmax, imax)
!      REAL PAR_m
!      REAL dens_ocn
!      REAL water_z_pp
!
!      REAL weathering_odic
!      REAL weathering_oalk
!      REAL surface_ocean
!
!      REAL, dimension(JT) :: CO3sat_ar
!      REAL, dimension(JT) :: CO3sat_ca

!output
!      REAL total_area_coral_an
!      REAL total_prod_coral_an
!      REAL total_mass_coral_an
!      REAL, dimension(LT,JT,NOC_CBR) :: coral_prod
!      REAL, dimension(LT,JT,NOC_CBR) :: coral_area
!      REAL, dimension(LT,JT,NOC_CBR) :: coral_cum_mass
!      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: coral_prod_out
!      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: coral_area_out

!nb comprend pas pourquoi ca n est pas de la dimension de scal...
!      REAL, dimension(LT,JT,NOC_CBR) :: corproda, corareaa
!
!      REAL coral_CO2
!
!      character*256 outfilename



      contains

! fonctions and subroutines
!----------------------------------------------------------------------
! initialisation
      subroutine ini_coastal

      use ncio,      only: nc_read

! local
      integer :: n,m
      integer :: i,j
      REAL :: thrld
      REAL :: somme_surface

!area (surface and hypsometry) from gebco
      thrld=-150. ! coastal ocean is defined above 150m

      print*, 'read hypsometry cotier'
      !!call nc_read_attr("GEBCO/hypsometry.nc", "title", testchar)
      !!write(*,*) "Title: ", trim(testchar)
      call nc_read("inputdata/hypsometry_cotier.nc","hypso_coast",      &
!tbd      hypso_coast_read) !(:,:,:)
      hypso_coast) !(:,:,:)
!tbd      print*, 'hypso size', ubound(hypso_coast_read)
      print*, 'hypso size', ubound(hypso_coast)
      call nc_read("inputdata/hypsometry_cotier.nc","level",level_hypso)
      print*,'level_hypso', level_hypso
      !nb tableau a bord replie :valeur(1)=valeur(nmax-1) et
      !valeur(nmax)=valeur(2)
!nb below not needed anymore tbd
!      do n=2,NOC_CBR-1
!        hypso_coast(:,:,n)=hypso_coast_read(:,:,n-1)
!      enddo
!      hypso_coast(:,:,1)=hypso_coast(:,:,NOC_CBR-1)
!      hypso_coast(:,:,NOC_CBR)=hypso_coast(:,:,2)


     ! 1rst loop on lon lat JT to read coastal surface
      surface_coast(:,:,:)=0.
      somme_surface=0.
      do j=1,JT 
       do n=1,NOC_CBR
        do i=1,LT
          if (MGT(i,j,n).eq.1) then ! pour tester si on est ds locean
           do m = 1,kmax_hypso
            if ((level_hypso(m) .LE. ZX(j)*(-1)) .AND. (level_hypso(m)  &
             .GT. ZX(j+1)*(-1)) .AND. (level_hypso(m) .GT. thrld)) then
             surface_coast(i,j,n)=surface_coast(i,j,n)                  &
             +hypso_coast(i,m,n)
             if (hypso_coast(i,m,n) .GT. 0) then
              write(*,*) 'level_hypso', level_hypso(m), 'ZX(j)*(-1)',   &
              ZX(j)*(-1), 'ZX(j+1)*(-1)', ZX(j+1)*(-1)
              write(*,*) 'surface_coast', surface_coast(i,j,n)
              !somme_surface=somme_surface+surface_coast(i,j,n)
              !write(*,*) 'somme_surface =', somme_surface
             endif
            endif 
           enddo
           somme_surface=somme_surface+surface_coast(i,j,n)
           write(*,*) 'somme_surface =', somme_surface
          endif
        enddo
       enddo
      enddo 

      end subroutine ini_coastal

!----------------------------------------------------------------------

!----------------------------------------------------------------------
      end module coastal_mod
#endif
