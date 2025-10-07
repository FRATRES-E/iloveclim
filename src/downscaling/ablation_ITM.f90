!-----|--1--------2---------3---------4---------5---------6---------7-|
!      This module implements the Insolation Temperature Melt (ITM)
!       parametrisation, here taken from Robinson, 2011, PhD Dissert.
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche, Emiel Spanier
!      Date   : 04 mars 2014
!      Derniere modification : 04 mars 2014
!-----|--1--------2---------3---------4---------5---------6---------7-|

       MODULE ablation_ITM


       IMPLICIT NONE

! --- dmr&es Shortwave radiation and sensible heat flux constant
! ---         Values can be taken between -40 and -60 W.m^-2
! ---        Robinson (2011) takes it at -55 to be close to the PDD
       REAL, PRIVATE, PARAMETER :: c_rad_heat_cst = -55.0

! --- dmr&es Latent heat of ice melting (J.kg^-1)
       REAL, PRIVATE, PARAMETER :: l_ice_melt = 334.0E3

! --- dmr&es Longwave radiation coefficient (W.m^-2.K^-1)
       REAL, PRIVATE, PARAMETER :: lwa_rad_coef = 10.0

! --- dmr&es Density of water (kg.m^-3)
       REAL, PRIVATE, PARAMETER :: wat_dens = 1000.0

! --- dmr&es Day length in seconds
       REAL, PRIVATE, PARAMETER :: day_length = 24.*60.*60.

! --- dmr&es Constant used in the computation afterwards ...
       REAL, PRIVATE, PARAMETER :: const_pre = day_length/(wat_dens*l_ice_melt)

       CONTAINS

!-----|--1--------2---------3---------4---------5---------6---------7-|
! --- dmr&es Computes transmissivity to yield surface shortwave from top
!      of the atmosphere shortwave incoming flux. Computed in Robinson (2011)
!      from a regression on Greenland data
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE transmissivity_Greenland(transmit,nx,ny    &
                                                   ,surf_height)

!$ USE OMP_LIB

       INTEGER, INTENT(IN) :: nx, ny

! --- dmr&es Surface height (m)
       REAL,DIMENSION(nx,ny),INTENT(IN) :: surf_height

! --- dmr&es Transmissivity ()
       REAL,DIMENSION(nx,ny),INTENT(OUT) :: transmit


       REAL, PARAMETER :: intercept=0.46, slope=6.E-5

       INTEGER :: i,j

!$OMP PARALLEL WORKSHARE
       FORALL (i=1:nx,j=1:ny)
         transmit(i,j) = intercept + slope*surf_height(i,j)
       END FORALL
!$OMP END PARALLEL WORKSHARE       

       RETURN

       END SUBROUTINE transmissivity_Greenland

!-----|--1--------2---------3---------4---------5---------6---------7-|
!-----|--1--------2---------3---------4---------5---------6---------7-|
       SUBROUTINE Potential_daily_melt(melt_rate,shortw_surf,surf_air_temp)

!$ USE OMP_LIB

! --- dmr&es Shortwave surface radiation (W.m^-2), Surface air
!             temperature (°C)
!     NOTA BENE: shortw_surf already includes the albedo at surface
!     effect, so it is a NET shortwave surface flux
       REAL,DIMENSION(:,:,:),INTENT(IN) :: shortw_surf, surf_air_temp
! --- dmr&es Melt rate (m.day^-1)
       REAL,DIMENSION(:,:,:),INTENT(OUT) :: melt_rate

! --- dmr&es Local variables
       INTEGER :: i,j,k

! --- dmr&es ITM formula implementation
!$OMP PARALLEL WORKSHARE
       FORALL (i=LBOUND(melt_rate,dim=1):UBOUND(melt_rate,dim=1)        &
              ,j=LBOUND(melt_rate,dim=2):UBOUND(melt_rate,dim=2)        &
              ,k=LBOUND(melt_rate,dim=3):UBOUND(melt_rate,dim=3))
         melt_rate(i,j,k) = max(const_pre*(shortw_surf(i,j,k)           &
                     +c_rad_heat_cst+lwa_rad_coef*surf_air_temp(i,j,k)) &
                            ,0.0d0)
       END FORALL
!$OMP END PARALLEL WORKSHARE      

       RETURN

       END SUBROUTINE Potential_daily_melt

!-----|--1--------2---------3---------4---------5---------6---------7-|
!-----|--1--------2---------3---------4---------5---------6---------7-|
       SUBROUTINE Potential_daily_melt_biascorr(melt_rate,shortw_surf,surf_air_temp,tannbias, biascoef)

!$ USE OMP_LIB

! --- dmr&es Shortwave surface radiation (W.m^-2), Surface air
!             temperature (°C)
!     NOTA BENE: shortw_surf already includes the albedo at surface
!     effect, so it is a NET shortwave surface flux
       REAL,DIMENSION(:,:,:),INTENT(IN) :: shortw_surf, surf_air_temp   &
                                          , tannbias
       REAL,                 INTENT(IN) :: biascoef
! --- dmr&es Melt rate (m.day^-1)
       REAL,DIMENSION(:,:,:),INTENT(OUT) :: melt_rate

! --- dmr&es Local variables
       INTEGER :: i,j,k

! --- dmr&es ITM formula implementation
!$OMP PARALLEL WORKSHARE
       FORALL (i=LBOUND(melt_rate,dim=1):UBOUND(melt_rate,dim=1)        &
              ,j=LBOUND(melt_rate,dim=2):UBOUND(melt_rate,dim=2)        &
              ,k=LBOUND(melt_rate,dim=3):UBOUND(melt_rate,dim=3))
         melt_rate(i,j,k) = max(const_pre*(shortw_surf(i,j,k)           &
                     +c_rad_heat_cst*(1.+tannbias(i,j,k)*biascoef)      &
                     +lwa_rad_coef*surf_air_temp(i,j,k))                &
                            ,0.0d0)
       END FORALL
!$OMP END PARALLEL WORKSHARE      

       RETURN

       END SUBROUTINE Potential_daily_melt_biascorr

!-----|--1--------2---------3---------4---------5---------6---------7-|
!-----|--1--------2---------3---------4---------5---------6---------7-|
       SUBROUTINE top_to_surf_shortwave(surf_swave,nx,ny,top_swave      &
                                       ,sheight)

!$ USE OMP_LIB

       INTEGER, INTENT(IN) :: nx, ny
! --- dmr&es Surface height (m), top of the atmosphere shortwave (W.m^-2)
       REAL,DIMENSION(nx,ny),INTENT(IN) :: sheight, top_swave

! --- dmr&es Surface incoming shortwave (W.m^-2)
       REAL,DIMENSION(nx,ny),INTENT(OUT) :: surf_swave


! --- dmr&es Local variables
       REAL, DIMENSION(nx,ny) :: trans
       INTEGER :: i,j

       CALL transmissivity_Greenland(trans,nx,ny,sheight)
       
!$OMP PARALLEL WORKSHARE
       FORALL (i=1:nx,j=1:ny)
         surf_swave(i,j) = top_swave(i,j)*trans(i,j)
       END FORALL
!$OMP END PARALLEL WORKSHARE

       RETURN

       END SUBROUTINE top_to_surf_shortwave

       END MODULE ablation_ITM
