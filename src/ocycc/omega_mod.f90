! Omega computation
! Written by N. Bouttes, G. Munhoven...
!
!----------------------------------------------------------------------
#include "choixcomposantes.h"

      module omega_mod

! use section
      use declars_mod
      use global_constants_mod, only: dblp=>dp, ip

! declaration section
       implicit none

!input output
      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: omega_arag3D
      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: omega_calc3D
      REAL(kind=8), dimension(LT,JT,NOC_CBR), target :: PH

      contains
!----------------------------------------------------------------------
       subroutine calc_omega_ar(i,j,n)

       use const_mod, only: gpes, rho0
       use carbonate_speciation_mod, only: calc_co3sat_arag, incche
       use loveclim_transfer_mod, only: TM, SM, mid_level
       use marine_bio_mod, only: OCO3, OALK, ODIC
       !use coral_mod, only: omega_arag3D

       integer(kind=ip), intent(in) :: i,j,n

       real(kind=dblp) :: p_bar
       real(kind=dblp) :: sCO2, xpCO2, xCO2, xHCO3, xCO3, z_h
       real(kind=dblp) :: omega_arag, CO3sat_ar

!nb       computes CO32- concentration at saturation
          ! needs temperature in Kelvin
          ! pressure in bar=rho * g * z = Kg/m3 * N/kg * m = N/m2 = Pa
          ! 1 bar = 1e5 Pa
          ! rho = b * rho0 /g ?
          ! en attendant use rho0
          ! beware mid_level is negative...
          !write(*,*) 'computation of p_bar', rho0, gpes,mid_level(j+1)
          p_bar= rho0 * gpes * (-1) * mid_level(j+1)*1e-5
          !p_bar=2.0
          !write(*,*) 'p_bar', p_bar
          CO3sat_ar=calc_co3sat_arag(TM(i,j,n)+273.15, SM(i,j,n), p_bar)
          !CO3sat_ar=calc_co3sat_arag(30+273.15, 35.0, p_bar)
          !CO3sat_ar=CO3sat_ar ! in mol/kg

          !write(*,*) 'ODIC et OALK', ODIC(i,j,n), OALK(,j,n)
          call incche(TM(i,j,n)+273.15,SM(i,j,n),p_bar,ODIC(i,j,n),     &
          !call incche(30+273.15,35.0 ,p_bar,ODIC(i,j,n),               &
            OALK(i,j,n),sCO2,xpCO2,xCO2,xHCO3,xCO3,z_h)

          OCO3(i,j,n)=xCO3

          omega_arag=OCO3(i,j,n)/CO3sat_ar
          omega_arag3D(i,j,n)=omega_arag
          !write(*,*),'j',j, OCO3(i,j,n)
          !if (OCO3(i,j,n).gt.0) then
          ! print*, 'omega_arag ', OCO3(i,j,n)*1e6, CO3sat_ar*1e6,
          ! omega_arag
          !endif

          !nb keep ph=-log[h]
          PH(i,j,n)=-log10(z_h)

       end subroutine calc_omega_ar
!----------------------------------------------------------------------
!----------------------------------------------------------------------
       subroutine calc_omega_ca(i,j,n)

       use const_mod, only: gpes, rho0
       use carbonate_speciation_mod, only: calc_co3sat_calc, incche
       use loveclim_transfer_mod, only: TM, SM, mid_level
       use marine_bio_mod, only: OCO3, OALK, ODIC
       !use coral_mod, only: omega_arag3D

       integer(kind=ip), intent(in) :: i,j,n

       real(kind=dblp) :: p_bar
       real(kind=dblp) :: sCO2, xpCO2, xCO2, xHCO3, xCO3, z_h
       real(kind=dblp) :: omega_calc, CO3sat_ca

!nb       computes CO32- concentration at saturation
          ! needs temperature in Kelvin
          ! pressure in bar=rho * g * z = Kg/m3 * N/kg * m = N/m2 = Pa
          ! 1 bar = 1e5 Pa
          ! rho = b * rho0 /g ?
          ! en attendant use rho0
          ! beware mid_level is negative...
          !write(*,*) 'computation of p_bar', rho0, gpes,mid_level(j+1)
          p_bar= rho0 * gpes * (-1) * mid_level(j+1)*1e-5
          !p_bar=2.0
          !write(*,*) 'p_bar', p_bar
          CO3sat_ca=calc_co3sat_calc(TM(i,j,n)+273.15, SM(i,j,n), p_bar)
          !CO3sat_ar=calc_co3sat_arag(30+273.15, 35.0, p_bar)
          !CO3sat_ar=CO3sat_ar ! in mol/kg

          !write(*,*) 'ODIC et OALK', ODIC(i,j,n), OALK(,j,n)
          call incche(TM(i,j,n)+273.15,SM(i,j,n),p_bar,ODIC(i,j,n),     &
          !call incche(30+273.15,35.0 ,p_bar,ODIC(i,j,n),               &
            OALK(i,j,n),sCO2,xpCO2,xCO2,xHCO3,xCO3,z_h)

          OCO3(i,j,n)=xCO3

          omega_calc=OCO3(i,j,n)/CO3sat_ca
          omega_calc3D(i,j,n)=omega_calc
          !write(*,*),'j',j, OCO3(i,j,n)
          !if (OCO3(i,j,n).gt.0) then
          ! print*, 'omega_arag ', OCO3(i,j,n)*1e6, CO3sat_ar*1e6,
          ! omega_arag
          !endif

          !nb keep ph=-log[h]
          PH(i,j,n)=-log10(z_h)

       end subroutine calc_omega_ca
!----------------------------------------------------------------------

                                                                                 
      end module omega_mod
