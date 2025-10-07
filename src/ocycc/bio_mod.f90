! module de declaration de variables
! bio
! remplace bio.inc
! C. Dumas 29-01-2007
!  Verifie pour C2.5 07-12-2007 C. Dumas & D. Roche

module bio_mod

! USE declar_mod

implicit none

!dmr --- Pour test de compilation
INTEGER, PARAMETER :: IT=32
INTEGER, PARAMETER :: NS=64
!dmr --- Pour test de compilation

!dmr incoherence en couple LCM - VECCARB avec veget.h

!dmr %%% REAL :: a,bet,gamm,fmax,avecube,tmin,npp,nppmax,v1,v2,v3
!dmr %%% REAL :: c1t,c2t,c3t,c1g,c2g,c3g
!dmr %%% REAL :: d1t,d2t,d3t,d1g,d2g,d3g
!dmr %%% REAL :: e1t,e2t,e3t,e1g,e2g,e3g
!dmr %%% REAL :: f1t,f2t,f3t,f1g,f2g,f3g
!dmr %%% REAL :: k1t,k2t,k3t,k1g,k2g,k3g
!dmr %%% REAL :: t1t,t2t,t3t,t4t,t1g,t2g,t3g,t4g
!dmr %%% REAL :: ps1,ps2,ps3,ps4,ps5,soilt
!dmr %%% REAL :: forshare_st,t1tn,t1td, desshare_st,nlshare_st,
!dmr %%% REAL :: g4share_st
!dmr %%% REAL :: deng,dentd,dentn,laig,lait
!dmr %%% REAL :: ave_t,ave_pr,ave_pr05,desmin,desmax
!dmr %%% REAL :: ades, acr, k0t, k0g, k4g
!dmr %%% REAL :: c14init,c14tdec,c13init,c13frac,c13ratio,c13frac4,c14ratio
REAL :: lon, lat
! ajout C2.5 dmr & dc
!dmr %%% REAL :: gdd0, gdd0_min,gdd0_max,
REAL :: co2,co2_max
REAL :: tempor1_nb, tempor2_nb
!dmr&vm %%% INTEGER :: INI_STEP

!dmr %%% REAL, dimension(IT,NS) :: b1t, b2t, b3t, b4t
!dmr %%% REAL, dimension(IT,NS) :: b1g, b2g, b3g, b4g
!dmr %%% REAL, dimension(IT,NS) :: b1t14, b2t14, b3t14, b4t14
!dmr %%% REAL, dimension(IT,NS) :: b1g14, b2g14, b3g14, b4g14
!dmr %%% REAL, dimension(IT,NS) :: b1t13, b2t13, b3t13, b4t13
!dmr %%% REAL, dimension(IT,NS) :: b1g13, b2g13, b3g13, b4g13
REAL, dimension(IT,NS) :: carea, init_flag

! ajout reprise C2.5C dmr & nb
INTEGER :: test_veget

end module bio_mod
