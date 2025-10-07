! Declaration of sediment parameters
! replaces parameter.sediment.include
! N. Bouttes 17.02.2009

module parameter_sediment_mod

!USE declar_mod

implicit none

!***************  Sediments *******************************************
!INTEGER, PARAMETER :: ip_sedmain_max=LT*JT*NOC
INTEGER, PARAMETER :: ip_sedmain_max=72*20*3
INTEGER, PARAMETER :: nz_sedmain_max=8
INTEGER, PARAMETER :: i_bur_rec_max=20

end module parameter_sediment_mod

