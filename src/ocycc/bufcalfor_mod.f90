! module de declaration de variables
! calottes forcees
! remplace buffercalottesforcees.inc
! N. Bouttes 13.03.2008
module bufcalfor_mod

! USE declar_mod

implicit none
!dmr --- Pour essai de compilation
INTEGER, PARAMETER :: IT=32
INTEGER, PARAMETER :: NS=64
!dmr --- Pour essai de compilation

! NFT est le nombre de lignes dans le fichier de niveau marin
INTEGER, PARAMETER :: NFT=287
REAL, dimension(IT,NS) :: FRGLC_ORIG, HORO_ORIG, HOROCTL, FRGLC_CTL
REAL, dimension(IT,NS) :: FROCN_CTL,NOA_CTL
REAL, dimension(NFT) :: SEALEVLU, TDATE
REAL :: sealevel,slconv,sealevel_prev

end module bufcalfor_mod
