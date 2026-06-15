!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!      write_isoatm_mod : writes isotopic atmospheric restart file
!
!      Auteur : Didier M. Roche, Thibaut Caley
!      Date   : 02 mai 2012
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

module write_isoatm_mod

  use global_constants_mod, only: dblp=>dp, ip

  implicit none
  private

  public :: write_isoatm

contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
  subroutine write_isoatm
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

#if ( ISOATM >= 1 )

    use iso_param_mod, only: ieau, neauiso
    use comphys,       only: rmoisg, torain, tosnow

    implicit none

    integer(kind=ip) :: wisoatm_restartdat_id

    open(newunit=wisoatm_restartdat_id, file='wisoatm_restart.dat', &
         status='unknown', form='unformatted')

    write(wisoatm_restartdat_id) rmoisg(:,:,ieau+1:neauiso), &
                                 torain(:,:,ieau+1:neauiso), &
                                 tosnow(:,:,ieau+1:neauiso)

    close(wisoatm_restartdat_id)

#endif

  end subroutine write_isoatm

end module write_isoatm_mod
