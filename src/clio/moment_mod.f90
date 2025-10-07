!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module est un portage en FORTRAN 90 du module CLIO initial
!       cree en FORTRAN 77.
!      Il contient une routine suplementaire pour simuler les
!      equivalences du code FORTRAN 77 via des POINTER en FORTRAN 90
!      (dans l'environnement logiciel LUDUS)
!
!      Auteur : ??, Didier M. Roche
!      Date   : ??
!      Derniere modification : 19 Aout 2014
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module moment_mod

      use para0_mod, only: imax, jmax

      implicit none


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       vicom was initially designed as an equivalence with the variables below. It is now replaced by pointers / targets.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(imax,jmax,35) :: vicmom

      real(kind=8), dimension(imax,jmax)    :: sxg, syg, sxxg, syyg, sxyg, sxn, syn, sxxn, syyn, sxyn, sxa, sya, sxxa, syya, sxya  &
                                             , sxc0, syc0, sxxc0, syyc0, sxyc0, sxc1, syc1, sxxc1, syyc1, sxyc1, sxc2, syc2, sxxc2 &
                                             , syyc2, sxyc2, sxst, syst, sxxst, syyst, sxyst
      contains


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Procedure copy_from_vicmon is designed to copy vicmon data into the appropriate data structures
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine copy_from_vicmon()

        sxg(:,:)   = vicmom(:,:,1)
        syg(:,:)   = vicmom(:,:,2)
        sxxg(:,:)  = vicmom(:,:,3)
        syyg(:,:)  = vicmom(:,:,4)
        sxyg(:,:)  = vicmom(:,:,5)
        sxn(:,:)   = vicmom(:,:,6)
        syn(:,:)   = vicmom(:,:,7)
        sxxn(:,:)  = vicmom(:,:,8)
        syyn(:,:)  = vicmom(:,:,9)
        sxyn(:,:)  = vicmom(:,:,10)
        sxa(:,:)   = vicmom(:,:,11)
        sya(:,:)   = vicmom(:,:,12)
        sxxa(:,:)  = vicmom(:,:,13)
        syya(:,:)  = vicmom(:,:,14)
        sxya(:,:)  = vicmom(:,:,15)
        sxc0(:,:)  = vicmom(:,:,16)
        syc0(:,:)  = vicmom(:,:,17)
        sxxc0(:,:) = vicmom(:,:,18)
        syyc0(:,:) = vicmom(:,:,19)
        sxyc0(:,:) = vicmom(:,:,20)
        sxc1(:,:)  = vicmom(:,:,21)
        syc1(:,:)  = vicmom(:,:,22)
        sxxc1(:,:) = vicmom(:,:,23)
        syyc1(:,:) = vicmom(:,:,24)
        sxyc1(:,:) = vicmom(:,:,25)
        sxc2(:,:)  = vicmom(:,:,26)
        syc2(:,:)  = vicmom(:,:,27)
        sxxc2(:,:) = vicmom(:,:,28)
        syyc2(:,:) = vicmom(:,:,29)
        sxyc2(:,:) = vicmom(:,:,30)
        sxst(:,:)  = vicmom(:,:,31)
        syst(:,:)  = vicmom(:,:,32)
        sxxst(:,:) = vicmom(:,:,33)
        syyst(:,:) = vicmom(:,:,34)
        sxyst(:,:) = vicmom(:,:,35)

        return

      end subroutine copy_from_vicmon

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Procedure copy_to_vicmon is designed to copy appropriate data structures into the vicmon container
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine copy_to_vicmon()

        vicmom(:,:,1)  = sxg(:,:)
        vicmom(:,:,2)  = syg(:,:)
        vicmom(:,:,3)  = sxxg(:,:)
        vicmom(:,:,4)  = syyg(:,:)
        vicmom(:,:,5)  = sxyg(:,:)
        vicmom(:,:,6)  = sxn(:,:)
        vicmom(:,:,7)  = syn(:,:)
        vicmom(:,:,8)  = sxxn(:,:)
        vicmom(:,:,9)  = syyn(:,:)
        vicmom(:,:,10) = sxyn(:,:)
        vicmom(:,:,11) = sxa(:,:)
        vicmom(:,:,12) = sya(:,:)
        vicmom(:,:,13) = sxxa(:,:)
        vicmom(:,:,14) = syya(:,:)
        vicmom(:,:,15) = sxya(:,:)
        vicmom(:,:,16) = sxc0(:,:)
        vicmom(:,:,17) = syc0(:,:)
        vicmom(:,:,18) = sxxc0(:,:)
        vicmom(:,:,19) = syyc0(:,:)
        vicmom(:,:,20) = sxyc0(:,:)
        vicmom(:,:,21) = sxc1(:,:)
        vicmom(:,:,22) = syc1(:,:)
        vicmom(:,:,23) = sxxc1(:,:)
        vicmom(:,:,24) = syyc1(:,:)
        vicmom(:,:,25) = sxyc1(:,:)
        vicmom(:,:,26) = sxc2(:,:)
        vicmom(:,:,27) = syc2(:,:)
        vicmom(:,:,28) = sxxc2(:,:)
        vicmom(:,:,29) = syyc2(:,:)
        vicmom(:,:,30) = sxyc2(:,:)
        vicmom(:,:,31) = sxst(:,:)
        vicmom(:,:,32) = syst(:,:)
        vicmom(:,:,33) = sxxst(:,:)
        vicmom(:,:,34) = syyst(:,:)
        vicmom(:,:,35) = sxyst(:,:)

        return

      end subroutine copy_to_vicmon

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      end module moment_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
