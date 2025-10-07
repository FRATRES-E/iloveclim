!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      Ce module est un portage en FORTRAN 90 du module CLIO initial
!       cree en FORTRAN 77.
!      (dans l'environnement logiciel LUDUS)
!
!      Auteur : ??, Didier M. Roche
!      Date   : ??
!      Derniere modification : 19 Aout 2014
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module thermo_mod

      use para0_mod, only: nbpt
      use para_mod,  only: nkb0

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! originally written with IMPLICIT INTEGER (I-N), REAL (A-H, O-Z)
! transferred as such in fortran90
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Thermo.com is incorparated by an instruction include in thersf.f ,
!      fontbc.f and acrlbq.f. It comprises the commons associated to
!      thermodynamic ice computation
!
! Correspondance between the variables
! qlbqb   qlbq
! qcmbqb  qcmbq
! thcmb   thcm
! fstbqb  fstrbq
! fltbqb  ffltbq
! fscbqb  fscmbq
! fsolgb  fsolg
! ratbqb  ratbqg
! psbqb   psbq
! tabqb   tabq
! qabqb   qabq
! vabqb   vabq
! qfvbqb  qfvbq
! tsb     ts
! tfub    tfu
! hnpbqb  zhnpbq
! hnbqb   hnbq
! hgbqb   hgbq
! albqb   albq
! qstbqb  qstobq
! fbbqb   fbbq
! tbqb    tbq
! dmgbqb  dmgbq
! dmnbqb  dmnbq
! qlbbqb  qlbsbq
! cldqb   cloud
! dmgwib  dmgwi
! npb     number of points where computations has to be done
! npac    correspondance between the points
! fratsb  firg
! fcsb    fcsg
! fleb    fleg
! dvsbqb  dvosbq
! dvbbqb  dvobbq
! dvlbqb  dvolbq
! dvnbqb  dvonbq
! hgcolb  hgcol (applies only when Cvhg is activated)
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "combq"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      real(kind=8), dimension(nbpt) :: qlbqb, qcmbqb, thcmb, fstbqb, fltbqb, fscbqb, fsolgb, ratbqb, psbqb, tabqb, qabqb, vabqb    &
                  , qfvbqb, tsb, tfub, hnpbqb, hnbqb, hgbqb, albqb, qstbqb, fbbqb, dmgbqb, dmnbqb, qlbbqb, cldqb, dmgwib, hgcolb
      integer,      dimension(nbpt) :: npb, npac

      real(kind=8), dimension(nbpt,nkb0):: tbqb

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr Bloc "comdbq"
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8), dimension(nbpt):: fratsb, fcsb, fleb, dvsbqb, dvbbqb&
                  , dvlbqb, dvnbqb


      end module thermo_mod
!
!cp2  common/fdericeb/fderb(nbpt)
!-- End of file 'thermo.com'
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr --- The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
