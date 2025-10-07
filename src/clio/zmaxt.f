!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:58 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:58 CET 2009

!
!     ***************************
      FUNCTION zmaxt(ndeb,nfin,a)
!     ***************************
!
!                CALCULATE MAXIMUM VALUE OF THE COMPONENTS OF A VECTOR.

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

!! END_OF_USE_SECTION

!
      real(kind=dblp) :: zmaxt
      real(kind=dblp), dimension(*) :: a
      integer(kind=ip) :: ji, ndeb, nfin

!
      zmaxt = a(ndeb)
      do 10 ji=ndeb,nfin
        zmaxt = max(zmaxt,a(ji))
10    continue
!
      return
      end
