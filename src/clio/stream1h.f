!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:54 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:54 CET 2009

      SUBROUTINE stream1h(psid,ishsf,iehsf,jshsf,jehsf,ndhsf)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
! Pour chaque detroit nvd (nvd=1,ndhsf) :
!  calcule le transport Horizontal qui traverse la section
!  meridienne de longitude ishsf, entre les latitudes  jshsf et jehsf
!   ou zonale de latitude  jshsf, entre les longitudes ishsf et iehsf
! en separant Flux > 0 et Flux < 0.
!  modif : 01/06/95

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
!! END_OF_USE_SECTION

!--dummy argument :
      real(kind=dblp), dimension(2,nhsfmx):: psid
      integer(kind=ip), dimension(nhsfmx) :: ishsf, iehsf, jshsf, jehsf


!---- more locales
      integer(kind=ip):: i, ii, j, jj, k, ndhsf, nv
      real(kind=dblp) :: u1dyz, v2dxz


!--facteur de conversion (en Sverdrup) (= svrdrp) .

      do nv=1,ndhsf
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Boucle sur l'ensembles des detroits

      psid(1,nv) = 0.0
      psid(2,nv) = 0.0

      if (iehsf(nv).eq.0) then
!--integration de (-u).(-dy) :
        ii = ishsf(nv)
        do j=jshsf(nv),jehsf(nv)
         do k=ks1,ks2
          u1dyz = tms(ii,j,k) * tms(ii-1,j,k) * cmy(ii,j,1) *
     &            dz(k) * ( u(ii,j,k) + u(ii,j+1,k) )
          psid(1,nv) = psid(1,nv) + min(zero, u1dyz)
          psid(2,nv) = psid(2,nv) + max(zero, u1dyz)
         enddo
        enddo
        psid(1,nv) = 0.5 * dy * psid(1,nv) * svrdrp
        psid(2,nv) = 0.5 * dy * psid(2,nv) * svrdrp
      endif

      if (jehsf(nv).eq.0) then
!--integration de vb.dx :
        jj = jshsf(nv)
        do i=ishsf(nv),iehsf(nv)
         do k=ks1,ks2
          v2dxz = tms(i,jj,k) * tms(i,jj-1,k) * cmx(i,jj,2) *
     &            dz(k) * ( v(i,jj,k) + v(i+1,jj,k) )
          psid(1,nv) = psid(1,nv) + min(zero, v2dxz)
          psid(2,nv) = psid(2,nv) + max(zero, v2dxz)
         enddo
        enddo
        psid(1,nv) = 0.5 * dx * psid(1,nv) * svrdrp
        psid(2,nv) = 0.5 * dx * psid(2,nv) * svrdrp
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      enddo

      return
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine stream1h -
      end
