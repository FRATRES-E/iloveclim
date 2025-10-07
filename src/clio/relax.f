!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:52 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:52 CET 2009

      SUBROUTINE relax(ttdif,residu)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  calcul des temperatures de rappel/diffusion de S.Rahmstorf
!--Convention : fflx, fss : + vers le haut, si refroidit l'ocean.
! ahrap <-> Diffusivite (Unite : m^2/s) <-> mu/(Dz*Rho*Cp) (<- S.Rahmstorf)
!  ahu  <-> 1 / tps_rappel (Unite : /s) <-> gamma/(Dz*Rho*Cp) (<- S.Rahmstorf)
!            et tps_rappel = 1_an / Fraction_d'annee(=rapp0(1)*unstyr)
!  ahe  <-> Coef de sous/sur_relaxation (sans unite).
!  pseudo flux (fss/dz, fflx) (Unite : flux/(Dz.rho.Cp) <-> deg_C.s-1 )
!-----
!  modif : 02/09/97

!! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip           

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      
      use newunit_clio_mod, only: clio3_out_id
      
!! END_OF_USE_SECTION

!--dummy variables :
      real(kind=dblp), dimension(imax,jmax) :: ttdif, residu


!--variables conservees d'un appel a l'autre :
      real(kind=dblp), dimension(imax,jmax), save :: cphix, cphiy
     &                                             , unsdia

!--variables locales :
      real(kind=dblp), dimension(imax,jmax) :: phix, phiy
      
      real(kind=dblp) :: ccdif, ddiag, ddmx, ddsc, epsil2, fflx, rsdm
     &                 , rsdmn, rsdmx, ttm1ob, ttm2ob
      
      integer(kind=ip) :: i, j, k

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  1 ) Initialisation et Calculs preliminaires :
!-----------------------------------------------------------------------

      if (numit.eq.nstart) then
!--done only at the 1rst iteration.

      do j=1,jmax
       do i=1,imax
         cphix(i,j)  = 0.
         cphiy(i,j)  = 0.
         unsdia(i,j) = 0.
       enddo
      enddo

      k = ks2
      ddsc = 0.
      do j=js1,js2
       do i=is1(j),is2(j)
         ddsc = ddsc + ctmi(i,j,k,0)*abs(ttdif(i,j))
       enddo
      enddo
      ddsc = ddsc * zsurf
      write(clio3_out_id,2222) ' numit, ddsc(debut) :', numit, ddsc

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!--Calcul des coeffs de flux :
      ccdif = ahrap * unsdx
      do j=js1,js2
       do i=is1(j),1+is2(j)
        if ( tms(i,j,k).eq.one .and. scalr(i,j,k,1).eq.spvr )
     &    write(clio3_out_id,'(A,2I4,1P2E18.10)') 'X: i, j,tms,scalr :',
     &       i,  j, tms(i,j,k), scalr(i,j,k,1)
        if (tms(i-1,j,k).eq.one.and.scalr(i-1,j,k,1).eq.spvr)
     &    write(clio3_out_id,'(A,2I4,1P2E18.10)') 'X:i-1,j,tms,scalr :',
     &       i-1, j, tms(i-1,j,k), scalr(i-1,j,k,1)
!---
        ttm1ob = min(tms(i-1,j,k), (scalr(i-1,j,k,1)-spvr),
     &               tms(i,j,k),   (scalr(i,j,k,1)-spvr)   )
!       cphix(i,j) = ccdif * smxy(i,j,1) * ttm1ob
        cphix(i,j) = ccdif * smxy(i,j,1) * tms(i-1,j,k) * tms(i,j,k)
       enddo
      enddo
      ccdif = ahrap * unsdy
      do j=js1,1+js2
       do i=isf1(j),isf2(j)
        if ( tms(i,j,k).eq.one .and. scalr(i,j,k,1).eq.spvr )
     &    write(clio3_out_id,'(A,2I4,1P2E18.10)') 'Y:i, j, tms,scalr :',
     &       i,  j, tms(i,j,k), scalr(i,j,k,1)
        if (tms(i,j-1,k).eq.one.and.scalr(i,j-1,k,1).eq.spvr)
     &    write(clio3_out_id,'(A,2I4,1P2E18.10)') 'Y:i,j-1,tms,scalr :',
     &       i, j-1, tms(i,j-1,k), scalr(i,j-1,k,1)
!---
        ttm2ob = min(tms(i,j-1,k), (scalr(i,j-1,k,1)-spvr),
     &               tms(i,j,k),   (scalr(i,j,k,1)-spvr)   )
!       cphiy(i,j) = ccdif * cmxy(i,j,2) * ttm2ob
        cphiy(i,j) = ccdif * cmxy(i,j,2) * tms(i,j-1,k) * tms(i,j,k)
       enddo
      enddo

!--Calcul des coeffs diagonaux :
      epsil2 = epsil * unsdx * unsdy
      ddmx = cstmin
      do j=js1,js2
       do i=is1(j),is2(j)
         ddiag = ahu + smxy(i,j,0) *
     &     ( unsdx * (cphix(i,j) + cphix(i+1,j))
     &     + unsdy * (cphiy(i,j) + cphiy(i,j+1)) )
         unsdia(i,j) = tms(i,j,k) / max(epsil2, ddiag)
         ddmx = max(unsdia(i,j), ddmx)
       enddo
      enddo
      write(clio3_out_id,*) 'Max. 1/Coef_Diagonaux =', ddmx

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- end of the "only 1rst iteration" executed part.
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  2 ) calcul de l'elevation, cas avec Filtre (ahe non nul) :          |
!-----------------------------------------------------------------------

!- Calcul de phiy & phix :
      do j=js1,js2
       do i=is1(j),1+is2(j)
        phix(i,j) = cphix(i,j) * ( ttdif(i-1,j) - ttdif(i,j) )
       enddo
      enddo
      do j=js1,1+js2
       do i=isf1(j),isf2(j)
        phiy(i,j) = cphiy(i,j) * ( ttdif(i,j-1) - ttdif(i,j) )
       enddo
      enddo

!--Calcul du residu (unite : degre_C) :
      do j=js1,js2
       do i=is1(j),is2(j)
         fflx = -ahu * ttdif(i,j) + smxy(i,j,0) *
     &        ( unsdx * (phix(i,j) - phix(i+1,j))
     &        + unsdy * (phiy(i,j) - phiy(i,j+1)) )
         residu(i,j) = ( fflx - fss(i,j,1)*unsdz(ks2) ) * unsdia(i,j)
         ttdif(i,j) = ttdif(i,j) + ahe * residu(i,j)
       enddo
      enddo

!-- evaluation globale du residu :
      if ( numit.eq.nstart .or. mod(numit,ninfo).eq.0) then
        ddsc = 0.
        rsdm = 0.
        rsdmn = cstmax
        rsdmx = cstmin
        do j=js1,js2
         do i=is1(j),is2(j)
           ddsc = ddsc + ctmi(i,j,ks2,0)*abs(ttdif(i,j))
           rsdm = rsdm + ctmi(i,j,ks2,0)*abs(residu(i,j))
           rsdmn = min(residu(i,j), rsdmn)
           rsdmx = max(residu(i,j), rsdmx)
         enddo
        enddo
        ddsc = ddsc * zsurf
        rsdm = rsdm * zsurf
        write(clio3_out_id,2222) ' it, ddsc, rsdm,n,x :',
     &               numit, ddsc, rsdm, rsdmn, rsdmx
      endif

!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!  3 ) raccords cycliques et autres .                                  |
!-----------------------------------------------------------------------

      call raccord(ttdif, zero, 1, 8)

      return
 2222 format(A,I6,1P4E14.6)
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!- fin de la routine relax -
      end
