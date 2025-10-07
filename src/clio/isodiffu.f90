!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:48 CET 2009

      MODULE isodiffu_mod

!! START_OF_USE_SECTION
      
      use para0_mod, only: kmax
      
!! END_OF_USE_SECTION

      IMPLICIT NONE

      real, dimension(kmax), save, private :: dtsdx, dtsdy, dtsdz !--variables locales conservees d'un appel a l'autre :

      contains

      SUBROUTINE isodiffu_init()

!! START_OF_USE_SECTION

      use bloc_mod, only: numit, nstart
      use bloc0_mod, only: ahs, ai, unsdz, unsdx, unsdy, dts

      use newunit_clio_mod, only: clio3_out_id
      
!! END_OF_USE_SECTION

      IMPLICIT NONE

      !dmr --- local loop variable
      integer :: k

!dmr [UPDATED]      if (numit.eq.nstart .and. ns.eq.1) then
      if (numit.eq.nstart) then
      
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ATTENTION : Modif de AHS !!!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        write(clio3_out_id,'(A,3(A,F6.0))') &
          ' *** isodiffu ; Isopycn. *** :', &
          ' ahs=', ahs(kmax),' +', ai(kmax),'(=ai)'
        
        do k=1,kmax
          ahs(k) = ahs(k)+ai(k)
          dtsdz(k) = unsdz(k)*dts(k)
          dtsdx(k) = unsdx*dts(k)
          dtsdy(k) = unsdy*dts(k)
        enddo
 
!- fin du traitement 1ere iter.
      endif
    
      END SUBROUTINE isodiffu_init


      SUBROUTINE isodiffu(scalat,scal_kopie)
!~       SUBROUTINE isodiffu(scalat,scal_kopie, ns)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!=============================================================
!  Adaptation of Cox small-slope Isopycnal Tensor
!  References:
!  -Cox, M.D, 1987, Isopycnal diffusion in a z-coordinate
!        model, Ocean Modelling, 74, 1-5
!  -Redi, M.H., 1982, Oceanic isopycnal mixing by coordinate
!        rotation, JPO, 12, 1154-1158
!=============================================================
!  modif : 25/05/99
!  last updated version: 2021/12/02 -- dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use const_mod, only: epsil
      use para0_mod, only: imax, jmax, kmax

!dmr [NOTUSED]      use para_mod,  only:
      
      use bloc0_mod, only: ai, unsdx, unsdy, js1, js2, ks1, ks2, is1, is2, tms, unsdzw, ims1, ims2, iberpm, jberp, ibera &
                   , jberam, iberp, iberam, jbera, jberpm
      use bloc_mod, only: cmy, cmx, smy, smx, smxy
      use isoslope_mod, only: ttm1, c1x, c1x, c4x, jsdom1, jsdom2, ttm2, c4y, c2y

!dmr [DELETED]      use isoslope_mod, only: dscalz    ! local variable now
!dmr [DELETED]      use isoslope_mod, only: fisoz     ! local variable now
!dmr [NOTUSED]      use comunit, only: iuo


      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!- by reference variables :
!dmr --- scalat is scal(i,j,k,ns) for the tracer considered at the partial
!         step we are in, see comment in scale.f
!dmr ---        in scale.f scalat(i,j,k) = scal(i,j,k,ns)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!~       integer                        , intent(in)    :: ns

      real, dimension(imax,jmax,kmax), intent(inout) :: scalat
      real, dimension(imax,jmax,kmax), intent(in) :: scal_kopie

!- local variables :
      
! dmr      real, dimension(imax,kmax)      :: dscalx
      real, dimension(imax,jmax,kmax)   :: dscalx_new = 0.0d0
      real, dimension(imax,jmax,kmax)   :: dscalz = 0.0d0
      
      real, dimension(jmax,kmax)        :: dscaly = 0.0d0, fisoy = 0.0d0
      
      real, dimension(imax,jmax,kmax)   :: fisox = 0.0d0
      real, dimension(imax,jmax,kmax+1) :: fisoz = 0.0d0

      !$OMP THREADPRIVATE(dscalx_new,dscalz,dscaly,fisoy,fisox,fisoz)

!- loop control variables
      integer :: i,j,k
      integer :: km0, km1, kp0, kp1

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ****** GRADIENTS OF SCALAR scal *******
!- NB1: for Gradient or Fluxes: in X => isf2(j)+1 & in Y => js2+1
!- NB2: isf1 is prefered to is1 to take into account
!        the Northern Margin (including Bering).
!- NB3: We take Minus flux * metric coefficient !!
!- NB4: BERING Raccord cmx&cmy coherent ONLY at point 2
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  2 ) Merid. Section : Compute Isopycn.Flx in X & Z(part.1) Directions
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!- Vertical Gradient of scalar "ns" :

!dmr switched k and j to be in FORTRAN natural ordering

      do k=ks1+1,ks2
       do j=js1,js2
        do i=is1(j),is2(j)+1
          dscalz(i,j,k) = tms(i,j,k-1) * unsdzw(k) * ( scal_kopie(i,j,k) - scal_kopie(i,j,k-1) )
        enddo
       enddo ! on k, ks1+1 -> ks2       
      enddo  ! on j, js1   -> js2


!- Cyclic Conditions :

#if ( L_TEST >= 1 )
!dmr switched k and j to be in FORTRAN natural ordering


       do k=ks1+1,ks2
         do j =js1,js2

        dscalz(1,j,k) = dscalz(ims2,j,k)
        dscalz(imax,j,k) = dscalz(ims1,j,k)
        
         enddo ! on j
       enddo ! on k


#endif        


!dmr switched k and j to be in FORTRAN natural ordering



       do k=ks1,ks2
       
        km1 = max(k-1, ks1)
        km0 = km1 + 1
        kp1 = min(k+1, ks2)
        kp0 = kp1 - 1
        
       do j =js1,js2        
        do i=is1(j),is2(j)+1
!- Grad.X of scalar "ns" :
          dscalx_new(i,j,k) = ttm1(i,j,k) * unsdx * smx(i,j,1) * ( scal_kopie(i,j,k) - scal_kopie(i-1,j,k) )
!- X.Flx of scalar "ns" :
          fisox(i,j,k) = cmy(i,j,1) * ai(k) * c1x(i,j,k) * ( dscalz(i,j,kp1) + dscalz(i-1,j,kp1)                                &            
            + dscalz(i,j,km0) + dscalz(i-1,j,km0) ) / ( tms(i,j,kp0) + tms(i-1,j,kp0) + tms(i,j,km1) + tms(i-1,j,km1) + epsil )
        enddo ! on is1 -> is2+1
        enddo ! on j       
        
       enddo  ! on ks1 -> ks2


!dmr switched k and j to be in FORTRAN natural ordering
!dmr [NOTA] Following loop is inefficient, since calling k and k-1

      do k=ks1+1,ks2
        do j =js1,js2        
!- Z.Flx (1rst part) of scalar "ns" :
         do i=is1(j),is2(j)
          fisoz(i,j,k) = c4x(i,j,k) * ( dscalx_new(i,j,k) + dscalx_new(i+1,j,k-1) + dscalx_new(i,j,k-1) + dscalx_new(i+1,j,k) ) &
                                    / ( ttm1(i,j,k) + ttm1(i+1,j,k-1) + ttm1(i,j,k-1) + ttm1(i+1,j,k) + epsil )
         enddo ! on is1   -> is2
        enddo ! on j
       enddo  ! on ks1+1 -> ks2

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!- End of external loop on all Meridional Section, index "j".



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  3 ) Zonal Section : Compute Isopycn.Flx in Y & Z(part.2) Directions
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!- Raccord a Bering :

#if ( L_TEST == 3 )

        do k=ks1+1,ks2
          dscalz(iberpm,jberp,k) = dscalz(ibera, jberam,k)
          dscalz(iberp, jberp,k) = dscalz(iberam,jberam,k)
          dscalz(iberam,jbera,k) = dscalz(iberp, jberpm,k)
          dscalz(ibera, jbera,k) = dscalz(iberpm,jberpm,k)
        enddo

#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!- Start  external loop on all Zonal Section, index "i".

      do i=ims1,ims2
       do k=ks1,ks2
        km1 = max(k-1, ks1)
        km0 = km1 + 1
        kp1 = min(k+1, ks2)
        kp0 = kp1 - 1
        do j=jsdom1(i),jsdom2(i)+1
!- Grad.Y of scalar "ns" :
          dscaly(j,k) = ttm2(i,j,k) * unsdy * smy(i,j,2) * ( scal_kopie(i,j,k) - scal_kopie(i,j-1,k) )
        enddo
       enddo

       do k=ks1+1,ks2
        do j=jsdom1(i),jsdom2(i)
!- Z.Flx (2nd  part) of scalar "ns" :
          fisoz(i,j,k) = ai(k) * ( fisoz(i,j,k) + c4y(i,j,k) * ( dscaly(j,k) + dscaly(j+1,k-1)                                 &
              + dscaly(j,k-1) + dscaly(j+1,k) ) / ( ttm2(i,j,k) + ttm2(i,j+1,k-1) + ttm2(i,j,k-1) + ttm2(i,j+1,k) + epsil ) )
        enddo 
       enddo
      enddo ! on i

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  4 ) Bilan des Flux Incorpore dans "scalat".
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!- BALANCE OF FLUXES
!- !! We take Minus the gradient of fluxes we took minus the flux
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      do i=ims1,ims2
       do k=ks1,ks2
         km1 = max(k-1, ks1)
         km0 = km1 + 1
         kp1 = min(k+1, ks2)
         kp0 = kp1 - 1
         do j=jsdom1(i),jsdom2(i)+1
!- Y.Flx of scalar "ns" :
            fisoy(j,k) = cmx(i,j,2) * ai(k) * c2y(i,j,k) * ( dscalz(i,j,kp1) + dscalz(i,j-1,kp1)                                 &
             + dscalz(i,j,km0) + dscalz(i,j-1,km0) ) / ( tms(i,j,kp0) + tms(i,j-1,kp0) + tms(i,j,km1) + tms(i,j-1,km1) + epsil )
          enddo
        enddo

      do k=ks1,ks2
       do j=jsdom1(i),jsdom2(i)
        scalat(i,j,k) = scalat(i,j,k) + smxy(i,j,0)*(dtsdx(k)*(fisox(i,j,k)-fisox(i+1,j,k))+dtsdy(k)*(fisoy(j,k)-fisoy(j+1,k)))  &
                      + dtsdz(k) * (fisoz(i,j,k)-fisoz(i,j,k+1))
        enddo
       enddo 
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!- End of external loop on all Zonal Section, index "i".
      enddo
      
      return
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!- fin de la routine isodiffu -
      end subroutine isodiffu


      END MODULE isodiffu_mod
