c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
#include "choixcomposantes.h"
!dmr -- Added optional components choice - Tue Jun 14 15:59:05 CEST 2011
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      module fait_pointer

      contains

      SUBROUTINE fait_pointer_CC_OCN(KOD)

!! START_OF_USE_SECTION

       use global_constants_mod, only: dblp=>dp, ip, tzero=>tK_zero_C      

#if ( PATH >= 1 )
       use path_mod, only: sea_mask, epais, k_fond, k_surf
#endif

      USE declars_mod, only: LT, NOC_CBR, JT
      use para0_mod, only: kmax
      use dynami_mod, only: area
      use bloc0_mod, only: dz, tms, scal, dz, zw, z
      use bloc_mod, only: normUV_ERA5,zmix_cc,LFe_PISCES 
      use ice_mod, only: albq, tairox, tairoy, fsolcn

#if ( BATHY >= 1 )
      use update_clio_bathy_tools, only: tms_prev, diff_tms
      use C_res_mod, only: alk_oc_rest
#endif

      use loveclim_transfer_mod, only: MGT, DVOL, TM, SM, TM_surface,
     &        ZZ, ZX, mid_level, FRICE, WS_OC, SABST_O, SQRO2, ZMIX_OCN
     &      , OVOL, WIND_ERA5, IRON_LIM, total_area

      use mbiota_mod, only: oc_bottom_cell

!! END_OF_USE_SECTION

       IMPLICIT NONE

       integer(kind=ip) :: k, n, i
       integer(kind=ip), intent(in) :: KOD

!dmr !!!
!dmr --- Commentaire
!dmr --- Les tableaux CLIO et CLIMBER ne sont pas ecrits dans le meme ordre.
!dmr       Pour CLIO : VarOcean(imax,jmax,kmax)   => Lon, Lat, Levels
!dmr       Pour CLIMBER : VarOcean(LT,JT,NOC_CBR) => Lat, Levels, Lon
!dmr --- Il s'ensuit donc naturellement que pour passer de l'un a l'autre a niveau fixe : 
!dmr       VarCLIMBER (:,level=level0,:) = TRANSPOSE(VarCLIO(:,:,level=level0)
!dmr !!!


       do n=1,NOC_CBR
         do k=1,kmax
          do i=1,LT

         DVOL(i,(kmax+1-k),n) = area(n+1,i)*tms(n+1,i,k)*dz(k)
         MGT(i,(kmax+1-k),n) = tms(n+1,i,k)
#if ( BATHY >= 1)
         MGT_prev(i,(kmax+1-k),n) = tms_prev(n+1,i,k)
         DVOL_prev(i,(kmax+1-k),n) = area(n+1,i)*tms_prev(n+1,i,k)*dz(k)
         diff_MGT(i,(kmax+1-k),n) = diff_tms(n+1,i,k)
#endif

         TM(i,(kmax+1-k),n) = scal(n+1,i,k,1) - tzero

         SM(i,(kmax+1-k),n) = scal(n+1,i,k,2)
         zz(kmax+1-k) = dz(k) ! defined in defgrid, grid cell thickness
         IRON_LIM(i,(kmax+1-k),n) = LFe_PISCES(n+1,i,k)

          enddo
         enddo
        enddo
        
       do i=1,LT
         do n=1,NOC_CBR
             TM_surface(i,n) = TM(i,1,n)
         enddo
       enddo


       do k=1,kmax+1
         ZX(kmax+2-k) = zw(k)*(-1.0)
!#if ( CORAL == 1 )
         mid_level(kmax+2-k) = z(k) ! negative values
!#endif
       enddo

       do n=1,NOC_CBR
        do i=1,LT
         FRICE(i,n) = 1.0-albq(n+1,i)
         WS_OC(i,n) = sqrt(tairox(n+1,i)**2+tairoy(n+1,i)**2)
         SABST_O(i,n) = fsolcn(n+1,i)
         SQRO2(i,n) = area(n+1,i)
         ZMIX_OCN(i,n) = zmix_CC(n+1,i)
         WIND_ERA5(i,n,:) = normUV_ERA5(n+1,i,:)
        enddo
       enddo

       !initialisation at the beginning of ocean volume and ocean area
       if (KOD.eq.0) then
          OVOL = SUM(DVOL, MASK=(MGT.EQ.1))
          write(*,*) 'OVOL=', OVOL
          total_area = SUM(SQRO2, MASK=(MGT(:,1,:).EQ.1))
          write(*,*) 'total_area=', total_area
#if ( BATHY >=1 )
          OVOL_prev = SUM(DVOL_prev, MASK=(MGT_prev.EQ.1))
          write(*,*) 'OVOL_prev=', OVOL_prev
!         OALK_ini=alk_oc_rest*OVOL_prev/OVOL
!        write(*,*) 'OALK_ini dans fait_pointer_CC_OCN ', OALK_ini
#endif


!dmr --- [NOTA] "Is this cell a bottom cell? Defined once here for good.

        do n=1,UBOUND(oc_bottom_cell,DIM=3)
        do k=1,UBOUND(oc_bottom_cell,DIM=2)
        do i=1,UBOUND(oc_bottom_cell,DIM=1)
        if ( k.eq.JT ) then
          oc_bottom_cell(i,k,n) = .true.
        else
          if (MGT(i,k+1,n).ne.1) then
            oc_bottom_cell(i,k,n) = .true.
          else
            oc_bottom_cell(i,k,n) = .false.
          endif
        endif

        enddo
        enddo
        enddo


       endif

#if ( PATH >= 1 )

       sea_mask => tms
       epais => dz
       k_surf(:,:) = kmax
       k_fond(:,:) = kfs(:,:)
#endif

       END SUBROUTINE fait_pointer_CC_OCN

      end module fait_pointer
