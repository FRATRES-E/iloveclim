!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:53 CET 2009

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      module ec_convec_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      

      use global_constants_mod, only: dblp=>dp, ip
      use newunit_mod, only: error_id

      implicit none


      contains
      
#if ( ISOATM == 0 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! ***   Moist convective adjustment
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      SUBROUTINE ec_convec_two

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Imports ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      use comatm, only: dtime, grav, rgas, plevel, nlat, nlon, iwater
      use comdyn, only: geopg
      use comphys,only: cpair, dp1, dp2, gamd, ivavm, potfac1, potfac2, rainmax, relhmax, rkappa, rlatsub, rlatvap, rowat, tzero  &
                 , rmoisg, temp2g, temp4g, temp2g, rmoisg, rmoisg, rmoisg, rmoisg, dyrain, dyrain, dysnow, dysnow, corain, corain &
                 , cosnow, cosnow, rmoisg, temp4g, corain, cosnow, thforg2, temp2g, temp4g, vhforg1, vhforg2, rmoisg, temp4g      &
                 , tetacr, tetacr, temp2g, temp4g, dyrain, corain, dysnow, cosnow, temp2g, temp4g, gams, teta, torain, tosnow
      use comrunlabel_mod, only: irunlabelf
      use comsurf_mod
      use comunit



#if ( DOWNSTS == 1 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr added for the precipitation downscaling
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use ecbilt_topography,       only: nb_levls, qmount_virt
      use atmphys_d,               only: ec_ptmoisgp_d, ec_ptmoisgp_dreduce
      use vertDownsc_mod,          only: tsurf_d

#if ( DOWNSCALING == 2 )
! dmr&afq : sub_grid_not_flat == mask : is 0.0 when sub-grid is flat and nothing is to be done

      use input_subgrid2L,         only: sub_grid_notflat, max_nb_points, relhum_sg, dyrain_sg, dysnow_sg, corain_sg, cosnow_sg   &
                         , torain_sg, tosnow_sg, tsurf_sg, nbpointssg, weights_low_sg, index_low_sg, area_sg, area_sg_onecb       &
                         , nneigh, index_interpL2G, weights_interpL2G, sumweights_interpL2G, cluster_var_subgrid_in_ecbilt_dble
      use taillesGrilles,          only: sgnx,sgny,sgnxm,sgnym,sgd
      use vertDownsc_mod,          only: tsurf_d_interp_sg
      use interpolate_mod,         only: interpolate
#endif /* DOWNSCALING == 2 */

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


#endif /* DOWNSTS == 1 */

#if ( DOWNSTS == 1 )
      use vertDownsc_mod, only: dyrain_d, dysnow_d, corain_d, cosnow_d, torain_d, tosnow_d
#endif

!!    USE OMP_LIB


      implicit none

      integer ncmax,iconvn,i,j

      double precision  :: dcmoisg


      real*8  qsatcr,tsatcr,pref,t500,qsat500,pot2g,pot4g
      real*8  fachulp,facteta,factemv,factems,pmount,tmount
      real*8  qmax,ec_qsat,hulp,redrain
      real*8  temp2go,temp4go,ec_levtempgp
      real*8  ec_detqmax,drainm,crainm,dqmdt

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  New variables for vertical downscaling on nb_ipoints virtual surfaces
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( DOWNSTS == 1 )

      logical                                 :: success = .false.
      double precision, dimension(nb_levls)   :: tmount_d,qmax_d,dqmdt_d,dcmoisg_d

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Series of variables for looping over the pertaining sub-grid points
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      double precision                          :: weight_low
      integer                                   :: ind_low, I_did_something

#if ( DOWNSCALING == 2 )
      integer :: iconvn_grisli

! dmr [NOTA] Have to keep the (dy,co)(rain,snow)_grisli variables on all points since there is a data SUM array in the end
! dmr        Same story for temp2g_grisli, temp4g_grisli, thforg2_grisli, vhforg1_grisli, vhforg2_grisli

! dmr The cluster below hence represents all variables providing a feedback from small to large scale
      double precision, dimension(max_nb_points):: dyrain_grisli , dysnow_grisli , corain_grisli , cosnow_grisli, thforg2_grisli  &
                     , vhforg1_grisli, vhforg2_grisli, temp2g_grisli , temp4g_grisli, rmoisg_grisli
! dmr Those two below are used only once as a transfer matrix. Not much I can do for them
      double precision, dimension(max_nb_points):: temp2g_grisli_diag, temp4g_grisli_diag

      double precision :: qmax_hr, dqmdt_hr, tsurf_hr, dcmoisg_hr, rmoisg_hr

      double precision :: dcmoisg_accum
      double precision :: drainm_grisli, crainm_grisli
      double precision :: teta_grisli, tetacr_grisli,gams_grisli
      double precision, parameter :: fact_humid = 1.0

      double precision,dimension(sgnxm,sgnym,sgd)   :: rmoisg_interp,temp2g_interp,temp4g_interp,tsurf_d_interp

      double precision,dimension(nlat,nlon,max_nb_points) :: rmoisg_interp_sg,temp2g_interp_sg,temp4g_interp_sg

!     ,tsurf_interp_sg_tmp  !! now in vertDownsc_mod
!     double precision, dimension(nlat,nlon,nb_levls,max_nb_points) ::
!     &  tsurf_d_interp_sg  !! now in vertDownsc_mod
      logical :: stability_coarse

#endif

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Added for performance monitoring
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

      integer(kind=8) :: count_strt, count_stop

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

      integer :: n_point, lokaal_npoints, nbd

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
      double precision :: rain_cg, snow_cg, rain_part

#endif /* DOWNSTS */

#define SYS_CLOCK 0

#if ( DOWNSTS == 1 )
       integer :: nl
#endif /* DOWNSTS */

#if ( SYS_CLOCK > 0 )

      call system_clock(count=count_strt)

#endif

! dmr [TODO] those computations could be moved to a common, no need to do this at every
!      call !
! This is to be checked in the STD SVN VERSION ...
      fachulp=0.622d0*(rlatvap**2)/(cpair*rgas)
      facteta=0.6d0*rgas*(2.d0**rkappa)/grav
      factemv=rlatvap*grav*rowat/cpair
      factems=rlatsub*grav*rowat/cpair
      ncmax=0
      crainm=0.5*rainmax
      drainm=rainmax-crainm
! [TODO]

#if ( DOWNSCALING == 2 )
! afq, it is possible to rain more/less in subgrid:
      crainm_grisli=crainm*1.d0
      drainm_grisli=drainm*1.d0
#endif


!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

#if (DOWNSCALING == 2 )

      do nbd=1,sgd
         
!~ !$OMP PARALLEL SECTIONS

!~ !$OMP SECTION
         success = interpolate(index_interpL2G(:,:,1:sgnx(nbd),1:sgny(nbd),nbd),            &
              weights_interpL2G(:,1:sgnx(nbd),1:sgny(nbd),nbd),                           &
              sumweights_interpL2G(1:sgnx(nbd),1:sgny(nbd),nbd),                        &
              transpose(rmoisg(:,:,iwater)),rmoisg_interp(1:sgnx(nbd),1:sgny(nbd),nbd),     &
              sgnx(nbd), sgny(nbd), 3, nneigh, nlon, nlat)

!~ !$OMP SECTION
         success = interpolate(index_interpL2G(:,:,1:sgnx(nbd),1:sgny(nbd),nbd),            &
              weights_interpL2G(:,1:sgnx(nbd),1:sgny(nbd),nbd),                           &
              sumweights_interpL2G(1:sgnx(nbd),1:sgny(nbd),nbd),                        &
              transpose(temp2g),temp2g_interp(1:sgnx(nbd),1:sgny(nbd),nbd),                 &
              sgnx(nbd), sgny(nbd), 3, nneigh, nlon, nlat)

!~ !$OMP SECTION
         success = interpolate(index_interpL2G(:,:,1:sgnx(nbd),1:sgny(nbd),nbd),            &
              weights_interpL2G(:,1:sgnx(nbd),1:sgny(nbd),nbd),                           &
              sumweights_interpL2G(1:sgnx(nbd),1:sgny(nbd),nbd),                        &
              transpose(temp4g),temp4g_interp(1:sgnx(nbd),1:sgny(nbd),nbd),                 &
              sgnx(nbd), sgny(nbd), 3, nneigh, nlon, nlat)

!~ !$OMP END PARALLEL SECTIONS

      enddo ! on the number of subgrids

      call cluster_var_subgrid_in_ecbilt_dble(rmoisg_interp, rmoisg_interp_sg)
      call cluster_var_subgrid_in_ecbilt_dble(temp2g_interp, temp2g_interp_sg)
      call cluster_var_subgrid_in_ecbilt_dble(temp4g_interp, temp4g_interp_sg)

#endif

#if ( SYS_CLOCK > 0 )
      call system_clock(count=count_stop)
      write(*,*) '("Time/1 = ",i6," counts")',count_stop - count_strt
#endif

      do j=1,nlon
        do i=1,nlat

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
! dmr   STD concistency check
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

          iconvn=0 ! global to the grid cell ...
#if ( DOWNSCALING == 2 )
          stability_coarse = .False.
#endif

          do while (iconvn.le.10) !afq -- 10 continue

#if ( DOWNSTS == 1 )

#if ( SYS_CLOCK > 0 )
      call system_clock(count=count_strt)
#endif

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
! dmr  Vertical downscaling counter part of ec_ptmoisg
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

         success = ec_ptmoisgp_d(temp2g(i,j),temp4g(i,j),geopg(i,j,2),tmount_d,qmax_d,dqmdt_d,qmount_virt(:))

         qmax_d(:) = relhmax*qmax_d(:)

         dcmoisg_d(:)=0d0

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
! dmr&afq  Vertical downscaling on vertical levs for precip
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

         do nl=1, nb_levls

            if (rmoisg(i,j,iwater).gt.qmax_d(nl)) then

! redrain is dependent on dqmdt ...
              redrain=1d0+dqmdt_d(nl)*relhmax*rowat*rlatvap*grav/(cpair*dp2)

              dcmoisg_d(nl) = (rmoisg(i,j,iwater)-qmax_d(nl))/(redrain*dtime)

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
! ***         if the air is supersaturated initially, this is due to
! ***         large scale convergence of moisture and large scale
! ***         temperature changes. Excessive moisture
! ***         is then removed as dynamic rain
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

              if (iconvn.eq.0) then

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Temperatures are positive, we are dealing with rain ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

                if (tsurf_d(i,j,nl).ge.tzero + 0.1d0 ) then
                  dyrain_d(i,j,nl)= dyrain_d(i,j,nl)+ dcmoisg_d(nl)

                  if (dyrain_d(i,j,nl).gt.drainm) then
                    dcmoisg_d(nl)=drainm-dyrain_d(i,j,nl)+dcmoisg_d(nl)
                    dyrain_d(i,j,nl)=drainm
                  endif

                else ! tsurf < tzero

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Temperatures are negative, we are dealing with dynamic snow ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

                  dysnow_d(i,j,nl)=dysnow_d(i,j,nl) + dcmoisg_d(nl)

                  if (dysnow_d(i,j,nl).gt.drainm) then
                    dcmoisg_d(nl)=drainm-dysnow_d(i,j,nl) + dcmoisg_d(nl)
                    dysnow_d(i,j,nl)=drainm
                  endif

                endif  ! on tsurf.ge.tzero, dynamic part

                else ! on iconvn == 0
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Temperatures are positive, we are dealing with rain ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

                if (tsurf_d(i,j,nl).ge.tzero + 0.1d0 ) then

                  corain_d(i,j,nl)=corain_d(i,j,nl) + dcmoisg_d(nl)

                  if (corain_d(i,j,nl).gt.crainm) then
                    dcmoisg_d(nl)=crainm-corain_d(i,j,nl) + dcmoisg_d(nl)
                    corain_d(i,j,nl)=crainm
                  endif

                else ! on tsurf > tzero
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Temperatures are negative, we are dealing with snow ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

                  cosnow_d(i,j,nl)=cosnow_d(i,j,nl) + dcmoisg_d(nl)
                  if (cosnow_d(i,j,nl).gt.crainm) then
                    dcmoisg_d(nl)=crainm-cosnow_d(i,j,nl) + dcmoisg_d(nl)
                    cosnow_d(i,j,nl)=crainm
                  endif

                endif ! on tsurf > tzero
              endif ! on iconvn == 0

             endif ! on rmoisg > qmax

           enddo ! on nb_levls

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
! dmr&afq  End vertical downscaling on vertical levs for precip
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

#if ( SYS_CLOCK > 0 )
      call system_clock(count=count_stop)
      write(*,*) '("Time/2 = ",i6," counts")',count_stop-count_strt
#endif

#endif

#if ( DOWNSCALING == 2 )

         if (sub_grid_notflat(i,j) <= 0.0d0 .or. stability_coarse) then ! grid is flat ... normal procedure!

#if ( SYS_CLOCK > 0 )
           call system_clock(count=count_strt)
#endif

#endif /* on DOWNSCALING == 2 */

! ***     calculate pressure and temperature at the ground
! ***     and the maximum water content

            call ec_ptmoisgp(tmount,qmax,i,j,dqmdt)

! ***     relhmax defines the relative humidity at which oversaturation
! ***     occurs

            qmax=relhmax*qmax

            pot2g=temp2g(i,j)/potfac1
            pot4g=temp4g(i,j)/potfac2
            teta(i,j)=0.5d0*(pot2g-pot4g)
! ***       dry adiabatic lapse rate
            tetacr(i,j)=0d0

            dcmoisg=0d0

            if (rmoisg(i,j,iwater).gt.qmax) then

! ***     calculate rain reduction factor to account for increased
! ***     moisture capacity due to latent heat release

              redrain=1d0+dqmdt*relhmax*rowat*rlatvap*grav/(cpair*dp2)

              dcmoisg=(rmoisg(i,j,iwater)-qmax)/(redrain*dtime)

! ***         if the air is supersaturated initially, this is due to
! ***         large scale convergence of moisture and large scale
! ***         temperature changes. Excessive moisture
! ***         is then removed as dynamic rain

              if (iconvn.eq.0) then

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Temperatures are positive, we are dealing with rain ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

                if (tsurf(i,j).ge.tzero + 0.1d0 ) then

                  dyrain(i,j,iwater)=dyrain(i,j,iwater) + dcmoisg
                  if (dyrain(i,j,iwater).gt.drainm) then
                    dcmoisg=drainm-dyrain(i,j,iwater)+dcmoisg
                    dyrain(i,j,iwater)=drainm
                  endif
                  thforg2(i,j)=thforg2(i,j) + factemv*dcmoisg/dp2
                  vhforg2(i,j)=vhforg2(i,j) + factemv*dcmoisg/dp2

                else ! tsurf < tzero
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Temperatures are negative, we are dealing with dynamic snow ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

                  dysnow(i,j,iwater)=dysnow(i,j,iwater) + dcmoisg
                  if (dysnow(i,j,iwater).gt.drainm) then
                    dcmoisg=drainm-dysnow(i,j,iwater)+dcmoisg
                    dysnow(i,j,iwater)=drainm
                  endif
                  thforg2(i,j)=thforg2(i,j) + factems*dcmoisg/dp2
                  vhforg2(i,j)=vhforg2(i,j) + factems*dcmoisg/dp2

                endif  ! on tsurf.ge.tzero, dynamic part
              else ! convective precipitation, iconvn != 0

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Temperatures are positive, we are dealing with rain ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
                if (tsurf(i,j).ge.tzero + 0.1d0 ) then

                  corain(i,j,iwater)=corain(i,j,iwater)+ dcmoisg
                  if (corain(i,j,iwater).gt.crainm) then
                    dcmoisg=crainm-corain(i,j,iwater)+dcmoisg
                    corain(i,j,iwater)=crainm
                  endif
                  temp4g(i,j) =temp4g(i,j) + factemv*dcmoisg*dtime/dp2
                  vhforg2(i,j)=vhforg2(i,j)+ factemv*dcmoisg/dp2

                else

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Temperatures are negative, we are dealing with snow ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

                  cosnow(i,j,iwater)=cosnow(i,j,iwater)+ dcmoisg
                  if (cosnow(i,j,iwater).gt.crainm) then
                    dcmoisg=crainm-cosnow(i,j,iwater)+dcmoisg
                    cosnow(i,j,iwater)=crainm
                  endif
                  temp4g(i,j) =temp4g(i,j) + factems*dcmoisg*dtime/dp2
                  vhforg2(i,j)=vhforg2(i,j)+ factems*dcmoisg/dp2

                endif ! on tsurf > tzero
              endif ! on iconvn == 0

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!   From the above code, corain, cosnow and dcmoisg are determined
!    we are ready to modify the moisture content in the atmosphere
!    (rmoisg)
!
! ***         moisture changes due to dynamic and convective rain
!
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

            rmoisg(i,j,iwater)=rmoisg(i,j,iwater)-ivavm*dcmoisg*dtime

              if (rmoisg(i,j,iwater).lt.0d0) then
                rmoisg(i,j,iwater)= 0d0
              endif
              ! [TODO] -> pas conservatif ajout à cormois?

! ***         calculate moist adiabatic lapse rate

              pot2g=temp2g(i,j)/potfac1
              pot4g=temp4g(i,j)/potfac2
              teta(i,j)=0.5d0*(pot2g-pot4g)

!              t500 = ec_levtempgp(plevel(2),i,j)
              t500 = 0.5d0 * (temp2g(i,j) + temp4g(i,j))

              qsat500=ec_qsat(plevel(2),t500)

              hulp=1d0 + fachulp*qsat500/(t500**2)
              gams(i,j)=gamd*(1+(rlatvap*qsat500)/(rgas*t500))/hulp
              tetacr(i,j)=0.5d0*t500*facteta*(gamd-gams(i,j))

#if ( DOWNSCALING == 2 )
              if (sub_grid_notflat(i,j) > 0.0d0) then ! or stability_coarse==T
!!                 write(*,*) "stability_coarse", i,j,corain(i,j)
!!     &            , corain_sg(i,j,1), corain_sg(i,j,2)
!!     &            , corain_sg(i,j,3), corain_sg(i,j,4)
                 corain_sg(i,j,:)=corain_sg(i,j,:)+corain(i,j,iwater)
                 cosnow_sg(i,j,:)=cosnow_sg(i,j,:)+cosnow(i,j,iwater)
              endif

#endif

            endif ! on rmoisg > qmax

#if ( SYS_CLOCK > 0 )
            call system_clock(count=count_stop)
            write(*,*) '("Time/2bis = ",i6," counts")', count_stop-count_strt
#endif

#if ( DOWNSCALING == 2 )
            else ! (sub_grid is not flat and iconvn=whatever, do the downscaling procedure ...)

#if ( SYS_CLOCK > 0 )
            call system_clock(count=count_strt)
#endif

!dmr&afq --- [???] Are the following lines really useful? (also for the non-downscaled version)
          pot2g=temp2g(i,j)/potfac1
          pot4g=temp4g(i,j)/potfac2
          teta(i,j)=0.5d0*(pot2g-pot4g)

! ***       dry adiabatic lapse rate
          tetacr(i,j)=0d0

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
! dmr   Reset work variables ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

          dcmoisg = 0d0
          dcmoisg_accum = 0d0

! afq -- sub-grid cells initialised to coarse-grid values
! afq -- at the end we do a sum, so we need to fill up only the effective subgrid


! dmr I have modified the sum so that it operates only on the effective one ...
!     [TODO] --- remove this init?
          thforg2_grisli(:) = 0d0
          temp2g_grisli(:) = 0d0
          temp4g_grisli(:) = 0d0
          vhforg1_grisli(:) = 0d0
          vhforg2_grisli(:) = 0d0
! dmr    --- until here

          lokaal_npoints = nbpointssg(i,j)

          thforg2_grisli(1:lokaal_npoints) = thforg2(i,j)
          temp2g_grisli(1:lokaal_npoints) = temp2g(i,j)
          temp4g_grisli(1:lokaal_npoints) = temp4g(i,j)
          vhforg1_grisli(1:lokaal_npoints) = vhforg1(i,j)
          vhforg2_grisli(1:lokaal_npoints) = vhforg2(i,j)

          rmoisg_grisli(1:lokaal_npoints) = rmoisg_interp_sg(i,j,1:lokaal_npoints)

          temp2g_grisli_diag(1:lokaal_npoints) = temp2g_interp_sg(i,j,1:lokaal_npoints)
          temp4g_grisli_diag(1:lokaal_npoints) = temp4g_interp_sg(i,j,1:lokaal_npoints)

          dyrain_grisli(:)  = 0.0d0
          dysnow_grisli(:)  = 0.0d0
          corain_grisli(:)  = 0.0d0
          cosnow_grisli(:)  = 0.0d0

          rain_cg = 0.0d0
          snow_cg = 0.0d0

          I_did_something = 0

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
! dmr   Start loop on sub-grid cells ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

          do n_point = 1, lokaal_npoints

            iconvn_grisli = 0

            do while (iconvn_grisli.lt.10) ! afq -- 30 continue

! afq, need to be recomputed for each convective loop:
            pot2g=temp2g_grisli(n_point)/potfac1
            pot4g=temp4g_grisli(n_point)/potfac2
            teta_grisli=0.5d0*(pot2g-pot4g)
! ***       dry adiabatic lapse rate
            tetacr_grisli=0d0

            dcmoisg_hr = 0.0d0

            ind_low = index_low_sg(i,j,n_point)
            weight_low = weights_low_sg(i,j,n_point)

! NOTA (dmr) use of reduce function allows to spare the computation of 9 out of 11 vertical levels
            success = ec_ptmoisgp_dreduce(temp2g_grisli_diag(n_point),temp4g_grisli_diag(n_point), geopg(i,j,2),tmount_d,      &  
                                          qmax_d,dqmdt_d,qmount_virt(:),ind_low,ind_low+1)

            qmax_d(ind_low:ind_low+1) = relhmax * qmax_d(ind_low:ind_low+1)

! [DELETE] VARIABLES qmax_grisli, dqmdt_grisli, tsurf_grisli
! [REPLACE] with:    qmax_hr    , dqmdt_hr    , tsurf_hr

            ! ... but we only use two of them ... can spare 9 levels computation!
            qmax_hr =  qmax_d(ind_low)*weight_low + qmax_d(ind_low+1) * (1.d0 - weight_low)
            dqmdt_hr = dqmdt_d(ind_low)*weight_low +  dqmdt_d(ind_low+1)*(1.d0 - weight_low)

            tsurf_hr = tsurf_d_interp_sg(ind_low,i,j,n_point)*weight_low                                                        &
                     + tsurf_d_interp_sg(ind_low+1,i,j,n_point)*(1.d0 - weight_low)

            tsurf_sg(i,j,n_point) = tsurf_hr

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
! ***         if the air is supersaturated initially, this is due to
! ***         large scale convergence of moisture and large scale
! ***         temperature changes. Excessive moisture
! ***         is then removed as dynamic rain
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

               if (rmoisg_grisli(n_point).gt.qmax_hr.and.(dcmoisg_accum.lt.(rmoisg(i,j,iwater)/(ivavm*dtime)))) then

                  I_did_something = I_did_something + 1

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
! ***     calculate rain reduction factor to account for increased
! ***     moisture capacity due to latent heat release
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

! redrain is dependent on dqmdt ...
              redrain=1d0+dqmdt_hr*relhmax*rowat*rlatvap*grav/(cpair*dp2)

              dcmoisg_hr = (rmoisg_grisli(n_point)-qmax_hr)/(redrain*dtime)

              dcmoisg_accum = dcmoisg_accum + dcmoisg_hr * area_sg(i,j,n_point)/area_sg_onecb(i,j)

              if (dcmoisg_accum.ge.(rmoisg(i,j,iwater)/(ivavm*dtime))) then

                 dcmoisg_hr=( rmoisg(i,j,iwater)/(ivavm*dtime) - dcmoisg_accum ) * ( area_sg_onecb(i,j) / area_sg(i,j,n_point) )
                 dcmoisg_accum = 1.e9
              endif

              if (iconvn_grisli.eq.0) then

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Temperatures are positive, we are dealing with rain ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

                if (tsurf_hr.ge.tzero + 0.1d0 ) then
                  dyrain_grisli(n_point)= dyrain_grisli(n_point) + dcmoisg_hr

                  if (dyrain_grisli(n_point).gt.drainm_grisli) then
                    dcmoisg_hr=drainm_grisli-dyrain_grisli(n_point) +dcmoisg_hr
                    dyrain_grisli(n_point)=drainm_grisli
                  endif

                  thforg2_grisli(n_point)=thforg2_grisli(n_point)+factemv*dcmoisg_hr/dp2
                  vhforg2_grisli(n_point)=vhforg2_grisli(n_point)+factemv*dcmoisg_hr/dp2

                else ! tsurf < tzero

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Temperatures are negative, we are dealing with dynamic snow ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

                  dysnow_grisli(n_point)=dysnow_grisli(n_point) + dcmoisg_hr

                  if (dysnow_grisli(n_point).gt.drainm_grisli) then
                    dcmoisg_hr=drainm_grisli-dysnow_grisli(n_point) + dcmoisg_hr
                    dysnow_grisli(n_point)=drainm_grisli
                  endif

                  thforg2_grisli(n_point)=thforg2_grisli(n_point) + factems*dcmoisg_hr/dp2
                  vhforg2_grisli(n_point)=vhforg2_grisli(n_point) + factems*dcmoisg_hr/dp2

                endif  ! on tsurf.ge.tzero, dynamic part

         else  ! dyn precipitation done (iconvn = 0), do convective:

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Temperatures are positive, we are dealing with rain ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

                if (tsurf_hr.ge.tzero + 0.1d0 ) then


                  corain_grisli(n_point)=corain_grisli(n_point) + dcmoisg_hr

                  if (corain_grisli(n_point).gt.crainm_grisli) then
                    dcmoisg_hr=crainm_grisli-corain_grisli(n_point)+dcmoisg_hr
                    corain_grisli(n_point)=crainm_grisli
                  endif

                  temp4g_grisli(n_point)=temp4g_grisli(n_point) + factemv*dcmoisg_hr*dtime/dp2
                  vhforg2_grisli(n_point)=vhforg2_grisli(n_point) + factemv*dcmoisg_hr/dp2

                else
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!       Temperatures are negative, we are dealing with snow ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

                  cosnow_grisli(n_point)=cosnow_grisli(n_point) + dcmoisg_hr
                  
                  if (cosnow_grisli(n_point).gt.crainm_grisli) then
                    dcmoisg_hr=crainm_grisli-cosnow_grisli(n_point) + dcmoisg_hr
                    cosnow_grisli(n_point)=crainm_grisli
                  endif
                  temp4g_grisli(n_point)=temp4g_grisli(n_point)  + factems*dcmoisg_hr*dtime/dp2
                  vhforg2_grisli(n_point)=vhforg2_grisli(n_point)+ factems*dcmoisg_hr/dp2

                endif ! on tsurf > tzero

             endif ! on iconvn_grisli == 0

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
!   From the above code, corain, cosnow and dcmoisg are determined
!    we are ready to modify the moisture content in the sub-grid
!    (rmoisg_grisli)
!
! ***         moisture changes due to sub-grid convective rain
!
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

              rmoisg_grisli(n_point)=rmoisg_grisli(n_point)-ivavm*dcmoisg_hr*dtime

              if (rmoisg_grisli(n_point).lt.0d0) then

                 rmoisg_grisli(n_point) = 0d0
! afq -- not a problem: the var. used for water cons. is rmoisg, not rmoisg_grisli
              endif

              pot2g=temp2g_grisli(n_point)/potfac1
              pot4g=temp4g_grisli(n_point)/potfac2
              teta_grisli=0.5d0*(pot2g-pot4g)

              t500 = 0.5d0 * (temp2g_grisli(n_point) + temp4g_grisli(n_point))

              qsat500=ec_qsat(plevel(2),t500)

              hulp=1d0 + fachulp*qsat500/(t500**2)
              gams_grisli=gamd*(1+(rlatvap*qsat500)/(rgas*t500))/hulp
              tetacr_grisli=0.5d0*t500*facteta*(gamd-gams_grisli)

           endif    ! on rmoisg_grisli > qmax_grisli


           if (teta_grisli .lt. (tetacr_grisli-0.1)) then

               pot2g=(dp1*temp2g_grisli(n_point)+dp2*temp4g_grisli(n_point)+dp2*potfac2*2.*tetacr_grisli)/(potfac1*dp1+potfac2*dp2)

               pot4g=pot2g - 2.d0*tetacr_grisli
               temp2go=temp2g_grisli(n_point)
               temp4go=temp4g_grisli(n_point)
               temp2g_grisli(n_point)=pot2g*potfac1
               temp4g_grisli(n_point)=pot4g*potfac2
               vhforg1_grisli(n_point)=vhforg1_grisli(n_point) + (temp2g_grisli(n_point)-temp2go)/dtime
               vhforg2_grisli(n_point)=vhforg2_grisli(n_point) + (temp4g_grisli(n_point)-temp4go)/dtime

               iconvn_grisli=iconvn_grisli+1

            else ! if stable: we do nothing and move on to next sub-grid point
               exit
            endif

         enddo ! on do while iconvn_grisli.lt.10
            
         if (iconvn_grisli.ge.10) then
            write(*,*) "SUB-GRID CONVEC DOES NOT REACH EQUIL"
            write(*,*) i,j,n_point,rmoisg_grisli(n_point), rmoisg(i,j,iwater),qmax_hr
            !read(*,*)
         endif

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
! dmr  Here I am still looping on the grisli cells ...
!        This part is parsed only once per grisli cell
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

             dyrain_sg(i,j,n_point) = dyrain_sg(i,j,n_point) + dyrain_grisli(n_point)
             dysnow_sg(i,j,n_point) = dysnow_sg(i,j,n_point) + dysnow_grisli(n_point)

             corain_sg(i,j,n_point) = corain_sg(i,j,n_point) + corain_grisli(n_point)
             cosnow_sg(i,j,n_point) = cosnow_sg(i,j,n_point) + cosnow_grisli(n_point)

! afq -- We add the rel. humidity output (normally done in atmmois):
             qmax_hr=qmax_hr/relhmax ! corrected for relhmax for consistency with atmmois
             if(qmax_hr.gt.0.) then
                relhum_sg(i,j,n_point) = min(rmoisg_grisli(n_point)/qmax_hr,1.)
             else
                relhum_sg(i,j,n_point) = 0. !as in atmmois_mod.f...
             endif
          enddo ! on nbISMpoints ...

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
! dmr  Back on the grid cell of ECbilt ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

          if ( I_did_something > 0 ) then ! useless to compute if none of the cell is modified

              ! here define the ECBilt variables as being the weighted sum of the sub-grid cell
              ! variables ...

            rain_cg=(SUM(dyrain_grisli(1:lokaal_npoints)*area_sg(i,j,1:lokaal_npoints),dim=1)+SUM(corain_grisli(1:lokaal_npoints) &
                       *area_sg(i,j,1:lokaal_npoints),dim=1) ) / area_sg_onecb(i,j)
            snow_cg=(SUM(dysnow_grisli(1:lokaal_npoints)*area_sg(i,j,1:lokaal_npoints),dim=1)+SUM(cosnow_grisli(1:lokaal_npoints) &
                       *area_sg(i,j,1:lokaal_npoints),dim=1) ) / area_sg_onecb(i,j)

            dcmoisg    = rain_cg + snow_cg

            dyrain(i,j,iwater) = dyrain(i,j,iwater)+SUM(dyrain_grisli(1:lokaal_npoints)*area_sg(i,j,1:lokaal_npoints),dim=1) & 
                    /area_sg_onecb(i,j)
            dysnow(i,j,iwater) = dysnow(i,j,iwater)+SUM(dysnow_grisli(1:lokaal_npoints)*area_sg(i,j,1:lokaal_npoints),dim=1) & 
                    /area_sg_onecb(i,j)

            corain(i,j,iwater) = corain(i,j,iwater)+SUM(corain_grisli(1:lokaal_npoints)*area_sg(i,j,1:lokaal_npoints),dim=1) & 
                    /area_sg_onecb(i,j)
            cosnow(i,j,iwater) = cosnow(i,j,iwater)+SUM(cosnow_grisli(1:lokaal_npoints)*area_sg(i,j,1:lokaal_npoints),dim=1) & 
                    /area_sg_onecb(i,j)

            temp2g(i,j) = SUM(temp2g_grisli(1:lokaal_npoints)*area_sg(i,j,1:lokaal_npoints),dim=1) /area_sg_onecb(i,j)

            temp4g(i,j) = SUM(temp4g_grisli(1:lokaal_npoints)*area_sg(i,j,1:lokaal_npoints),dim=1) /area_sg_onecb(i,j)

            thforg2(i,j) = SUM(thforg2_grisli(1:lokaal_npoints)*area_sg(i,j,1:lokaal_npoints),dim=1) /area_sg_onecb(i,j)

            vhforg1(i,j) = SUM(vhforg1_grisli(1:lokaal_npoints)*area_sg(i,j,1:lokaal_npoints),dim=1) /area_sg_onecb(i,j)

            vhforg2(i,j) = SUM(vhforg2_grisli(1:lokaal_npoints)*area_sg(i,j,1:lokaal_npoints),dim=1) /area_sg_onecb(i,j)

            rmoisg(i,j,iwater)=rmoisg(i,j,iwater)-ivavm*dcmoisg*dtime

! ***         calculate moist adiabatic lapse rate

            pot2g=temp2g(i,j)/potfac1
            pot4g=temp4g(i,j)/potfac2
            teta(i,j)=0.5d0*(pot2g-pot4g)

            t500 = 0.5d0 * (temp2g(i,j) + temp4g(i,j))

            qsat500=ec_qsat(plevel(2),t500)

            hulp=1d0 + fachulp*qsat500/(t500**2)
            gams(i,j)=gamd*(1+(rlatvap*qsat500)/(rgas*t500))/hulp
            tetacr(i,j)=0.5d0*t500*facteta*(gamd-gams(i,j))

! dmr should be moved above (sub grid) !!!             endif ! on rmoisg > qmax
            endif ! on I_did_something

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|
! dmr  End of the downscaling procedure ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2---------3-|

#if ( SYS_CLOCK > 0 )
            call system_clock(count=count_stop)
            write(*,*) '("Time/2third = ",i6," counts")', count_stop - count_strt
#endif

            endif ! on sub_grid being flat or not ...
#endif /* on DOWNSCALING == 2 */

            if (teta(i,j) .lt. (tetacr(i,j)-0.1)) then


              pot2g=(dp1*temp2g(i,j)+dp2*temp4g(i,j)+dp2*potfac2*2.*tetacr(i,j))/(potfac1*dp1+potfac2*dp2)

              pot4g=pot2g - 2.d0*tetacr(i,j)
              temp2go=temp2g(i,j)
              temp4go=temp4g(i,j)
              temp2g(i,j)=pot2g*potfac1
              temp4g(i,j)=pot4g*potfac2
              vhforg1(i,j)=vhforg1(i,j) + (temp2g(i,j)-temp2go)/dtime
              vhforg2(i,j)=vhforg2(i,j) + (temp4g(i,j)-temp4go)/dtime
              iconvn=iconvn+1

              ! afq -- if (iconvn.lt.10) then

                if (dyrain(i,j,iwater).eq.drainm) then
                  write(error_id,*) 'in latlon ',i,j, ' dyrain'
                  call ec_error(122)
                  exit !afq -- goto 20

                else if (corain(i,j,iwater).eq.crainm) then
                  write(error_id,*) 'in latlon ',i,j, ' corain'
                  call ec_error(122)
                  exit !afq -- goto 20

                else if (dysnow(i,j,iwater).eq.drainm) then
                  write(error_id,*) 'in latlon ',i,j, ' dysnow'
                  call ec_error(122)
                  exit !afq -- goto 20

                else if (cosnow(i,j,iwater).eq.crainm) then
                  write(error_id,*) 'in latlon ',i,j, ' cosnow'
                  call ec_error(122)
                  exit !afq -- goto 20
                endif

#if ( DOWNSCALING == 2 )
! afq -- with the downscaling we can not improve large scale stability...
              if ( (sub_grid_notflat(i,j)>0.0d0).and.(I_did_something.gt.0) ) then
                 stability_coarse = .True.
              endif
#endif
              ! afq --goto 10
           else ! teta(i,j) >= (tetacr(i,j)-0.1)
              exit
           endif

        enddo !on do while iconvn<=10


        if (iconvn.ge.10) then
           write(error_id,*) 'warning in lat-lon point: ',i,j
           write(error_id,*) temp2g(i,j),temp4g(i,j)
           call ec_error(120)
           ! afq -- goto 20
        endif
             

! afq --  20      continue ! sortir de cette là ...

          if (iconvn.gt.ncmax) ncmax=iconvn

          torain(i,j,iwater)= dyrain(i,j,iwater) + corain(i,j,iwater)
          tosnow(i,j,iwater)= dysnow(i,j,iwater) + cosnow(i,j,iwater)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

        enddo ! spatial loop
      enddo   ! spatial loop


#if ( DOWNSCALING == 2 )
     do n_point=1,max_nb_points
     do j = 1,nlon
     do i = 1,nlat
          torain_sg(i,j,n_point) = dyrain_sg(i,j,n_point) + corain_sg(i,j,n_point)
          tosnow_sg(i,j,n_point) = dysnow_sg(i,j,n_point) + cosnow_sg(i,j,n_point)
      enddo
      enddo
      enddo
          
#endif


#if ( DOWNSTS == 1 )
     do n_point=1,nb_levls
     do j = 1,nlon
     do i = 1,nlat
          torain_d(i,j,n_point) = dyrain_d(i,j,n_point) + corain_d(i,j,n_point)
          tosnow_d(i,j,n_point) = dysnow_d(i,j,n_point) + cosnow_d(i,j,n_point)
      enddo
      enddo
      enddo          
#endif


!
      return
      end subroutine ec_convec_two
#endif
   
      end module ec_convec_mod
