!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of iLOVECLIM/CLIO
!!      iLOVECLIM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
!!      as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!!
!!      Foobar is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!!      of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!!
!!      You should have received a copy of the GNU General Public License along with Foobar.
!!      If not, see <http://www.gnu.org/licenses/>.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"
!dmr -- Adding the choice of components through the pre-processing options

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: atmos_composition_mod
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module is a wrapper of engtur to allow exportation of needed variables in the context of TIDEMIX
!
!>     @date Creation date: April, 2020, 20th
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module engtur_mod

       use global_constants_mod, only: dblp=>dp
       use para0_mod, only: imax, jmax, kmax
       
       use newunit_clio_mod, only: clio3_out_id
      
       implicit none

       public :: engtur

       private

      ! NOTE_avoid_public_variables_if_possible


#if ( TIDEMIX == 1 )
      real(kind=dblp), dimension(imax,jmax)     , save         :: ensq_tmx, esho_tmx
      real(kind=dblp), dimension(imax,jmax,kmax), save, public :: emix_ini, emix_tmx, avs_tide
#endif

      contains
      

      SUBROUTINE engtur
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  evolution of Turbulent Kinetic Energy  (diffusion is completly implicit)
!  modif : 2019-12-17, dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod, only: zero, one, epsil, rho0
      use para0_mod, only: imax, jmax, kmax
      use para_mod,  only:
      use bloc0_mod, only: avnu0, avnub, bvf, cdbot, is1, iu1, js1, js2, ju1, ju2, ks1, ks2, q2turb, tms, unsdzw, is2, dts &
                         , dzw, tmu, kfs, iu2, avsdz, ust2b, ust2s, phifv, avudz, fqajc, kfu, phifu, u, v                  
#if ( TIDEMIX == 1 )
      use bloc0_mod, only: unsdz, z, zw, hs
#endif      

      use bloc_mod,  only: kajul, lstab, nstart, numit, q2tmin, tm2tur, vlturb, avqdz, avuloc
      use ice_mod,   only:

#if ( TIDEMIX == 1 )
      use ncio,      only: nc_read
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      implicit none


      real(kind=dblp), dimension(kmax), save :: avs2mx
      real(kind=dblp)                 , save :: bvfmix

      real(kind=dblp) :: avsmix,ccdif,ccdt, demi,des,epsil2,prod,tm4u,zinsto,zinztz

      integer(kind=ip) :: i,j,k,km2,kstab

      real(kind=dblp), dimension(imax,kmax) :: ff = 0.0_dblp, aa = 0.0_dblp, bb, cc
      real(kind=dblp), dimension(kmax)      :: cczdt
      real(kind=dblp), dimension(imax,jmax,kmax) :: zaju = 0.0_dblp

#if ( TIDEMIX == 1 )
!~       real(kind=dblp), dimension(imax,jmax)     , save :: ensq_tmx, esho_tmx
!~       real(kind=dblp), dimension(imax,jmax,kmax), save :: emix_ini, emix_tmx, avs_tide
      real(kind=dblp), dimension(imax,jmax) :: zfact, ebot_tmx, ecri_tmx
      real(kind=dblp), dimension(imax,jmax) :: hbot_tmx, hcri_tmx
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!    CCW* coeff relatif a l'interface (= place des W).
!    CCZ* coeff relatif au centre des boites (= place de T,S)
!    CCWI : tient compte du taux Implicite.
!    CCWABS(k) intervient avec |w| ;
!              CCWUP(k) avec scal(k) ; CCWDowN(k) avec scal(k-1)
!  pour le taux de Decentrement calcule a partir du Nb de courant : CCWDUP(k)
!    avec dts(k), CCWDDowN(k) avec dts(k-1) ey CCWDMY(k) avec la moyenne des 2
!    CCZDT**(k) intervient pour le bilan des flux V, boite "k" ;
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      if (numit.le.nstart) then
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        avsmix = avnub(1)
        bvfmix = avnu0(1)
             if (lstab.eq.0) write(clio3_out_id,'(2(A,1PE10.3))') '          AVS +', avsmix, ' si bvf <', bvfmix
        do k=ks1+1,ks2
         avs2mx(k) = unsdzw(k) * avsmix * 0.5
        enddo

#if ( TIDEMIX == 1 )    
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  Initialization of variables needed for new tidal mixing parameterization
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!!
!!              - Read the input data in 6 files :
!!              1. dissipation scaling with squared buoyancy frequency (mixing_power_nsq.nc)
!!              2. dissipation due to shoaling internal tides (mixing_power_sho.nc)
!!              3. bottom-intensified dissipation at critical slopes (mixing_power_cri.nc)
!!              4. bottom-intensified dissipation above abyssal hills (mixing_power_bot.nc)
!!              5. decay scale for critical slope dissipation (decay_scale_cri.nc)
!!              6. decay scale for abyssal hill dissipation (decay_scale_bot.nc)

      call nc_read("inputdata/clio/tidal_mixing.nc","ensq_tmx",ensq_tmx)
      call nc_read("inputdata/clio/tidal_mixing.nc","esho_tmx",esho_tmx)
      call nc_read("inputdata/clio/tidal_mixing.nc","ecri_tmx",ecri_tmx)
      call nc_read("inputdata/clio/tidal_mixing.nc","ebot_tmx",ebot_tmx)
      call nc_read("inputdata/clio/tidal_mixing.nc","hcri_tmx",hcri_tmx)
      call nc_read("inputdata/clio/tidal_mixing.nc","hbot_tmx",hbot_tmx)

!!               Initialize fields with zeros.
      emix_ini(:,:,:) = 0.
      emix_tmx(:,:,:) = 0.
      avs_tide(:,:,:) = 0.

!!               Distribute energy from the cri and bot components in the vertical.
!!               These two vertical distributions are static.
!!
!!               'cri' component: distribute energy over the water column
!!               using an exponential decay above the seafloor.
      zfact(:,:) = 0.
      do j=js1,js2                ! part independent of the level
         do i=is1(j),is2(j)
            zfact(i,j) = rho0 * (  1. - EXP( -hs(i,j) / hcri_tmx(i,j) )  )
            if( zfact(i,j).gt.0 ) zfact(i,j) = ecri_tmx(i,j) / zfact(i,j)
         enddo
      enddo

      do k=ks1+1,ks2              ! complete with the level-dependent part
         emix_ini(:,:,k) = zfact(:,:) *                                             &
            (  EXP( ( -z(k-1) - hs(:,:) ) / hcri_tmx(:,:) )                         &
             - EXP( ( -z(k  ) - hs(:,:) ) / hcri_tmx(:,:) )  )                      &
            * tms(:,:,k-1) / dzw(k)
      enddo

!!                'bot' component: distribute energy over the water column
!!                using an algebraic decay above the seafloor.
      zfact(:,:) = 0.
      do j=js1,js2                ! part independent of the level
         do i=is1(j),is2(j)
            if( hs(i,j).gt.0 )  zfact(i,j) = ebot_tmx(i,j) *                        &
                (  1. +  hbot_tmx(i,j) / hs(i,j)  ) / rho0
         enddo
      enddo

      do k=ks1+1,ks2              ! complete with the level-dependent part
         emix_ini(:,:,k) = emix_ini(:,:,k) + zfact(:,:) *                           &
            (   1. / ( 1. + ( hs(:,:) + z(k-1) ) / hbot_tmx(:,:) )                  &
             -  1. / ( 1. + ( hs(:,:) + z(k  ) ) / hbot_tmx(:,:) )   )              &
            * tms(:,:,k-1) / dzw(k)
      enddo
!!              The other two components depend on the evolving stratification.
!!              Hence they are added to emix_tmx further down, at each iteration.
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! End of specific treatment of first iteration.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      endif

#if ( TIDEMIX == 1 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! Add to emix_tmx the energy from 'nsq' and 'sho' components.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!!               'nsq' component: distribute energy over the water column
!!               as proportional to N^2 (bvf).
      zfact(:,:) = 0.
      do k=ks1+1,ks2              ! part independent of the level
         zfact(:,:) = zfact(:,:) + dzw(k) * max(zero,bvf(:,:,k)) * tms(:,:,k-1)
      enddo

      do j=js1,js2
         do i=is1(j),is2(j)
            if( zfact(i,j).gt.0 ) zfact(i,j) = ensq_tmx(i,j) / ( rho0 * zfact(i,j) )
         enddo
      enddo

      do k=ks1+1,ks2              ! complete with the level-dependent part
         emix_tmx(:,:,k) = emix_ini(:,:,k) +                                       &
            zfact(:,:) * max(zero,bvf(:,:,k)) * tms(:,:,k-1)
      enddo

!!               'sho' component: distribute energy over the water column
!!               as proportional to N (sqrt(bvf)).
      zfact(:,:) = 0.
      do k=ks1+1,ks2              ! part independent of the level
         zfact(:,:) = zfact(:,:) + dzw(k) * sqrt(  max(zero,bvf(:,:,k))  ) * tms(:,:,k-1)
      enddo

      do j=js1,js2
         do i=is1(j),is2(j)
            if( zfact(i,j).gt.0 ) zfact(i,j) = esho_tmx(i,j) / ( rho0 * zfact(i,j) )
         enddo
      enddo

      do k=ks1+1,ks2              ! complete with the level-dependent part
         emix_tmx(:,:,k) = emix_tmx(:,:,k) +                                       &
            zfact(:,:) * sqrt(  max(zero,bvf(:,:,k))  ) * tms(:,:,k-1)
      enddo

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! Derive tidal diffusivity from the tidal mixing energy: avs_tide = (1/6) * emix_tmx / bvf  ,
! and bound the tidal diffusivity by 1.4e-7 and 1e-2.
! Next, update avqdz to incorporate the effect of tidal mixing.
! Update of avudz and avsdz is performed at the end of this subroutine.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      avs_tide(:,:,:) = emix_tmx(:,:,:) / (6. * bvf(:,:,:) + epsil*epsil)
      avs_tide(:,:,:) = min( max(1.4d-7,avs_tide(:,:,:)), 1.d-2 )

      do k=ks1,ks2          
        avqdz(:,:,k) = avqdz(:,:,k) + unsdz(k) * unsdz(k) *                        &
        	( avs_tide(:,:,k+1)*(zw(k+1) - z(k)) +                                 &
        	  avs_tide(:,:,k  )*( z(k)  - zw(k)) )
      enddo
#endif

      if (cdbot.ne.zero) then
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!--Tension de fond : TauX,Y = Cdrag * |V(bot.)| * u,v(bot.) = Coeff * u,v(bot.)
!  calcul du Coeff. (dans avudz(-,-,1))(pour resolution implicite).
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      do j=ju1,ju2
       do i=iu1(j),iu2(j)
         avudz(i,j,1) = sqrt( u(i,j,kfu(i,j))*u(i,j,kfu(i,j))                      &
                            + v(i,j,kfu(i,j))*v(i,j,kfu(i,j)) ) * cdbot
         phifu(i,j) = -avudz(i,j,1)*u(i,j,kfu(i,j))
         phifv(i,j) = -avudz(i,j,1)*v(i,j,kfu(i,j))
       enddo
      enddo

      do j=js1,js2
       do i=is1(j),is2(j)
        tm4u = tmu(i,j,ks2) + tmu(i+1,j+1,ks2)                                     &
             + tmu(i,j+1,ks2) + tmu(i+1,j,ks2) + epsil
        ust2b(i,j)=(sqrt(phifu(i,j)*phifu(i,j)+phifv(i,j)*phifv(i,j))              &
          + sqrt(phifu(i+1,j)*phifu(i+1,j)+phifv(i+1,j)*phifv(i+1,j))              &
          + sqrt(phifu(i,j+1)*phifu(i,j+1)+phifv(i,j+1)*phifv(i,j+1))              &
          + sqrt(phifu(i+1,j+1)*phifu(i+1,j+1)                                     &
                  + phifv(i+1,j+1)*phifv(i+1,j+1))  ) /tm4u
       enddo
      enddo

!-----
      endif


      do k=ks1,ks2
        cczdt(k) = dts(ks2) * unsdzw(k)
      enddo

      if (lstab.eq.0) then
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  1 ) Ajust.Conv. par augmentation de la Diffusion Verticale .
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!- Ajoute avsmix(/dzw) si bvf < bvfmix :
      demi = 0.5
!ajo  do k=1+ks1,ks2
!ajo   do j=js1,js2
!ajo    do i=is1(j),is2(j)
!ajo     avsdz(i,j,k) = avsdz(i,j,k)
!ajo &         + ( avs2mx(k) + sign(avs2mx(k), bvfmix-bvf(i,j,k)) )
!ajo     fqajc(i,j,k) = fqajc(i,j,k)
!ajo &         + ( demi + sign(demi, bvfmix-bvf(i,j,k)) )
!ajo    enddo
!ajo   enddo
!ajo  enddo

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!- Convective adjustment if instability > kajul levels
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      do j=js1,js2
       do i=is1(j),is2(j)

        zinsto = 0.0
        kstab=kajul
        do k=kajul+1,ks2
          zinztz=max(zero,sign(one,bvf(i,j,k) ) )
          kstab=int(zinztz)*k+(1-int(zinztz))*kstab
          zinsto=min(zinsto+zinztz,one)
        enddo
        do k=ks2,kajul+1,-1
!         avsdz(i,j,k) = avsdz(i,j,k)
          zaju(i,j,k) = + ( avs2mx(k) + sign(avs2mx(k), bvfmix-bvf(i,j,k)) ) * (1.0-zinsto)
         fqajc(i,j,k) = fqajc(i,j,k)  + ( demi + sign(demi, bvfmix-bvf(i,j,k)) ) * (1.0-zinsto)
        enddo
!       do k=1+ks1,kajul
        do k=1+ks1,kstab
!        avsdz(i,j,k) = avsdz(i,j,k)
         zaju(i,j,k) = + ( avs2mx(k) + sign(avs2mx(k), bvfmix-bvf(i,j,k)) )
         fqajc(i,j,k) = fqajc(i,j,k) + ( demi + sign(demi, bvfmix-bvf(i,j,k)) ) * (1.0-zinsto)
!--Note :subsurface conv aju are not taken into account in fqajc
        enddo

       enddo
      enddo

      do j=js1,js2
       do i=is1(j),is2(j)
        do k=ks1,ks2
         avqdz(i,j,k)=avqdz(i,j,k)+zaju(i,j,k)
        enddo
       enddo
      enddo

      endif
!--Debut de la boucle externe sur l'indice de latitude j :
      do j=js1,js2
!-----

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  2 ) Building the matrix                                          |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      epsil2=epsil*epsil

!--Boundary condition
      do i=is1(j),is2(j)
       q2turb(i,j,kfs(i,j)) = max(q2tmin,6.51d0*ust2b(i,j))
       q2turb(i,j,ks2+1)    = max(q2tmin,6.51d0*ust2s(i,j))
!      q2turb(i,j,ks2+1)    = max(q2tmin,6.51d0*
!    &                   ust2s(i,j)*(one+sdvt(i,j)))
!    &                   ust2s(i,j)*(one+varfor*.375*sdvt(i,j)))
      enddo

!--

      do k=ks1+1,ks2
       do i=is1(j),is2(j)

!       prod   =  tms(i,j,k-1)* 2 * dts(ks2) * (avudz(i,j,k)*
        prod   =  tms(i,j,k-1)* 2 * dts(ks2) * (avuloc(i,j,k)* dzw(k)*tm2tur(i,j,k)           &
                + max(zero,-avsdz(i,j,k)*dzw(k)*bvf(i,j,k)) )

        des    =  2*dts(ks2)*(sqrt(q2turb(i,j,k))/(16.6*vlturb(i,j,k))                        &
                + max(zero,avsdz(i,j,k)*dzw(k)*bvf(i,j,k)/q2turb(i,j,k)) )

        aa(i,k) = 1.0 + tms(i,j,k-1)*des

        ff(i,k) = q2turb(i,j,k) + prod

       enddo
      enddo

!-N.B.-: avsdz =  Dif.Vert. / dz !!
      do k=ks1+1,ks2
       do i=is1(j),is2(j)
        ccdif = avqdz(i,j,k-1)

!- effet du flux phiz(k) sur q2(k) : aa(k)*q2(k) + bb(k)*q2(k-1)
        ccdt = cczdt(k) * tms(i,j,k-1)
        bb(i,k) = ccdt * (-ccdif )
        aa(i,k) = aa(i,k) + ccdt * (ccdif )
!- effet du flux phiz(k) sur q2(k-1) : aa(k-1)*q2(k-1) + cc(k-1)*q2(k)
        km2=max(k-2,ks1)
        ccdt = cczdt(k-1) * tms(i,j,km2)
        aa(i,k-1) = aa(i,k-1) + ccdt * (ccdif )
        cc(i,k-1) = ccdt * (-ccdif)
       enddo
      enddo

      do i=is1(j),is2(j)
        ff(i,ks2) = ff(i,ks2) + tms(i,j,ks2)*cczdt(ks2)*avqdz(i,j,ks2)*q2turb(i,j,ks2+1)
        aa(i,ks2) = aa(i,ks2) + tms(i,j,ks2)*cczdt(ks2)*avqdz(i,j,ks2)
        aa(i,ks1) = 1.0
      enddo



!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  3 ) decomposition LU - (methode et notation : cf Linear Algebra pp165-167)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!--calcul de 1/alpha(k) dans aa(i,k) et beta(k) dans bb(i,k)
      do i=is1(j),is2(j)
        aa(i,ks1) = 1.0 / aa(i,ks1)
      enddo
      do k=ks1+1,ks2
       do i=is1(j),is2(j)
        bb(i,k) = bb(i,k) * aa(i,k-1)
        aa(i,k) = 1.0 / ( aa(i,k) - bb(i,k) * cc(i,k-1) )
       enddo
      enddo

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  4 ) substitutions avant et arriere .                                |
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!--calcul de g(k) dans ff(i,k) :
      do k=ks1+1,ks2
       do i=is1(j),is2(j)
        ff(i,k) = ff(i,k) - bb(i,k) * ff(i,k-1)
       enddo
      enddo
!--calcul de x(k) dans scal(i,j,k,ns) :
      do i=is1(j),is2(j)
       q2turb(i,j,ks2) = max(q2turb(i,j,ks2+1),ff(i,ks2) * aa(i,ks2))
      enddo

      do k=ks2-1,ks1,-1
       do i=is1(j),is2(j)
        q2turb(i,j,k) =  max(q2tmin,(ff(i,k) - cc(i,k) * q2turb(i,j,k+1)) * aa(i,k)*tms(i,j,k) )
       enddo
      enddo


!--Fin de la boucle externe sur l'indice de latitude j .
      enddo


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!  5 ) Ajust.Conv. par augmentation de la Diffusion Verticale (scal)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      if (lstab.eq.0) then
       do j=js1,js2
        do i=is1(j),is2(j)
         do k=ks1,ks2
          avsdz(i,j,k)=avsdz(i,j,k)+zaju(i,j,k)
         enddo
        enddo
       enddo
      endif

#if ( TIDEMIX == 1 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! Update avsdz and avudz to incorporate the effect of tidal mixing.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!!               Add avs_tide to avsdz and avudz
      do k=ks1+1,ks2
        avsdz(:,:,k) = avsdz(:,:,k) + unsdzw(k) * avs_tide(:,:,k)
        avudz(:,:,k) = avudz(:,:,k) + unsdzw(k) * avs_tide(:,:,k)
      enddo
#endif

      return
      end subroutine engtur
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the subroutine here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module engtur_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

















