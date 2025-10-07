!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:47 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:47 CET 2009

      SUBROUTINE icdadv(xjour)
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!  This routine advcect sea ice properties.
!
!  modif : 29/11/99

      use global_constants_mod, only: dblp=>dp, ip           

      use const_mod, only: one, zero

      use para0_mod, only: imax, jmax
      use moment_mod, only: sxg, syg, sxxg, syyg, sxyg, sxn, syn, sxxn
     &                    , syyn, sxyn, sxa, sya, sxxa, syya, sxya, sxc0
     &                    , syc0, sxxc0, syyc0, sxyc0, sxc1, syc1, sxxc1
     &                    , syyc1, sxyc1, sxc2, syc2, sxxc2, syyc2
     &                    , sxyc2, sxst, syst, sxxst, syyst, sxyst

      use bloc0_mod, only: tmu, ku2, dts, ks2, js1, js2, is1, is2, jeq
     &             , jcl1, jcl2, ims1, ims2
      use bloc_mod, only: idyn
      use ice_mod, only: hnbq, hgbq, albq, alcd, tbq, qstobq, tfsn
     &           , tfsg, xlg, alcr, acrit, hndif, tfu
      use dynami_mod, only: area, ug, vg, bound, dxc1, dxc2, hnm, hgm
     &              , iameth, dfhu, dfhv
      use newunit_clio_mod, only: clio3_out_id      


      implicit none
!
      real(kind=dblp), dimension(imax,jmax) :: ut, vt, sm, s1, amskx, amsky
     &        , difhx, difhy, fld0, fld1, s0g, s0n, s0a, s0c0, s0c1
     &        , s0c2, s0st
     
!age &         ,s0an(imax,jmax),s0ag(imax,jmax)

!--- loads of locales

       integer(kind=ip):: i, ip1, j, jour, jp1, k, nitad, jj
       real(kind=dblp) :: acrith, albq_prev, cfl, dta, onelo, usnit
     &                 , vbord, xjour, zeps0, zeps1, zerolo, zignm
     &                 , zindb, zindg, zindhe, zindn, zusvog, zusvon



!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!--1. Initialisation                                                   |
!-----------------------------------------------------------------------
!
      jour=int(xjour)
      zeps0 = 1.0d-16
      zeps1 = 1.0d-20
      onelo =  one
      zerolo = zero
      do j=1,jmax
       do i=1,imax
        sm(i,j)=area(i,j)
       enddo
      enddo

      if (idyn.eq.1) then
!
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!--2.Advection of sea ice properties.                                  |
!-----------------------------------------------------------------------
!
!--2.1. Computation of velocities for advection.
!-----------------------------------------------
!
! vbord factor between 1 and 2 to take into account slip
! or no-slip boundary conditions.
!
        vbord = 1.0+(1.0-bound)
        do 80 j=1,jmax-1
          do 70 i=1,imax
            ip1     =  i+1-(imax-2)*(i/imax)
!           ut(i,j) =  0.5*(ug(ip1,j)+ug(ip1,j+1))
!           vt(i,j) =  0.5*(vg(i,j+1)+vg(ip1,j+1))
            ut(i,j) = (ug(ip1,j)+ug(ip1,j+1))/
     &          (max(tmu(ip1,j,ku2)+tmu(ip1,j+1,ku2),vbord))
            vt(i,j) =  (vg(i,j+1)+vg(ip1,j+1))/
     &          (max(tmu(i,j+1,ku2)+tmu(ip1,j+1,ku2),vbord))
70        continue
80        continue

!
!  2.1.1. Boundary conditions at the top of the grids.
!
        do 90 i=1,imax
          ut(i,jmax) = 0.0
          vt(i,jmax) = 0.0
90      continue
!
!  2.1.2. CFL test for stability.
!
         cfl    = 0.0
         do 110 j=1,jmax
           do 100 i=1,imax
            ip1  = (i+1)-(imax-2)*(i/imax)
             cfl  = max(cfl,(abs(ut(i,j))*dts(ks2))/
     &             (0.5*(dxc1(i,j)+dxc1(ip1,j))))
100          continue
110        continue
         do 130 j=1,jmax-1
          jp1 = (j+1)
           do 120 i=1,imax
             cfl = max(cfl,(abs(vt(i,j))*dts(ks2))/
     &            (0.5*(dxc2(i,j)+dxc2(i,jp1))))
120          continue
130        continue
         if (cfl.gt.0.5) then
          write(clio3_out_id,'(a31,i3,a14,f10.6)')
     &          'violation of cfl criterion the ',jour,'th day, cfl = ',cfl
        endif
!
!--2.2. Transported properties.
!-------------------------------
!
        do 150 j=1,jmax
          do 140 i=1,imax
!
!  Snow volume.
!
            s0n(i,j)  = hnm(i,j)*area(i,j)
!
!  Ice volume.
!
            s0g(i,j)  = hgm(i,j)*area(i,j)
!
!  Surface covered by ice.
!
            s0a(i,j)  = (1.0-albq(i,j))*area(i,j)
!
!  Heat content of the snow layer.
!
            s0c0(i,j) = (tbq(i,j,1)/tfsn)*s0n(i,j)
!
!  Heat content of the first ice layer.
!
            s0c1(i,j) = (tbq(i,j,2)/tfsg)*s0g(i,j)
!
!  Heat content of the second ice layer.
!
            s0c2(i,j) = (tbq(i,j,3)/tfsg)*s0g(i,j)
!
!  Heat reservoir for brine pockets.
!
            s0st(i,j) = (qstobq(i,j)/xlg)*s0a(i,j)
!
!  Age of snow
!
!age        s0an(i,j) = agen(i,j)*s0n(i,j)
!
!  Age of ice
!
!age        s0ag(i,j) = ageg(i,j)*s0g(i,j)
!
140       continue
150        continue
!
!--2.3. Calls to advection and diffusion routines.
!-------------------------------------------------
!
!  Advection.
!
!        If ice drift field is too fast, use an
!        appropriate time step for advection.
!
        nitad = 1+int(max(zero,sign(one,cfl-0.5)))
        usnit = 1.0/real(nitad)
        dta   = dts(ks2)
        if (iameth.eq.1) then
          do k=1,nitad
!
            do j=1,jmax
              do i=1,imax
                s1(i,j) = s0g(i,j)
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0g,s1)
            do j=1,jmax
              do i=1,imax
               s1(i,j) = 0.5*(s0g(i,j)+s1(i,j))
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0g,s1)
            do j=1,jmax
              do i=1,imax
                s0g(i,j) = s1(i,j)
             enddo
            enddo
!
            do j=1,jmax
              do i=1,imax
                s1(i,j) = s0n(i,j)
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0n,s1)
            do j=1,jmax
              do i=1,imax
                s1(i,j) = 0.5*(s0n(i,j)+s1(i,j))
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0n,s1)
            do j=1,jmax
              do i=1,imax
                s0n(i,j) = s1(i,j)
              enddo
            enddo
!
            do j=1,jmax
              do i=1,imax
                s1(i,j) = s0a(i,j)
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0a,s1)
            do j=1,jmax
              do i=1,imax
                s1(i,j) = 0.5*(s0a(i,j)+s1(i,j))
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0a,s1)
            do j=1,jmax
              do i=1,imax
                s0a(i,j) = s1(i,j)
              enddo
            enddo
!
            do j=1,jmax
              do i=1,imax
                s1(i,j) = s0c0(i,j)
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0c0,s1)
            do j=1,jmax
              do i=1,imax
                s1(i,j) = 0.5*(s0c0(i,j)+s1(i,j))
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0c0,s1)
            do j=1,jmax
              do i=1,imax
                s0c0(i,j) = s1(i,j)
              enddo
            enddo
!
            do j=1,jmax
              do i=1,imax
                s1(i,j) = s0c1(i,j)
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0c1,s1)
            do j=1,jmax
              do i=1,imax
                s1(i,j) = 0.5*(s0c1(i,j)+s1(i,j))
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0c1,s1)
            do j=1,jmax
              do i=1,imax
                s0c1(i,j) = s1(i,j)
              enddo
            enddo
!
            do j=1,jmax
              do i=1,imax
                s1(i,j) = s0c2(i,j)
                    enddo
            enddo
            call adv(dta*usnit,ut,vt,s0c2,s1)
            do j=1,jmax
              do i=1,imax
                s1(i,j) = 0.5*(s0c2(i,j)+s1(i,j))
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0c2,s1)
            do j=1,jmax
              do i=1,imax
                s0c2(i,j) = s1(i,j)
              enddo
            enddo
!
            do j=1,jmax
              do i=1,imax
                s1(i,j) = s0st(i,j)
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0st,s1)
            do j=1,jmax
              do i=1,imax
                s1(i,j) = 0.5*(s0st(i,j)+s1(i,j))
              enddo
            enddo
            call adv(dta*usnit,ut,vt,s0st,s1)
            do j=1,jmax
              do i=1,imax
                s0st(i,j) = s1(i,j)
              enddo
            enddo
!
!age        do j=1,jmax
!age          do i=1,imax
!age            s1(i,j) = s0an(i,j)
!age          enddo
!age        enddo
!age        call adv(dta*usnit,ut,vt,s0an,s1)
!age        do j=1,jmax
!age          do i=1,imax
!age            s1(i,j) = 0.5*(s0an(i,j)+s1(i,j))
!age          enddo
!age        enddo
!age        call adv(dta*usnit,ut,vt,s0an,s1)
!age        do j=1,jmax
!age          do i=1,imax
!age            s0an(i,j) = s1(i,j)
!age          enddo
!age        enddo
!
!age         do j=1,jmax
!age          do i=1,imax
!age            s1(i,j) = s0ag(i,j)
!age          enddo
!age        enddo
!age        call adv(dta*usnit,ut,vt,s0ag,s1)
!age        do j=1,jmax
!age          do i=1,imax
!age            s1(i,j) = 0.5*(s0ag(i,j)+s1(i,j))
!age          enddo
!age        enddo
!age        call adv(dta*usnit,ut,vt,s0ag,s1)
!age        do j=1,jmax
!age          do i=1,imax
!age            s0ag(i,j) = s1(i,j)
!age          enddo
!age        enddo
!
          enddo
        else
          if (mod(jour,2).eq.0) then
            do k=1,nitad
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0g,sxg,sxxg,syg,syyg,sxyg)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0g,sxg,sxxg,syg,syyg,sxyg)
!
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0n,sxn,sxxn,syn,syyn,sxyn)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0n,sxn,sxxn,syn,syyn,sxyn)
!
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0a,sxa,sxxa,sya,syya,sxya)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0a,sxa,sxxa,sya,syya,sxya)
!
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0c0,sxc0,sxxc0,syc0,syyc0,sxyc0)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0c0,sxc0,sxxc0,syc0,syyc0,sxyc0)
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0c1,sxc1,sxxc1,syc1,syyc1,sxyc1)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0c1,sxc1,sxxc1,syc1,syyc1,sxyc1)
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0c2,sxc2,sxxc2,syc2,syyc2,sxyc2)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0c2,sxc2,sxxc2,syc2,syyc2,sxyc2)
              call advx(dta*usnit,ut,onelo,
     &                  sm,s0st,sxst,sxxst,syst,syyst,sxyst)
              call advy(dta*usnit,vt,zerolo,
     &                  sm,s0st,sxst,sxxst,syst,syyst,sxyst)
!age          call advx(dta*usnit,ut,onelo,
!age &                  sm,s0an,sxagn,sxxagn,syagn,syyagn,sxyagn)
!age          call advy(dta*usnit,vt,zerolo,
!age &                  sm,s0an,sxagn,sxxagn,syagn,syyagn,sxyagn)
!age          call advx(dta*usnit,ut,onelo,
!age &                  sm,s0ag,sxagg,sxxagg,syagg,syyagg,sxyagg)
!age          call advy(dta*usnit,vt,zerolo,
!age &                  sm,s0ag,sxagg,sxxagg,syagg,syyagg,sxyagg)
            enddo
          else
            do k=1,nitad
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0g,sxg,sxxg,syg,syyg,sxyg)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0g,sxg,sxxg,syg,syyg,sxyg)
!
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0n,sxn,sxxn,syn,syyn,sxyn)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0n,sxn,sxxn,syn,syyn,sxyn)
!
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0a,sxa,sxxa,sya,syya,sxya)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0a,sxa,sxxa,sya,syya,sxya)
!
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0c0,sxc0,sxxc0,syc0,syyc0,sxyc0)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0c0,sxc0,sxxc0,syc0,syyc0,sxyc0)
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0c1,sxc1,sxxc1,syc1,syyc1,sxyc1)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0c1,sxc1,sxxc1,syc1,syyc1,sxyc1)
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0c2,sxc2,sxxc2,syc2,syyc2,sxyc2)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0c2,sxc2,sxxc2,syc2,syyc2,sxyc2)
              call advy(dta*usnit,vt,onelo,
     &                  sm,s0st,sxst,sxxst,syst,syyst,sxyst)
              call advx(dta*usnit,ut,zerolo,
     &                  sm,s0st,sxst,sxxst,syst,syyst,sxyst)
!age          call advy(dta*usnit,vt,onelo,
!age &                  sm,s0an,sxagn,sxxagn,syagn,syyagn,sxyagn)
!age          call advx(dta*usnit,ut,zerolo,
!age &                  sm,s0an,sxagn,sxxagn,syagn,syyagn,sxyagn)
!age          call advy(dta*usnit,vt,onelo,
!age &                  sm,s0ag,sxagg,sxxagg,syagg,syyagg,sxyagg)
!age          call advx(dta*usnit,ut,zerolo,
!age &                  sm,s0ag,sxagg,sxxagg,syagg,syyagg,sxyagg)
            enddo
          endif
        endif

!
!     extra amount of lead opening/closing due to shearing
!     deformation (hibler, 1980, 1984; flato and hibler, 1991,1995;
!     harder and lemke, 1994) and ridging (schulkes, 1995).
!     also stern et al. (1995).
!
      do j=1,jmax
        do i=1,imax
          s0a(i,j) = s0a(i,j)-alcr(i,j)*dta*area(i,j)
        enddo
      enddo
!
!  Diffusion.
!
      do j=js1-1,js2
        jp1 = j+1
        do i=is1(j)-1,is2(j)
          ip1        = (i+1)-(imax-2)*(i/imax)
          amskx(i,j) = (1.0-max(zero,sign(one,-s0a(i,j))))*
     &                         (1.0-max(zero,sign(one,-s0a(ip1,j))))
          amsky(i,j) = (1.0-max(zero,sign(one,-s0a(i,j))))*
     &                         (1.0-max(zero,sign(one,-s0a(i,jp1))))
!
        enddo
      enddo
      do j=1,jmax
         do i=1,imax
            difhx(i,j) = 0.
            difhy(i,j) = 0.
         enddo
      enddo
      do j=js1-1,js2-1
        jp1 = (j+1)
        do i=is1(j)-1,is2(j)
          ip1          = (i+1)
          difhx(i,j)   = amskx(i,j)*dfhu(i,j)
     &                   /(0.5*(dxc1(i,j)+dxc1(ip1,j)))
          difhy(i,j)   = amsky(i,j)*dfhv(i,j)
     &                   /(0.5*(dxc2(i,j)+dxc2(i,jp1)))
        enddo
      enddo
!
      do j=1,jmax
        do i=1,imax
          fld0(i,j) = s0g(i,j)/area(i,j)
        enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
        do i=1,imax
          s0g(i,j) = max(zero,fld1(i,j)*area(i,j))
        enddo
      enddo
!
      do j=1,jmax
        do i=1,imax
          fld0(i,j) = s0n(i,j)/area(i,j)
        enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
        do i=1,imax
          s0n(i,j) = max(zero,fld1(i,j)*area(i,j))
        enddo
      enddo
!
      do j=1,jmax
        do i=1,imax
          fld0(i,j) = s0a(i,j)/area(i,j)
        enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
        do i=1,imax
          s0a(i,j) = max(zero,fld1(i,j)*area(i,j))
        enddo
      enddo
!
      do j=1,jmax
        do i=1,imax
          fld0(i,j) = s0c0(i,j)/area(i,j)
        enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
        do i=1,imax
          s0c0(i,j) = max(zero,fld1(i,j)*area(i,j))
        enddo
      enddo
!
      do j=1,jmax
        do i=1,imax
          fld0(i,j) = s0c1(i,j)/area(i,j)
        enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
        do i=1,imax
          s0c1(i,j) = max(zero,fld1(i,j)*area(i,j))
        enddo
      enddo
!
      do j=1,jmax
        do i=1,imax
          fld0(i,j) = s0c2(i,j)/area(i,j)
        enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
        do i=1,imax
          s0c2(i,j) = max(zero,fld1(i,j)*area(i,j))
        enddo
      enddo
!
      do j=1,jmax
        do i=1,imax
          fld0(i,j) = s0st(i,j)/area(i,j)
        enddo
      enddo
      call diffus(dta,difhx,difhy,fld0,fld1)
      do j=1,jmax
        do i=1,imax
          s0st(i,j) = max(zero,fld1(i,j)*area(i,j))
        enddo
      enddo
!
!age  do j=1,jmax
!age    do i=1,imax
!age      fld0(i,j) = s0an(i,j)/area(i,j)
!age    enddo
!age  enddo
!age  call diffus(dta,difhx,difhy,fld0,fld1)
!age  do j=1,jmax
!age    do i=1,imax
!age      s0an(i,j) = fld1(i,j)*area(i,j)
!age    enddo
!age  enddo
!
!age  do j=1,jmax
!age    do i=1,imax
!age      fld0(i,j) = s0ag(i,j)/area(i,j)
!age    enddo
!age  enddo
!age  call diffus(dta,difhx,difhy,fld0,fld1)
!age  do j=1,jmax
!age    do i=1,imax
!age      s0ag(i,j) = fld1(i,j)*area(i,j)
!age    enddo
!age  enddo

!
!--2.5. Up-dating and limitation of sea ice properties
!       after transport.
!------------------------------------------------------
!
        do 190 j=js1,js2
          zindhe = real(max(0,isign(1,j-jeq)))
          do 180 i=is1(j),is2(j)
!
!  2.5.1. Recover mean values over the grid squares.
!
            s0n(i,j)    = max(zero,s0n(i,j)/area(i,j))
            s0g(i,j)    = max(zero,s0g(i,j)/area(i,j))
            s0a(i,j)    = max(zero,s0a(i,j)/area(i,j))
            s0c0(i,j)   = max(zero,s0c0(i,j)/area(i,j))
            s0c1(i,j)   = max(zero,s0c1(i,j)/area(i,j))
            s0c2(i,j)   = max(zero,s0c2(i,j)/area(i,j))
            s0st(i,j)   = max(zero,s0st(i,j)/area(i,j))
!age        s0an(i,j)   = max(zero,s0an(i,j)/area(i,j))
!age        s0ag(i,j)   = max(zero,s0ag(i,j)/area(i,j))
!
!  2.5.2. Recover in situ values.
!
            zindb       = max(zero,sign(one,s0a(i,j)-1.0e-06))
            acrith      = 1.0-(zindhe*acrit(1)+(1.0-zindhe)*acrit(2))
            s0a(i,j)    = zindb*min(s0a(i,j),acrith)
            hnbq(i,j)   = zindb*(s0n(i,j)/max(s0a(i,j),zeps0))
            hgbq(i,j)   = zindb*(s0g(i,j)/max(s0a(i,j),zeps0))
            zindn       = max(zero,sign(one,hnbq(i,j)-1.0e-06))
            zindg       = max(zero,sign(one,hgbq(i,j)-1.0e-03))
            zindb       = max(zindn,zindg)
            s0a(i,j)    = zindb*s0a(i,j)
            albq_prev   = albq(i,j)
            albq(i,j)   = 1.0-s0a(i,j)
            alcd(i,j)   = (albq(i,j)-albq_prev)/dta
     &                   -alcr(i,j)
            hnbq(i,j)   = zindn*hnbq(i,j)
            hgbq(i,j)   = zindg*hgbq(i,j)
            zusvon      = 1.0/max(hnbq(i,j)*s0a(i,j),zeps0)
            zusvog      = 1.0/max(hgbq(i,j)*s0a(i,j),zeps0)
            zignm       = max(zero,sign(one,hndif-hnbq(i,j)))
            tbq(i,j,1)  = zindn*(zignm*tbq(i,j,1)+(1.0-zignm)*
     &                    min(max(173.15*one,
     &                    tfsn*zusvon*s0c0(i,j)),tfu(i,j)))+
     &                          (1.0-zindn)*tfu(i,j)
            tbq(i,j,2)  = zindg*min(max(173.15*one,tfsg*zusvog*
     &                    s0c1(i,j)),tfu(i,j))+
     &                          (1.0-zindg)*tfu(i,j)
            tbq(i,j,3)  = zindg*min(max(173.15*one,tfsg*zusvog*
     &                    s0c2(i,j)),tfu(i,j))+
     &                          (1.0-zindg)*tfu(i,j)
            qstobq(i,j) = zindb*xlg*s0st(i,j)/max(s0a(i,j),zeps0)
!age        agen(i,j)   = zindn*s0an(i,j)*zusvon
!age        ageg(i,j)   = zindg*s0ag(i,j)*zusvog
180       continue
190        continue
!
!  2.5.3. Raccord cyclique (not neccesary here)
!
        do 200 jj=jcl1,jcl2
           hnbq(ims1-1,jj) = hnbq(ims2,jj)
           hnbq(ims2+1,jj) = hnbq(ims1,jj)
           hgbq(ims1-1,jj) = hgbq(ims2,jj)
           hgbq(ims2+1,jj) = hgbq(ims1,jj)
           albq(ims1-1,jj) = albq(ims2,jj)
           albq(ims2+1,jj) = albq(ims1,jj)
           alcd(ims1-1,jj) = alcd(ims2,jj)
           alcd(ims2+1,jj) = alcd(ims1,jj)
           tbq(ims1-1,jj,1) = tbq(ims2,jj,1)
           tbq(ims2+1,jj,1) = tbq(ims1,jj,1)
           tbq(ims1-1,jj,2) = tbq(ims2,jj,2)
           tbq(ims2+1,jj,2) = tbq(ims1,jj,2)
           tbq(ims1-1,jj,3) = tbq(ims2,jj,3)
           tbq(ims2+1,jj,3) = tbq(ims1,jj,3)
           qstobq(ims1-1,jj) = qstobq(ims2,jj)
           qstobq(ims2+1,jj) = qstobq(ims1,jj)
!age       agen(ims1-1,jj) = agen(ims2,jj)
!age       agen(ims2+1,jj) = agen(ims1,jj)
!age       ageg(ims1-1,jj) = ageg(ims2,jj)
!age       ageg(ims2+1,jj) = ageg(ims1,jj)
 200    continue
!
      endif

      return
!---5----|----5----|----5----|----5----|----5----|----5----|----5----|72--5----|
!- fin de la routine icdadv -
      end
