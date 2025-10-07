!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:55 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:55 CET 2009

      SUBROUTINE streamfunc(nupdate,usnjm,nstreamout)
!
! this routine is called in SUBROUTINE outave (see: outave.f)
! and computes the meridional overturning streamFUNCTIONs for the 
! MERIDIONAL VELOCITY in 4 basins
!
! also the MERIDIONAL HEAT and MERIDIONAL SALT FLUX are computed in 4 basins
! 
! if GM option is used ==> 'eddy/bolus' (viso) contribution is nonzero 
! diffusion terms are included in calculations
!
! because the fresh water flux is related to vertical velocity at the surface
! the meridional streamFUNCTION uuu(k=kmax) is nonzero
! this has consequences for the meridional heat flux and salt flux
! ===> compute these via an anomaly: T_ref = 1.3 K and S_ref = 34.7 p.s.u.
!
! INPUT:
! 1) nupdate (-1=initialize; 0=compute and update; 1=prepare and write output) 
! 2) normalization factor usnjm(=1 for snapshot, =1/30 for monthly aver. etc.)
! 3) controls output mode in outave.f and is used in SUBROUTINE streaminitout
!
! OUTPUT: 
! output is written for monthly or annual outputs (see outave.f: nstreamout)
! meridional overturning streamFUNCTIONs vs latitude and depth
! j=1 corresponds to true latitude -81.0 degree South, increment: 3 degrees
! 22(=20+extra bottom and surface) depth levels for overturning streamFUNCTION
! 1 depth level for meridional heat flux and salt flux 
! 4 basins (global, Indian, Pacific, Atlantic)
!
! NOTE (see SUBROUTINE streaminitout):
! unit = 38: stream.dat for streamFUNCTION meridional velocity v         
! unit = 39: heatsaltflux.dat for meridional heat flux and salt flux 
!
! ! GrADS str*.ctl files are provided for in streaminitout 
!
! uu0/uuu: v               ==> meridional overturning
! vv0/vvv: v * scal(...,1) ==> meridional heat flux
! ww0/vvv: v * scal(...,2) ==> meridional salt flux
!
! *** NOT included: Bering Strait contributions! ***
!

!! START_OF_USE_SECTION

      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use isoslope_mod
      use ice_mod
      use reper_mod
      
      use stream_mod, only: uuu, vvv, i1zon, i2zon, jmtt, kmtt
     >                    , heatsaltfluxdat_id, streamdat_id
!! END_OF_USE_SECTION

      
c~       real*4 uuu(jmtt,0:kmtt,0:nbsmax)
c~       real*4 vvv(jmtt,0:nbsmax+4)
c~       common/streamflu/ uuu,vvv


c~       integer i1zon(jmax,0:nbsmax),i2zon(jmax,0:nbsmax)
c~       common/streamzon/ i1zon,i2zon


      integer i,j,k,nb,ii,jj,nupdate,nrecstr0,nrecstr1
      real*8 avd,ccydif1,ccydif2,ttm,v00
      real*8 dxsvrdrp,dxsvrdrpT,dxsvrdrpS,usnjm
      real*8 uu0(jmax),vv0(jmax),ww0(jmax)
      save nrecstr0,nrecstr1

      integer(kind=ip) :: jcr, jcr1, jj1, jj2, jsigne, nstreamout



! ------------------------ c nupdate = -1 ==> initialize
      if (nupdate .eq. -1) then



       call streaminit
       call streaminitout(nstreamout)
       nrecstr0 = 0
       nrecstr1 = 0
       do nb=0,nbsmax
        do j=1,jmtt
         do k=0,kmtt
          uuu(j,k,nb) = 0.0
         enddo
        enddo
       enddo
       do nb=0,nbsmax+4
        do j=1,jmtt
         vvv(j,nb) = 0.0
        enddo
       enddo


! -------------------------c nupdate = 0 ==> compute streamFUNCTIONs and update
      elseif (nupdate .eq.  0) then
!
! compute streamFUNCTION for three(four) basins; tms = mask (zero or one)
! uu0/uuu: v               ==> meridional overturning
! vv0/vvv: v * scal(...,1) ==> meridional heat flux
! ww0/www: v * scal(...,2) ==> meridional fresh water flux
!
! begin nb-loop
      do nb=1,nbsmax
! begin k-loop
       do k=1,kmax
        ccydif1 = unsdy * dts(k)
        ccydif2 = 2.0 * unsdy * ahs(k)
        do j=1,jmax
         uu0(j) = 0.d0
         vv0(j) = 0.d0
         ww0(j) = 0.d0
        enddo
        do j=ju1,jeq-1
         do i=i1zon(j,nb),i2zon(j,nb)
          ii = icl(i)
          ttm = tms(ii,j-1,k)*tms(ii,j,k)
          v00 = 0.5*(v(ii,j,k)+v(icl(i+1),j,k)) + viso(ii,j,k)
          avd = min(1.0,abs(ccydif1*smy(ii,j,2)*v00)+alphmi(k)) 
          v00 = ttm * cmx(ii,j,2) * v00
          avd = avd * abs(v00) + ttm * ccydif2 * cmxy(ii,j,2) 
!          avd = 0.0
          uu0(j) = uu0(j) + v00  
          vv0(j) = vv0(j) + v00 * (scal(ii,j-1,k,1)+scal(ii,j,k,1))
     &                    + avd * (scal(ii,j-1,k,1)-scal(ii,j,k,1)) 
          ww0(j) = ww0(j) + v00 * (scal(ii,j-1,k,2)+scal(ii,j,k,2))
     &                    + avd * (scal(ii,j-1,k,2)-scal(ii,j,k,2)) 
         enddo
        enddo
        do j=jeq,jmvlat
         do i=i1zon(j,nb),i2zon(j,nb)
          ii = icl(i)
          jcr = jmaxl(ii,j)
          ttm = tms(ii,jcr-1,k)*tms(ii,jcr,k)
          v00 = 0.5*(v(ii,jcr,k)+v(icl(i+1),jcr,k)) + viso(ii,jcr,k)
          avd = min(1.0,abs(ccydif1*smy(ii,jcr,2)*v00)+alphmi(k)) 
          v00 = ttm * cmx(ii,jcr,2) * v00
          avd = avd * abs(v00) + ttm * ccydif2 * cmxy(ii,jcr,2)
!          avd = 0.0
          uu0(j) = uu0(j) + v00  
          vv0(j) = vv0(j) + v00 * (scal(ii,jcr-1,k,1)+scal(ii,jcr,k,1))
     &                    + avd * (scal(ii,jcr-1,k,1)-scal(ii,jcr,k,1)) 
          ww0(j) = ww0(j) + v00 * (scal(ii,jcr-1,k,2)+scal(ii,jcr,k,2))
     &                    + avd * (scal(ii,jcr-1,k,2)-scal(ii,jcr,k,2)) 
         enddo
         do i=1+i1zon(j,nb),i2zon(j,nb)
          ii = icl(i)
          jcr=jmaxl(ii,j)
          jcr1=jmaxl(icl(i-1),j)
          jsigne = sign(1,jcr1-jcr)
          jj1 = min(jcr,jcr1)
          jj2 = max(jcr,jcr1) - 1
          do jj=jj1,jj2
           ttm = tms(icl(i-1),jj,k)*tms(ii,jj,k) * DFLOAT(jsigne)
           v00 = 0.5*(u(ii,jj,k)+u(ii,jj+1,k)) + uiso(ii,jj,k)
           avd = min(1.0,abs(ccydif1*smx(ii,jj,1)*v00)+alphmi(k)) 
           v00 = ttm * cmy(ii,jj,1) * v00
           avd = avd * abs(v00) + ttm * ccydif2 * cmxy(ii,jj,1)  
!           avd = 0.0
           uu0(j) = uu0(j) + v00  
           vv0(j) = vv0(j) + v00 * (scal(icl(i-1),jj,k,1)+scal(ii,jj,k,1))
     &                     + avd * (scal(icl(i-1),jj,k,1)-scal(ii,jj,k,1)) 
           ww0(j) = ww0(j) + v00 * (scal(icl(i-1),jj,k,2)+scal(ii,jj,k,2))
     &                     + avd * (scal(icl(i-1),jj,k,2)-scal(ii,jj,k,2)) 
          enddo
         enddo
        enddo
        do j=1,jmtt
         uuu(j,k,nb)= uuu(j,k,nb) + uu0(j)*dz(k)
         vvv(j,nb)=   vvv(j,nb)   + vv0(j)*dz(k)
         vvv(j,nb+4)= vvv(j,nb+4) + ww0(j)*dz(k)   
        enddo
! end k-loop
       enddo
! end nb-loop
      enddo



! ------------------------ c nupdate = 1 ==> prepare and write output
      elseif (nupdate .eq.  1) then



! normalization factor usnjm(=1 for snapshot, =1/30 for monthly aver. etc.)
       if (usnjm.eq.0.d0) usnjm = 1.0d0
       dxsvrdrp = -dx*svrdrp*usnjm
       dxsvrdrpT=  0.5*dx*usnjm*rho0*cpo*1.d-15 
       dxsvrdrpS=  0.5*dx*svrdrp*usnjm 


! summing from bottom up
       do nb=1,nbsmax
        do k=2,kmax
         do j=1,jmtt
           uuu(j,k,nb)= uuu(j,k,nb) + uuu(j,k-1,nb)
         enddo
        enddo
       enddo


! add 3 basins to obtain global basin
       do nb=1,nbsmax
        do k=1,kmax
         do j=1,jmtt
          uuu(j,k,0)=uuu(j,k,0)+uuu(j,k,nb)
         enddo
        enddo        
       enddo
       do nb=1,nbsmax
        do j=1,jmtt
         vvv(j,0)=vvv(j,0)+vvv(j,nb)
         vvv(j,4)=vvv(j,4)+vvv(j,nb+4)
        enddo 
       enddo


! put proper scale factor: dxsvrdrp(T/S)
       do nb=0,nbsmax
        do k=1,kmax
         do j=1,jmtt
          uuu(j,k,nb)=uuu(j,k,nb)*dxsvrdrp
         enddo
        enddo
       enddo
       do nb=0,nbsmax
        do j=1,jmtt
         vvv(j,nb)  =vvv(j,nb)  *dxsvrdrpT
         vvv(j,nb+4)=vvv(j,nb+4)*dxsvrdrpS
        enddo
       enddo


! the meridional heat flux and salt flux
! ===> compute these via an anomaly: T_ref = 1.3 K and S_ref = 34.7 psu
       k=kmax
       do nb=0,nbsmax
        do j=1,jmtt
         vvv(j,nb)  =vvv(j,nb)  +uuu(j,k,nb)*rho0*cpo*1.d-9*(273.15+1.3)  
         vvv(j,nb+4)=vvv(j,nb+4)+uuu(j,k,nb)*34.7
        enddo
       enddo


! save for ouput in evolu
       do nb=0,nbsmax
        do k=2,kmax
         do j=1,jmtt
          vwx(j,k,nb)=uuu(j,k,nb)
         enddo
        enddo
       enddo


! remove southern ocean part for partial basins: for GrADS undef = -1.0e20
! j=15 corresponds to approximately 32S
       do nb=1,nbsmax
        do k=1,kmtt
         do j=1,16
          uuu(j,k,nb)= -1.0e20
         enddo
        enddo
       enddo
       do nb=1,nbsmax
        do j=1,16
         vvv(j,nb)  = -1.0e20
         vvv(j,nb+4)= -1.0e20
        enddo
       enddo


! write arrays to output files
       nrecstr0=nrecstr0+1           
       nrecstr1=nrecstr1+1  
       write(streamdat_id,rec=nrecstr0) uuu
       write(heatsaltfluxdat_id,rec=nrecstr1) vvv


! re-initialize arrays
!      do nb=0,nbsmax
!       do j=1,jmtt
!         do k=1,kmtt
!          uuu(j,k,nb) = 0.0
!         enddo
!       enddo
!      enddo
!      do nb=0,nbsmax+4
!       do j=1,jmtt
!        vvv(j,nb) = 0.0
!       enddo
!      enddo


      endif



      return
      end
