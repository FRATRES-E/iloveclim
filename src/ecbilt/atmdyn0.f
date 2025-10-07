!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:31 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:31 CET 2009

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_iatmdyn
!-----------------------------------------------------------------------
! *** initialise parameters and operators and read initial state
!-----------------------------------------------------------------------

      USE comatm
      USE comdyn
      USE comphys
      use comsurf_mod, only: fractn, nld, epss ! afq, topo=0 over the oceans
#if ( F_PALAEO_FWF == 2 )
     >                       , thi_chge
#endif
      use newunit_mod, only: coef_dat_id, berg_dat_id, win_dat_id, sum_dat_id
      use global_constants_mod, only: dblp=>dp, ip
      use comemic_mod, only: fini, irunlabel
      use comcoup_mod
      use comunit
#if ( NC_BERG >= 1 )      
      use ncio, only: nc_read
#endif      
#if ( NC_BERG == 2 )
      use update_clio_bathy_tools, only: la_date
#endif

#if ( ISM == 2 )
! dmr FLAG AJOUT GRISLI
      USE input_flagsGRIS
      USE output_ECBilt
! dmr FLAG AJOUT GRISLI
#endif

      implicit none

#if ( ISM == 1 )
#include "ismecv.com"
#endif

      integer i,j,k1,k2,k,l,m,n,ifail,ii,jj,i1,j1,nn
      real*8  pigr4,dis,dif,rll,ininag(nlat,nlon),asum
      real*8  r1,a,b,c,d,e,sqn,rsqn
      real*8  rnorm,rh0,dd
      real*8  agg(nlat,nlon), agg1(nlat,nlon), agg2(nlat,nlon)
      real*8  fw(nsh2),fs(nsh2),fors(nsh2,nvl), fmu(nlat,2)
      real*8  forw(nsh2,nvl),wsx(nsh2),areafac
      real*8  spv
      REAL*4  outdata(nlon,nlat)

      integer(kind=ip):: topography_read_ctl_id, topography_read_dat_id
     &                   , inatdyn_dat_id

#if ( NC_BERG >= 1 )      
      real*8 aggT(nlon,nlat)
#endif
#if ( NC_BERG == 2 )
      character*30 name_file
      character(len=5) :: charI
#endif 

      read (coef_dat_id) nshm, ll

! *** real parameters


      pigr4=4.d0*pi
      rl1=1.0d0/rrdef1**2
      rl2=1.0d0/rrdef2**2
      relt1=max(0.0d0,rl1/(trel*pigr4))
      relt2=max(0.0d0,rl2/(trel*pigr4))
      dis=max(0.0d0,1.0d0/(tdis*pigr4))
      rll=dble(ll(nsh))
      dif=max(0.0d0,1.0d0/(tdif*pigr4*(rll*(rll+1))**idif))

! *** zonal derivative operator

      k2=0
      do m=0,nm
        k1=k2+1
        k2=k2+nshm(m)
!PB        do k=k1,k2
        do k=k1,253
        
          rm(k)=dble(m)
        enddo
      enddo

! *** laplace/helmholtz direct and inverse operators

      do j=0,5
        rinhel(1,j)=0.0d0
      enddo

      diss(1,1)=0.0d0
      diss(1,2)=0.0d0

      do k=2,nsh
        r1=dble(ll(k)*(ll(k)+1))
        a=-r1-3.0d0*rl1
        b=-r1-3.0d0*rl2
        c=-r1-rl1
        d=-r1-rl2
        e=a*d+b*c
        rinhel(k,0)=-r1
        rinhel(k,1)=-1.0d0/r1
        rinhel(k,2)= d/e
        rinhel(k,3)= b/e
        rinhel(k,4)=-c/e
        rinhel(k,5)= a/e
        diss(k,2)=dis*r1
        diss(k,1)=-dif*r1**idif
      enddo

      do j=0,5
        do k=1,nsh
          rinhel(k+nsh,j)=rinhel(k,j)
        enddo
      enddo

      do j=1,2
        do k=1,nsh
          diss(k+nsh,j)=diss(k,j)
        enddo
      enddo

! *** legendre associated FUNCTIONs and derivatives

      read (coef_dat_id) pp
      read (coef_dat_id) pd
      read (coef_dat_id) pw

! *** compensation for normalization in nag fft routines

      sqn=sqrt(dble(nlon))
      rsqn=1d0/sqn
      do k=1,nsh
        do i=1,nlat
          pp(i,k)=pp(i,k)*sqn
          pd(i,k)=pd(i,k)*sqn
          pw(i,k)=pw(i,k)*rsqn
        enddo
      enddo


! *** computation of the weight needed in OASIS
!     include 'weight.h'

! *** initialization of coefficients for fft


      do j=1,nlon
        do i=1,nlat
          ininag(i,j)=1.0d0
        enddo
      enddo

      ifail=0
      call c06fpf (nlat,nlon,ininag,'i',trigd,wgg,ifail)

      ifail=0
      call c06fqf (nlat,nlon,ininag,'i',trigi,wgg,ifail)

! *** orography and dissipation terms

! *** fmu(i,1): sin(phi(i))
! *** fmu(i,2): 1-sin**2(phi(i))

      rnorm=1.0d0/sqrt(3.0d0*nlon)
      do i=1,nlat
        fmu(i,1)=rnorm*pp(i,2)
        fmu(i,2)=1.d0-fmu(i,1)**2
      enddo

! *** height of orography in meters

#if ( NC_BERG == 1 )
      call nc_read("inputdata/berg.nc","agg1",aggT)
      agg1 = transpose(aggT)
#elif ( NC_BERG == 2 )
      write(*,*) 'date for berg.nc', la_date
      write(charI,'(I5.5)'), la_date
      name_file='inputdata/berg_'//trim(charI)//'k.nc'
      call nc_read(name_file,"agg1",aggT)
      agg1 = transpose(aggT)

#if ( F_PALAEO_FWF == 2 )
      call nc_read(name_file,"thi_chge",aggT)
      thi_chge = -1 * transpose(aggT) !flux is positive if ice loss
      !write(*,*) '!!!!!!!!!!!!!!!!!!!!!'
      !write(*,*) 'thi change atmdyn0', thi_chge
#endif

#else      
      read (berg_dat_id) agg1
#endif
            
! afq -- 0 over the oceans instead of the wavy spectral signal
      where (fractn(:,:,nld).LT.epss)
          agg1(:,:)=0.d0
      endwhere


      rh0=max(0.0d0,0.001d0/h0)

#if ( ISM == 2 )
! dmr FLAG AJOUT GRISLI
! dmr lignes copiees sur la version AGISM ci-dessous
! dmr topoECB est la topo GRISLI interpolee ECBilt
      if ((nord_GRIS.GE.1).OR.(sud_GRIS.GE.1)) then
      ! initialisation de la topo type icesheet
        DO i=LBOUND(topoECB,1),UBOUND(topoECB,1)
          DO j=LBOUND(topoECB,2),UBOUND(topoECB,2)
           if (masqueECB(i,j).GT.0.0) agg1(i,j)=topoECB(i,j)
          ENDDO
        ENDDO
      endif ! pas d initialisation

! dmr FLAG AJOUT GRISLI
#endif

#if ( ISM == 1 )
      do j=1,nlon
        do i=1,nlat
           spv=0.
          if ((flgisma.AND.(i.lt.15)).OR.(flgismg.AND.(i.gt.15))) then
          if (rmount_ism(i,j).GT.spv) agg1(i,j)=rmount_ism(i,j)
          endif
        enddo
      enddo
#endif

         do j=1,nlon
          do i=1,nlat
!          agg(i,j)=fmu(i,1)*agg1(i,j)*rh0
          agg(i,j) = agg1(i,j)*rh0
          rmount(i,j)=agg1(i,j)
          if (rmount(i,j).lt.0d0) rmount(i,j)=0d0
         enddo
        enddo


#if (1)
      outdata = 0.0

      open(newunit=topography_read_ctl_id,
     &     file='outputdata/atmos/topography-read.ctl')
      write(topography_read_ctl_id,fmt="('dset   ^topography-read.dat')")
      write(topography_read_ctl_id,fmt="('options big_endian')")
      write(topography_read_ctl_id,fmt="('undef ',1p,e12.4)") -1.0e20
      write(topography_read_ctl_id,fmt="('title ECBILT orography')")
      write(topography_read_ctl_id,
     &      fmt="('xdef ',i3,' linear ',2f7.2)") 64,0.00,5.625
      write(topography_read_ctl_id,fmt="('ydef ',i3,' levels')") 32
      write(topography_read_ctl_id,
     & fmt="(' -85.7606 -80.2688 -74.7445 -69.2130 -63.6786')")
      write(topography_read_ctl_id,
     & fmt="(' -58.1430 -52.6065 -47.0696 -41.5325 -35.9951')")
      write(topography_read_ctl_id,
     & fmt="(' -30.4576 -24.9199 -19.3822 -13.8445 -8.30670')")
      write(topography_read_ctl_id,
     & fmt="(' -2.76890 2.76890 8.30670 13.8445 19.3822')")
      write(topography_read_ctl_id,
     & fmt="('  24.9199 30.4576 35.9951 41.5325 47.0696')")
      write(topography_read_ctl_id,
     & fmt="('  52.6065 58.1430 63.6786 69.2130 74.7445')")
      write(topography_read_ctl_id,fmt="('  80.2688 85.7606')")
      write(topography_read_ctl_id,
     & fmt="('zdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
      write(topography_read_ctl_id,
     & fmt="('tdef ',i4,' linear 1jan0001  1YR')") 1
      write(topography_read_ctl_id,fmt="('vars 1')")
      write(topography_read_ctl_id,
     & fmt="('var1       1  99 orography ECBILT')")
      write(topography_read_ctl_id,fmt="('endvars')")
      close(topography_read_ctl_id)
      do i=1,nlon
        do j=1,nlat
          outdata(i,j)=rmount(j,i)
        enddo
      enddo
      open(newunit=topography_read_dat_id,CONVERT='BIG_ENDIAN',
     &      file='outputdata/atmos/topography-read.dat'
     &         ,form='unformatted',
     &         access='direct',recl=Size(outdata)*Kind(outdata(1,1)))
      write(topography_read_dat_id,REC=1) outdata
      close(topography_read_dat_id)

#endif



! *** surface dependent friction

      lgdiss=((addisl.gt.0.0).or.(addish.gt.0.0))

      call ec_ggtosp (agg,orog)
      call ec_ddl (orog,ws)
      call ec_sptogg (ws,dorodl,pp)
      call ec_sptogg (orog,dorodm,pd)

      if (lgdiss) then

#if ( NC_BERG == 1 )     
      call nc_read("inputdata/berg.nc","agg2",aggT)
      agg2 = transpose(aggT)
#elif ( NC_BERG == 2 )
      write(charI,'(I5.5)'), la_date
      name_file='inputdata/berg_'//trim(charI)//'k.nc'
      call nc_read(name_file,"agg2",aggT)
      agg2 = transpose(aggT)
#else      
        read (berg_dat_id) agg2
#endif

        do j=1,nlon
          do i=1,nlat
            agg(i,j)=1.0d0+addisl*agg2(i,j)+
     &                addish*(1.0d0-exp(-0.001d0*agg1(i,j)))
          enddo
        enddo

        call ec_ggtosp (agg,ws)
        call ec_ddl (ws,wsx)

        call ec_sptogg (ws,rdiss,pp)
        call ec_sptogg (wsx,ddisdx,pp)
        call ec_sptogg (ws,ddisdy,pd)

        dd=0.5d0*diss(2,2)
        do j=1,nlon
          do i=1,nlat
            ddisdx(i,j)=dd*ddisdx(i,j)/fmu(i,2)
            ddisdy(i,j)=dd*ddisdy(i,j)*fmu(i,2)
          enddo
        enddo

      endif

! *** forcing term

      do l=1,nvl
        do k=1,nsh2
          forw(k,l)=0.0d0
          fors(k,l)=0.0d0
        enddo
      enddo

      do j=1,nlon
        do i=1,nlat
          forcggw1(i,j)=0d0
          forcggs1(i,j)=0d0
        enddo
      enddo

      if (iartif.eq.1) then

! ***   forcing data from liu, for is winter forcing, fors is summer

        read(win_dat_id,910)   ((forw(k,l),k=1,nsh2),l=1,nvl)
        read(sum_dat_id,910)   ((fors(k,l),k=1,nsh2),l=1,nvl)

! ***   liu's forcing at 200mb  in grid point

        do k=1,nsh2
          fw(k)=forw(k,1)
          fs(k)=fors(k,1)
        enddo

        call ec_sptogg(fw,forcggw1,pp)
        call ec_sptogg(fs,forcggs1,pp)

      endif

! *** input initial qprime and for

      if (irunlabel.eq.0) then
        do k=1,nsh2
          do l=1,nvl
            qprime(k,l)=0d0
            psi(k,l)=0d0
            for(k,l)=0d0
          enddo
          do l=1,ntl
            psit(k,l)=0d0
          enddo
        enddo
        do j=1,nlon
          do i=1,nlat
            u800(i,j)=0d0
            u500(i,j)=0d0
            u200(i,j)=0d0
            uv10(i,j)=0d0
            uvw10(i,j)=0d0

!dmr @-@ iceb0
! JONO march 2004 init/reset surface windvector to drive icebergs
            utot10(i,j)=0d0
            vtot10(i,j)=0d0
! JONO end
!dmr @-@ iceb0

            do l=1,nvl
              utot(i,j,l)=0d0
              udivg(i,j,l)=0d0
              divg(i,j,l)=0d0
            enddo
          enddo
        enddo
      else
        open(newunit=inatdyn_dat_id,file='startdata/inatdyn'//fini//'.dat',
     *        form='unformatted')
        read(inatdyn_dat_id) qprime,for
        close(inatdyn_dat_id)
      endif

      write(6,*) 'ensemble',iens,numens
      if (iens.eq.1) then
       write(6,*) 'ensemble',iens,numens
       write(6,*) 1-log(float(numens))/23
       qprime=qprime*(1-log(float(numens))/23)
      endif

      if (ipert.ne.0) call ec_addperturb

  910 format(10e12.5)

      return
      end

!1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
!  (C) Copr. 1986-92 Numerical Recipes Software +.-).
      double precision function ran1(idum)

      implicit none
      integer idum,ia,im,iq,ir,ntab,ndiv
      real am,eps,rnmx
      parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,
     *ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer j,k,iv(ntab),iy
      save iv,iy
      data iv /ntab*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do j=ntab+8,1,-1
          k=idum/iq
          idum=ia*(idum-k*iq)-ir*k
          if (idum.lt.0) idum=idum+im
          if (j.le.ntab) iv(j)=idum
        enddo
        iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=min(am*iy,rnmx)

      return
      end function ran1

!1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
      SUBROUTINE ec_addperturb

#if ( COMATM == 1 )
      USE comatm
      USE comdyn
      USE comphys
      use comemic_mod, only:
#endif

#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#include "comphys.h"
#include "comemic.h"
#endif

      real*8 qpgg1(nlat,nlon),qpgg2(nlat,nlon),qpgg3(nlat,nlon)

      integer :: i,j
      double precision :: ran1

      write(*,*)'ipert=',ipert

      call ec_sptogg(qprime(1,1), qpgg1,pp)
      call ec_sptogg(qprime(1,2), qpgg2,pp)
      call ec_sptogg(qprime(1,3), qpgg3,pp)


      do i=1,nlat
        do j=1,nlon
          qpgg1(i,j)=qpgg1(i,j)*(1.025-0.05*ran1(ipert))
        enddo
      enddo

      do i=1,nlat
        do j=1,nlon
          qpgg2(i,j)=qpgg2(i,j)*(1.025-0.05*ran1(ipert))
        enddo
      enddo

      do i=1,nlat
        do j=1,nlon
          qpgg3(i,j)=qpgg3(i,j)*(1.025-0.05*ran1(ipert))
        enddo
      enddo

      call ec_ggtosp(qpgg1,qprime(1,1))
      call ec_ggtosp(qpgg2,qprime(1,2))
      call ec_ggtosp(qpgg3,qprime(1,3))

      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_ddt
!----------------------------------------------------------------------
! *** computation of time derivative of the potential vorticity fields
! *** input qprime, psi, psit
! *** output dqprdt
! *** NOTE psit is destroyed
!----------------------------------------------------------------------



#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer k,l
      real*8  dum1,dum2

! *** advection of potential vorticity at upper level

      call ec_jacob (psi(1,1),qprime(1,1),dqprdt(1,1))

! *** advection of potential vorticity at middle level

      call ec_jacob (psi(1,2),qprime(1,2),dqprdt(1,2))

! *** advection of potential vorticity and dissipation at lower level

      call ec_jacobd (psi(1,3),qprime(1,3),dqprdt(1,3))

! *** relaxation of temperature and forcing

      do k=1,nsh2
        dum1=relt1*psit(k,1)
        dum2=relt2*psit(k,2)
        dqprdt(k,1)=dqprdt(k,1)+dum1              +for(k,1)
        dqprdt(k,2)=dqprdt(k,2)-dum1+dum2         +for(k,2)
        dqprdt(k,3)=dqprdt(k,3)          -dum2    +for(k,3)
      enddo

! *** explicit horizontal diffusion

      do l=1,3
        do k=1,nsh2
          dqprdt(k,l)=dqprdt(k,l)+diss(k,1)*qprime(k,l)
        enddo
      enddo
      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_jacob (psiloc,pvor,sjacob)
!----------------------------------------------------------------------
! *** advection of potential vorticity
! *** input psiloc, pvor
! *** output sjacob
!----------------------------------------------------------------------


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif


      integer i,j,k
      real*8  psiloc(nsh2), pvor(nsh2), sjacob(nsh2),vv(nsh2)
      real*8  dpsidl(nlat,nlon), dpsidm(nlat,nlon), dvordl(nlat,nlon),
     *        dvordm(nlat,nlon), gjacob(nlat,nlon), dpsidls(nsh2)

! *** space derivatives of potential vorticity

      call ec_ddl (pvor,vv)
      call ec_sptogg (vv,dvordl,pp)
      call ec_sptogg (pvor,dvordm,pd)

! *** space derivatives of streamFUNCTION

      call ec_ddl (psiloc,dpsidls)
      call ec_sptogg (dpsidls,dpsidl,pp)
      call ec_sptogg (psiloc,dpsidm,pd)

! *** jacobian term

      do j=1,nlon
        do i=1,nlat
          gjacob(i,j)=dpsidm(i,j)*dvordl(i,j)-dpsidl(i,j)*dvordm(i,j)
        enddo
      enddo

      call ec_ggtosp (gjacob,sjacob)

! *** planetary vorticity advection

      do k=1,nsh2
        sjacob(k)=sjacob(k)-dpsidls(k)
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_jacobd (psiloc,pvor,sjacob)
!----------------------------------------------------------------------
! *** advection of potential vorticity and dissipation on gaussian grid
! *** input psiloc, pvor
! *** output sjacob
!----------------------------------------------------------------------


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer i,j,k
      real*8  psiloc(nsh2), pvor(nsh2), sjacob(nsh2)
      real*8  dpsidl(nlat,nlon), dpsidm(nlat,nlon), dvordl(nlat,nlon),
     *        dvordm(nlat,nlon), gjacob(nlat,nlon), vv(nsh2),
     *        azeta(nlat,nlon),dpsidls(nsh2)

! *** space derivatives of potential vorticity

      call ec_ddl (pvor,vv)
      call ec_sptogg (vv,dvordl,pp)
      call ec_sptogg (pvor,dvordm,pd)

! *** space derivatives of streamFUNCTION

      call ec_ddl (psiloc,dpsidls)
      call ec_sptogg (dpsidls,dpsidl,pp)
      call ec_sptogg (psiloc,dpsidm,pd)

! *** jacobian term + orographic forcing

      do j=1,nlon
        do i=1,nlat
          gjacob(i,j)=dpsidm(i,j)*(dvordl(i,j)+sinfi(i)*dorodl(i,j))-
     *                dpsidl(i,j)*(dvordm(i,j)+sinfi(i)*dorodm(i,j))
        enddo
      enddo

! *** dissipation


      if (lgdiss) then

! ***   spatially varying dissipation

        do k=1,nsh2
          vv(k)=diss(k,2)*psiloc(k)
        enddo

        call ec_sptogg (vv,azeta,pp)

        do j=1,nlon
          do i=1,nlat
            gekdis(i,j)=-dpsidm(i,j)*ddisdy(i,j)
     *              -dpsidl(i,j)*ddisdx(i,j)
     *              +rdiss(i,j)*azeta(i,j)
            gjacob(i,j)=gjacob(i,j) + gekdis(i,j)
          enddo
        enddo

        call ec_ggtosp (gjacob,sjacob)

      else

! ***   uniform dissipation

        call ec_ggtosp (gjacob,sjacob)

        do k=1,nsh2
          sjacob(k)=sjacob(k)+diss(k,2)*psi(k,3)
        enddo

      endif

! *** planetary vorticity advection

      do k=1,nsh2
        sjacob(k)=sjacob(k)-dpsidls(k)
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_jacobr (psiloc,pvor,sjacob)
!-----------------------------------------------------------------------
! *** computation of jacobian without planetary vorticity
! *** input psiloc, pvor
! *** output sjacob
!-----------------------------------------------------------------------


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer i,j,k
      real*8  psiloc(nsh2), pvor(nsh2), sjacob(nsh2),vv(nsh2)
      real*8  dpsidl(nlat,nlon), dpsidm(nlat,nlon), dvordl(nlat,nlon),
     *        dvordm(nlat,nlon), gjacob(nlat,nlon), dpsidls(nsh2)

! *** space derivatives of potential vorticity

      call ec_ddl (pvor,vv)
      call ec_sptogg (vv,dvordl,pp)
      call ec_sptogg (pvor,dvordm,pd)

! *** space derivatives of streamFUNCTION

      call ec_ddl (psiloc,dpsidls)
      call ec_sptogg (dpsidls,dpsidl,pp)
      call ec_sptogg (psiloc,dpsidm,pd)

! *** jacobian term

      do j=1,nlon
        do i=1,nlat
          gjacob(i,j)=dpsidm(i,j)*dvordl(i,j)-dpsidl(i,j)*dvordm(i,j)
        enddo
      enddo

      call ec_ggtosp (gjacob,sjacob)


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_omoro (psiloc,sjacob)
!-----------------------------------------------------------------------
! *** computation of jacobian without planetary vorticity
! *** input psiloc, pvor
! *** output sjacob
!-----------------------------------------------------------------------


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer i,j,k
      real*8  psiloc(nsh2), sjacob(nsh2),vv(nsh2)
      real*8  dpsidl(nlat,nlon), dpsidm(nlat,nlon),
     *        gjacob(nlat,nlon), dpsidls(nsh2)

! *** space derivatives of streamFUNCTION

      call ec_ddl (psiloc,dpsidls)
      call ec_sptogg (dpsidls,dpsidl,pp)
      call ec_sptogg (psiloc,dpsidm,pd)

! *** jacobian term

      do j=1,nlon
        do i=1,nlat
          gjacob(i,j)=dpsidm(i,j)*dorodl(i,j)-dpsidl(i,j)*dorodm(i,j)
          gjacob(i,j)=sinfi(i)*gjacob(i,j)
        enddo
      enddo

      call ec_ggtosp (gjacob,sjacob)


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_ddl (as,dadl)
!-----------------------------------------------------------------------
! *** zonal derivative in spectral space
! *** input spectral field as
! *** output spectral field dadl which is as differentiated wrt lambda
!-----------------------------------------------------------------------



#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer k
      real*8 as(nsh,2), dadl(nsh,2)

      do k=1,nsh
        dadl(k,1)=-rm(k)*as(k,2)
        dadl(k,2)= rm(k)*as(k,1)
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_sptogg (as,agg,pploc)

!-----------------------------------------------------------------------
! *** conversion from spectral coefficients to gaussian grid
! *** input  spectral field as, legendre polynomials pploc (pp or pd)
! ***        where pp are legendre polynomials and pd derivatives with
! ***        respect to sin(fi)
! *** output gaussian grid agg
!-----------------------------------------------------------------------



#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer i,ifail,j,k,k1,k2,m,mi,mr,nlon1
      real*8  as(nsh,2), agg(nlat,nlon), pploc(nlat,nsh)

! *** inverse legendre transform

      do j=1,nlon
        do i=1,nlat
          agg(i,j)=0.0d0
        enddo
      enddo

      nlon1=nlon+1
      k2=nshm(0)

      do k=1,k2
        do i=1,nlat
          agg(i,1)=agg(i,1)+as(k,1)*pploc(i,k)
        enddo
      enddo

      do m=1,nm
        mr=m+1
        mi=nlon1-m
        k1=k2+1
        k2=k2+nshm(m)
        do k=k1,k2
          do i=1,nlat
            agg(i,mr)=agg(i,mr)+as(k,1)*pploc(i,k)
          enddo
          do i=1,nlat
            agg(i,mi)=agg(i,mi)-as(k,2)*pploc(i,k)
          enddo
        enddo
      enddo

! *** inverse fourier transform

      ifail=0
      call c06fqf (nlat,nlon,agg,'r',trigi,wgg,ifail)

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_ggtosp (agg,as)
!-----------------------------------------------------------------------
! *** conversion from gaussian grid (agg) to spectral coefficients (as)
! *** input array agg is destroyed
! *** output as contains spectral coefficients
!-----------------------------------------------------------------------



#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer ir,ifail,j,k,k1,k2,m,mi,mr,nlon1,i
      real*8  as(nsh,2), agg(nlat,nlon)
!
! *** fourier transform
!
      ifail=0
      call c06fpf (nlat,nlon,agg,'r',trigd,wgg,ifail)
!
! *** legendre transform
!
      do ir=1,2
        do k=1,nsh
          as(k,ir)=0.0d0
        enddo
      enddo

      nlon1=nlon+1

      k2=nshm(0)

      do k=1,k2
        do i=1,nlat
          as(k,1)=as(k,1)+agg(i,1)*pw(i,k)
        enddo
      enddo

      do m=1,nm
        mr=m+1
        mi=nlon1-m
        k1=k2+1
        k2=k2+nshm(m)
        do k=k1,k2
          do i=1,nlat
            as(k,1)=as(k,1)+agg(i,mr)*pw(i,k)
            as(k,2)=as(k,2)+agg(i,mi)*pw(i,k)
          enddo
        enddo
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_rggtosp (agg,as)
!-----------------------------------------------------------------------
! *** conversion from gaussian grid (agg) to spectral coefficients (as)
! *** input array agg is conserved
! *** output as contains spectral coefficients
!-----------------------------------------------------------------------



#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer i,ifail,ir,j,k,k1,k2,m,mi,mr,nlon1
      real*8 as(nsh,2), agg(nlat,nlon)
      real*8 store(nlat,nlon)

      do j=1,nlon
        do i=1,nlat
          store(i,j)=agg(i,j)
        enddo
      enddo

! *** fourier transform

      ifail=0
      call c06fpf (nlat,nlon,store,'r',trigd,wgg,ifail)

! *** legendre transform

      do ir=1,2
        do k=1,nsh
          as(k,ir)=0.0d0
        enddo
      enddo

      nlon1=nlon+1

      k2=nshm(0)

      do k=1,k2
        do i=1,nlat
          as(k,1)=as(k,1)+store(i,1)*pw(i,k)
        enddo
      enddo

      do m=1,nm
        mr=m+1
        mi=nlon1-m
        k1=k2+1
        k2=k2+nshm(m)
        do k=k1,k2
          do i=1,nlat
            as(k,1)=as(k,1)+store(i,mr)*pw(i,k)
            as(k,2)=as(k,2)+store(i,mi)*pw(i,k)
          enddo
        enddo
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_qtopsi
!-----------------------------------------------------------------------
! *** computation of streamFUNCTION from potential vorticity
! *** input  qprime which is potential vorticity field
! *** output psi, the streamFUNCTION and psit, the layer thicknesses
!-----------------------------------------------------------------------



#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer k
      real*8  r3

      do k=1,nsh2
        ws(k)=qprime(k,1)+qprime(k,3)
        psi(k,1)=rinhel(k,1)*(ws(k)+qprime(k,2))
        psi(k,2)=ws(k)-2.d0*qprime(k,2)
        psi(k,3)=qprime(k,1)-qprime(k,3)
      enddo

      do k=1,nsh2
        psit(k,1)=rinhel(k,2)*psi(k,2)+rinhel(k,3)*psi(k,3)
        psit(k,2)=rinhel(k,4)*psi(k,2)+rinhel(k,5)*psi(k,3)
      enddo

      r3=1./3
      do k=1,nsh2
        psi(k,2)=r3*(psi(k,1)-psit(k,1)+psit(k,2))
        psi(k,1)=psi(k,2)+psit(k,1)
        psi(k,3)=psi(k,2)-psit(k,2)
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_psitoq
!-----------------------------------------------------------------------
! *** computation of potential vorticity from stream FUNCTION
! *** input psi streamFUNCTION
! *** output qprime, the potential vorticity and psit, the layer thick.
!-----------------------------------------------------------------------



#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none

#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif


      integer k

      do k=1,nsh2
        psit(k,1)=psi(k,1)-psi(k,2)
        psit(k,2)=psi(k,2)-psi(k,3)
        qprime(k,1)=rinhel(k,0)*psi(k,1)-rl1*psit(k,1)
        qprime(k,2)=rinhel(k,0)*psi(k,2)+rl1*psit(k,1)-rl2*psit(k,2)
        qprime(k,3)=rinhel(k,0)*psi(k,3)+rl2*psit(k,2)
      enddo
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_psiq(sfin,qout)
!-----------------------------------------------------------------------
! ***  computation of potential vorticity qout from stream FUNCTION sfin
!-----------------------------------------------------------------------



#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none

#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer k
      real*8  sfin(nsh2,nvl),qout(nsh2,nvl),tus(nsh2)

      do k=1,nsh2
        tus(k)=rl1*sfin(k,1)-rl1*sfin(k,2)
      enddo

      do k=1,nsh2
        qout(k,1)=rinhel(k,0)*sfin(k,1)-tus(k)
        qout(k,2)=rinhel(k,0)*sfin(k,2)+tus(k)
      enddo

      do k=1,nsh2
        tus(k)=rl2*sfin(k,2)-rl2*sfin(k,3)
      enddo

      do k=1,nsh2
        qout(k,2)=qout(k,2)-tus(k)
        qout(k,3)=rinhel(k,0)*sfin(k,3)+tus(k)
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_qpsi(qin,sfout)
!-----------------------------------------------------------------------
! *** computation of streamFUNCTION bb from potential vorticity qin
!-----------------------------------------------------------------------



#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none

#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      real*8  qin(nsh2,nvl),sfout(nsh2,nvl), tus(nsh2,ntl), r3
      integer k

      do k=1,nsh2
        ws(k)=qin(k,1)+qin(k,3)
        sfout(k,1)=rinhel(k,1)*(ws(k)+qin(k,2))
        sfout(k,2)=ws(k)-2.*qin(k,2)
        sfout(k,3)=qin(k,1)-qin(k,3)
      enddo

      do k=1,nsh2
        tus(k,1)=rinhel(k,2)*sfout(k,2)+rinhel(k,3)*sfout(k,3)
        tus(k,2)=rinhel(k,4)*sfout(k,2)+rinhel(k,5)*sfout(k,3)
      enddo

      r3=1./3
      do k=1,nsh2
        sfout(k,2)=r3*(sfout(k,1)-tus(k,1)+tus(k,2))
        sfout(k,1)=sfout(k,2)+tus(k,1)
        sfout(k,3)=sfout(k,2)-tus(k,2)
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_qpsit(qin,tus)
!-----------------------------------------------------------------------
! *** computation of streamFUNCTION bb from potential vorticity qin
!-----------------------------------------------------------------------



#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none

#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif


      real*8  qin(nsh2,nvl),tus(nsh2,ntl), r3,sfout(nsh2,nvl)
      integer k

      do k=1,nsh2
        ws(k)=qin(k,1)+qin(k,3)
        sfout(k,1)=rinhel(k,1)*(ws(k)+qin(k,2))
        sfout(k,2)=ws(k)-2.*qin(k,2)
        sfout(k,3)=qin(k,1)-qin(k,3)
      enddo

      do k=1,nsh2
        tus(k,1)=rinhel(k,2)*sfout(k,2)+rinhel(k,3)*sfout(k,3)
        tus(k,2)=rinhel(k,4)*sfout(k,2)+rinhel(k,5)*sfout(k,3)
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_fmtofs (y,z)
!-----------------------------------------------------------------------
! *** transforms franco's format to the french format for global fields
! *** input  y spectral coefficients in franco's format
! *** output z spectral coefficients in french format
!-----------------------------------------------------------------------


!script >>> Les declarations suivantes ont ete faite par un script
#if ( COMATM == 1 )
      USE comatm
#endif
!sript <<<
      implicit none

!script >>> Les declarations suivantes ont ete faite par un script
#if ( COMATM == 0 )
#include "comatm.h"
#endif
!sript <<<

      integer   m,n,k,indx,l
      real*8    y(nsh2,nvl),z(nsh2,nvl)

      do l=1,nvl
        k=1
        do m=0,nm
          do n=max(m,1),nm
            k=k+1
            if (m.eq.0) then
              indx=n**2
            else
              indx=n**2+2*m-1
            endif
            z(indx,l)=y(k,l)
            if (m.ne.0) z(indx+1,l)=y(k+nsh,l)
          enddo
        enddo
      enddo
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_fstofm (y,z,ntr)
!-----------------------------------------------------------------------
! *** transforms the french format to franco's format for global fields
! *** input  y spectral coef. in french format, ntr is truncation limit
! *** output z spectral coefficients in franco's format
!-----------------------------------------------------------------------


!script >>> Les declarations suivantes ont ete faite par un script
#if ( COMATM == 1 )
      USE comatm
#endif
!sript <<<
      implicit none

!script >>> Les declarations suivantes ont ete faite par un script
#if ( COMATM == 0 )
#include "comatm.h"
#endif
!sript <<<

      integer   m,n,k,indx,i,l,ntr
      real*8    y(nsh2,nvl),z(nsh2,nvl)

      do l=1,nvl
        do i=1,nsh2
          z(i,l)=0d0
        enddo
        k=1
        do m=0,nm
          do n=max(m,1),nm
            k=k+1
            if ((m.le.ntr).and.(n.le.ntr)) then
              if (m.eq.0) then
                indx=n**2
              else
                indx=n**2+2*m-1
              endif
              z(k,l)=y(indx,l)
              if (m.ne.0) z(k+nsh,l)=y(indx+1,l)
            endif
          enddo
        enddo
      enddo
      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_forward
!-----------------------------------------------------------------------
! *** performs a fourth order runge kutta time step at truncation nm
! *** with time step dt
! *** dqdt calculates the time derivative
! *** input  qprime at current time
! *** output qprime at current time plus dt
!-----------------------------------------------------------------------


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer  k,l,nvar
      real*8   dt2,dt6
      real*8   y(nsh2,nvl),dydt(nsh2,nvl),yt(nsh2,nvl)
      real*8   dyt(nsh2,nvl),dym(nsh2,nvl)

      nvar=(nm+2)*nm
      dt2=dtt*0.5d0
      dt6=dtt/6d0
      call ec_fmtofs(qprime,y)
      call ec_dqdt(y,dydt)
      do l=1,nvl
        do k=1,nvar
          yt(k,l)=y(k,l)+dt2*dydt(k,l)
        enddo
      enddo
      call ec_dqdt(yt,dyt)
      do l=1,nvl
        do k=1,nvar
          yt(k,l)=y(k,l)+dt2*dyt(k,l)
        enddo
      enddo
      call ec_dqdt(yt,dym)
      do l=1,nvl
        do k=1,nvar
          yt(k,l)=y(k,l)+dtt*dym(k,l)
          dym(k,l)=dyt(k,l)+dym(k,l)
        enddo
      enddo
      call ec_dqdt(yt,dyt)
      do l=1,nvl
        do k=1,nvar
          y(k,l)=y(k,l)+dt6*(dydt(k,l)+dyt(k,l)+2.*dym(k,l))
        enddo
      enddo
      call ec_fstofm(y,qprime,nm)
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_dqdt(y,dydt)
!-----------------------------------------------------------------------
! *** computation of time derivative of the potential vorticity field
! *** input  y potential vorticity in french format
! *** output dydt time derivative of y in french format
! *** values of qprime, psi and psit are changed
!-----------------------------------------------------------------------


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none

#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif


      real*8  y(nsh2,nvl),dydt(nsh2,nvl)

      call ec_fstofm(y,qprime,nm)
      call ec_qtopsi
      call ec_ddt
      call ec_fmtofs(dqprdt,dydt)
      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_psitogeo
!-----------------------------------------------------------------------
! *** computes geopotential in [m2/s2] from the streamFUNCTION
! *** by solving the linear balance equation:
! *** del phi = (1 - mu**2 ) d psi/dmu + mu del psi
! *** the global mean value is not determined and set to zero
! *** input:  psi
! *** output: grpsi1,grpsi2,grpsi3,geopg(nlat,nlon,nvl)
!-----------------------------------------------------------------------


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
      USE comphys
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#include "comphys.h"
#endif




      integer i,j,l
      real*8  dmu(nlat),cdim,tempfac(ntl)
      real*8  delpsis(nsh2),delpsig(nlat,nlon)
      real*8  dmupsig(nlat,nlon),delgeog(nlat,nlon)
      real*8  delgeos(nsh2),geos(nsh2)

#if ( ABEL == 1 )
      real*8  x1d(nsh2), sdim
      real*8  chi_scaled(nsh2)
!      real*8  chigtest(nlat,nlon,nvl)
      real*8  dchidls(nsh2)
      real*8  dchidlg(nlat,nlon)
#endif

      do i=1,nlat
        dmu(i)=1-sinfi(i)**2
      enddo

      cdim=(om**2)*(radius**2)
#if ( ABEL == 1 )
      sdim=(om)*(radius**2)

      do l=1,nvl
         do k=1,nsh2
            x1d(k)=(radius**2)*divs(k,l)*rinhel(k,1)
         enddo
      enddo

      do k=1,nsh2
         chi_scaled(k)=x1d(k)/sdim
      enddo

!      do l=1,nvl
!         call ec_sptogg(chi(1,l),chigtest(1,1,l),pp)
!         do j=1,nlon
!            do i=1,nlat
!               write(980,*) chigtest(i,j,l)
!            enddo
!         enddo
!      enddo
#endif

      do l=1,nvl

! *** solve linear balance equation

        call ec_lap(psi(1,l),delpsis)
        call ec_sptogg(delpsis,delpsig,pp)
        call ec_sptogg(psi(1,l),dmupsig,pd)
#if ( ABEL == 1 )
         call ec_ddl(chi_scaled,dchidls)
         call ec_sptogg(dchidls,dchidlg,pp)
#endif

        do j=1,nlon
          do i=1,nlat
            delgeog(i,j)=dmu(i)*dmupsig(i,j)+
     *                      sinfi(i)*delpsig(i,j)
#if ( ABEL == 1 )
     *              -1*dchidlg(i,j)
#endif
          enddo
        enddo
        call ec_ggtosp(delgeog,delgeos)
        call ec_lapinv(delgeos,geos)
        geos(1)=0.d0
        call ec_sptogg(geos,geopg(1,1,l),pp)
        do j=1,nlon
          do i=1,nlat
            geopg(i,j,l)=cdim*geopg(i,j,l)
          enddo
        enddo

      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_omega3
!-----------------------------------------------------------------------
! *** computes the vertical velocity at the two temperature levels
! *** and the surface in pa/sec
! *** input psi,psit
! *** output omegs is vertical velocity at three levels in spectral form
! ***        omegg at gaussian grid
!-----------------------------------------------------------------------

#if ( COMATM == 1 )
      USE comatm, only: fzero,om,pi,radius,tlevel,rgas,dp
      USE comdyn
      use comemic_mod, only:
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#include "comemic.h"
#endif

      integer i,j,k,l
      real*8  adoro(nsh2),adpsit(nsh2,ntl),dpsitdt(nsh2,ntl)
      real*8  facom(nvl),facoc
      real*8  facd1,facd2
      real*8  facekm
      real*8  omegsd(nsh2)

      facom(1)=(dp*om)/(rrdef1**2*fzero)
      facom(2)=(dp*om)/(rrdef2**2*fzero)
      facom(3)=(dp*om)

      call ec_jacobr(psi(1,1),psit(1,1),adpsit(1,1))
      call ec_jacobr(psi(1,2),psit(1,2),adpsit(1,2))
      call ec_omoro(psi(1,3),adoro)

      call ec_ddt

      call ec_qpsit(dqprdt,dpsitdt)

      facoc=trel*4.d0*pi


! *** ekman dissipation induced omega at ground level
! *** rein

      do j=1,nlon
        do i=1,nlat
!          if (i.gt.14.and.i.lt.19) then
!            gekdis(i,j)=gekdis(i,j)/abs(sinfi(14))
!          else
!            gekdis(i,j)=gekdis(i,j)/abs(sinfi(i))
!          endif
          gekdis(i,j)=gekdis(i,j)/fzero
        enddo
      enddo

      call ec_ggtosp(gekdis,omegsd)

! *** first contributions that are odd in the equator

      do k=1,nsh2

! ***   level1  (350 hpa)

        omegs(k,1)= dpsitdt(k,1) -
     *              adpsit(k,1) + psit(k,1)/facoc

! ***   level2  (650 hpa)

        omegs(k,2)= dpsitdt(k,2) -
     *              adpsit(k,2) + psit(k,2)/facoc

! ***   level 3 (surface)

!        omegs(k,3)= omegsd(k) + adoro(k)
        omegs(k,3)= 0d0

      enddo

      do l=1,nvl
        do k=1,nsh2
          omegs(k,l)=omegs(k,l)*facom(l)
        enddo
      enddo

      do l=1,nvl
        call ec_sptogg(omegs(1,l),omegg(1,1,l),pp)
      enddo

      do l=1,nvl
        do j=1,nlon
          do i=1,nlat/2
            omegg(i,j,l)=-omegg(i,j,l)
          enddo
        enddo
      enddo

      do l=1,nvl
        call ec_ggtosp(omegg(1,1,l),omegs(1,l))
      enddo


! *** second contributions that are even in the equator
! *** diabatic forcing and lift due to orography

      facd1=(rgas*(dp**2))/((rrdef1*radius*om*fzero)**2*tlevel(1))
      facd2=(rgas*(dp**2))/((rrdef2*radius*om*fzero)**2*tlevel(2))

      do k=1,nsh2
        omegs(k,1)=omegs(k,1) - dfor1(k)*facd1
        omegs(k,2)=omegs(k,2) - dfor2(k)*facd2
!        omegs(k,3)=omegs(k,3) + adoro(k)*facom(3)
      enddo

      do l=1,nvl
        omegs(1,l)=0.
        call ec_sptogg(omegs(1,l),omegg(1,1,l),pp)
      enddo

      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_diver
!-----------------------------------------------------------------------
! *** computes divergence from omega using conservation of mass
! *** input omegs is vertical velocity in pa/s
! *** output divs is divergence in spectral form
! ***        divg at gaussian grid
!-----------------------------------------------------------------------



#if ( COMATM == 1 )
      USE comatm, only: nsh2, dp
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif


      integer k,l

      do k=1,nsh2
        divs(k,1)=-omegs(k,1)/dp
        divs(k,2)=(omegs(k,1)-omegs(k,2))/dp
        divs(k,3)=(omegs(k,2)-omegs(k,3))/dp
      enddo

      do l=1,nvl
        call ec_sptogg(divs(1,l),divg(1,1,l),pp)
      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_divwin
!-----------------------------------------------------------------------
! *** computes divergent wind from the divergence
! *** input  divs is divergence at the three pressure levels
! *** output udivg zonal divergent wind at gaussian grid
! ***        vdivg meridional divergent wind at gaussian grid
!-----------------------------------------------------------------------


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer i,j,k,l
      real*8  r1,r2
      real*8  x(nsh2),xhelp(nsh2),dxdl(nlat,nlon),dxdm(nlat,nlon)

      r1=radius
      r2=radius**2

      do l=1,nvl

         do k=1,nsh2
           x(k)=r2*divs(k,l)*rinhel(k,1)
           chi(k,l)=x(k)
         enddo

         call ec_sptogg(x,chig(1,1,l),pp)

         call ec_ddl(x,xhelp)
         call ec_sptogg(xhelp,dxdl,pp)
         call ec_sptogg(x,dxdm,pd)

         do j=1,nlon
           do i=1,nlat
             udivg(i,j,l)=dxdl(i,j)/(r1*cosfi(i))
             vdivg(i,j,l)=dxdm(i,j)*cosfi(i)/r1
           enddo
         enddo

      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_geowin
!-----------------------------------------------------------------------
! *** computation of geostrophic winds at all levels
! *** input psi streamFUNCTION in spectral form
! *** output u200,v200 wind components at 200 hPa
! ***        u500,v500 wind components at 500 hPa
! ***        u800,v800 wind components at 800 hPa
!
!-----------------------------------------------------------------------


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif

      integer i,j,k,l
      real*8  facwin2
      real*8  dpsdl(nlat,nlon),dpsdm(nlat,nlon),psik(nsh2),vv(nsh2)

! *** space derivatives of streamFUNCTION

      facwin2=radius*om

      do l=1,nvl

        do k=1,nsh2
          psik(k)=psi(k,l)
        enddo

        call ec_ddl (psik,vv)
        call ec_sptogg (vv,dpsdl,pp)
        call ec_sptogg (psik,dpsdm,pd)

        if (l.eq.1) then
          do j=1,nlon
            do i=1,nlat
              u200(i,j)=-facwin2*dpsdm(i,j)*cosfi(i)
              v200(i,j)=+facwin2*dpsdl(i,j)/cosfi(i)
            enddo
          enddo
        endif

        if (l.eq.2) then
          do j=1,nlon
            do i=1,nlat
              u500(i,j)=-facwin2*dpsdm(i,j)*cosfi(i)
              v500(i,j)=+facwin2*dpsdl(i,j)/cosfi(i)
            enddo
          enddo
        endif

        if (l.eq.3) then
          do j=1,nlon
            do i=1,nlat
              u800(i,j)=-facwin2*dpsdm(i,j)*cosfi(i)
              v800(i,j)=+facwin2*dpsdl(i,j)/cosfi(i)
            enddo
          enddo
        endif

      enddo


      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_totwind
!-----------------------------------------------------------------------
! *** computation of total wind at all levels
!
!-----------------------------------------------------------------------


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer i,j,k,l

      do l=1,nvl

        if (l.eq.1) then
          do j=1,nlon
            do i=1,nlat
              utot(i,j,l)=u200(i,j) + udivg(i,j,l)
              vtot(i,j,l)=v200(i,j) + vdivg(i,j,l)
            enddo
          enddo
        endif

        if (l.eq.2) then
          do j=1,nlon
            do i=1,nlat
              utot(i,j,l)=u500(i,j) + udivg(i,j,l)
              vtot(i,j,l)=v500(i,j) + vdivg(i,j,l)
            enddo
          enddo
        endif

        if (l.eq.3) then
          do j=1,nlon
            do i=1,nlat
              utot(i,j,l)=u800(i,j) + udivg(i,j,l)
              vtot(i,j,l)=v800(i,j) + vdivg(i,j,l)
            enddo
          enddo
        endif

      enddo


      return
      end


!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_lap(xs,xsl)
!-----------------------------------------------------------------------
! *** computation of laplace operator in spectral domain
! *** input  xs  field in spectral form
! *** output xsl laplace of xs in spectral form
!-----------------------------------------------------------------------


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer k
      real*8  xs(nsh2),xsl(nsh2)

      do k=1,nsh2
        xsl(k)=xs(k)*rinhel(k,0)
      enddo

      return
      end



!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_lapinv(xsl,xs)
!-----------------------------------------------------------------------
! *** computation of laplace operator in spectral domain
! *** input  xsl field in spectral form
! *** output xs  inverse laplace of xs in spectral form
!-----------------------------------------------------------------------


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif

      implicit none

#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif



      integer k
      real*8  xs(nsh2),xsl(nsh2)

      do k=1,nsh2
        xs(k)=xsl(k)*rinhel(k,1)
      enddo

      return
      end

!23456789012345678901234567890123456789012345678901234567890123456789012
      SUBROUTINE ec_forcdaily
!-----------------------------------------------------------------------
! *** calculates artificial forcing as a FUNCTION of the time
! *** of the year
!-----------------------------------------------------------------------

#if ( COMATM == 1 )
      USE comatm
      USE comdyn
      USE comphys
      use comemic_mod, only: day
#endif

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#include "comphys.h"
#include "comemic.h"
#endif


      integer i,j

      do i=1,nlat
        do j=1,nlon
          forcgg1(i,j)=dabs(180.d0-day)*forcggw1(i,j)/180.d0 +
     *             forcggs1(i,j) - dabs(180.d0-day)*forcggs1(i,j)/180.d0
        enddo
      enddo

      return
      end


