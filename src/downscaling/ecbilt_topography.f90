!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:35 CET 2009

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!      MODULE: ecbilt_topography
!
!>     @author  Didier M. Roche (dmr)
!>     @author  Aurélien Quiquet (aurel)
!
!>     @brief This module [ecbilt_topography] is handling all operations related
!              to the topography (orography) for the ECBilt model.
!             Initial topography read is still done for now in atmdyn ...
!
!>     @date Creation date: January, 26th, 2016
!>     @date Last modification: SLastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

      module ecbilt_topography

       implicit none
       private

       public :: ec_topo

       integer, parameter, public :: nb_levls = 11
#if ( DOWNSTS == 1)
       double precision, dimension(nb_levls), public :: rmount_virt, qmount_virt
#endif


      contains


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!      ROUTINE: define_virt_levels
!
!>     @brief This function simply defined the each level to downscale on
!              , in meters
!
!      DESCRIPTION:
!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

      function define_virt_levels(vert_virt_levls) result(returnValue)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   By reference variables ...
!>    @param[in]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

       double precision, dimension(nb_levls), intent(out) :: vert_virt_levls

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  Through commons variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

       logical :: returnValue

       vert_virt_levls( 1) =   10.0d0
       vert_virt_levls( 2) =  250.0d0
       vert_virt_levls( 3) =  500.0d0
       vert_virt_levls( 4) =  750.0d0
       vert_virt_levls( 5) = 1000.0d0
       vert_virt_levls( 6) = 1250.0d0
       vert_virt_levls( 7) = 1500.0d0
       vert_virt_levls( 8) = 2000.0d0
       vert_virt_levls( 9) = 3000.0d0
       vert_virt_levls(10) = 4000.0d0
       vert_virt_levls(11) = 5000.0d0

       returnValue = .true.

      end function define_virt_levels

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!      ROUTINE: rmount_to_qmount
!
!>     @brief This function converts the rmounts virtual layers to qmount ones
!
!      DESCRIPTION:
!
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

      function rmount_to_qmount(rmount_virt_levls,qmount_virt_levls)              &
     & result(returnValue)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr   By reference variables ...
!>    @param[in]
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

       double precision, dimension(nb_levls), intent(in)  :: rmount_virt_levls
       double precision, dimension(nb_levls), intent(out) :: qmount_virt_levls

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  Through commons variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

       logical :: returnValue

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|
! dmr repeat procedure for the min and max of the cells
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|

       ! From what I read, qmount is just a 9 points smoothed version of rmount
       !  no reason to exist here?

       qmount_virt_levls(:) = rmount_virt_levls(:)

       returnValue = .true.

      end function rmount_to_qmount

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!




!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!      Cette routine ec_topo est une reprise du couplage ECBilt
!       -- AGISM (Ph. Huybrecht), dans le code ECBilt
!      Le but de cette routine est de mettre à jour la topo dans ECBilt
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche (pour la reprise GRISLI)
!      Date   : 19 Novembre 2008
!      Derniere modification : 26 mai 2009, 29 mai 2013
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

! ----------------------------------------------------------------------
! *** replace topography by ISM topography
! ----------------------------------------------------------------------

       subroutine ec_topo

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  :
!       Variables de sortie :
!-----|--1--------2---------3---------4---------5---------6---------7-|

#if ( COMATM == 1 )
      use comatm,         only: nlat, nlon, nsh2, nvl, nsh, nm, ntl, pi
      use comdyn
      use comphys
      use comsurf_mod, only: rmountn, fractn, nld, epss
!afq -- deprecated      use comunit,     only: iuo
#endif

#if ( WRAP_EVOL == 1 )
      use comemic_mod,    only:
#endif

      use newunit_mod, only: berg_dat_id
      use output_ECBilt,  only: topoECB, modiffractn, masqueECB

#if ( ISM == 2 || ISM == 3 || DOWNSCALING >= 1 )
      use input_subgrid2L, only: nbpointssg
#endif

#if ( F_PALAEO == 1 )
      use output_ECBilt,  only: where_update
#endif

![DELETE] #if ( DOWNSTS == 1)
![DELETE] ! dmr Beware in spite of its name topdifECB are absolute orographies, not differences.
![DELETE]       use output_ECBilt,  only: qmount_dif, topdifECB
![DELETE] #endif

       implicit none

#if ( COMATM == 0 )
         ! comdyn provides
         ! a massive amount of variables ...
#include "comdyn.h"
         ! comphys.h ptovides
         ! qmount(nlat,nlon)
#include "comphys.h"
         ! comsurf.h ptovides
         ! rmountn(nlat,nlon,ntyps)
#include "comsurf.h"
#endif

#if ( COMATM == 0 )
#include "comemic.h"
          ! comunit.h provides
          ! iuo
#include "comunit.h"
#endif

          ! comcoup.h provides ... nothing?
! #include "comcoup.h"


!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|


       integer :: i,j,i1,j1,ii,jj,k1,k2,k,l,m,n,ifail,nn

       double precision ::   pigr4,dis,dif,rll
       double precision ::   r1,a,b,c,d,e,sqn,rsqn
       double precision ::   rnorm,rh0,dd
       double precision ::  asum,spv,areafac

       double precision, dimension(nlat,nlon) :: ininag,agg, agg1, agg2
       double precision, dimension(nsh2) :: fw,fs,wsx
       double precision, dimension(nsh2,nvl) :: fors,forw
       double precision, dimension(nlat,2) :: fmu

       real*4  outdata(64,32)

       logical          :: success

!-pour changer topo (dyn+thermo), mettre le nouveau champ dans rmount_ism:
        spv=0.

! dmr version en syntaxe tableau du 26 mai 2009
        modiffractn = 0.0

#if ( F_PALAEO == 1 )
        where ((where_update.GT.0).AND.(fractn(:,:,nld).GT.epss))
          rmount = topoECB
          rmountn(:,:,nld) = topoECB
        end where
#endif

#if ( DOWNSCALING >= 1 )
        where (nbpointssg.GT.0)
          rmount = topoECB
          rmountn(:,:,nld) = topoECB
        end where
#endif

        where (rmountn(:,:,nld).LT.0.0)
           rmountn(:,:,nld) = 0.0
        end where

#if ( ISM == 2 || ISM == 3 || DOWNSCALING == 2 )
! dmr ci-dessous le cas ou l on a une case GRISLI (ou sous grille!) avec de la glace
! dmr mais pas de terre dans ECBilt = on rajoute un morceau tout ptit
! dmr de continent !!
        where ((masqueECB.GT.0.0).AND.(fractn(:,:,nld).LT.epss))
           fractn(:,:,nld) = epss + 1.0D-3
           modiffractn = 1.0
        end where
#endif
        qmount(:,:) = rmountn(:,:,nld)

        do j=1,nlon
        do i=1,nlat
          asum=0d0
          do i1=-1,1
            do j1=-1,1
              ii=i+i1
              jj=j+j1
              if (ii.lt.1) then
                ii=1
                jj=jj+nlon/2
              endif
              if (ii.gt.nlat) then
                ii=nlat
                jj=jj+nlon/2
              endif
              if (jj.lt.1) jj=jj+nlon
              if (jj.gt.nlon) jj=jj-nlon
                asum=asum+rmountn(ii,jj,nld)
              enddo
            enddo
            qmount(i,j)=asum/9d0
          enddo
        enddo

#if ( DOWNSTS == 1)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|
! dmr   If a sub-grid is available, I do not care which, I just proceed to
!        compute the adequate virtual levels!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8|

        success = define_virt_levels(rmount_virt)
        success = rmount_to_qmount(rmount_virt,qmount_virt)

! [WARNING]
! dmr --- BEWARE rmount_d and qmount_d are not computed anymore !!!
! dmr --- Next step: move the computations done in atmphys0.f on the vertical
!          levels to module atmphys_d also on the vertical ...
! [WARNING]

#endif

       open(2080,file='outputdata/globals/topographie',status='unknown')
         do i=1,nlat
           write(2080,'(1P64E10.3)') (rmountn(i,j,nld),j=1,nlon)
         enddo
         close(2080)

!        read (iuo+1) nshm, ll

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
        do k=k1,k2
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

       rewind (berg_dat_id)
       read (berg_dat_id) agg1
       rh0=max(0.0d0,0.001d0/h0)

! afq -- we
       where (fractn(:,:,nld).LT.epss)
         agg1(:,:) = 0.d0
       endwhere
         
       do j=1,nlon
         do i=1,nlat
#if ( F_PALAEO == 1 )
           if ((where_update(i,j).GT.0).AND.(fractn(i,j,nld).GT.epss)) then
             agg1(i,j)=topoECB(i,j)
           endif
#endif
#if ( DOWNSCALING >= 1 )
           if (nbpointssg(i,j).GT.0.0) then
             agg1(i,j)=topoECB(i,j)
           endif
#endif
         enddo
       enddo

         do j=1,nlon
          do i=1,nlat
!          agg(i,j)=fmu(i,1)*agg1(i,j)*rh0
          agg(i,j) = agg1(i,j)*rh0
          rmount(i,j)=agg1(i,j)
          if (rmount(i,j).lt.0d0) rmount(i,j)=0d0
         enddo
        enddo

! *** surface dependent friction

      lgdiss=((addisl.gt.0.0).or.(addish.gt.0.0))

      call ec_ggtosp (agg,orog)
      call ec_ddl (orog,ws)
      call ec_sptogg (ws,dorodl,pp)
      call ec_sptogg (orog,dorodm,pd)

      if (lgdiss) then

        read (berg_dat_id) agg2
        do j=1,nlon
          do i=1,nlat
            agg(i,j)=1.0d0+addisl*agg2(i,j)+                                      &
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

      endif ! ldiss

#if ( ROUTEAU == 1 ) /* Routing of surface fwf is dependent on topography, change it thus */
      CALL ec_inirunoff_2
#endif

      open(777,file='outputdata/atmos/orographu.ctl')
      write(777,fmt="('dset   ^orographu.dat')")
      write(777,fmt="('options big_endian')")
      write(777,fmt="('undef ',1p,e12.4)") -1.0e20
      write(777,fmt="('title ECBILT orography')")
      write(777,fmt="('xdef ',i3,' linear ',2f7.3)") 64,0.00,5.625
      write(777,fmt="('ydef ',i3,' levels')") 32
      write(777,fmt="(' -85.7606 -80.2688 -74.7445 -69.2130 -63.6786')")
      write(777,fmt="(' -58.1430 -52.6065 -47.0696 -41.5325 -35.9951')")
      write(777,fmt="(' -30.4576 -24.9199 -19.3822 -13.8445 -8.30670')")
      write(777,fmt="(' -2.76890 2.76890 8.30670 13.8445 19.3822')")
      write(777,fmt="('  24.9199 30.4576 35.9951 41.5325 47.0696')")
      write(777,fmt="('  52.6065 58.1430 63.6786 69.2130 74.7445')")
      write(777,fmt="('  80.2688 85.7606')")
      write(777,fmt="('zdef ',i3,' linear ',2f7.2)") 1,0.00,1.00
      write(777,fmt="('tdef ',i4,' linear 1jan0001  1YR')") 1
      write(777,fmt="('vars 1')")
      write(777,fmt="('orog       1  99 orography ECBILT')")
      write(777,fmt="('endvars')")
      close(777)
      do i=1,nlon
        do j=1,nlat
          outdata(i,j)=rmount(j,i)
        enddo
      enddo
      open(777,CONVERT='BIG_ENDIAN',file='outputdata/atmos/orographu.dat          &
     &',form='unformatted', access='direct',recl=Size(outdata)*Kind(outdata(1,1)))
      write(777,REC=1) outdata
      close(777)

      return

      end subroutine ec_topo

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!      End of the function here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

       end module ecbilt_topography

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
