!dmr -- Ajout du choix optionnel des composantes - Tue Dec 15 16:33:22 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Tue Dec 15 16:33:22 CET 2009

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Cette routine ec_masq sert mettre a jour le masque calotte
!       dans ECBilt
!       (dans l'environnement logiciel LUDUS)
!
!      Auteur : Didier M. Roche
!      Date   : 18 decembre 2008
!      Derniere modification : 16 novembre 2010, D.M. Roche, C. Dumas
!-----|--1--------2---------3---------4---------5---------6---------7-|

       SUBROUTINE ec_masq(init)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables d'entree  : init = est on a l'initialisation ?
!       Variables de sortie :
!-----|--1--------2---------3---------4---------5---------6---------7-|

       USE input_icemask
#if ( DOWNSCALING == 2 || ISM >= 2 )
       USE output_ecbilt
#endif
#if ( DOWNSCALING == 2 )
       use input_subgrid2L, only : nbpointssg
#endif
#if ( ISM >= 2 )
       USE input_flagsgris
#endif
#if ( F_PALAEO == 1 )
       USE output_ecbilt, ONLY: masqueECB, where_update
#endif
#if ( NC_IMSK >= 1 )      
      use ncio, only: nc_read
#endif
#if ( NC_IMSK == 2 )
      use update_clio_bathy_tools, only: la_date
#endif
       IMPLICIT NONE

       LOGICAL, INTENT(IN) :: init
       REAL(KIND=4), DIMENSION(32,64) :: localmasq
#if ( NC_IMSK >= 1 )       
       REAL(KIND=4), DIMENSION(64,32) :: localmasqT
#endif       
#if ( NC_IMSK == 2 )
      character*30 name_file
      character(len=5) :: charI
#endif

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Variables locales
!-----|--1--------2---------3---------4---------5---------6---------7-|

       INTEGER i,j

       if (init) then
#if ( NC_IMSK == 1 )
      call nc_read("inputdata/icemask.nc","icemask",localmasqT)
      localmasq = transpose(localmasqT)
#elif ( NC_IMSK == 2 )
      write(charI,'(I5.5)'), la_date
      name_file='inputdata/icemask_'//trim(charI)//'k.nc'
      call nc_read(name_file,"icemask",localmasqT)
      localmasq = transpose(localmasqT)
#else
       open(251,CONVERT='BIG_ENDIAN',file='inputdata/icemask.dat'
     &  ,form='unformatted'
     &  ,access='direct',
     &   recl=32*64*4)
       read(251, REC=1) ((localmasq(i,j),j=1,64),i=1,32)
       close(251)
#endif       

       icemask = localmasq

       else
#if ( DOWNSCALING == 2 )
      ! afq -- we are updating icemask only over the subgrid
       where((nbpointssg(:,:).GT.0).and.(masqueECB(:,:).GT.0.3))
         icemask(:,:)=1.0
       elsewhere((nbpointssg(:,:).GT.0).and.(masqueECB(:,:).LE.0.3))
         icemask(:,:)=0.0
       endwhere
#endif
#if ( ISM >= 2 )

        if (nord_GRIS.GT.1) then
          WHERE(masqueECB(16:32,:).GT.0.3)
            icemask(16:32,:) = 1.0
          ELSEWHERE
            icemask(16:32,:) = 0.0
          ENDWHERE
        endif

        if (sud_GRIS.GT.1) THEN
          WHERE(masqueECB(1:15,:).GT.0.3)
            icemask(1:15,:) = 1.0
          ELSEWHERE
            icemask(1:15,:) = 0.0
          ENDWHERE
        endif

        WRITE(*,*) "Icemask updated"
#endif
#if ( F_PALAEO == 1 )
       WHERE(where_update.EQ.1)
          icemask = masqueECB
       ENDWHERE
#endif
       endif

       CALL update_masq

       END SUBROUTINE ec_masq

! =================================

      SUBROUTINE update_masq

#if ( COMATM == 1 )
      USE comatm
      USE comdyn
      use comphys
      use comemic_mod, only: fracto
      use comrunlabel_mod, only: irunlabelf
      use comsurf_mod, only:
      use comunit
#endif

      USE input_icemask

      implicit none


#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#include "comphys.h"
#include "comemic.h"
#include "comrunlabel.h"
#include "comsurf.h"
#include "comunit.h"
#endif

      INTEGER i,j
      REAL(KIND=4) ::  outdata(64,32)

! afq -- icemask cannot be 1 over the ocean (model crashes), so
!        we keep the info, only for the albedo (used in ec_co2oc)
      icemaskalb(:,:)=icemask(:,:)

      do j=1,nlon
       do i=1,nlat

cdmr --- lgmv9.1 : le masque glaciaire ne doit pas etre si c est une case oean
         if (fracto(i,j).gt.0.99) icemask(i,j) = 0.0
cdmr --- lgmv9.1 : le masque glaciaire ne doit pas etre si c est une case oean
         if ((icemask(i,j).gt.0.9).and.(i.gt.6)) then
          IF ((j.gt.10.and.j.lt.24).and.(i.gt.18.and.i.lt.25)) then
           ! do nothing this is Thibetan Plateau, might cause trouble
          else
            irn(i,j,2) = 27 ! Greenland type
          endif
         else if ((icemask(i,j).gt.0.1).and.(i.le.6)) then
          irn(i,j,2) = 23 ! Antarctica type
        endif

       enddo
      enddo

      open(777,file='outputdata/atmos/icemasku.ctl')
      write(777,fmt="('dset   ^icemasku.dat')")
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
      write(777,fmt="('masq       1  99 orography ECBILT')")
      write(777,fmt="('endvars')")
      close(777)

      do i=1,nlon
        do j=1,nlat
          outdata(i,j)=icemask(j,i)
        enddo
      enddo
      open(777,CONVERT='BIG_ENDIAN',file='outputdata/atmos/icemasku.dat
     &         ',form='unformatted',
     &         access='direct',recl=Size(outdata)*Kind(outdata(1,1)))
      write(777,REC=1) outdata
      close(777)

       END SUBROUTINE update_masq


#if ( F_PALAEO == 1 )

!-----|--1--------2---------3---------4---------5---------6---------7-|
       SUBROUTINE load_topo_masq(init)
!-----|--1--------2---------3---------4---------5---------6---------7-|

       USE output_ecbilt, ONLY : masqueECB, masq_trans, where_update
     &              , topoECB, topo_trans, topoECB_ti, time_max


#if ( COMATM == 1 )
      USE comatm
      USE comdyn
#endif
! fl added indx_topo to disentangle the effect of ice sheets topo & albedo
      use palaeo_timer_mod, only: indx_masq, indx_topo

       IMPLICIT NONE

#if ( COMATM == 0 )
#include "comatm.h"
#include "comdyn.h"
#endif

! dmr deprecated replaced by indx_masq       INTEGER, INTENT(IN) :: timer
       LOGICAL, INTENT(IN) :: init

!-----|--1--------2---------3---------4---------5---------6---------7-|
!-----|--1--------2---------3---------4---------5---------6---------7-|

       INTEGER :: indx, ios, t, i, j
       REAL(KIND=4) :: dummy
       REAL(KIND=4) ::  outdata(64,32)


       INTEGER, PARAMETER :: lim_sud_sico_ECB = 20
       INTEGER, PARAMETER :: update_time = 250 ! in yrs
! Current files contains 801,000 years of data with a 250 yrs step.
!         files start at 800,000 kyrs BP and stops at -750 yrs BP
!
! Start at 150,000 yrs B.P. is thus: 3200.-(150,000./250.-1.) = 2601
       INTEGER, PARAMETER :: indx_first = 2601 ! in indicies ...
       INTEGER, PARAMETER :: irunlabel_ref = 3000

       CHARACTER*22 :: filin


       IF (init) then

!-----|--1--------2---------3---------4---------5---------6---------7-|
!        Read masks
!-----|--1--------2---------3---------4---------5---------6---------7-|
         filin="masks-cycles.txt"

         OPEN(2013,file=filin,iostat=ios)

         DO t=1, time_max
           READ(2013,'(F8.0)') dummy

           DO i=nlat,lim_sud_sico_ECB,-1
             READ(2013,'(64I1)') (masq_trans(i,j,t),j=1,nlon)
           ENDDO

         ENDDO
         CLOSE(2013)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!        Read orographies
!-----|--1--------2---------3---------4---------5---------6---------7-|
         filin="orography-cycles.txt"

         OPEN(2013,file=filin,iostat=ios)

         DO t=1, time_max
           READ(2013,'(F8.0)') dummy

           DO i=nlat,lim_sud_sico_ECB,-1
            READ(2013,'(64F8.0)') (topo_trans(i,j,t),j=1,nlon)
         ENDDO
       ENDDO

       CLOSE(2013)


!-----|--1--------2---------3---------4---------5---------6---------7-|
!        Build where array
!-----|--1--------2---------3---------4---------5---------6---------7-|

       where_update(lim_sud_sico_ECB:nlat,:) = 1
       where_update(1:lim_sud_sico_ECB,:) = 0

!-----|--1--------2---------3---------4---------5---------6---------7-|
!        Save initial topography
!-----|--1--------2---------3---------4---------5---------6---------7-|

         topoECB_ti = rmount
       ENDIF


!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Calculate the right indx ... (timer for topo/mask)
!-----|--1--------2---------3---------4---------5---------6---------7-|

!       indx = MIN(indx_first+(timer-irunlabel_ref)/250,time_max)
! dmr --- Better version: retreive the indx from the global palaeo timer
       indx = indx_masq

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Update the topo / mask fields
!-----|--1--------2---------3---------4---------5---------6---------7-|

       WHERE((masq_trans(:,:,indx).EQ.0).AND.(where_update(:,:).EQ.1))
         masqueECB(:,:) = 1.0d0
       ELSEWHERE
         masqueECB(:,:) = 0.0d0
       ENDWHERE

!       where_update(:,:) = where_trans(:,:,indx)

! fl changed indx to indx_topo
       topoECB(:,:) = topoECB_ti(:,:) + topo_trans(:,:,indx_topo)

      open(777,file='outputdata/atmos/icemaski.ctl')
      write(777,fmt="('dset   ^icemaski.dat')")
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
      write(777,fmt="('masq       1  99 orography ECBILT')")
      write(777,fmt="('endvars')")
      close(777)

      do i=1,nlon
        do j=1,nlat
          outdata(i,j)=masqueECB(j,i)
        enddo
      enddo
      open(777,CONVERT='BIG_ENDIAN',file='outputdata/atmos/icemaski.dat
     &         ',form='unformatted',
     &         access='direct',recl=Size(outdata)*Kind(outdata(1,1)))
      write(777,REC=1) outdata
      close(777)
       END SUBROUTINE load_topo_masq
#endif /* on F_PALAEO */


#if ( F_PALAEO == 3 )


!-----|--1--------2---------3---------4---------5---------6---------7-|
       SUBROUTINE load_topo_masq_hr(init)
!-----|--1--------2---------3---------4---------5---------6---------7-|

       use input_subgrid2L, only: forcedISM, forcedMSK
#if ( F_PALAEO_FWF == 1 )
       use input_subgrid2L, only: forcedTHI
#endif
       use output_ecbilt, only : topoECB
       use comdyn, only: rmount
       use ncio

       implicit none

       logical, intent(in) :: init


       CHARACTER*26 :: filin

       IF (init) then

!-----|--1--------2---------3---------4---------5---------6---------7-|
!        Read orography and ice mask
!-----|--1--------2---------3---------4---------5---------6---------7-|
         filin="Gano_40k-0k_hemin40.nc"
         call nc_read(filin,"sur",forcedISM(:,:,:))
         call nc_read(filin,"msk",forcedMSK(:,:,:))
#if ( F_PALAEO_FWF == 1 )
         call nc_read(filin,"thi",forcedTHI(:,:,:))
#endif
! afq -- similarly to what is done for PALAEO==1, we need to initialise topoECB         
         topoECB = rmount
       ENDIF

      END SUBROUTINE load_topo_masq_hr


#endif /* on F_PALAEO */

