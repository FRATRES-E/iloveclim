!dmr -- Ajout du choix optionnel des composantes - Tue Dec 15 16:33:22 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Tue Dec 15 16:33:22 CET 2009

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!      ec_masq_mod : mise a jour du masque calotte dans ECBilt
!
!      Auteur : Didier M. Roche
!      Date   : 18 decembre 2008
!      Derniere modification : 16 novembre 2010, D.M. Roche, C. Dumas
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

module ec_masq_mod

  use global_constants_mod, only: dblp=>dp, ip

  implicit none
  private

  public :: ec_masq, update_masq
#if ( F_PALAEO == 1 )
  public :: load_topo_masq
#endif
#if ( F_PALAEO == 3 )
  public :: load_topo_masq_hr
#endif

contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
  subroutine ec_masq(init)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
!  Variables d'entree  : init = est on a l'initialisation ?
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

    use input_icemask, only: icemask, nlat, nlon

#if ( DOWNSCALING == 2 || ISM >= 2 )
    use output_ecbilt                        ! external — only: not available
#endif
#if ( DOWNSCALING == 2 )
    use input_subgrid2L, only: nbpointssg
#endif
#if ( ISM >= 2 )
    use input_flagsgris, only: nord_GRIS, sud_GRIS, masqueECB
#endif
#if ( F_PALAEO == 1 )
    use output_ecbilt,   only: masqueECB, where_update
#endif
#if ( NC_IMSK >= 1 )
    use ncio,            only: nc_read
#endif
#if ( NC_IMSK == 2 )
    use update_clio_bathy_tools, only: la_date
#endif

    implicit none

    logical, intent(in)                  :: init
    real(kind=4), dimension(nlat,nlon)   :: localmasq
#if ( NC_IMSK >= 1 )
    real(kind=4), dimension(nlon,nlat)   :: localmasqT
#endif
#if ( NC_IMSK == 2 )
    character(len=30) :: name_file
    character(len=5)  :: charI
#endif

    integer(kind=ip) :: i, j

    if (init) then
#if ( NC_IMSK == 1 )
      call nc_read("inputdata/icemask.nc", "icemask", localmasqT)
      localmasq = transpose(localmasqT)
#elif ( NC_IMSK == 2 )
      write(charI,'(I5.5)') la_date
      name_file = 'inputdata/icemask_'//trim(charI)//'k.nc'
      call nc_read(name_file, "icemask", localmasqT)
      localmasq = transpose(localmasqT)
#else
      open(251, convert='BIG_ENDIAN', file='inputdata/icemask.dat', &
           form='unformatted', access='direct', &
           recl=nlat*nlon*kind(localmasq(1,1)))
      read(251, rec=1) ((localmasq(i,j), j=1,64), i=1,32)
      close(251)
#endif

      icemask = localmasq

    else

#if ( DOWNSCALING == 2 )
      ! afq -- we are updating icemask only over the subgrid
      where ((nbpointssg(:,:) .gt. 0) .and. (masqueECB(:,:) .gt. 0.3))
        icemask(:,:) = 1.0
      elsewhere ((nbpointssg(:,:) .gt. 0) .and. (masqueECB(:,:) .le. 0.3))
        icemask(:,:) = 0.0
      endwhere
#endif

#if ( ISM >= 2 )
      if (nord_GRIS .gt. 1) then
        where (masqueECB(16:32,:) .gt. 0.3)
          icemask(16:32,:) = 1.0
        elsewhere
          icemask(16:32,:) = 0.0
        endwhere
      endif

      if (sud_GRIS .gt. 1) then
        where (masqueECB(1:15,:) .gt. 0.3)
          icemask(1:15,:) = 1.0
        elsewhere
          icemask(1:15,:) = 0.0
        endwhere
      endif

      write(*,*) "Icemask updated"
#endif

#if ( F_PALAEO == 1 )
      where (where_update .eq. 1)
        icemask = masqueECB
      endwhere
#endif

    endif

    call update_masq

  end subroutine ec_masq


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
  subroutine update_masq
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

    use comatm,          only: nlat, nlon
    use comphys,         only: irn
    use comemic_mod,     only: fracto
    use comrunlabel_mod, only: irunlabelf
    use input_icemask,   only: icemask, icemaskalb

    implicit none

    integer(kind=ip)             :: i, j
    real(kind=4)                 :: outdata(nlon,nlat)

    ! dmr lgmv9.1 : le masque glaciaire ne doit pas etre si c est une case ocean
    icemaskalb(:,:) = icemask(:,:)

    do j = 1, nlon
      do i = 1, nlat
        if (fracto(i,j) .gt. 0.99) icemask(i,j) = 0.0
        if ((icemask(i,j) .gt. 0.9) .and. (i .gt. 6)) then
          if ((j .gt. 10 .and. j .lt. 24) .and. &
              (i .gt. 18 .and. i .lt. 25)) then
            ! do nothing — Thibetan Plateau, might cause trouble
          else
            irn(i,j,2) = 27  ! Greenland type
          endif
        else if ((icemask(i,j) .gt. 0.1) .and. (i .le. 6)) then
          irn(i,j,2) = 23   ! Antarctica type
        endif
      enddo
    enddo

    open(777, file='outputdata/atmos/icemasku.ctl')
    write(777,fmt="('dset   ^icemasku.dat')")
    write(777,fmt="('options big_endian')")
    write(777,fmt="('undef ',1p,e12.4)") -1.0e20
    write(777,fmt="('title ECBILT orography')")
    write(777,fmt="('xdef ',i3,' linear ',2f7.3)") 64, 0.00, 5.625
    write(777,fmt="('ydef ',i3,' levels')")        32
    write(777,fmt="(' -85.7606 -80.2688 -74.7445 -69.2130 -63.6786')")
    write(777,fmt="(' -58.1430 -52.6065 -47.0696 -41.5325 -35.9951')")
    write(777,fmt="(' -30.4576 -24.9199 -19.3822 -13.8445 -8.30670')")
    write(777,fmt="(' -2.76890 2.76890 8.30670 13.8445 19.3822')")
    write(777,fmt="('  24.9199 30.4576 35.9951 41.5325 47.0696')")
    write(777,fmt="('  52.6065 58.1430 63.6786 69.2130 74.7445')")
    write(777,fmt="('  80.2688 85.7606')")
    write(777,fmt="('zdef ',i3,' linear ',2f7.2)") 1, 0.00, 1.00
    write(777,fmt="('tdef ',i4,' linear 1jan0001  1YR')") 1
    write(777,fmt="('vars 1')")
    write(777,fmt="('masq       1  99 orography ECBILT')")
    write(777,fmt="('endvars')")
    close(777)

    do i = 1, nlon
      do j = 1, nlat
        outdata(i,j) = icemask(j,i)
      enddo
    enddo

    open(777, convert='BIG_ENDIAN', &
         file='outputdata/atmos/icemasku.dat', &
         form='unformatted', access='direct', &
         recl=size(outdata)*kind(outdata(1,1)))
    write(777, rec=1) outdata
    close(777)

  end subroutine update_masq


#if ( F_PALAEO == 1 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
  subroutine load_topo_masq(init)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

    use output_ecbilt,    only: masqueECB, masq_trans, where_update, &
                                topoECB, topo_trans, topoECB_ti, time_max
    use comatm,           only: nlat, nlon
    use comdyn,           only: rmount
    use palaeo_timer_mod, only: indx_masq, indx_topo

    implicit none

    logical, intent(in)  :: init

    integer(kind=ip) :: indx, ios, t, i, j
    real(kind=4)     :: dummy
    real(kind=4)     :: outdata(nlon,nlat)

    integer(kind=ip), parameter :: lim_sud_sico_ECB = 20
    integer(kind=ip), parameter :: update_time  = 250
    integer(kind=ip), parameter :: indx_first   = 2601
    integer(kind=ip), parameter :: irunlabel_ref = 3000

    character(len=22) :: filin

    if (init) then

      ! Read masks
      filin = "masks-cycles.txt"
      open(2013, file=filin, iostat=ios)
      do t = 1, time_max
        read(2013,'(F8.0)') dummy
        do i = nlat, lim_sud_sico_ECB, -1
          read(2013,'(64I1)') (masq_trans(i,j,t), j=1,nlon)
        enddo
      enddo
      close(2013)

      ! Read orographies
      filin = "orography-cycles.txt"
      open(2013, file=filin, iostat=ios)
      do t = 1, time_max
        read(2013,'(F8.0)') dummy
        do i = nlat, lim_sud_sico_ECB, -1
          read(2013,'(64F8.0)') (topo_trans(i,j,t), j=1,nlon)
        enddo
      enddo
      close(2013)

      ! Build where array
      where_update(lim_sud_sico_ECB:nlat,:) = 1
      where_update(1:lim_sud_sico_ECB,:)    = 0

      ! Save initial topography
      topoECB_ti = rmount

    endif

    ! dmr --- retrieve indx from the global palaeo timer
    indx = indx_masq

    ! Update topo / mask fields
    where ((masq_trans(:,:,indx) .eq. 0) .and. &
           (where_update(:,:) .eq. 1))
      masqueECB(:,:) = 1.0_dblp
    elsewhere
      masqueECB(:,:) = 0.0_dblp
    endwhere

    ! fl changed indx to indx_topo
    topoECB(:,:) = topoECB_ti(:,:) + topo_trans(:,:,indx_topo)

    open(777, file='outputdata/atmos/icemaski.ctl')
    write(777,fmt="('dset   ^icemaski.dat')")
    write(777,fmt="('options big_endian')")
    write(777,fmt="('undef ',1p,e12.4)") -1.0e20
    write(777,fmt="('title ECBILT orography')")
    write(777,fmt="('xdef ',i3,' linear ',2f7.3)") 64, 0.00, 5.625
    write(777,fmt="('ydef ',i3,' levels')")        32
    write(777,fmt="(' -85.7606 -80.2688 -74.7445 -69.2130 -63.6786')")
    write(777,fmt="(' -58.1430 -52.6065 -47.0696 -41.5325 -35.9951')")
    write(777,fmt="(' -30.4576 -24.9199 -19.3822 -13.8445 -8.30670')")
    write(777,fmt="(' -2.76890 2.76890 8.30670 13.8445 19.3822')")
    write(777,fmt="('  24.9199 30.4576 35.9951 41.5325 47.0696')")
    write(777,fmt="('  52.6065 58.1430 63.6786 69.2130 74.7445')")
    write(777,fmt="('  80.2688 85.7606')")
    write(777,fmt="('zdef ',i3,' linear ',2f7.2)") 1, 0.00, 1.00
    write(777,fmt="('tdef ',i4,' linear 1jan0001  1YR')") 1
    write(777,fmt="('vars 1')")
    write(777,fmt="('masq       1  99 orography ECBILT')")
    write(777,fmt="('endvars')")
    close(777)

    do i = 1, nlon
      do j = 1, nlat
        outdata(i,j) = masqueECB(j,i)
      enddo
    enddo

    open(777, convert='BIG_ENDIAN', &
         file='outputdata/atmos/icemaski.dat', &
         form='unformatted', access='direct', &
         recl=size(outdata)*kind(outdata(1,1)))
    write(777, rec=1) outdata
    close(777)

  end subroutine load_topo_masq
#endif


#if ( F_PALAEO == 3 )
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!
  subroutine load_topo_masq_hr(init)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8!

    use input_subgrid2L, only: forcedISM, forcedMSK
#if ( F_PALAEO_FWF == 1 )
    use input_subgrid2L, only: forcedTHI
#endif
    use output_ecbilt,   only: topoECB
    use comdyn,          only: rmount
    use ncio,            only: nc_read

    implicit none

    logical, intent(in) :: init

    character(len=26) :: filin

    if (init) then
      filin = "Gano_40k-0k_hemin40.nc"
      call nc_read(filin, "sur", forcedISM(:,:,:))
      call nc_read(filin, "msk", forcedMSK(:,:,:))
#if ( F_PALAEO_FWF == 1 )
      call nc_read(filin, "thi", forcedTHI(:,:,:))
#endif
      ! similarly to F_PALAEO==1, initialise topoECB
      topoECB = rmount
    endif

  end subroutine load_topo_masq_hr
#endif


end module ec_masq_mod
