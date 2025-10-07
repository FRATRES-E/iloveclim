!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:40 CET 2009
#include "choixcomposantes.h"
!dmr -- Ajout du choix optionnel des composantes - Thu Dec 17 11:57:40 CET 2009

      SUBROUTINE ec_inierror

      implicit none

      integer nwarns

      common /cerror/nwarns

      nwarns = 0

      return
      end

      SUBROUTINE ec_error(ierr)

      USE comatm
      USE comdyn
      USE comphys
      use comemic_mod, only:
      use comcoup_mod
      use comsurf_mod
      use comunit
      use global_constants_mod, only: dblp=>dp, ip
      use newunit_mod, only: error_id
      
#if ( ISOATM >= 1 )
       USE iso_param_mod, ONLY : ieau
#endif

      implicit none


      integer ierr,nwarns,i,j
      integer(kind=ip):: state_grads_id
      common /cerror/nwarns

      if (ierr.lt.100) then

      if (ierr.eq.1) then
        write(error_id,*) 'error in routine surftl of landmodel'
        write(error_id,*) 'too many iterations (> 20)'
      endif

      if (ierr.eq.2) then
        write(error_id,*) 'error in routine surftl of landmodel'
        write(error_id,*) 'divergence of solution'
      endif

      if (ierr.eq.3) then
        write(error_id,*) 'error in routine test of atmdiag0.f'
        write(error_id,*) 'surface temperature out of range'
      endif

      if (ierr.eq.4) then
        write(error_id,*) 'error in routine iniglobal of initial0.f'
        write(error_id,*) 'error in reading gausspoints'
      endif

      if (ierr.eq.5) then
        write(error_id,*) 'error in routine rooster of oceandyn0.f'
        write(error_id,*) 'error in reading landseamask'
      endif

      if (ierr.eq.6) then
        write(error_id,*) 'error in routine seaice of icemodel0.f'
        write(error_id,*) 'tijs undefined'
      endif

      if (ierr.eq.7) then
        write(error_id,*) 'error in routine surftl of landmodel0.f'
        write(error_id,*) 'too many iterations (> 100) in zbrac'
      endif

      if (ierr.eq.8) then
        write(error_id,*) 'error in routine surftl of landmodel0.f'
        write(error_id,*) 'too many iterations (> 100) in zbrent'
      endif

      if (ierr.eq.9) then
        write(error_id,*) 'error in routine surfacetemp of icemodel0.f'
        write(error_id,*) 'too many iterations (> 100)'
      endif

      if (ierr.eq.10) then
        write(error_id,*) 'error in routine roostl of lakemodel0.f'
        write(error_id,*) 'error in reading lakemask'
      endif

      if (ierr.eq.11) then
        write(error_id,*) 'error in routine surfacetemp of icemodel0.f'
        write(error_id,*) 'too many iterations (> 100) in zbrac'
      endif

      if (ierr.eq.12) then
        write(error_id,*) 'error in routine surfacetemp of icemodel0.f'
        write(error_id,*) 'too many iterations (> 100) in zbrent'
      endif

      if (ierr.eq.13) then
        write(error_id,*) 'error in routine zbrent of root0.f'
        write(error_id,*) 'root not in specified interval'
      endif

      if (ierr.eq.15) then
        write(error_id,*) 'lake salinity out of range (<15)'
      endif

      if (ierr.eq.16) then
        write(error_id,*) 'error in routine iniland of landmodel0.f'
        write(error_id,*) 'land point not in a specified landbasin'
      endif

      if (ierr.eq.17) then
        write(error_id,*) 'error in reading ice mask in oceanfixed0.f'
      endif

      if (ierr.eq.18) then
        write(error_id,*) 'unrealistic ground pressure'
      endif

      if (ierr.eq.19) then
        write(error_id,*) 'too many ocean-lake neighbours in inioclacp'
      endif

      if (ierr.eq.20) then
        write(error_id,*) 'failure in expint called by detqmax'
      endif

      if (ierr.eq.99) then
        write(error_id,*) 'not a number detected in surfacetemperature'
      endif

      else

        if (ierr.eq.117) then
          write(error_id,*) ' longwave par. out of range too low temp'
        endif

        if (ierr.eq.118) then
          write(error_id,*) ' longwave par. out of range too high temp'
        endif

        if (ierr.eq.120) then
          write(error_id,*) 'convec: 10 iterations, still unstable'
        endif

        if (ierr.eq.121) then
          write(error_id,*) 'qmax out of range'
        endif

        if (ierr.eq.122) then
          write(error_id,*) 'rain larger than maximum set by rainmax'
        endif

        if (ierr.eq.123) then
          write(error_id,*) 'failure in netcdf export'
        endif

        if (ierr.eq.999) then
          write(error_id,*) 'integration succesfully completed'
        endif

        call flush(error_id)
        nwarns=nwarns+1
        if (nwarns.gt.100000) goto 10

        return

      endif

 10   continue

      write(error_id,*) 'program aborted due to soft error'
      write(error_id,*) 'see state.grads for state of the system'

      open(newunit=state_grads_id,file='state.grads',form='unformatted')

      write(state_grads_id) ((real(u200(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(v200(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(u500(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(v500(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(u800(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(v800(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(udivg(i,j,1)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(vdivg(i,j,1)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(udivg(i,j,2)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(vdivg(i,j,2)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(udivg(i,j,3)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(vdivg(i,j,3)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(temp0g(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(temp2g(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(temp4g(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(tempsg(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(tsurf(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(torain(i,j,iwater)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(tosnow(i,j,iwater)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(rmoisg(i,j,iwater)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(cormois(i,j,iwater)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(adsnow(i,j,iwater)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(ahic(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(eflux(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(hflux(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(hesws(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(dlrads(i,j)),j=1,nlon),i=1,nlat)
      write(state_grads_id) ((real(ulrads(i,j)),j=1,nlon),i=1,nlat)
      close(state_grads_id)
      stop

      return
      end
