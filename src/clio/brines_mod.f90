#include "choixcomposantes.h"

      module brines_mod

! use section

      use para0_mod, only: imax, jmax
      use global_constants_mod, only: str_len, dblp=>dp, ip
      use newunit_clio_mod, only: mouchard_id

      implicit none

      real, dimension(imax, jmax) :: frac_2D
      real :: frac_1D
      integer la_date_brines
      integer(ip),parameter, public ::  frac_nb = 43
      real(dblp), dimension(frac_nb), public :: frac_date
      real(dblp), dimension(frac_nb), public :: frac_value

      contains

!-----------------------------------------------------------------------
! Pour la plongee des brines
!-----------------------------------------------------------------------
      subroutine brines

      ! START_OF_USE_SECTION
      use global_constants_mod, only: dblp=>dp, ip

      use const_mod

      use para0_mod
      use para_mod
      use bloc0_mod
      use bloc_mod
      use thermo_mod, only:
      use ice_mod, only:
      use update_clio_bathy_tools, only: salglobal


!cfc  use trace_mod
!! END_OF_USE_SECTION

!---
! Ccpl [Ccp0] => ligne specifique a la version avec [sans] couplage .
!---

!cfc  include 'trace.com'
!fl modified nb_tracer to include PaTh

      integer(kind=ip):: l, j, i, li, k
      real(kind=dblp) :: vcell, dif


      real frac
      real*8 totsalts, totsaltf
      real*8 totsalt1, totsalt2, totsalt3
      real*8 vol, vol_global
      integer, parameter :: nb_tracer=16 !12 without PaTh, 16 with
      integer indices(nb_tracer)
      real*8 scal_surface(nb_tracer)
      real*8 sum_brines(nb_tracer),somme_test(nb_tracer)
      real*8 coef
      real*8 vfonte
!      real*8 vbrines pour la temperature
      real*8 fluxbrines(imax,jmax)
      common / common_brines / fluxbrines
      real*8 salinite_ref
      integer nmax

      nmax=nsmax-1 ! all tracers except temperature
!      nmax=12
!      nmax=1
      salinite_ref=salglobal
!      salinite_ref=34.7031774516157 !!!!
!      salinite_ref=34.703177 !!!!
!      salinite_ref=35.703177 !!!! at LGM
!      print*, salinite_ref
!cccc Parametre a modifier ccccccccccccccccccccccc
#if ( BRINES == 1 )
      frac=0.8 !brines parameter : between 0 and 1
#elif ( BRINES == 4 )
      call get_brines_frac
      frac=frac_1D
#endif
!ccccccccccccccccccccccccccccccccccccccccccccccccc

      totsalts=0.0
      totsaltf=0.0
      vfonte=0.0

      indices(1)=2 ! salinite
      indices(2)=3 !=> ODOC
      indices(3)=4 !=> ODOCS
      indices(4)=5 !=> ODIC
      indices(5)=6 !=> OPO4
      indices(6)=7 !=> ONO3
      indices(7)=8 !=> OALK
      indices(8)=9 !=> OO2
      indices(9)=10 !=> OC13
      indices(10)=11 !=> ODOC13
      indices(11)=12 !=> ODOCS13
      indices(12)=13 !=> OC14
      indices(13)=14 !fl padiss
      indices(14)=15 !fl papart
      indices(15)=16 !fl thdiss

!      temperature...

      do l=1,nb_tracer
        sum_brines(l)=0.0
        somme_test(l)=0.0
      enddo

      totsalt1=0.0
      vol=0.0
      do i=1,imax
        do j=1,jmax
          do k=1,kmax
          vol=vol+aire(i,j)*dz(k)*tms(i,j,k)
          totsalt1=totsalt1+scal(i,j,k,2)*aire(i,j)*dz(k)*tms(i,j,k)
          enddo
        enddo
      enddo
      totsalt1=totsalt1/vol
!      write(mouchard_id,*) 'salinity global 1', totsalt1
      salinite_ref=totsalt1



!nb phiss(i,j,2) flux de sel rejete: convention positif vers le haut
!c si dans le sud
      do j=js1,js2
        if (j.lt.32) then !dans le sud
          do i=is1(j),is2(j)
            if (fluxbrines(i,j).gt.0.0) then ! si flux de sel rejete              
!            write(mouchard_id,*) 'fluxbrines brines',fluxbrines(i,j)
#if ( BRINES == 3 )
             frac=frac_2D(i,j)
#endif
             !print*, 'frac in brines_mod ', frac
             !cefficient d'enrichissement
             coef=fluxbrines(i,j)/scal(i,j,ks2,2)
!l              do l=1,1
              do l=1,nmax
               li=indices(l)
!               print*, l, li
               ! on garde en memoire les valeures de surface
               scal_surface(l)=scal(i,j,ks2,li)

               !en surface on enleve du sel
               k=ks2
               ! attention cas special sel pas d enrichissement (deja
               ! inclus)
               if (li.eq.2) then
                scal(i,j,k,li)=scal(i,j,k,li)                           &
!     >              -frac*fluxbrines(i,j)*tms(i,j,k)
                   -frac*coef*scal_surface(l)*tms(i,j,k)
!              totsalts=totsalts-frac*fluxbrines(i,j)*tms(i,j,k)
               totsalts=totsalts-frac*coef*scal_surface(l)*tms(i,j,k)   &
                    *dz(k)*aire(i,j)
               else ! pour les autres l]
                vcell=dz(k)*aire(i,j)
                call surface_brines(scal(i,j,k,li)                      &
                      ,coef,frac,sum_brines(l),tms(i,j,k),vcell)
               endif

               ! au fond on ajoute du sel
               k=kfs(i,j)
!               print*,'fond1', li, scal(i,j,k,li)
               call fond_brines(scal(i,j,k,li)                          &
                 ,coef,frac,scal_surface(l),dz(k),dz(ks2),tms(i,j,k))
!               print*,'fond2', li, scal(i,j,k,li)
!              totsaltf=totsaltf+frac*fluxbrines(i,j)*dz(ks2)/dz(k)
               totsaltf=totsaltf+frac*coef*scal_surface(l)              &
                      *dz(ks2)/dz(k)*tms(i,j,k)*dz(k)*aire(i,j)
               ! pour la temperature
               !vbrines=frac * FSIO(i,n) * TSTOC * SQRO2(i,n)           &
               !          / ( sb-SM(i,j,n) )
               !TM(i,j,n)=(TM(i,j,n)*Vcell+TM_ori*vbrines)/Vcell
              enddo

            elseif(fluxbrines(i,j).lt.0.0) then !si fonte de glace
              !calcul du volume de la case
              vcell=dz(ks2)*aire(i,j)
              vfonte=vfonte+vcell 
            endif
          enddo
        endif
      enddo

! correction pour les variables sauf sel et temperature : rejet lors de
! la fonte

      do j=js1,js2
        if (j.lt.32) then !dans le sud
          do i=is1(j),is2(j)
            if (fluxbrines(i,j).lt.0.0) then !si fonte de glace
              vcell=dz(ks2)*aire(i,j)
!l              do l=1,1
              do l=1,nmax
               li=indices(l)
               if (li.gt.2) then ! correction sauf pour le sel
                call correction_brines(scal(i,j,ks2,li)                 &
                    ,sum_brines(l),vfonte,vcell,tms(i,j,ks2))
                somme_test(l)=somme_test(l)                             &
                 +sum_brines(l) * vcell / vfonte
               endif
              enddo
            endif
          enddo
        endif
      enddo
 
!      do l=1,nb_tracer
!        print*,'test somme',l,sum_brines(l),somme_test(l)
!        print*,'test somme dif',sum_brines(l)-somme_test(l)
!      enddo

!      write(mouchard_id,*) 'totsalts, totsaltf',totsalts, totsaltf
!      print*, 'totsalts, totsaltf',totsalts, totsaltf
!      if ((totsalts+totsaltf).ne.0) write(mouchard_id,*)
!     >   'totsalts+totsaltf different de zero !!',totsalts+totsaltf

! correction de la salinite globale : calcul de la salinite moyenne
! globale
      totsalt2=0.0
      vol=0.0
      do i=1,imax
        do j=1,jmax
          do k=1,kmax
          vol=vol+aire(i,j)*dz(k)*tms(i,j,k)
          totsalt2=totsalt2+scal(i,j,k,2)*aire(i,j)*dz(k)*tms(i,j,k)
          enddo
        enddo
      enddo
      vol_global=vol
      totsalt2=totsalt2/vol
!      write(mouchard_id,*) 'salinity global 2', totsalt2
!      print*, 'salinity global 2', totsalt2

! correction de la salinite partout
!      if (abs(totsalt1-totsalt2).ge.0) then
      if (abs(salinite_ref-totsalt2).ge.0) then
        dif=(salinite_ref-totsalt2)
!        write(mouchard_id,*) 'difference salinity', dif
!        print*, 'difference salinity', dif
        do i=1,imax
         do j=1,jmax
          do k=1,kmax
           scal(i,j,k,2)=scal(i,j,k,2)+dif
          enddo
         enddo
        enddo
      endif

! verification apres correction
      totsalt3=0.0
      vol=0.0
      do i=1,imax
        do j=1,jmax
          do k=1,kmax
          vol=vol+aire(i,j)*dz(k)*tms(i,j,k)
          totsalt3=totsalt3+scal(i,j,k,2)*aire(i,j)*dz(k)*tms(i,j,k)
          enddo
        enddo
      enddo
      totsalt3=totsalt3/vol
      write(mouchard_id,*) 'salinity global', totsalt3
!      print*, 'salinity global 3', totsalt2, totsalt3
!     >    ,salinite_ref-totsalt3



      return

      end subroutine brines
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE surface_brines(Xsurf,k,frac,sum_X,tms,Vsurf)

      implicit none

      REAL Xsurf,k,frac,sum_X,tms,Vsurf

      sum_X = sum_X + k * Xsurf * Vsurf * tms

      Xsurf = Xsurf                                                     &
            + ( 1- frac ) * k * Xsurf * tms ! * Vsurf / Vsurf !enrichissement


      RETURN

      end subroutine surface_brines
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE fond_brines(Xbot,k,frac,Xsurf,Vbot,Vsurf,tms)

      implicit none

      REAL Xbot,k,frac,Xsurf,Vbot,Vsurf,tms

      Xbot = Xbot                                                       &
            + frac * k * Xsurf * Vsurf / Vbot * tms ! enrichissement

      RETURN

      end subroutine fond_brines
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      SUBROUTINE correction_brines(Xsurf,sum_X,Vfonte,Vsurf,tms)

      implicit none

      REAL Xsurf,sum_X,Vfonte,Vsurf,tms

!      Xsurf = Xsurf - sum_X * Vsurf / Vfonte * tms
      Xsurf = Xsurf - sum_X / Vfonte * tms ! * Vsurf / Vsurf
!      Xsurf = Xsurf - sum_X / Vfonte

      RETURN

      end subroutine correction_brines
!-----------------------------------------------------------------------
      subroutine read_frac2D(frac_filename)

      use ncio, only: nc_read

      implicit none

      integer i
      real, dimension(imax-2, jmax) :: frac_2Dtemp
      !real, dimension( jmax, imax-2) :: frac_2Dtemp
      !print*, 'frac_2Dtemp', frac_2Dtemp
      character(str_len) :: frac_filename

      call nc_read(frac_filename,"frac",frac_2Dtemp)
      !print*, 'frac_2Dtemp', frac_2Dtemp

      !nb tableau a bord replie :on recopie les colonnes croisees
      ! valeur(1)=valeur(nmax-1) et
      ! valeur(nmax)=valeur(2)
      do i=2,imax-1
        frac_2D(i,:)=frac_2Dtemp(i-1,:)
      enddo
      frac_2D(1,:)=frac_2D(imax-1,:)
      frac_2D(imax,:)=frac_2D(2,:)

      return

      end subroutine read_frac2D

!-----------------------------------------------------------------------
      subroutine read_frac1D
      
       use comunit, only: iuo
       use global_constants_mod, only: ip

       implicit none

       integer i
       integer(ip) :: Frac_scenario_txt_id

        !write(*,*) 'frac_nb= ', frac_nb
        open(newunit=Frac_scenario_txt_id,                              &
             file='inputdata/clio/Frac_scenario.txt')
        do i=1,frac_nb
         read(Frac_scenario_txt_id,*) frac_date(i), frac_value(i)
         write(mouchard_id,*) 'frac BRINES ' &
                            , frac_date(i), frac_value(i), i
        enddo

        close(Frac_scenario_txt_id)


      return

      end subroutine read_frac1D
!-----------------------------------------------------------------------
      subroutine get_brines_frac

      use comemic_mod, only: iyear
      use palaeo_timer_mod, only: palaeo_year
      
      implicit none

      integer i
      real, dimension(imax-2, jmax) :: frac_2Dtemp
     
      la_date_brines = palaeo_year - iyear

      if (la_date_brines.ge.frac_date(1)) then
        frac_1D = frac_value(1)
      elseif (la_date_brines.le.frac_date(frac_nb)) then
        frac_1D = frac_value(frac_nb)
      else
         do i=2,frac_nb
           if (la_date_brines.le.frac_date(i-1).and. &
               la_date_brines.gt.frac_date(i)) then
             frac_1D = frac_value(i-1) + (frac_value(i)-frac_value(i-1)) &
                       /(frac_date(i-1)-frac_date(i)) * (frac_date(i-1)-la_date_brines)
           endif
         enddo
      endif

      write(mouchard_id,*) 'frac for brines', frac_1D, la_date_brines &
                          ,palaeo_year, iyear

      return

      end subroutine get_brines_frac
!-----------------------------------------------------------------------



      end module brines_mod
