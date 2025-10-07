!dmr -- Adding the choice of components through the pre-processing options
#include "choixcomposantes.h"
!dmr -- Adding the choice of components through the pre-processing options

       module palaeo_timer_mod

#if ( F_PALAEO == 3 )
         use input_subgrid2L, only: nb_steps_fism, update_time_fism, indx_fism
! afq --- these variables should ideally be here instead in input_GRISLI2L
!         but this is not possible (yet?) because of compilation order
#endif
         
!-dmr to Be moved to global module data, private
       integer(kind=4), parameter, private :: time_max = 3197           &
               , update_time = 250, nb_steps=3200, timeghg_max = 798000
#if ( F_PALAEO == 3 )
       integer(kind=4), parameter, private :: time_max_fism = 161
#endif

!-dmr --- to be moved as global, saved variables, public
       integer,save, public :: ref_irunlabel, accel_rt, palaeo_year     &
               , equil_flag

       integer(kind=4), parameter, public  ::                           &
                 reload_topo = update_time * 360

! dmr these are the outpout computations of the subroutine timer ...
! fl added indx_topo to disentangle the effect of the ice sheets topo & albedo
       integer(kind=4), save, public :: indx_masq, indx_topo, indx_ghg, celest_year
       contains

       subroutine palaeo_timer(init,timer,n_days)

#if ( UNCORFLUX >= 1 )
       use zonesormenfwf_mod, only: set_datesBP
#endif

       implicit none

       integer, intent(in) :: init, timer
       integer, intent(in), optional :: n_days


!-dmr --- these are the local variables ...
       integer :: indx_first, indx_ghgfirst
#if ( F_PALAEO == 3 )
       integer :: indx_first_fism
#endif

       if ((init.EQ.1).or.(init.eq.2)) then

        OPEN(1234,file="palaeo_transient_parameters.dat")

        READ(1234,*)
        READ(1234,*) palaeo_year

        READ(1234,*)
        READ(1234,*)
        READ(1234,*)
        READ(1234,*) ref_irunlabel


        READ(1234,*)
        READ(1234,*)
        READ(1234,*)
        READ(1234,*)
        READ(1234,*) accel_rt

        READ(1234,*)
        READ(1234,*)
        READ(1234,*)
        READ(1234,*) equil_flag

        CLOSE(1234)

       endif

       if (init.ne.2) then
!dmr --- to get the righ masqs ...

! Current files contains 801,000 years of data with a 250 yrs step.▫
!         files start at 800,000 kyrs BP and stops at -750 yrs BP
! BEWARE: only 3197 datapoints make sense ...
! Start at 150,000 yrs B.P. is thus: 3200.-(150,000./250.-1.) = 2601

!-dmr --- indx_first is the index at the ref_irunlabel
       indx_first = MIN(MAX(nb_steps-(palaeo_year/update_time-1),1),time_max)
!-dmr --- The GHG.dat file contains 798000 lines of data, the first one being 798000 B.P.
       indx_ghgfirst = MIN(MAX(timeghg_max-(palaeo_year-1),1),timeghg_max)

!-dmr --- In the following timer is meant to be irunlabelf+iyear
       indx_masq = MIN(  indx_first +                                   &
             ((timer-ref_irunlabel)*accel_rt/update_time)*(1-equil_flag)&
                 , time_max )
! fl tested indx_masq/topo = time_max to assess the effect of the ice sheets topo or albedo only
! fl These index are taken into account in ec_masq.f to update the icemask or the topography
! fl Currently indx_topo is not used in ec_masq.f (Ganopolski only, currently Peltier)
! fl By default they are equivalent
       indx_topo = indx_masq
!       indx_masq = time_max
!       indx_topo = time_max

#if ( F_PALAEO == 3 )
       indx_first_fism = MIN(MAX(nb_steps_fism-(palaeo_year/update_time_fism),1),time_max_fism)
       indx_fISM = MIN(  indx_first_fism +                              &
             ((timer-ref_irunlabel)*accel_rt/update_time_fism)*(1-equil_flag)&
                 , time_max_fism )
#endif

       celest_year = palaeo_year - (timer-ref_irunlabel)*accel_rt*(1-equil_flag)


       indx_ghg = MIN(  indx_ghgfirst +                                 &
             ((timer-ref_irunlabel)*accel_rt)*(1-equil_flag)            &
                 , timeghg_max )

       !write(*,*) timer, indx_masq, celest_year, indx_ghg

#if ( UNCORFLUX >= 1 )
       if (present(n_days).and.(init.eq.1)) then
         call set_datesBP(celest_year,celest_year-n_days/360)
       endif
#endif

       endif !init ne 2

       end subroutine palaeo_timer

       end module palaeo_timer_mod
