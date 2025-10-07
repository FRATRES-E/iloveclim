      MODULE WINDFORC_CLIOERA5


      use global_constants_mod, only: str_len, dblp=>dp, months_year_i, days_year360d_i, days_year365d_i
      use para0_mod,            only: imax, jmax

      IMPLICIT NONE

      PUBLIC :: get_daily_ERA5_CLIO_WINDSSTRENGTH

      CONTAINS

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      FUNCTION get_daily_ERA5_CLIO_WINDS() result(winds_CLIO_uv)

        use ncio, only: nc_read

        character(len=str_len), parameter :: filename_NC_uwinds="../compil/sources/climato_Ucomponent_1930_1960-CLIO-T-masked.nc"
        character(len=str_len), parameter :: filename_NC_vwinds="../compil/sources/climato_Vcomponent_1930_1960-CLIO-T-masked.nc"

        real(kind=dblp), dimension(imax,jmax,days_year360d_i,2) :: winds_CLIO_uv, fixed_winds_CLIO_uv


        integer, parameter  :: imaxm2 = imax-2, imaxm1 = imax-1
        real(kind=dblp), dimension(imaxm2,jmax,days_year365d_i+1) :: daily_winds_to_read

        integer :: checked_mask

        logical :: f_exists

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! Verifcation
!        character(len=str_len) :: absolute_path
!        call getcwd(absolute_path)
!        write(*,*) "Fotran cherche fichier :", trim(absolute_path)
! ----------

        inquire(file=filename_NC_uwinds,exist=f_exists)

        if (f_exists) then
               ! proceed with opening and reading

               call nc_read(filename_NC_uwinds,"u10m",daily_winds_to_read)

               winds_CLIO_uv(2:imaxm1,:,:,1) = daily_winds_to_read(:,:,:days_year360d_i)
               winds_CLIO_uv(1,:,:,1) = daily_winds_to_read(imaxm2,:,:days_year360d_i)
               winds_CLIO_uv(imax,:,:,1) = daily_winds_to_read(1,:,:days_year360d_i)

        else
               WRITE(*,*) "File "//TRIM(filename_NC_uwinds)//" does not exist! [STOP]"
               stop
        endif

        inquire(file=filename_NC_vwinds,exist=f_exists)

        if (f_exists) then
               ! proceed with opening and reading

               call nc_read(filename_NC_vwinds,"v10m",daily_winds_to_read)

               winds_CLIO_uv(2:imaxm1,:,:,2) = daily_winds_to_read(:,:,:days_year360d_i)
               winds_CLIO_uv(1,:,:,2) = daily_winds_to_read(imaxm2,:,:days_year360d_i)
               winds_CLIO_uv(imax,:,:,2) = daily_winds_to_read(1,:,:days_year360d_i)

        else
               WRITE(*,*) "File "//TRIM(filename_NC_vwinds)//" does not exist! [STOP]"
               stop
        endif

        checked_mask = check_ERA5_CLIO_WINDS(winds_CLIO_uv,quick_fix=.TRUE.,fixed_windfield=fixed_winds_CLIO_uv)
        ! Modify the input field if desired
        winds_CLIO_uv(:,:,:,:) = fixed_winds_CLIO_uv(:,:,:,:)

      END FUNCTION get_daily_ERA5_CLIO_WINDS
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       FUNCTION get_daily_ERA5_Ucomponent() result(winds_Ucomponent)

        real(kind=dblp), dimension(imax,jmax,days_year360d_i) :: winds_Ucomponent
        real(kind=dblp), dimension(imax,jmax,days_year360d_i,2) :: variable_uv

         variable_uv(:,:,:,:) = get_daily_ERA5_CLIO_WINDS()
         winds_Ucomponent(:,:,:) = variable_uv(:,:,:,1)

       END FUNCTION get_daily_ERA5_Ucomponent

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       FUNCTION get_daily_ERA5_Vcomponent() result(winds_Vcomponent)

        real(kind=dblp), dimension(imax,jmax,days_year360d_i) :: winds_Vcomponent
        real(kind=dblp), dimension(imax,jmax,days_year360d_i,2) :: variable_uv

         variable_uv(:,:,:,:) = get_daily_ERA5_CLIO_WINDS()
         winds_Vcomponent(:,:,:) = variable_uv(:,:,:,2)

       END FUNCTION get_daily_ERA5_Vcomponent

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       FUNCTION get_daily_ERA5_CLIO_WINDSSTRENGTH() result(winds_CLIO_normuv)

        real(kind=dblp), dimension(imax,jmax,days_year360d_i) :: winds_CLIO_normuv
        real(kind=dblp), dimension(imax,jmax,days_year360d_i,2) :: variable_uv

        integer :: i,j,k

         variable_uv(:,:,:,:) = get_daily_ERA5_CLIO_WINDS()
         do i=1,imax
         do j=1,jmax
         do k=1,days_year360d_i
           winds_CLIO_normuv(i,j,k) = SQRT(variable_uv(i,j,k,1)**2+variable_uv(i,j,k,2)**2)
           !WRITE(*,*) "fonctions vent",winds_CLIO_normuv(i,j,k)
         enddo
         enddo
         enddo
         WRITE(*,*) "MINMAX :::", MINVAL(winds_CLIO_normuv), MAXVAL(winds_CLIO_normuv)

       END FUNCTION get_daily_ERA5_CLIO_WINDSSTRENGTH

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
       FUNCTION check_ERA5_CLIO_WINDS(windsCLIOuv,quick_fix,fixed_windfield) result(nb_probs)
          use bloc0_mod, only: tms

          real(kind=dblp), dimension(imax,jmax,days_year360d_i,2), intent(in)            :: windsCLIOuv
          real(kind=dblp), dimension(imax,jmax,days_year360d_i,2), intent(out), optional :: fixed_windfield
          logical                                                , intent(in),  optional :: quick_fix

          real(kind=dblp), parameter :: acceptable_limit=1.E10_dblp
          real(kind=dblp), dimension(days_year360d_i,2) :: summed_winds

          integer :: i,j,k,l,m,n, summed_cells
          integer :: nb_probs
          integer, dimension(:,:), allocatable :: probs_location
          integer, parameter :: search_range = 3

          nb_probs = 0

          k = UBOUND(tms,dim=3) ! This means I take the surface only

          do j=1,UBOUND(tms,dim=2)
            do i=1,UBOUND(tms,dim=1)
              if ((tms(i,j,k).gt.0.).and.(windsCLIOuv(i,j,1,1).ge.acceptable_limit)) then
                 nb_probs = nb_probs + 1
              endif
            enddo
          enddo

          if (nb_probs.gt.0) then
            WRITE(*,*) "SUMMARY PROBLEM GRID WINDS ERA5", nb_probs
          endif

          if (present(quick_fix).and.present(fixed_windfield)) then

          QUICKFIX: if ((nb_probs.lt.100).and.quick_fix) then

            fixed_windfield(:,:,:,:) = windsCLIOuv(:,:,:,:)

            allocate(probs_location(2,nb_probs))
            nb_probs = 0

            do j=1,UBOUND(tms,dim=2)
              do i=1,UBOUND(tms,dim=1)
                if ((tms(i,j,k).gt.0.).and.(windsCLIOuv(i,j,1,1).ge.acceptable_limit)) then
                   nb_probs = nb_probs + 1
                   probs_location(1,nb_probs) = i
                   probs_location(2,nb_probs) = j
                endif
              enddo
            enddo


            PROBS: do n=1,nb_probs

              summed_cells = 0
              summed_winds(:,:) = 0.0_dblp

              do m=probs_location(2,n)-search_range,probs_location(2,n)+search_range

                if ((m.LT.LBOUND(tms,dim=2)).OR.(m.GT.UBOUND(tms,2))) CYCLE

                do l=probs_location(1,n)-search_range,probs_location(1,n)+search_range

                  if ((l.LT.LBOUND(tms,dim=1)).OR.(l.GT.UBOUND(tms,1))) CYCLE

                  if (windsCLIOuv(l,m,1,1).lt.acceptable_limit) then
                    summed_winds(:,1) = summed_winds(:,1) + windsCLIOuv(l,m,:,1)
                    summed_winds(:,2) = summed_winds(:,2) + windsCLIOuv(l,m,:,2)
                    summed_cells = summed_cells + 1
                  endif
                enddo
              enddo

             if (summed_cells.gt.0) then
                  fixed_windfield(probs_location(1,n),probs_location(2,n),:,:) = summed_winds(:,:) / summed_cells
             else
                  WRITE(*,*) "Unfixable error in ERA5 @", probs_location(1,n),probs_location(2,n) &
                            , windsCLIOuv(probs_location(1,n),probs_location(2,n),1,1)
             endif

            enddo PROBS

          endif QUICKFIX

          endif
       END FUNCTION check_ERA5_CLIO_WINDS

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      END MODULE WINDFORC_CLIOERA5
