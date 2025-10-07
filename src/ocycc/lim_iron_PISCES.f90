      MODULE IRONLIM_PISCES


      use global_constants_mod, only: str_len, ip, dblp=>dp
      use para0_mod,            only: imax, jmax, kmax

      IMPLICIT NONE

      PUBLIC :: get_PISCES_CLIO_IRONLIM

      INTEGER(kind=ip), PARAMETER :: prof = 6

      CONTAINS

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      FUNCTION read_iron_lim_CLIO() result(ironLim_CLIO)

        use ncio, only: nc_read

        character(len=str_len), parameter :: filename_NC="../compil/sources/LFe.nc"

        real(kind=dblp), dimension(imax,jmax,prof) :: ironLim_CLIO, fixed_iron_CLIO
        integer, parameter  :: imaxm2 = imax-2, imaxm1 = imax-1
        real(kind=dblp), dimension(imaxm2,jmax) :: ironLim_to_read

        integer :: i,checked_mask
        logical :: f_exists
        character(len=6), dimension(prof) :: var_names =["LFe5m ","LFe15m","LFe29m", "LFe45m","LFe64m","LFe89m"]

        inquire(file=filename_NC,exist=f_exists)

        if (f_exists) then

            ! Boucle sur les variables 
            do i = 1, size(var_names) 
 
               ! proceed with opening and reading
               call nc_read(filename_NC,TRIM(var_names(i)),ironLim_to_read)

               ironLim_CLIO(2:imaxm1,:,i) = ironLim_to_read(:,:)
               ironLim_CLIO(1,:,i) = ironLim_to_read(imaxm2,:)
               ironLim_CLIO(imax,:,i) = ironLim_to_read(1,:)

            enddo

        else
               WRITE(*,*) "File "//TRIM(filename_NC)//" does not exist! [STOP]"
               stop
        endif
 

        checked_mask = check_ERA5_CLIO_IRON(ironLim_CLIO,quick_fix=.TRUE.,fixed_ironfield=fixed_iron_CLIO)

        ! Modify the input field if desired
        ironLim_CLIO(:,:,:) = fixed_iron_CLIO(:,:,:)

      END FUNCTION read_iron_lim_CLIO
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       FUNCTION get_PISCES_CLIO_IRONLIM() result(ironLim_CLIO_depth)

        real(kind=dblp), dimension(imax,jmax,prof) :: ironLim_CLIO
        real(kind=dblp), dimension(imax,jmax,kmax) :: ironLim_CLIO_depth

        integer :: i,j,k,l

         ironLim_CLIO(:,:,:) = read_iron_lim_CLIO()
         
         do i=1,imax
           do j=1,jmax
             do k=1,kmax

             ironLim_CLIO_depth(i,j,k) = 1.0 ! pas de lim en profondeur             

             if ( k >= 15 ) then ! Limitation en surface
                 l = kmax-k+1 
                 ironLim_CLIO_depth(i,j,k) = ironLim_CLIO(i,j,l)
             endif

             !WRITE(*,*) "ironLIm", ironLim_CLIO_depth(i,j,20)
             !WRITE(*,*) "ironLIm bef", ironLim_CLIO(i,j,1)

             enddo
           enddo
         enddo
         
         WRITE(*,*) "MINMAX :::", MINVAL(ironLim_CLIO_depth), MAXVAL(ironLim_CLIO_depth)

       END FUNCTION get_PISCES_CLIO_IRONLIM
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----

       FUNCTION check_ERA5_CLIO_IRON(ironLim_CLIO,quick_fix,fixed_ironfield) result(nb_probs)
          use bloc0_mod, only: tms

          real(kind=dblp), dimension(imax,jmax,prof), intent(in) :: ironLim_CLIO
          real(kind=dblp), dimension(imax,jmax,prof), intent(out), optional :: fixed_ironfield
          logical                                                , intent(in),  optional :: quick_fix

          real(kind=dblp), parameter :: acceptable_limit=1.E10_dblp
          real(kind=dblp) :: summed_iron

          integer :: p,i,j,k,l,m,n,t, summed_cells
          integer :: nb_probs
          integer, dimension(:,:), allocatable :: probs_location
          integer, parameter :: search_range = 3

          nb_probs = 0
  
          do p=1,prof
           !k=UBOUND(tms,dim=3)-prof+1
           k=UBOUND(tms,dim=3)-p+1

           do j=1,UBOUND(tms,dim=2)
            do i=1,UBOUND(tms,dim=1)
              
                if ((tms(i,j,k).gt.0.).and.(ironLim_CLIO(i,j,p).ge.acceptable_limit)) then
                   nb_probs = nb_probs + 1
                endif

               enddo   
            enddo
          enddo

          if (nb_probs.gt.0) then
            WRITE(*,*) "SUMMARY PROBLEM GRID LIM IRON", nb_probs
          endif


          ! Application de la correction : 
          
          if (present(quick_fix).and.present(fixed_ironfield)) then

          !QUICKFIX: if ((nb_probs.lt.100).and.quick_fix) then
          QUICKFIX: if ((nb_probs.lt.1000).and.quick_fix) then

            fixed_ironfield(:,:,:) = ironLim_CLIO(:,:,:)
            allocate(probs_location(2,nb_probs))
            
           DEPTH : do p=1,prof
           ! k=UBOUND(tms,dim=3)-prof+1
            k=UBOUND(tms,dim=3)-p+1

            nb_probs = 0

            do j=1,UBOUND(tms,dim=2)
              do i=1,UBOUND(tms,dim=1)
                
                  if ((tms(i,j,k).gt.0.).and.(ironLim_CLIO(i,j,p).ge.acceptable_limit)) then
                     nb_probs = nb_probs + 1
                     probs_location(1,nb_probs) = i
                     probs_location(2,nb_probs) = j
                  endif

              enddo
            enddo

            PROBS: do n=1,nb_probs
              summed_cells = 0
              summed_iron = 0.0_dblp

              do m=probs_location(2,n)-search_range,probs_location(2,n)+search_range
                if ((m.LT.LBOUND(tms,dim=2)).OR.(m.GT.UBOUND(tms,2))) CYCLE

                do l=probs_location(1,n)-search_range,probs_location(1,n)+search_range
                  if ((l.LT.LBOUND(tms,dim=1)).OR.(l.GT.UBOUND(tms,1))) CYCLE

                     if (ironLim_CLIO(l,m,p).lt.acceptable_limit) then
                      summed_iron = summed_iron + ironLim_CLIO(l,m,p)
                      summed_cells = summed_cells + 1
                     endif

                enddo
              enddo
              
             if (summed_cells.gt.0) then
                  fixed_ironfield(probs_location(1,n),probs_location(2,n), p) = &
                        summed_iron / summed_cells
             else
                  WRITE(*,*) "Unfixable error in ERA5 @", probs_location(1,n),probs_location(2,n) &
                            , k, ironLim_CLIO(probs_location(1,n),probs_location(2,n),p)
             endif

            enddo PROBS
            
            enddo DEPTH

          endif QUICKFIX

          endif

      END FUNCTION check_ERA5_CLIO_IRON
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      END MODULE IRONLIM_PISCES
