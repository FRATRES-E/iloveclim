subroutine linear_regression(nbt,nmin,nmax,x,y,intercept,slope)

   implicit none

   integer                        ,intent(in)  :: nbt       ! the length of the input arrays
   integer                        ,intent(in)  :: nmin      ! we might want to compute the reg. on a part of x,y, nmin
   integer                        ,intent(in)  :: nmax      ! we might want to compute the reg. on a part of x,y, nmax
   double precision,dimension(nbt),intent(in)  :: x         !
   double precision,dimension(nbt),intent(in)  :: y         !
   double precision               ,intent(out) :: intercept ! y-intercept of least-squares best fit line
   double precision               ,intent(out) :: slope     ! slope of least-squares best fit line
   
   integer          ::  nbp            ! number of considered points
   double precision ::  r              ! squared correlation coefficient
   double precision ::  sumx           ! sum of x
   double precision ::  sumx2          ! sum of x**2
   double precision ::  sumxy          ! sum of x * y
   double precision ::  sumy           ! sum of y
   double precision ::  sumy2          ! sum of y**2

   integer          ::  nb             ! loop integer
   
   nbp = nmax - nmin + 1

   if ((nbt.lt.nbp).or.(nbp.le.0.)) then
      write (*,*) "Linear regression: X and Y are shorter than required length!"
      STOP
   endif

   sumx = 0.d0
   sumx2 = 0.0d0
   sumxy = 0.0d0
   sumy  = 0.0d0
   sumy2 = 0.0d0
   
   do nb=nmin,nmax                  ! loop for all data points
      sumx  = sumx + x(nb)          ! compute sum of x
      sumx2 = sumx2 + x(nb) * x(nb) ! compute sum of x**2
      sumxy = sumxy + x(nb) * y(nb) ! compute sum of x * y
      sumy  = sumy + y(nb)          ! compute sum of y
      sumy2 = sumy2 + y(nb) * y(nb) ! compute sum of y**2
   end do

   if ((nbp * sumx2 - sumx**2).lt.1.d-9) then
      
      write (*,*) "Linear regression: X does not vary!"
      !STOP
      slope = 0
      intercept = sumx/nbp
      r = 1
      
   else

      slope = (nbp * sumxy  -  sumx * sumy) / (nbp * sumx2 - sumx**2)            ! compute slope
      intercept = (sumy * sumx2  -  sumx * sumxy) / (nbp * sumx2  -  sumx**2)    ! compute y-intercept
      r = (sumxy - sumx * sumy / nbp) /                                     &    ! compute correlation coefficient
                     sqrt((sumx2 - sumx**2/nbp) * (sumy2 - sumy**2/nbp))
      
   endif
   
   return
end subroutine linear_regression
