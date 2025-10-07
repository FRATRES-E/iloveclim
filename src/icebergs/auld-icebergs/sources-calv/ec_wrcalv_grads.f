!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to convert mass that has not been used to produce icebergs into  !
! format that can be read by grads to see where the calving takes place       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE ec_wrcalv_grads

#include "para.com"

      integer :: ii,jj,irec,ijm
      real (kind = 4) ::mr(imax,jmax)
      real (kind = 8) :: mass(imax,jmax)
      logical ::  lexist


      ijm=Size(mr)*Kind(mr(1,1))

!      write(*,*)'ec_wrcal_grads: imax,jmax',imax,jmax
      inquire(file='outputdata/icebergs/mass_grads.dat',exist=lexist)

      if (lexist) then
      open(15,file='outputdata/icebergs/mass_grads.dat'
     &,form='unformatted',status='old')
      read(15) mass
      close(15)

      do ii=1,imax
        do jj= 1,jmax
           mr(ii,jj)=mass(ii,jj)
        enddo
      enddo


      open(8,file='outputdata/icebergs/massbin_grads.dat'
     &,form='unformatted',access='direct',recl=ijm,status='replace')
      open(20,file='wrcal_grads.dat',form='formatted',status='replace')
      irec=1
      write(8,rec=irec) ((mr(ii,jj),ii=1,imax),jj=1,jmax)
      write(20,*) ((mr(ii,jj),ii=1,imax),jj=1,jmax)

      close(8)
      close(20)
      
      else
      print*, 'no mass_grads.dat!!!'
      endif

      return 
         
      END SUBROUTINE ec_wrcalv_grads
