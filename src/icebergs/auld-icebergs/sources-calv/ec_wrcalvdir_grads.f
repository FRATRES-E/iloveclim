!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to convert mass that has not been used to produce icebergs into  !
! format that can be read by grads to see where the calving takes place       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE ec_wrcalvdir_grads

#include "para.com"

      integer :: ii,jj,irec,ijm
      real (kind = 4) ::mr(imax,jmax)
      real (kind = 8) :: mass(imax,jmax)


      ijm=Size(mr)*Kind(mr(1,1))

!      write(*,*)'ec_wrcal_grads: imax,jmax',imax,jmax
      
      open(15,file='outputdata/icebergs/calvdir_grads.dat'
     &,form='unformatted',status='old')
      read(15) mass
      close(15)

      do ii=1,imax
        do jj= 1,jmax
           mr(ii,jj)=mass(ii,jj)
        enddo
      enddo


      open(8,file='outputdata/icebergs/calvdir_bin.dat'
     &,form='unformatted',access='direct',recl=ijm,status='replace')
      open(20,file='calvdir.dat',form='formatted',status='replace')
      irec=1
      write(8,rec=irec) ((mr(ii,jj),ii=1,imax),jj=1,jmax)
      write(20,*) ((mr(ii,jj),ii=1,imax),jj=1,jmax)

      close(8)
      close(20)
      
      return 
         
      END SUBROUTINE ec_wrcalvdir_grads
