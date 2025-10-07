!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        module ec_detqmax_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         use global_constants_mod, only: dblp=>dp, ip

         implicit none


         contains
       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        function ec_detqmax(tmount,i,j,dqmdt)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! *** determines the maximum water content in latlon point
! *** i,j for given ground- and 650 and 350 hPa temperature
! *** by linear interpolation in qmtabel
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      use comatm, only: alogpl2tl2, nlat, nlon, nsh2, nvl, alogtl12, alogtl1pl2, grav, nsh, nm, ntl, rgas, rlogtl12
      use comdyn, only: geopg
      use comphys, only: tqmi, tqmj, tqmk, iqmtab, jqmtab, kqmtab, tqmimin, tqmjmin, tqmkmin, qmtabel, rdtqmi, rdtqmj, rdtqmk   &
                 , gpm500, hmoisr, temp4g, temp2g, qmount
      
      use comunit

      use OMP_LIB

      implicit none

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=dblp), intent(in)  :: tmount
      integer(kind=ip),intent(in)  :: i,j
      real(kind=dblp), intent(out) :: dqmdt

      integer(kind=ip)             :: ii,jj,kk
      real(kind=dblp)              :: ti,tj,tk
      real(kind=dblp)              :: qmax,ec_detqmax,dqmdi,dqmdj,dqmdk
      real(kind=dblp)              :: hmount,z500,t500,alpha,dtgdt

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


      ti=temp4g(i,j)
      tj=tmount-temp4g(i,j)
      tk=temp4g(i,j)-temp2g(i,j)


!~ c~       if (ti.lt.tqmi(0)) then
!~ c~         ti=tqmi(0)
!~ c~ !        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!~ c~ !        call ec_error(121)
!~ c~       endif
!~ c~       if (ti.gt.tqmi(iqmtab)) then
!~ c~         ti=tqmi(iqmtab)
!~ c~ !        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!~ c~ !        call ec_error(121)
!~ c~       endif

      ti = min(max(ti,tqmi(0)),tqmi(iqmtab))


!~ c~       if (tj.lt.tqmj(0)) then
!~ c~         tj=tqmj(0)
!~ c~ !        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!~ c~ !        call ec_error(121)
!~ c~       endif
!~ c~       if (tj.gt.tqmj(jqmtab)) then
!~ c~         tj=tqmj(jqmtab)
!~ c~ !        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!~ c~ !        call ec_error(121)
!~ c~       endif


      tj = min(max(tj,tqmj(0)),tqmj(jqmtab))


!~ c~       if (tk.lt.tqmk(0)) then
!~ c~         tk=tqmk(0)
!~ c~ !        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!~ c~ !        call ec_error(121)
!~ c~       endif
!~ c~       if (tk.gt.tqmk(kqmtab)) then
!~ c~         tk=tqmk(kqmtab)
!~ c~ !        write(29,*) 'in latlon ',i,j,' t500 ',t500,' tmount ',tmount
!~ c~ !        call ec_error(121)
!~ c~       endif


      tk = min(max(tk,tqmk(0)),tqmk(kqmtab))


      ii=min(iqmtab-1,int((ti-tqmimin)*rdtqmi))
      jj=min(jqmtab-1,int((tj-tqmjmin)*rdtqmj))
      kk=min(kqmtab-1,int((tk-tqmkmin)*rdtqmk))

      dqmdi=(qmtabel(ii+1,jj,kk)-qmtabel(ii,jj,kk))*rdtqmi
      dqmdj=(qmtabel(ii,jj+1,kk)-qmtabel(ii,jj,kk))*rdtqmj
      dqmdk=(qmtabel(ii,jj,kk+1)-qmtabel(ii,jj,kk))*rdtqmk

    
      qmax = qmtabel(ii,jj,kk) + (ti-tqmi(ii))*dqmdi + (tj-tqmj(jj))*dqmdj + (tk-tqmk(kk))*dqmdk

      qmax = min(max(qmax,0d0),0.2d0)

!~ c~       if (qmax.lt.0d0) qmax=0d0


!~ c~       if (qmax.gt.0.2) then
!~ c~ !dmr --- Tentative fix to try and stop the crash of model
!~ c~         qmax = 0.2
!~ c~ !        write(iuo+29,*) 'in latlon ',i,j,' qmax ',qmax, 'tentative fix'
!~ c~         write(iuo+29,*) 'in latlon ',i,j,' qmax ',qmax, 'tentative fix'
!~ c~ !dmr --- Tentative fix to try and stop the crash of model
!~ c~ !        call ec_error(121)
!~ c~       endif


      alpha=(temp2g(i,j)-temp4g(i,j))*rlogtl12
      t500=temp4g(i,j)+alpha*alogpl2tl2
      z500=gpm500*grav
      hmount=qmount(i,j)*hmoisr*grav


      dtgdt=(rgas*t500*alogtl1pl2 + (hmount-geopg(i,j,2)-z500))/(rgas*tmount*alogtl12)


      dqmdt=dqmdi + dqmdj * (dtgdt - 1d0) + dqmdk


      ec_detqmax=0.9*qmax
!     ec_detqmax=0.85*qmax

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        end function ec_detqmax
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|     
       end module ec_detqmax_mod
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
