       SUBROUTINE fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
       IMPLICIT NONE
       INTEGER mwt,ndata
       REAL a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
! USES gammq
! Given a set of data points x(1:ndata),y(1:ndata) with individual standard deviations
! sig(1:ndata), fit them to a straight line y = a + bx by minimizing χ2. Returned are
! a,b and their respective probable uncertainties siga and sigb, the chi-square chi2, and
! the goodness-of-fit probability q (that the fit would have χ2 this large or larger). If mwt=0
! on input, then the standard deviations are assumed to be unavailable: q is returned as 1.0
! and the normalization of chi2 is to unit standard deviation on all points.
       INTEGER i
       REAL sigdat,ss,st2,sx,sxoss,sy,t,wt,gammq
       sx=0. ! Initialize sums to zero.
       sy=0.
       st2=0.
       b=0.

       if(mwt.ne.0) then ! Accumulate sums ...
       ss=0.
       do i=1,ndata ! ...with weights
       wt=1./(sig(i)**2)
       ss=ss+wt
       sx=sx+x(i)*wt
       sy=sy+y(i)*wt
       enddo 
       else
       do i=1,ndata ! ...or without weights.
       sx=sx+x(i)
       sy=sy+y(i)
       enddo
       ss=float(ndata)
       endif
       sxoss=sx/ss
       if(mwt.ne.0) then
       do i=1,ndata
       t=(x(i)-sxoss)/sig(i)
       st2=st2+t*t
       b=b+t*y(i)/sig(i)
       enddo
       else
       do i=1,ndata
       t=x(i)-sxoss
       st2=st2+t*t
       b=b+t*y(i)
       enddo
       endif
       b=b/st2 ! Solve for a, b, σa, and σb.
       a=(sy-sx*b)/ss
       siga=sqrt((1.+sx*sx/(ss*st2))/ss)
       sigb=sqrt(1./st2)
       chi2=0. ! Calculate χ2
       q=1.
       if(mwt.eq.0) then
       do i=1,ndata
       chi2=chi2+(y(i)-a-b*x(i))**2
       enddo
       sigdat=sqrt(chi2/(ndata-2)) ! For unweighted data evaluate typical sig using
! chi2, and adjust the standard deviations.
       siga=siga*sigdat
       sigb=sigb*sigdat
       else
       do i=1,ndata
       chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
       enddo
       if(ndata.gt.2) q=gammq(0.5*(ndata-2),0.5*chi2) ! Equation (15.2.12).
       endif
       return
       END
