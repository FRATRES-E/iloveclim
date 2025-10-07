       SUBROUTINE gser(gamser,a,x,gln)
       INTEGER ITMAX
       REAL a,gamser,gln,x,EPS
       PARAMETER (ITMAX=100,EPS=3.e-7)
C USES gammln
       ! Returns the incomplete gamma function P(a, x) evaluated by its series representation as
       ! gamser. Also returns lnΓ(a) as gln.
       INTEGER n
       REAL ap,del,sum,gammln
       gln=gammln(a)
       if(x.le.0.)then
       if(x.lt.0.) then
         write(*,*) "x < 0 in gser"
         read(*,*)
       endif
       gamser=0.
       return
       endif
       ap=a
       sum=1./a
       del=sum
       do n=1,ITMAX
       ap=ap+1.
       del=del*x/ap
       sum=sum+del
       if(abs(del).lt.abs(sum)*EPS)goto 1
       enddo
       write(*,*) "a too large, ITMAX too small in gser"
       read(*,*)
    1  gamser=sum*exp(-x+a*log(x)-gln)
       return
       END
