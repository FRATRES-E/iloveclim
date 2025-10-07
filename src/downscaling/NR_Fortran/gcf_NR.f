       SUBROUTINE gcf(gammcf,a,x,gln)
       INTEGER ITMAX
       REAL a,gammcf,gln,x,EPS,FPMIN
C       PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
       PARAMETER (ITMAX=1000,EPS=3.e-5,FPMIN=1.e-30)
C USES gammln
       ! Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation
       ! as gammcf. Also returns ln Γ(a) as gln.
       ! Parameters: ITMAX is the maximum allowed number of iterations; EPS is the relative accuracy;
       ! FPMIN is a number near the smallest representable floating-point number.
       INTEGER i
       REAL an,b,c,d,del,h,gammln
       gln=gammln(a)
       b=x+1.-a ! Set up for evaluating continued fraction by modified
       c=1./FPMIN ! Lentz’s method (§5.2) with b0 = 0.
       d=1./b
       h=d
       do i=1,ITMAX ! Iterate to convergence.
       an=-i*(i-a)
       b=b+2.
       d=an*d+b
       if(abs(d).lt.FPMIN)d=FPMIN
       c=b+an/c
       if(abs(c).lt.FPMIN)c=FPMIN
       d=1./d
       del=d*c
       h=h*del
       if(abs(del-1.).lt.EPS)goto 1
       enddo
       write(*,*) "a too large, ITMAX too small in gcf"
       read(*,*)
    1  gammcf=exp(-x+a*log(x)-gln)*h ! Put factors in front.
       return
       END
