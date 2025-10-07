      FUNCTION gammq(a,x)
      REAL a,gammq,x
C USES gcf,gser
      ! Returns the incomplete gamma function Q(a, x) ≡ 1 − P(a, x).
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.) then
        write(*,*) "bad arguments in gammq"
        read(*,*)
      endif
      if(x.lt.a+1.)then ! Use the series representation
      call gser(gamser,a,x,gln)
      gammq=1.-gamser ! and take its complement.
      else ! Use the continued fraction representation.
      call gcf(gammcf,a,x,gln)
      gammq=gammcf
      endif
      return
      END
