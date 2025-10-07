      SUBROUTINE moment(data,n,ave,adev,sdev,var,skew,curt)
      INTEGER n
      REAL adev,ave,curt,sdev,skew,var,data(n)
      ! Given an array of data(1:n), this routine returns its mean ave, average deviation adev,
      ! standard deviation sdev, variance var, skewness skew, and kurtosis curt.
      INTEGER j
      REAL p,s,ep
      if(n.le.1) then
        write(*,*) 'n must be at least 2 in moment'
        read(*,*)
      endif
      s=0. ! First pass to get the mean.
      do j=1,n
      s=s+data(j)
      enddo
      ave=s/n
      adev=0. ! Second pass to get the first (absolute), second, third, and fourth
      var=0. ! moments of the deviation from the mean.
      skew=0.
      curt=0.
      ep=0.
      do j=1,n
      s=data(j)-ave
      ep=ep+s
      adev=adev+abs(s)
      p=s*s
      var=var+p
      p=p*s
      skew=skew+p
      p=p*s
      curt=curt+p
      enddo
      adev=adev/n ! Put the pieces together according to the conventional definitions.
      var=(var-ep**2/n)/(n-1) ! Corrected two-pass formula.
      sdev=sqrt(var)
      if(var.ne.0.)then
      skew=skew/(n*sdev**3)
      curt=curt/(n*var**2)-3.
      else
!dmr&cd      pause 'no skew or kurtosis when zero variance in moment'
      skew=9999.
      curt=9999.
!dmr&cd
      endif
      return
      END
