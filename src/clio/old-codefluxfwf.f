









#if ( 0 )
      lawrflux = 0.05   ! St. Lawrence Outlet
      hudsonflux = 0.04 ! Hudson Strait+River
      waisflux = 0.0    ! West Antarc. Ice Sheet
!     HUDSON Strait+River
      do i=99,100
         do j=48,51
            fluxarea=fluxarea+(area(i,j)*tms(i,j,ks2))
         enddo
      enddo
      do i=99,100
         do j=48,51
!     [m/step] = [m3/s]/[m2]*[s/step]
            fluxinmeters(i,j) =
     &           ((hudsonflux*1.E6/fluxarea)*ddtb)
         enddo
      enddo

!     St. Lawrence Outlet
      fluxarea=0.0
      do i=97,97
         do j=44,44
            fluxarea=fluxarea+(area(i,j)*tms(i,j,ks2))
         enddo
      enddo
      do i=97,97
         do j=44,44
!     [m/step] = [m3/s]/[m2]*[s/step]
            fluxinmeters(i,j) =
     &           ((lawrflux*1.E6/fluxarea)*ddtb)
         enddo
      enddo

!     West Antarctic Ice Sheet
!     not all cells are ocean!!
      fluxarea=0.0
      do i=54,110
         do j=5,9
            fluxarea=fluxarea+(area(i,j)*tms(i,j,ks2))
         enddo
      enddo
      do i=54,110
         do j=5,9
!     [m/step] = [m3/s]/[m2]*[s/step]
            fluxinmeters(i,j) =
     &           ((waisflux*1.E6/fluxarea)*ddtb)*tms(i,j,ks2)
         enddo
      enddo
#endif /* on if above to switch on / off */
