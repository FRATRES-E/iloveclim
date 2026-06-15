!==============================================================================
! MODULE atm_thermo_mod
! Pure thermodynamic functions for the ECBilt QG3L atmospheric model.
!
! Contains functions with no side effects, usable by any atmospheric module
! without creating circular dependencies:
!   ec_qsat       — saturation mixing ratio
!   ec_levtempgp  — temperature at a given pressure level (grid-point)
!   ec_expint     — exponential integral E_n(x)
!
! Extracted from atmphys_mod (atmphys0.f) to break the circular dependency
! between atmphys_mod and atmmois_mod.
!==============================================================================
module atm_thermo_mod

  use global_constants_mod, only: dblp=>dp, silp=>sp, ip

  implicit none
  private

  public :: ec_qsat, ec_levtempgp, ec_expint, ec_surfmois

contains

!------------------------------------------------------------------------------
!  FUNCTION ec_qsat
!> @brief Saturation mixing ratio
!> @param[in]  press  pressure [Pa]
!> @param[in]  temp   temperature [K]
!> @return     saturation mixing ratio
!------------------------------------------------------------------------------
      function ec_qsat(press,temp) result(qsat)

      use comphys, only: tzero,cc1,cc2,cc3

      implicit none

      real(kind=dblp), intent(in) :: press, temp
      real(kind=dblp)             :: qsat

      qsat = cc1*exp(cc2*(temp-tzero)/(temp-cc3)) &
           / press

      end function ec_qsat


!------------------------------------------------------------------------------
!  FUNCTION ec_levtempgp
!> @brief Temperature at pressure level plev [Pa] at grid point (i,j)
!> @details Assumes temperature varies linearly with log(p)
!> @param[in]  plev  pressure level [Pa]
!> @param[in]  i     latitude index
!> @param[in]  j     longitude index
!> @return     temperature [K]
!------------------------------------------------------------------------------
      function ec_levtempgp(plev,i,j) result(tlev)

      use comatm,  only: rlogtl12
      use comphys, only: temp2g, temp4g

      implicit none

      real(kind=dblp), intent(in) :: plev
      integer,         intent(in) :: i, j
      real(kind=dblp)             :: tlev, r

      r    = log(plev/65000.e0_dblp)*rlogtl12
      tlev = temp4g(i,j) + r*(temp2g(i,j)-temp4g(i,j))

      end function ec_levtempgp


!------------------------------------------------------------------------------
!  FUNCTION ec_expint
!> @brief Exponential integral E_n(x)
!> @param[in]  n  order (>= 0)
!> @param[in]  x  argument (>= 0)
!> @return     E_n(x)
!------------------------------------------------------------------------------
      function ec_expint(n,x) result(res)

      use newunit_mod, only: error_id
      use error0_mod, only: ec_error

      implicit none

      integer,         intent(in) :: n
      real(kind=dblp), intent(in) :: x
      real(kind=dblp)             :: res

      integer,         parameter :: maxit = 100
      real(kind=dblp), parameter :: eps   = 1.e-10_dblp
      real(kind=dblp), parameter :: fpmin = 1.e-30_dblp
      real(kind=dblp), parameter :: euler = 0.5772156649_dblp

      integer         :: i, ii, nm1
      real(kind=dblp) :: a, b, c, d, del, fact, h, psi

      nm1 = n - 1

      if (n.lt.0 .or. x.lt.0. .or. &
          (x.eq.0. .and. (n.eq.0 .or. n.eq.1))) then
        write(error_id,*) 'bad arguments in ec_expint'
        call ec_error(20)

      else if (n.eq.0) then
        res = exp(-x)/x

      else if (x.eq.0.) then
        res = 1./nm1

      else if (x.gt.1.) then
        b = x + n
        c = 1./fpmin
        d = 1./b
        h = d
        do i = 1, maxit
          a   = -i*(nm1+i)
          b   = b + 2.
          d   = 1./(a*d+b)
          c   = b + a/c
          del = c*d
          h   = h*del
          if (abs(del-1.) .lt. eps) then
            res = h*exp(-x)
            return
          endif
        enddo
        call ec_error(20)

      else
        if (nm1.ne.0) then
          res = 1./nm1
        else
          res = -log(x) - euler
        endif
        fact = 1.
        do i = 1, maxit
          fact = -fact*x/i
          if (i.ne.nm1) then
            del = -fact/(i-nm1)
          else
            psi = -euler
            do ii = 1, nm1
              psi = psi + 1./ii
            enddo
            del = fact*(-log(x)+psi)
          endif
          res = res + del
          if (abs(del) .lt. abs(res)*eps) return
        enddo
        call ec_error(20)
      endif

      end function ec_expint



!------------------------------------------------------------------------------
!  FUNCTION ec_surfmois
!> @brief Computes surface specific humidity qsurfn from ground pressure
!>        and surface temperature via the saturation mixing ratio.
!>        Includes downscaling branch (#if DOWNSTS).
!> @param[in]  nn  surface type index
!> @return     .true. on success
!------------------------------------------------------------------------------
      function ec_surfmois(nn) result(returnValue)

      use comatm,          only: nlat, nlon
      use comsurf_mod,     only: qsurfn, pgroundn, tsurfn

#if ( DOWNSTS == 1 )
      use comsurf_mod,     only: nld
      use vertDownsc_mod,  only: pground_d, tsurfn_d, qsurfn_d
      use ecbilt_topography, only: nb_levls
#endif

      implicit none

      integer(kind=ip), intent(in) :: nn

      logical          :: returnValue
      integer(kind=ip) :: i, j
#if ( DOWNSTS == 1 )
      integer(kind=ip) :: nb_down
#endif

      do j = 1, nlon
        do i = 1, nlat
          qsurfn(i,j,nn) = ec_qsat(pgroundn(i,j,nn), tsurfn(i,j,nn))
#if ( DOWNSTS == 1 )
          if (nn.eq.nld) then
            do nb_down = 1, nb_levls
              qsurfn_d(i,j,nn,nb_down) = &
                ec_qsat(pground_d(i,j,nb_down), tsurfn_d(i,j,nn,nb_down))
            enddo
          endif
#endif
        enddo
      enddo

      returnValue = .true.

      end function ec_surfmois


end module atm_thermo_mod
