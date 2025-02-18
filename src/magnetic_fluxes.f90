!> \file
!> \brief Compute toroidal and poloidal magnetic flux profiles.

!> \brief Compute the radial derivative of the enclosed toroidal magnetic flux.
!>
!> @param x evaluation location; radial flux variable (=TOROIDAL FLUX ONLY IF APHI=1)
FUNCTION torflux_deriv (x)
  USE stel_kinds
  USE vmec_input, ONLY: tf => aphi

  REAL(rprec), INTENT(IN) :: x
  REAL(rprec) :: torflux_deriv
  INTEGER     :: i

  ! This evaluates the aphi polynomial (possibly modified by the user in the input file).

  ! TOKAMAK/STELLARATOR (default is tf(1) = 1)
  torflux_deriv = 0
  DO i = UBOUND(tf,1), LBOUND(tf,1), -1
     torflux_deriv = x*torflux_deriv + i*tf(i)
  END DO

END FUNCTION torflux_deriv

!> \brief Compute the radial derivative of the enclosed poloidal magnetic flux.
!>
!> polflux_deriv == d(chi)/dx = iota(TF(x)) * torflux_deriv(x)
!>
!> @param x evaluation location; radial flux variable (=TOROIDAL FLUX ONLY IF APHI=1)
FUNCTION polflux_deriv (x)
  USE stel_kinds

  REAL(rprec), INTENT(IN) :: x
  REAL(rprec) :: tf
  REAL(rprec) :: polflux_deriv
  REAL(rprec), EXTERNAL :: torflux, torflux_deriv, piota

  ! TOKAMAK/STELLARATOR: dchi/ds = iota * dphi/ds
  ! piota is assumed to be a function of the TF(x) on input
  tf = torflux(x)
  tf = MIN(tf, 1.0_dp)
  polflux_deriv = piota(tf)*torflux_deriv(x)

  ! TODO: how does above code work if the iota profile is not specified,
  ! i.e., in constrained-toroidal-current mode?
  ! or is it simply not used...?

END FUNCTION polflux_deriv

!> \brief Compute the enclosed toroidal magnetic flux.
!>
!> @param x evaluation location; radial flux variable (=TOROIDAL FLUX ONLY IF APHI=1)
FUNCTION torflux (x)
  USE stel_kinds

  REAL(rprec), INTENT(IN)  :: x
  REAL(rprec) :: torflux, h, xi
  REAL(rprec), EXTERNAL :: torflux_deriv
  INTEGER     :: i

  ! some ad-hoc trapezoidal integration radially outward to given x
  ! with a fixed number (100) of quadrature nodes?
  h = 1.E-2_dp*x
  torflux = 0.0_dp
  DO i=1,101
     xi = (i-1)*h
     torflux = torflux + torflux_deriv(xi)
  END DO
  torflux = torflux-0.5_dp*(torflux_deriv(0._dp)+torflux_deriv(x))
  torflux = h*torflux

END FUNCTION torflux

!> \brief Compute the enclosed poloidal magnetic flux.
!>
!> @param x evaluation location; radial flux variable (=TOROIDAL FLUX ONLY IF APHI=1)
FUNCTION polflux (x)
  USE stel_kinds

  REAL(rprec), INTENT(IN) :: x
  REAL(rprec) :: polflux, h, xi
  REAL(rprec), EXTERNAL :: polflux_deriv
  INTEGER     :: i

  ! some ad-hoc trapezoidal integration radially outward to given x
  ! with a fixed number (100) of quadrature nodes?
  h = 1.E-2_dp*x
  polflux = 0.0_dp
  DO i=1,101
     xi = (i-1)*h
     polflux = polflux + polflux_deriv(xi)
  END DO
  polflux = polflux-0.5_dp*(polflux_deriv(0._dp)+polflux_deriv(x))
  polflux = h*polflux

END FUNCTION polflux
