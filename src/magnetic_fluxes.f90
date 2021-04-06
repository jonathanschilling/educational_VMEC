!> \file
FUNCTION torflux_deriv (x)
  USE stel_kinds
  USE vmec_main, ONLY: zero
  USE vmec_input, ONLY: tf => aphi

  REAL(rprec), INTENT(IN) :: x !< radial flux variable (=TOROIDAL FLUX ONLY IF APHI=1)
  REAL(rprec) :: torflux_deriv
  REAL(rprec), EXTERNAL :: polflux_deriv, piota
  INTEGER     :: i

  ! TOKAMAK/STELLARATOR (default is tf(1) = 1)
  torflux_deriv = 0
  DO i = UBOUND(tf,1), LBOUND(tf,1), -1
     torflux_deriv = x*torflux_deriv + i*tf(i)
  END DO

END FUNCTION torflux_deriv


FUNCTION polflux_deriv (x)
  USE stel_kinds

  REAL(rprec), INTENT(IN) :: x !< radial flux variable (=TOROIDAL FLUX ONLY IF APHI=1)
  REAL(rprec) :: tf
  REAL(rprec) :: polflux_deriv !< polflux_deriv == d(chi)/dx = iota(TF(x)) * torflux_deriv(x)
  REAL(rprec), EXTERNAL :: torflux, torflux_deriv, piota

  ! TOKAMAK/STELLARATOR: dchi/ds = iota * dphi/ds
  ! piota is assumed to be a function of the TF(x) on input
  tf = torflux(x)
  tf = MIN(tf, 1.0_dp)
  polflux_deriv = piota(tf)*torflux_deriv(x)

END FUNCTION polflux_deriv


FUNCTION torflux (x)
  USE stel_kinds

  REAL(rprec), INTENT(IN)  :: x !< radial flux variable (=TOROIDAL FLUX ONLY IF APHI=1)
  REAL(rprec) :: torflux, h, xi
  REAL(rprec), EXTERNAL :: torflux_deriv
  INTEGER     :: i

  h = 1.E-2_dp*x
  torflux = 0
  DO i=1,101
     xi = (i-1)*h
     torflux = torflux + torflux_deriv(xi)
  END DO
  torflux = torflux-0.5_dp*(torflux_deriv(0._dp)+torflux_deriv(x))
  torflux = h*torflux

END FUNCTION torflux


FUNCTION polflux (x)
  USE stel_kinds

  REAL(rprec), INTENT(IN) :: x !< radial flux variable (=TOROIDAL FLUX ONLY IF APHI=1)
  REAL(rprec) :: polflux, h, xi
  REAL(rprec), EXTERNAL :: polflux_deriv
  INTEGER     :: i

  h = 1.E-2_dp*x
  polflux = 0
  DO i=1,101
     xi = (i-1)*h
     polflux = polflux + polflux_deriv(xi)
  END DO
  polflux = polflux-0.5_dp*(polflux_deriv(0._dp)+polflux_deriv(x))
  polflux = h*polflux

END FUNCTION polflux
