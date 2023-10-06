!> \file
MODULE stel_constants

  USE stel_kinds, ONLY: rprec, dp

  implicit none

!----------------------------------------------------------------------
!  Mathematical constants
!----------------------------------------------------------------------

  REAL(dp), PARAMETER :: pi = 4.0_dp * datan(1.0_dp) ! 3.14159265358979323846264338328_dp
  REAL(dp), PARAMETER :: pio2 = 2.0_dp * datan(1.0_dp) ! pi/2
  REAL(dp), PARAMETER :: twopi = 8.0_dp * datan(1.0_dp) ! 2.0_dp*pi
  REAL(dp), PARAMETER :: sqrt2 = 1.41421356237309504880168872_dp
  REAL(dp), PARAMETER :: degree = twopi / 360
  REAL(dp), PARAMETER :: one = 1.0_dp
  REAL(dp), PARAMETER :: zero = 0.0_dp

!----------------------------------------------------------------------
!  Physical constants
!------------------------------------------------------------------

  REAL(dp), PARAMETER :: mu0 = 2.0e-7_dp * twopi ! * 1.0e-7_dp

END MODULE stel_constants
