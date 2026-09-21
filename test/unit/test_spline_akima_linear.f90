!> \file
!> \brief Unit test for the Akima spline on data that is linear over part of its span.
!>
!> Akima's slope at a knot weights the two neighbouring secants by the change of
!> slope on the far side of each, so where consecutive knots lie on a line the
!> slope at them is the slope of that line and the interpolant reproduces the
!> line between them. The data below follow y = x up to x = 0.5 and then bend,
!> so on [0, 0.5] the interpolant must equal x and its integral from 0 must
!> equal x**2/2.
program test_spline_akima_linear
  use stel_kinds, only: rprec
  implicit none

  integer, parameter :: n = 5
  real(rprec), parameter :: tol = 1.0e-12_rprec

  real(rprec) :: knots(n)  = (/ 0.0_rprec, 0.25_rprec, 0.5_rprec, &
                                0.75_rprec, 1.0_rprec /)
  real(rprec) :: values(n) = (/ 0.0_rprec, 0.25_rprec, 0.5_rprec, &
                                0.9_rprec, 0.2_rprec /)
  real(rprec) :: xq(5) = (/ 0.1_rprec, 0.2_rprec, 0.3_rprec, &
                            0.4_rprec, 0.45_rprec /)

  integer :: i, iflag
  real(rprec) :: y, y_integral, max_deviation

  max_deviation = 0.0_rprec
  do i = 1, size(xq)
    call spline_akima(xq(i), y, knots, values, n, iflag)
    call spline_akima_int(xq(i), y_integral, knots, values, n, iflag)
    max_deviation = max(max_deviation, abs(y - xq(i)), &
                        abs(y_integral - 0.5_rprec*xq(i)**2))
  end do

  write (*, '(A, ES12.4)') 'max deviation from the line = ', max_deviation
  if (max_deviation > tol) then
    write (*, '(A)') 'FAIL: Akima interpolant does not reproduce the linear data'
    error stop 1
  end if
  write (*, '(A)') 'PASS: Akima interpolant reproduces the linear data'
end program test_spline_akima_linear
