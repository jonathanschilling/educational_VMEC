!> \file
!> \brief Unit test for the Akima spline right-edge curvature.
!>
!> The corrected right edge mirrors the left edge, so interpolating the mirror
!> image of a dataset at mirrored query points reproduces the original
!> interpolant. This reflection symmetry holds only when the right-edge phantom
!> points use the right-edge curvature cr; the earlier code reused the left
!> curvature cl on the right edge, which broke the symmetry. The test therefore
!> fails against the pre-fix routine and passes against the fixed one.
program test_spline_akima
  use stel_kinds, only: rprec
  implicit none

  integer, parameter :: n = 5
  real(rprec), parameter :: tol = 1.0e-12_rprec

  real(rprec) :: knots(n)  = (/ 0.0_rprec, 0.2_rprec, 0.45_rprec, &
                                0.7_rprec, 1.0_rprec /)
  real(rprec) :: values(n) = (/ 0.3_rprec, 0.7_rprec, 0.55_rprec, &
                                0.9_rprec, 0.1_rprec /)
  real(rprec) :: mirrored_knots(n), mirrored_values(n)
  real(rprec) :: xq(5) = (/ 0.05_rprec, 0.18_rprec, 0.50_rprec, &
                            0.83_rprec, 0.96_rprec /)

  integer :: i, iflag
  real(rprec) :: forward, reflected, max_asymmetry

  ! mirror about the knot span [0, 1]: reverse the nodes and map x -> 1 - x
  do i = 1, n
    mirrored_knots(i)  = 1.0_rprec - knots(n + 1 - i)
    mirrored_values(i) = values(n + 1 - i)
  end do

  max_asymmetry = 0.0_rprec
  do i = 1, size(xq)
    call spline_akima(xq(i),            forward,   knots,          values,          n, iflag)
    call spline_akima(1.0_rprec-xq(i),  reflected, mirrored_knots, mirrored_values, n, iflag)
    max_asymmetry = max(max_asymmetry, abs(forward - reflected))
  end do

  write (*, '(A, ES12.4)') 'max reflection asymmetry = ', max_asymmetry
  if (max_asymmetry > tol) then
    write (*, '(A)') 'FAIL: Akima interpolant is not reflection-symmetric'
    error stop 1
  end if
  write (*, '(A)') 'PASS: Akima interpolant is reflection-symmetric'
end program test_spline_akima
