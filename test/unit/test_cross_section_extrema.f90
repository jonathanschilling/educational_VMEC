!> \file
!> \brief Unit test for the cross-section extrema of the geometry table in eqfor.
!>
!> An ellipse \f$R = R_0 + a \cos\theta\f$, \f$Z = Z_0 + \kappa a \sin\theta\f$
!> has its Z extrema at \f$Z_0 \pm \kappa a\f$, its R extrema at
!> \f$R_0 \pm a\f$ and \f$R = R_0\f$ at the largest Z, wherever \f$Z_0\f$ is.
!> The ellipse is stored as a lasym run stores it, centered at \f$Z_0 = 0\f$
!> and at \f$Z_0 = -2\f$, and as a stellarator-symmetric run stores it,
!> centered at \f$Z_0 = 0\f$.
program test_cross_section_extrema
  use stel_kinds, only: rprec
  use stel_constants, only: pi
  use cross_section, only: cross_section_extrema
  implicit none

  integer, parameter :: ns = 1, nzeta = 1, ntheta3 = 16
  integer, parameter :: ntheta2 = ntheta3/2 + 1
  real(rprec), parameter :: r0 = 3.0_rprec, a = 1.0_rprec, kappa = 1.5_rprec
  real(rprec), parameter :: tol = 1.0e-14_rprec

  real(rprec) :: r1(ns*nzeta*ntheta3, 0:1), z1(ns*nzeta*ntheta3, 0:1)
  real(rprec) :: z0(3), theta, xmin, xmax, zmin, zmax, rzmax, err, worst
  logical :: lasym(3)
  integer :: i, k

  z0 = (/ 0.0_rprec, -2.0_rprec, 0.0_rprec /)
  lasym = (/ .true., .true., .false. /)

  worst = 0.0_rprec
  do k = 1, 3
    r1 = 0.0_rprec
    z1 = 0.0_rprec
    do i = 1, ntheta3
      theta = 2.0_rprec*pi*(i - 1)/ntheta3
      r1(i, 0) = r0 + a*cos(theta)
      z1(i, 0) = z0(k) + kappa*a*sin(theta)
    end do

    call cross_section_extrema(r1, z1, 1.0_rprec, 1, 1, ns, nzeta, ntheta2, &
                               ntheta3, lasym(k), xmin, xmax, zmin, zmax, rzmax)

    err = max(abs(xmin - (r0 - a)), abs(xmax - (r0 + a)),          &
              abs(zmin - (z0(k) - kappa*a)), abs(zmax - (z0(k) + kappa*a)), &
              abs(rzmax - r0))
    write (*, '(A, L2, A, F5.1, 5(A, F8.4))') 'lasym', lasym(k), ', Z0 =', &
         z0(k), ': xmin', xmin, ', xmax', xmax, ', zmin', zmin, ', zmax', &
         zmax, ', rzmax', rzmax
    worst = max(worst, err)
  end do

  if (worst > tol) then
    write (*, '(A, ES12.4)') 'FAIL: extrema differ from the ellipse by ', worst
    error stop 1
  end if
  write (*, '(A)') 'PASS: cross-section extrema follow the ellipse'
end program test_cross_section_extrema
