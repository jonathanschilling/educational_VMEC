!> \file
!> \brief Unit test for the value of the sum_atan profiles at s = 1.
!>
!> A term \f$(2/\pi) \arctan(b x^p / (1 - x)^q)\f$ of a sum_atan profile tends
!> to sign(b) for q > 0, to \f$(2/\pi) \arctan(b)\f$ for q = 0 and to zero for
!> q < 0 as x -> 1. pcurr and piota at x = 1 must return that limit; pcurr(1)
!> is the value profil1d scales the enclosed current by to match curtor.
program test_sum_atan
  use stel_kinds, only: rprec
  use stel_constants, only: pi
  use vmec_input, only: ac, ai, bloat, pcurr_type, piota_type
  implicit none

  real(rprec), external :: pcurr, piota
  real(rprec), parameter :: tol = 1.0e-14_rprec

  ! amplitude, b, p, q of a single term, one case per column
  real(rprec), parameter :: terms(4, 4) = reshape((/ &
       1.0_rprec,  2.0_rprec, 2.0_rprec,  1.0_rprec, &
       1.0_rprec, -2.0_rprec, 2.0_rprec,  1.0_rprec, &
       1.0_rprec,  2.0_rprec, 2.0_rprec,  0.0_rprec, &
       1.0_rprec,  2.0_rprec, 2.0_rprec, -1.0_rprec /), (/ 4, 4 /))

  real(rprec) :: expected(4), current_edge, iota_edge, worst
  integer :: k

  expected = (/ 1.0_rprec, -1.0_rprec, 2.0_rprec / pi * atan(2.0_rprec), &
                0.0_rprec /)

  bloat = 1.0_rprec
  pcurr_type = 'sum_atan'
  piota_type = 'sum_atan'

  worst = 0.0_rprec
  do k = 1, 4
    ac = 0.0_rprec
    ai = 0.0_rprec
    ac(1:4) = terms(:, k)
    ai(1:4) = terms(:, k)
    current_edge = pcurr(1.0_rprec)
    iota_edge = piota(1.0_rprec)
    write (*, '(A, I2, 3(A, ES14.6))') 'case', k, ': pcurr(1) = ', &
         current_edge, ', piota(1) = ', iota_edge, ', limit = ', expected(k)
    worst = max(worst, abs(current_edge - expected(k)), &
                abs(iota_edge - expected(k)))
  end do

  if (worst > tol) then
    write (*, '(A, ES12.4)') 'FAIL: s = 1 value differs from the limit by ', &
         worst
    error stop 1
  end if
  write (*, '(A)') 'PASS: sum_atan profiles take their limit at s = 1'
end program test_sum_atan
