!> \file
SUBROUTINE belicu(bx, by, bz, cos1, sin1, rp, zp)
  USE vacmod, vm_bz => bz
  USE biotsavart
  IMPLICIT NONE

  REAL(rprec), DIMENSION(nuv2), INTENT(in)  :: cos1, sin1, rp, zp
  REAL(rprec), DIMENSION(nuv2), INTENT(out) :: bx, by, bz

  INTEGER :: j
  REAL(rprec), DIMENSION(3) :: xpt, bvec

  DO j = 1,nuv2
     xpt(1) = rp(j) * cos1(j)
     xpt(2) = rp(j) * sin1(j)
     xpt(3) = zp(j)

     CALL bsc_b (single_coil, xpt, bvec)

     bx(j) = bvec(1)
     by(j) = bvec(2)
     bz(j) = bvec(3)
  END DO

  CALL cleanup_biotsavart

END SUBROUTINE belicu
