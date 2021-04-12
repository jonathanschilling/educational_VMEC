!> \file
SUBROUTINE belicu(torcur, bx, by, bz, cos1, sin1, rp, zp)
  USE vacmod, vm_bz => bz
  USE biotsavart
  IMPLICIT NONE

  REAL(rprec), INTENT(in) :: torcur
  REAL(rprec), DIMENSION(nuv2), INTENT(in)  :: cos1, sin1, rp, zp
  REAL(rprec), DIMENSION(nuv2), INTENT(out) :: bx, by, bz

  INTEGER :: i, j, kper, kv
  REAL(rprec) :: current(1)
  REAL(rprec), DIMENSION(3) :: xpt, bvec

  ! net toroidal plasma current in A (?)
  current = torcur/mu0

  print *, "net toroidal current: ", current

  ! COMPUTE WIRE SEGMENTS (DO NOT CLOSE LOOP, CLOSURE DONE IN biotsavart ROUTINES)
  i = 1
  DO kper = 1, nvper
     DO kv = 1, nv
        ! xpts == xpt of _s_ource (current filament)
        xpts(1,i) = raxis_nestor(kv)*(cosper(kper)*cosuv(kv) - sinper(kper)*sinuv(kv))
        xpts(2,i) = raxis_nestor(kv)*(sinper(kper)*cosuv(kv) + cosper(kper)*sinuv(kv))
        xpts(3,i) = zaxis_nestor(kv)
        i = i + 1
     END DO
  END DO

  print *, "points along axis: ", i-1

  ! INITIALIZE COIL-RELATED QUANTITIES
  CALL initialize_biotsavart (current, xpt=xpts)

  ! evaluate magnetic field
  ! due to full net toroidal plasma current along magnetic axis
  ! at all grid points on the boundary
  DO j = 1,nuv2
     ! xpt == evaluation position
     xpt(1) = rp(j) * cos1(j)
     xpt(2) = rp(j) * sin1(j)
     xpt(3) = zp(j)

     CALL bsc_b (single_coil, xpt, bvec)

     bx(j) = bvec(1)
     by(j) = bvec(2)
     bz(j) = bvec(3)
  END DO

  print *, "eval at nuv2 points: ", nuv2

  ! cleanup internal state of biotsavart module
  CALL cleanup_biotsavart

END SUBROUTINE belicu
