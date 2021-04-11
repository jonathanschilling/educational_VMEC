!> \file
SUBROUTINE tolicu(torcur)
  USE vparams, ONLY: mu0
  USE vacmod
  USE biotsavart
  IMPLICIT NONE

  REAL(rprec), INTENT(in) :: torcur

  INTEGER :: i, kper, kv
  REAL(rprec) :: current(1)

  ! COMPUTE WIRE SEGMENTS (DO NOT CLOSE LOOP, CLOSURE DONE IN biotsavart ROUTINES)

  current = torcur/mu0

  i = 1
  DO kper = 1, nvper
     DO kv = 1, nv
        xpts(1,i) = raxis_nestor(kv)*(cosper(kper)*cosuv(kv) - sinper(kper)*sinuv(kv))
        xpts(2,i) = raxis_nestor(kv)*(sinper(kper)*cosuv(kv) + cosper(kper)*sinuv(kv))
        xpts(3,i) = zaxis_nestor(kv)
        i = i + 1
     END DO
  END DO

  ! INITIALIZE COIL-RELATED QUANTITIES
  CALL initialize_biotsavart (current, xpt=xpts)

END SUBROUTINE tolicu
