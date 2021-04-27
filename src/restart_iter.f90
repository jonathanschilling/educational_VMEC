!> \file
SUBROUTINE restart_iter(time_step)
  USE vmec_main
  USE xstuff
  IMPLICIT NONE

  REAL(rprec), intent(inout) :: time_step

  REAL(rprec), PARAMETER :: c1p03 = 1.03_dp
  REAL(rprec), PARAMETER :: cp90  = 0.90_dp

  SELECT CASE (irst)
  CASE (2:3)

     ! restore previous good state
     xcdot(:neqs) = zero
     xc(:neqs) = xstore(:neqs)

     ! ---- reduce time step ----
     ! this is only executed when irst ==2 or ==3
     ! first case: irst == 2:
     !  => irst-2 == 0
     !  => 3-irst == 1
     !  => resulting operation: time_step *= 0.9  (reduce time step)
     ! second case: irst == 3:
     !  => irst-2 == 1
     !  => 3-irst == 0
     !  => resulting operation: time_step /= 1.03 (reduce time step)
     time_step = time_step*(  (irst-2)/c1p03  &
                            + (3-irst)*cp90  )

     IF (irst .eq. 2) THEN
        ijacob = ijacob + 1
        iter1 = iter2
     END IF

     irst = 1

     RETURN
  CASE DEFAULT
     ! save current state vector
     xstore(:neqs) = xc(:neqs)
     RETURN
  END SELECT

END SUBROUTINE restart_iter
