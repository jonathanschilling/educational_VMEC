!> \file
!> \brief Save current or restore previous good state vector and reduce time step.

!> \brief Save current or restore previous good state vector and reduce time step.
!>
!> @param time_step time step to be modified if convergence is bad
SUBROUTINE restart_iter(time_step)

  USE vmec_main
  USE xstuff

  IMPLICIT NONE

  REAL(rprec), intent(inout) :: time_step

  REAL(rprec), PARAMETER :: c1p03 = 1.03_dp
  REAL(rprec), PARAMETER :: cp90  = 0.90_dp

  SELECT CASE (first)
  CASE (2:3)

     ! restore previous good state
     xcdot(:neqs) = zero
     xc(:neqs) = xstore(:neqs)

     ! ---- reduce time step ----
     ! this is only executed when first ==2 or ==3
     ! first case: first == 2:
     !  => first-2 == 0
     !  => 3-first == 1
     !  => resulting operation: time_step *= 0.9  (reduce time step)
     ! second case: first == 3:
     !  => first-2 == 1
     !  => 3-first == 0
     !  => resulting operation: time_step /= 1.03 (reduce time step)
     time_step = time_step*(  (first-2)/c1p03  &
                            + (3-first)*cp90  )

     IF (first .eq. 2) THEN

        print *, "bad jacobian --> restart_iter"

        ijacob = ijacob + 1
        iter1 = iter2
     END IF

     first = 1

     RETURN
  CASE DEFAULT
     ! save current state vector, e.g. first=1
     xstore(:neqs) = xc(:neqs)
     RETURN
  END SELECT

END SUBROUTINE restart_iter
