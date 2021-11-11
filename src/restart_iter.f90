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

  if (first .ne. 1) then
    print *, "bad jacobian --> restart_iter (first = ", first, ")"
  end if
  
  SELECT CASE (first)
  case(2)
     ! restore previous good state
     xcdot(:neqs) = zero
     xc(:neqs) = xstore(:neqs)

     ! ---- reduce time step ----
     time_step = time_step * cp90

     ijacob = ijacob + 1
     iter1 = iter2

     first = 1

     RETURN
  CASE (3)

     ! restore previous good state
     xcdot(:neqs) = zero
     xc(:neqs) = xstore(:neqs)

     ! ---- reduce time step ----
     time_step = time_step / c1p03

     first = 1

     RETURN
  CASE DEFAULT ! 1, 4
     ! save current state vector, e.g. first=1
     xstore(:neqs) = xc(:neqs)
     RETURN
  END SELECT

END SUBROUTINE restart_iter
