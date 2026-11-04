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

  !> reverts without progress at which a multi-grid step goes over to the PARVMEC time-step control
  INTEGER, PARAMETER :: max_reverts_without_progress = 8

!   if (first .ne. 1) then
!     print *, "bad jacobian --> restart_iter (first = ", first, ")"
!   end if

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

     CALL count_revert

     RETURN
  CASE (3)

     ! restore previous good state
     xcdot(:neqs) = zero
     xc(:neqs) = xstore(:neqs)

     ! ---- reduce time step ----
     time_step = time_step / c1p03

     first = 1

     CALL count_revert

     RETURN
  CASE DEFAULT ! 1, 4
     ! save current state vector, e.g. first=1
     xstore(:neqs) = xc(:neqs)
     RETURN
  END SELECT

CONTAINS

  !> The VMEC 8.52 control reverts without limit to a stored state.
  !> That state can be worse than the one just reached.
  SUBROUTINE count_revert

    reverts_without_progress = reverts_without_progress + 1

    IF (.not.parvmec_control .and. reverts_without_progress.ge.max_reverts_without_progress) THEN
       parvmec_control = .true.
       res0 = -1
       res1 = -1
       PRINT 100, reverts_without_progress
    END IF

100 FORMAT(2x,'NO PROGRESS IN ',i0,' RESTARTS: CONTINUING WITH THE PARVMEC TIME-STEP CONTROL')

  END SUBROUTINE count_revert

END SUBROUTINE restart_iter
