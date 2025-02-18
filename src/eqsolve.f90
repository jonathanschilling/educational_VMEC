!> \file
!> \brief Iteratively evolve the Fourier coefficients that specify the equilibrium.

!> \brief Iteratively evolve the Fourier coefficients that specify the equilibrium.
!>
!> @param ier_flag error flag
SUBROUTINE eqsolve(ier_flag)

  USE vmec_main
  USE vmec_params, ONLY: ns4, jac75_flag, norm_term_flag,        &
                         bad_jacobian_flag, successful_term_flag
  USE realspace
  USE xstuff

  use dbgout, only: skip_dbgout_collison

  IMPLICIT NONE

  INTEGER, intent(inout) :: ier_flag

  REAL(rprec), PARAMETER :: p98 = 0.98_dp
  REAL(rprec), PARAMETER :: p96 = 0.96_dp

  REAL(rprec) :: w0
  LOGICAL :: liter_flag

  ! TODO: initial value of lreset_internal was undefined! --> now set to false
  LOGICAL :: lreset_internal

!   print *, " eqsolve"

  ! need to do this before jump label 20 to allow restart_iter to do its magic
  lreset_internal = .false.

  liter_flag = iter2 .eq. 1 ! true at start of current multi-grid iteration

  ! COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
20 CONTINUE ! try again

  ! !!! THIS must be the ONLY place where this gets incremented !!!
  num_eqsolve_retries = num_eqsolve_retries + 1

  ! print *, "goto 20, num_eqsolve_retries = ", num_eqsolve_retries

  ! RECOMPUTE INITIAL PROFILE, BUT WITH IMPROVED AXIS
  ! OR
  ! RESTART FROM INITIAL PROFILE, BUT WITH A SMALLER TIME-STEP
  IF (first .EQ. 2) THEN
     xc = 0.0_dp

     CALL profil3d (xc(1), xc(1+irzloff), lreset_internal)

     first = 1 ! tells restart_iter to store current xc in xstore
     IF (liter_flag) then
        ! Note that at this point, liter_flag could also simply contain (iter2 .eq. 1) (see above).
        ! (OFF IN v8.50)
        CALL restart_iter(delt0r)
     end if
  END IF

  ! start normal iterations
  liter_flag = .true.

  ! reset error flag
  ier_flag = norm_term_flag

  ! FORCE ITERATION LOOP
  iter_loop: DO WHILE (liter_flag)

     ! ADVANCE FOURIER AMPLITUDES OF R, Z, AND LAMBDA
     CALL evolve (delt0r, ier_flag, liter_flag)

     ! check for bad jacobian and bad initial guess for axis
     IF (ijacob.eq.0 .and.                                               &
         (ier_flag.eq.bad_jacobian_flag .or. first.eq.4) .and.           &
         ns.ge.3) THEN

        IF (ier_flag .eq. bad_jacobian_flag) THEN
           PRINT *, ' INITIAL JACOBIAN CHANGED SIGN!'
        END IF

        PRINT *, ' TRYING TO IMPROVE INITIAL MAGNETIC AXIS GUESS'

        CALL guess_axis (r1, z1, ru0, zu0)
        ijacob = 1

        ! prepare parameters to functions that get called due to lreset_internal and first=2
        lreset_internal = .true.
        first = 2
        GOTO 20 ! try again
     ELSE IF (ier_flag.ne.norm_term_flag .and.                          &
              ier_flag.ne.successful_term_flag) THEN
        ! if something went totally wrong even in this initial stuff,
        ! so do not continue
        RETURN
     ENDIF

     ! compute MHD energy
     w0 = wb + wp/(gamma - one)





     ! ADDITIONAL STOPPING CRITERION (set liter_flag to FALSE)

     ! the blocks for ijacob=25 or 50 are equal up to the point
     ! that for 25, delt0r is reset to 0.98*delt (delt given by user)
     ! and  for 50, delt0r is reset to 0.96*delt (delt given by user)
     IF (ijacob .eq. 25) THEN
        ! jacobian changed sign 25 times: hmmm? :-/
        first = 2
        CALL restart_iter(delt0r)
        delt0r = p98*delt
        PRINT 120, delt0r
        first = 1 ! done by restart_iter already...
        GOTO 20 ! try again
     ELSE IF (ijacob .eq. 50) THEN
        ! jacobian changed sign 50 times: what the hell? :-S
        first = 2
        CALL restart_iter(delt0r)
        delt0r = p96*delt
        PRINT 120, delt0r
        first = 1 ! done by restart_iter already...
        GOTO 20 ! try again
     ELSE IF (ijacob .ge. 75) THEN
        ! jacobian changed sign at least 75 times: time to give up :-(
        ier_flag = jac75_flag ! 'MORE THAN 75 JACOBIAN ITERATIONS (DECREASE DELT)'
        liter_flag = .false.
     ELSE IF (iter2.ge.niterv .and. liter_flag) THEN
        ! allowed number of iterations exceeded
        liter_flag = .false.
     ENDIF

  120 FORMAT(2x,'HAVING A CONVERGENCE PROBLEM: RESETTING DELT TO ',f8.3,&
        /,2x,'If this does NOT resolve the problem, try changing ',     &
             '(decrease OR increase) the value of DELT')





     ! TIME STEP CONTROL

     IF (iter2.eq.iter1 .or. res0.eq.-1) then
        ! if res0 has never been assigned (-1), give it the current value of fsq
        res0 = fsq
     end if

     ! res0 is the best force residual we got so far
     res0 = MIN(res0,fsq)

     IF (fsq.le.res0 .and. iter2-iter1.gt.10) THEN
        ! Store current state (first=1)
        ! --> was able to reduce force consistenly over at least 10 iterations
        CALL restart_iter(delt0r)
     ELSE IF (fsq.gt.100*res0 .and. iter2.gt.iter1) THEN
        ! Residuals are growing in time, reduce time step
        first = 2
     ELSE IF (iter2-iter1 .gt. ns4/2 .and.                              &
              iter2       .gt. 2*ns4 .and.                              &
              fsqr+fsqz   .gt. c1pm2       ) THEN ! 1.0e-2

        ! quite some iterations and quite large forces
        ! --> restart with different timestep

        ! TODO: maybe the threshold 0.01 is too large nowadays
        ! --> this could help fix the cases where VMEC gets stuck immediately at ~2e-3
        ! --> lower threshold, e.g. 1e-4 ?

        first = 3
     ENDIF





     IF (first .ne. 1) THEN
        ! Retrieve previous good state
        CALL restart_iter(delt0r)
        iter1 = iter2

        ! This breaks the dbgout logic, in that the same iter2
        ! and num_eqsolve_retries combination is touched twice.
        ! Hence, allow ONCE to overwrite the output file.
        skip_dbgout_collison = .true.
     ELSE
        ! Increment time step and printout every nstep iterations
        ! status report due or
        ! first iteration or
        ! iterations cancelled already (last iteration)
        IF (MOD(iter2, nstep) .eq. 0 .or.                               &
            iter2             .eq. 1 .or.                               &
            .not.liter_flag) then

           CALL printout(iter2, delt0r, w0)
        end if

        ! count iterations
        iter2 = iter2 + 1

        ! Disable again the temporary overwrite that was allowed above.
        ! Trust on VMEC to reach this branch in the next iteration.
        skip_dbgout_collison = .false.
     ENDIF




     ! ivac gets set to 1 in vacuum() of NESTOR
     IF (ivac .eq. 1) THEN
        ! vacuum pressure turned on at iter2 iterations (here)
        PRINT 110, iter2
        WRITE (nthreed, 110) iter2
110 FORMAT(/,2x,'VACUUM PRESSURE TURNED ON AT ',i4,' ITERATIONS'/)

        ! this makes ivac=1 --> ivac=2
        ivac = ivac + 1
     ENDIF

  END DO iter_loop

  ! write MHD energy at end of iterations for current number of surfaces
  WRITE (nthreed, 60) w0*twopi**2.0_dp
60 FORMAT(/,' MHD Energy = ',1p,e12.6)

END SUBROUTINE eqsolve
