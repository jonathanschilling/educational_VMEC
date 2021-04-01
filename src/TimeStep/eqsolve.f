!> \file
      SUBROUTINE eqsolve(ier_flag, lscreen)
      USE vmec_main
      USE vmec_params, ONLY: ntmax, ns4, jac75_flag, norm_term_flag,
     1                       bad_jacobian_flag,
     2                       successful_term_flag
      USE realspace
      USE xstuff
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ier_flag
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p98 = 0.98_dp, p96 = 0.96_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec) :: w0
      LOGICAL :: liter_flag, lreset_internal
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!        iequi   counter used to call -EQFOR- at end of run
!        ijacob  counter for number of times jacobian changes sign
!        irst    "counter" monitoring sign of jacobian;
!                resets R, Z, and Lambda when jacobian changes sign
!                and decreases time step
!        signgs  sign of Jacobian : must be =1 (right-handed) or =-1 (left-handed)

!        iterj   stores position in main iteration loop (j=1,2)
!        ivac    counts number of free-boundary iterations
!        ndamp   number of iterations over which damping is averaged
!        meven   parity selection label for even poloidal modes of R and Z
!        modd    parity selection label for odd poloidal modes of R and
!        gc      stacked array of        R, Z, Lambda Spectral force coefficients (see readin for stack order)
!        xc      stacked array of scaled R, Z, Lambda Fourier        coefficients


      print *, " eqsolve"

      liter_flag = iter2 .eq. 1 ! true at startup of program

 1000 CONTINUE

!
!     COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
!
   20 CONTINUE
!
!     RECOMPUTE INITIAL PROFILE, BUT WITH IMPROVED AXIS
!     OR
!     RESTART FROM INITIAL PROFILE, BUT WITH A SMALLER TIME-STEP
      IF (irst .EQ. 2) THEN
         xc = 0
         CALL profil3d (xc(1), xc(1+irzloff), lreset_internal, .FALSE.)
         irst = 1
         IF (liter_flag) then
            CALL restart_iter(delt0r)    !(OFF IN v8.50)
         end if
      END IF

      ! start normal iterations
      liter_flag = .true.

      ! reset error flag
      ier_flag = norm_term_flag

!
!     FORCE ITERATION LOOP
!
      iter_loop: DO WHILE (liter_flag)



         ! ADVANCE FOURIER AMPLITUDES OF R, Z, AND LAMBDA
         CALL evolve (delt0r, ier_flag, liter_flag, lscreen)




         ! check for bad jacobian and bad initial guess for axis
         IF (ijacob.eq.0 .and.
     1       (ier_flag.eq.bad_jacobian_flag .or. irst.eq.4) .and.
     1       ns.ge.3) THEN

            IF (lscreen) THEN
               IF (ier_flag .eq. bad_jacobian_flag) THEN
                  PRINT *, ' INITIAL JACOBIAN CHANGED SIGN!'
               END IF

               PRINT *, ' TRYING TO IMPROVE INITIAL MAGNETIC AXIS GUESS'
            END IF

            CALL guess_axis (r1, z1, ru0, zu0)
            lreset_internal = .true.
            ijacob = 1
            irst = 2
            GOTO 20
         ELSE IF (ier_flag.ne.norm_term_flag .and.
     1            ier_flag.ne.successful_term_flag) THEN
            ! if something went totally wrong even in this initial stuff,
            ! so do not continue
            RETURN
         ENDIF

         ! compute MHD energy
         w0 = wb + wp/(gamma - one)
!
!     ADDITIONAL STOPPING CRITERION (set liter_flag to FALSE)
!
         ! the blocks for ijacob=25 or 50 are equal up to the point
         ! that for 25, delt0r is reset to 0.98*delt (delt given by user)
         ! and  for 50, delt0r is reset to 0.96*delt (delt given by user)
         IF (ijacob .eq. 25) THEN
            ! jacobian changed sign 25 times: hmmm? :-/
            irst = 2
            CALL restart_iter(delt0r)
            delt0r = p98*delt
            IF (lscreen) PRINT 120, delt0r
            irst = 1
            GOTO 1000
         ELSE IF (ijacob .eq. 50) THEN
            ! jacobian changed sign 50 times: what the hell? :-S
            irst = 2
            CALL restart_iter(delt0r)
            delt0r = p96*delt
            IF (lscreen) PRINT 120, delt0r
            irst = 1
            GOTO 1000
         ELSE IF (ijacob .ge. 75) THEN
            ! jacobian changed sign at least 75 times: time to give up :-(
            ier_flag = jac75_flag ! 'MORE THAN 75 JACOBIAN ITERATIONS (DECREASE DELT)'
            liter_flag = .false.
         ELSE IF (iter2.ge.niter .and. liter_flag) THEN
            ! allowed number of iterations exceeded
            liter_flag = .false.
         ENDIF

  120 FORMAT(2x,'HAVING A CONVERGENCE PROBLEM: RESETTING DELT TO ',f8.3,
     1  /,2x,'If this does NOT resolve the problem, try changing ',
     2       '(decrease OR increase) the value of DELT')
!
!       TIME STEP CONTROL
!
         IF (iter2.eq.iter1 .or. res0.eq.-1) then
            res0 = fsq
         end if

         res0 = MIN(res0,fsq)

         IF (fsq.le.res0 .and. iter2-iter1.gt.10) THEN
            ! Store current state (irst=1)
            CALL restart_iter(delt0r)
         ELSE IF (fsq.gt.100*res0 .and. iter2.gt.iter1) THEN
            ! Residuals are growing in time, reduce time step
            irst = 2
         ELSE IF (iter2-iter1 .gt. ns4/2 .and.
     1            iter2       .gt. 2*ns4 .and.
     1            fsqr+fsqz   .gt. c1pm2       ) THEN ! 1.0e-2

            irst = 3
         ENDIF

         IF (irst .ne. 1) THEN
!           Retrieve previous good state
            CALL restart_iter(delt0r)
            iter1 = iter2
         ELSE
!           Increment time step and printout every nstep iterations
            IF (MOD(iter2, nstep) .eq. 0 .or. ! status report due
     1          iter2             .eq. 1 .or. ! first iteration
     1          .not.liter_flag) then         ! iterations cancelled already

               CALL printout(iter2, delt0r, w0, lscreen)
            end if

            ! count iterations
            iter2 = iter2 + 1
         ENDIF

         ! ivac gets set to 1 in vacuum() of NESTOR
         IF (ivac .eq. 1) THEN
            ! vacuum pressure turned on at iter2 iterations (here)
            IF (lscreen) PRINT 110, iter2
            WRITE (nthreed, 110) iter2
  110 FORMAT(/,2x,'VACUUM PRESSURE TURNED ON AT ',i4,' ITERATIONS'/)

            ivac = ivac + 1
         ENDIF



      END DO iter_loop

      WRITE (nthreed, 60) w0*twopi**2
   60 FORMAT(/,' MHD Energy = ',1p,e12.6)


      END SUBROUTINE eqsolve
