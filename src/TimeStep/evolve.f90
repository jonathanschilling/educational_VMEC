!> \file
SUBROUTINE evolve(time_step, ier_flag, liter_flag, lscreen)
  USE vmec_main
  USE vmec_params, ONLY: bad_jacobian_flag, successful_term_flag,       &
                         norm_term_flag
  USE xstuff
  USE timer_sub

  IMPLICIT NONE

  REAL(rprec), intent(in) :: time_step
  INTEGER, INTENT(inout) :: ier_flag
  LOGICAL, INTENT(inout) :: liter_flag
  LOGICAL, INTENT(in)  :: lscreen


  CHARACTER(LEN=*), PARAMETER :: fcn_message = "External calls to FUNCT3D: "
  REAL(rprec) :: fsq1, dtau, b1, fac

  ! COMPUTE MHD FORCES
  CALL funct3d (lscreen, ier_flag)

  ! COMPUTE ABSOLUTE STOPPING CRITERION
  IF (iter2.eq.1 .and. irst.eq.2) THEN
     ! first iteration, jacobian was not computed correctly
     ier_flag = bad_jacobian_flag
  ELSE IF (fsqr.le.ftolv .and. fsqz.le.ftolv .and. fsql.le.ftolv) THEN
     ! converged to desired tolerance
     liter_flag = .false.
     ier_flag = successful_term_flag
  ENDIF

  IF (ier_flag.ne.norm_term_flag .or. .not.liter_flag) then
     ! errornous iteration or shall not iterate further
     RETURN
  end if
  ! no error and not converged --> keep going...


  ! COMPUTE DAMPING PARAMETER (DTAU) AND
  ! EVOLVE R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE

  fsq1 = fsqr1 + fsqz1 + fsql1

  IF (iter2 .eq. iter1) then
     ! initialize all entries in otau to 0.15/time_step --> required for averaging
     ! otau: "over" tau --> 1/tau ???
     otau(:ndamp) = cp15/time_step
  end if


  IF (iter2 .gt. iter1 .and. fsq1 .ne. zero) then
     ! fsq is 1 (first iteration) or fsq1 from previous iteration

     ! fsq1/fsq is y_n assuming monotonic decrease of energy
     dtau = MIN(ABS(LOG(fsq1/fsq)), cp15)
  end if

  ! update backup copy of fsq1
  fsq = fsq1

  ! shift array for averaging to the left to make space at end for new entry
  otau(1:ndamp-1) = otau(2:ndamp)

  IF (iter2 .gt. iter1) then
     ! insert new 1/tau value at end of otau array
     otau(ndamp) = dtau/time_step
  end if

  ! averaging over ndamp entries : 1/ndamp*sum(otau)
  otav = SUM(otau(:ndamp))/ndamp

  !
  dtau = time_step*otav/2

  b1  = one - dtau
  fac = one/(one + dtau)

  ! THIS IS THE TIME-STEP ALGORITHM. IT IS ESSENTIALLY A CONJUGATE
  ! GRADIENT METHOD, WITHOUT THE LINE SEARCHES (FLETCHER-REEVES),
  ! BASED ON A METHOD GIVEN BY P. GARABEDIAN
  xcdot = fac*(b1*xcdot + time_step*gc) ! update velocity
  xc    = xc + time_step*xcdot          ! advance xc by velocity given in xcdot

END SUBROUTINE evolve