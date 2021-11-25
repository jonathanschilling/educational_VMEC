!> \file
!> \brief Take a single time step in Fourier space to evolve the Fourier coefficients
!>        describing the equilibrium towards force balance.

!> \brief Take a single time step in Fourier space to evolve the Fourier coefficients
!>        describing the equilibrium towards force balance.
!>
!> @param time_step step length in parameter space to take
!> @param ier_flag error flag
!> @param liter_flag keep running?
SUBROUTINE evolve(time_step, ier_flag, liter_flag)
  USE vmec_main
  USE vmec_params, ONLY: bad_jacobian_flag, successful_term_flag,       &
                         norm_term_flag, ntmax
  USE xstuff

  use dbgout

  IMPLICIT NONE

  REAL(rprec), intent(in) :: time_step
  INTEGER, INTENT(inout) :: ier_flag
  LOGICAL, INTENT(inout) :: liter_flag

  integer :: i
  CHARACTER(LEN=*), PARAMETER :: fcn_message = "External calls to FUNCT3D: "
  REAL(rprec) :: fsq1, dtau, b1, fac

  ! COMPUTE MHD FORCES
  CALL funct3d (ier_flag)

  ! COMPUTE ABSOLUTE STOPPING CRITERION
  IF (iter2.eq.1 .and. first.eq.2) THEN
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
     ! dtau is temporarily re-used for something else here...
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

  dtau = time_step*otav/2

  b1  = one - dtau
  fac = one/(one + dtau)

  ! debugging output: xc, xcdot, gc before time step; xc and xcdot also after time step
  if (dump_evolve .and. iter2.le.nDump) then
    write(dump_filename, 999) ns, iter2, trim(input_extension)
999 format('evolve_',i5.5,'_',i6.6,'.',a,'.json')

    call open_dbg_out(trim(dump_filename))

    call add_real_5d("xc_before",    ns, ntor1, mpol, ntmax, 2, &
                      reshape(xc, (/ ns, ntor1, mpol, ntmax, 2 /), order=(/ 1, 3, 4, 5, 2 /) ) )
    call add_real_5d("xcdot_before", ns, ntor1, mpol, ntmax, 2, &
                   reshape(xcdot, (/ ns, ntor1, mpol, ntmax, 2 /), order=(/ 1, 3, 4, 5, 2 /) ) )
    call add_real_5d("gc",           ns, ntor1, mpol, ntmax, 2, &
                      reshape(gc, (/ ns, ntor1, mpol, ntmax, 2 /), order=(/ 1, 3, 4, 5, 2 /) ) )
  end if

  ! THIS IS THE TIME-STEP ALGORITHM. IT IS ESSENTIALLY A CONJUGATE
  ! GRADIENT METHOD, WITHOUT THE LINE SEARCHES (FLETCHER-REEVES),
  ! BASED ON A METHOD GIVEN BY P. GARABEDIAN
  xcdot = fac*(b1*xcdot + time_step*gc) ! update velocity
  xc    = xc + time_step*xcdot          ! advance xc by velocity given in xcdot

  if (dump_evolve .and. iter2.le.nDump) then
    call add_real_5d("xc_after",    ns, ntor1, mpol, ntmax, 2, &
                     reshape(xc, (/ ns, ntor1, mpol, ntmax, 2 /), order=(/ 1, 3, 4, 5, 2 /) ) )
    call add_real_5d("xcdot_after", ns, ntor1, mpol, ntmax, 2, &
                  reshape(xcdot, (/ ns, ntor1, mpol, ntmax, 2 /), order=(/ 1, 3, 4, 5, 2 /) ) )

    call close_dbg_out()
  end if

END SUBROUTINE evolve
