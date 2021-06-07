!> \file
!> \brief Write the output files.

!> \brief Write the output files.
!>
!> @param ier_flag error flag
SUBROUTINE fileout(ier_flag)
  USE vmec_main
  USE vac_persistent
  USE realspace
  USE vmec_params, ONLY: signgs, norm_term_flag, successful_term_flag
  USE vforces
  USE xstuff, ONLY: xc, gc, xsave
  IMPLICIT NONE

  INTEGER, INTENT(inout) :: ier_flag

  INTEGER :: istat
  INTEGER :: loc_ier_flag
  INTEGER :: js
  INTEGER :: first0
  LOGICAL :: lterm

  REAL(rprec), ALLOCATABLE :: br_out(:), bz_out(:)

  CHARACTER(LEN=*), PARAMETER, DIMENSION(0:10) :: werror = (/       &
     'EXECUTION TERMINATED NORMALLY                            ',   &
     'INITIAL JACOBIAN CHANGED SIGN (IMPROVE INITIAL GUESS)    ',   &
     'FORCE RESIDUALS EXCEED FTOL: MORE ITERATIONS REQUIRED    ',   &
     'VMEC INDATA ERROR: NCURR.ne.1 but BLOAT.ne.1.            ',   &
     'MORE THAN 75 JACOBIAN ITERATIONS (DECREASE DELT)         ',   &
     'ERROR READING INPUT FILE OR NAMELIST                     ',   &
     'NEW AXIS GUESS STILL FAILED TO GIVE GOOD JACOBIAN        ',   &
     'PHIEDGE HAS WRONG SIGN IN VACUUM SUBROUTINE              ',   &
     'NS ARRAY MUST NOT BE ALL ZEROES                          ',   &
     'ERROR READING MGRID FILE                                 ',   &
     'VAC-VMEC I_TOR MISMATCH : BOUNDARY MAY ENCLOSE EXT. COIL ' /)

  ! COMPUTE REMAINING COVARIANT COMPONENT OF B (BSUBS),
  ! CYLINDRICAL COMPONENTS OF B (BR, BPHI, BZ), AND
  ! AVERAGE EQUILIBRIUM PROPERTIES AT END OF RUN
  iequi = 1
  lterm = (ier_flag.eq.norm_term_flag  .or. ier_flag.eq.successful_term_flag)

  loc_ier_flag = ier_flag
  if (ier_flag .eq. successful_term_flag) then
     loc_ier_flag = norm_term_flag
  end if

  IF (lterm) THEN
     ! Must save first value if in "restart" mode
     first0 = first
     CALL funct3d (istat)

     ! The sign of the jacobian MUST multiply phi to get the physically correct toroidal flux
     ! Below few lines compute the toroidal flux profile from phip by quadrature
     phi(1) = zero
     DO js = 2, ns
        phi(js) = phi(js-1) + phip(js)
     END DO
     phi = (signgs*twopi*hs)*phi

     first = first0

     ALLOCATE(br_out(nrzt), bz_out(nrzt), stat=istat)

     gc = xc
     CALL eqfor (br_out, bz_out, clmn, blmn, rcon(1,1), gc, ier_flag)

  END IF

  IF (ASSOCIATED(bzmn_o)) THEN

     ! Call WROUT to write output
     ! or
     ! error message if lterm = false
     CALL wrout (bzmn_o, azmn_o, clmn, blmn, crmn_o, czmn_e,            &
                 crmn_e, xsave, gc, loc_ier_flag, lterm)

     ! These are the last few lines that appear on screen / in the threed1 file
     PRINT 120, TRIM(werror(loc_ier_flag))
     IF (lterm) PRINT 10,  TRIM(input_extension), ijacob
     IF (nthreed .gt. 0) THEN
        WRITE (nthreed,120) TRIM(werror(loc_ier_flag))
        IF (lterm) then
           WRITE (nthreed, 10) TRIM(input_extension), ijacob
        end if
     END IF

  END IF

 10 FORMAT(' FILE : ',a,/,' NUMBER OF JACOBIAN RESETS = ',i4,/)
120 FORMAT(/1x,a,/)

  IF (ALLOCATED(br_out)) THEN
     DEALLOCATE (br_out, bz_out)
  END IF

END SUBROUTINE fileout
