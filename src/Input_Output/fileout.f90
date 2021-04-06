!> \file
SUBROUTINE fileout(iseq, ictrl_flag, ier_flag, lscreen)
  USE vmec_main
  USE vac_persistent
  USE realspace
  USE vmec_params, ONLY: mscale, nscale, signgs, uminus,            &
        norm_term_flag, output_flag, cleanup_flag, successful_term_flag
  USE vforces
  USE xstuff, ONLY: xc, gc, xsave, scalxc
  USE timer_sub
  IMPLICIT NONE

  INTEGER, INTENT(in) :: iseq, ictrl_flag
  INTEGER, INTENT(inout) :: ier_flag
  LOGICAL :: lscreen

  INTEGER :: istat, loc_ier_flag

  INTEGER :: js, istat1=0, irst0
  REAL(rprec), DIMENSION(:), POINTER :: lu, lv
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
  CHARACTER(LEN=*), PARAMETER :: Warning = " Error deallocating global memory FILEOUT"
  LOGICAL :: log_open, lwrite, loutput, lterm

  lu => czmn
  lv => crmn

  ! COMPUTE REMAINING COVARIANT COMPONENT OF B (BSUBS),
  ! CYLINDRICAL COMPONENTS OF B (BR, BPHI, BZ), AND
  ! AVERAGE EQUILIBRIUM PROPERTIES AT END OF RUN
  iequi = 1
  lterm = (ier_flag.eq.norm_term_flag  .or. ier_flag.eq.successful_term_flag)
  lwrite = lterm
  loutput = (IAND(ictrl_flag, output_flag) .ne. 0)
  loc_ier_flag = ier_flag
  if (ier_flag .eq. successful_term_flag) then
     loc_ier_flag = norm_term_flag
  end if

  IF (lwrite .and. loutput) THEN

     ! Must save irst value if in "restart" mode
     irst0 = irst
     CALL funct3d (lscreen, istat)

     ! The sign of the jacobian MUST multiply phi to get the physically correct toroidal flux
     phi(1) = zero
     DO js = 2, ns
        phi(js) = phi(js-1) + phip(js)
     END DO
     phi = (signgs*twopi*hs)*phi

     irst = irst0

     CALL second0 (teqfon)
     ALLOCATE(br_out(nrzt), bz_out(nrzt), stat=istat)
     gc = xc
     CALL eqfor (br_out, bz_out, clmn, blmn, rcon(1,1), gc, ier_flag)
     CALL second0 (teqfoff)
     timer(teqf) = timer(teqf) + teqfoff - teqfon
  END IF

  ! Call WROUT to write output or error message if lwrite = false
  IF (loutput .and. ASSOCIATED(bzmn_o)) THEN
     CALL second0 (twouton)

     CALL wrout (bzmn_o, azmn_o, clmn, blmn, crmn_o, czmn_e,            &
                 crmn_e, xsave, gc, loc_ier_flag, lwrite)
     CALL second0 (twoutoff)

     timer(twout) = timer(twout) + twoutoff - twouton

     IF (lscreen)             PRINT 120, TRIM(werror(loc_ier_flag))
     IF (lscreen .and. lterm) PRINT 10,  TRIM(input_extension), ijacob

     IF (nthreed .gt. 0) THEN
        WRITE (nthreed,120) TRIM(werror(loc_ier_flag))
        IF (.not. lterm) GOTO 1000
        WRITE (nthreed, 10) TRIM(input_extension), ijacob
        CALL write_times(nthreed, lscreen, lfreeb, .false.)
     END IF
  END IF

 10 FORMAT(' FILE : ',a,/,' NUMBER OF JACOBIAN RESETS = ',i4,/)
120 FORMAT(/1x,a,/)


  IF (ALLOCATED(br_out)) THEN
     DEALLOCATE (br_out, bz_out)
  END IF

1000 CONTINUE

  ! DEALLOCATE GLOBAL MEMORY AND CLOSE FILES
  IF (IAND(ictrl_flag, cleanup_flag).eq.0) RETURN

  IF (ALLOCATED(cosmu)) then
     DEALLOCATE(cosmu, sinmu, cosmum, sinmum, cosmui, cosmumi,          &
                sinmui, sinmumi, cosnv, sinnv, cosnvn, sinnvn,          &
                cosmui3, cosmumi3, cos01, sin01, stat=istat1)
     IF (istat1 .ne. 0) PRINT *, Warning // "#1"
  end if

  IF (ALLOCATED(xm)) then
     DEALLOCATE (xm, xn, ixm, xm_nyq, xn_nyq,                           &
                 jmin3, mscale, nscale, uminus, stat=istat1)
     IF (istat1 .ne. 0) PRINT *, Warning // "#2"
  end if

  IF (ALLOCATED(tanu)) then
     DEALLOCATE(tanu, tanv, sinper, cosper, sinuv, cosuv, cmns,         &
                sinu, cosu, sinv, cosv, sinui, cosui, csign, sinu1,     &
                cosu1, sinv1, cosv1, imirr, xmpot, xnpot, stat=istat1)
     IF (istat1 .ne. 0) PRINT *, Warning // "#3"
  end if

END SUBROUTINE fileout
