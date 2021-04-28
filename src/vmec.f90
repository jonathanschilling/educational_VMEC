!> \file
!> \brief Main program of VMEC

!> \brief Main program of VMEC
PROGRAM xvmec

  call run_vmec

end program xvmec

!> \brief Main driver for VMEC
subroutine vmec

  USE vmec_input
  USE safe_open_mod
  USE vparams, ONLY: nthreed
  USE vmec_params
  USE vmec_main
  use mgrid_mod, only: free_mgrid

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: bad_jacobian = "The jacobian was non-definite!"
  CHARACTER(LEN=*), PARAMETER :: Warning      = "Error deallocating global memory!"

  CHARACTER(LEN=120), DIMENSION(10) :: command_arg
  CHARACTER(LEN=120) :: arg
  CHARACTER(LEN=120) :: input_file
  INTEGER :: numargs
  INTEGER :: ier_flag
  INTEGER :: index_end
  INTEGER :: iseq
  INTEGER :: index_dat
  INTEGER :: istat1
  INTEGER :: ns_min
  INTEGER :: nsval
  INTEGER :: ns_old=0
  INTEGER :: igrid
  INTEGER :: jacob_off

  ! Read in command-line arguments to get input file or sequence file,
  ! screen display information, and restart information
  numargs = iargc()
  DO iseq = 1, numargs
     CALL getarg(iseq, command_arg(iseq))
  END DO

  IF (numargs .lt. 1) THEN
     STOP 'Invalid command line'
  ELSE IF (command_arg(1).eq.'-h' .or. command_arg(1).eq.'/h') THEN
     ! print help text
     PRINT *,' ENTER INPUT FILE NAME OR INPUT-FILE SUFFIX ON COMMAND LINE'
     PRINT *
     PRINT *,' For example: '
     PRINT *,'    xvmec input.tftr OR xvmec tftr OR xvmec ../input.tftr'
     PRINT *
     PRINT *,' Additional (optional) command arguments are allowed:'
     PRINT *
     PRINT *,'  xvmec <filename> [noscreen]'
     PRINT *
     PRINT *,' noscreen: supresses all output to screen (default, or "screen", displays output)'

     STOP
  END IF

  ! PARSE input_file into path/input.ext
  arg = command_arg(1)
  index_dat = INDEX(arg, '.') ! find position of '.'
  index_end = LEN_TRIM(arg)
  IF (index_dat .gt. 0) THEN
     ! found '.' in arg
     input_extension = arg(index_dat+1:index_end) ! extension is rest of arg after first '.'
     input_file = TRIM(arg)
  ELSE
     ! no '.' found --> extension given directly
     input_extension = TRIM(arg)
     input_file = 'input.'//TRIM(input_extension)
  END IF




  ier_flag = norm_term_flag


  ! INITIALIZE PARAMETERS
  CALL reset_params ! no further calls

  ! READ INPUT FILE (INDATA NAMELIST), MGRID_FILE (VACUUM FIELD DATA)
  CALL readin (input_file, ier_flag) ! calls read_indata, heading, read_mgrid, allocate_nunv (,flip_theta if required)

  IF (ier_flag .eq. 0) then
     ! reading of input file was successful, so start computation

     ! COMPUTE NS-INVARIANT ARRAYS
     CALL fixaray

     ns_old = 0
     delt0r = delt

     WRITE (nthreed, 30)
30 FORMAT(' FSQR, FSQZ = Normalized Physical Force Residuals',/,        &
          ' fsqr, fsqz = Preconditioned Force Residuals',/,1x,23('-'),/,&
          ' BEGIN FORCE ITERATIONS',/,1x,23('-'),/)

     do jacob_off = 0, 1
        ! jacob_off=1 indicates that an initial run with ns=3 shall be inserted
        ! before the user-provided ns values from ns_array are processed
        ! in the multi-grid run

        ! convergence flag: initially not converged yet
        iequi = 0

        IF (lfreeb .and. jacob_off.eq.1) then
           ! jacob_off=1 indicates that in the previous iteration jacobian was bad --> also need to restart vacuum calculations
           ivac = 1
        end if

        ns_min = 3

        ! multi-grid iterations: loop over ns_array
        ! jacob_off=0,1 is required to insert one ns=3 run before
        ! starting to work with the user-provided ns_array
        ! if the first ns value from ns_array gave a bad jacobian
        ITERATIONS: DO igrid = 1-jacob_off, multi_ns_grid

           IF (igrid .lt. 1) THEN
              ! igrid .lt. 1 can only happen when jacob_off == 1 (then igrid==0)

              ! TRY TO GET NON-SINGULAR JACOBIAN ON A 3 PT RADIAL MESH
              ! COMPUTE INITIAL SOLUTION ON COARSE GRID
              ! IF PREVIOUS SEQUENCE DID NOT CONVERGE WELL
              nsval = 3
              ftolv = 1.e-4_dp

              ! fully restart vacuum (why then assign ivac=1 then above???)
              ivac = -1

           ELSE
              ! proceed regularly with ns values from ns_array
              nsval = ns_array(igrid)
              IF (nsval .lt. ns_min) then
                 ! skip entries that have less flux surfaces than previous iteration
                 CYCLE
              end if

              ! update ns_min --> reduction in number of flux surfaces not allowed
              ns_min = nsval

              ftolv  = ftol_array(igrid)
              niterv = niter_array(igrid)
           END IF

           IF (ns_old .le. nsval) then
              ! initialize ns-dependent arrays
              ! and (if previous solution is available) interpolate to current ns value
              CALL initialize_radial(nsval, ns_old, delt0r)
           end if

           ! *HERE* is the *ACTUAL* call to the equilibrium solver !
           CALL eqsolve (ier_flag)

           ! break the multi-grid sequence if current number of flux surfaced
           ! did not reach convergence
           IF (ier_flag.ne.norm_term_flag .and. ier_flag.ne.successful_term_flag) then
              EXIT
           end if

        END DO ITERATIONS


        ! if did not converge only because jacobian was bad
        ! and the intermediate ns=3 run was not performed yet (jacob_off is still == 0),
        ! retry the whole thing again
        IF (ier_flag.ne.bad_jacobian_flag) THEN
           exit ! jacob_off loop
           ! otherwise, retry with initial ns=3 to fix bad jacobian
        END IF

     ! if ier_flag .eq. bad_jacobian_flag, repeat once again with ns=3 before
     end do ! jacob_off = 0, 1

  end if ! (ier_flag .eq. 0) after readin



  ! write output files
  CALL fileout (ier_flag)






  ! free memory
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



  CALL free_mem_funct3d
  CALL free_mem_ns
  CALL free_mem_nunv

  CALL free_mgrid (istat1)
  IF (istat1.ne.0) THEN
      PRINT *,'problem in free_mgrid'
      PRINT *,' istat1 = ',istat1
  ENDIF

  ! verbose error message
  SELECT CASE (ier_flag)
  CASE (bad_jacobian_flag)
     ! Bad jacobian even after axis reset and ns->3
     WRITE (      6, '(/,1x,a)') bad_jacobian
     WRITE (nthreed, '(/,1x,a)') bad_jacobian
  CASE DEFAULT
  END SELECT

    ! close threed1 file
  IF (nthreed .gt. 0) CLOSE (nthreed)

END subroutine vmec
