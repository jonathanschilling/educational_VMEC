      PROGRAM vmec
      USE vmec_input
      USE safe_open_mod
      USE vparams, ONLY: nthreed
      USE vmec_params, ONLY: bad_jacobian_flag,
     2    restart_flag, readin_flag, timestep_flag,
     3    output_flag, cleanup_flag, ns_error_flag,
     4    norm_term_flag, successful_term_flag, reset_jacdt_flag
      USE vmec_main
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      CHARACTER(LEN=*), PARAMETER ::
     2    bad_jacobian = "The jacobian was non-definite!"
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: numargs, ierr_vmec, index_end,
     1   iopen, isnml, iread, iseq, index_seq,
     2   index_dat, iunit, ncount, nsteps, i
      INTEGER, TARGET :: ictrl(5)
      CHARACTER(LEN=120) :: input_file, seq_ext, reset_file_name, arg
      CHARACTER(LEN=120) :: input_file0
      CHARACTER(LEN=120), DIMENSION(10) :: command_arg
      LOGICAL :: lscreen
!     Local variables
!
!     ictrl:   array(5) of control variables for running "runvmec" routine
!              see "runvmec" for a description
!

      INTEGER, POINTER :: ier_flag
      INTEGER :: ictrl_flag, iseq_count
      INTEGER :: ns_index, ns_min, nsval, ns_old=0, numsteps
      INTEGER :: igrid,
     1           jacob_off, niter_store
      INTEGER, SAVE :: igrid0
      LOGICAL :: lreset


!
!     Read in command-line arguments to get input file or sequence file,
!     screen display information, and restart information
!
!       INTERFACE
!          SUBROUTINE runvmec(ictrl_array, input_file0,
!      1                      lscreen, reset_file_name)
!          IMPLICIT NONE
!          INTEGER, INTENT(inout), TARGET :: ictrl_array(5)
!          LOGICAL, INTENT(in) :: lscreen
!          CHARACTER(LEN=*), INTENT(in) :: input_file0
!          CHARACTER(LEN=*), OPTIONAL :: reset_file_name
!          END SUBROUTINE runvmec
!       END INTERFACE

      INTERFACE
         SUBROUTINE initialize_radial(nsval, ns_old, delt0,
     1                                lscreen, reset_file_name)
         USE vmec_main
         IMPLICIT NONE
         INTEGER, INTENT(in) :: nsval
         INTEGER, INTENT(inout) :: ns_old
         CHARACTER(LEN=*), OPTIONAL :: reset_file_name
         LOGICAL, INTENT(in) :: lscreen
         REAL(rprec), INTENT(out) :: delt0
         END SUBROUTINE initialize_radial
      END INTERFACE


      CALL getcarg(1, command_arg(1), numargs)
      DO iseq = 2, numargs
         CALL getcarg(iseq, command_arg(iseq), numargs)
      END DO

      lscreen = .true.
      reset_file_name = " "

      IF (numargs .lt. 1) THEN
         STOP 'Invalid command line'
      ELSE IF (command_arg(1).eq.'-h' .or. command_arg(1).eq.'/h') THEN
         PRINT *,
     1   ' ENTER INPUT FILE NAME OR INPUT-FILE SUFFIX ON COMMAND LINE'
         PRINT *
         PRINT *,' For example: '
         PRINT *,'    xvmec input.tftr OR xvmec tftr ',
     1           'OR xvmec ../input.tftr'
         PRINT *
         PRINT *,' Additional (optional) command arguments are',
     1           ' allowed:'
         PRINT *
         PRINT *,'  xvmec <filename> [noscreen] [reset=reset_wout_file]'
         PRINT *
         PRINT *,' noscreen: supresses all output to screen ',
     1           ' (default, or "screen", displays output)'
         PRINT *,' name of reset wout file (defaults to none)'

         STOP
      ELSE
         DO iseq = 2, numargs
            arg = command_arg(iseq)
            IF (TRIM(arg).eq.'noscreen' .or. TRIM(arg).eq.'NOSCREEN')
     1         lscreen = .false.
            index_end = INDEX(arg, "reset=")
            index_seq = MAX(INDEX(arg, "RESET="), index_end)
            IF (index_seq .gt. 0) reset_file_name = arg(index_seq+6:)
         END DO
      END IF

!
!     Determine type of file opened (sequential or input-data)
!     ARG1 (char var)
!          By DEFAULT, ARG1 obtained from the command
!          line is parsed as follows to determine the input data file(s):
!               a. Attempt to OPEN file ARG1 (full path + file name).
!                  Look for the VSEQ NAMELIST to obtain nseq, nseq_select, and
!                  extension array. If they exist and nseq>0, VMEC will run
!                  sequentially using input determined from the array EXTENSION[i]
!                  or input.EXTENSION[i]
!               b. If the command argument is not a sequence NAMELIST, THEN the data file
!                  ARG1 or input.ARG1 is READ directly, with NSEQ=1.
!
      arg = command_arg(1)
      index_dat = INDEX(arg,'.')
      index_end = LEN_TRIM(arg)
      IF (index_dat .gt. 0) THEN
         seq_ext  = arg(index_dat+1:index_end)
         input_file = TRIM(arg)
      ELSE
         seq_ext = TRIM(arg)
         input_file = 'input.'//TRIM(seq_ext)
      END IF

!
!     CALL EQUILIBRIUM SOLVER
!
      ictrl = 0
      ictrl(1) =   restart_flag+readin_flag +timestep_flag
     1           + output_flag +cleanup_flag                !Sets all flags



!
!
!       CALL runvmec (ictrl, input_file, lscreen, reset_file_name)
!
!
      input_file0 = input_file
      ictrl_flag =  ictrl(1)
      numsteps   =  ictrl(3)
      ier_flag   => ictrl(2)
      ns_index   =  ictrl(4)
      iseq_count =  ictrl(5)
      CALL second0 (timeon)
!
!     PARSE input_file into path/input.ext
!
      index_dat = INDEX(input_file0,'input.')
      index_end = LEN_TRIM(input_file0)
      IF (index_dat .gt. 0) THEN
         input_file = TRIM(input_file0)
         input_extension  = input_file0(index_dat+6:index_end)
      ELSE
         input_extension = input_file0(1:index_end)
         input_file = 'input.'//TRIM(input_extension)
      END IF

!
!     INITIALIZE PARAMETERS
!
      lreset = (IAND(ictrl_flag, restart_flag) .ne. 0)

      IF (lreset) CALL reset_params

      IF (IAND(ictrl_flag, reset_jacdt_flag) .ne. 0) THEN
         ijacob = 0
         delt0r = delt
      END IF

      IF (IAND(ictrl_flag, readin_flag) .ne. 0) THEN
!
!        READ INPUT FILE (INDATA NAMELIST), MGRID_FILE (VACUUM FIELD DATA)
!
         CALL vsetup (iseq_count)
         CALL readin (input_file, iseq_count, ier_flag, lscreen)
         IF (ier_flag .ne. 0) GOTO 1000
!
!        COMPUTE NS-INVARIANT ARRAYS
!
         CALL fixaray

      END IF

      IF (lreset) THEN
!
!        COMPUTE INITIAL SOLUTION ON COARSE GRID
!        IF PREVIOUS SEQUENCE DID NOT CONVERGE WELL
!
!        IF (lreseta) THEN    !NOTE: where externally, lreseta = T, set restart_flag bit
!                                    (ictrl_flag = IOR(ictrl_flag,restart_flag))
         igrid0 = 1
         ns_old = 0
         IF (LEN_TRIM(reset_file_name) .ne. 0) igrid0 = multi_ns_grid
         WRITE (nthreed, 30)
         delt0r = delt
      ENDIF

  30  FORMAT(' FSQR, FSQZ = Normalized Physical Force Residuals',/,
     1   ' fsqr, fsqz = Preconditioned Force Residuals',/,1x,23('-'),/,
     2   ' BEGIN FORCE ITERATIONS',/,1x,23('-'),/)

      IF (ALL(ns_array.eq.0) .and. ns_index.le.0) THEN
         ier_flag = ns_error_flag
         GOTO 1000
      END IF

      jacob_off = 0

      IF (IAND(ictrl_flag, timestep_flag) == 0) GOTO 1000

  50  CONTINUE
      iequi = 0
      IF (lfreeb .and. jacob_off.eq.1) ivac = 1    !!restart vacuum calculations

      ns_min = 3

      ITERATIONS: DO igrid = igrid0-jacob_off, multi_ns_grid
         IF (igrid .lt. igrid0) THEN
!           TRY TO GET NON-SINGULAR JACOBIAN ON A 3 PT RADIAL MESH
            nsval = 3; ivac = -1
            ftolv = 1.e-4_dp
         ELSE IF (ns_index .gt. 0) THEN
            IF (ns_index .gt. SIZE(ns_array)) THEN
               ier_flag = ns_error_flag
               RETURN
            END IF
            nsval = ns_array(ns_index)
            IF (nsval .le. 0) STOP 'NSVAL <= 0: WRONG INDEX VALUE'
            ftolv = ftol_array(ns_index)
            niter = niter_array(ns_index)
         ELSE
            nsval = ns_array(igrid)
            IF (nsval .lt. ns_min) CYCLE
            ns_min = nsval
            ictrl(4) = igrid
            ftolv = ftol_array(igrid)
            niter = niter_array(igrid)
         END IF

         IF (ns_old .le. nsval)
     1      CALL initialize_radial(nsval, ns_old, delt0r,
     2                             lscreen, reset_file_name)

!     CONTROL NUMBER OF STEPS
         IF (numsteps > 0) THEN
            niter_store = niter
            niter = numsteps+iter2-1
         END IF

         CALL eqsolve (ier_flag, lscreen)

         IF (numsteps > 0) THEN
            niter = niter_store
         END IF

         IF (ier_flag.ne.norm_term_flag .and.
     1       ier_flag.ne.successful_term_flag) EXIT
         IF (numsteps>0 .or. ns_index>0) EXIT

      END DO ITERATIONS

  100 CONTINUE

      IF (ier_flag.eq.bad_jacobian_flag .and. jacob_off.eq.0) THEN
         jacob_off = 1
         GO TO 50
      END IF

      CALL second0 (timeoff)
      timer(tsum) = timer(tsum) + timeoff - timeon

 1000 CALL fileout (iseq_count, ictrl_flag, ier_flag, lscreen)















      ierr_vmec = ictrl(2)
      SELECT CASE (ierr_vmec)
      CASE (bad_jacobian_flag)                             !Bad jacobian even after axis reset and ns->3
         IF (lscreen) WRITE (6, '(/,1x,a)') bad_jacobian
         WRITE (nthreed, '(/,1x,a)') bad_jacobian
      CASE DEFAULT
      END SELECT

      END PROGRAM vmec
