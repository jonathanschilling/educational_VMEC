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
      INTEGER :: numargs, ier_flag, index_end,
     1   iopen, isnml, iread, iseq, index_seq,
     2   index_dat, iunit, ncount, nsteps, i
      INTEGER, TARGET :: ictrl(5)
      CHARACTER(LEN=120) :: input_file, seq_ext, arg
      CHARACTER(LEN=120) :: input_file0
      CHARACTER(LEN=120), DIMENSION(10) :: command_arg
      LOGICAL :: lscreen, lreset
      INTEGER :: ictrl_flag
      INTEGER :: ns_min, nsval, ns_old=0
      INTEGER :: igrid,
     1           jacob_off, niter_store
      INTEGER, SAVE :: igrid0


!     Read in command-line arguments to get input file or sequence file,
!     screen display information, and restart information
!
      CALL getcarg(1, command_arg(1), numargs)
      DO iseq = 2, numargs
         CALL getcarg(iseq, command_arg(iseq), numargs)
      END DO

      lscreen = .true.

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
         END DO
      END IF

!     PARSE input_file into path/input.ext
      arg = command_arg(1)
      index_dat = INDEX(arg,'.')
      index_end = LEN_TRIM(arg)
      IF (index_dat .gt. 0) THEN
         input_extension  = arg(index_dat+1:index_end)
         input_file = TRIM(arg)
      ELSE
         input_extension = TRIM(arg)
         input_file = 'input.'//TRIM(input_extension)
      END IF

      ictrl_flag =  restart_flag+readin_flag +timestep_flag
     1           + output_flag +cleanup_flag     ! Sets all flags

      CALL second0 (timeon)

!     INITIALIZE PARAMETERS
      CALL reset_params

!     READ INPUT FILE (INDATA NAMELIST), MGRID_FILE (VACUUM FIELD DATA)
      CALL vsetup
      CALL readin (input_file, ier_flag, lscreen)
      IF (ier_flag .ne. 0) GOTO 1000

!     COMPUTE NS-INVARIANT ARRAYS
      CALL fixaray

!     COMPUTE INITIAL SOLUTION ON COARSE GRID
!     IF PREVIOUS SEQUENCE DID NOT CONVERGE WELL
      igrid0 = 1
      ns_old = 0
      WRITE (nthreed, 30)
      delt0r = delt

  30  FORMAT(' FSQR, FSQZ = Normalized Physical Force Residuals',/,
     1   ' fsqr, fsqz = Preconditioned Force Residuals',/,1x,23('-'),/,
     2   ' BEGIN FORCE ITERATIONS',/,1x,23('-'),/)

      IF (ALL(ns_array.eq.0)) THEN
         ier_flag = ns_error_flag
         GOTO 1000
      END IF

      jacob_off = 0

  50  CONTINUE
      iequi = 0
      IF (lfreeb .and. jacob_off.eq.1) ivac = 1    !!restart vacuum calculations

      ns_min = 3

      ITERATIONS: DO igrid = igrid0-jacob_off, multi_ns_grid
         IF (igrid .lt. igrid0) THEN
!           TRY TO GET NON-SINGULAR JACOBIAN ON A 3 PT RADIAL MESH
            nsval = 3; ivac = -1
            ftolv = 1.e-4_dp
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
     2                             lscreen)

         ! *HERE* is the *ACTUAL* call to the equilibrium solver !
         CALL eqsolve (ier_flag, lscreen)

         IF (ier_flag.ne.norm_term_flag .and.
     1       ier_flag.ne.successful_term_flag) EXIT

      END DO ITERATIONS

      IF (ier_flag.eq.bad_jacobian_flag .and. jacob_off.eq.0) THEN
         jacob_off = 1
         GO TO 50
      END IF

      CALL second0 (timeoff)
      timer(tsum) = timer(tsum) + timeoff - timeon

 1000 CALL fileout (0, ictrl_flag, ier_flag, lscreen)

      SELECT CASE (ier_flag)
      CASE (bad_jacobian_flag)
         ! Bad jacobian even after axis reset and ns->3
         IF (lscreen) WRITE (6, '(/,1x,a)') bad_jacobian
         WRITE (nthreed, '(/,1x,a)') bad_jacobian
      CASE DEFAULT
      END SELECT

      END PROGRAM vmec
