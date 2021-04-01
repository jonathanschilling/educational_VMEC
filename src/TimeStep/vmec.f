!> \file
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
      CHARACTER(LEN=120) :: input_file, seq_ext, arg
      CHARACTER(LEN=120) :: input_file0
      CHARACTER(LEN=120), DIMENSION(10) :: command_arg
      LOGICAL :: lscreen, lreset
      INTEGER :: ictrl_flag
      INTEGER :: ns_min, nsval, ns_old=0
      INTEGER :: igrid,
     1           jacob_off, niter_store
      INTEGER, SAVE :: igrid0
      LOGICAL, PARAMETER :: lreset_xc = .false.

      print *, "vmec"


!     Read in command-line arguments to get input file or sequence file,
!     screen display information, and restart information
!
      CALL getcarg(1, command_arg(1), numargs)
      DO iseq = 2, numargs
         CALL getcarg(iseq, command_arg(iseq), numargs)
      END DO

      ! default: enable screen output
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
         PRINT *,'  xvmec <filename> [noscreen]'
         PRINT *
         PRINT *,' noscreen: supresses all output to screen ',
     1           ' (default, or "screen", displays output)'

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

      ! Sets all flags --> full feature set enabled!
      ictrl_flag =  restart_flag+readin_flag +timestep_flag
     1           + output_flag +cleanup_flag

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

      ! consistency check on requested number of flux surfaces
      IF (ALL(ns_array.eq.0)) THEN
         ier_flag = ns_error_flag ! 'NS ARRAY MUST NOT BE ALL ZEROES'
         GOTO 1000
      END IF

      ! jacob_off=1 indicates that an initial run with ns=3 shall be inserted
      ! before the user-provided ns values from ns_array are processed
      ! in the multi-grid run
      jacob_off = 0

  50  CONTINUE

      ! convergence flag: initially not converged yet
      iequi = 0

      ! jacob_off=1 indicates that in the previous interation (got back here by GOTO 50)
      ! jacobian was bad --> also need to
      ! restart vacuum calculations
      IF (lfreeb .and. jacob_off.eq.1) ivac = 1

      ns_min = 3

      ! multi-grid iterations: loop over ns_array
      ! jacob_off=0,1 is required to insert one ns=3 run before
      ! starting to work with the user-provided ns_array
      ! if the first ns value from ns_array gave a bad jacobian
      ITERATIONS: DO igrid = igrid0-jacob_off, multi_ns_grid


         IF (igrid .lt. igrid0) THEN
!           TRY TO GET NON-SINGULAR JACOBIAN ON A 3 PT RADIAL MESH
            nsval = 3; ivac = -1
            ftolv = 1.e-4_dp
         ELSE
            ! proceed regularly with ns values from ns_array
            nsval = ns_array(igrid)
            IF (nsval .lt. ns_min) then
               ! skip entries that have less than ns_min flux surfaces
               CYCLE
            end if

            ! update ns_min --> reduction in number of flux surfaces not allowed
            ns_min = nsval

            ftolv = ftol_array(igrid)
            niter = niter_array(igrid)
         END IF

         IF (ns_old .le. nsval) then
            ! initialize ns-dependent arrays
            ! and (if previous solution is available) interpolate to current ns value
            CALL initialize_radial(nsval, ns_old, delt0r,
     2                             lscreen)
         end if

         ! *HERE* is the *ACTUAL* call to the equilibrium solver !
         CALL eqsolve (ier_flag, lscreen)

         ! break the multi-grid sequence if current number of flux surfaced
         ! did not reach convergence
         IF (ier_flag.ne.norm_term_flag .and.
     1       ier_flag.ne.successful_term_flag) EXIT

      END DO ITERATIONS

      ! if did not converge only because jacobian was bad
      ! and the intermediate ns=3 run was not performed yet (jacob_off is still == 0),
      ! retry the whole thing again
      IF (ier_flag.eq.bad_jacobian_flag .and. jacob_off.eq.0) THEN
         jacob_off = 1
         GO TO 50
      END IF

      CALL second0 (timeoff)
      timer(tsum) = timer(tsum) + timeoff - timeon

      ! write output files
 1000 CALL fileout (0, ictrl_flag, ier_flag, lscreen)

      ! free memory
      CALL free_mem_funct3d
      CALL free_mem_ns (lreset_xc)
      CALL free_mem_nunv
      CALL free_persistent_mem

      ! close output files
      CALL close_all_files

      ! verbose error message
      SELECT CASE (ier_flag)
      CASE (bad_jacobian_flag)
         ! Bad jacobian even after axis reset and ns->3
         IF (lscreen) WRITE (6, '(/,1x,a)') bad_jacobian
         WRITE (nthreed, '(/,1x,a)') bad_jacobian
      CASE DEFAULT
      END SELECT

      END PROGRAM vmec
