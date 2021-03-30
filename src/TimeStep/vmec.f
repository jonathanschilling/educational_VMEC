      PROGRAM vmec
      USE vmec_input
      USE safe_open_mod
      USE vparams, ONLY: nthreed
      USE vmec_params, ONLY: bad_jacobian_flag,
     2    restart_flag, readin_flag, timestep_flag,
     3    output_flag, cleanup_flag,
     4    norm_term_flag, successful_term_flag
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
      INTEGER :: ictrl(5)
      CHARACTER(LEN=120) :: input_file, seq_ext, reset_file_name, arg
      CHARACTER(LEN=120), DIMENSION(10) :: command_arg
      LOGICAL :: lscreen
!     Local variables
!
!     ictrl:   array(5) of control variables for running "runvmec" routine
!              see "runvmec" for a description
!

!
!     Read in command-line arguments to get input file or sequence file,
!     screen display information, and restart information
!
      INTERFACE
         SUBROUTINE runvmec(ictrl_array, input_file0,
     1                      lscreen, reset_file_name)
         IMPLICIT NONE
         INTEGER, INTENT(inout), TARGET :: ictrl_array(5)
         LOGICAL, INTENT(in) :: lscreen
         CHARACTER(LEN=*), INTENT(in) :: input_file0
         CHARACTER(LEN=*), OPTIONAL :: reset_file_name
         END SUBROUTINE runvmec
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

      CALL runvmec (ictrl, input_file, lscreen, reset_file_name)

      ierr_vmec = ictrl(2)
      SELECT CASE (ierr_vmec)
      CASE (bad_jacobian_flag)                             !Bad jacobian even after axis reset and ns->3
         IF (lscreen) WRITE (6, '(/,1x,a)') bad_jacobian
         WRITE (nthreed, '(/,1x,a)') bad_jacobian
      CASE DEFAULT
      END SELECT

      END PROGRAM vmec
