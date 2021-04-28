!> \file
!> \brief Open output files.

!> \brief Open output files.
!>
!> @param extension input file "extension": part after \c 'input.' .
!> @param lfirst flag to indicate if this is the first call to this routine or not
SUBROUTINE open_output_files (extension, lfirst)
  USE safe_open_mod
  USE vparams, ONLY: nthreed, nthreed0
  IMPLICIT NONE

  CHARACTER(LEN=*) :: extension
  LOGICAL :: lfirst

  INTEGER :: iread, inthreed=0
  CHARACTER(LEN=120) :: threed1_file

  ! OPEN FILES FOR READING, WRITING
  inthreed = 0
  threed1_file = 'threed1.'//extension

  INQUIRE(FILE=threed1_file, OPENED=lfirst)
  lfirst = .not.lfirst
  IF (.not.lfirst) RETURN

  WRITE (*, '(33('' -''))')
  nthreed = nthreed0
  CALL safe_open(nthreed, iread, threed1_file, 'new', 'formatted')
  IF (iread .ne. 0) THEN
     PRINT *,' VMEC OUTPUT FILES ALREADY EXIST: OVERWRITING THEM ...'
     CALL safe_open(nthreed, inthreed, threed1_file, 'replace', 'formatted')
  ENDIF

END SUBROUTINE open_output_files
