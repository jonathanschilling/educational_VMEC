!> \file
SUBROUTINE open_output_files (extension, lscreen, lfirst)
  USE safe_open_mod
  USE vparams, ONLY: nthreed, nthreed0
  IMPLICIT NONE

  CHARACTER(LEN=*) :: extension
  LOGICAL :: lscreen, lfirst

  INTEGER :: iread, inthreed=0
  CHARACTER(LEN=120) :: threed1_file

  ! OPEN FILES FOR READING, WRITING
  inthreed = 0
  threed1_file = 'threed1.'//extension

  INQUIRE(FILE=threed1_file, OPENED=lfirst)
  lfirst = .not.lfirst
  IF (.not.lfirst) RETURN

  IF (lscreen) WRITE (*, '(33('' -''))')
  nthreed = nthreed0
  CALL safe_open(nthreed, iread, threed1_file, 'new', 'formatted')
  IF (iread .ne. 0) THEN
     IF (lscreen) PRINT *,' VMEC OUTPUT FILES ALREADY EXIST: OVERWRITING THEM ...'
     CALL safe_open(nthreed, inthreed, threed1_file, 'replace', 'formatted')
  ENDIF

END SUBROUTINE open_output_files