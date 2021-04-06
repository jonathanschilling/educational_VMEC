!> \file
SUBROUTINE vmec_getenv(ename, evalue)
  IMPLICIT NONE
  CHARACTER(LEN=*) :: ename, evalue

  CALL getenv(ename, evalue)

END SUBROUTINE vmec_getenv
