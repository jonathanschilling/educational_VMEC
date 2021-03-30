      SUBROUTINE getcarg(narg, arg, numargs)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in)  :: narg
      INTEGER, INTENT(out) :: numargs
      CHARACTER(LEN=*), INTENT(out) :: arg
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: numchars
C-----------------------------------------------

      INTEGER iargc
      numargs = iargc()
      CALL getarg(narg, arg)

      END SUBROUTINE getcarg
