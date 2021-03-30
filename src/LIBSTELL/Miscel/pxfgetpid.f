      SUBROUTINE vmec_getpid(ipid, ierror)
      IMPLICIT NONE
      INTEGER :: ipid, ierror
C      INTEGER, EXTERNAL :: getpid
      INTEGER :: getpid

      ierror = 0

      ipid = getpid()

      IF (ipid < 0) ierror = -ipid
      
      END SUBROUTINE vmec_getpid
