!> \file
      SUBROUTINE pxffork_g (ipid, ierror)
      IMPLICIT NONE
      INTEGER :: ipid, ierror

      ierror = 0
      IF (ipid < 0) ierror = -ipid
      
      END SUBROUTINE pxffork_g
