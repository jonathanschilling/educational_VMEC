!> \file
SUBROUTINE free_persistent_mem
  USE vmec_main
  USE xstuff
  USE mgrid_mod, ONLY: free_mgrid
  IMPLICIT NONE

  INTEGER :: istat1 = 0, istat2 = 0

  IF (ALLOCATED(xc)) DEALLOCATE (xc, scalxc, stat=istat1)
  CALL free_mgrid (istat2)

  IF (istat1.ne.0 .or. istat2.ne.0) THEN
      PRINT *,'problem in free_persistent_mem'
      PRINT *,' istat1 = ',istat1,' istat2 = ',istat2
  ENDIF

END SUBROUTINE free_persistent_mem
