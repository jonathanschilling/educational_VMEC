!> \file
SUBROUTINE free_mem_nunv
  USE vmec_main
  USE vacmod
  IMPLICIT NONE

  INTEGER :: istat1 = 0, istat2 = 0, istat3 = 0

  IF (ALLOCATED(bsubu0)) &
      DEALLOCATE (bsubu0, rbsq, dbsq, stat=istat1)
  IF (ALLOCATED(rmn_bdy)) &
      DEALLOCATE (rmn_bdy, zmn_bdy, stat=istat2)

  IF (ALLOCATED(amatsav)) then
      DEALLOCATE (amatsav, bvecsav, bsqsav, potvac,  &
                  raxis_nestor, zaxis_nestor, stat=istat3)
  end if

  IF (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0) THEN
      PRINT *,' deallocation problem in free_mem_nunv'
      PRINT *,' istat1 = ',istat1
      PRINT *,' istat2 = ',istat2
      PRINT *,' istat3 = ',istat3
  ENDIF

END SUBROUTINE free_mem_nunv
