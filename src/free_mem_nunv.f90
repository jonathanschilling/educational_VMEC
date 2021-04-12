!> \file
SUBROUTINE free_mem_nunv
  USE vmec_main
  USE vacmod, only: free_mem_nestor
  IMPLICIT NONE

  INTEGER :: istat1 = 0, istat2 = 0

  IF (ALLOCATED(bsubu0)) &
      DEALLOCATE (bsubu0, rbsq, dbsq, stat=istat1)
  IF (ALLOCATED(rmn_bdy)) &
      DEALLOCATE (rmn_bdy, zmn_bdy, stat=istat2)

  call free_mem_nestor
  if (allocated(bsqsav)) then
     deallocate(bsqsav)
  end if

  IF (istat1.ne.0 .or. istat2.ne.0) THEN
      PRINT *,' deallocation problem in free_mem_nunv'
      PRINT *,' istat1 = ',istat1
      PRINT *,' istat2 = ',istat2
  ENDIF

END SUBROUTINE free_mem_nunv
