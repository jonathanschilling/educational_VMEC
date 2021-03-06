!> \file
!> \brief allocate arrays depending on the number of Fourier coefficients \c nunv

!> \brief allocate arrays depending on the number of Fourier coefficients \c nunv
!>
SUBROUTINE allocate_nunv
  USE vmec_main
  USE vmec_params, ONLY: ntmax
  USE vacmod, only: allocate_nestor, nuv2
  IMPLICIT NONE

  INTEGER :: istat1

  CALL free_mem_nunv

  ALLOCATE (bsubu0(nznt), rbsq(nznt), dbsq(nznt), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #1 in allocate_nunv'

  ALLOCATE (rmn_bdy(0:ntor,0:mpol1,ntmax), zmn_bdy(0:ntor,0:mpol1,ntmax), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #2 in allocate_nunv'

  ! PERSISTENT ARRAYS (DURATION OF PROGRAM) for NESTOR
  IF (lfreeb) then
     call allocate_nestor

     ALLOCATE(bsqsav(nuv2,3))
  end if

END SUBROUTINE allocate_nunv
