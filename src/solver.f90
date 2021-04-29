!> \file
!> Solve a linear system of equations using \c dgesv

!> Solve a linear system of equations using \c dgesv
!>
!> @param amat coefficient matrix, size (m x m)
!> @param b right-hand side(s), size (m x nrhs)
!> @param m number of free parameters (=number of rows, cols in amat)
!> @param nrhs number of right-hand sides to solve simultaneously
!> @param info error flag
SUBROUTINE solver(amat, b, m, nrhs, info)

  USE stel_kinds

  IMPLICIT NONE

  INTEGER, INTENT(in)  :: m, nrhs
  INTEGER, INTENT(out) :: info
  REAL(rprec), INTENT(inout) :: amat(m,m)
  REAL(rprec), INTENT(inout) :: b(m,nrhs)

  INTEGER, ALLOCATABLE :: ipiv(:)

  info = 0
  ALLOCATE (ipiv(m))

  ! Compute the solution to a REAL system of linear equations
  !   AMAT * X = B,
  ! WHERE AMAT is an M-by-M matrix and X and B are N-by-NRHS matrices.

  ! FACTOR AMATRIX INTO LU FORM AND SOLVE BY GAUSSIAN ELIMINATION
  CALL dgesv (m, nrhs, amat, m, ipiv, b, m, info)

  ! IF (info .ne. 0) PRINT *, ' Condition No. = 0 in Solver'

  DEALLOCATE (ipiv)

END SUBROUTINE solver
