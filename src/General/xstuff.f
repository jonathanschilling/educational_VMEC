!> \file
      MODULE xstuff
      USE stel_kinds, ONLY: rprec
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: gc         !< stacked array of        R, Z, Lambda Spectral force coefficients (see readin for stack order)
      REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: xc !< stacked array of scaled R, Z, Lambda Fourier        coefficients (see readin for stack order)

      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xcdot      !< "velocity": change of Fourier coefficients per time step
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xsave
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: xstore     !< backup copy of last-known-good xc
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: scalxc
C-----------------------------------------------
      END MODULE xstuff
