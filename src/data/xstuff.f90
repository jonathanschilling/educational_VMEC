!> \file
MODULE xstuff
  USE stel_kinds, ONLY: rprec
  IMPLICIT NONE

  !> stacked array of R, Z, Lambda Spectral force coefficients (see readin for stack order)
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: gc

  !> stacked array of scaled R, Z, Lambda Fourier coefficients (see readin for stack order)
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: xc

  !> "velocity": change of Fourier coefficients per time step
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: xcdot

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: xsave

  !> backup copy of last-known-good xc
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: xstore

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: scalxc

END MODULE xstuff
