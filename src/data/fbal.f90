!> \file
MODULE fbal

  USE stel_kinds, ONLY: dp

  implicit none

  REAL(dp), DIMENSION(:), ALLOCATABLE :: rzu_fac
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rru_fac
  REAL(dp), DIMENSION(:), ALLOCATABLE :: frcc_fac
  REAL(dp), DIMENSION(:), ALLOCATABLE :: fzsc_fac

END MODULE fbal
