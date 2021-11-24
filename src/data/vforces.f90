!> \file
MODULE vforces

  USE stel_kinds, ONLY: rprec

  IMPLICIT NONE

  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: armn
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: brmn
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: crmn

  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: azmn
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: bzmn
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: czmn

  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: blmn
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: clmn

  REAL(rprec), POINTER, DIMENSION(:) :: armn_e
  REAL(rprec), POINTER, DIMENSION(:) :: armn_o
  REAL(rprec), POINTER, DIMENSION(:) :: brmn_e
  REAL(rprec), POINTER, DIMENSION(:) :: brmn_o
  REAL(rprec), POINTER, DIMENSION(:) :: crmn_e
  REAL(rprec), POINTER, DIMENSION(:) :: crmn_o

  REAL(rprec), POINTER, DIMENSION(:) :: azmn_e
  REAL(rprec), POINTER, DIMENSION(:) :: azmn_o
  REAL(rprec), POINTER, DIMENSION(:) :: bzmn_e
  REAL(rprec), POINTER, DIMENSION(:) :: bzmn_o
  REAL(rprec), POINTER, DIMENSION(:) :: czmn_e
  REAL(rprec), POINTER, DIMENSION(:) :: czmn_o

  REAL(rprec), POINTER, DIMENSION(:) :: blmn_e
  REAL(rprec), POINTER, DIMENSION(:) :: blmn_o
  REAL(rprec), POINTER, DIMENSION(:) :: clmn_e
  REAL(rprec), POINTER, DIMENSION(:) :: clmn_o

  ! B-type forces with contributions only from constraint
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: brmn_con
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: bzmn_con
  REAL(rprec), POINTER, DIMENSION(:) :: brmn_e_con
  REAL(rprec), POINTER, DIMENSION(:) :: brmn_o_con
  REAL(rprec), POINTER, DIMENSION(:) :: bzmn_e_con
  REAL(rprec), POINTER, DIMENSION(:) :: bzmn_o_con

END MODULE vforces
