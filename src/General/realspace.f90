!> \file
MODULE realspace
  USE stel_kinds
  IMPLICIT NONE

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: r1, ru, rv, zu, zv, rcon, zcon
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: z1
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: guu, guv, gvv,              &
     ru0, zu0, gcon, rcon0, zcon0, phip, chip, shalf, sqrts, wint
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: extra1, extra2, extra3, extra4

END MODULE realspace
