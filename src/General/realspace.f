!> \file
      MODULE realspace
      USE stel_kinds
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1   r1, ru, rv, zu, zv, rcon, zcon
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: z1
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: guu, guv, gvv, sigma_an,
     1   ru0, zu0, gcon, rcon0, zcon0, phip, chip, shalf, sqrts, wint
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET ::
     1   extra1, extra2, extra3, extra4
C-----------------------------------------------
      END MODULE realspace
