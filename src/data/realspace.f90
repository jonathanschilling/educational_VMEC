!> \file
MODULE realspace

  USE stel_kinds

  IMPLICIT NONE

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: r1
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ru
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: rv
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: z1
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: zu
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: zv
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: rcon
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: zcon

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: guu
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: guv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: gvv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: ru0
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: zu0
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: gcon
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rcon0
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: zcon0
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: phip  !< radial derivative of phi/(2*pi) on half-grid
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: chip  !< radial derivative of chi/(2*pi) on half-grid
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: shalf !< sqrt(s) ,two-dimensional array on half-grid
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: sqrts !< sqrt(s), two-dimensional array on full-grid
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: wint  !< two-dimensional array for normalizing angle integrations

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: extra1
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: extra2
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: extra3
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE, TARGET :: extra4

END MODULE realspace
