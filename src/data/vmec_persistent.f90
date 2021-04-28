!> \file
MODULE vmec_persistent

  USE stel_kinds, ONLY: rprec

  IMPLICIT NONE

  INTEGER, DIMENSION(:), ALLOCATABLE :: ixm
  INTEGER, DIMENSION(:), ALLOCATABLE :: jmin3

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmu
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: sinmu
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmum
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: sinmum
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmumi
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: sinmumi
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosnv
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: sinnv
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosnvn
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: sinnvn

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmui
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: sinmui
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmui3
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmumi3

  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: xm
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: xn
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: xm_nyq
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: xn_nyq

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: cos01
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: sin01

END MODULE vmec_persistent
