!> \file
MODULE vmec_persistent
  USE stel_kinds, ONLY: rprec
  IMPLICIT NONE

  INTEGER, DIMENSION(:), ALLOCATABLE :: ixm
  INTEGER, DIMENSION(:), ALLOCATABLE :: jmin3
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmu, sinmu
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmum, sinmum
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmumi, sinmumi
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosnv, sinnv
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosnvn, sinnvn

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmui, sinmui
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: cosmui3, cosmumi3
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: xm, xn, xm_nyq, xn_nyq
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: cos01, sin01

END MODULE vmec_persistent
