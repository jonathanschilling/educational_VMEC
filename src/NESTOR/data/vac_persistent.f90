!> \file
MODULE vac_persistent
  USE stel_kinds
  IMPLICIT NONE

  INTEGER,     DIMENSION(:),     ALLOCATABLE :: imirr

  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: sinper
  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: cosper
  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: sinuv
  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: cosuv
  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: tanu
  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: tanv
  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: xmpot
  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: xnpot
  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: csign

  REAL(rprec), DIMENSION(:,:),   ALLOCATABLE :: sinu
  REAL(rprec), DIMENSION(:,:),   ALLOCATABLE :: cosu
  REAL(rprec), DIMENSION(:,:),   ALLOCATABLE :: sinv
  REAL(rprec), DIMENSION(:,:),   ALLOCATABLE :: cosv
  REAL(rprec), DIMENSION(:,:),   ALLOCATABLE :: sinui
  REAL(rprec), DIMENSION(:,:),   ALLOCATABLE :: cosui
  REAL(rprec), DIMENSION(:,:),   ALLOCATABLE :: sinu1
  REAL(rprec), DIMENSION(:,:),   ALLOCATABLE :: cosu1
  REAL(rprec), DIMENSION(:,:),   ALLOCATABLE :: sinv1
  REAL(rprec), DIMENSION(:,:),   ALLOCATABLE :: cosv1

  REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: cmns

  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: bsubu_sur
  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: bsubv_sur
  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: bsupu_sur
  REAL(rprec), DIMENSION(:),     ALLOCATABLE :: bsupv_sur

END MODULE vac_persistent
