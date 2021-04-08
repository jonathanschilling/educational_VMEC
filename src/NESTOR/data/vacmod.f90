!> \file
MODULE vacmod
  USE vacmod0
  USE vac_persistent
  USE vparams, ONLY: zero, one, c2p0, cp5

  IMPLICIT NONE

  REAL(rprec), PARAMETER :: p5 = cp5
  REAL(rprec), PARAMETER :: two = c2p0

  INTEGER :: nfper
  INTEGER :: nvper

  REAL(rprec) :: bsubvvac
  REAL(rprec) :: pi2
  REAL(rprec) :: pi3
  REAL(rprec) :: pi4
  REAL(rprec) :: alp
  REAL(rprec) :: alu
  REAL(rprec) :: alv
  REAL(rprec) :: alvp
  REAL(rprec) :: onp
  REAL(rprec) :: onp2

  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: potvac

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bvecsav
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: amatsav

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bexni

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: brv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bphiv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bzv

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsqvac

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: r1b
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rub
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rvb
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: z1b
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: zub
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: zvb

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bexu
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bexv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bexn

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: auu
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: auv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: avv

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: snr
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: snv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: snz

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: drv

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: guu_b
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: guv_b
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: gvv_b

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rzb2

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rcosuv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rsinuv

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bredge
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bpedge
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bzedge

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: raxis_nestor
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: zaxis_nestor

END MODULE vacmod
