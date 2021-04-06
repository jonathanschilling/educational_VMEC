!> \file
MODULE vmercier
  USE vparams, ONLY: nsd, rprec, dp
  IMPLICIT NONE

  REAL(rprec), DIMENSION(nsd) :: Dshear
  REAL(rprec), DIMENSION(nsd) :: Dwell
  REAL(rprec), DIMENSION(nsd) :: Dcurr
  REAL(rprec), DIMENSION(nsd) :: Dmerc
  REAL(rprec), DIMENSION(nsd) :: Dgeod
END MODULE vmercier
