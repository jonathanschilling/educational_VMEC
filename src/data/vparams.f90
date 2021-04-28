!> \file
MODULE vparams

  USE stel_kinds
  USE stel_constants, ONLY: zero, twopi, mu0, one

  IMPLICIT NONE

  ! MAXIMUM PARAMETERS FOR VMEC CODE (FOR READING INPUT)
  ! USER SHOULD NOT ALTER THESE

  INTEGER, PARAMETER :: nsd = 10001 !< maximum number of radial nodes
  INTEGER, PARAMETER :: mpold = 101 !< maximum number of poloidal harmonics (in r,z,lam fourier series)
  INTEGER, PARAMETER :: ntord = 101 !< maximum number of toroidal harmonics
  INTEGER, PARAMETER :: ndatafmax  = 101
  INTEGER, PARAMETER :: nstore_seq = 100

  ! DERIVED (FROM FUNDAMENTAL) PARAMETERS FOR VMEC CODE
  INTEGER, PARAMETER :: mpol1d = mpold - 1
  INTEGER, PARAMETER :: ntor1d = ntord + 1

  ! file units
  INTEGER, PARAMETER :: nthreed0  = 9
  INTEGER, PARAMETER :: indata0   = nthreed0 + 2
  INTEGER, PARAMETER :: nwout0    = nthreed0 + 3
  INTEGER, PARAMETER :: jxbout0   = nthreed0 + 4
  INTEGER, PARAMETER :: nfort18   = 18
  INTEGER, PARAMETER :: nmercier0 = 52
  INTEGER            :: nthreed

  ! MISCELLANEOUS PARAMETERS
  REAL(rprec), PARAMETER :: c1pm2  =  1.e-2_dp
  REAL(rprec), PARAMETER :: cp15   =  0.15_dp
  REAL(rprec), PARAMETER :: cp25   = 0.25_dp
  REAL(rprec), PARAMETER :: cp5    = 0.50_dp
  REAL(rprec), PARAMETER :: c1pm8  = 1.0e-8_dp
  REAL(rprec), PARAMETER :: cbig   = 0.9e30_dp
  REAL(rprec), PARAMETER :: c2p0   = 2
  REAL(rprec), PARAMETER :: c3p0   = 3
  REAL(rprec), PARAMETER :: cp05   = 0.05_dp
  REAL(rprec), PARAMETER :: c1pm13 = 1.0e-13_dp
  REAL(rprec), PARAMETER :: osqrt2 = 0.707106781186547462_dp

END MODULE vparams
