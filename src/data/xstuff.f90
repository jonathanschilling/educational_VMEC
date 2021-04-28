!> \file
MODULE xstuff

  USE stel_kinds, ONLY: rprec

  IMPLICIT NONE

! LOCAL VARIABLES
!
! rbcc,rbss,rbcs,rbsc
!          boundary Fourier coefficient arrays for R (of cosu*cosv, etc)
! zbcc,zbss,zbcs,zbsc
!          boundary Fourier coefficient arrays for Z
!
! XCC*COS(MU)COS(NV), XCS*COS(MU)SIN(NV), ETC
!
! STACKING ORDER DEPENDS ON LASYM AND LTHREED.
! EACH COMPONENT XCC, XSS, XSC, XCS HAS SIZE = mns.
!
!   LTHREED=F,      LTHREED=F,      LTHREED=T,      LTHREED=T
!   LASYM=F         LASYM=T         LASYM=F         LASYM=T
!
!    rmncc           rmncc           rmncc           rmncc
!    zmnsc           rmnsc           rmnss           rmnss
!    lmnsc           zmnsc           zmnsc           rmnsc
!                    zmncc           zmncs           rmncs
!                    lmnsc           lmnsc           zmnsc
!                    lmncc           lmncs           zmncs
!                                                    zmncc
!                                                    zmnss
!                                                    lmnsc
!                                                    lmncs
!                                                    lmncc
!                                                    lmnss



  !> stacked array of R, Z, Lambda Spectral force coefficients (see above for stack order)
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: gc

  !> stacked array of scaled R, Z, Lambda Fourier coefficients (see above for stack order)
  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: xc

  !> "velocity": change of Fourier coefficients per time step
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: xcdot

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: xsave

  !> backup copy of last-known-good xc
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: xstore

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: scalxc

END MODULE xstuff
