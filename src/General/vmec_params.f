!> \file
      MODULE vmec_params
      USE stel_kinds, ONLY: rprec, dp
      USE vparams, ONLY: mpold
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: meven = 0, modd = 1
      INTEGER, PARAMETER :: ndamp = 10
      INTEGER, PARAMETER :: ns4 = 25

      INTEGER, PRIVATE :: ink
      INTEGER, PARAMETER, DIMENSION(0:mpold) ::
     1  jmin1 = (/ 1,1,(2,ink=2,mpold) /),        !starting js(m) values where R,Z are non-zero
     2  jmin2 = (/ 1,2,(2,ink=2,mpold) /),        !starting js(m) values for which R,Z are evolved
     3  jlam  = (/ 2,2,(2,ink=2,mpold) /)         !starting js(m) values for which Lambda is evolved

      INTEGER, PARAMETER :: norm_term_flag=0, bad_jacobian_flag=1,
     2                      jac75_flag=4, input_error_flag=5,
     3                      phiedge_error_flag=7,
     4                      ns_error_flag=8,
     5                      misc_error_flag=9,
     6                      successful_term_flag=11 !ftol force criterion has been met
      INTEGER, PARAMETER :: restart_flag=1, readin_flag=2,
     1                      timestep_flag=4,output_flag=8,
     2                      cleanup_flag=16, reset_jacdt_flag=32

      REAL(rprec), PARAMETER :: pdamp = 0.05_dp
      CHARACTER(LEN=*), PARAMETER :: version_ = '8.52'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ntmax, rcc, rss, rsc, rcs, zsc, zcs, zcc, zss
      INTEGER :: mnyq, nnyq
      INTEGER, ALLOCATABLE :: uminus(:)
      REAL(rprec), ALLOCATABLE :: mscale(:), nscale(:)
      REAL(rprec) :: signgs, lamscale=1

      END MODULE vmec_params
