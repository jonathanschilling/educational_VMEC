!> \file
MODULE vmec_params

  USE stel_kinds, ONLY: rprec, dp
  USE vparams, ONLY: mpold

  implicit none

  INTEGER, PARAMETER :: meven = 0  !< parity selection label for even poloidal modes of R and Z
  INTEGER, PARAMETER :: modd = 1   !< parity selection label for  odd poloidal modes of R and Z
  INTEGER, PARAMETER :: ndamp = 10 !< number of iterations over which damping is averaged
  INTEGER, PARAMETER :: ns4 = 25

  INTEGER, PRIVATE :: ink                            ! 0,1, 2,3,...,mpold
  INTEGER, PARAMETER, DIMENSION(0:mpold) :: jmin1 = (/ 1,1,(2,ink=2,mpold) /) !< starting js(m) values where R,Z are non-zero
  INTEGER, PARAMETER, DIMENSION(0:mpold) :: jmin2 = (/ 1,2,(2,ink=2,mpold) /) !< starting js(m) values for which R,Z are evolved
  INTEGER, PARAMETER, DIMENSION(0:mpold) :: jlam  = (/ 2,2,(2,ink=2,mpold) /) !< starting js(m) values for which Lambda is evolved

  INTEGER, PARAMETER :: norm_term_flag       =  0
  INTEGER, PARAMETER :: bad_jacobian_flag    =  1
  ! 2 not used anymore?
  ! 3: 'VMEC INDATA ERROR: NCURR.ne.1 but BLOAT.ne.1.'
  INTEGER, PARAMETER :: jac75_flag           =  4
  INTEGER, PARAMETER :: input_error_flag     =  5
  ! 6 not used anymore?
  INTEGER, PARAMETER :: phiedge_error_flag   =  7
  INTEGER, PARAMETER :: ns_error_flag        =  8 ! 'NS ARRAY MUST NOT BE ALL ZEROES'
  INTEGER, PARAMETER :: misc_error_flag      =  9 ! (can happen in mgrid_mod)
  ! 10: 'VAC-VMEC I_TOR MISMATCH : BOUNDARY MAY ENCLOSE EXT. COIL'
  INTEGER, PARAMETER :: successful_term_flag = 11 ! ftol force criterion has been met

  ! These denote(d) bits in the call to runvmec().
  INTEGER, PARAMETER :: restart_flag     =  1 ! (1 << 0)
  INTEGER, PARAMETER :: readin_flag      =  2 ! (1 << 1)
  INTEGER, PARAMETER :: timestep_flag    =  4 ! (1 << 2)
  INTEGER, PARAMETER :: output_flag      =  8 ! (1 << 3)
  INTEGER, PARAMETER :: cleanup_flag     = 16 ! (1 << 4)
  INTEGER, PARAMETER :: reset_jacdt_flag = 32 ! (1 << 5)

  REAL(rprec), PARAMETER :: pdamp = 0.05_dp
  CHARACTER(LEN=*), PARAMETER :: version_ = '8.52'

  INTEGER :: ntmax !< number of contributing Fourier basis function (can be 1, 2 or 4); assigned in read_indata()

  INTEGER :: rcc
  INTEGER :: rss
  INTEGER :: rsc
  INTEGER :: rcs
  INTEGER :: zsc
  INTEGER :: zcs
  INTEGER :: zcc
  INTEGER :: zss ! stacking indices for xc, gc, ...

  INTEGER :: mnyq
  INTEGER :: nnyq

  INTEGER, ALLOCATABLE :: uminus(:)
  REAL(rprec), ALLOCATABLE :: mscale(:) !< array for norming theta-trig functions (internal use only)
                                        !< so that the discrete SUM[cos(mu)*cos(m'u)] = .5 delta(m,m')
  REAL(rprec), ALLOCATABLE :: nscale(:) !< array for norming zeta -trig functions (internal use only)

  REAL(rprec)              :: signgs    !< sign of Jacobian : must be =1 (right-handed) or =-1 (left-handed)

  REAL(rprec) :: lamscale=1

  INTEGER, PARAMETER :: m0=0 !< from totzsp
  INTEGER, PARAMETER :: m1=1 !< from totzsp
  INTEGER, PARAMETER :: n0=0 !< from totzsp


END MODULE vmec_params
