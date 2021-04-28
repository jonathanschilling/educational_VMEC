!> \file
!> \brief Get current at midplane (?)

!> \brief Get current at midplane (?)
!>
!> @param curmid current at midplane (?)
!> @param izeta index in toroidal direction
!> @param gsqrt Jacobian
!> @param r12 \f$R^2\f$
SUBROUTINE getcurmid (curmid, izeta, gsqrt, r12)

  USE vmec_input, ONLY: rprec, dp, nzeta
  USE vmec_dim, ONLY: ns, ns1, ntheta2

  implicit none

  REAL(rprec) :: curmid(2*ns)
  REAL(rprec) :: izeta(ns,nzeta,*), gsqrt(ns,nzeta,*), r12(ns,nzeta,*)

  REAL(rprec) :: midcur(ns)

  ! THETA = pi, PHI = 0
  midcur(2:ns) = r12(2:ns,1,ntheta2)/gsqrt(2:ns,1,ntheta2)

  curmid(1) = izeta(ns,1,ntheta2)*midcur(ns)
  curmid(2:ns1) = 0.5_dp*izeta(ns1:2:-1,1,ntheta2)*(midcur(ns1:2:-1) + midcur(ns:3:-1))

  ! THETA = 0, PHI = 0
  midcur(2:ns) = r12(2:ns,1,1)/gsqrt(2:ns,1,1)

  curmid(ns+1:2*ns-1) = 0.5_dp*izeta(2:ns1,1,1)*(midcur(2:ns1) + midcur(3:ns))

  curmid(ns) = 0.5_dp*(curmid(ns-1) + curmid(ns+1))
  curmid(2*ns) = 2*curmid(2*ns-1) - curmid(2*ns-2)

END SUBROUTINE getcurmid
