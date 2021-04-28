!> \file
!> \brief Extends function from \c ntheta2 to \c ntheta3 range.

!> \brief Extends function from \c ntheta2 to \c ntheta3 range.
!>
!> @param bsubsu tangential derivative of covariant magnetic field component \f$\partial B_s / \partial \theta\f$
!> @param bsubsv tangential derivative of covariant magnetic field component \f$\partial B_s / \partial \zeta\f$
SUBROUTINE fsym_invfft (bsubsu, bsubsv)
  USE vmec_main, ONLY: rprec, ns, nzeta, ntheta1, ntheta2, ntheta3, ireflect
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(inout) :: bsubsu, bsubsv

  INTEGER :: ir, i, jkz, jkr

  ! EXTENDS FUNCTION FROM ntheta2 to ntheta3 range
  ! ASSUMES bsubsu,v(0) ~ cos(mu-nv)   bsubsu,v(1) ~ sin(mu-nv)

  DO i = 1+ntheta2, ntheta1
     ! -theta
     ir = ntheta1+2-i
     DO jkz= 1, ns*nzeta
        ! -zeta
        jkr = ireflect(jkz)
        bsubsu(jkz,i,0) = bsubsu(jkr,ir,0) - bsubsu(jkr,ir,1)
        bsubsv(jkz,i,0) = bsubsv(jkr,ir,0) - bsubsv(jkr,ir,1)
     END DO
  END DO

  bsubsu(:,:ntheta2,0)=bsubsu(:,:ntheta2,0) + bsubsu(:,:ntheta2,1)
  bsubsv(:,:ntheta2,0)=bsubsv(:,:ntheta2,0) + bsubsv(:,:ntheta2,1)

END SUBROUTINE fsym_invfft
