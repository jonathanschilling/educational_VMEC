!> \file
!> \brief Fourier transforms

!> \brief Extends \f$B_s\f$ from \c ntheta2 interval to full \c ntheta3 interval in angle \f$\theta\f$.
!>
!> @param bout output \f$B_s\f$
!> @param bs_s symmetric part of \f$B_s\f$
!> @param bs_a anti-symmetric part of \f$B_s\f$
SUBROUTINE fext_fft (bout, bs_s, bs_a)
  USE vmec_main
  IMPLICIT NONE

  REAL(rprec), DIMENSION(nzeta,ntheta3), INTENT(out) :: bout
  REAL(rprec), DIMENSION(nzeta,ntheta2), INTENT(in) ::  bs_s, bs_a

  INTEGER :: ir, i, kz, kzr

  ! Extends bs from ntheta2 interval to full ntheta3 interval in angle u
  ! bs_s ~ cos(mu-nv) (     symmetric)
  ! bs_a ~ sin(mu-nv) (anti-symmetric)
  ! ntheta2 = pi
  bout(:,1:ntheta2) = bs_s(:,1:ntheta2) + bs_a(:,1:ntheta2)
  DO i = 1+ntheta2, ntheta3
     ! -theta
     ir = ntheta1+2-i
     DO kz= 1, nzeta
        ! -zeta
        ! kzr = ireflect(kz*ns)/ns
        kzr = nzeta+2-kz
        IF (kz .eq. 1) kzr = 1

        bout(kz,i) = bs_s(kzr,ir) - bs_a(kzr,ir)
     END DO
  END DO
END SUBROUTINE fext_fft
