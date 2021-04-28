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
  ! bs_s ~ cos(mu-nv) (symmetric);   bs_a ~ sin(mu-nv) (anti-symmetric)
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

!> \brief Contract bs,bu,bv from full \c nu interval to half-u interval
!>        so cos, sin integrals can be performed on half-u interval.
!>
!> @param bs output \f$B_s\f$
!> @param bu output \f$B_\theta\f$
!> @param bv output \f$B_zeta\f$
!> @param bs_s symmetric part of \f$B_s\f$
!> @param bu_s symmetric part of \f$B_\theta\f$
!> @param bv_s symmetric part of \f$B_\zeta\f$
!> @param bs_a anti-symmetric part of \f$B_s\f$
!> @param bu_a anti-symmetric part of \f$B_\theta\f$
!> @param bv_a anti-symmetric part of \f$B_\zeta\f$
SUBROUTINE fsym_fft (bs, bu, bv, bs_s, bu_s, bv_s, bs_a, bu_a, bv_a)

  USE vmec_main

  IMPLICIT NONE

  REAL(rprec), DIMENSION(nzeta,ntheta3), INTENT(in)      :: bs
  REAL(rprec), DIMENSION(nzeta,ntheta3,0:1), INTENT(in)  :: bu, bv
  REAL(rprec), DIMENSION(nzeta,ntheta2,0:1), INTENT(out) :: bu_s, bv_s, bu_a, bv_a
  REAL(rprec), DIMENSION(nzeta,ntheta2) :: bs_s, bs_a

  INTEGER :: ir, i, kz, kzr

  ! CONTRACTS bs,bu,bv FROM FULL nu INTERVAL TO HALF-U INTERVAL
  ! SO COS,SIN INTEGRALS CAN BE PERFORMED ON HALF-U INTERVAL
  !
  ! bs_s(v,u) = .5*( bs(v,u) - bs(-v,-u) )     ! * SIN(mu - nv)
  ! bs_a(v,u) = .5*( bs(v,u) + bs(-v,-u) )     ! * COS(mu - nv)
  !
  ! bu, bv have opposite parity

  DO i = 1, ntheta2
     ! -theta
     ir = ntheta1+2-i
     IF (i == 1) ir = 1
     DO kz = 1, nzeta
        ! -zeta
        ! kzr = ireflect(ns*kz)/ns
        kzr = nzeta+2-kz
        IF (kz .eq. 1) kzr = 1
        bs_a(kz,i)      = cp5*(bs(kz,i)+bs(kzr,ir))
        bs_s(kz,i)      = cp5*(bs(kz,i)-bs(kzr,ir))
        bu_a(kz,i,:) = cp5*(bu(kz,i,:)-bu(kzr,ir,:))
        bu_s(kz,i,:) = cp5*(bu(kz,i,:)+bu(kzr,ir,:))
        bv_a(kz,i,:) = cp5*(bv(kz,i,:)-bv(kzr,ir,:))
        bv_s(kz,i,:) = cp5*(bv(kz,i,:)+bv(kzr,ir,:))
        END DO
     END DO

END SUBROUTINE fsym_fft
