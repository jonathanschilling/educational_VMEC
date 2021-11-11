!> \file
!> \brief Convert internal mode representation to standard form for output
!>        (coefficients of cos(mu-nv), sin(mu-nv) without internal \c mscale , \c nscale norms).

!> \brief Convert internal mode representation to standard form for output
!>        (coefficients of cos(mu-nv), sin(mu-nv) without internal \c mscale , \c nscale norms).
!>
!> @param rmnc stellarator-symmetric Fourier coefficients of \f$R\f$
!> @param zmns stellarator-symmetric Fourier coefficients of \f$Z\f$
!> @param lmns stellarator-symmetric Fourier coefficients of \f$\lambda\f$
!> @param rmns non-stellarator-symmetric Fourier coefficients of \f$R\f$
!> @param zmnc non-stellarator-symmetric Fourier coefficients of \f$Z\f$
!> @param lmnc non-stellarator-symmetric Fourier coefficients of \f$\lambda\f$
!> @param rzl_array state vector (all Fourier coefficients) of VMEC
!> @param js index of flux surface at which to do the conversion
SUBROUTINE convert(rmnc, zmns, lmns, rmns, zmnc, lmnc, rzl_array, js)

  USE vmec_main
  USE vmec_params

  IMPLICIT NONE

  REAL(rprec), DIMENSION(mnmax), INTENT(out) :: rmnc
  REAL(rprec), DIMENSION(mnmax), INTENT(out) :: zmns
  REAL(rprec), DIMENSION(mnmax), INTENT(out) :: lmns
  REAL(rprec), DIMENSION(mnmax), INTENT(out) :: rmns
  REAL(rprec), DIMENSION(mnmax), INTENT(out) :: zmnc
  REAL(rprec), DIMENSION(mnmax), INTENT(out) :: lmnc
  REAL(rprec), DIMENSION(ns, 0:ntor, 0:mpol1, 3*ntmax), INTENT(in) :: rzl_array
  INTEGER, INTENT(in) :: js

  REAL(rprec), PARAMETER :: p5 = 0.5_dp

  INTEGER :: rmncc, rmnss, rmncs, rmnsc, zmncs, zmnsc, &
             zmncc, zmnss, lmncs, lmnsc, lmncc, lmnss
  INTEGER :: mn, m, n, n1
  REAL(rprec) :: t1, sign0

  ! CONVERTS INTERNAL MODE REPRESENTATION TO STANDARD FORM FOR OUTPUT
  ! (COEFFICIENTS OF COS(mu-nv), SIN(mu-nv) WITHOUT internal mscale,nscale norms)

  rmncc = rcc
  rmnss = rss
  rmnsc = rsc
  rmncs = rcs
  zmnsc = zsc + ntmax
  zmncc = zcc + ntmax
  zmncs = zcs + ntmax
  zmnss = zss + ntmax
  lmnsc = zsc + 2*ntmax
  lmncc = zcc + 2*ntmax
  lmncs = zcs + 2*ntmax
  lmnss = zss + 2*ntmax

  ! init to zero, since not assigned if .not. lthreed
  zmns(1:ntor+1) = 0
  lmns(1:ntor+1) = 0

  ! DO M = 0 MODES SEPARATELY (ONLY KEEP N >= 0 HERE: COS(-NV), SIN(-NV))
  mn = 0
  m = 0
  DO n = 0, ntor
     t1 = mscale(m)*nscale(n)
     mn = mn + 1
     rmnc(mn) = t1*rzl_array(js,n,m,rmncc)
     IF (.not. lthreed) CYCLE
     zmns(mn) =-t1*rzl_array(js,n,m,zmncs)
     lmns(mn) =-t1*rzl_array(js,n,m,lmncs)
  END DO

  IF (lthreed .and. js.eq.1) THEN
     ! extrapolate to axis if 3D
     mn = 0
     DO n = 0, ntor
        t1 = mscale(m)*nscale(n)
        mn = mn + 1
        lmns(mn) =-t1*(2*rzl_array(2,n,m,lmncs) - rzl_array(3,n,m,lmncs))
     END DO
  END IF

  ! safeguard against spurious DC component in lambda (?)
  ! may have been used for storing iota variation...
  lmns(1) = 0

  ! now come the m>0, n=-ntor, ..., ntor entries
  DO m = 1, mpol1
     DO n = -ntor, ntor
        n1 = ABS(n)
        t1 = mscale(m)*nscale(n1)
        mn = mn + 1
        IF (n .eq. 0) THEN
           rmnc(mn) = t1*rzl_array(js,n,m,rmncc)
           zmns(mn) = t1*rzl_array(js,n,m,zmnsc)
           lmns(mn) = t1*rzl_array(js,n,m,lmnsc)
        ELSE IF (js .gt. 1) THEN
           sign0 = n/n1
           IF (.not.lthreed) sign0 = 0
           rmnc(mn) = p5*t1*(rzl_array(js,n1,m,rmncc)+sign0*rzl_array(js,n1,m,rmnss))
           zmns(mn) = p5*t1*(rzl_array(js,n1,m,zmnsc)-sign0*rzl_array(js,n1,m,zmncs))
           lmns(mn) = p5*t1*(rzl_array(js,n1,m,lmnsc)-sign0*rzl_array(js,n1,m,lmncs))
        ELSE IF (js .eq. 1) THEN
           ! no m>=1 component in magnetic axis (js==1)
           rmnc(mn) = 0
           zmns(mn) = 0
           lmns(mn) = 0
        END IF
     END DO
  END DO

  IF (mn .ne. mnmax) STOP 'Error in Convert!'

  IF (.not.lasym) THEN
     rmns = 0
     zmnc = 0
     lmnc = 0
     RETURN
  END IF

  ! add non-stellarator-symmetric terms now

  mn = 0; m = 0
  rmns(1:ntor+1) = 0
  DO n = 0, ntor
     t1 = mscale(m)*nscale(n)
     mn = mn + 1
     zmnc(mn) = t1*rzl_array(js,n,m,zmncc)
     lmnc(mn) = t1*rzl_array(js,n,m,lmncc)
     IF (.not.lthreed) CYCLE
     rmns(mn) =-t1*rzl_array(js,n,m,rmncs)                           !ers-fixed sign
  END DO

  DO m = 1, mpol1
     DO n = -ntor, ntor
        n1 = ABS(n)
        t1 = mscale(m)*nscale(n1)
        mn = mn + 1
        IF (n .eq. 0) THEN
           rmns(mn) = t1*rzl_array(js,n,m,rmnsc)
           zmnc(mn) = t1*rzl_array(js,n,m,zmncc)
           lmnc(mn) = t1*rzl_array(js,n,m,lmncc)
        ELSE IF (js .gt. 1) THEN
           sign0 = n/n1
           ! ers-corrected rmnsc <-> rmncs
           rmns(mn) = p5*t1*(rzl_array(js,n1,m,rmnsc)-sign0*rzl_array(js,n1,m,rmncs))
           zmnc(mn) = p5*t1*(rzl_array(js,n1,m,zmncc)+sign0*rzl_array(js,n1,m,zmnss))
           lmnc(mn) = p5*t1*(rzl_array(js,n1,m,lmncc)+sign0*rzl_array(js,n1,m,lmnss))
        ELSE IF (js .eq. 1) THEN
           rmns(mn) = 0
           zmnc(mn) = 0
           lmnc(mn) = 0
        END IF
     END DO
  END DO

END SUBROUTINE convert
