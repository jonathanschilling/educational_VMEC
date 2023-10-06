!> \file
!> \brief Compute the spectral width of the surface geometry Fourier coefficients.

!> \brief Compute the spectral width of the surface geometry Fourier coefficients.
!>
!> @param rmn Fourier coefficients of \f$R\f$
!> @param zmn Fourier coefficients of \f$Z\f$
SUBROUTINE spectrum(rmn, zmn)

  USE vmec_main
  USE vmec_params, ONLY: mscale, nscale, ntmax, rss, zcs, rsc, zcc

  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: rmn
  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: zmn

  INTEGER :: ntype
  INTEGER :: n
  INTEGER :: m
  INTEGER :: js
  REAL(rprec) :: scale_fac
  REAL(rprec), DIMENSION(ns) :: t1
  REAL(rprec), DIMENSION(ns) :: dnumer
  REAL(rprec), DIMENSION(ns) :: denom

  ! CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES,
  ! R+(at rsc) = .5(rsc + zcc),
  ! R-(at zcc) = .5(rsc - zcc),
  ! TO REQUIRED rsc, zcc FORMS

! #ifndef _HBANGLE
  IF (lthreed) CALL convert_sym  (rmn(:,:,:,rss), zmn(:,:,:,zcs))
  IF (lasym)   CALL convert_asym (rmn(:,:,:,rsc), zmn(:,:,:,zcc))
! #end /* ndef _HBANGLE */

  dnumer(2:ns) = zero
  denom(2:ns) = zero

  DO ntype = 1,ntmax
    DO n = 0,ntor
      DO m = 1,mpol1

         scale_fac = (mscale(m)*nscale(n))**2.0_dp

         DO js = 2,ns
            t1(js) =(rmn(js,n,m,ntype)**2.0_dp + zmn(js,n,m,ntype)**2.0_dp)*scale_fac
         END DO

         dnumer(2:ns) = dnumer(2:ns) + t1(2:ns)*xmpq(m,3)
         denom (2:ns) = denom (2:ns) + t1(2:ns)*xmpq(m,2)

      END DO
    END DO
  ENDDO

  specw(2:ns) = dnumer(2:ns)/denom(2:ns)

END SUBROUTINE spectrum
