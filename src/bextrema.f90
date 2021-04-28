!> \file
!> \brief Computes minimum and maximum \f$|\mathbf{B}|\f$ along \f$\zeta\f$ between two angle lines (\f$\theta = 0, \pi\f$).

!> \brief Computes minimum and maximum \f$|\mathbf{B}|\f$ along \f$\zeta\f$ between two angle lines (\f$\theta = 0, \pi\f$).
!>
!> @param modb magnitude of magnetic field \f$|\mathbf{B}|\f$
!> @param bmin minimum value of \f$|\mathbf{B}|\f$
!> @param bmax maximum value of \f$|\mathbf{B}|\f$
!> @param nzeta number of grid points in toroidal direction
!> @param ntheta number of grid points in poloidal direction
SUBROUTINE bextrema(modb, bmin, bmax, nzeta, ntheta)

  USE stel_kinds

  IMPLICIT NONE

  INTEGER, intent(in) :: nzeta, ntheta
  REAL(rprec), INTENT(in)  :: modb(nzeta,ntheta)
  REAL(rprec), INTENT(out) :: bmin(ntheta), bmax(ntheta)

  INTEGER :: ku

  DO ku = 1,ntheta
     bmin(ku)  = MINVAL(modb(:,ku))
     bmax(ku)  = MAXVAL(modb(:,ku))
  ENDDO

END SUBROUTINE bextrema
