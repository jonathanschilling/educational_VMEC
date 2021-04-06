!> \file
SUBROUTINE bextrema(modb, bmin, bmax, nzeta, ntheta)
  USE stel_kinds
  IMPLICIT NONE

  INTEGER, intent(in) :: nzeta, ntheta
  REAL(rprec), INTENT(in)  :: modb(nzeta,ntheta)
  REAL(rprec), INTENT(out) :: bmin(ntheta), bmax(ntheta)

  INTEGER :: ku

  ! Computes MAX, MIN of |B| along v (zeta) between two angle lines (theta = 0,pi)
  DO ku = 1,ntheta
     bmin(ku)  = MINVAL(modb(:,ku))
     bmax(ku)  = MAXVAL(modb(:,ku))
  ENDDO

END SUBROUTINE bextrema
