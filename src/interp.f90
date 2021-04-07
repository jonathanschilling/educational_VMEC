!> \file
SUBROUTINE interp(xnew, xold, scalxc, nsnew, nsold)
  USE vmec_main, ONLY: dp, rprec, mnsize
  USE vmec_params, ONLY: ntmax
  USE vmec_persistent, ONLY: ixm
  IMPLICIT NONE

  INTEGER, intent(in) :: nsnew, nsold
  REAL(rprec), DIMENSION(nsnew,mnsize,3*ntmax), INTENT(out)   :: xnew
  REAL(rprec), DIMENSION(nsnew,mnsize,3*ntmax), INTENT(in)    :: scalxc
  REAL(rprec), DIMENSION(nsold,mnsize,3*ntmax), intent(inout) :: xold

  REAL(rprec), PARAMETER :: zero=0, one=1

  INTEGER :: ntype, js, js1, js2
  REAL(rprec) :: hsold, sj, s1, xint

  IF (nsold .le. 0) RETURN

  hsold = one/(nsold - 1)

  ! INTERPOLATE R,Z AND LAMBDA ON FULL GRID
  ! (EXTRAPOLATE M=1 MODES,OVER SQRT(S), TO ORIGIN)
  ! ON ENTRY, XOLD = X(COARSE MESH) * SCALXC(COARSE MESH)
  ! ON EXIT,  XNEW = X(NEW MESH)   [ NOT SCALED BY 1/SQRTS ]
  DO ntype = 1, 3*ntmax

     ! extrapolation to axis for odd-m modes (?)
     WHERE (MOD(ixm(:mnsize), 2) .eq. 1) &
       xold(1,:,ntype) = 2*xold(2,:,ntype) - xold(3,:,ntype)

     DO js = 1, nsnew
        sj = REAL(js - 1, rprec)/(nsnew - 1)
        js1 = 1 + ((js - 1)*(nsold - 1))/(nsnew - 1)
        js2 = MIN(js1 + 1, nsold)
        s1 = (js1 - 1)*hsold
        xint = (sj - s1)/hsold
        xint = MIN(one,xint)
        xint = MAX(zero,xint)
        xnew(js,:,ntype) = (   (one - xint)*xold(js1,:,ntype) &
                             +        xint *xold(js2,:,ntype)   )/scalxc(js,:,1)
     END DO

     ! Zero M=1 modes at origin
     WHERE (MOD(ixm(:mnsize),2) .eq. 1) &
        xnew(1,:,ntype) = 0
  END DO

  END SUBROUTINE interp
