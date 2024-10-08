!> \file
!> \brief Interpolate \f$R\f$, \f$Z\f$ and \f$\lambda\f$ on full grid.

!> \brief Interpolate \f$R\f$, \f$Z\f$ and \f$\lambda\f$ on full grid.
!>
!> @param xnew interpolated state vector (nsnew surfaces)
!> @param xold interpolation basis: old state vector (nsold surfaces)
!> @param scalxc scaling factors to normalize the new state vector to
!> @param nsnew new number of flux surfaces
!> @param nsold old number of flux surfaces
SUBROUTINE interp(xnew, xold, scalxc, nsnew, nsold)
  USE vmec_main, ONLY: dp, rprec, mnsize
  USE vmec_params, ONLY: ntmax
  USE vmec_persistent, ONLY: ixm
  use vmec_dim
  use vmec_input

  use dbgout

  IMPLICIT NONE

  INTEGER, intent(in) :: nsnew, nsold
  REAL(rprec), DIMENSION(nsnew,mnsize,3*ntmax), INTENT(out)   :: xnew
  REAL(rprec), DIMENSION(nsnew,mnsize,3*ntmax), INTENT(in)    :: scalxc
  REAL(rprec), DIMENSION(nsold,mnsize,3*ntmax), intent(inout) :: xold

  REAL(rprec), PARAMETER :: zero=0.0_dp, one=1.0_dp

  INTEGER :: ntype, js
  integer, dimension(nsnew) :: js1, js2
  REAL(rprec) :: hsold
  real(rprec), dimension(nsnew) :: sj, s1, xint

  IF (nsold .le. 0) RETURN

  hsold = one/(nsold - 1.0_dp)

  ! INTERPOLATE R,Z AND LAMBDA ON FULL GRID
  ! (EXTRAPOLATE M=1 MODES,OVER SQRT(S), TO ORIGIN)
  ! ON ENTRY, XOLD = X(COARSE MESH) * SCALXC(COARSE MESH)
  ! ON EXIT,  XNEW = X(NEW MESH)   [ NOT SCALED BY 1/SQRTS ]
  DO ntype = 1, 3*ntmax

     ! extrapolation to axis for odd-m modes (?)
     WHERE (MOD(ixm(:mnsize), 2) .eq. 1) &
       xold(1,:,ntype) = 2.0_dp*xold(2,:,ntype) - xold(3,:,ntype)

     ! radial interpolation from old, coarse state vector to new, finer state vector
     DO js = 1, nsnew

        sj(js) = REAL(js - 1, rprec)/(nsnew - 1.0_dp)

        js1(js) = 1 + ((js - 1)*(nsold - 1))/(nsnew - 1)
        js2(js) = MIN(js1(js) + 1, nsold)

        s1(js) = (js1(js) - 1.0_dp)*hsold

        xint(js) = (sj(js) - s1(js))/hsold
        xint(js) = MIN(one, xint(js))
        xint(js) = MAX(zero,xint(js))

        xnew(js,:,ntype) = (   (one - xint(js))*xold(js1(js),:,ntype) &
                             +        xint(js) *xold(js2(js),:,ntype)   )/scalxc(js,:,1)
     END DO

     ! Zero M=1 modes at origin
     ! Actually, all odd-m modes are zeroed!
     WHERE (MOD(ixm(:mnsize), 2) .eq. 1) &
        xnew(1,:,ntype) = 0.0_dp

  END DO

  if (open_dbg_context("interp")) then

    call add_real_1d("sj", nsnew, sj)
    call add_int_1d("js1", nsnew, js1-1) ! index shift: Fortran starts at 1, Java at 0
    call add_int_1d("js2", nsnew, js2-1) ! index shift: Fortran starts at 1, Java at 0
    call add_real_1d("s1", nsnew, s1)
    call add_real_1d("xint", nsnew, xint)

    call add_real_5d("xold",   3, ntmax, nsold, ntor1, mpol, xold,   order=(/ 3, 4, 5, 2, 1 /) )
    call add_real_5d("xnew",   3, ntmax, nsnew, ntor1, mpol, xnew,   order=(/ 3, 4, 5, 2, 1 /) )
    call add_real_5d("scalxc", 3, ntmax, nsnew, ntor1, mpol, scalxc, order=(/ 3, 4, 5, 2, 1 /) )

    call close_dbg_out()
  end if

END SUBROUTINE interp
