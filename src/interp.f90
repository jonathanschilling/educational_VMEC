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
  USE vmec_main, ONLY: dp, rprec, mnsize, input_extension
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

  REAL(rprec), PARAMETER :: zero=0, one=1

  INTEGER :: ntype, js, mn
  integer, dimension(nsnew) :: js1, js2
  REAL(rprec) :: hsold
  real(rprec), dimension(nsnew) :: sj, s1, xint

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

     ! radial interpolation from old, coarse state vector to new, finer state vector
     DO js = 1, nsnew
        sj(js) = REAL(js - 1, rprec)/(nsnew - 1)
        js1(js) = 1 + ((js - 1)*(nsold - 1))/(nsnew - 1)
        js2(js) = MIN(js1(js) + 1, nsold)
        s1(js) = (js1(js) - 1)*hsold
        xint(js) = (sj(js) - s1(js))/hsold
        xint(js) = MIN(one, xint(js))
        xint(js) = MAX(zero,xint(js))
        xnew(js,:,ntype) = (   (one - xint(js))*xold(js1(js),:,ntype) &
                             +        xint(js) *xold(js2(js),:,ntype)   )/scalxc(js,:,1)
     END DO

     ! Zero M=1 modes at origin
     WHERE (MOD(ixm(:mnsize), 2) .eq. 1) &
        xnew(1,:,ntype) = 0
  END DO

  if (dump_interp) then
    write(dump_filename, 999) nsold, nsnew, trim(input_extension)
999 format('interp_',i5.5,'_',i5.5,'.',a,'.json')

    call open_dbg_out(dump_filename)

    call add_real_1d("sj", nsnew, sj)
    call add_int_1d("js1", nsnew, js1)
    call add_int_1d("js2", nsnew, js2)
    call add_real_1d("s1", nsnew, s1)
    call add_real_1d("xint", nsnew, xint)

    call add_real_5d("xold",   nsold, ntor1, mpol, ntmax, 2, xold,   order=(/ 1, 3, 4, 5, 2 /) )
    call add_real_5d("xnew",   nsnew, ntor1, mpol, ntmax, 2, xnew,   order=(/ 1, 3, 4, 5, 2 /) )
    call add_real_5d("scalxc", nsnew, ntor1, mpol, ntmax, 2, scalxc, order=(/ 1, 3, 4, 5, 2 /) ) 

    call close_dbg_out()
  end if

END SUBROUTINE interp
