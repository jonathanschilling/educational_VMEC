!> \file
!> \brief Compute preconditioning matrix elements for \f$R\f$, \f$Z\f$ force.

!> \brief Compute preconditioning matrix elements for \f$R\f$, \f$Z\f$ force.
!>
!> Note that in all parameter names, x=(r,z)
!>
!> @param lu1
!> @param bsq
!> @param gsqrt
!> @param r12
!> @param xs
!> @param xu12
!> @param xue
!> @param xuo
!> @param xodd
!> @param axm
!> @param axd
!> @param bxm
!> @param bxd
!> @param cx
!> @param eqfactor
!> @param trigmult
SUBROUTINE precondn(lu1, bsq, gsqrt, r12, &
                    xs, xu12, xue, xuo, xodd, &
                    axm, axd, bxm, bxd, cx, trigmult)
  USE vmec_main
  USE vmec_params, ONLY: signgs
  USE realspace
  IMPLICIT NONE

  REAL(rprec), DIMENSION(nrzt),   INTENT(in)  :: lu1, bsq, gsqrt, r12, xs, xu12, xue, xuo, xodd
  REAL(rprec), DIMENSION(ns+1,2), INTENT(out) :: axm, axd, bxm, bxd
  REAL(rprec), DIMENSION(ns+1),   INTENT(out) :: cx
  REAL(rprec), DIMENSION(nznt),   INTENT(in)  :: trigmult

  INTEGER :: js, l, lk
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ax, bx
  REAL(rprec) :: temp(ns+1)
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: ptau, ptau2
  REAL(rprec) :: t1, t2, t3, pfactor

  ! COMPUTE PRECONDITIONING MATRIX ELEMENTS FOR R,Z FORCE.
  ! NOTE THAT THE NORMALIZATION IS:
  !
  ! AX(off-diag) ~ <(cosmui cosmu cosnv cosnv) 2(R**2 * Xu**2 * bsq/gsqrt)>
  ! Factor of 2 arising from 1/gsqrt**2 in bsq
  !
  ! Now, cosmui cosmu ~ mscale(0)**2, cosnv**2 ~ nscale(0)**2
  ! Therefore, AX ~ (mscale(0)*nscale(0))**2 2<R**2*Xu**2*bsq/gsqrt>
  !               ~ 2*r0scale**2 <...>
  ALLOCATE (ax(ns+1,4), bx(ns+1,4), ptau(nznt), ptau2(nznt))

  ax = 0
  bx = 0
  cx = 0
  temp = 0

  ! pfactor = -2*r0scale**2       !v8.50
  pfactor = -4*r0scale**2        !restored in v8.51

  DO js = 2,ns
    ! COMPUTE DOMINANT (1/DELTA-S)**2 PRECONDITIONING MATRIX ELEMENTS
    lk = 0
    DO l = js,nrzt,ns
      lk = lk + 1
      t1 = pfactor*r12(l)*bsq(l)
      ptau2(lk) = r12(l)*t1/gsqrt(l)
      t1 = t1*wint(l)
      temp(js) = temp(js) + t1*trigmult(lk)*xu12(l)
      ptau(lk) = r12(l)*t1/gsqrt(l)
      t1 = xu12(l)*ohs
      t2 = cp25*(xue(l)/shalf(js) + xuo(l))/shalf(js)
      t3 = cp25*(xue(l-1)/shalf(js) + xuo(l-1))/shalf(js)
      ax(js,1) = ax(js,1) + ptau(lk)*t1*t1
      ax(js,2) = ax(js,2) + ptau(lk)*(-t1+t3)*(t1+t2)
      ax(js,3) = ax(js,3) + ptau(lk)*(t1+t2)*(t1+t2)
      ax(js,4) = ax(js,4) + ptau(lk)*(-t1+t3)*(-t1+t3)
    END DO

    ! COMPUTE PRECONDITIONING MATRIX ELEMENTS FOR M**2, N**2 TERMS
    lk = 0
    DO l = js,nrzt,ns
      lk = lk+1
      t1 = cp5*(xs(l) + cp5*xodd(l)/shalf(js))
      t2 = cp5*(xs(l) + cp5*xodd(l-1)/shalf(js))
      bx(js,1) = bx(js,1) + ptau(lk)*t1*t2
      bx(js,2) = bx(js,2) + ptau(lk)*t1*t1
      bx(js,3) = bx(js,3) + ptau(lk)*t2*t2
      cx(js) = cx(js) + cp25*pfactor*lu1(l)**2*gsqrt(l)*wint(l)
    END DO
  end do

  temp(1)=0
  temp(2:ns)=temp(2:ns)/vp(2:ns)
  temp(ns+1)=0
  DO js = 1,ns
    axm(js,1) =-ax(js,1)
    axd(js,1) = ax(js,1) + ax(js+1,1)
    axm(js,2) = ax(js,2) * sm(js) * sp(js-1)
    axd(js,2) = ax(js,3)*sm(js)**2 + ax(js+1,4)*sp(js)**2
    bxm(js,1) = bx(js,1)
    bxm(js,2) = bx(js,1) * sm(js) * sp(js-1)
    bxd(js,1) = bx(js,2) + bx(js+1,3)
    bxd(js,2) = bx(js,2)*sm(js)**2 + bx(js+1,3)*sp(js)**2
    cx(js)    = cx(js) + cx(js+1)
    temp(js) = signgs*(temp(js) + temp(js+1))
  END DO

  axm(ns+1,:) = 0
  axd(ns+1,:) = 0
  bxm(ns+1,:) = 0
  bxd(ns+1,:) = 0

  DEALLOCATE (ax, bx, ptau, ptau2)

END SUBROUTINE precondn
