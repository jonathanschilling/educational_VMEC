!> \file
!> \brief Normalization parameters for \f$\lambda\f$

!> \brief Normalization parameters for \f$\lambda\f$
!>
!> @param overg inverse of Jacobian \f$1/\sqrt{g}\f$
!> @param guu metric element \f$g_{\theta \theta}\f$
!> @param guv metric element \f$g_{\theta \zeta}\f$
!> @param gvv metric element \f$g_{\zeta \zeta}\f$
SUBROUTINE lamcal(overg, guu, guv, gvv)
  USE vmec_main
  USE vmec_params, ONLY: ntmax, jlam, lamscale
  USE realspace, ONLY: sqrts

  use dbgout

  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,nznt), INTENT(in) :: overg, guu, guv, gvv

  REAL(rprec), PARAMETER :: damping_fac=2.0_dp

  INTEGER :: m,n,js
  REAL(rprec) :: tnn, tnm, tmm, power, pfactor0, pfactor

  blam(:ns) = SUM(guu*overg, dim=2) ! over surface
  clam(:ns) = SUM(gvv*overg, dim=2) ! over surface
  dlam(:ns) = SUM(guv*overg, dim=2) ! over surface
  blam(1) = blam(2) ! constant extrapolation to axis
  clam(1) = clam(2) ! constant extrapolation to axis
  dlam(1) = dlam(2) ! constant extrapolation to axis
  blam(ns+1) =  0.0_dp   ! virtual "ghost" point beyond LCFS
  clam(ns+1) =  0.0_dp   ! virtual "ghost" point beyond LCFS
  dlam(ns+1) =  0.0_dp   ! virtual "ghost" point beyond LCFS
  DO js = 2, ns
     blam(js) = cp5*(blam(js) + blam(js+1))
     clam(js) = cp5*(clam(js) + clam(js+1))
     dlam(js) = cp5*(dlam(js) + dlam(js+1))
  END DO

  faclam = 0.0_dp
  pfactor0 = damping_fac/(2.0_dp*r0scale*lamscale)**2.0_dp

  DO m = 0, mpol1
     tmm = m*m
     power = MIN(tmm/256.0_dp, 8._dp)

     pfactor = pfactor0
     DO n = 0, ntor
        IF (m.eq.0 .and. n.eq.0) CYCLE

        ! sometimes helps convergence
        ! IF (n .gt. 1) pfactor = pfactor0/4

        tnn = (n*nfp)**2.0_dp
        tnm = 2.0_dp*m*n*nfp
        DO js = jlam(m), ns

           ! b: coupling between n and n ?
           ! d: coupling between n and m ?
           ! c: coupling between m and m ?

           faclam(js,n,m,1) = (blam(js)*tnn + SIGN(dlam(js),blam(js))*tnm + clam(js)*tmm)
           IF (faclam(js,n,m,1) .eq. zero) then
               faclam(js,n,m,1) = -1.E-10_dp
           end if

           ! Damps m > 16 modes
           faclam(js,n,m,1) = (pfactor/faclam(js,n,m,1)) * sqrts(js)**power

        END DO
     END DO
  END DO

  ! extend to rest of Fourier basis
  DO n = 2, ntmax
     faclam(:ns,0:ntor,0:mpol1,n) = faclam(:ns,0:ntor,0:mpol1,1)
  END DO

  ! ADD NORM FOR CHIP (PREVIOUSLY IOTA) FORCE, STORED IN lmnsc(m=0,n=0) COMPONENT
  DO js = 1, ns
     faclam(js,0,0,1) = (pfactor0*lamscale**2)/blam(js)
  END DO

  ! check lamcal output
  if (open_dbg_context("lamcal", num_eqsolve_retries)) then

    call add_real("pfactor0", pfactor0)

    call add_real_1d("blam", ns+1, blam)
    call add_real_1d("clam", ns+1, clam)
    call add_real_1d("dlam", ns+1, dlam)

    call add_real_3d("faclam", ns, ntor+1, mpol, faclam(:ns,0:ntor,0:mpol1,1)) !, order=(/ 2, 3, 1 /) ) ! TODO: check order !

    call close_dbg_out()
  end if

END SUBROUTINE lamcal
