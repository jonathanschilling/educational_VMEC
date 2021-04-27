!> \file
SUBROUTINE lamcal(overg, guu, guv, gvv)
  USE vmec_main
  USE vmec_params, ONLY: ntmax, jlam, lamscale
  USE realspace, ONLY: sqrts
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,nznt), INTENT(in) :: overg, guu, guv, gvv

  REAL(rprec), PARAMETER :: damping_fac=2

  INTEGER :: m,n,js
  REAL(rprec) :: tnn, tnm, tmm, power, pfactor0, pfactor

  blam(:ns) = SUM(guu*overg, dim=2)
  clam(:ns) = SUM(gvv*overg, dim=2)
  dlam(:ns) = SUM(guv*overg, dim=2)
  blam(1) = blam(2)
  clam(1) = clam(2)
  dlam(1) = dlam(2)
  blam(ns+1) =  0
  clam(ns+1) =  0
  dlam(ns+1) =  0
  DO js = 2, ns
     blam(js) = cp5*(blam(js) + blam(js+1))
     clam(js) = cp5*(clam(js) + clam(js+1))
     dlam(js) = cp5*(dlam(js) + dlam(js+1))
  END DO

  faclam = 0
  pfactor0 = damping_fac/(2*r0scale*lamscale)**2

  DO m = 0, mpol1
     tmm = m*m
     power = MIN(tmm/256, 8._dp)

     pfactor = pfactor0
     DO n = 0, ntor
        IF (m.eq.0 .and. n.eq.0) CYCLE

        ! sometimes helps convergence
        ! IF (n .gt. 1) pfactor = pfactor0/4

        tnn = (n*nfp)**2
        tnm = 2*m*n*nfp
        DO js = jlam(m), ns
           faclam(js,n,m,1) = (blam(js)*tnn + SIGN(dlam(js),blam(js))*tnm + clam(js)*tmm)
           IF (faclam(js,n,m,1) .eq. zero) then
               faclam(js,n,m,1) = -1.E-10_dp
           end if

           ! Damps m > 16 modes
           faclam(js,n,m,1) = (pfactor/faclam(js,n,m,1)) * sqrts(js)**power

        END DO
     END DO
  END DO

  DO n = 2, ntmax
     faclam(:ns,0:ntor,0:mpol1,n) = faclam(:ns,0:ntor,0:mpol1,1)
  END DO

  ! ADD NORM FOR CHIP (PREVIOUSLY IOTA) FORCE, STORED IN lmnsc(m=0,n=0) COMPONENT
  DO js = 1, ns
     faclam(js,0,0,1) = (pfactor0*lamscale**2)/blam(js)
  END DO

END SUBROUTINE lamcal
