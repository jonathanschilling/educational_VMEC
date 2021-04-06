!> \file
SUBROUTINE residue (gcr, gcz, gcl)
  USE vmec_main, p5 => cp5
  USE vmec_params, ONLY: rss, zcs, rsc, zcc,                        &
                         meven, modd, ntmax, signgs
  USE realspace, ONLY: phip
  USE xstuff
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: gcr, gcz, gcl

  INTEGER, PARAMETER :: n0=0, m0=0, m1=1
  INTEGER, PARAMETER :: n3d=0, nasym=1

  INTEGER :: nsfix, jedge, delIter
  REAL(rprec) :: r1, tnorm, fac

  ! IMPOSE M=1 MODE CONSTRAINT TO MAKE THETA ANGLE
  ! INVARIANT TO PHI-SHIFTS (AND THETA SHIFTS FOR ASYMMETRIC CASE)
  ! (ZCS = RSS, ZSS = RCS ARE THE CORRECT POLAR RELATIONS)
  !
  ! SYMMETRIC PERTURBATIONS (BASED ON POLAR RELATIONS):
  !    RSS(n) = ZCS(n), n != 0
  ! ASYMMETRIC PERTURBATIONS:
  !    RSC(n) = ZCC(n), ALL n
  !
  ! INTERNALLY:
  !    XC(rss) = .5*(Rss + Zcs), XC(zcs) = .5*(Rss - Zcs) -> 0
  !    XC(rsc) = .5*(Rsc + Zcc), XC(zcc) = .5*(Rsc - Zcc) -> 0
  ! THIS IMPLIES THE CONSTRAINT
  !    3D ONLY : GC(zcs) = 0;
  !    ASYM:     GC(zcc) = 0
  IF (lthreed) CALL constrain_m1(gcr(:,:,m1,rss), gcz(:,:,m1,zcs))
  IF (lasym)   CALL constrain_m1(gcr(:,:,m1,rsc), gcz(:,:,m1,zcc))

  ! COMPUTE INVARIANT RESIDUALS
  r1 = one/(2*r0scale)**2
  jedge = 0

  ! SPH-JAH013108: MUST INCLUDE EDGE FORCE (INITIALLY) FOR V3FITA TO WORK
  ! ADD A V3FIT RELATED FLAG? ADD fsq criterion first
  delIter = iter2-iter1

  ! Coding for VMEC2000 run stand-alone
  IF (delIter.lt.50 .and. (fsqr+fsqz).lt.1.E-6_dp) then
     ! include edge contribution only under certain circumstances ?
     jedge = 1
  end if

  CALL getfsq (gcr, gcz, fsqr, fsqz, r1*fnorm, jedge)

  fsql = fnormL*SUM(gcl*gcl)
  fedge = r1*fnorm*SUM(gcr(ns,:,:,:)**2 + gcz(ns,:,:,:)**2)

  ! PERFORM PRECONDITIONING AND COMPUTE RESIDUES

  ! m = 1 constraint scaling
  IF (lthreed) CALL scale_m1(gcr(:,:,1,rss), gcz(:,:,1,zcs))
  IF (lasym)   CALL scale_m1(gcr(:,:,1,rsc), gcz(:,:,1,zcc))
  jedge = 0
  CALL scalfor (gcr, arm, brm, ard, brd, crd, jedge)
  jedge = 1
  CALL scalfor (gcz, azm, bzm, azd, bzd, crd, jedge)

  !SPH: add fnorm1 ~ 1/R**2, since preconditioned forces gcr,gcz ~ Rmn or Zmn
  CALL getfsq (gcr, gcz, fsqr1, fsqz1, fnorm1, m1)
  !SPH: THIS IS NOT INVARIANT UNDER PHIP->A*PHIP, AM->A**2*AM IN PROFIL1D
  !     (EXTCUR -> A*EXTCUR for FREE BOUNDARY)
  gcl = faclam*gcl
  fsql1 = hs*SUM(gcl*gcl)
  !030514      fsql1 = hs*lamscale**2*SUM(gcl*gcl)

END SUBROUTINE residue


SUBROUTINE constrain_m1(gcr, gcz)
  USE vmec_main, p5 => cp5
  IMPLICIT NONE

  REAL(dp), DIMENSION(ns,0:ntor), INTENT(inout) :: gcr, gcz

  REAL(dp), PARAMETER :: FThreshold = 1.E-6_dp
  REAL(dp) :: temp(ns,0:ntor)

  ! COMPUTE INTERNAL gr, gz
  ! NOTE: internal gz => 0 for both values of lconm1 (although gz is different)
  ! FOR lconm1=T, gcr(internal) = gcr+gcz, gcz(internal) = gcr-gcz->0
  IF (lconm1) THEN
     temp = gcr
     gcr = osqrt2*(gcr + gcz)
     gcz = osqrt2*(temp - gcz)
  END IF

  !v8.50: ADD iter2<2 so reset=<WOUT_FILE> works
  IF (fsqz.LT.FThreshold .OR. iter2.LT.2) gcz = 0

END SUBROUTINE constrain_m1


SUBROUTINE scale_m1(gcr, gcz)
  USE vmec_main
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor), INTENT(inout) :: gcr, gcz

  INTEGER, PARAMETER :: nodd=2
  INTEGER :: n
  REAL(rprec) :: fac(ns)

  IF (.not.lconm1) RETURN

  fac = (ard(:ns,nodd)+brd(:ns,nodd))/                                        &
        (ard(:ns,nodd)+brd(:ns,nodd)+azd(:ns,nodd)+bzd(:ns,nodd))
  DO n = 0, ntor
     gcr(:,n) = fac*gcr(:,n)
  END DO

  fac = (azd(:ns,nodd)+bzd(:ns,nodd))/                                        &
        (ard(:ns,nodd)+brd(:ns,nodd)+azd(:ns,nodd)+bzd(:ns,nodd))
  DO n = 0, ntor
     gcz(:,n) = fac*gcz(:,n)
  END DO

END SUBROUTINE scale_m1
