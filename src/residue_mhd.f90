!> \file
!> \brief Compute invariant residuals

!> \brief Compute invariant residuals
!>
!> @param gcr \f$R\f$-component of forces
!> @param gcz \f$Z\f$-component of forces
!> @param gcl \f$\lambda\f$-component of forces
SUBROUTINE residue_mhd (gcr, gcz, gcl, fsqrz, old_fsqz)

  USE vmec_main, p5 => cp5
  USE vmec_params, ONLY: rss, zcs, rsc, zcc, meven, modd, ntmax
  USE xstuff

  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: gcr
  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: gcz
  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: gcl
  real(rprec), intent(in) :: fsqrz, old_fsqz

  INTEGER, PARAMETER :: m1=1

  INTEGER :: jedge, delIter
  REAL(rprec) :: r1
  logical, parameter :: skip_scalfor_dbg = .true.

  ! IMPOSE M=1 MODE CONSTRAINT TO MAKE THETA ANGLE
  ! INVARIANT TO PHI-SHIFTS (AND THETA SHIFTS FOR ASYMMETRIC CASE)
  ! (ZCS = RSS, ZSS = RCS ARE THE CORRECT POLAR RELATIONS)

! #ifndef _HBANGLE
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
  IF (lthreed) CALL constrain_m1(gcr(:,:,m1,rss), gcz(:,:,m1,zcs), old_fsqz)
  IF (lasym)   CALL constrain_m1(gcr(:,:,m1,rsc), gcz(:,:,m1,zcc), old_fsqz)
! #end /* ndef _HBANGLE */

  ! COMPUTE INVARIANT RESIDUALS
  r1 = one/(2*r0scale)**2 ! --> actually look at r1*fnorm --> scaling factor for forces (?)
  jedge = 0

  ! SPH-JAH013108: MUST INCLUDE EDGE FORCE (INITIALLY) FOR V3FITA TO WORK
  ! ADD A V3FIT RELATED FLAG? ADD fsq criterion first
  delIter = iter2-iter1

  ! Coding for VMEC2000 run stand-alone
  IF (delIter.lt.50 .and. fsqrz.lt.1.E-6_dp) then
     ! include edge contribution only if converged well enough fast enough (?)
!      print *, "include edge force in residue"
     jedge = 1
  end if

  CALL getfsq (gcr, gcz, fsqr_mhd, fsqz_mhd, r1*fnorm, jedge)

  !fsql = fnormL*SUM(gcl*gcl)
  !fedge = r1*fnorm * SUM(gcr(ns,:,:,:)**2 + gcz(ns,:,:,:)**2)

  ! PERFORM PRECONDITIONING AND COMPUTE RESIDUES

! #ifndef _HBANGLE
  ! m = 1 constraint scaling
  IF (lthreed) CALL scale_m1(gcr(:,:,1,rss), gcz(:,:,1,zcs))
  IF (lasym)   CALL scale_m1(gcr(:,:,1,rsc), gcz(:,:,1,zcc))

  jedge = 0
  CALL scalfor (gcr, arm, brm, ard, brd, crd, jedge, skip_scalfor_dbg)
  jedge = 1
  CALL scalfor (gcz, azm, bzm, azd, bzd, crd, jedge, skip_scalfor_dbg)
! #end /* ndef _HBANGLE */

  !SPH: add fnorm1 ~ 1/R**2, since preconditioned forces gcr,gcz ~ Rmn or Zmn
  CALL getfsq (gcr, gcz, fsqr1_mhd, fsqz1_mhd, fnorm1, m1) ! m1 is simply == 1 --> include edge

  !SPH: THIS IS NOT INVARIANT UNDER PHIP->A*PHIP, AM->A**2*AM IN PROFIL1D
  !     (EXTCUR -> A*EXTCUR for FREE BOUNDARY)
  gcl = faclam*gcl
  !fsql1 = hs*SUM(gcl*gcl)
  !030514      fsql1 = hs*lamscale**2*SUM(gcl*gcl)

END SUBROUTINE residue_mhd
