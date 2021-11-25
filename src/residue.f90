!> \file
!> \brief Compute invariant residuals

!> \brief Compute invariant residuals
!>
!> @param gcr \f$R\f$-component of forces
!> @param gcz \f$Z\f$-component of forces
!> @param gcl \f$\lambda\f$-component of forces
SUBROUTINE residue (gcr, gcz, gcl, fsqrz, old_fsqz)

  USE vmec_main, p5 => cp5
  USE vmec_params, ONLY: rss, zcs, rsc, zcc, meven, modd, ntmax
  USE xstuff

  use dbgout

  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: gcr
  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: gcz
  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: gcl
  real(rprec), intent(in) :: fsqrz, old_fsqz

  INTEGER, PARAMETER :: n0=0
  INTEGER, PARAMETER :: m0=0
  INTEGER, PARAMETER :: m1=1
  INTEGER, PARAMETER :: n3d=0
  INTEGER, PARAMETER :: nasym=1

  INTEGER :: jedge, j, n, m, i
  INTEGER :: delIter
  REAL(rprec) :: r1

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

  ! dump physical forces
  if (dump_physical_gc .and. iter2.le.max_dump) then
    write(dump_filename, 998) ns, iter2, trim(input_extension)
998 format('phys_gc_',i5.5,'_',i6.6,'.',a,'.json')

    call open_dbg_out(dump_filename)

    call add_real_4d("gcr", ns, ntmax, ntor1, mpol, &
            reshape(gcr, (/ ns, ntmax, ntor1, mpol /), order=(/ 1, 3, 4, 2 /) ) )
    call add_real_4d("gcz", ns, ntmax, ntor1, mpol, &
            reshape(gcz, (/ ns, ntmax, ntor1, mpol /), order=(/ 1, 3, 4, 2 /) ) )
    call add_real_4d("gcl", ns, ntmax, ntor1, mpol, &
            reshape(gcl, (/ ns, ntmax, ntor1, mpol /), order=(/ 1, 3, 4, 2 /) ) )

    call close_dbg_out()
  end if

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

  CALL getfsq (gcr, gcz, fsqr, fsqz, r1*fnorm, jedge)

  fsql = fnormL*SUM(gcl*gcl)
  fedge = r1*fnorm * SUM(gcr(ns,:,:,:)**2 + gcz(ns,:,:,:)**2)

  if (dump_fsq .and. iter2.le.max_dump) then
    write(dump_filename, 997) ns, iter2, trim(input_extension)
997 format('fsq_'i5.5,'_',i6.6,'.',a,'.json')

    call open_dbg_out(dump_filename)

    call add_real("r0scale", r0scale)
    call add_real("r1", r1)
    call add_real("fnorm", fnorm)
    call add_real("fnormL", fnormL)
    call add_int("jedge", jedge)
    call add_real("fsqr", fsqr)
    call add_real("fsqz", fsqz)
    call add_real("fsql", fsql)
    call add_real("fedge", fedge)

    call close_dbg_out()
  end if

  ! PERFORM PRECONDITIONING AND COMPUTE RESIDUES

! #ifndef _HBANGLE
  ! m = 1 constraint scaling
  IF (lthreed) CALL scale_m1(gcr(:,:,1,rss), gcz(:,:,1,zcs))
  IF (lasym)   CALL scale_m1(gcr(:,:,1,rsc), gcz(:,:,1,zcc))

  ! dump forces after scale_m1 has been applied
  if (dump_scale_m1 .and. iter2.le.max_dump) then
    write(dump_filename, 996) ns, iter2, trim(input_extension)
996 format('scale_m1_',i5.5,'_',i6.6,'.',a,'.json')

    call open_dbg_out(dump_filename)

    call add_real_4d("gcr", ns, ntmax, ntor1, mpol, &
            reshape(gcr, (/ ns, ntmax, ntor1, mpol /), order=(/ 1, 3, 4, 2 /) ) )
    call add_real_4d("gcz", ns, ntmax, ntor1, mpol, &
            reshape(gcz, (/ ns, ntmax, ntor1, mpol /), order=(/ 1, 3, 4, 2 /) ) )

    call close_dbg_out()
  end if

  jedge = 0
  CALL scalfor (gcr, arm, brm, ard, brd, crd, jedge)
  jedge = 1
  CALL scalfor (gcz, azm, bzm, azd, bzd, crd, jedge)
! #end /* ndef _HBANGLE */

  ! dump forces after scalfor has been applied
  if (dump_scalfor_out .and. iter2.le.max_dump) then
    write(dump_filename, 995) ns, iter2, trim(input_extension)
995 format('scalfor_out_',i5.5,'_',i6.6,'.',a,'.json')

    call open_dbg_out(dump_filename)

    call add_read_2d("arm", ns+1, 2, arm)
    call add_read_2d("ard", ns+1, 2, ard)
    call add_read_2d("brm", ns+1, 2, brm)
    call add_read_2d("brd", ns+1, 2, brd)
    call add_read_1d("crd", ns+1,    crd)
    call add_read_2d("azm", ns+1, 2, azm)
    call add_read_2d("azd", ns+1, 2, azd)
    call add_read_2d("bzm", ns+1, 2, bzm)
    call add_read_2d("bzd", ns+1, 2, bzd)

    call add_real_4d("gcr", ns, ntmax, ntor1, mpol, &
            reshape(gcr, (/ ns, ntmax, ntor1, mpol /), order=(/ 1, 3, 4, 2 /) ) )
    call add_real_4d("gcz", ns, ntmax, ntor1, mpol, &
            reshape(gcz, (/ ns, ntmax, ntor1, mpol /), order=(/ 1, 3, 4, 2 /) ) )

    call close_dbg_out()
  end if

  !SPH: add fnorm1 ~ 1/R**2, since preconditioned forces gcr,gcz ~ Rmn or Zmn
  CALL getfsq (gcr, gcz, fsqr1, fsqz1, fnorm1, m1) ! m1 is simply == 1 --> include edge

  !SPH: THIS IS NOT INVARIANT UNDER PHIP->A*PHIP, AM->A**2*AM IN PROFIL1D
  !     (EXTCUR -> A*EXTCUR for FREE BOUNDARY)
  gcl = faclam*gcl
  fsql1 = hs*SUM(gcl*gcl)
  !030514      fsql1 = hs*lamscale**2*SUM(gcl*gcl)

  if (dump_fsq1 .and. iter2.le.max_dump) then
    write(dump_filename, 994) ns, iter2, trim(input_extension)
994 format('fsq1_',i5.5,'_',i6.6,'.',a,'.json')

    call open_dbg_out(dump_filename)

    call add_real("fnorm1", fnorm1)
    call add_real("fsqr1", fsqr1)
    call add_real("fsqz1", fsqz1)
    call add_real("fsql1", fsql1)

    call close_dbg_out()
  end if

END SUBROUTINE residue

!> \brief Compute internal \c gr , \c gz required for \f$m=1\f$ constraint
!>
!> @param gcr \f$R\f$-component of forces
!> @param gcz \f$Z\f$-component of forces
SUBROUTINE constrain_m1(gcr, gcz, old_fsqz)

  USE vmec_main, p5 => cp5

  IMPLICIT NONE

  REAL(dp), DIMENSION(ns,0:ntor), INTENT(inout) :: gcr, gcz
  real(dp), intent(in) :: old_fsqz

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
  IF (old_fsqz.LT.FThreshold .OR. iter2.LT.2) then
     !write(*,*) "zero z-force in constrain_m1"

     ! ensure that the m=1 constraint is satisfied exactly
     ! --> the corresponding m=1 coeffs of R,Z are constrained to be zero
     !     and thus must not be "forced" (by the time evol using gc) away from zero
     gcz = 0
  end if

END SUBROUTINE constrain_m1

!> \brief Compute internal \c gr , \c gz required for \f$m=1\f$ constraint
!>
!> @param gcr \f$R\f$-component of forces
!> @param gcz \f$Z\f$-component of forces
SUBROUTINE scale_m1(gcr, gcz)

  USE vmec_main

  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor), INTENT(inout) :: gcr, gcz

  INTEGER, PARAMETER :: nodd=2
  INTEGER :: n
  REAL(rprec) :: fac(ns)

  IF (lconm1) then

     fac =  (ard(:ns,nodd)+brd(:ns,nodd)                            )      &
           /(ard(:ns,nodd)+brd(:ns,nodd)+azd(:ns,nodd)+bzd(:ns,nodd))
     DO n = 0, ntor
        gcr(:,n) = fac*gcr(:,n)
     END DO

     fac =  (                            azd(:ns,nodd)+bzd(:ns,nodd))      &
           /(ard(:ns,nodd)+brd(:ns,nodd)+azd(:ns,nodd)+bzd(:ns,nodd))
     DO n = 0, ntor
        gcz(:,n) = fac*gcz(:,n)
     END DO

  end if

END SUBROUTINE scale_m1
