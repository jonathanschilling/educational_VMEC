!> \file
!> \brief Compute invariant residuals

!> \brief Compute invariant residuals
!>
!> @param gcr \f$R\f$-component of forces
!> @param gcz \f$Z\f$-component of forces
!> @param gcl \f$\lambda\f$-component of forces
SUBROUTINE residue (gcr, gcz, gcl)

  USE vmec_main, p5 => cp5
  USE vmec_params, ONLY: rss, zcs, rsc, zcc, meven, modd, ntmax
  USE xstuff

  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: gcr
  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: gcz
  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: gcl

  INTEGER, PARAMETER :: n0=0
  INTEGER, PARAMETER :: m0=0
  INTEGER, PARAMETER :: m1=1
  INTEGER, PARAMETER :: n3d=0
  INTEGER, PARAMETER :: nasym=1

  INTEGER :: jedge, j, n, m, i
  INTEGER :: delIter
  REAL(rprec) :: r1

  character(len=255) :: dump_filename
  logical            :: dump_physical_gc = .false.
  logical            :: dump_fsq = .false.
  logical            :: dump_scale_m1 = .false.

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
  IF (lthreed) CALL constrain_m1(gcr(:,:,m1,rss), gcz(:,:,m1,zcs))
  IF (lasym)   CALL constrain_m1(gcr(:,:,m1,rsc), gcz(:,:,m1,zcc))
! #end /* ndef _HBANGLE */

  ! dump physical forces
  if (dump_physical_gc) then
    write(dump_filename, 998) trim(input_extension)
998 format('phys_gc.',a)

    open(unit=42, file=trim(dump_filename), status="unknown")

    write(42, *) "# ns ntor mpol1 ntmax"
    write(42, *) ns, ntor, mpol1, ntmax

    write(42, *) "# j n m ntmax gcr gcz gcl"
    do j=1, ns
      do n=0, ntor
        do m=0, mpol1
          do i=1, ntmax
            write(42, *) j, n, m, i, &
              gcr(j, n, m, i), &
              gcz(j, n, m, i), &
              gcl(j, n, m, i)
          end do
        end do
      end do
    end do

    print *, "dumped physical gc to '"//trim(dump_filename)//"'"
    close(42)

    stop
  end if

  ! COMPUTE INVARIANT RESIDUALS
  r1 = one/(2*r0scale)**2 ! --> actually look at r1*fnorm --> scaling factor for forces (?)
  jedge = 0

  ! SPH-JAH013108: MUST INCLUDE EDGE FORCE (INITIALLY) FOR V3FITA TO WORK
  ! ADD A V3FIT RELATED FLAG? ADD fsq criterion first
  delIter = iter2-iter1

  ! Coding for VMEC2000 run stand-alone
  IF (delIter.lt.50 .and. (fsqr+fsqz).lt.1.E-6_dp) then
     ! include edge contribution only if converged well enough fast enough (?)
!      print *, "include edge force in residue"
     jedge = 1
  end if

  CALL getfsq (gcr, gcz, fsqr, fsqz, r1*fnorm, jedge)

  fsql = fnormL*SUM(gcl*gcl)
  fedge = r1*fnorm * SUM(gcr(ns,:,:,:)**2 + gcz(ns,:,:,:)**2)

  if (dump_fsq) then
    write(dump_filename, 997) trim(input_extension)
997 format('fsq.',a)

    open(unit=42, file=trim(dump_filename), status="unknown")

    write(42, *) "# r0scale r1 fnorm fnormL jedge fsqr fsqz fsql fedge"
    write(42, *) r0scale, r1, fnorm, fnormL, jedge, fsqr, fsqz, fsql, fedge

    print *, "dumped fsq to '"//trim(dump_filename)//"'"
    close(42)

    stop
  end if

  ! PERFORM PRECONDITIONING AND COMPUTE RESIDUES

! #ifndef _HBANGLE
  ! m = 1 constraint scaling
  IF (lthreed) CALL scale_m1(gcr(:,:,1,rss), gcz(:,:,1,zcs))
  IF (lasym)   CALL scale_m1(gcr(:,:,1,rsc), gcz(:,:,1,zcc))

  ! dump forces after scale_m1 has been applied
  if ((lthreed .or. lasym) .and. dump_scale_m1) then
    write(dump_filename, 996) trim(input_extension)
996 format('scale_m1.',a)

    open(unit=42, file=trim(dump_filename), status="unknown")

    write(42, *) "# ns ntor mpol1 ntmax"
    write(42, *) ns, ntor, mpol1, ntmax

    write(42, *) "# j n m ntmax gcr gcz"
    do j=1, ns
      do n=0, ntor
        do m=0, mpol1
          do i=1, ntmax
            write(42, *) j, n, m, i, &
              gcr(j, n, m, i), &
              gcz(j, n, m, i)
          end do
        end do
      end do
    end do

    print *, "dumped scale_m1 output to '"//trim(dump_filename)//"'"
    close(42)

    stop
  end if





  jedge = 0
  CALL scalfor (gcr, arm, brm, ard, brd, crd, jedge)
  jedge = 1
  CALL scalfor (gcz, azm, bzm, azd, bzd, crd, jedge)
! #end /* ndef _HBANGLE */

  !SPH: add fnorm1 ~ 1/R**2, since preconditioned forces gcr,gcz ~ Rmn or Zmn
  CALL getfsq (gcr, gcz, fsqr1, fsqz1, fnorm1, m1) ! m1 is simply == 1 --> include edge

  !SPH: THIS IS NOT INVARIANT UNDER PHIP->A*PHIP, AM->A**2*AM IN PROFIL1D
  !     (EXTCUR -> A*EXTCUR for FREE BOUNDARY)
  gcl = faclam*gcl
  fsql1 = hs*SUM(gcl*gcl)
  !030514      fsql1 = hs*lamscale**2*SUM(gcl*gcl)

END SUBROUTINE residue

!> \brief Compute internal \c gr , \c gz required for \f$m=1\f$ constraint
!>
!> @param gcr \f$R\f$-component of forces
!> @param gcz \f$Z\f$-component of forces
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
  IF (fsqz.LT.FThreshold .OR. iter2.LT.2) then
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
