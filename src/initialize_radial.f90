!> \file
!> \brief Allocates memory for radial arrays and initializes radial profiles.

!> \brief Allocates memory for radial arrays and initializes radial profiles.
!>
!> @param nsval new number of flux surfaces
!> @param ns_old old number of flux surfaces (from previous multi-grid iteration)
!> @param delt0 time step to be used in the new multi-grid iteration
SUBROUTINE initialize_radial(nsval, ns_old, delt0)
  USE vmec_main
  USE vmec_params, ONLY: ntmax
  USE realspace
  USE xstuff
  IMPLICIT NONE

  INTEGER, INTENT(in)      :: nsval
  INTEGER, INTENT(inout)   :: ns_old
  REAL(rprec), INTENT(out) :: delt0

  INTEGER :: neqs_old = 0 ! note that this sets neqs_old to zero only once at program start!
  LOGICAL :: lreset_internal, linterp

  ! Allocates memory for radial arrays and initializes radial profiles
  ! Loads data (if available) from a reset file

  ! print *, "initialize_radial"

  ! Set timestep control parameters
  fsq    = one

  iter2  = 1
  iter1  = iter2

  ijacob = 0
  first  = 1
  res0   = -1
  delt0  = delt

  ! start: INITIALIZE MESH-DEPENDENT SCALARS

  ! radial grid: only depends on ns
  ns = nsval
  ns1 = ns-1
  hs = one/ns1
  ohs = one/hs ! == ns1, but real-valued variant to avoid some kind of roundoff error ?

  ! real-space grid on surfaces
  nrzt = ns*nznt

  ! Fourier-space resolution
  mns = ns*mnsize ! number of flux surfaces * number of Fourier coeffs per surface --> total size of Fourier basis (n>=0)
  irzloff = ntmax*mns ! total number of Fourier coeffs for each of R, Z and Lambda --> including (lasym, lthreed)-dependent ntmax
  neqs = 3*irzloff ! degrees of freedom == total number of (R,Z,Lambda) Fourier coefficients

  ! end: INITIALIZE MESH-DEPENDENT SCALARS

  WRITE (nthreed, 10) ns, mnmax, ftolv, niterv
  PRINT 10, ns, mnmax, ftolv, niterv
10 FORMAT(/'  NS = ',i4,' NO. FOURIER MODES = ',i4,' FTOLV = ',1p,e10.3,' NITER = ',i6)

  lreset_internal = .true.

  ! check that interpolating from coarse to fine mesh
  ! and that old solution is available
  linterp = (ns_old.lt.ns .and. ns_old.ne.0)

  IF (ns_old .ne. ns) then

     ! ALLOCATE NS-DEPENDENT ARRAYS
     CALL allocate_ns(linterp, neqs_old)

     ! SAVE THIS FOR INTERPOLATION
     ! gc is re-used here for old, scaled xc
     IF (neqs_old.gt.0 .and. linterp) THEN
        gc(1:neqs_old)=scalxc(1:neqs_old)*xstore(1:neqs_old)
     END IF

     ! reset Fourier coefficients vector if lreset was specified
     xcdot = 0
     IF (lreset_internal) THEN
       xc = 0
     END IF

     ! COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
     CALL profil1d()

     ! TODO: lreset .and. .not.linter?
     ! If xc is overwritten by interp() anyway, why bother to initialize it in profil3d()?
     CALL profil3d(xc(1), xc(1+irzloff), lreset_internal)

     ! first.eq.1 at entry of restart_iter means to store xc in xstore
     first = 1
     CALL restart_iter(delt)

     ! INTERPOLATE FROM COARSE (ns_old) TO NEXT FINER (ns) RADIAL GRID
     IF (linterp) THEN
        ! print *, "interpolate from previous solution"

        CALL interp (xc, gc, scalxc, ns, ns_old)
     END IF

     ns_old = ns
     neqs_old = neqs
  end if

END SUBROUTINE initialize_radial
