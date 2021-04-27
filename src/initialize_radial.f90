!> \file
SUBROUTINE initialize_radial(nsval, ns_old, delt0)
  USE vmec_main
  USE vmec_params, ONLY: ntmax
  USE realspace
  USE xstuff
  IMPLICIT NONE

  INTEGER, INTENT(in)      :: nsval
  INTEGER, INTENT(inout)   :: ns_old
  REAL(rprec), INTENT(out) :: delt0

  INTEGER :: neqs_old=0
  LOGICAL :: lreset_internal, linterp

  ! Allocates memory for radial arrays and initializes radial profiles
  ! Loads data (if available) from a reset file

  print *, "initialize_radial"

  ! Set timestep control parameters
  fsq    = one
  iter2  = 1
  iter1  = iter2
  ijacob = 0
  irst   = 1
  res0   = -1

  ! INITIALIZE MESH-DEPENDENT SCALARS
  ns = nsval
  ns1 = ns-1
  delt0 = delt
  hs = one/ns1
  ohs = one/hs ! == ns1 ?
  mns = ns*mnsize ! number of flux surfaces * number of Fourier coeffs per surface --> total size of Fourier basis (n>=0)
  irzloff = ntmax*mns ! total number of Fourier coeffs for each of R, Z and Lambda --> including (lasym, lthreed)-dependent ntmax
  nrzt = nznt*ns
  neqs = 3*irzloff ! degrees of freedom == total number of (R,Z,Lambda) Fourier coefficients

  WRITE (nthreed, 10) ns, mnmax, ftolv, niterv
  PRINT 10, ns, mnmax, ftolv, niterv
10 FORMAT(/'  NS = ',i4,' NO. FOURIER MODES = ',i4,' FTOLV = ',1p,e10.3,' NITER = ',i6)

  ! ALLOCATE NS-DEPENDENT ARRAYS
  lreset_internal = .true.
  linterp = (ns_old.lt.ns .and. ns_old.ne.0)

  IF (ns_old .ne. ns) then
     CALL allocate_ns(linterp, neqs_old)

     ! SAVE THIS FOR INTERPOLATION
     IF (neqs_old.gt.0 .and. linterp) THEN
        gc(1:neqs_old)=scalxc(1:neqs_old)*xstore(1:neqs_old)
     END IF

     ! COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
     CALL profil1d (xc, xcdot, lreset_internal)
     CALL profil3d (xc(1), xc(1+irzloff), lreset_internal)

     irst = 1
     CALL restart_iter(delt)

     ! INTERPOLATE FROM COARSE (ns_old) TO NEXT FINER (ns) RADIAL GRID
     IF (linterp) THEN
        CALL interp (xc, gc, scalxc, ns, ns_old)
     END IF

     ns_old = ns
     neqs_old = neqs
  end if

END SUBROUTINE initialize_radial
