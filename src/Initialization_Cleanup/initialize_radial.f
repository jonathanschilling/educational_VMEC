      SUBROUTINE initialize_radial(nsval, ns_old, delt0, lscreen)
      USE vmec_main
      USE vmec_params, ONLY: ntmax
      USE realspace
      USE xstuff
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: nsval
      INTEGER, INTENT(inout) :: ns_old
      REAL(rprec), INTENT(out) :: delt0
      LOGICAL, INTENT(in) :: lscreen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: neqs2_old=0
      LOGICAL :: lreset_internal, linterp
C-----------------------------------------------
!
!     Allocates memory for radial arrays and initializes radial profiles
!     Loads data (if available) from a reset file
!
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!        hs      radial mesh size increment
!        irzloff offset in xc array between R,Z,L components
!        neqs    total number of equations to evolve (size of xc)
C-----------------------------------------------
!     Set timestep control parameters
      fsq     = one
      iter2 = 1
      iter1 = iter2
      ijacob = 0
      irst = 1
      res0 = -1

!
!       INITIALIZE MESH-DEPENDENT SCALARS
!
      ns = nsval
      ns1 = ns-1
      delt0 = delt
      hs = one/ns1
      ohs = one/hs
      mns = ns*mnsize
      irzloff = ntmax*mns
      nrzt = nznt*ns
      neqs = 3*irzloff
      neqs1 = neqs + 1
      neqs2 = neqs1 + 1

      WRITE (nthreed, 10) ns, mnmax, ftolv, niter
       IF (lscreen) PRINT 10, ns, mnmax, ftolv, niter
   10 FORMAT(/'  NS = ',i4,' NO. FOURIER MODES = ',i4,' FTOLV = ',
     1   1p,e10.3,' NITER = ',i6)

!
!     ALLOCATE NS-DEPENDENT ARRAYS
!
      lreset_internal = .true.
      linterp = (ns_old.lt.ns .and. ns_old.ne.0)
      IF (ns_old .eq. ns) RETURN
      CALL allocate_ns(linterp, neqs2_old)

!
!     SAVE THIS FOR INTERPOLATION
!
      IF (neqs2_old.gt.0 .and. linterp) THEN
         gc(1:neqs2_old)=scalxc(1:neqs2_old)*xstore(1:neqs2_old)
      END IF
!
!     COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
!
      CALL profil1d (xc, xcdot, lreset_internal)
      CALL profil3d (xc(1), xc(1+irzloff), lreset_internal, linterp)

      irst = 1
      CALL restart_iter(delt)

!
!     INTERPOLATE FROM COARSE (ns_old) TO NEXT FINER (ns) RADIAL GRID
!
      IF (linterp) THEN
         CALL interp (xc, gc, scalxc, ns, ns_old)
      END IF
      ns_old = ns
      neqs2_old = neqs2

      END SUBROUTINE initialize_radial
