!> \file
!> \brief allocate arrays depending on the number of flux surfaces \c ns

!> \brief allocate arrays depending on the number of flux surfaces \c ns
!>
!> @param linterp interpolate from coars to finer mesh?
!> @param neqs_old previous number of degrees-of-freedom, i.e., Fourier coefficients for \f$R\f$, \f$Z\f$ and \f$\lambda\f$
SUBROUTINE allocate_ns (linterp, neqs_old)
  USE vmec_main
  USE vmec_params, ONLY: ntmax
  USE realspace
  USE vforces
  USE xstuff
  USE mgrid_mod

  IMPLICIT NONE

  LOGICAL, INTENT(in) :: linterp
  INTEGER, INTENT(in) :: neqs_old

  INTEGER :: ndim, nsp1, istat1
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: xc_old, scalxc_old

  ndim  = 1 + nrzt ! TODO: why +1? some magical hidden storage at the end of the array ?
  nsp1  = 1 + ns   ! TODO: why +1? some magical hidden storage at the end of the array ?

  IF (neqs_old .gt. 0 .and. ALLOCATED(scalxc) .and. linterp) THEN
     ! Save old xc, scalxc for possible interpolation or IF iterations restarted on same mesh...
     ALLOCATE(xc_old(neqs_old), scalxc_old(neqs_old), stat=istat1)
     IF (istat1.ne.0) STOP 'allocation error #1 in allocate_ns'
         xc_old(:neqs_old) =     xc(:neqs_old)
     scalxc_old(:neqs_old) = scalxc(:neqs_old)
  END IF

  ! ALLOCATES MEMORY FOR NS-DEPENDENT ARRAYS
  ! FIRST BE SURE TO FREE MEMORY PREVIOUSLY ALLOCATED
  CALL free_mem_ns

  ALLOCATE (phip(ndim), chip(ndim), shalf(ndim), sqrts(ndim), wint(ndim), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #2 in allocate_ns'
  phip=0; chip=0; shalf=0; sqrts=0; wint=0

  ALLOCATE( ireflect(ns*nzeta), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #3 in allocate_ns'

  ALLOCATE( ard(nsp1,2),arm(nsp1,2),brd(nsp1,2),brm(nsp1,2),  &
            azd(nsp1,2),azm(nsp1,2),bzd(nsp1,2),bzm(nsp1,2), &
            sm(ns), sp(0:ns), bmin(ntheta2,ns), bmax(ntheta2,ns), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #6 in allocate_ns'

  ALLOCATE( iotaf(nsp1), crd(nsp1), mass(ns), phi(ns), presf(ns),   &
            jcuru(ns), jcurv(ns), jdotb(ns), buco(ns), bvco(ns),    &
            bucof(ns), bvcof(ns), chi(ns),                          &
            bdotgradv(ns), equif(ns), specw(ns), tcon(ns),          &
            psi(ns),yellip(2,ns),yinden(2,ns), ytrian(2,ns),yshift(2,ns),   &
            ygeo(2,ns),overr(ns), faclam(ns,0:ntor,0:mpol1,ntmax),    &
            iotas(nsp1), phips(nsp1), chips(nsp1), pres(nsp1),      &
            beta_vol(ns), jperp2(ns), jpar2(ns), bdotb(ns),         &
            phipf(ns), chipf(ns), blam(nsp1), clam(nsp1),           &
            dlam(nsp1), icurv(ns+1), vpphi(ns), bdamp(ns),          &
            presgrad(ns), vp(nsp1), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #7 in allocate_ns'

  iotaf(nsp1) = 0 ! TODO: why explicitly zero out only the last entry? hidden storage?

  ALLOCATE (gc(neqs), gc_con(neqs), gc_mhd(neqs), xcdot(neqs), xsave(neqs), xstore(neqs), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #9 in allocate_ns'
  xstore = zero

  IF (.not.ALLOCATED(xc)) THEN
     ALLOCATE (xc(neqs), scalxc(neqs), stat=istat1)
     IF (istat1.ne.0) STOP 'allocation error #10 in allocate_ns'
     xc = zero
  END IF

  ! FIRST STORE COARSE-MESH XC FOR INTERPOLATION
  IF (ALLOCATED(xc_old)) THEN
     xstore(1:neqs_old) =     xc_old(1:neqs_old)
     scalxc(1:neqs_old) = scalxc_old(1:neqs_old)
     DEALLOCATE (xc_old, scalxc_old)
  END IF

  ! Allocate nrzt-dependent arrays (persistent) for funct3d
  CALL allocate_funct3d

END SUBROUTINE allocate_ns
