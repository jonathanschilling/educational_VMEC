!> \file
SUBROUTINE allocate_ns (linterp, neqs2_old)
  USE vmec_main
  USE vmec_params, ONLY: ntmax
  USE realspace
  USE vforces
  USE xstuff
  USE mgrid_mod
  USE fbal
  IMPLICIT NONE

  INTEGER, INTENT(in) :: neqs2_old
  LOGICAL, INTENT(inout) :: linterp

  INTEGER :: ndim, nsp1, istat1
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: xc_old, scalxc_old

  ! FIRST STORE COARSE-MESH XC FOR INTERPOLATION
  ndim  = 1 + nrzt
  nsp1  = 1 + ns

  IF (neqs2_old .gt. 0 .and. ALLOCATED(scalxc) .and. linterp) THEN
     ! Save old xc, scalxc for possible interpolation or IF iterations restarted on same mesh...
     ALLOCATE(xc_old(neqs2_old), scalxc_old(neqs2_old), stat=istat1)
     IF (istat1.ne.0) STOP 'allocation error #1 in allocate_ns'
         xc_old(:neqs2_old) =     xc(:neqs2_old)
     scalxc_old(:neqs2_old) = scalxc(:neqs2_old)
  END IF

  ! ALLOCATES MEMORY FOR NS-DEPENDENT ARRAYS
  ! FIRST BE SURE TO FREE MEMORY PREVIOUSLY ALLOCATED
  CALL free_mem_ns (.true.)

  ALLOCATE (phip(ndim), chip(ndim), shalf(ndim), sqrts(ndim), wint(ndim), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #2 in allocate_ns'
  phip=0; chip=0; shalf=0; sqrts=0; wint=0

  ALLOCATE( ireflect(ns*nzeta), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #3 in allocate_ns'

  ALLOCATE( ard(nsp1,2),arm(nsp1,2),brd(nsp1,2),brm(nsp1,2),  &
            azd(nsp1,2),azm(nsp1,2),bzd(nsp1,2), bzm(nsp1,2), &
            sm(ns), sp(0:ns), bmin(ntheta2,ns), bmax(ntheta2,ns), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #6 in allocate_ns'

  ALLOCATE( iotaf(nsp1), crd(nsp1), mass(ns), phi(ns), presf(ns),   &
            jcuru(ns), jcurv(ns), jdotb(ns), buco(ns), bvco(ns),    &
            bucof(ns), bvcof(ns), chi(ns),                          &
            bdotgradv(ns), equif(ns), specw(ns), tcon(ns),          &
            psi(ns),yellip(ns),yinden(ns), ytrian(ns),yshift(ns),   &
            ygeo(ns),overr(ns), faclam(ns,0:ntor,0:mpol1,ntmax),    &
            iotas(nsp1), phips(nsp1), chips(nsp1), pres(nsp1),      &
            beta_vol(ns), jperp2(ns), jpar2(ns), bdotb(ns),         &
            phipf(ns), chipf(ns), blam(nsp1), clam(nsp1),           &
            dlam(nsp1), rru_fac(ns), rzu_fac(ns), frcc_fac(ns),     &
            fzsc_fac(ns), icurv(ns+1), vpphi(ns), bdamp(ns),        &
            presgrad(ns), vp(nsp1), r01(ns), z01(ns), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #7 in allocate_ns'
  frcc_fac = 0
  fzsc_fac = 0

  iotaf(nsp1) = 0

  ALLOCATE (gc(neqs2), xcdot(neqs2), xsave(neqs2), xstore(neqs2), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #9 in allocate_ns'
  xstore = zero

  IF (.not.ALLOCATED(xc)) THEN
     ALLOCATE (xc(neqs2), scalxc(neqs2), stat=istat1)
     IF (istat1.ne.0) STOP 'allocation error #10 in allocate_ns'
     xc(:neqs2) = zero
  END IF

  IF (ALLOCATED(xc_old)) THEN
     xstore(1:neqs2_old) =     xc_old(1:neqs2_old)
     scalxc(1:neqs2_old) = scalxc_old(1:neqs2_old)
     DEALLOCATE (xc_old, scalxc_old)
  END IF

  ! Allocate nrzt-dependent arrays (persistent) for funct3d
  CALL allocate_funct3d

END SUBROUTINE allocate_ns
