!> \file
MODULE vmec_main
  USE vmec_dim
  USE vmec_input
  USE vmec_persistent
  USE vmec_params, ONLY: ndamp
  USE vparams
  IMPLICIT NONE

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: &
      ard, arm, brd, brm, azd, azm, bzd, bzm, bmin, bmax
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: &
      crd, iotaf, phipf, chipf, phi, beta_vol, &
      jcuru, jcurv, jdotb, &
      buco, bvco, bdotgradv, equif, specw, tcon, &
      psi, yellip, yinden, ytrian, yshift, ygeo, overr, &
      sm, sp, pres, vp, jpar2, jperp2, bdotb, &
      blam, clam, dlam, vpphi, presgrad, &
      bdamp, bucof, bvcof, chi
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: presf !< pressure profile on full-grid, mass/phip**gamma
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: chips !< poloidal flux (same as chip), one-dimensional array
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: phips !< toroidal flux (same as phip), one-dimensional array
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: iotas !< rotational transform , on half radial mesh
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: icurv !< (-)toroidal current inside flux surface (vanishes like s)
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: mass  !< mass profile on half-grid
  REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: faclam
  REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: faclam0

  REAL(rprec), ALLOCATABLE :: xcl0(:)

  REAL(rprec), DIMENSION(0:mpol1d,3) :: xmpq
  REAL(rprec), DIMENSION(0:mpol1d) :: faccon
  REAL(rprec) :: hs !< radial mesh size increment
  REAL(rprec) :: currv, aspect, ohs, voli, &
     r00, r0scale, z00, fsqsum0, &
     fnorm, fsqr=1, fsqz=1, fsql=1, fnorm1, fnorml, &
     fsqr1, fsqz1, fsql1, fsq, fedge, wb, wp

  REAL(rprec) :: ftolv

  !> time-step algorithm
  REAL(rprec) :: otav
  REAL(rprec), DIMENSION(ndamp) :: otau

  REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: rmn_bdy, zmn_bdy
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bsqsav
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsubu0, dbsq, rbsq
  REAL(rprec) :: rbtor, rbtor0, ctor, delbsq, res0, delt0r
  REAL(rprec), DIMENSION(ndatafmax) :: spfa, spfa2, hp, sifa, sifa2, hi

  LOGICAL :: lthreed
  LOGICAL :: lconm1

  LOGICAL :: lflip !< from init_geometry

  INTEGER, DIMENSION(:), ALLOCATABLE :: ireflect !< two-dimensional array for computing 2pi-v angle
  INTEGER :: multi_ns_grid
  INTEGER :: itfsq
  INTEGER :: ndatap
  INTEGER :: ndatai
  integer :: niterv  !< max iterations for current multi-grid iteration

  integer :: neqs    !< total number of equations to evolve (size of xc)
  integer :: irzloff !< offset in xc array between R,Z,L components
  integer :: iequi   !< counter used to call -EQFOR- at end of run
  integer :: ijacob  !< counter for number of times jacobian changes sign
  integer :: irst    !< "counter" monitoring sign of jacobian;
                     !< resets R, Z, and Lambda when jacobian changes sign
                     !< and decreases time step
  integer :: iter1   !< stores position in main iteration loop
  integer :: iter2   !< stores position in main iteration loop
  integer :: ivac    !< counts number of free-boundary iterations


END MODULE vmec_main
