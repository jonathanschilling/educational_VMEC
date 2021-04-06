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
      crd, iotaf, phipf, chipf, mass, phi, presf, beta_vol, &
      jcuru, jcurv, jdotb, &
      buco, bvco, bdotgradv, equif, specw, tcon, &
      psi, yellip, yinden, ytrian, yshift, ygeo, overr, &
      sm, sp, iotas, phips, chips, pres, vp, jpar2, jperp2, bdotb, &
      blam, clam, dlam, icurv, vpphi, presgrad, &
      r01, z01, bdamp, bucof, bvcof, chi
  REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: faclam, faclam0

  REAL(rprec), ALLOCATABLE :: xcl0(:)

  REAL(rprec), DIMENSION(0:mpol1d,3) :: xmpq
  REAL(rprec), DIMENSION(0:mpol1d) :: faccon
  REAL(rprec) :: dcon, currv, aspect, hs, ohs, voli, &
     signiota, rc0mse, r00, r0scale, z00, dkappa, fsqsum0, &
     pressum0, fnorm, fsqr=1, fsqz=1, fsql=1, fnorm1, fnorml, &
     fsqr1, fsqz1, fsql1, fsq, fedge, wb, wp, r00b, z00b, fz00_edge
  REAL(rprec) :: ftolv, otav, alphaR, alphaZ
  REAL(rprec), DIMENSION(ndamp) :: otau
  REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: rmn_bdy, zmn_bdy
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bsqsav
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsubu0, dbsq, rbsq
  REAL(rprec) :: rbtor, rbtor0, ctor, delbsq, res0, delt0r
  REAL(rprec), DIMENSION(ndatafmax) :: spfa, spfa2, hp, sifa, sifa2, hi
  LOGICAL :: lthreed, lconm1
  INTEGER, DIMENSION(:), ALLOCATABLE :: ireflect
  INTEGER :: multi_ns_grid, itfsq, iresidue, neqs, neqs1, neqs2, irzloff, ndatap, ndatai

  integer :: iequi  !< counter used to call -EQFOR- at end of run
  integer :: ijacob !< counter for number of times jacobian changes sign
  integer :: irst   !< "counter" monitoring sign of jacobian;
                    !< resets R, Z, and Lambda when jacobian changes sign
                    !< and decreases time step
  integer :: iter1  !< stores position in main iteration loop
  integer :: iter2  !< stores position in main iteration loop
  integer :: ivac   !< counts number of free-boundary iterations


END MODULE vmec_main
