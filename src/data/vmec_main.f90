!> \file
MODULE vmec_main

  USE vmec_dim
  USE vmec_input
  USE vmec_persistent
  USE vmec_params, ONLY: ndamp
  USE vparams

  IMPLICIT NONE

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ard
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: arm
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: brd
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: brm
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: azd
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: azm
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bzd
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bzm
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bmin
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bmax

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: crd
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: iotaf
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: phipf
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: chipf
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: phi
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: beta_vol
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: jcuru
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: jcurv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: jdotb
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: buco
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bvco
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bdotgradv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: equif
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: specw
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: tcon
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: psi
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: yellip
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: yinden
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: ytrian
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: yshift
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: ygeo
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: overr
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: sm
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: sp
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: pres
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: vp
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: jpar2
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: jperp2
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bdotb
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: blam
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: clam
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: dlam
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: vpphi
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: presgrad
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bdamp
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bucof
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bvcof
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: chi

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: presf !< pressure profile on full-grid, mass/phip**gamma
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: chips !< poloidal flux (same as chip), one-dimensional array
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: phips !< toroidal flux (same as phip), one-dimensional array
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: iotas !< rotational transform , on half radial mesh
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: icurv !< (-)toroidal current inside flux surface (vanishes like s)
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: mass  !< mass profile on half-grid

  REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: faclam
  REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: faclam0

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bsqsav

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bredge
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bpedge
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bzedge

  REAL(rprec), ALLOCATABLE :: xcl0(:)

  REAL(rprec), DIMENSION(0:mpol1d,3) :: xmpq
  REAL(rprec), DIMENSION(0:mpol1d) :: faccon

  REAL(rprec) :: hs !< radial mesh size increment
  REAL(rprec) :: currv
  REAL(rprec) :: aspect
  REAL(rprec) :: ohs
  REAL(rprec) :: voli
  REAL(rprec) :: r00
  REAL(rprec) :: r0scale
  REAL(rprec) :: z00
  REAL(rprec) :: fsqsum0
  REAL(rprec) :: fnorm
  REAL(rprec) :: fsqr=1
  REAL(rprec) :: fsqz=1
  REAL(rprec) :: fsql=1
  REAL(rprec) :: fnorm1
  REAL(rprec) :: fnorml
  REAL(rprec) :: fsqr1
  REAL(rprec) :: fsqz1
  REAL(rprec) :: fsql1
  REAL(rprec) :: fsq
  REAL(rprec) :: fedge
  REAL(rprec) :: wb
  REAL(rprec) :: wp

  REAL(rprec) :: router
  REAL(rprec) :: rinner

  REAL(rprec) :: ftolv

  !> time-step algorithm
  REAL(rprec) :: otav
  REAL(rprec), DIMENSION(ndamp) :: otau

  REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: rmn_bdy
  REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: zmn_bdy

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsubu0
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: dbsq
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rbsq

  REAL(rprec) :: rbtor
  REAL(rprec) :: rbtor0
  REAL(rprec) :: ctor
  REAL(rprec) :: delbsq
  REAL(rprec) :: res0
  REAL(rprec) :: delt0r

  REAL(rprec), DIMENSION(ndatafmax) :: spfa
  REAL(rprec), DIMENSION(ndatafmax) :: spfa2
  REAL(rprec), DIMENSION(ndatafmax) :: hp
  REAL(rprec), DIMENSION(ndatafmax) :: sifa
  REAL(rprec), DIMENSION(ndatafmax) :: sifa2
  REAL(rprec), DIMENSION(ndatafmax) :: hi

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
  integer :: iter1   !< number of iterations at which the currently active evolution was branched off from
  integer :: iter2   !< total number of iterations
  integer :: ivac    !< counts number of free-boundary iterations

  integer :: vacuum_calls = 0
  integer :: profil3d_calls = 0

END MODULE vmec_main
