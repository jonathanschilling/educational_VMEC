      MODULE angle_constraints
      USE vmec_main, ONLY: ns, mpol, ntor, dp, mpol1, lthreed, lasym
      USE vmec_params, ONLY: signgs, ntmax, rcc, rss, zsc, zcs, rsc, rcs, zss, zcc
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
      INTEGER, PARAMETER :: pexp=4, m0=0, m1=1, m2=2, m3=3
      LOGICAL, PARAMETER :: lorigin=.FALSE.
      INTEGER            :: mrho, m, istat
      REAL, PARAMETER    :: p5=0.5_dp, zero=0
      REAL(dp), ALLOCATABLE :: t1m(:), t2m(:), cos_HB(:), sin_HB(:)
      REAL(dp), ALLOCATABLE :: rz_array0(:,:,:), xtempa(:)
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: arhod, arhom, brhod, brhom,   &
         ard2, arm2, azd2, azm2, arhod2, arhom2,                             &
         brd2, brm2, bzd2, bzm2, brhod2, brhom2
      REAL(dp), DIMENSION(:), ALLOCATABLE   :: crhod, sin2u, cos2u, sfact
      REAL(dp) :: sqp5
      END MODULE angle_constraints
