!> Solve the vacuum field of NESTOR on a circular torus of major radius 1 m and
!> minor radius 0.3 m around two circular current rings in Z = 0, at R = 0.94 m
!> and R = 1.06 m, with opposite currents of 10 kA. Outside the torus the field
!> of the rings is the gradient of a single-valued potential, so the total
!> vacuum field on the boundary is zero and the field NESTOR returns there is
!> its error. The command line gives nfp, ntor, nzeta, ntheta, lasym and the
!> tolerance, 5e-4 if not given. The test stops with an error unless the rms
!> of that field over the rms of the field of the rings is below the tolerance.
program test_nestor_two_rings
  use stel_kinds, only: dp
  use vacmod0, only: set_nestor_sizes, nu3, nuv2
  use vacmod, only: allocate_nestor, free_mem_nestor, brv, bphiv, bzv, brad, bphi, bz
  use mgrid_mod, only: nr0b, nz0b, np0b, rminb, rmaxb, zminb, zmaxb, delrb, delzb, bvac
  implicit none

  real(dp), parameter :: pi = 4*atan(1.0_dp), mu0 = 4.0e-7_dp*pi
  real(dp), parameter :: major = 1.0_dp, minor = 0.3_dp
  real(dp), parameter :: ring_radius(2) = [0.94_dp, 1.06_dp]
  real(dp), parameter :: ring_current(2) = [1.0e4_dp, -1.0e4_dp]
  integer, parameter :: mpol = 8, mnmax = 2

  ! The field of the rings fills one plane of the mgrid arrays, which becoil
  ! uses for every toroidal plane. The even number of points in Z keeps the
  ! grid off the plane of the rings.
  integer, parameter :: nr = 701, nz = 700
  real(dp), parameter :: spacing = 1.0e-3_dp

  real(dp) :: tolerance = 5.0e-4_dp

  ! The rings carry no net toroidal current, rbtor is only compared with the
  ! R*BTOR of the vacuum field, and signgs is the sign of the Jacobian of this
  ! boundary.
  real(dp), parameter :: plascur = 0, rbtor = 1, signgs = -1

  character(len=32) :: argument
  integer :: nfp, ntor, nzeta, ntheta, ivac, ier_flag, i, k, l
  logical :: lasym
  real(dp) :: xm(mnmax), xn(mnmax), rmnc(mnmax), zmns(mnmax), zero_mn(mnmax)
  real(dp) :: r, z, residual
  real(dp), allocatable :: wint(:), raxis(:), zaxis(:)

  call get_command_argument(1, argument)
  read (argument, *) nfp
  call get_command_argument(2, argument)
  read (argument, *) ntor
  call get_command_argument(3, argument)
  read (argument, *) nzeta
  call get_command_argument(4, argument)
  read (argument, *) ntheta
  call get_command_argument(5, argument)
  read (argument, *) lasym
  if (command_argument_count() > 5) then
    call get_command_argument(6, argument)
    read (argument, *) tolerance
  end if

  nr0b = nr
  nz0b = nz
  np0b = 1
  rminb = 0.65_dp
  rmaxb = rminb + (nr - 1)*spacing
  zmaxb = (nz - 1)*spacing/2
  zminb = -zmaxb
  delrb = spacing
  delzb = spacing
  allocate(bvac(nr*nz, 3))
  bvac = 0
  do k = 1, nz
    z = zminb + (k - 1)*spacing
    do i = 1, nr
      r = rminb + (i - 1)*spacing
      l = i + nr*(k - 1)
      call add_ring_field(ring_radius(1), ring_current(1), r, z, bvac(l, 1), bvac(l, 3))
      call add_ring_field(ring_radius(2), ring_current(2), r, z, bvac(l, 1), bvac(l, 3))
    end do
  end do

  call set_nestor_sizes(nfp, ntor, mpol, nzeta, ntheta, lasym)
  call allocate_nestor

  xm = [0, 1]
  xn = 0
  rmnc = [major, minor]
  zmns = [0.0_dp, minor]
  zero_mn = 0

  ! The weights of VMEC: uniform on the full poloidal grid of lasym, and on the
  ! poloidal half grid with half weight at its ends otherwise.
  allocate(wint(nuv2), raxis(nzeta), zaxis(nzeta))
  if (lasym) then
    wint = 1.0_dp/nuv2
  else
    wint = 1.0_dp/(nzeta*(nu3 - 1))
    wint(:nzeta) = wint(:nzeta)/2
    wint(nuv2 - nzeta + 1:) = wint(nuv2 - nzeta + 1:)/2
  end if
  raxis = major
  zaxis = 0

  ivac = 0
  ier_flag = 0
  call vacuum(rmnc, zero_mn, zmns, zero_mn, xm, xn, plascur, rbtor, wint, 0, ivac, &
              mnmax, ier_flag, lasym, signgs, raxis, zaxis)

  residual = sqrt(sum(brv**2 + bphiv**2 + bzv**2)/sum(brad**2 + bphi**2 + bz**2))
  print '(a,i0,a,i0,a,i0,a,i0,a,l1,a,es10.3)', 'nfp = ', nfp, ', ntor = ', ntor, ', nzeta = ', nzeta, &
        ', ntheta = ', ntheta, ', lasym = ', lasym, ': rms of the vacuum field over that of the rings = ', residual
  call free_mem_nestor
  if (residual > tolerance) error stop 'the vacuum field on the boundary exceeds the tolerance'

contains

  !> Add the field of a circular ring of radius a and current I in Z = 0 at
  !> (r, z) to b_r and b_z.
  subroutine add_ring_field(a, current, r, z, b_r, b_z)
    real(dp), intent(in) :: a, current, r, z
    real(dp), intent(inout) :: b_r, b_z
    real(dp) :: alpha2, beta2, beta, c, k, e

    alpha2 = (r - a)**2 + z**2
    beta2 = (r + a)**2 + z**2
    beta = sqrt(beta2)
    call elliptic_integrals(4*a*r/beta2, k, e)
    c = mu0*current/pi
    b_r = b_r + c*z/(2*alpha2*beta*r)*((a**2 + r**2 + z**2)*e - alpha2*k)
    b_z = b_z + c/(2*alpha2*beta)*((a**2 - r**2 - z**2)*e + alpha2*k)
  end subroutine add_ring_field

  !> Complete elliptic integrals K and E of parameter m, by the
  !> arithmetic-geometric mean.
  subroutine elliptic_integrals(m, k, e)
    real(dp), intent(in) :: m
    real(dp), intent(out) :: k, e
    real(dp) :: a, b, c, a_next, weight, sum_c2

    a = 1
    b = sqrt(1 - m)
    c = sqrt(m)
    weight = 0.5_dp
    sum_c2 = weight*c**2
    do while (c > epsilon(a)*a)
      c = (a - b)/2
      a_next = (a + b)/2
      b = sqrt(a*b)
      a = a_next
      weight = 2*weight
      sum_c2 = sum_c2 + weight*c**2
    end do
    k = pi/(2*a)
    e = k*(1 - sum_c2)
  end subroutine elliptic_integrals

end program test_nestor_two_rings
