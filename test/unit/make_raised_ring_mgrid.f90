!> Write the mgrid file named on the command line: a vertical field of decay index
!> 0.5, in which a ring of major radius 2.5 m is vertically stable, and a 1/R toroidal field.
program make_raised_ring_mgrid
  use netcdf
  use stel_kinds, only: rprec
  implicit none

  integer, parameter :: ir = 41, jz = 41, stringsize = 30
  real(rprec), parameter :: rmin = 1.5_rprec, rmax = 3.5_rprec, zmin = -1.0_rprec, zmax = 1.0_rprec
  real(rprec), parameter :: r0 = 2.5_rprec, minor = 0.5_rprec, b0 = 1.0_rprec, current = 250.0e3_rprec
  real(rprec), parameter :: decay_index = 0.5_rprec
  real(rprec), parameter :: pi = 4.0_rprec*atan(1.0_rprec), mu0 = 4.0e-7_rprec*pi

  character(len=512) :: filename
  character(len=stringsize) :: coil_group
  integer :: ncid, varid, i, j
  integer :: stringsize_dim, groups_dim, mode_dim, coils_dim, rad_dim, zee_dim, phi_dim
  integer :: br_id, bp_id, bz_id, group_id, mode_id, raw_id
  real(rprec) :: vertical, r, z
  real(rprec), dimension(ir, jz, 1) :: br, bp, bz

  ! the vertical field that balances the hoop force of a thin ring
  vertical = -mu0*current/(4.0_rprec*pi*r0)*(log(8.0_rprec*r0/minor) - 1.25_rprec)
  do j = 1, jz
    do i = 1, ir
      r = rmin + (rmax - rmin)*(i - 1)/(ir - 1)
      z = zmin + (zmax - zmin)*(j - 1)/(jz - 1)
      br(i, j, 1) = -vertical*decay_index*r*z/r0**2
      bp(i, j, 1) = b0*r0/r
      bz(i, j, 1) = vertical*(1.0_rprec + 0.5_rprec*decay_index*(1.0_rprec - (r**2 - 2.0_rprec*z**2)/r0**2))
    end do
  end do

  call get_command_argument(1, filename)
  call checked(nf90_create(trim(filename), nf90_clobber, ncid))
  call checked(nf90_def_dim(ncid, 'stringsize', stringsize, stringsize_dim))
  call checked(nf90_def_dim(ncid, 'external_coil_groups', 1, groups_dim))
  call checked(nf90_def_dim(ncid, 'dim_00001', 1, mode_dim))
  call checked(nf90_def_dim(ncid, 'external_coils', 1, coils_dim))
  call checked(nf90_def_dim(ncid, 'rad', ir, rad_dim))
  call checked(nf90_def_dim(ncid, 'zee', jz, zee_dim))
  call checked(nf90_def_dim(ncid, 'phi', 1, phi_dim))
  call checked(nf90_def_var(ncid, 'coil_group', nf90_char, [stringsize_dim, groups_dim], group_id))
  call checked(nf90_def_var(ncid, 'mgrid_mode', nf90_char, [mode_dim], mode_id))
  call checked(nf90_def_var(ncid, 'raw_coil_cur', nf90_double, [coils_dim], raw_id))
  call checked(nf90_def_var(ncid, 'br_001', nf90_double, [rad_dim, zee_dim, phi_dim], br_id))
  call checked(nf90_def_var(ncid, 'bp_001', nf90_double, [rad_dim, zee_dim, phi_dim], bp_id))
  call checked(nf90_def_var(ncid, 'bz_001', nf90_double, [rad_dim, zee_dim, phi_dim], bz_id))
  call define_scalars()
  call checked(nf90_enddef(ncid))

  call put_integer('ir', ir)
  call put_integer('jz', jz)
  call put_integer('kp', 1)
  call put_integer('nfp', 1)
  call put_integer('nextcur', 1)
  call put_real('rmin', rmin)
  call put_real('rmax', rmax)
  call put_real('zmin', zmin)
  call put_real('zmax', zmax)
  coil_group = 'vertical_and_toroidal'
  call checked(nf90_put_var(ncid, group_id, coil_group))
  call checked(nf90_put_var(ncid, mode_id, 'R'))
  call checked(nf90_put_var(ncid, raw_id, [1.0_rprec]))
  call checked(nf90_put_var(ncid, br_id, br))
  call checked(nf90_put_var(ncid, bp_id, bp))
  call checked(nf90_put_var(ncid, bz_id, bz))
  call checked(nf90_close(ncid))

contains

  subroutine define_scalars()
    integer :: k
    character(len=7), parameter :: integers(5) = [character(len=7) :: 'ir', 'jz', 'kp', 'nfp', 'nextcur']
    character(len=4), parameter :: reals(4) = [character(len=4) :: 'rmin', 'rmax', 'zmin', 'zmax']
    do k = 1, size(integers)
      call checked(nf90_def_var(ncid, trim(integers(k)), nf90_int, varid))
    end do
    do k = 1, size(reals)
      call checked(nf90_def_var(ncid, trim(reals(k)), nf90_double, varid))
    end do
  end subroutine define_scalars

  subroutine put_integer(name, value)
    character(len=*), intent(in) :: name
    integer, intent(in) :: value
    call checked(nf90_inq_varid(ncid, name, varid))
    call checked(nf90_put_var(ncid, varid, value))
  end subroutine put_integer

  subroutine put_real(name, value)
    character(len=*), intent(in) :: name
    real(rprec), intent(in) :: value
    call checked(nf90_inq_varid(ncid, name, varid))
    call checked(nf90_put_var(ncid, varid, value))
  end subroutine put_real

  subroutine checked(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      error stop 'cannot write the mgrid file'
    end if
  end subroutine checked

end program make_raised_ring_mgrid
