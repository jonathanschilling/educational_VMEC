program check_jacobian_recovery
  use netcdf
  use stel_kinds, only: rprec
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  integer :: ncid, varid, ns, ier, i
  real(rprec) :: residuals(3)
  character(len=4), parameter :: names(3) = [character(len=4) :: 'fsqr', 'fsqz', 'fsql']
  character(len=512) :: filename

  call get_command_argument(1, filename)
  call checked(nf90_open(trim(filename), nf90_nowrite, ncid))
  call checked(nf90_inq_varid(ncid, 'ns', varid))
  call checked(nf90_get_var(ncid, varid, ns))
  call checked(nf90_inq_varid(ncid, 'ier_flag', varid))
  call checked(nf90_get_var(ncid, varid, ier))
  do i = 1, 3
    call checked(nf90_inq_varid(ncid, names(i), varid))
    call checked(nf90_get_var(ncid, varid, residuals(i)))
  end do
  call checked(nf90_close(ncid))

  if (ns /= 31 .or. ier /= 0) error stop 'did not reach the requested equilibrium'
  if (.not. all(ieee_is_finite(residuals))) error stop 'non-finite force residual'
  if (any(residuals > 1.0e-9_rprec)) error stop 'requested force tolerance was not met'

contains

  subroutine checked(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      error stop 'cannot read equilibrium output'
    end if
  end subroutine checked

end program check_jacobian_recovery
