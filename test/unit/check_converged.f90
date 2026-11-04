!> Stop with an error unless the wout file named on the command line holds a run
!> that terminated normally with its three force residuals at or below ftolv.
program check_converged
  use netcdf
  use stel_kinds, only: rprec
  implicit none

  integer :: ncid, varid, ier_flag, i
  real(rprec) :: ftolv, residuals(3)
  character(len=4), parameter :: names(3) = [character(len=4) :: 'fsqr', 'fsqz', 'fsql']
  character(len=512) :: filename

  call get_command_argument(1, filename)
  call checked(nf90_open(trim(filename), nf90_nowrite, ncid))
  call checked(nf90_inq_varid(ncid, 'ier_flag', varid))
  call checked(nf90_get_var(ncid, varid, ier_flag))
  call checked(nf90_inq_varid(ncid, 'ftolv', varid))
  call checked(nf90_get_var(ncid, varid, ftolv))
  do i = 1, 3
    call checked(nf90_inq_varid(ncid, names(i), varid))
    call checked(nf90_get_var(ncid, varid, residuals(i)))
  end do
  call checked(nf90_close(ncid))

  print '(a,i0,a,3es10.2,a,es10.2)', 'ier_flag = ', ier_flag, &
        ', fsqr fsqz fsql =', residuals, ', ftolv =', ftolv
  if (ier_flag /= 0) error stop 'the run did not terminate normally'
  if (any(residuals > ftolv)) error stop 'a force residual is above ftolv'

contains

  subroutine checked(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      error stop 'cannot read the wout file'
    end if
  end subroutine checked

end program check_converged
