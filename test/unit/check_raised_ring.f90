!> Stop with an error unless the run in the wout file named on the command line
!> terminated normally with its boundary below the height it started from.
program check_raised_ring
  use netcdf
  use stel_kinds, only: rprec
  implicit none

  real(rprec), parameter :: raised = 1.0e-3_rprec
  character(len=512) :: filename
  integer :: ncid, varid, dimid, ier_flag, ns, mnmax
  real(rprec), allocatable :: zmnc(:,:)

  call get_command_argument(1, filename)
  call checked(nf90_open(trim(filename), nf90_nowrite, ncid))
  call checked(nf90_inq_varid(ncid, 'ier_flag', varid))
  call checked(nf90_get_var(ncid, varid, ier_flag))
  call checked(nf90_inq_dimid(ncid, 'radius', dimid))
  call checked(nf90_inquire_dimension(ncid, dimid, len=ns))
  call checked(nf90_inq_dimid(ncid, 'mn_mode', dimid))
  call checked(nf90_inquire_dimension(ncid, dimid, len=mnmax))
  allocate(zmnc(mnmax, ns))
  call checked(nf90_inq_varid(ncid, 'zmnc', varid))
  call checked(nf90_get_var(ncid, varid, zmnc))
  call checked(nf90_close(ncid))

  print '(a,i0,a,es11.3)', 'ier_flag = ', ier_flag, ', boundary height = ', zmnc(1, ns)
  if (ier_flag /= 0) error stop 'the run did not terminate normally'
  if (abs(zmnc(1, ns)) >= raised) error stop 'the boundary did not descend'

contains

  subroutine checked(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      error stop 'cannot read the wout file'
    end if
  end subroutine checked

end program check_raised_ring
