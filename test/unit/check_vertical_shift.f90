!> Stop with an error unless the two wout files named on the command line hold the
!> same equilibrium apart from its vertical offset.
program check_vertical_shift
  use netcdf
  use stel_kinds, only: rprec
  implicit none

  real(rprec), parameter :: tolerance = 1.0e-9_rprec
  character(len=4), parameter :: names(6) = [character(len=4) :: &
       'rmnc', 'zmns', 'lmns', 'rmns', 'zmnc', 'lmnc']
  character(len=512) :: filename(2)
  integer :: ncid(2), i, k, varid, dimid, ns, mnmax
  real(rprec), allocatable :: a(:,:), b(:,:)
  real(rprec) :: largest

  do k = 1, 2
    call get_command_argument(k, filename(k))
    call checked(nf90_open(trim(filename(k)), nf90_nowrite, ncid(k)))
  end do
  call checked(nf90_inq_dimid(ncid(1), 'radius', dimid))
  call checked(nf90_inquire_dimension(ncid(1), dimid, len=ns))
  call checked(nf90_inq_dimid(ncid(1), 'mn_mode', dimid))
  call checked(nf90_inquire_dimension(ncid(1), dimid, len=mnmax))
  allocate(a(mnmax, ns), b(mnmax, ns))

  largest = 0.0_rprec
  do i = 1, size(names)
    call checked(nf90_inq_varid(ncid(1), names(i), varid))
    call checked(nf90_get_var(ncid(1), varid, a))
    call checked(nf90_inq_varid(ncid(2), names(i), varid))
    call checked(nf90_get_var(ncid(2), varid, b))
    ! the m = 0, n = 0 coefficient of Z is the offset itself
    if (names(i) == 'zmnc') then
      a(1, :) = 0.0_rprec
      b(1, :) = 0.0_rprec
    end if
    print '(a,a,es10.2)', names(i), ': largest difference', maxval(abs(a - b))
    largest = max(largest, maxval(abs(a - b)))
  end do

  do k = 1, 2
    call checked(nf90_close(ncid(k)))
  end do
  if (largest > tolerance) error stop 'a vertical shift changes the equilibrium'

contains

  subroutine checked(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      error stop 'cannot read the wout file'
    end if
  end subroutine checked

end program check_vertical_shift
