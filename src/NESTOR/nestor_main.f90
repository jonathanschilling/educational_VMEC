program nestor
  use stel_kinds, only: dp
  use vacmod0, only: set_nestor_sizes
  use vacmod, only: allocate_nestor
  use nestor_io

  implicit none

  character(len=255) :: command_arg, vac_file
  INTEGER :: numargs, istat1

  print *, "Hello from NESTOR"

  numargs = iargc()
  if (numargs.ne.1) then
     stop "useage: xnestor vacin.nc"
  end if
  call getarg(1, command_arg)

  vac_file = trim(command_arg)

  print *, "  input filename: '",trim(vac_file),"'"

  call read_nestor_inputs(trim(vac_file))

  call set_nestor_sizes(nfp, ntor, mpol, nzeta, ntheta, lasym)

  call allocate_nestor

  CALL vacuum (rmnc, rmns, zmns, zmnc, xm, xn,                         &
               ctor, rbtor, wint, ivacskip, ivac, mnmax, ier_flag, &
               lasym, signgs, raxis, zaxis)


end program
