!> \file
!> \brief Main program of stand-alone version of NESTOR

!> \brief Main program of stand-alone version of NESTOR
program nestor
  use stel_kinds, only: dp
  use vacmod0, only: set_nestor_sizes
  use vacmod, only: allocate_nestor, free_mem_nestor, amatsav, bvecsav, bsubvvac
  use nestor_io

  use vmec_input, only: &
    dump_vac1n_vacuum  , &
    dump_vac1n_precal  , &
    dump_vac1n_surface , &
    dump_vac1n_bextern , &
    dump_vac1n_analyt  , &
    dump_vac1n_greenf  , &
    dump_vac1n_fourp   , &
    dump_vac1n_fouri   , &
    dump_vac1n_solver  , &
    dump_vac1n_bsqvac, &
    vi_ie => input_extension

  implicit none

  character(len=255) :: command_arg, vac_file, vacout_file
  INTEGER :: numargs, istat1

  print *, "Hello from NESTOR"

  vi_ie = "test.vmec"
  dump_vac1n_vacuum  = .true.
  dump_vac1n_precal  = .true.
  dump_vac1n_surface = .true.
  dump_vac1n_bextern = .true.
  dump_vac1n_analyt  = .true.
  dump_vac1n_greenf  = .true.
  dump_vac1n_fourp   = .true.
  dump_vac1n_fouri   = .true.
  dump_vac1n_solver  = .true.
  dump_vac1n_bsqvac  = .true.

  numargs = iargc()
  if (numargs.ne.1) then
     stop "usage: xnestor vacin.nc"
  end if
  call getarg(1, command_arg)

  vac_file = trim(command_arg)

  print *, "  input filename: '",trim(vac_file),"'"

  call read_nestor_inputs(trim(vac_file))

  print *, "set nestor sizes: ", nfp, ntor, mpol, nzeta, ntheta, lasym
  call set_nestor_sizes(nfp, ntor, mpol, nzeta, ntheta, lasym)

  call allocate_nestor

  ! copy over from input file, not though vacuum call
  amatsav = amatsav_nestor
  bvecsav = bvecsav_nestor
  bsubvvac = bsubvvac_nestor

  CALL vacuum (rmnc, rmns, zmns, zmnc, xm, xn,                         &
               ctor, rbtor, wint, ivacskip, ivac, mnmax, ier_flag, &
               lasym, signgs, raxis, zaxis)

  ! contruct filename into which to write output of stand-alone NESTOR
  write(vac_file, "(A,I6.6,A)") "vac/vacout_"//TRIM(input_extension)//"_", &
     vacuum_calls, ".nc"
  print *, "write NESTOR output to '"//trim(vac_file)//"'"
  call write_nestor_outputs(trim(vac_file), lasym, ivac, ier_flag)

  call free_mem_nestor

end program
