program scale_vmec

  use vmec_input
  use stel_kinds

  implicit none

  integer :: iunit = 11, istat, numargs
  real(dp) :: B_scale, R_scale
  character(len=200) :: inputFilename, arg

  print *,"Usage: scale_vmec input.* B_scale R_scale"

  numargs = command_argument_count() ! This function is a fortran 2003 intrinsic.
  if (numargs .ne. 3) stop "Wrong number of arguments"

  call get_command_argument(1, inputFilename)

  call get_command_argument(2, arg)
  read (arg,*) B_scale
  call get_command_argument(3, arg)
  read (arg,*) R_scale

  print *,"B_scale =",B_scale, "  R_scale =",R_scale

  open(unit=iunit, file=inputFilename, action="read", status="old", iostat=istat)
  if (istat /= 0) then
     print *,"Error opening input file ", trim(inputFilename)
     stop
  end if

  call read_indata_namelist(iunit, istat)
  if (istat .ne. 0) then
     print *,"Error reading namelist! istat=",istat
     stop
  else
     print *,"Successfully read indata namelist."
  end if

  close(unit = iunit)

  !------------------------------------------
  ! Now do the scaling.

  rbc = rbc * R_scale
  rbs = rbs * R_scale
  zbc = zbc * R_scale
  zbs = zbs * R_scale
  raxis_cc = raxis_cc * R_scale
  raxis_cs = raxis_cs * R_scale
  zaxis_cc = zaxis_cc * R_scale
  zaxis_cs = zaxis_cs * R_scale
  phiedge = phiedge * B_scale * R_scale * R_scale
  pres_scale = pres_scale * B_scale * B_scale
  curtor = curtor * B_scale * R_scale

  !------------------------------------------

  open(unit=iunit, file=trim(inputFilename)//"_scaled", action="write", iostat=istat)
  if (istat /= 0) then
     print *,"Error opening output file."
     stop
  end if

  call write_indata_namelist(iunit, istat)
  if (istat .ne. 0) then
     print *,"Error writing scaled namelist! istat=",istat
     stop
  else
     print *,"Successfully wrote scaled indata namelist."
  end if

  close(unit = iunit)
  print *,"Good bye!"

end program scale_vmec
