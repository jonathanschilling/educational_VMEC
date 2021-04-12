module nestor_io
  use stel_kinds, only: dp
IMPLICIT NONE

  character(len=255) :: input_extension
  character(len=255) :: mgrid_file
  real(dp), dimension(:), ALLOCATABLE :: extcur
  real(dp), dimension(:), ALLOCATABLE :: raxis, zaxis
  real(dp), dimension(:), ALLOCATABLE :: xm, xn, rmnc, zmns, rmns, zmnc
  real(dp), dimension(:), ALLOCATABLE :: wint

  integer :: nfp, ntor, mpol, ntheta, nzeta, nextcur
  logical :: lasym
  integer  :: ier_flag, ivac, ivacskip, mnmax
  real(dp) :: ctor, rbtor, signgs



  CHARACTER(LEN=*), PARAMETER, DIMENSION(1) :: &
     mn1dim = (/'mn_mode'/), &
     mnpotdim = (/'mn_mode_pot'/), &
     nzntdim = (/'nznt'/), &
     nzetadim = (/'nzeta'/), &
     nextcurim = (/'nextcur'/)

  CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: &
     r2dim = (/'mn_mode','radius '/)

  character(len=*), parameter :: &
     vn_vacuum_calls = 'vacuum_calls', &
     vn_ier_flag = "ier_flag", &
     vn_mgrid = "mgrid_file", &
     vn_inputext = "input_extension", &
     vn_ivacskip = "ivacskip", &
     vn_ivac = "ivac", &
     vn_ns = "ns", &
     vn_nfp = "nfp", &
     vn_ntor = "ntor", &
     vn_mpol = "mpol", &
     vn_nzeta = "nzeta", &
     vn_ntheta = "ntheta", &
     vn_mnmax = "mnmax", &
     vn_pmod = "xm", &
     vn_tmod = "xn", &
     vn_rmnc = "rmnc", &
     vn_zmns = "zmns", &
     vn_rmns = "rmns", &
     vn_zmnc = "zmnc", &
     vn_rbtor = "rbtor", &
     vn_ctor = "ctor", &
     vn_lasym = "lasym", &
     vn_signgs = "signgs", &
     vn_extcur = "extcur", &
     vn_raxis_nestor = "raxis_nestor", &
     vn_zaxis_nestor = "zaxis_nestor", &
     vn_wint = "wint"

  CONTAINS


subroutine write_nestor_inputs(vac_file,                               &
                  vacuum_calls, ier_flag, mgrid_file, input_extension, &
                  ivacskip, ivac, ns, nfp, ntor, mpol, nzeta, ntheta,  &
                  mnmax, xm, xn, rmnc, zmns, rmns, zmnc,               &
                  rbtor, ctor, lasym, signgs, extcur_nestor,           &
                  raxis_nestor, zaxis_nestor, wint, nznt)

  USE ezcdf
  use stel_kinds, only: dp
  USE mgrid_mod, ONLY: nextcur

  character(len=*), intent(in) :: vac_file, mgrid_file, input_extension
  integer, intent(in) :: &
     vacuum_calls, ier_flag, ivacskip, ivac, ns, nfp, ntor, mpol, nzeta, ntheta, mnmax, nznt
  REAL(dp), intent(in) :: rbtor, ctor, signgs
  REAL(dp), DIMENSION(mnmax), INTENT(in) :: xm, xn, rmnc, rmns, zmns, zmnc
  LOGICAL, intent(in) :: lasym
  real(dp), dimension(nextcur), intent(in) :: extcur_nestor
  REAL(dp), DIMENSION(nzeta), INTENT(in) :: raxis_nestor, zaxis_nestor
  REAL(dp), DIMENSION(nznt), INTENT(in) :: wint

  integer :: istat_vac, nvac

  istat_vac = 0
  CALL cdf_open(nvac, trim(vac_file), 'w', istat_vac)
  IF (istat_vac .ne. 0) STOP 'Error writing input file for NESTOR'

  ! definitions for variables
  CALL cdf_define(nvac, vn_vacuum_calls, vacuum_calls)
  CALL cdf_define(nvac, vn_ier_flag,     ier_flag)
  CALL cdf_define(nvac, vn_mgrid,        trim(mgrid_file))
  CALL cdf_define(nvac, vn_inputext,     trim(input_extension))
  CALL cdf_define(nvac, vn_ivacskip,     ivacskip)
  CALL cdf_define(nvac, vn_ivac,         ivac)
  CALL cdf_define(nvac, vn_ns,           ns)
  CALL cdf_define(nvac, vn_nfp,          nfp)
  CALL cdf_define(nvac, vn_ntor,         ntor)
  CALL cdf_define(nvac, vn_mpol,         mpol)
  CALL cdf_define(nvac, vn_nzeta,        nzeta)
  CALL cdf_define(nvac, vn_ntheta,       ntheta)
  CALL cdf_define(nvac, vn_mnmax,        mnmax)
  CALL cdf_define(nvac, vn_pmod,         xm,   dimname=mn1dim)
  CALL cdf_define(nvac, vn_tmod,         xn,   dimname=mn1dim)
  CALL cdf_define(nvac, vn_rmnc,         rmnc, dimname=r2dim)
  CALL cdf_define(nvac, vn_zmns,         zmns, dimname=r2dim)
  if (lasym) then
     CALL cdf_define(nvac, vn_rmns,      rmns, dimname=r2dim)
     CALL cdf_define(nvac, vn_zmnc,      zmnc, dimname=r2dim)
  end if ! lasym
  CALL cdf_define(nvac, vn_rbtor,        rbtor)
  CALL cdf_define(nvac, vn_ctor,         ctor)
  CALL cdf_define(nvac, vn_lasym,        lasym)
  CALL cdf_define(nvac, vn_signgs,       signgs)
  CALL cdf_define(nvac, vn_extcur,       extcur_nestor, dimname=nextcurim)
  CALL cdf_define(nvac, vn_raxis_nestor, raxis_nestor, dimname=nzetadim)
  CALL cdf_define(nvac, vn_zaxis_nestor, zaxis_nestor, dimname=nzetadim)
  CALL cdf_define(nvac, vn_wint,         wint, dimname=nzntdim)

  ! actually write data
  CALL cdf_write(nvac, vn_vacuum_calls,  vacuum_calls)
  CALL cdf_write(nvac, vn_ier_flag,      ier_flag)
  CALL cdf_write(nvac, vn_mgrid,         trim(mgrid_file))
  CALL cdf_write(nvac, vn_inputext,      trim(input_extension))
  CALL cdf_write(nvac, vn_ivacskip,      ivacskip)
  CALL cdf_write(nvac, vn_ivac,          ivac)
  CALL cdf_write(nvac, vn_ns,            ns)
  CALL cdf_write(nvac, vn_nfp,           nfp)
  CALL cdf_write(nvac, vn_ntor,          ntor)
  CALL cdf_write(nvac, vn_mpol,          mpol)
  CALL cdf_write(nvac, vn_nzeta,         nzeta)
  CALL cdf_write(nvac, vn_ntheta,        ntheta)
  CALL cdf_write(nvac, vn_mnmax,         mnmax)
  CALL cdf_write(nvac, vn_pmod,          xm)
  CALL cdf_write(nvac, vn_tmod,          xn)
  CALL cdf_write(nvac, vn_rmnc,          rmnc)
  CALL cdf_write(nvac, vn_zmns,          zmns)
  if (lasym) then
     CALL cdf_write(nvac, vn_rmns,       rmns)
     CALL cdf_write(nvac, vn_zmnc,       zmnc)
  end if ! lasym
  CALL cdf_write(nvac, vn_rbtor,         rbtor)
  CALL cdf_write(nvac, vn_ctor,          ctor)
  CALL cdf_write(nvac, vn_lasym,         lasym)
  CALL cdf_write(nvac, vn_signgs,        signgs)
  CALL cdf_write(nvac, vn_extcur,        extcur_nestor)
  CALL cdf_write(nvac, vn_raxis_nestor,  raxis_nestor)
  CALL cdf_write(nvac, vn_zaxis_nestor,  zaxis_nestor)
  CALL cdf_write(nvac, vn_wint,          wint)

  CALL cdf_close(nvac)
end subroutine


subroutine read_nestor_inputs(vac_file)
  use ezcdf
  use stel_kinds, only: dp
  use mgrid_mod, only: read_mgrid

  character(len=*), intent(in) :: vac_file

  integer :: istat_vac, nvac, ierror
  INTEGER, DIMENSION(3)   :: dimlens, mgrid_name_lens

  istat_vac = 0
  CALL cdf_open(nvac, trim(vac_file), 'r', istat_vac)
  IF (istat_vac .ne. 0) STOP 'Error reading input file of NESTOR'

  CALL cdf_read(nvac, vn_ier_flag,      ier_flag)

  CALL cdf_inquire(nvac, vn_mgrid, mgrid_name_lens)
  CALL cdf_read(nvac, vn_mgrid,    mgrid_file(1:mgrid_name_lens(1)))

  CALL cdf_inquire(nvac, vn_inputext, dimlens)
  CALL cdf_read(nvac, vn_inputext,    input_extension(1:dimlens(1)))

  CALL cdf_read(nvac, vn_ivacskip,      ivacskip)
  CALL cdf_read(nvac, vn_ivac,          ivac)
  CALL cdf_read(nvac, vn_nfp,           nfp)
  CALL cdf_read(nvac, vn_ntor,          ntor)
  CALL cdf_read(nvac, vn_mpol,          mpol)
  CALL cdf_read(nvac, vn_nzeta,         nzeta)
  CALL cdf_read(nvac, vn_ntheta,        ntheta)
  CALL cdf_read(nvac, vn_mnmax,         mnmax)
  CALL cdf_read(nvac, vn_lasym,         lasym)
  CALL cdf_read(nvac, vn_rbtor,         rbtor)
  CALL cdf_read(nvac, vn_ctor,          ctor)
  CALL cdf_read(nvac, vn_signgs,        signgs)

  CALL cdf_inquire(nvac, vn_extcur, dimlens)
  ALLOCATE (extcur(dimlens(1)), stat = ierror)

  CALL cdf_inquire(nvac, vn_raxis_nestor, dimlens)
  ALLOCATE (raxis(dimlens(1)), stat = ierror)
  CALL cdf_inquire(nvac, vn_zaxis_nestor, dimlens)
  ALLOCATE (zaxis(dimlens(1)), stat = ierror)

  CALL cdf_inquire(nvac, vn_pmod, dimlens)
  ALLOCATE (xm(dimlens(1)), stat = ierror)
  CALL cdf_inquire(nvac, vn_tmod, dimlens)
  ALLOCATE (xn(dimlens(1)), stat = ierror)

  CALL cdf_inquire(nvac, vn_rmnc, dimlens)
  ALLOCATE (rmnc(dimlens(1)), stat = ierror)
  CALL cdf_inquire(nvac, vn_zmns, dimlens)
  ALLOCATE (zmns(dimlens(1)), stat = ierror)

  CALL cdf_inquire(nvac, vn_wint, dimlens)
  ALLOCATE (wint(dimlens(1)), stat = ierror)

  CALL cdf_read(nvac, vn_extcur,        extcur)
  CALL cdf_read(nvac, vn_raxis_nestor,  raxis)
  CALL cdf_read(nvac, vn_zaxis_nestor,  zaxis)

  CALL cdf_read(nvac, vn_pmod,          xm)
  CALL cdf_read(nvac, vn_tmod,          xn)
  CALL cdf_read(nvac, vn_rmnc,          rmnc)
  CALL cdf_read(nvac, vn_zmns,          zmns)

  if (lasym) then
     CALL cdf_inquire(nvac, vn_rmns, dimlens)
     ALLOCATE (rmns(dimlens(1)), stat = ierror)
     CALL cdf_inquire(nvac, vn_zmnc, dimlens)
     ALLOCATE (zmnc(dimlens(1)), stat = ierror)

     CALL cdf_read(nvac, vn_rmns,       rmns)
     CALL cdf_read(nvac, vn_zmnc,       zmnc)
  end if ! lasym

  CALL cdf_read(nvac, vn_wint,          wint)

  call cdf_close(nvac)

  if (ierror .ne. 0) then
     stop "problem reading NESTOR input file"
  end if

  ierror = 0
  CALL read_mgrid (mgrid_file(1:mgrid_name_lens(1)), extcur, nzeta, nfp, .true., ierror)
  if (ierror .ne. 0) then
     stop "could not read mgrid file"
  end if

end subroutine



subroutine write_nestor_outputs(vac_file, vacuum_calls, lasym, &
    ivac, ier_flag)

  USE ezcdf
  use stel_kinds, only: dp

  character(len=*), intent(in) :: vac_file
  integer, intent(in) :: vacuum_calls
  logical, intent(in) :: lasym

  integer, intent(in) :: ivac, ier_flag

  integer :: istat_vac, nvac

  istat_vac = 0
  CALL cdf_open(nvac, trim(vac_file), 'w', istat_vac)
  IF (istat_vac .ne. 0) STOP 'Error writing output file of NESTOR'

  call cdf_define(nvac, vn_ivac, ivac)
  call cdf_define(nvac, vn_ier_flag, ier_flag)


  ! actually write data
  call cdf_write(nvac, vn_ivac, ivac)
  call cdf_write(nvac, vn_ier_flag, ier_flag)




  CALL cdf_close(nvac)

end subroutine



end module
