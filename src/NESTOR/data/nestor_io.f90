module nestor_io
  use stel_kinds, only: dp
IMPLICIT NONE

  character(len=255) :: input_extension
  character(len=255) :: mgrid_file
  real(dp), dimension(:), ALLOCATABLE :: extcur
  real(dp), dimension(:), ALLOCATABLE :: raxis
  real(dp), dimension(:), ALLOCATABLE :: zaxis
  real(dp), dimension(:), ALLOCATABLE :: xm
  real(dp), dimension(:), ALLOCATABLE :: xn
  real(dp), dimension(:), ALLOCATABLE :: rmnc
  real(dp), dimension(:), ALLOCATABLE :: zmns
  real(dp), dimension(:), ALLOCATABLE :: rmns
  real(dp), dimension(:), ALLOCATABLE :: zmnc
  real(dp), dimension(:), ALLOCATABLE :: wint

  integer :: nfp
  integer :: ntor
  integer :: mpol
  integer :: ntheta
  integer :: nzeta
  integer :: nextcur
  integer :: ier_flag
  integer :: ivac
  integer :: ivacskip
  integer :: mnmax
  integer :: vacuum_calls

  logical :: lasym

  real(dp) :: ctor
  real(dp) :: rbtor
  real(dp) :: signgs

  integer :: mnpd2_nestor
  real(dp), dimension(:), ALLOCATABLE :: amatsav_nestor
  real(dp), dimension(:), ALLOCATABLE :: bvecsav_nestor
  real(dp) :: bsubvvac_nestor

  CHARACTER(LEN=*), PARAMETER, DIMENSION(1) :: mn1dim = (/'mn_mode'/)
  CHARACTER(LEN=*), PARAMETER, DIMENSION(1) :: mnpotdim = (/'mn_mode_pot'/)
  CHARACTER(LEN=*), PARAMETER, DIMENSION(1) :: nzntdim = (/'nznt'/)
  CHARACTER(LEN=*), PARAMETER, DIMENSION(1) :: nzetadim = (/'nzeta'/)
  CHARACTER(LEN=*), PARAMETER, DIMENSION(1) :: nextcurim = (/'nextcur'/)
  CHARACTER(LEN=*), PARAMETER, DIMENSION(1) :: bvecsavdim =(/'mnpd2'/)
  CHARACTER(LEN=*), PARAMETER, DIMENSION(1) :: amatsavdim =(/'mnpd2_times_mnpd2'/)

  CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: r2dim = (/'mn_mode','radius '/)

  character(len=*), parameter :: &
     vn_vacuum_calls = 'vacuum_calls', &
     vn_ier_flag     = "ier_flag", &
     vn_mgrid        = "mgrid_file", &
     vn_inputext     = "input_extension", &
     vn_ivacskip     = "ivacskip", &
     vn_ivac         = "ivac", &
     vn_nfp          = "nfp", &
     vn_ntor         = "ntor", &
     vn_mpol         = "mpol", &
     vn_nzeta        = "nzeta", &
     vn_ntheta       = "ntheta", &
     vn_mnmax        = "mnmax", &
     vn_pmod         = "xm", &
     vn_tmod         = "xn", &
     vn_rmnc         = "rmnc", &
     vn_zmns         = "zmns", &
     vn_rmns         = "rmns", &
     vn_zmnc         = "zmnc", &
     vn_rbtor        = "rbtor", &
     vn_ctor         = "ctor", &
     vn_lasym        = "lasym", &
     vn_signgs       = "signgs", &
     vn_extcur       = "extcur", &
     vn_raxis_nestor = "raxis_nestor", &
     vn_zaxis_nestor = "zaxis_nestor", &
     vn_wint         = "wint"

  character(len=*), parameter :: &
     vn_bsqvac   = "bsqvac", &
     vn_mnpd     = "mnpd", &
     vn_xmpot    = "xmpot", &
     vn_xnpot    = "xnpot", &
     vn_potvac   = "potvac", &
     vn_brv      = "brv", &
     vn_bphiv    = "bphiv", &
     vn_bzv      = "bzv", &
     vn_bsubvvac = "bsubvvac", &
     vn_amatsav = "amatsav", &
     vn_bvecsav = "bvecsav", &
     vn_mnpd2   = "mnpd2"

  ! for debugging
  character(len=*), parameter :: &
     vn_r1b    = "r1b"   , &
     vn_rub    = "rub"   , &
     vn_rvb    = "rvb"   , &
     vn_z1b    = "z1b"   , &
     vn_zub    = "zub"   , &
     vn_zvb    = "zvb"   , &
     vn_ruu    = "ruu"   , &
     vn_ruv    = "ruv"   , &
     vn_rvv    = "rvv"   , &
     vn_zuu    = "zuu"   , &
     vn_zuv    = "zuv"   , &
     vn_zvv    = "zvv"   , &
     vn_guu_b  = "guu_b" , &
     vn_guv_b  = "guv_b" , &
     vn_gvv_b  = "gvv_b" , &
     vn_rzb2   = "rzb2"  , &
     vn_snr    = "snr"   , &
     vn_snv    = "snv"   , &
     vn_snz    = "snz"   , &
     vn_drv    = "drv"   , &
     vn_auu    = "auu"   , &
     vn_auv    = "auv"   , &
     vn_avv    = "avv"   , &
     vn_rcosuv = "rcosuv", &
     vn_rsinuv = "rsinuv", &
     vn_brad   = "brad", &
     vn_bphi   = "bphi", &
     vn_bz     = "bz", &
     vn_bexu   = "bexu", &
     vn_bexv   = "bexv", &
     vn_bexn   = "bexn", &
     vn_bexni  = "bexni", &
     vn_grpmn  = "grpmn", &
     vn_adp    = "adp"   , &
     vn_adm    = "adm"   , &
     vn_cma    = "cma"   , &
     vn_sqrtc  = "sqrtc" , &
     vn_sqrta  = "sqrta" , &
     vn_delt1u = "delt1u", &
     vn_azp1u  = "azp1u" , &
     vn_azm1u  = "azm1u" , &
     vn_cma11u = "cma11u", &
     vn_r1p    = "r1p"   , &
     vn_r1m    = "r1m"   , &
     vn_r0p    = "r0p"   , &
     vn_r0m    = "r0m"   , &
     vn_ra1p   = "ra1p"  , &
     vn_ra1m   = "ra1m"  , &
     vn_sqad1u = "sqad1u", &
     vn_sqad2u = "sqad2u", &
     vn_all_tlp = "all_tlp", &
     vn_all_tlm = "all_tlm", &
     vn_all_slp = "all_slp", &
     vn_all_slm = "all_slm", &
     vn_m_map = "m_map", &
     vn_n_map = "n_map", &
     vn_green = "green", &
     vn_greenp = "greenp", &
     vn_tanu = "tanu", &
     vn_tanv = "tanv", &
     vn_gstore = "gstore", &
     vn_grpmn_m_map = "grpmn_m_map", &
     vn_grpmn_n_map = "grpmn_n_map", &
     vn_imirr = "imirr", &
     vn_amatrix = "amatrix", &
     vn_potu = "potu", &
     vn_potv = "potv", &
     vn_bsubu = "bsubu", &
     vn_bsubv = "bsubv"

  CONTAINS


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
  mgrid_file = mgrid_file(1:mgrid_name_lens(1))

  CALL cdf_inquire(nvac, vn_inputext, dimlens)
  CALL cdf_read(nvac, vn_inputext,    input_extension(1:dimlens(1)))
  input_extension = input_extension(1:dimlens(1))

  CALL cdf_read(nvac, vn_vacuum_calls,  vacuum_calls)
  CALL cdf_read(nvac, vn_ivacskip,      ivacskip)
  CALL cdf_read(nvac, vn_ivac,          ivac)
  !print *, "read in nestor ivac=",ivac
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

  CALL cdf_read(nvac, vn_mnpd2, mnpd2_nestor)

  CALL cdf_inquire(nvac, vn_amatsav, dimlens)
  if (dimlens(1) .ne. mnpd2_nestor*mnpd2_nestor) then
     print *, "dimension mismatch in amatsav: shall be mnpd2*mnpd2=", &
        mnpd2_nestor*mnpd2_nestor," but is ",dimlens(1)
  end if
  ALLOCATE (amatsav_nestor(dimlens(1)), stat = ierror)

  CALL cdf_inquire(nvac, vn_bvecsav, dimlens)
  if (dimlens(1) .ne. mnpd2_nestor) then
     print *, "dimension mismatch in bvecsav: shall be mnpd2=", &
        mnpd2_nestor," but is ",dimlens(1)
  end if
  ALLOCATE (bvecsav_nestor(dimlens(1)), stat = ierror)

  CALL cdf_read(nvac, vn_amatsav,       amatsav_nestor)
  CALL cdf_read(nvac, vn_bvecsav,       bvecsav_nestor)

  CALL cdf_read(nvac, vn_bsubvvac,      bsubvvac_nestor)

  call cdf_close(nvac)

  if (ierror .ne. 0) then
     stop "problem reading NESTOR input file"
  end if

  ierror = 0
  CALL read_mgrid (trim(mgrid_file), extcur, nzeta, nfp, .false., ierror)
  if (ierror .ne. 0) then
     stop "could not read mgrid file"
  end if

end subroutine



subroutine write_nestor_outputs(vac_file, lasym, ivac, ier_flag)
  USE ezcdf
  use stel_kinds, only: dp
  !use vacmod, only: brv, bphiv, bzv, bsqvac, mnpd, xmpot, xnpot, potvac, mnpd2, bsubvvac, &
  !   amatsav, bvecsav
  use vacmod

  character(len=*), intent(in) :: vac_file
  logical, intent(in) :: lasym

  integer, intent(in) :: ivac, ier_flag

  integer :: istat_vac, nvac

  istat_vac = 0
  CALL cdf_open(nvac, trim(vac_file), 'w', istat_vac)
  IF (istat_vac .ne. 0) STOP 'Error writing output file of NESTOR'

  call cdf_define(nvac, vn_ivac, ivac)
  call cdf_define(nvac, vn_ier_flag, ier_flag)
  call cdf_define(nvac, vn_bsqvac  , bsqvac  )
  call cdf_define(nvac, vn_mnpd    , mnpd    )
  call cdf_define(nvac, vn_mnpd2    , mnpd2    )
  call cdf_define(nvac, vn_xmpot   , xmpot   )
  call cdf_define(nvac, vn_xnpot   , xnpot   )
  call cdf_define(nvac, vn_potvac  , potvac(1:mnpd2)  )
  call cdf_define(nvac, vn_brv     , brv     )
  call cdf_define(nvac, vn_bphiv   , bphiv   )
  call cdf_define(nvac, vn_bzv     , bzv     )
  call cdf_define(nvac, vn_bsubvvac, bsubvvac)
  call cdf_define(nvac, vn_amatsav, amatsav)
  call cdf_define(nvac, vn_bvecsav, bvecsav)

  ! below quantities are only written for debugging
  call cdf_define(nvac, vn_r1b   , r1b   )
  call cdf_define(nvac, vn_rub   , rub   )
  call cdf_define(nvac, vn_rvb   , rvb   )
  call cdf_define(nvac, vn_z1b   , z1b   )
  call cdf_define(nvac, vn_zub   , zub   )
  call cdf_define(nvac, vn_zvb   , zvb   )
  call cdf_define(nvac, vn_ruu   , ruu   )
  call cdf_define(nvac, vn_ruv   , ruv   )
  call cdf_define(nvac, vn_rvv   , rvv   )
  call cdf_define(nvac, vn_zuu   , zuu   )
  call cdf_define(nvac, vn_zuv   , zuv   )
  call cdf_define(nvac, vn_zvv   , zvv   )
  call cdf_define(nvac, vn_guu_b , guu_b )
  call cdf_define(nvac, vn_guv_b , guv_b )
  call cdf_define(nvac, vn_gvv_b , gvv_b )
  call cdf_define(nvac, vn_rzb2  , rzb2  )
  call cdf_define(nvac, vn_snr   , snr   )
  call cdf_define(nvac, vn_snv   , snv   )
  call cdf_define(nvac, vn_snz   , snz   )
  call cdf_define(nvac, vn_drv   , drv   )
  call cdf_define(nvac, vn_auu   , auu   )
  call cdf_define(nvac, vn_auv   , auv   )
  call cdf_define(nvac, vn_avv   , avv   )
  call cdf_define(nvac, vn_rcosuv, rcosuv)
  call cdf_define(nvac, vn_rsinuv, rsinuv)
  call cdf_define(nvac, vn_brad, brad)
  call cdf_define(nvac, vn_bphi, bphi)
  call cdf_define(nvac, vn_bz  , bz  )
  call cdf_define(nvac, vn_bexu, bexu)
  call cdf_define(nvac, vn_bexv, bexv)
  call cdf_define(nvac, vn_bexn, bexn)
  call cdf_define(nvac, vn_bexni, bexni)
  call cdf_define(nvac, vn_grpmn , grpmn )
  call cdf_define(nvac, vn_adp   , adp   )
  call cdf_define(nvac, vn_adm   , adm   )
  call cdf_define(nvac, vn_cma   , cma   )
  call cdf_define(nvac, vn_sqrtc , sqrtc )
  call cdf_define(nvac, vn_sqrta , sqrta )
  call cdf_define(nvac, vn_delt1u, delt1u)
  call cdf_define(nvac, vn_azp1u , azp1u )
  call cdf_define(nvac, vn_azm1u , azm1u )
  call cdf_define(nvac, vn_cma11u, cma11u)
  call cdf_define(nvac, vn_r1p   , r1p   )
  call cdf_define(nvac, vn_r1m   , r1m   )
  call cdf_define(nvac, vn_r0p   , r0p   )
  call cdf_define(nvac, vn_r0m   , r0m   )
  call cdf_define(nvac, vn_ra1p  , ra1p  )
  call cdf_define(nvac, vn_ra1m  , ra1m  )
  call cdf_define(nvac, vn_sqad1u, sqad1u)
  call cdf_define(nvac, vn_sqad2u, sqad2u)
  call cdf_define(nvac, vn_all_tlp, all_tlp)
  call cdf_define(nvac, vn_all_tlm, all_tlm)
  call cdf_define(nvac, vn_all_slp, all_slp)
  call cdf_define(nvac, vn_all_slm, all_slm)
  call cdf_define(nvac, vn_m_map, m_map_wrt)
  call cdf_define(nvac, vn_n_map, n_map_wrt)
  call cdf_define(nvac, vn_green, green)
  call cdf_define(nvac, vn_greenp, greenp)
  call cdf_define(nvac, vn_tanu, tanu)
  call cdf_define(nvac, vn_tanv, tanv)
  call cdf_define(nvac, vn_gstore, gstore)
  call cdf_define(nvac, vn_grpmn_m_map, grpmn_m_map_wrt)
  call cdf_define(nvac, vn_grpmn_n_map, grpmn_n_map_wrt)
  call cdf_define(nvac, vn_imirr, imirr)
  call cdf_define(nvac, vn_amatrix, amatrix)
  call cdf_define(nvac, vn_potu, potu)
  call cdf_define(nvac, vn_potv, potv)
  call cdf_define(nvac, vn_bsubu, bsubu)
  call cdf_define(nvac, vn_bsubv, bsubv)



  ! actually write data
  !print *, "write ivac=",ivac
  call cdf_write(nvac, vn_ivac,     ivac)
  call cdf_write(nvac, vn_ier_flag, ier_flag)
  call cdf_write(nvac, vn_bsqvac  , bsqvac  )
  call cdf_write(nvac, vn_mnpd    , mnpd    )
  call cdf_write(nvac, vn_mnpd2    , mnpd2    )
  call cdf_write(nvac, vn_xmpot   , xmpot   )
  call cdf_write(nvac, vn_xnpot   , xnpot   )
  call cdf_write(nvac, vn_potvac  , potvac(1:mnpd2)  )
  call cdf_write(nvac, vn_brv     , brv     )
  call cdf_write(nvac, vn_bphiv   , bphiv   )
  call cdf_write(nvac, vn_bzv     , bzv     )
  call cdf_write(nvac, vn_bsubvvac, bsubvvac)
  call cdf_write(nvac, vn_amatsav, amatsav)
  call cdf_write(nvac, vn_bvecsav, bvecsav)
  ! below quantities are only written for debugging
  call cdf_write(nvac, vn_r1b   , r1b   )
  call cdf_write(nvac, vn_rub   , rub   )
  call cdf_write(nvac, vn_rvb   , rvb   )
  call cdf_write(nvac, vn_z1b   , z1b   )
  call cdf_write(nvac, vn_zub   , zub   )
  call cdf_write(nvac, vn_zvb   , zvb   )
  call cdf_write(nvac, vn_ruu   , ruu   )
  call cdf_write(nvac, vn_ruv   , ruv   )
  call cdf_write(nvac, vn_rvv   , rvv   )
  call cdf_write(nvac, vn_zuu   , zuu   )
  call cdf_write(nvac, vn_zuv   , zuv   )
  call cdf_write(nvac, vn_zvv   , zvv   )
  call cdf_write(nvac, vn_guu_b , guu_b )
  call cdf_write(nvac, vn_guv_b , guv_b )
  call cdf_write(nvac, vn_gvv_b , gvv_b )
  call cdf_write(nvac, vn_rzb2  , rzb2  )
  call cdf_write(nvac, vn_snr   , snr   )
  call cdf_write(nvac, vn_snv   , snv   )
  call cdf_write(nvac, vn_snz   , snz   )
  call cdf_write(nvac, vn_drv   , drv   )
  call cdf_write(nvac, vn_auu   , auu   )
  call cdf_write(nvac, vn_auv   , auv   )
  call cdf_write(nvac, vn_avv   , avv   )
  call cdf_write(nvac, vn_rcosuv, rcosuv)
  call cdf_write(nvac, vn_rsinuv, rsinuv)
  call cdf_write(nvac, vn_brad, brad)
  call cdf_write(nvac, vn_bphi, bphi)
  call cdf_write(nvac, vn_bz  , bz  )
  call cdf_write(nvac, vn_bexu, bexu)
  call cdf_write(nvac, vn_bexv, bexv)
  call cdf_write(nvac, vn_bexn, bexn)
  call cdf_write(nvac, vn_bexni, bexni)
  call cdf_write(nvac, vn_grpmn , grpmn )
  call cdf_write(nvac, vn_adp   , adp   )
  call cdf_write(nvac, vn_adm   , adm   )
  call cdf_write(nvac, vn_cma   , cma   )
  call cdf_write(nvac, vn_sqrtc , sqrtc )
  call cdf_write(nvac, vn_sqrta , sqrta )
  call cdf_write(nvac, vn_delt1u, delt1u)
  call cdf_write(nvac, vn_azp1u , azp1u )
  call cdf_write(nvac, vn_azm1u , azm1u )
  call cdf_write(nvac, vn_cma11u, cma11u)
  call cdf_write(nvac, vn_r1p   , r1p   )
  call cdf_write(nvac, vn_r1m   , r1m   )
  call cdf_write(nvac, vn_r0p   , r0p   )
  call cdf_write(nvac, vn_r0m   , r0m   )
  call cdf_write(nvac, vn_ra1p  , ra1p  )
  call cdf_write(nvac, vn_ra1m  , ra1m  )
  call cdf_write(nvac, vn_sqad1u, sqad1u)
  call cdf_write(nvac, vn_sqad2u, sqad2u)
  call cdf_write(nvac, vn_all_tlp, all_tlp)
  call cdf_write(nvac, vn_all_tlm, all_tlm)
  call cdf_write(nvac, vn_all_slp, all_slp)
  call cdf_write(nvac, vn_all_slm, all_slm)
  call cdf_write(nvac, vn_m_map, m_map_wrt)
  call cdf_write(nvac, vn_n_map, n_map_wrt)
  call cdf_write(nvac, vn_green, green)
  call cdf_write(nvac, vn_greenp, greenp)
  call cdf_write(nvac, vn_tanu, tanu)
  call cdf_write(nvac, vn_tanv, tanv)
  call cdf_write(nvac, vn_gstore, gstore)
  call cdf_write(nvac, vn_grpmn_m_map, grpmn_m_map_wrt)
  call cdf_write(nvac, vn_grpmn_n_map, grpmn_n_map_wrt)
  call cdf_write(nvac, vn_imirr, imirr)
  call cdf_write(nvac, vn_amatrix, amatrix)
  call cdf_write(nvac, vn_potu, potu)
  call cdf_write(nvac, vn_potv, potv)
  call cdf_write(nvac, vn_bsubu, bsubu)
  call cdf_write(nvac, vn_bsubv, bsubv)


  CALL cdf_close(nvac)

  !print *, "write_nestor_outputs completed"
end subroutine






end module


subroutine write_nestor_inputs(vac_file,                               &
                  vacuum_calls, ier_flag, mgrid_file, input_extension, &
                  ivacskip, ivac, nfp, ntor, mpol, nzeta, ntheta,  &
                  mnmax, xm, xn, rmnc, zmns, rmns, zmnc,               &
                  rbtor, ctor, lasym, signgs, extcur_nestor,           &
                  raxis_nestor, zaxis_nestor, wint, nznt, amatsav, bvecsav, mnpd2, &
                  bsubvvac)

  USE ezcdf
  use stel_kinds, only: dp
  USE mgrid_mod, ONLY: nextcur
  use nestor_io, only: &
    mn1dim, mnpotdim, nzntdim, nzetadim, nextcurim, r2dim, &
    bvecsavdim, amatsavdim, &
    vn_vacuum_calls, &
    vn_ier_flag    , &
    vn_mgrid       , &
    vn_inputext    , &
    vn_ivacskip    , &
    vn_ivac        , &
    vn_nfp         , &
    vn_ntor        , &
    vn_mpol        , &
    vn_nzeta       , &
    vn_ntheta      , &
    vn_mnmax       , &
    vn_pmod        , &
    vn_tmod        , &
    vn_rmnc        , &
    vn_zmns        , &
    vn_rmns        , &
    vn_zmnc        , &
    vn_rbtor       , &
    vn_ctor        , &
    vn_lasym       , &
    vn_signgs      , &
    vn_extcur      , &
    vn_raxis_nestor, &
    vn_zaxis_nestor, &
    vn_wint, &
    vn_amatsav, &
    vn_bvecsav, &
    vn_mnpd2, &
    vn_bsubvvac

  character(len=*), intent(in) :: vac_file, mgrid_file, input_extension
  integer, intent(in) :: &
     vacuum_calls, ier_flag, ivacskip, ivac, nfp, ntor, mpol, nzeta, ntheta, mnmax, nznt, &
     mnpd2
  REAL(dp), intent(in) :: rbtor, ctor, signgs
  REAL(dp), DIMENSION(mnmax), INTENT(in) :: xm, xn, rmnc, rmns, zmns, zmnc
  LOGICAL, intent(in) :: lasym
  real(dp), dimension(nextcur), intent(in) :: extcur_nestor
  REAL(dp), DIMENSION(nzeta), INTENT(in) :: raxis_nestor, zaxis_nestor
  REAL(dp), DIMENSION(nznt), INTENT(in) :: wint
  real(dp), dimension(mnpd2), intent(in) :: bvecsav
  real(dp), dimension(mnpd2*mnpd2), intent(in) :: amatsav
  real(dp), intent(in) :: bsubvvac

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
  CALL cdf_define(nvac, vn_bvecsav,      bvecsav, dimname=bvecsavdim)
  CALL cdf_define(nvac, vn_amatsav,      amatsav, dimname=amatsavdim)
  CALL cdf_define(nvac, vn_mnpd2,         mnpd2)
  CALL cdf_define(nvac, vn_bsubvvac,      bsubvvac)

  ! actually write data
  CALL cdf_write(nvac, vn_vacuum_calls,  vacuum_calls)
  CALL cdf_write(nvac, vn_ier_flag,      ier_flag)
  CALL cdf_write(nvac, vn_mgrid,         trim(mgrid_file))
  CALL cdf_write(nvac, vn_inputext,      trim(input_extension))
  CALL cdf_write(nvac, vn_ivacskip,      ivacskip)
  CALL cdf_write(nvac, vn_ivac,          ivac)
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
  CALL cdf_write(nvac, vn_bvecsav,      bvecsav)
  CALL cdf_write(nvac, vn_amatsav,      amatsav)
  CALL cdf_write(nvac, vn_mnpd2,         mnpd2)
  CALL cdf_write(nvac, vn_bsubvvac,      bsubvvac)

  CALL cdf_close(nvac)
end subroutine

subroutine read_nestor_outputs(vac_file, ier_flag, ivac)
  USE ezcdf
  use stel_kinds, only: dp
  use vacmod, only: brv, bphiv, bzv, bsqvac, mnpd, xmpot, xnpot, potvac, mnpd2, bsubvvac, &
     amatsav, bvecsav
  use nestor_io, only: &
    vn_ier_flag    , &
    vn_ivac        , &
    vn_bsqvac  , &
    vn_mnpd    , &
    vn_xmpot   , &
    vn_xnpot   , &
    vn_potvac  , &
    vn_brv     , &
    vn_bphiv   , &
    vn_bzv     , &
    vn_bsubvvac, &
    vn_amatsav , &
    vn_bvecsav, &
    vn_mnpd2

  character(len=*), intent(in) :: vac_file
  integer, intent(out) :: ier_flag
  integer, intent(out) :: ivac

  integer :: istat_vac, nvac

  istat_vac = 0
  CALL cdf_open(nvac, trim(vac_file), 'r', istat_vac)
  IF (istat_vac .ne. 0) STOP 'Error writing output file of NESTOR'

  ! read data; assume allocated since allocate_nestor was also called in main VMEC executable
  call cdf_read(nvac, vn_ivac,     ivac)
  call cdf_read(nvac, vn_ier_flag, ier_flag)
  call cdf_read(nvac, vn_bsqvac  , bsqvac  )
  call cdf_read(nvac, vn_mnpd    , mnpd    )
  call cdf_read(nvac, vn_mnpd2    , mnpd2    )
  call cdf_read(nvac, vn_xmpot   , xmpot   )
  call cdf_read(nvac, vn_xnpot   , xnpot   )
  call cdf_read(nvac, vn_potvac  , potvac(1:mnpd2)  )
  call cdf_read(nvac, vn_brv     , brv     )
  call cdf_read(nvac, vn_bphiv   , bphiv   )
  call cdf_read(nvac, vn_bzv     , bzv     )
  call cdf_read(nvac, vn_bsubvvac, bsubvvac)
  call cdf_read(nvac, vn_amatsav, amatsav)
  call cdf_read(nvac, vn_bvecsav, bvecsav)

  CALL cdf_close(nvac)

end subroutine
