!> \file
MODULE vmec_input

  USE vparams, ONLY: rprec, dp, mpol1d, ntord, ndatafmax
  USE vsvd0

  IMPLICIT NONE

  INTEGER, PARAMETER :: mpol_default = 6
  INTEGER, PARAMETER :: ntor_default = 0
  INTEGER, PARAMETER :: ns_default   = 31
  INTEGER, PARAMETER :: niter_default   = 100
  REAL(rprec), PARAMETER :: ftol_default = 1.E-10_dp

  INTEGER                           :: nfp
  INTEGER                           :: ncurr
  INTEGER                           :: nstep
  INTEGER                           :: nvacskip
  INTEGER                           :: mpol
  INTEGER                           :: ntor
  INTEGER                           :: ntheta
  INTEGER                           :: nzeta
  INTEGER                           :: mfilter_fbdy
  INTEGER                           :: nfilter_fbdy

  INTEGER,     DIMENSION(100)       :: ns_array
  INTEGER,     DIMENSION(100)       :: niter_array
  REAL(rprec), DIMENSION(100)       :: ftol_array

  REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) :: rbc
  REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) :: zbs
  REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) :: rbs
  REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) :: zbc
  REAL(rprec)                       :: curtor
  REAL(rprec)                       :: delt
  REAL(rprec)                       :: tcon0
  REAL(rprec)                       :: gamma
  REAL(rprec)                       :: bloat
  REAL(rprec)                       :: pres_scale
  REAL(rprec)                       :: spres_ped !< value of s beyond which pressure profile is flat (pedestal)
  REAL(rprec)                       :: phiedge   !< value of real toroidal flux at plasma edge (s=1)
  REAL(rprec), DIMENSION(0:20)      :: am        !< array of coefficients in phi-series for mass (NWT/m**2)
  REAL(rprec), DIMENSION(0:20)      :: ai        !< array of coefficients in phi-series for iota (ncurr=0)
  REAL(rprec), DIMENSION(0:20)      :: ac        !< array of coefficients in phi-series for the quantity d(Icurv)/ds = toroidal
                                                 !< current density * Vprime, so Icurv(s) = Itor(s) (used for ncurr=1)
  REAL(rprec), DIMENSION(1:20)      :: aphi
  CHARACTER(len=20)                 :: pcurr_type
  CHARACTER(len=20)                 :: piota_type
  CHARACTER(len=20)                 :: pmass_type
  REAL(rprec), DIMENSION(ndatafmax) :: am_aux_s
  REAL(rprec), DIMENSION(ndatafmax) :: am_aux_f
  REAL(rprec), DIMENSION(ndatafmax) :: ai_aux_s
  REAL(rprec), DIMENSION(ndatafmax) :: ai_aux_f
  REAL(rprec), DIMENSION(ndatafmax) :: ac_aux_s
  REAL(rprec), DIMENSION(ndatafmax) :: ac_aux_f

  REAL(rprec), DIMENSION(0:ntord)   :: raxis_cc
  REAL(rprec), DIMENSION(0:ntord)   :: raxis_cs
  REAL(rprec), DIMENSION(0:ntord)   :: zaxis_cc
  REAL(rprec), DIMENSION(0:ntord)   :: zaxis_cs
  REAL(rprec), DIMENSION(nigroup)   :: extcur
  LOGICAL                           :: lfreeb
  LOGICAL                           :: lasym
  
  ! switch between implementations of NESTOR:
  ! vac1 (magnetic scalar potential, both Stellarator and Tokamak)
  ! vac2/vac3 (surface current density, Stellarator/Tokamak)
  ! Only vac1 is available in educational_VMEC.
  ! Therefore, this flag is ignored here.
  integer                           :: vac_1_2

  ! RESET lbsubs DEFAULT FLAG TO FALSE TO CAPTURE CURRENT SHEETS!
  ! LOGICAL, PARAMETER :: lbsubs = .false. ! False to use (correct)  bsubs calculation (from metrics)
  !                                        ! True  to use (modified) bsubs calculation (from mag. diff. eq.)
  ! lbsubs is now a namelist input variable, so user can change.
  ! LOGICAL, PARAMETER :: lbsubs = .true.  ! True  to use NEW bsubs calculation (from mag. diff. eq.)
  !                                        ! False to use OLD bsubs calculation (from metrics)
  LOGICAL                           :: lbsubs

  CHARACTER(len=200)                :: mgrid_file
  CHARACTER(len=100)                :: input_extension

  ! debugging output

  !> maximum number of iterations for which to dump data
  integer :: max_dump                 =    2

  !> individual flags to control debug output loosely related to
  !> similarly-named routines (checkpoints along iterations)
  logical :: dump_add_fluxes          = .false.
  logical :: dump_metric              = .false.
  logical :: dump_volume              = .false.
  logical :: dump_bcontrav            = .false.
  logical :: dump_bcov                = .false.
  logical :: dump_lambda_forces       = .false.
  logical :: dump_bcov_full           = .false.
  logical :: dump_precondn            = .false.
  logical :: dump_forceNorms_tcon     = .false.
  logical :: dump_lulv_comb           = .false.
  logical :: dump_calc_fbal           = .false.
  logical :: dump_evolve              = .false.
  logical :: dump_fixaray             = .false.
  logical :: dump_spectral_constraint = .false.
  logical :: dump_forces              = .false.
  logical :: dump_funct3d_geometry    = .false.
  logical :: dump_constraint_force    = .false.
  logical :: dump_guess_axis          = .false.
  logical :: dump_interp              = .false.
  logical :: dump_jacobian            = .false.
  logical :: dump_lamcal              = .false.
  logical :: dump_profil1d            = .false.
  logical :: dump_profil3d            = .false.
  logical :: dump_readin_boundary     = .false.
  logical :: dump_phys_gc             = .false.
  logical :: dump_fsq                 = .false.
  logical :: dump_scale_m1            = .false.
  logical :: dump_scalfor_out         = .false.
  logical :: dump_fsq1                = .false.
  logical :: dump_scalfor_R           = .false.
  logical :: dump_scalfor_Z           = .false.
  logical :: dump_symforce            = .false.
  logical :: dump_tomnsps             = .false.
  logical :: dump_tomnspa             = .false.
  logical :: dump_multigrid_result    = .false.
  logical :: dump_rbsq                = .false.

  ! debugging output flags for NESTOR
  logical :: dump_vac1n_vacuum  = .false.
  logical :: dump_vac1n_precal  = .false.
  logical :: dump_vac1n_surface = .false.
  logical :: dump_vac1n_bextern = .false.
  logical :: dump_vac1n_analyt  = .false.
  logical :: dump_vac1n_greenf  = .false.
  logical :: dump_vac1n_fourp   = .false.
  logical :: dump_vac1n_fouri   = .false.
  logical :: dump_vac1n_bsqvac  = .false.
  


  NAMELIST /indata/ &
     mgrid_file,    &
     nfp,           &
     ncurr,         &
     nstep,         &
     nvacskip,      &
     delt,          &
     gamma,         &
     am,            &
     ai,            &
     ac,            &
     aphi,          &
     pcurr_type,    &
     pmass_type,    &
     piota_type,    &
     bloat,         &
     am_aux_s,      &
     am_aux_f,      &
     ai_aux_s,      &
     ai_aux_f,      &
     ac_aux_s,      &
     ac_aux_f,      &
     rbc,           &
     zbs,           &
     rbs,           &
     zbc,           &
     spres_ped,     &
     pres_scale,    &
     raxis_cc,      &
     zaxis_cs,      &
     raxis_cs,      &
     zaxis_cc,      &
     mpol,          &
     ntor,          &
     ntheta,        &
     nzeta,         &
     mfilter_fbdy,  &
     nfilter_fbdy,  &
     niter_array,   &
     ns_array,      &
     ftol_array,    &
     tcon0,         &
     curtor,        &
     extcur,        &
     phiedge,       &
     lfreeb,        &
     lasym,         &
     lbsubs,        &
     max_dump                , &
     dump_add_fluxes         , &
     dump_metric             , &
     dump_volume             , &
     dump_bcontrav           , &
     dump_bcov               , &
     dump_lambda_forces      , &
     dump_bcov_full          , &
     dump_precondn           , &
     dump_forceNorms_tcon    , &
     dump_lulv_comb          , &
     dump_calc_fbal          , &
     dump_evolve             , &
     dump_fixaray            , &
     dump_spectral_constraint, &
     dump_forces             , &
     dump_funct3d_geometry   , &
     dump_constraint_force   , &
     dump_guess_axis         , &
     dump_interp             , &
     dump_jacobian           , &
     dump_lamcal             , &
     dump_profil1d           , &
     dump_profil3d           , &
     dump_readin_boundary    , &
     dump_phys_gc            , &
     dump_fsq                , &
     dump_scale_m1           , &
     dump_scalfor_out        , &
     dump_fsq1               , &
     dump_scalfor_R          , &
     dump_scalfor_Z          , &
     dump_symforce           , &
     dump_tomnsps            , &
     dump_tomnspa            , &
     dump_multigrid_result   , &
     dump_rbsq               , &
     vac_1_2, &
     dump_vac1n_vacuum,  &
     dump_vac1n_precal,  &
     dump_vac1n_surface, &
     dump_vac1n_bextern, &
     dump_vac1n_analyt,  &
     dump_vac1n_greenf,  &
     dump_vac1n_fourp,   &
     dump_vac1n_fouri,   &
     dump_vac1n_bsqvac
     

CONTAINS

SUBROUTINE read_indata_namelist (iunit, istat)
  INTEGER, intent(in)    :: iunit
  INTEGER, intent(inout) :: istat

  character(len=1000) :: line

  ! INITIALIZATIONS
  gamma = 0
  spres_ped = 1
  mpol = mpol_default
  ntor = ntor_default
  ntheta = 0
  nzeta = 0

     ns_array =  0;    ns_array(1) =    ns_default
   ftol_array =  0;  ftol_array(1) =  ftol_default
  niter_array = -1; niter_array(1) = niter_default

  bloat = 1
  rbc = 0
  rbs = 0
  zbs = 0
  zbc = 0
  nfp = 1
  ncurr = 0
  nstep = 10
  nvacskip = 1
  delt = 1

  am = 0
  ai = 0
  ac = 0

  aphi = 0
  aphi(1) = 1

  pres_scale = 1

  raxis_cc = 0
  zaxis_cs = 0
  raxis_cs = 0
  zaxis_cc = 0

  mfilter_fbdy = -1
  nfilter_fbdy = -1

  tcon0 = 1
  curtor = 0
  extcur = 0
  phiedge = 1

  mgrid_file = 'NONE'

  lfreeb = .true.
  lasym = .false.
  lbsubs = .false.

  pcurr_type = 'power_series'
  piota_type = 'power_series'
  pmass_type = 'power_series'

  am_aux_s(:) = -1
  ac_aux_s(:) = -1
  ai_aux_s(:) = -1

  READ (iunit, nml=indata, iostat=istat)

  if (istat .ne. 0) then
    ! help to debug invalid inputs:
    ! re-read last line that lead to error and print it to screen
    backspace(iunit)
    read(iunit,fmt='(A)') line
    write(*,'(A)') 'Invalid line in namelist: '//trim(line)
  end if

END SUBROUTINE read_indata_namelist


SUBROUTINE write_indata_namelist (iunit, istat)
  IMPLICIT NONE

  INTEGER, INTENT(in) :: iunit
  INTEGER, INTENT(inout) :: istat
  INTEGER :: iftol,i,n,m
  INTEGER, DIMENSION(1) :: ins
  CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
  CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
  CHARACTER(LEN=*), PARAMETER :: outint1 = "(2X,A,1X,'=',1X,I1.1)"
  CHARACTER(LEN=*), PARAMETER :: outint2 = "(2X,A,1X,'=',1X,I2.2)"
  CHARACTER(LEN=*), PARAMETER :: outint3 = "(2X,A,1X,'=',1X,I3.3)"
  CHARACTER(LEN=*), PARAMETER :: outint4 = "(2X,A,1X,'=',1X,I4.4)"
  CHARACTER(LEN=*), PARAMETER :: outint5 = "(2X,A,1X,'=',1X,I5.5)"
  CHARACTER(LEN=*), PARAMETER :: outint6 = "(2X,A,1X,'=',1X,I6.6)"
  CHARACTER(LEN=*), PARAMETER :: outflt="(2X,A,1X,'=',1X,ES22.12E3)"
  CHARACTER(LEN=*), PARAMETER :: outexp="(2X,A,1X,'=',1X,ES22.12E3)"
  IF (istat < 0) RETURN
  WRITE(iunit,'(A)') '!----- Runtime Parameters -----'
  WRITE(iunit,'(A)') '&INDATA'
  WRITE(iunit,outflt) 'DELT',delt
  WRITE(iunit,outint) 'NSTEP',nstep
  WRITE(iunit,outflt) 'TCON0',tcon0
  ins = MAXLOC(ns_array)
  WRITE(iunit,'(a,(1p,4i14))') '  NS_ARRAY =    ',(ns_array(i), i=1,ins(1))
  iftol = 1
  DO WHILE(ftol_array(iftol).ne.0 .and. iftol.lt.100)
     iftol = iftol + 1
  END DO
  WRITE(iunit,'(a,(1p,4e14.6))')'  FTOL_ARRAY =  ',(ftol_array(i), i=1,iftol - 1)
  ins = MINLOC(niter_array)
  IF (ins(1) > 1) WRITE(iunit,'(a,(1p,4i14))') '  NITER_ARRAY = ',(niter_array(i), i=1,ins(1)-1)
  WRITE(iunit,'(A)') '!----- Grid Parameters -----'
  WRITE(iunit,outboo) 'LASYM',lasym
  WRITE(iunit,outint4) 'NFP',nfp
  WRITE(iunit,outint4) 'MPOL',mpol
  WRITE(iunit,outint4) 'NTOR',ntor
  WRITE(iunit,outflt) 'PHIEDGE',phiedge
  WRITE(iunit,'(A)') '!----- Free Boundary Parameters -----'
  WRITE(iunit,outboo) 'LFREEB',lfreeb
  IF (lfreeb) THEN
     WRITE (iunit, '(2x,3a)') "MGRID_FILE = '",TRIM(mgrid_file),"'"
     WRITE(iunit,outint4) 'NZETA',nzeta
     DO n=1,SIZE(extcur)
        IF (extcur(n) == 0) CYCLE
        WRITE(iunit,'(2X,A,I3.3,A,ES22.12E3)') 'EXTCUR(',n,') = ',extcur(n)
     END DO
     WRITE(iunit,outint4) 'NVACSKIP',nvacskip
  END IF
  WRITE(iunit,'(A)') '!----- Pressure Parameters -----'
  WRITE(iunit,outflt) 'GAMMA',gamma
  WRITE(iunit,outflt) 'BLOAT',bloat
  WRITE(iunit,outflt) 'SPRES_PED',spres_ped
  WRITE(iunit,outflt) 'PRES_SCALE',pres_scale
  WRITE(iunit,'(2x,3a)') "PMASS_TYPE = '",TRIM(pmass_type),"'"
  WRITE(iunit,'(a,(1p,4e22.14))')'  AM = ', (am(i-1), i=1,SIZE(am))
  i = minloc(am_aux_s(2:),DIM=1)
  IF (i > 4) THEN
     WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AM_AUX_S = ', (am_aux_s(n), n=1,i)
     WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AM_AUX_F = ', (am_aux_f(n), n=1,i)
  END IF

  WRITE(iunit,'(A)') '!----- Current/Iota Parameters -----'
  WRITE(iunit,outexp) 'CURTOR',curtor
  WRITE(iunit,outint) 'NCURR',ncurr
  WRITE (iunit, '(2x,3a)') "PIOTA_TYPE = '",TRIM(piota_type),"'"
  WRITE (iunit,'(a,(1p,4e22.14))') '  AI = ',(ai(n-1), n=1,SIZE(ai))
  i = minloc(ai_aux_s(2:),DIM=1)
  IF (i > 4) THEN
     WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AI_AUX_S = ', (ai_aux_s(n), n=1,i)
     WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AI_AUX_F = ', (ai_aux_f(n), n=1,i)
  END IF
  WRITE (iunit, '(2x,3a)') "PCURR_TYPE = '",TRIM(pcurr_type),"'"
  WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AC = ',(ac(n-1), n=1,SIZE(ac))
  i = minloc(ac_aux_s(2:),DIM=1)
  IF (i > 4) THEN
     WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AC_AUX_S = ', (ac_aux_s(n), n=1,i)
     WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AC_AUX_F = ', (ac_aux_f(n), n=1,i)
  END IF

  WRITE(iunit,'(A)') '!----- Axis Parameters ----- '
  WRITE (iunit,'(a,(1p,4e22.14))') '  RAXIS_CC = ',(raxis_cc(n), n=0,ntor)
  IF (lasym) then
     WRITE (iunit,'(a,(1p,4ES22.12E3))') '  RAXIS_CS = ',(raxis_cs(n), n=0,ntor)
     WRITE (iunit,'(a,(1p,4ES22.12E3))') '  ZAXIS_CC = ',(zaxis_cc(n), n=0,ntor)
  end if
  WRITE (iunit,'(a,(1p,4ES22.12E3))') '  ZAXIS_CS = ',(zaxis_cs(n), n=0,ntor)

  WRITE(iunit,'(A)') '!----- Boundary Parameters -----'
  DO m = 0, mpol - 1
     DO n = -ntor, ntor
        IF ((rbc(n,m).ne.0) .or. (zbs(n,m).ne.0)) THEN
           WRITE(iunit,'(2(A,I4.3,A,I3.3,A,ES22.12E3))') &
              '  RBC(',n,',',m,') = ',rbc(n,m),'    ZBS(',n,',',m,') = ',zbs(n,m)
           IF (lasym) then
              WRITE(iunit,'(2(A,I4.3,A,I3.3,A,ES22.12E3))') &
                 '  RBS(',n,',',m,') = ',rbs(n,m),'    ZBC(',n,',',m,') = ',zbc(n,m)
           end if
        END IF
     END DO
  END DO

  WRITE(iunit,'(A)') '/'

  RETURN

END SUBROUTINE write_indata_namelist

END MODULE vmec_input
