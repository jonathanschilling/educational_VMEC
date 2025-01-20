!> \file
!> \brief Reading of \c wout VMEC output file.

!> \brief Reading of \c wout VMEC output file.
      MODULE read_wout_mod
!
!     USE READ_WOUT_MOD to include variables dynamically allocated
!     in this module
!     Call DEALLOCATE_READ_WOUT to free this memory when it is no longer needed
!
!     Reads in output from VMEC equilibrium code(s), contained in wout file
!
!     Contained subroutines:
!
!     read_wout_file      wrapper alias called to read/open wout file
!     read_wout_text      called by read_wout_file to read text file wout
!     read_wout_nc        called by read_wout_file to read netcdf file wout
!
!     Post-processing routines
!
!     mse_pitch           user-callable function to compute mse pitch angle
!                         for the computed equilibrium
!
      USE stel_kinds
      USE mgrid_mod

      IMPLICIT NONE
C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
! Variable names (vn_...) : put eventually into library, used by read_wout too...
      CHARACTER(LEN=*), PARAMETER ::
     1  vn_version = 'version_',
     2  vn_extension = 'input_extension', vn_mgrid = 'mgrid_file',
     3  vn_magen = 'wb', vn_therm = 'wp', vn_gam = 'gamma',
     4  vn_maxr = 'rmax_surf', vn_minr = 'rmin_surf',
     5  vn_maxz = 'zmax_surf', vn_fp = 'nfp',
     6  vn_radnod = 'ns', vn_polmod = 'mpol', vn_tormod = 'ntor',
     7  vn_maxmod = 'mnmax', vn_maxit = 'niter',
     8  vn_asym = 'lasym', vn_free = 'lfreeb',
     9  vn_error = 'ier_flag', vn_aspect = 'aspect',
     A  vn_maxmod_nyq = 'mnmax_nyq',
     B  vn_beta = 'betatotal', vn_pbeta = 'betapol',
     C  vn_tbeta = 'betator', vn_abeta = 'betaxis',
     D  vn_b0 = 'b0', vn_rbt0 = 'rbtor0', vn_rbt1 = 'rbtor',
     E  vn_sgs = 'signgs', vn_lar = 'IonLarmor', vn_modB = 'volavgB',
     F  vn_ctor = 'ctor', vn_amin = 'Aminor_p', vn_Rmaj = 'Rmajor_p',
     G  vn_vol = 'volume_p', vn_am = 'am', vn_ai = 'ai', vn_ac = 'ac',
     G  vn_ah = 'hot particle fraction', vn_atuname = 'T-perp/T-par',
     H  vn_pmass_type = 'pmass_type', vn_piota_type = 'piota_type',
     I  vn_pcurr_type = 'pcurr_type',
     J  vn_am_aux_s = 'am_aux_s', vn_am_aux_f = 'am_aux_f',
     K  vn_ai_aux_s = 'ai_aux_s', vn_ai_aux_f = 'ai_aux_f',
     L  vn_ac_aux_s = 'ac_aux_s', vn_ac_aux_f = 'ac_aux_f',
     M  vn_mse = 'imse', vn_thom = 'itse',
     N  vn_pmod = 'xm', vn_tmod = 'xn', vn_pmod_nyq = 'xm_nyq',
     O  vn_tmod_nyq = 'xn_nyq',
     P  vn_racc = 'raxis_cc', vn_zacs = 'zaxis_cs',
     Q  vn_racs = 'raxis_cs', vn_zacc = 'zaxis_cc', vn_iotaf = 'iotaf',
     Q  vn_qfact='q-factor', vn_chi='chi', vn_chipf='chipf',
     R  vn_presf = 'presf', vn_phi = 'phi', vn_phipf = 'phipf',
     S  vn_jcuru = 'jcuru', vn_jcurv = 'jcurv', vn_iotah = 'iotas',
     T  vn_mass = 'mass', vn_presh = 'pres', vn_betah = 'beta_vol',
     U  vn_buco = 'buco', vn_bvco = 'bvco', vn_vp = 'vp',
     V  vn_specw = 'specw', vn_phip = 'phips', vn_jdotb = 'jdotb',
     W  vn_bdotb = 'bdotb', vn_overr = 'over_r',
     X  vn_bgrv = 'bdotgradv', vn_merc = 'DMerc', vn_mshear = 'DShear',
     Y  vn_mwell = 'DWell', vn_mcurr = 'DCurr', vn_mgeo = 'DGeod',
     Z  vn_equif = 'equif', vn_fsq = 'fsqt', vn_wdot = 'wdot',
     1  vn_ftolv = 'ftolv', vn_fsql= 'fsql', vn_fsqr = 'fsqr',
     2  vn_fsqz = 'fsqz',
     3  vn_extcur = 'extcur', vn_curlab = 'curlabel', vn_rmnc = 'rmnc',
     4  vn_zmns = 'zmns', vn_lmns = 'lmns', vn_gmnc = 'gmnc',
     5  vn_bmnc = 'bmnc', vn_bsubumnc = 'bsubumnc',
     6  vn_bsubvmnc = 'bsubvmnc', vn_bsubsmns = 'bsubsmns',
     7  vn_bsupumnc = 'bsupumnc', vn_bsupvmnc = 'bsupvmnc',
     8  vn_rmns = 'rmns', vn_zmnc = 'zmnc',
     9  vn_lmnc = 'lmnc', vn_gmns = 'gmns', vn_bmns = 'bmns',
     A  vn_bsubumns = 'bsubumns', vn_bsubvmns = 'bsubvmns',
     B  vn_bsubsmnc = 'bsubsmnc', vn_bsupumns = 'bsupumns',
     C  vn_bsupvmns = 'bsupvmns',
     D  vn_bsubumnc_sur = 'bsubumnc_sur',
     E  vn_bsubvmnc_sur = 'bsubvmnc_sur',
     F  vn_bsupumnc_sur = 'bsupumnc_sur',
     G  vn_bsupvmnc_sur = 'bsupvmnc_sur',
     H  vn_bsubumns_sur = 'bsubumns_sur',
     I  vn_bsubvmns_sur = 'bsubvmns_sur',
     J  vn_bsupumns_sur = 'bsupumns_sur',
     K  vn_bsupvmns_sur = 'bsupvmns_sur',
     D  vn_rbc = 'rbc', vn_zbs = 'zbs', vn_rbs = 'rbs', vn_zbc = 'zbc',
     E  vn_potvac = 'potvac'

! Long names (ln_...)
      CHARACTER(LEN=*), PARAMETER ::
     1  ln_version = 'VMEC Version',
     2  ln_extension = 'Input file extension',
     3  ln_mgrid = 'MGRID file',
     4  ln_magen = 'Magnetic Energy', ln_therm = 'Thermal Energy',
     5  ln_gam = 'Gamma', ln_maxr = 'Maximum R', ln_minr = 'Minimum R',
     6  ln_maxz = 'Maximum Z', ln_fp = 'Field Periods',
     7  ln_radnod = 'Radial nodes', ln_polmod = 'Poloidal modes',
     8  ln_tormod = 'Toroidal modes', ln_maxmod = 'Fourier modes',
     8  ln_maxmod_nyq = 'Fourier modes (Nyquist)',
     9  ln_maxit = 'Max iterations',
     1  ln_asym = 'Asymmetry', ln_recon = 'Reconstruction',
     1  ln_free = 'Free boundary',
     2  ln_error = 'Error flag', ln_aspect = 'Aspect ratio',
     3  ln_beta = 'Total beta', ln_pbeta = 'Poloidal beta',
     4  ln_tbeta = 'Toroidal beta', ln_abeta = 'Beta axis',
     5  ln_b0 = 'RB-t over R axis', ln_rbt0 = 'RB-t axis',
     6  ln_rbt1 = 'RB-t edge', ln_sgs = 'Sign jacobian',
     7  ln_lar = 'Ion Larmor radius', ln_modB = 'avg mod B',
     8  ln_ctor = 'Toroidal current', ln_amin = 'minor radius',
     9  ln_Rmaj = 'major radius', ln_vol = 'Plasma volume',
     1  ln_mse = 'Number of MSE points',
     1  ln_thom = 'Number of Thompson scattering points',
     1  ln_am = 'Specification parameters for mass(s)',
     1  ln_ac = 'Specification parameters for <J>(s)',
     1  ln_ai = 'Specification parameters for iota(s)',
     1  ln_pmass_type = 'Profile type specifier for mass(s)',
     1  ln_pcurr_type = 'Profile type specifier for <J>(s)',
     1  ln_piota_type = 'Profile type specifier for iota(s)',
     1  ln_am_aux_s = 'Auxiliary-s parameters for mass(s)',
     1  ln_am_aux_f = 'Auxiliary-f parameters for mass(s)',
     1  ln_ac_aux_s = 'Auxiliary-s parameters for <J>(s)',
     1  ln_ac_aux_f = 'Auxiliary-f parameters for <J>(s)',
     1  ln_ai_aux_s = 'Auxiliary-s parameters for iota(s)',
     1  ln_ai_aux_f = 'Auxiliary-f parameters for iota(s)',
     4  ln_pmod = 'Poloidal mode numbers',
     5  ln_tmod = 'Toroidal mode numbers',
     4  ln_pmod_nyq = 'Poloidal mode numbers (Nyquist)',
     5  ln_tmod_nyq = 'Toroidal mode numbers (Nyquist)',
     5  ln_racc = 'raxis (cosnv)', ln_racs = 'raxis (sinnv)',
     6  ln_zacs = 'zaxis (sinnv)', ln_zacc = 'zaxis (cosnv)',
     7  ln_iotaf = 'iota on full mesh',
     7  ln_qfact = 'q-factor on full mesh',
     8  ln_presf = 'pressure on full mesh',
     8  ln_phi = 'Toroidal flux on full mesh',
     9  ln_phipf = 'd(phi)/ds: Toroidal flux deriv on full mesh',
     9  ln_chi = 'Poloidal flux on full mesh',
     9  ln_chipf = 'd(chi)/ds: Poroidal flux deriv on full mesh',
     9  ln_jcuru = 'j dot gradu full',
     1  ln_jcurv = 'j dot gradv full', ln_iotah = 'iota half',
     2  ln_mass = 'mass half', ln_presh = 'pressure half',
     3  ln_betah = 'beta half', ln_buco = 'bsubu half',
     4  ln_bvco = 'bsubv half', ln_vp = 'volume deriv half',
     5  ln_specw = 'Spectral width half',
     6  ln_phip = 'tor flux deriv over 2pi half',
     7  ln_jdotb = 'J dot B', ln_bdotb = 'B dot B',
     7  ln_bgrv = 'B dot grad v',
     8  ln_merc = 'Mercier criterion', ln_mshear = 'Shear Mercier',
     9  ln_mwell = 'Well Mercier', ln_mcurr = 'Current Mercier',
     1  ln_mgeo = 'Geodesic Mercier', ln_equif='Average force balance',
     1  ln_fsq = 'Residual decay',
     2  ln_wdot = 'Wdot decay', ln_extcur = 'External coil currents',
     2  ln_fsqr = 'Residual decay - radial',
     2  ln_fsqz = 'Residual decay - vertical',
     2  ln_fsql = 'Residual decay - hoop',
     2  ln_ftolv = 'Residual decay - requested',
     3  ln_curlab = 'External current names',

     3  ln_rmnc = 'cosmn component of cylindrical R, full mesh',
     4  ln_zmns = 'sinmn component of cylindrical Z, full mesh',
     4  ln_lmns = 'sinmn component of lambda, half mesh',
     5  ln_gmnc = 'cosmn component of jacobian, half mesh',
     6  ln_bmnc = 'cosmn component of mod-B, half mesh',
     6  ln_bsubumnc = 'cosmn covariant u-component of B, half mesh',
     6  ln_bsubvmnc = 'cosmn covariant v-component of B, half mesh',
     7  ln_bsubsmns = 'sinmn covariant s-component of B, full mesh',

     8  ln_bsubumnc_sur = 'cosmn bsubu of B, surface',
     9  ln_bsubvmnc_sur = 'cosmn bsubv of B, surface',
     A  ln_bsupumnc_sur = 'cosmn bsupu of B, surface',
     B  ln_bsupvmnc_sur = 'cosmn bsupv of B, surface',

     7  ln_bsupumnc = 'BSUPUmnc half',
     8  ln_bsupvmnc = 'BSUPVmnc half',

     3  ln_rmns = 'sinmn component of cylindrical R, full mesh',
     4  ln_zmnc = 'cosmn component of cylindrical Z, full mesh',
     4  ln_lmnc = 'cosmn component of lambda, half mesh',
     5  ln_gmns = 'sinmn component of jacobian, half mesh',
     6  ln_bmns = 'sinmn component of mod-B, half mesh',
     6  ln_bsubumns = 'sinmn covariant u-component of B, half mesh',
     6  ln_bsubvmns = 'sinmn covariant v-component of B, half mesh',
     7  ln_bsubsmnc = 'cosmn covariant s-component of B, full mesh',

     8  ln_bsubumns_sur = 'sinmn bsubu of B, surface',
     9  ln_bsubvmns_sur = 'sinmn bsubv of B, surface',
     A  ln_bsupumns_sur = 'sinmn bsupu of B, surface',
     B  ln_bsupvmns_sur = 'sinmn bsupv of B, surface',

     4  ln_bsupumns = 'BSUPUmns half', ln_bsupvmns = 'BSUPVmns half',
     6  ln_rbc = 'Initial boundary R cos(mu-nv) coefficients',
     7  ln_zbs = 'Initial boundary Z sin(mu-nv) coefficients',
     8  ln_rbs = 'Initial boundary R sin(mu-nv) coefficients',
     9  ln_zbc = 'Initial boundary Z cos(mu-nv) coefficients',
     1  ln_potvac = 'Vacuum Potential on Boundary'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nfp, ns, mpol, ntor, mnmax, mnmax_nyq, niter,
     1    iasym, ierr_vmec, imse, itse,
     2    isnodes, ipnodes, imatch_phiedge, isigng, mnyq, nnyq, ntmax
      REAL(rprec) :: wb, wp, gamma, pfac, rmax_surf, rmin_surf,
     1    zmax_surf, aspect, betatot, betapol, betator, betaxis, b0,
     2    tswgt, msewgt, flmwgt, bcwgt, phidiam, version_,
     3    delphid, IonLarmor, VolAvgB,
     3    fsql, fsqr, fsqz, ftolv,
     4    Aminor, Rmajor, Volume, RBtor, RBtor0, Itor,
     5    machsq !SAL
      REAL(rprec), ALLOCATABLE :: rzl_local(:,:,:,:)
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    rmnc, zmns, lmns, rmns, zmnc, lmnc, bmnc, gmnc, bsubumnc,
     2    bsubvmnc, bsubsmns, bsupumnc, bsupvmnc, currvmnc,
     3    currumnc, bbc, raxis, zaxis
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    bmns, gmns, bsubumns, bsubvmns, bsubsmnc,
     2    bsupumns, bsupvmns, currumns, currvmns
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1   iotas, iotaf, presf, phipf, mass, pres, beta_vol, xm, xn,
     1   qfact, chipf, phi, chi,
     2   xm_nyq, xn_nyq, phip, buco, bvco, vp, overr, jcuru, jcurv,
     3   specw, jdotb, bdotb, bdotgradv, fsqt, wdot, am, ac, ai,
     3   am_aux_s, am_aux_f, ac_aux_s, ac_aux_f, ai_aux_s, ai_aux_f,
     3   Dmerc, Dshear, Dwell, Dcurr, Dgeod, equif, extcur,
     4   sknots, ystark, y2stark, pknots, ythom, y2thom,
     5   anglemse, rmid, qmid, shear, presmid, alfa, curmid, rstark,
     6   qmeas, datastark, rthom, datathom, dsiobt, potvac

      LOGICAL :: lasym, lthreed, lwout_opened=.false.
      CHARACTER :: mgrid_file*200, input_extension*100
      CHARACTER :: pmass_type*20, pcurr_type*20, piota_type*20

      INTEGER, PARAMETER :: norm_term_flag=0,
     1   bad_jacobian_flag=1, more_iter_flag=2, jac75_flag=4

!     OVERLOAD SUBROUTINE READ_WOUT_FILE TO ACCEPT BOTH UNIT NO. (OPENED EXTERNALLY)
!     OR FILENAME (HANDLE OPEN/CLOSE HERE)
      INTERFACE read_wout_file
          MODULE PROCEDURE readw_and_open
      END INTERFACE

      PRIVATE :: read_wout_nc
      PRIVATE :: norm_term_flag, bad_jacobian_flag,
     1           more_iter_flag, jac75_flag

      CONTAINS

      SUBROUTINE readw_and_open(file_or_extension, ierr, iopen)
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: ierr
      INTEGER, OPTIONAL :: iopen
      CHARACTER(LEN=*), INTENT(in) :: file_or_extension
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: iunit_init = 10
      LOGICAL :: isnc
      CHARACTER(len=LEN_TRIM(file_or_extension)+10) :: filename
C-----------------------------------------------
!
!     THIS SUBROUTINE READS THE WOUT FILE CREATED BY THE VMEC CODE
!     AND STORES THE DATA IN THE READ_WOUT MODULE
!
!     FIRST, CHECK IF THIS IS A FULLY-QUALIFIED PATH NAME
!     MAKE SURE wout IS NOT EMBEDDED IN THE NAME (PERVERSE USER...)
!
      filename = 'wout'
      CALL parse_extension(filename, file_or_extension, isnc)
      CALL flush(6)
      IF (isnc) THEN
         CALL read_wout_nc(filename, ierr)
      END IF

      IF (PRESENT(iopen)) iopen = ierr
      lwout_opened = (ierr .eq. 0)
      ! WHEN READING A NETCDF FILE, A BAD RUN MAY PREVENT XN FROM BEING
      ! READ, SUBSEQUENTLY WE MUST CHECK TO SEE IF XN HAS BEEN ALLOCATED
      ! BEFORE DOING ANYTHING WITH IT OTHERWISE WE DEFAULT LTHREED TO
      ! FALSE.  - SAL 09/07/11
      IF (ALLOCATED(XN)) THEN
         lthreed = ANY(NINT(xn) .ne. 0)
      ELSE
         lthreed = .FALSE.
      END IF

      END SUBROUTINE readw_and_open

      SUBROUTINE read_wout_nc(filename, ierr)
      USE ezcdf
      USE stel_constants, ONLY: mu0
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: ierr
      CHARACTER(LEN=*), INTENT(in) :: filename
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nwout, ierror
      INTEGER, DIMENSION(3)   :: dimlens
      REAL(rprec) :: ohs
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: raxis_cc, raxis_cs,
     1                                          zaxis_cs, zaxis_cc
C-----------------------------------------------
! Open cdf File
      CALL cdf_open(nwout,filename,'r', ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Error opening wout .nc file'
         RETURN
      END IF

! Be sure all arrays are deallocated
      CALL read_wout_deallocate

      ! reset iasym, needed when using the same instance multiple times -jons
      iasym = 0

! Read in scalar variables
      CALL cdf_read(nwout, vn_error, ierr_vmec)

      IF (ierr_vmec.ne.norm_term_flag .and. ierr_vmec.ne.more_iter_flag)
     1   GOTO 1000

      CALL cdf_read(nwout, vn_version, version_)
      CALL cdf_read(nwout, vn_extension, input_extension)
      CALL cdf_read(nwout, vn_mgrid, mgrid_file)
      CALL cdf_read(nwout, vn_magen, wb)
      CALL cdf_read(nwout, vn_therm, wp)
      CALL cdf_read(nwout, vn_gam, gamma)
      CALL cdf_read(nwout, vn_maxr, rmax_surf)
      CALL cdf_read(nwout, vn_minr, rmin_surf)
      CALL cdf_read(nwout, vn_maxz, zmax_surf)
      CALL cdf_read(nwout, vn_fp, nfp)
      CALL cdf_read(nwout, vn_radnod, ns)
      CALL cdf_read(nwout, vn_polmod, mpol)
      CALL cdf_read(nwout, vn_tormod, ntor)
      CALL cdf_read(nwout, vn_maxmod, mnmax)
      mnmax_nyq = -1
      CALL cdf_read(nwout, vn_maxmod_nyq, mnmax_nyq)
      CALL cdf_read(nwout, vn_maxit, niter)
      CALL cdf_read(nwout, vn_asym, lasym)
      IF (lasym) iasym = 1
      CALL cdf_read(nwout, vn_free, lfreeb)
      CALL cdf_read(nwout, vn_aspect, aspect)
      CALL cdf_read(nwout, vn_beta, betatot)
      CALL cdf_read(nwout, vn_pbeta, betapol)
      CALL cdf_read(nwout, vn_tbeta, betator)
      CALL cdf_read(nwout, vn_abeta, betaxis)
      CALL cdf_read(nwout, vn_b0, b0)
      CALL cdf_read(nwout, vn_rbt0, rbtor0)
      CALL cdf_read(nwout, vn_rbt1, rbtor)
      CALL cdf_read(nwout, vn_sgs, isigng)
      CALL cdf_read(nwout, vn_lar, IonLarmor)
      CALL cdf_read(nwout, vn_modB, volAvgB)
      CALL cdf_read(nwout, vn_ctor, Itor)
      CALL cdf_read(nwout, vn_amin, Aminor)
      CALL cdf_read(nwout, vn_rmaj, Rmajor)
      CALL cdf_read(nwout, vn_vol, volume)
      CALL cdf_read(nwout, vn_ftolv, ftolv)
      CALL cdf_read(nwout, vn_fsqr, fsqr)
      CALL cdf_read(nwout, vn_fsqz, fsqz)
      CALL cdf_read(nwout, vn_fsql, fsql)
      CALL cdf_read(nwout, vn_pcurr_type, pcurr_type)
      CALL cdf_read(nwout, vn_piota_type, piota_type)
      CALL cdf_read(nwout, vn_pmass_type, pmass_type)
      imse = -1
      CALL cdf_read(nwout, vn_nextcur, nextcur)

      mgrid_mode = 'N'
      CALL cdf_inquire(nwout, vn_mgmode, dimlens, ier=ierror)
      IF (ierror.eq.0) CALL cdf_read(nwout, vn_mgmode, mgrid_mode)

! Inquire existence, dimensions of arrays for allocation
! 1D Arrays
      CALL cdf_inquire(nwout, vn_pmod, dimlens)
      ALLOCATE (xm(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_tmod, dimlens)
      ALLOCATE (xn(dimlens(1)), stat = ierror)
      IF (mnmax_nyq .gt. 0) THEN
         CALL cdf_inquire(nwout, vn_pmod_nyq, dimlens)
         ALLOCATE (xm_nyq(dimlens(1)), stat = ierror)
         CALL cdf_inquire(nwout, vn_tmod_nyq, dimlens)
         ALLOCATE (xn_nyq(dimlens(1)), stat = ierror)
      END IF

      CALL cdf_inquire(nwout, vn_racc, dimlens)
      ALLOCATE (raxis_cc(0:dimlens(1)-1), stat = ierror)
      CALL cdf_inquire(nwout, vn_zacs, dimlens)
      ALLOCATE (zaxis_cs(0:dimlens(1)-1), stat = ierror)
      IF (lasym) THEN
         CALL cdf_inquire(nwout, vn_racs, dimlens)
         ALLOCATE (raxis_cs(0:dimlens(1)-1), stat = ierror)
         CALL cdf_inquire(nwout, vn_zacc, dimlens)
         ALLOCATE (zaxis_cc(0:dimlens(1)-1), stat = ierror)
      END IF

!  Profile coefficients, dimensioned from 0
      CALL cdf_inquire(nwout, vn_am, dimlens)
      ALLOCATE (am(0:dimlens(1)-1), stat = ierror)
      CALL cdf_inquire(nwout, vn_ac, dimlens)
      ALLOCATE (ac(0:dimlens(1)-1), stat = ierror)
      CALL cdf_inquire(nwout, vn_ai, dimlens)
      ALLOCATE (ai(0:dimlens(1)-1), stat = ierror)

      CALL cdf_inquire(nwout, vn_ac_aux_s, dimlens)
      ALLOCATE (ac_aux_s(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_ac_aux_f, dimlens)
      ALLOCATE (ac_aux_f(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_ai_aux_s, dimlens)
      ALLOCATE (ai_aux_s(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_ai_aux_f, dimlens)
      ALLOCATE (ai_aux_f(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_am_aux_s, dimlens)
      ALLOCATE (am_aux_s(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_am_aux_f, dimlens)
      ALLOCATE (am_aux_f(dimlens(1)), stat = ierror)

      CALL cdf_inquire(nwout, vn_iotaf, dimlens)
      ALLOCATE (iotaf(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_qfact, dimlens)
      ALLOCATE (qfact(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_presf, dimlens)
      ALLOCATE (presf(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_phi, dimlens)
      ALLOCATE (phi(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_chi, dimlens)
      ALLOCATE (chi(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_phipf, dimlens)
      ALLOCATE (phipf(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_chipf, dimlens)
      ALLOCATE (chipf(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_jcuru, dimlens)
      ALLOCATE (jcuru(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_jcurv, dimlens)
      ALLOCATE (jcurv(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_iotah, dimlens)
      ALLOCATE (iotas(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mass, dimlens)
      ALLOCATE (mass(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_presh, dimlens)
      ALLOCATE (pres(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_betah, dimlens)
      ALLOCATE (beta_vol(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_buco, dimlens)
      ALLOCATE (buco(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bvco, dimlens)
      ALLOCATE (bvco(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_vp, dimlens)
      ALLOCATE (vp(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_specw, dimlens)
      ALLOCATE (specw(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_phip, dimlens)
      ALLOCATE (phip(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_overr, dimlens)
      ALLOCATE (overr(dimlens(1)), stat = ierror)

      CALL cdf_inquire(nwout, vn_jdotb, dimlens)
      ALLOCATE (jdotb(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bdotb, dimlens)
      ALLOCATE (bdotb(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bgrv, dimlens)
      ALLOCATE (bdotgradv(dimlens(1)), stat = ierror)

      CALL cdf_inquire(nwout, vn_merc, dimlens)
      ALLOCATE (Dmerc(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mshear, dimlens)
      ALLOCATE (Dshear(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mwell, dimlens)
      ALLOCATE (Dwell(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mcurr, dimlens)
      ALLOCATE (Dcurr(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mgeo, dimlens)
      ALLOCATE (Dgeod(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_equif, dimlens)
      ALLOCATE (equif(dimlens(1)), stat = ierror)

      CALL cdf_inquire(nwout, vn_fsq, dimlens)
      ALLOCATE (fsqt(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_wdot, dimlens)
      ALLOCATE (wdot(dimlens(1)), stat = ierror)

      IF (nextcur .gt. 0) THEN
         CALL cdf_inquire(nwout, vn_extcur, dimlens)
         ALLOCATE (extcur(dimlens(1)), stat = ierror)
!NOTE: curlabel is an array of CHARACTER(30) strings - defined in mgrid_mod
!      so dimlens(1) == 30 (check this) and dimlens(2) is the number of strings in the array
         CALL cdf_inquire(nwout, vn_curlab, dimlens)
         ALLOCATE (curlabel(dimlens(2)), stat = ierror)
         ! SAL
         CALL cdf_inquire(nwout, vn_potvac, dimlens, ier = ierror)
         IF (ierror == 0) ALLOCATE (potvac(1:dimlens(1)), stat = ierror)
      ENDIF

! 2D Arrays
      CALL cdf_inquire(nwout, vn_rmnc, dimlens)
      ALLOCATE (rmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_zmns, dimlens)
      ALLOCATE (zmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_lmns, dimlens)
      ALLOCATE (lmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_gmnc, dimlens)
      ALLOCATE (gmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bmnc, dimlens)
      ALLOCATE (bmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubumnc, dimlens)
      ALLOCATE (bsubumnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubvmnc, dimlens)
      ALLOCATE (bsubvmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubsmns, dimlens)
      ALLOCATE (bsubsmns(dimlens(1),dimlens(2)), stat = ierror)

!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM
      CALL cdf_inquire(nwout, vn_bsupumnc, dimlens)
      ALLOCATE (bsupumnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsupvmnc, dimlens)
      ALLOCATE (bsupvmnc(dimlens(1),dimlens(2)), stat = ierror)

      IF (.NOT. lasym) GO TO 800

      CALL cdf_inquire(nwout, vn_rmns, dimlens)
      ALLOCATE (rmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_zmnc, dimlens)
      ALLOCATE (zmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_lmnc, dimlens)
      ALLOCATE (lmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_gmns, dimlens)
      ALLOCATE (gmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bmns, dimlens)
      ALLOCATE (bmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubumns, dimlens)
      ALLOCATE (bsubumns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubvmns, dimlens)
      ALLOCATE (bsubvmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubsmnc, dimlens)
      ALLOCATE (bsubsmnc(dimlens(1),dimlens(2)), stat = ierror)

!     ELIMINATE THESE EVENTUALLY: DO NOT NEED THEM
      CALL cdf_inquire(nwout, vn_bsupumns, dimlens)
      ALLOCATE (bsupumns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsupvmns, dimlens)
      ALLOCATE (bsupvmns(dimlens(1),dimlens(2)), stat = ierror)

 800  CONTINUE

! Read Arrays
      CALL cdf_read(nwout, vn_pmod, xm)
      CALL cdf_read(nwout, vn_tmod, xn)
      IF (mnmax_nyq .le. 0) THEN
         mnmax_nyq = mnmax
         ALLOCATE (xm_nyq(mnmax_nyq), xn_nyq(mnmax_nyq), stat=ierror)
         xm_nyq = xm;  xn_nyq = xn
      ELSE
         CALL cdf_read(nwout, vn_pmod_nyq, xm_nyq)
         CALL cdf_read(nwout, vn_tmod_nyq, xn_nyq)
      END IF

      mnyq = INT(MAXVAL(xm_nyq));  nnyq = INT(MAXVAL(ABS(xn_nyq)))/nfp

      CALL cdf_read(nwout, vn_racc, raxis_cc)
      CALL cdf_read(nwout, vn_zacs, zaxis_cs)

      IF (SIZE(raxis_cc) .ne. ntor+1)
     1   STOP 'WRONG SIZE(raxis_cc) in READ_WOUT_NC'
      ALLOCATE (raxis(0:ntor,2), zaxis(0:ntor,2), stat=ierror)
      raxis(:,1) = raxis_cc(0:ntor);   zaxis(:,1) = zaxis_cs(0:ntor)
      raxis(:,2) = 0;                  zaxis(:,2) = 0
      DEALLOCATE (raxis_cc, zaxis_cs, stat=ierror)

      CALL cdf_read(nwout, vn_rmnc, rmnc)
      CALL cdf_read(nwout, vn_zmns, zmns)
      CALL cdf_read(nwout, vn_lmns, lmns)
      CALL cdf_read(nwout, vn_gmnc, gmnc)              !Half mesh
      CALL cdf_read(nwout, vn_bmnc, bmnc)              !Half mesh
      CALL cdf_read(nwout, vn_bsubumnc, bsubumnc)      !Half mesh
      CALL cdf_read(nwout, vn_bsubvmnc, bsubvmnc)      !Half mesh
      CALL cdf_read(nwout, vn_bsubsmns, bsubsmns)      !Full mesh
!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM (can express in terms of lambdas)
      CALL cdf_read(nwout, vn_bsupumnc, bsupumnc)
      CALL cdf_read(nwout, vn_bsupvmnc, bsupvmnc)
      IF (lasym) THEN
         CALL cdf_read(nwout, vn_racs, raxis_cs)
         CALL cdf_read(nwout, vn_zacc, zaxis_cc)
         raxis(:,2) = raxis_cs;   zaxis(:,2) = zaxis_cc
         DEALLOCATE (raxis_cs, zaxis_cc, stat=ierror)
         CALL cdf_read(nwout, vn_rmns, rmns)
         CALL cdf_read(nwout, vn_zmnc, zmnc)
         CALL cdf_read(nwout, vn_lmnc, lmnc)
         CALL cdf_read(nwout, vn_gmns, gmns)
         CALL cdf_read(nwout, vn_bmns, bmns)
         CALL cdf_read(nwout, vn_bsubumns, bsubumns)
         CALL cdf_read(nwout, vn_bsubvmns, bsubvmns)
         CALL cdf_read(nwout, vn_bsubsmnc, bsubsmnc)
!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM
         CALL cdf_read(nwout, vn_bsupumns, bsupumns)
         CALL cdf_read(nwout, vn_bsupvmns, bsupvmns)
      END IF

      CALL cdf_read(nwout, vn_am, am)
      CALL cdf_read(nwout, vn_ac, ac)
      CALL cdf_read(nwout, vn_ai, ai)

      CALL cdf_read(nwout, vn_am_aux_s, am_aux_s)
      CALL cdf_read(nwout, vn_am_aux_f, am_aux_f)
      CALL cdf_read(nwout, vn_ac_aux_s, ac_aux_s)
      CALL cdf_read(nwout, vn_ac_aux_f, ac_aux_f)
      CALL cdf_read(nwout, vn_ai_aux_s, ai_aux_s)
      CALL cdf_read(nwout, vn_ai_aux_f, ai_aux_f)

      CALL cdf_read(nwout, vn_iotaf, iotaf)
      CALL cdf_read(nwout, vn_qfact, qfact)
      CALL cdf_read(nwout, vn_presf, presf)
      CALL cdf_read(nwout, vn_phi, phi)
      CALL cdf_read(nwout, vn_phipf, phipf)
      CALL cdf_read(nwout, vn_chi, chi)
      CALL cdf_read(nwout, vn_chipf, chipf)
      CALL cdf_read(nwout, vn_jcuru, jcuru)
      CALL cdf_read(nwout, vn_jcurv, jcurv)


!     HALF-MESH quantities
!     NOTE: jdotb is in units_of_A (1/mu0 incorporated in jxbforce...)
!     prior to version 6.00, this was output in internal VMEC units...
      CALL cdf_read(nwout, vn_iotah, iotas)
      CALL cdf_read(nwout, vn_mass, mass)
      CALL cdf_read(nwout, vn_presh, pres)
      CALL cdf_read(nwout, vn_betah, beta_vol)
      CALL cdf_read(nwout, vn_buco, buco)
      CALL cdf_read(nwout, vn_bvco, bvco)
      CALL cdf_read(nwout, vn_vp, vp)
      CALL cdf_read(nwout, vn_specw, specw)
      CALL cdf_read(nwout, vn_phip, phip)
      CALL cdf_read(nwout, vn_jdotb, jdotb)
      CALL cdf_read(nwout, vn_bdotb, bdotb)
      CALL cdf_read(nwout, vn_bgrv, bdotgradv)

!     MERCIER_CRITERION
      CALL cdf_read(nwout, vn_merc, Dmerc)
      CALL cdf_read(nwout, vn_mshear, Dshear)
      CALL cdf_read(nwout, vn_mwell, Dwell)
      CALL cdf_read(nwout, vn_mcurr, Dcurr)
      CALL cdf_read(nwout, vn_mgeo, Dgeod)
      CALL cdf_read(nwout, vn_equif, equif)

      CALL cdf_read(nwout, vn_fsq, fsqt)
      CALL cdf_read(nwout, vn_wdot, wdot)

      IF (nextcur .gt. 0) THEN
         CALL cdf_read(nwout, vn_extcur, extcur)
         CALL cdf_read(nwout, vn_curlab, curlabel)
      ENDIF

      !SAL Addition
      IF (ALLOCATED(potvac)) CALL cdf_read(nwout,vn_potvac, potvac)

 1000 CONTINUE

      CALL cdf_close(nwout, ierr)

      IF (.not.ALLOCATED(bsubumnc)) RETURN                              !Moved this here because ns may not be set. SAL -09/07/11
!
!     COMPUTE CONTRAVARIANT CURRENT COMPONENTS IN AMPS
!     ON THE FULL RADIAL MESH, WHERE JACOBIAN = SQRT(G)
!
!     CURRU = SQRT(G) * J dot grad(u)
!     CURRV = SQRT(G) * J dot grad(v)
!
      ohs = (ns-1)


      IF (ierror .eq. 0) CALL Compute_Currents(ierror)

      IF (ierr. ne. 0)   PRINT *,"in read_wout_nc ierr=",ierr
      IF (ierror. ne. 0) PRINT *,"in read_wout_nc ierror=",ierror

      END SUBROUTINE read_wout_nc

      SUBROUTINE Compute_Currents(ierror)
      USE stel_constants, ONLY: mu0
      IMPLICIT NONE
      INTEGER, INTENT(out) :: ierror
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: js
      REAL(rprec) :: ohs, hs, shalf(ns), sfull(ns)
      REAL(rprec), DIMENSION(mnmax_nyq) :: bu1, bu0, bv1, bv0, t1, t2,
     &                                      t3
!-----------------------------------------------
!
!     Computes current harmonics for currXmn == sqrt(g)*JsupX, X = u,v
!     [Corrected above "JsubX" to "JsupX", JDH 2010-08-16]

!     NOTE: bsub(s,u,v)mn are on HALF radial grid
!          (in earlier versions, bsubsmn was on FULL radial grid)

!
      ohs = (ns-1)
      hs  = 1._dp/ohs

      DO js = 2, ns
         shalf(js) = SQRT(hs*(js-1.5_dp))
         sfull(js) = SQRT(hs*(js-1))
      END DO

      ALLOCATE (currumnc(mnmax_nyq,ns), currvmnc(mnmax_nyq,ns),         &
     &          stat=ierror)
      IF (ierror .ne. 0) RETURN

      DO js = 2, ns-1
         WHERE (MOD(INT(xm_nyq),2) .EQ. 1)
            t1 = 0.5_dp*(shalf(js+1)*bsubsmns(:,js+1) +                  &
     &                   shalf(js)  *bsubsmns(:,js)) /sfull(js)
            bu0 = bsubumnc(:,js  )/shalf(js)
            bu1 = bsubumnc(:,js+1)/shalf(js+1)
            t2 = ohs*(bu1-bu0)*sfull(js)+0.25_dp*(bu0+bu1)/sfull(js)
            bv0 = bsubvmnc(:,js  )/shalf(js)
            bv1 = bsubvmnc(:,js+1)/shalf(js+1)
            t3 = ohs*(bv1-bv0)*sfull(js)+0.25_dp*(bv0+bv1)/sfull(js)
         ELSEWHERE
            t1 = 0.5_dp*(bsubsmns(:,js+1)+bsubsmns(:,js))
            t2 = ohs*(bsubumnc(:,js+1)-bsubumnc(:,js))
            t3 = ohs*(bsubvmnc(:,js+1)-bsubvmnc(:,js))
         ENDWHERE
         currumnc(:,js) = -xn_nyq(:)*t1 - t3
         currvmnc(:,js) = -xm_nyq(:)*t1 + t2
      END DO

      WHERE (xm_nyq .LE. 1)
         currvmnc(:,1) =  2*currvmnc(:,2) - currvmnc(:,3)
         currumnc(:,1) =  2*currumnc(:,2) - currumnc(:,3)
      ELSEWHERE
         currvmnc(:,1) = 0
         currumnc(:,1) = 0
      ENDWHERE

      currumnc(:,ns) = 2*currumnc(:,ns-1) - currumnc(:,ns-2)
      currvmnc(:,ns) = 2*currvmnc(:,ns-1) - currvmnc(:,ns-2)
      currumnc = currumnc/mu0;   currvmnc = currvmnc/mu0

      IF (.NOT.lasym) RETURN

      ALLOCATE (currumns(mnmax_nyq,ns), currvmns(mnmax_nyq,ns),         &
     &           stat=ierror)

      DO js = 2, ns-1
         WHERE (MOD(INT(xm_nyq),2) .EQ. 1)
            t1 = 0.5_dp*(shalf(js+1)*bsubsmnc(:,js+1)                   &
     &          +         shalf(js)  *bsubsmnc(:,js)) / sfull(js)
            bu0 = bsubumns(:,js  )/shalf(js+1)
            bu1 = bsubumns(:,js+1)/shalf(js+1)
            t2 = ohs*(bu1-bu0)*sfull(js) + 0.25_dp*(bu0+bu1)/sfull(js)
            bv0 = bsubvmns(:,js  )/shalf(js)
            bv1 = bsubvmns(:,js+1)/shalf(js+1)
            t3 = ohs*(bv1-bv0)*sfull(js)+0.25_dp*(bv0+bv1)/sfull(js)
         ELSEWHERE
            t1 = 0.5_dp*(bsubsmnc(:,js+1) + bsubsmnc(:,js))
            t2 = ohs*(bsubumns(:,js+1)-bsubumns(:,js))
            t3 = ohs*(bsubvmns(:,js+1)-bsubvmns(:,js))
         END WHERE
         currumns(:,js) =  xn_nyq(:)*t1 - t3
         currvmns(:,js) =  xm_nyq(:)*t1 + t2
      END DO

      WHERE (xm_nyq .LE. 1)
         currvmns(:,1) =  2*currvmns(:,2) - currvmns(:,3)
         currumns(:,1) =  2*currumns(:,2) - currumns(:,3)
      ELSEWHERE
         currvmns(:,1) = 0
         currumns(:,1) = 0
      END WHERE
      currumns(:,ns) = 2*currumns(:,ns-1) - currumns(:,ns-2)
      currvmns(:,ns) = 2*currvmns(:,ns-1) - currvmns(:,ns-2)
      currumns = currumns/mu0;   currvmns = currvmns/mu0

      END SUBROUTINE Compute_Currents

      SUBROUTINE read_wout_deallocate
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat(10)
!-----------------------------------------------
      istat = 0
      lwout_opened = .false.

      IF (ALLOCATED(extcur)) DEALLOCATE (extcur,
     1         stat = istat(1))
      IF (ALLOCATED(curlabel)) DEALLOCATE (curlabel,
     1         stat = istat(1))
      IF (ALLOCATED(overr)) DEALLOCATE (overr, stat = istat(2))

      IF (ALLOCATED(xm)) DEALLOCATE (xm, xn, xm_nyq, xn_nyq,
     1  rmnc, zmns, lmns, bmnc, gmnc, bsubumnc, iotaf, presf, phipf,
     2  bsubvmnc, bsubsmns, bsupumnc, bsupvmnc, currvmnc, iotas, mass,
     3  pres, beta_vol, phip, buco, bvco, phi, vp, jcuru, am, ac, ai,
     4  jcurv, specw, Dmerc, Dshear, Dwell, Dcurr, Dgeod, equif, jdotb,
     5  bdotb, bdotgradv, raxis, zaxis, fsqt, wdot, stat = istat(3))

      IF (ALLOCATED(chipf)) DEALLOCATE (chipf, chi)

      IF (ALLOCATED(am_aux_s)) DEALLOCATE (am_aux_s, am_aux_f, ac_aux_s,
     1  ac_aux_f, ai_aux_s, ai_aux_f, stat=istat(6))

      IF (ALLOCATED(rmns)) DEALLOCATE (rmns, zmnc, lmnc,
     1    bmns, gmns, bsubumns, bsubvmns, bsubsmnc,
     2    bsupumns, bsupvmns, stat=istat(5))

      IF (ALLOCATED(currumnc)) DEALLOCATE (currumnc)
      IF (ALLOCATED(currumns)) DEALLOCATE (currumns, currvmns)
      IF (ALLOCATED(rzl_local)) DEALLOCATE (rzl_local)

      ! SAL Addition
      IF (ALLOCATED(potvac)) DEALLOCATE(potvac)

      IF (ANY(istat .ne. 0)) THEN
        PRINT *,istat
        STOP 'Deallocation error in read_wout_deallocate'
      END IF

      END SUBROUTINE read_wout_deallocate

      SUBROUTINE tosuvspace (s_in, u_in, v_in, gsqrt,
     1     bsupu, bsupv, jsupu, jsupv, lam)
      USE stel_constants, ONLY: zero, one
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(in) :: s_in, u_in, v_in
      REAL(rprec), INTENT(out), OPTIONAL :: gsqrt, bsupu, bsupv,
     1    jsupu, jsupv, lam
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
      INTEGER :: m, n, n1, mn, ipresent, jslo, jshi
      REAL(rprec) :: hs1, wlo, whi, wlo_odd, whi_odd
      REAL(rprec), DIMENSION(mnmax_nyq) :: gmnc1, gmns1, bsupumnc1,
     1   bsupumns1, bsupvmnc1, bsupvmns1, jsupumnc1, jsupumns1,
     2   jsupvmnc1, jsupvmns1, wmins, wplus, lammns1, lammnc1
      REAL(rprec) :: cosu, sinu, cosv, sinv, tcosmn, tsinmn, sgn
      REAL(rprec) :: cosmu(0:mnyq), sinmu(0:mnyq),
     1               cosnv(0:nnyq), sinnv(0:nnyq)
      LOGICAL :: lgsqrt, lbsupu, lbsupv, ljsupu, ljsupv, llam
C-----------------------------------------------
!
!     COMPUTE VARIOUS HALF/FULL-RADIAL GRID QUANTITIES AT THE INPUT POINT
!     (S, U, V) , WHERE
!        S = normalized toroidal flux (0 - 1),
!        U = poloidal angle
!        V = N*phi = toroidal angle * no. field periods
!
!     HALF-RADIAL GRID QUANTITIES
!     gsqrt, bsupu, bsupv
!
!     FULL-RADIAL GRID QUANTITIES
!     dbsubuds, dbsubvds, dbsubsdu, dbsubsdv
!
C-----------------------------------------------
      IF (s_in.lt.zero .or. s_in.gt.one) THEN
         WRITE(6, *)
     1   ' In tosuvspace, s(flux) must be between 0 and 1'
         RETURN
      END IF

      IF (.not.lwout_opened) THEN
         WRITE(6, *)
     1   ' tosuvspace can only be called AFTER opening wout file!'
         RETURN
      END IF

!
!     SETUP TRIG ARRAYS
!
      cosu = COS(u_in);   sinu = SIN(u_in)
      cosv = COS(v_in);   sinv = SIN(v_in)

      cosmu(0) = 1;    sinmu(0) = 0
      cosnv(0) = 1;    sinnv(0) = 0
      DO m = 1, mnyq
         cosmu(m) = cosmu(m-1)*cosu - sinmu(m-1)*sinu
         sinmu(m) = sinmu(m-1)*cosu + cosmu(m-1)*sinu
      END DO

      DO n = 1, nnyq
         cosnv(n) = cosnv(n-1)*cosv - sinnv(n-1)*sinv
         sinnv(n) = sinnv(n-1)*cosv + cosnv(n-1)*sinv
      END DO


!
!     FIND INTERPOLATED s VALUE AND COMPUTE INTERPOLATION WEIGHTS wlo, whi
!     RECALL THAT THESE QUANTITIES ARE ON THE HALF-RADIAL GRID...
!     s-half(j) = (j-1.5)*hs, for j = 2,...ns
!
      hs1 = one/(ns-1)
      jslo = INT(c1p5 + s_in/hs1)
      jshi = jslo+1
      wlo = (hs1*(jshi-c1p5) - s_in)/hs1
      whi = 1 - wlo
      IF (jslo .eq. ns) THEN
!        USE Xhalf(ns+1) = 2*Xhalf(ns) - Xhalf(ns-1) FOR "GHOST" POINT VALUE 1/2hs OUTSIDE EDGE
!        THEN, X = wlo*Xhalf(ns) + whi*Xhalf(ns+1) == Xhalf(ns) + whi*(Xhalf(ns) - Xhalf(ns-1))
         jshi = jslo-1
         wlo = 1+whi; whi = -whi
      ELSE IF (jslo .eq. 1) THEN
         jslo = 2
      END IF

!
!     FOR ODD-m MODES X ~ SQRT(s), SO INTERPOLATE Xmn/SQRT(s)
!
      whi_odd = whi*SQRT(s_in/(hs1*(jshi-c1p5)))
      IF (jslo .ne. 1) THEN
         wlo_odd = wlo*SQRT(s_in/(hs1*(jslo-c1p5)))
      ELSE
         wlo_odd = 0
         whi_odd = SQRT(s_in/(hs1*(jshi-c1p5)))
      END IF

      WHERE (MOD(NINT(xm_nyq(:)),2) .eq. 0)
         wmins = wlo
         wplus = whi
      ELSEWHERE
         wmins = wlo_odd
         wplus = whi_odd
      END WHERE

      ipresent = 0
      lgsqrt = PRESENT(gsqrt)
      IF (lgsqrt) THEN
         gsqrt = 0 ;  ipresent = ipresent+1
         gmnc1 = wmins*gmnc(:,jslo) + wplus*gmnc(:,jshi)
         IF (lasym)
     1   gmns1 = wmins*gmns(:,jslo) + wplus*gmns(:,jshi)
      END IF
      lbsupu = PRESENT(bsupu)
      IF (lbsupu) THEN
         bsupu = 0 ;  ipresent = ipresent+1
         bsupumnc1 = wmins*bsupumnc(:,jslo) + wplus*bsupumnc(:,jshi)
         IF (lasym)
     1   bsupumns1 = wmins*bsupumns(:,jslo) + wplus*bsupumns(:,jshi)
      END IF
      lbsupv = PRESENT(bsupv)
      IF (lbsupv) THEN
         bsupv = 0 ;  ipresent = ipresent+1
         bsupvmnc1 = wmins*bsupvmnc(:,jslo) + wplus*bsupvmnc(:,jshi)
         IF (lasym)
     1   bsupvmns1 = wmins*bsupvmns(:,jslo) + wplus*bsupvmns(:,jshi)
      END IF
      llam = PRESENT(lam)
      IF (llam) THEN
         lam = 0 ;  ipresent = ipresent+1
         lammns1 = wmins*lmns(:,jslo) + wplus*lmns(:,jshi)
         IF (lasym)
     1   lammnc1 = wmins*lmnc(:,jslo) + wplus*lmnc(:,jshi)
      END IF

      IF (ipresent .eq. 0) GOTO 1000

!
!     COMPUTE GSQRT, ... IN REAL SPACE
!     tcosmn = cos(mu - nv);  tsinmn = sin(mu - nv)
!
      DO mn = 1, mnmax_nyq
         m = NINT(xm_nyq(mn));  n = NINT(xn_nyq(mn))/nfp
         n1 = ABS(n);   sgn = SIGN(1,n)
         tcosmn = cosmu(m)*cosnv(n1) + sgn*sinmu(m)*sinnv(n1)
         tsinmn = sinmu(m)*cosnv(n1) - sgn*cosmu(m)*sinnv(n1)
         IF (lgsqrt) gsqrt = gsqrt + gmnc1(mn)*tcosmn
         IF (lbsupu) bsupu = bsupu + bsupumnc1(mn)*tcosmn
         IF (lbsupv) bsupv = bsupv + bsupvmnc1(mn)*tcosmn
         IF (llam)   lam = lam + lammns1(mn)*tsinmn
      END DO

      IF (.not.lasym) GOTO 1000

      DO mn = 1, mnmax_nyq
         m = NINT(xm_nyq(mn));  n = NINT(xn_nyq(mn))/nfp
         n1 = ABS(n);   sgn = SIGN(1,n)
         tcosmn = cosmu(m)*cosnv(n1) + sgn*sinmu(m)*sinnv(n1)
         tsinmn = sinmu(m)*cosnv(n1) - sgn*cosmu(m)*sinnv(n1)
         IF (lgsqrt) gsqrt = gsqrt + gmns1(mn)*tsinmn
         IF (lbsupu) bsupu = bsupu + bsupumns1(mn)*tsinmn
         IF (lbsupv) bsupv = bsupv + bsupvmns1(mn)*tsinmn
         IF (llam)   lam = lam + lammnc1(mn)*tcosmn
      END DO

 1000 CONTINUE

!     FULL-MESH QUANTITIES
!
!     FIND INTERPOLATED s VALUE AND COMPUTE INTERPOLATION WEIGHTS wlo, whi
!     RECALL THAT THESE QUANTITIES ARE ON THE FULL-RADIAL GRID...
!     s-full(j) = (j-1)*hs, for j = 1,...ns
!
      hs1 = one/(ns-1)
      jslo = 1+INT(s_in/hs1)
      jshi = jslo+1
      IF (jslo .eq. ns) jshi = ns
      wlo = (hs1*(jshi-1) - s_in)/hs1
      whi = 1 - wlo
!
!     FOR ODD-m MODES X ~ SQRT(s), SO INTERPOLATE Xmn/SQRT(s)
!
      whi_odd = whi*SQRT(s_in/(hs1*(jshi-1)))
      IF (jslo .ne. 1) THEN
         wlo_odd = wlo*SQRT(s_in/(hs1*(jslo-1)))
      ELSE
         wlo_odd = 0
         whi_odd = SQRT(s_in/(hs1*(jshi-1)))
      END IF

      WHERE (MOD(NINT(xm_nyq(:)),2) .eq. 0)
         wmins = wlo
         wplus = whi
      ELSEWHERE
         wmins = wlo_odd
         wplus = whi_odd
      END WHERE

      ipresent = 0
      ljsupu = PRESENT(jsupu)
      IF (ljsupu) THEN
         IF (.not.lgsqrt) STOP 'MUST compute gsqrt for jsupu'
         jsupu = 0 ;  ipresent = ipresent+1
         jsupumnc1 = wmins*currumnc(:,jslo) + wplus*currumnc(:,jshi)
         IF (lasym)
     1   jsupumns1 = wmins*currumns(:,jslo) + wplus*currumns(:,jshi)
      END IF

      ljsupv = PRESENT(jsupv)
      IF (ljsupv) THEN
         IF (.not.lgsqrt) STOP 'MUST compute gsqrt for jsupv'
         jsupv = 0 ;  ipresent = ipresent+1
         jsupvmnc1 = wmins*currvmnc(:,jslo) + wplus*currvmnc(:,jshi)
         IF (lasym)
     1   jsupvmns1 = wmins*currvmns(:,jslo) + wplus*currvmns(:,jshi)
      END IF

      IF (ipresent .eq. 0) RETURN

      DO mn = 1, mnmax_nyq
         m = NINT(xm_nyq(mn));  n = NINT(xn_nyq(mn))/nfp
         n1 = ABS(n);   sgn = SIGN(1,n)
         tcosmn = cosmu(m)*cosnv(n1) + sgn*sinmu(m)*sinnv(n1)
         IF (ljsupu) jsupu = jsupu + jsupumnc1(mn)*tcosmn
         IF (ljsupv) jsupv = jsupv + jsupvmnc1(mn)*tcosmn
      END DO

      IF (.not.lasym) GOTO 2000

      DO mn = 1, mnmax_nyq
         m = NINT(xm_nyq(mn));  n = NINT(xn_nyq(mn))/nfp
         n1 = ABS(n);   sgn = SIGN(1,n)
         tsinmn = sinmu(m)*cosnv(n1) - sgn*cosmu(m)*sinnv(n1)
         IF (ljsupu) jsupu = jsupu + jsupumns1(mn)*tsinmn
         IF (ljsupv) jsupv = jsupv + jsupvmns1(mn)*tsinmn
      END DO

 2000 CONTINUE

      IF (ljsupu) jsupu = jsupu/gsqrt
      IF (ljsupv) jsupv = jsupv/gsqrt

      END SUBROUTINE tosuvspace

      SUBROUTINE LoadRZL
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER     :: rcc, rss, zsc, zcs, rsc, rcs, zcc, zss
      INTEGER     :: mpol1, mn, m, n, n1
      REAL(rprec) :: sgn
C-----------------------------------------------
!
!     Arrays must be stacked (and ns,ntor,mpol ordering imposed)
!     as coefficients of cos(mu)*cos(nv), etc
!     Only need R, Z components(not lambda, for now anyhow)
!
      IF (ALLOCATED(rzl_local)) RETURN

      mpol1 = mpol-1
      rcc = 1;  zsc = 1
      IF (.not.lasym) THEN
         IF (lthreed) THEN
            ntmax = 2
            rss = 2;  zcs = 2
         ELSE
            ntmax = 1
         END IF
      ELSE
         IF (lthreed) THEN
            ntmax = 4
            rss = 2;  rsc = 3;  rcs = 4
            zcs = 2;  zcc = 3;  zss = 4
         ELSE
            ntmax = 2
            rsc = 2;  zcc = 2
         END IF
      END IF

!     only ALLOCATE 2*ntmax, don't need lambdas
      zsc = 1+ntmax; zcs = zcs+ntmax; zcc = zcc+ntmax; zss = zss+ntmax
      ALLOCATE(rzl_local(ns,0:ntor,0:mpol1,2*ntmax), stat=n)
      IF (n .ne. 0) STOP 'Allocation error in LoadRZL'
      rzl_local = 0

      DO mn = 1, mnmax
         m = NINT(xm(mn));  n = NINT(xn(mn))/nfp; n1 = ABS(n)
         sgn = SIGN(1, n)
         rzl_local(:,n1,m,rcc) = rzl_local(:,n1,m,rcc) + rmnc(mn,:)
         rzl_local(:,n1,m,zsc) = rzl_local(:,n1,m,zsc) + zmns(mn,:)
         IF (lthreed) THEN
            rzl_local(:,n1,m,rss) = rzl_local(:,n1,m,rss)
     1                            + sgn*rmnc(mn,:)
            rzl_local(:,n1,m,zcs) = rzl_local(:,n1,m,zcs)
     1                            - sgn*zmns(mn,:)
         END IF
         IF (lasym) THEN
            rzl_local(:,n1,m,rsc) = rzl_local(:,n1,m,rsc)
     1                            + rmns(mn,:)
            rzl_local(:,n1,m,zcc) = rzl_local(:,n1,m,zcc)
     1                            + zmnc(mn,:)
            IF (lthreed) THEN
                rzl_local(:,n1,m,rcs) = rzl_local(:,n1,m,rcs)
     1                                - sgn*rmns(mn,:)
                rzl_local(:,n1,m,zss) = rzl_local(:,n1,m,zss)
     1                                + sgn*zmnc(mn,:)
            END IF
         END IF
      END DO

!     ADDED by SAL for Vecpot calc
      IF (.not. ALLOCATED(chi))  ALLOCATE(chi(1:ns))
      DO mn = 1, ns
        chi(mn) = SUM(iotaf(1:mn)*phipf(1:mn))
      END DO

      END SUBROUTINE LoadRZL

      END MODULE read_wout_mod

