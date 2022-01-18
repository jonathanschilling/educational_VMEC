!> \file
!> \brief Write the output files of VMEC

!> \brief Write the output files of VMEC
!>
!> @param bsq
!> @param gsqrt
!> @param bsubu
!> @param bsubv
!> @param bsubs
!> @param bsupv
!> @param bsupu
!> @param rzl_array
!> @param gc_array
!> @param ier_flag error flag
SUBROUTINE wrout(bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu, rzl_array, gc_array, ier_flag)
  USE vmec_main
  USE vparams, p5 => cp5, two => c2p0
  USE vmec_input, ONLY: ns_array, ftol_array
  USE vmec_params
  USE vmercier
  USE vmec_persistent
  USE xstuff
  USE vmec_io
  USE realspace, ONLY: phip, gsqrta=>z1, z1=>z1
  USE vforces, ONLY: bsupua=>brmn_e, bsupva=>czmn_o, bsqa=>bzmn_e,  &
                     bsubsa=>armn_e, bsubua=>azmn_e, bsubva=>armn_o
  USE vacmod, ONLY: potvac
  USE ezcdf
  USE read_wout_mod, ONLY: vn_version, vn_extension, vn_mgrid,      &
    vn_magen, vn_therm, vn_gam, vn_maxr, vn_minr, vn_maxz, vn_fp,   &
    vn_radnod, vn_polmod, vn_tormod, vn_maxmod, vn_maxit,           &
    vn_asym, vn_free, vn_error, vn_aspect, vn_beta,                 &
    vn_pbeta, vn_tbeta, vn_abeta, vn_b0, vn_rbt0, vn_maxmod_nyq,    &
    vn_rbt1, vn_sgs, vn_lar, vn_modB, vn_ctor, vn_amin, vn_Rmaj,    &
    vn_vol, vn_ac, vn_ai, vn_am,                                    &
    vn_pmass_type, vn_pcurr_type, vn_piota_type,                    &
    vn_am_aux_s, vn_am_aux_f, vn_ac_aux_s, vn_ac_aux_f,             &
    vn_ai_aux_s, vn_ai_aux_f,                                       &
    vn_ftolv, vn_fsqr, vn_fsqz, vn_fsql,                            &
    vn_pmod, vn_tmod, vn_pmod_nyq, vn_tmod_nyq,                     &
    vn_racc, vn_zacs, vn_racs, vn_zacc, vn_iotaf, vn_qfact,         &
    vn_presf, vn_phi, vn_phipf, vn_jcuru, vn_jcurv, vn_iotah,       &
    vn_chi, vn_chipf,                                               &
    vn_mass, vn_presh, vn_betah, vn_buco, vn_bvco, vn_vp, vn_specw, &
    vn_phip, vn_jdotb, vn_overr, vn_bgrv, vn_merc, vn_mshear,       &
    vn_mwell, vn_mcurr, vn_mgeo, vn_equif, vn_fsq, vn_wdot,         &
    vn_extcur, vn_curlab, vn_rmnc, vn_zmns, vn_lmns, vn_gmnc,       &
    vn_bmnc, vn_bsubumnc, vn_bsubvmnc, vn_bsubsmns,                 &
    vn_bsupumnc, vn_bsupvmnc, vn_rmns, vn_zmnc, vn_lmnc, vn_gmns,   &
    vn_bmns, vn_bsubumns, vn_bsubvmns, vn_bsubsmnc, vn_bsupumns,    &
    vn_bsupvmns, vn_rbc, vn_zbs, vn_rbs, vn_zbc, vn_potvac,         &
    ln_version, ln_extension, ln_mgrid,                             &
    ln_magen, ln_therm, ln_gam, ln_maxr, ln_minr, ln_maxz, ln_fp,   &
    ln_radnod, ln_polmod, ln_tormod, ln_maxmod, ln_maxit,           &
    ln_asym, ln_free, ln_error, ln_aspect, ln_beta,                 &
    ln_pbeta, ln_tbeta, ln_abeta, ln_b0, ln_rbt0, ln_maxmod_nyq,    &
    ln_rbt1, ln_sgs, ln_lar, ln_modB, ln_ctor, ln_amin, ln_Rmaj,    &
    ln_mse, ln_thom, ln_next,                                       &
    ln_pmod, ln_tmod, ln_pmod_nyq, ln_tmod_nyq, ln_racc, ln_zacs,   &
    ln_racs, ln_zacc, ln_iotaf, ln_qfact, ln_am, ln_ac, ln_ai,      &
    ln_pmass_type, ln_pcurr_type, ln_piota_type,                    &
    ln_am_aux_s, ln_am_aux_f, ln_ac_aux_s, ln_ac_aux_f,             &
    ln_ai_aux_s, ln_ai_aux_f, ln_chi, ln_chipf,                     &
    ln_presf, ln_phi, ln_phipf, ln_jcuru, ln_jcurv, ln_iotah,       &
    ln_mass, ln_presh, ln_betah, ln_buco, ln_bvco, ln_vp, ln_specw, &
    ln_vol, ln_phip, ln_jdotb, ln_bgrv, ln_merc, ln_mshear,         &
    ln_mwell, ln_mcurr, ln_mgeo, ln_equif, ln_fsq, ln_wdot,         &
    ln_extcur, ln_curlab, ln_rmnc, ln_zmns, ln_lmns, ln_gmnc,       &
    ln_bmnc, ln_bsubumnc, ln_bsubvmnc, ln_bsubsmns,                 &
    ln_bsupumnc, ln_bsupvmnc, ln_rmns, ln_zmnc, ln_lmnc, ln_gmns,   &
    ln_bmns, ln_bsubumns, ln_bsubvmns, ln_bsubsmnc, ln_bsupumns,    &
    ln_bsupvmns, ln_rbc, ln_zbs, ln_rbs, ln_zbc, ln_potvac

  USE safe_open_mod
  USE mgrid_mod
  IMPLICIT NONE

  INTEGER, INTENT(in) :: ier_flag

  ! reverse ns, mnmax for backwards compatibility
  REAL(rprec), DIMENSION(mnmax,ns,3*MAX(ntmax/2,1)), INTENT(inout), TARGET :: rzl_array, gc_array
  REAL(rprec), DIMENSION(ns,nznt), INTENT(inout) :: bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu

  REAL(rprec) :: qfact(ns)

  CHARACTER(LEN=*), PARAMETER, DIMENSION(1) :: &
               r1dim = (/'radius'/),           &
               mn1dim = (/'mn_mode'/),         &
               mn2dim = (/'mn_mode_nyq'/),     &
               mnpddim = (/'mnpd'/),           &
               currg = (/'ext_current'/),      &
               currl = (/'current_label'/)
  CHARACTER(LEN=*), PARAMETER, DIMENSION(2) ::          &
               r2dim = (/'mn_mode',    'radius '    /), &
               r3dim = (/'mn_mode_nyq','radius     '/)

  INTEGER :: j, js, jlk, mn, lk,                              &
             m, n, k, iwout0, n1, nwout, istat, indx1(1)
  REAL(rprec) :: dmult, tcosi, tsini, vversion, sgn, tmult,         &
                 ftolx1
  REAL(rprec), POINTER, DIMENSION(:,:) :: rmnc, rmns, zmns, zmnc, lmns, lmnc
  REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: gmnc, bmnc, gmns, bmns, &
    bsubumnc, bsubvmnc, bsubsmns, bsubumns, bsubvmns, bsubsmnc
  REAL(rprec), DIMENSION(mnmax) :: rmnc1, zmns1, lmns1,             &
     rmns1, zmnc1, lmnc1, bmodmn, bmodmn1
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: gmn, bmn,               &
     bsubumn, bsubvmn, bsubsmn, bsupumn, bsupvmn
  CHARACTER(LEN=120) :: wout_file
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: xfinal

  ! ELIMINATE THESE EVENTUALLY
  REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: bsupumnc, bsupumns, bsupvmnc, bsupvmns



  ! THIS SUBROUTINE CREATES THE FILE WOUT.
  ! IT CONTAINS THE CYLINDRICAL COORDINATE SPECTRAL COEFFICIENTS
  ! RMN,ZMN (full), LMN (half_mesh - CONVERTED FROM INTERNAL full REPRESENTATION),
  ! AS WELL AS COEFFICIENTS (ON NYQ MESH) FOR COMPUTED QUANTITIES:
  ! BSQ, BSUPU,V, BSUBU,V, GSQRT (HALF); BSUBS (FULL-CONVERTED IN JXBFORCE)



  ! Pointer assignments for storage arrays
  n1 = MAX(1,ntmax/2)
  rmnc => rzl_array(:,:,1)            ! store COS(mu-nv) components
  zmns => rzl_array(:,:,1+n1)         ! store SIN(mu-nv)
  lmns => rzl_array(:,:,1+2*n1)       ! store SIN(mu-nv)

  IF (lasym) THEN
     rmns => gc_array(:,:,1)          ! store SIN(mu-nv)
     zmnc => gc_array(:,:,1+n1)       ! store COS(mu-nv)
     lmnc => gc_array(:,:,1+2*n1)     ! store COS(mu-nv)
  END IF

  ALLOCATE (gmn(mnmax_nyq), bmn(mnmax_nyq),                       &
     bsubumn(mnmax_nyq), bsubvmn(mnmax_nyq), bsubsmn(mnmax_nyq), &
     bsupumn(mnmax_nyq), bsupvmn(mnmax_nyq), stat=istat)

  ALLOCATE (gmnc(mnmax_nyq,ns), bmnc(mnmax_nyq,ns),               &
            bsubumnc(mnmax_nyq,ns), bsubvmnc(mnmax_nyq,ns),       &
            bsubsmns(mnmax_nyq,ns), bsupumnc(mnmax_nyq,ns),       &
            bsupvmnc(mnmax_nyq,ns), stat=istat)
  IF (lasym) THEN
    ALLOCATE (gmns(mnmax_nyq,ns), bmns(mnmax_nyq,ns),             &
              bsubumns(mnmax_nyq,ns), bsubvmns(mnmax_nyq,ns),     &
              bsubsmnc(mnmax_nyq,ns), bsupumns(mnmax_nyq,ns),     &
              bsupvmns(mnmax_nyq,ns), stat=istat)
  END IF
  IF (istat .ne. 0) STOP 'Error allocating arrays in VMEC WROUT'

  ! IF (nextcur .eq. 0) THEN
  !    DO j = SIZE(extcur), 1, -1
  !      IF (extcur(j) .ne. zero) THEN
  !          nextcur = j
  !          EXIT
  !       END IF
  !    END DO
  ! END IF

  ! ftol info evaluated here!
  indx1=MAXLOC(ns_array)
  ftolx1=ftol_array(indx1(1))

  ! NYQUIST FREQUENCY REQUIRES FACTOR OF 1/2
  IF (mnyq .ne. 0) cosmui(:,mnyq) = p5*cosmui(:,mnyq)
  IF (nnyq .ne. 0) cosnv (:,nnyq) = p5*cosnv (:,nnyq)

  ! use wout_file as temporary storage for parsing version number to REAL
  wout_file = version_
  READ (wout_file, *) vversion

  wout_file = 'wout_' // TRIM(input_extension) // '.nc'
  CALL cdf_open(nwout, wout_file, 'w', iwout0)
  IF (iwout0 .ne. 0) STOP 'Error opening wout.nc file VMEC WROUT'

  !================================
  ! Define Variables
  !================================
  ! Scalars
  CALL cdf_define(nwout, vn_version, vversion)
  CALL cdf_define(nwout, vn_extension, input_extension)
  CALL cdf_define(nwout, vn_mgrid, mgrid_file)
  CALL cdf_define(nwout, vn_pcurr_type, pcurr_type)
  CALL cdf_define(nwout, vn_pmass_type, pmass_type)
  CALL cdf_define(nwout, vn_piota_type, piota_type)
  CALL cdf_define(nwout, vn_magen, wb)
  CALL cdf_define(nwout, vn_therm, wp)
  CALL cdf_define(nwout, vn_gam, gamma)
  CALL cdf_define(nwout, vn_maxr, rmax_surf)
  CALL cdf_define(nwout, vn_minr, rmin_surf)
  CALL cdf_define(nwout, vn_maxz, zmax_surf)
  CALL cdf_define(nwout, vn_fp, nfp)
  CALL cdf_define(nwout, vn_radnod, ns)
  CALL cdf_define(nwout, vn_polmod, mpol)
  CALL cdf_define(nwout, vn_tormod, ntor)
  CALL cdf_define(nwout, vn_maxmod, mnmax)
  CALL cdf_define(nwout, vn_maxmod_nyq, mnmax_nyq)
  CALL cdf_define(nwout, vn_maxit, iter2)
  CALL cdf_define(nwout, vn_asym, lasym)
  CALL cdf_define(nwout, vn_free, lfreeb)
  CALL cdf_define(nwout, vn_error, ier_flag)
  CALL cdf_define(nwout, vn_aspect, aspect)
  CALL cdf_define(nwout, vn_beta, betatot)
  CALL cdf_define(nwout, vn_pbeta, betapol)
  CALL cdf_define(nwout, vn_tbeta, betator)
  CALL cdf_define(nwout, vn_abeta, betaxis)
  CALL cdf_define(nwout, vn_b0, b0)
  CALL cdf_define(nwout, vn_rbt0, rbtor0)
  CALL cdf_define(nwout, vn_rbt1, rbtor)
  CALL cdf_define(nwout, vn_sgs, NINT(signgs))
  CALL cdf_define(nwout, vn_lar, IonLarmor)
  CALL cdf_define(nwout, vn_modB, volAvgB)
  CALL cdf_define(nwout, vn_ctor, ctor)
  CALL cdf_define(nwout, vn_amin, Aminor_p)
  CALL cdf_define(nwout, vn_Rmaj, Rmajor_p)
  CALL cdf_define(nwout, vn_vol, volume_p)
  CALL cdf_define(nwout, vn_ftolv, ftolx1)
  CALL cdf_define(nwout, vn_fsqr, fsqr)
  CALL cdf_define(nwout, vn_fsqz, fsqz)
  CALL cdf_define(nwout, vn_fsql, fsql)

  CALL cdf_define(nwout, vn_nextcur, nextcur)
  CALL cdf_define(nwout, vn_extcur, extcur(1:nextcur), dimname=currg)
  CALL cdf_define(nwout, vn_mgmode, mgrid_mode)

  ! 1D Arrays
  CALL cdf_define(nwout, vn_pmod, xm, dimname=mn1dim)
  CALL cdf_setatt(nwout, vn_pmod, ln_pmod)
  CALL cdf_define(nwout, vn_tmod, xn, dimname=mn1dim)
  CALL cdf_setatt(nwout, vn_tmod, ln_tmod)
  CALL cdf_define(nwout, vn_pmod_nyq, xm_nyq, dimname=mn2dim)
  CALL cdf_setatt(nwout, vn_pmod_nyq, ln_pmod_nyq)
  CALL cdf_define(nwout, vn_tmod_nyq, xn_nyq, dimname=mn2dim)
  CALL cdf_setatt(nwout, vn_tmod_nyq, ln_tmod_nyq)

  CALL cdf_define(nwout, vn_racc, raxis_cc(0:ntor), dimname=(/'n_tor'/))
  CALL cdf_setatt(nwout, vn_racc, ln_racc)
  CALL cdf_define(nwout, vn_zacs, zaxis_cs(0:ntor), dimname=(/'n_tor'/))
  CALL cdf_setatt(nwout, vn_zacs, ln_zacs)
  IF (lasym) THEN
     CALL cdf_define(nwout, vn_racs, raxis_cs(0:ntor), dimname=(/'n_tor'/))
     CALL cdf_setatt(nwout, vn_racs, ln_racs)
     CALL cdf_define(nwout, vn_zacc, zaxis_cc(0:ntor), dimname=(/'n_tor'/))
     CALL cdf_setatt(nwout, vn_zacc, ln_zacc)
  END IF

  j = SIZE(am)-1
  CALL cdf_define(nwout, vn_am, am(0:j), dimname=(/'preset'/))
  j = SIZE(ac)-1
  CALL cdf_define(nwout, vn_ac, ac(0:j), dimname=(/'preset'/))
  j = SIZE(ai)-1
  CALL cdf_define(nwout, vn_ai, ai(0:j), dimname=(/'preset'/))

  j = SIZE(am_aux_s)
  CALL cdf_define(nwout, vn_am_aux_s, am_aux_s(1:j), dimname=(/'ndfmax'/))
  j = SIZE(am_aux_f)
  CALL cdf_define(nwout, vn_am_aux_f, am_aux_f(1:j), dimname=(/'ndfmax'/))
  j = SIZE(ai_aux_s)
  CALL cdf_define(nwout, vn_ai_aux_s, ai_aux_s(1:j), dimname=(/'ndfmax'/))
  j = SIZE(ai_aux_f)
  CALL cdf_define(nwout, vn_ai_aux_f, ai_aux_f(1:j), dimname=(/'ndfmax'/))
  j = SIZE(ac_aux_s)
  CALL cdf_define(nwout, vn_ac_aux_s, ac_aux_s(1:j), dimname=(/'ndfmax'/))
  j = SIZE(ac_aux_f)
  CALL cdf_define(nwout, vn_ac_aux_f, ac_aux_f(1:j), dimname=(/'ndfmax'/))


  CALL cdf_define(nwout, vn_iotaf, iotaf(1:ns), dimname=r1dim)
  CALL cdf_setatt(nwout, vn_iotaf, ln_iotaf)

  qfact=HUGE(qfact)
  WHERE (iotaf(1:ns) .NE. zero) qfact=one/iotaf(1:ns)
  CALL cdf_define(nwout, vn_qfact, qfact(1:ns), dimname=r1dim)
  CALL cdf_setatt(nwout, vn_iotaf, ln_qfact)

  CALL cdf_define(nwout, vn_presf, presf,       dimname=r1dim)
  CALL cdf_setatt(nwout, vn_presf, ln_presf, units='Pa')
  CALL cdf_define(nwout, vn_phi,   phi,         dimname=r1dim)
  CALL cdf_setatt(nwout, vn_phi,   ln_phi, units='wb')
  CALL cdf_define(nwout, vn_phipf, phipf,       dimname=r1dim)
  CALL cdf_setatt(nwout, vn_phipf, ln_phipf)
  CALL cdf_define(nwout, vn_chi,   chi,         dimname=r1dim)
  CALL cdf_setatt(nwout, vn_chi,   ln_chi, units='wb')
  CALL cdf_define(nwout, vn_chipf, phipf,       dimname=r1dim) ! TODO: wrong quantity !!!
  CALL cdf_setatt(nwout, vn_chipf, ln_chipf)
  CALL cdf_define(nwout, vn_jcuru, jcuru,       dimname=r1dim)
  CALL cdf_define(nwout, vn_jcurv, jcurv,       dimname=r1dim)

  CALL cdf_define(nwout, vn_iotah, iotas(1:ns), dimname=r1dim)
  CALL cdf_setatt(nwout, vn_iotah, ln_iotah)
  CALL cdf_define(nwout, vn_mass,  mass,        dimname=r1dim)
  CALL cdf_setatt(nwout, vn_mass,  ln_mass)
  CALL cdf_define(nwout, vn_presh, pres(1:ns),  dimname=r1dim)
  CALL cdf_setatt(nwout, vn_presh, ln_presh, units='Pa')
  CALL cdf_define(nwout, vn_betah, beta_vol,    dimname=r1dim)
  CALL cdf_define(nwout, vn_buco,  buco,        dimname=r1dim)
  CALL cdf_define(nwout, vn_bvco,  bvco,        dimname=r1dim)
  CALL cdf_define(nwout, vn_vp,    vp(1:ns),    dimname=r1dim)
  CALL cdf_define(nwout, vn_specw, specw,       dimname=r1dim) ! TODO: specw on full grid ?
  CALL cdf_define(nwout, vn_phip,  phips(1:ns), dimname=r1dim)
  CALL cdf_define(nwout, vn_overr, overr(1:ns), dimname=r1dim)

  CALL cdf_define(nwout, vn_jdotb, jdotb, dimname=r1dim)
  CALL cdf_define(nwout, vn_bgrv, bdotgradv, dimname=r1dim)

  CALL cdf_define(nwout, vn_merc,   Dmerc,  dimname=r1dim)
  CALL cdf_define(nwout, vn_mshear, Dshear, dimname=r1dim)
  CALL cdf_define(nwout, vn_mwell,  Dwell,  dimname=r1dim)
  CALL cdf_define(nwout, vn_mcurr,  Dcurr,  dimname=r1dim)
  CALL cdf_define(nwout, vn_mgeo,   Dgeod,  dimname=r1dim)
  CALL cdf_define(nwout, vn_equif,  equif,  dimname=r1dim)

  IF (lfreeb .and. nextcur.gt.0 .and. ALLOCATED(curlabel)) THEN
     CALL cdf_define(nwout, vn_curlab, curlabel(1:nextcur), dimname=currl)
     ! SAL for potvac
     CALL cdf_define(nwout, vn_potvac, potvac, dimname=mnpddim)
  ENDIF

  ! 2D Arrays
  CALL cdf_define(nwout, vn_rmnc, rmnc, dimname=r2dim)
  CALL cdf_setatt(nwout, vn_rmnc, ln_rmnc, units='m')
  CALL cdf_define(nwout, vn_zmns, zmns, dimname=r2dim)
  CALL cdf_setatt(nwout, vn_zmns, ln_zmns, units='m')
  CALL cdf_define(nwout, vn_lmns, lmns, dimname=r2dim)
  CALL cdf_setatt(nwout, vn_lmns, ln_lmns)
  CALL cdf_define(nwout, vn_gmnc, gmnc, dimname=r3dim)
  CALL cdf_setatt(nwout, vn_gmnc, ln_gmnc)
  CALL cdf_define(nwout, vn_bmnc, bmnc, dimname=r3dim)
  CALL cdf_setatt(nwout, vn_bmnc, ln_bmnc)
  CALL cdf_define(nwout, vn_bsubumnc, bsubumnc, dimname=r3dim)
  CALL cdf_setatt(nwout, vn_bsubumnc, ln_bsubumnc)
  CALL cdf_define(nwout, vn_bsubvmnc, bsubvmnc, dimname=r3dim)
  CALL cdf_setatt(nwout, vn_bsubvmnc, ln_bsubvmnc)
  CALL cdf_define(nwout, vn_bsubsmns, bsubsmns, dimname=r3dim)
  CALL cdf_setatt(nwout, vn_bsubsmns, ln_bsubsmns)

  ! ELIMINATE THESE EVENTUALLY: DON'T NEED THEM - CAN COMPUTE FROM GSQRT
  CALL cdf_define(nwout, vn_bsupumnc, bsupumnc, dimname=r3dim)
  CALL cdf_define(nwout, vn_bsupvmnc, bsupvmnc, dimname=r3dim)

  IF (lasym) then
     CALL cdf_define(nwout, vn_rmns, rmns, dimname=r2dim)
     CALL cdf_setatt(nwout, vn_rmns, ln_rmns, units='m')
     CALL cdf_define(nwout, vn_zmnc, zmnc, dimname=r2dim)
     CALL cdf_setatt(nwout, vn_zmnc, ln_zmnc, units='m')
     CALL cdf_define(nwout, vn_lmnc, lmnc, dimname=r2dim)
     CALL cdf_setatt(nwout, vn_lmnc, ln_lmnc)
     CALL cdf_define(nwout, vn_gmns, gmns, dimname=r3dim)
     CALL cdf_setatt(nwout, vn_gmns, ln_gmns)
     CALL cdf_define(nwout, vn_bmns, bmns, dimname=r3dim)
     CALL cdf_setatt(nwout, vn_bmns, ln_bmns)
     CALL cdf_define(nwout, vn_bsubumns, bsubumns, dimname=r3dim)
     CALL cdf_setatt(nwout, vn_bsubumns, ln_bsubumns)
     CALL cdf_define(nwout, vn_bsubvmns, bsubvmns, dimname=r3dim)
     CALL cdf_setatt(nwout, vn_bsubvmns, ln_bsubvmns)
     CALL cdf_define(nwout, vn_bsubsmnc, bsubsmnc, dimname=r3dim)
     CALL cdf_setatt(nwout, vn_bsubsmnc, ln_bsubsmnc)

     ! ELIMINATE THESE EVENTUALLY: DON'T NEED THEM
     CALL cdf_define(nwout, vn_bsupumns, bsupumns, dimname=r3dim)
     CALL cdf_define(nwout, vn_bsupvmns, bsupvmns, dimname=r3dim)
  end if

  !================================
  ! Write Variables
  !================================
  ! Scalars
  CALL cdf_write(nwout, vn_version, vversion)
  CALL cdf_write(nwout, vn_extension, input_extension)
  CALL cdf_write(nwout, vn_mgrid, mgrid_file)
  CALL cdf_write(nwout, vn_pcurr_type, pcurr_type)
  CALL cdf_write(nwout, vn_piota_type, piota_type)
  CALL cdf_write(nwout, vn_pmass_type, pmass_type)
  CALL cdf_write(nwout, vn_magen, wb)
  CALL cdf_write(nwout, vn_therm, wp)
  CALL cdf_write(nwout, vn_gam, gamma)
  CALL cdf_write(nwout, vn_maxr, rmax_surf)
  CALL cdf_write(nwout, vn_minr, rmin_surf)
  CALL cdf_write(nwout, vn_maxz, zmax_surf)
  CALL cdf_write(nwout, vn_fp, nfp)
  CALL cdf_write(nwout, vn_radnod, ns)
  CALL cdf_write(nwout, vn_polmod, mpol)
  CALL cdf_write(nwout, vn_tormod, ntor)
  CALL cdf_write(nwout, vn_maxmod, mnmax)
  CALL cdf_write(nwout, vn_maxmod_nyq, mnmax_nyq)
  CALL cdf_write(nwout, vn_maxit, iter2)
  CALL cdf_write(nwout, vn_asym, lasym)
  CALL cdf_write(nwout, vn_free, lfreeb)
  CALL cdf_write(nwout, vn_error, ier_flag)

  CALL cdf_write(nwout, vn_aspect, aspect)
  CALL cdf_write(nwout, vn_beta, betatot)
  CALL cdf_write(nwout, vn_pbeta, betapol)
  CALL cdf_write(nwout, vn_tbeta, betator)
  CALL cdf_write(nwout, vn_abeta, betaxis)
  CALL cdf_write(nwout, vn_b0, b0)
  CALL cdf_write(nwout, vn_rbt0, rbtor0)
  CALL cdf_write(nwout, vn_rbt1, rbtor)
  CALL cdf_write(nwout, vn_sgs, NINT(signgs))
  CALL cdf_write(nwout, vn_lar, IonLarmor)
  CALL cdf_write(nwout, vn_modB, volAvgB)
  CALL cdf_write(nwout, vn_ctor, ctor/mu0)
  CALL cdf_write(nwout, vn_amin, Aminor_p)
  CALL cdf_write(nwout, vn_rmaj, Rmajor_p)
  CALL cdf_write(nwout, vn_vol, volume_p)
  CALL cdf_write(nwout, vn_ftolv, ftolx1)
  CALL cdf_write(nwout, vn_fsql, fsql)
  CALL cdf_write(nwout, vn_fsqr, fsqr)
  CALL cdf_write(nwout, vn_fsqz, fsqz)

  CALL cdf_write(nwout, vn_nextcur, nextcur)
  IF (nextcur .gt. 0) THEN
     CALL cdf_write(nwout, vn_extcur, extcur(1:nextcur))
     CALL cdf_write(nwout, vn_mgmode, mgrid_mode)
  ENDIF
  IF (lfreeb) THEN
     ! TODO: write current labels
     CALL cdf_write(nwout, vn_potvac, potvac)
  END IF

  ! 1D Arrays
  CALL cdf_write(nwout, vn_pmod, xm)
  CALL cdf_write(nwout, vn_tmod, xn)
  CALL cdf_write(nwout, vn_pmod_nyq, xm_nyq)
  CALL cdf_write(nwout, vn_tmod_nyq, xn_nyq)

  ALLOCATE (xfinal(neqs), stat=js) ! re-use js for allocation return code
  IF (js .ne. 0) STOP 'Allocation error for xfinal in WROUT!'
  xfinal = xc ! ignore passed rzl_array !!!
  ! --> rzl_array is re-used to store temporary  symmetric rmnc, ... !
  ! -->  gc_array is re-used to store temporary asymmetric rmns, ... ! (for lasym)

  ! MUST CONVERT m=1 MODES... FROM INTERNAL TO PHYSICAL FORM
  ! Extrapolation of m=0 Lambda (cs) modes, which are not evolved at j=1, done in CONVERT

  IF (lthreed) CALL convert_sym  (xfinal(1+mns*(rss-1)), xfinal(1+irzloff+mns*(zcs-1)))
  IF (lasym)   CALL convert_asym (xfinal(1+mns*(rsc-1)), xfinal(1+irzloff+mns*(zcc-1)))

  ! CONVERT TO rmnc, zmns, lmns, etc EXTERNAL representation (without internal mscale, nscale)
  ! IF B^v ~ phip + lamu, MUST DIVIDE BY phipf(js) below to maintain old-style format
  RADIUS1: DO js = 1, ns

     CALL convert (rmnc1, zmns1, lmns1, rmns1, zmnc1, lmnc1, xfinal, js)

     rmnc(:,js) = rmnc1(:)
     zmns(:,js) = zmns1(:)
     lmns(:,js) = (lmns1(:)/phipf(js)) * lamscale
     IF (lasym) THEN
        rmns(:,js) = rmns1(:)
        zmnc(:,js) = zmnc1(:)
        lmnc(:,js) = (lmnc1(:)/phipf(js)) * lamscale
     END IF

  END DO RADIUS1

  DEALLOCATE (xfinal)

  ! INTERPOLATE LAMBDA ONTO HALF-MESH FOR BACKWARDS CONSISTENCY WITH EARLIER VERSIONS OF VMEC
  ! AND SMOOTHS POSSIBLE UNPHYSICAL "WIGGLE" ON RADIAL MESH
  WHERE (NINT(xm) .le. 1) lmns(:,1) = lmns(:,2)
  DO js = ns,2,-1
     WHERE (MOD(NINT(xm),2) .eq. 0)
        lmns(:,js) = p5*(lmns(:,js) + lmns(:,js-1))
     ELSEWHERE
        lmns(:,js) = p5*(sm(js)*lmns(:,js) + sp(js-1)*lmns(:,js-1))
     END WHERE
  END DO

  lmns(:,1) = 0
  raxis_cc(0:ntor) = rmnc(1:ntor+1,1)
  zaxis_cs(0:ntor) = zmns(1:ntor+1,1)

  IF (lasym) then
     WHERE (NINT(xm) .le. 1) lmnc(:,1) = lmnc(:,2)
     DO js = ns,2,-1
        WHERE (MOD(NINT(xm),2) .eq. 0)
           lmnc(:,js) = p5*(lmnc(:,js) + lmnc(:,js-1))
        ELSEWHERE
           lmnc(:,js) = p5*(sm(js)*lmnc(:,js) + sp(js-1)*lmnc(:,js-1))
        END WHERE
     END DO

     lmnc(:,1) = 0;
     raxis_cs(0:ntor) = rmns(1:ntor+1,1)
     zaxis_cc(0:ntor) = zmnc(1:ntor+1,1)
  end if

  ! COMPUTE |B| = SQRT(|B|**2) and store in bsq, bsqa
  DO js = 2, ns
     bsq(js,:nznt) = SQRT(2*ABS(bsq(js,:nznt)-pres(js)))
  END DO

  tmult = p5/r0scale**2
  !SPH: FIXED THIS 03-05-07 TO CALL symmetrization routine
  IF (lasym) THEN
     ! Changed integration norm in fixaray, SPH012314
     tmult = 2*tmult
     bsubs(1,:) = 0
     CALL symoutput (bsq,   gsqrt,  bsubu,  bsubv,  bsupu,   &
                     bsupv,  bsubs,                          &
                     bsqa,  gsqrta, bsubua, bsubva, bsupua,  &
                     bsupva, bsubsa )
  END IF

  ! Fourier-transform derived quantities for each surface individually
  RADIUS2: DO js = 2, ns
     gmn = 0
     bmn = 0
     bsubumn = 0
     bsubvmn = 0
     bsubsmn = 0
     bsupumn = 0
     bsupvmn = 0

     MN2: DO mn = 1, mnmax_nyq
        n = NINT(xn_nyq(mn))/nfp
        m = NINT(xm_nyq(mn))
        n1 = ABS(n)
        dmult = mscale(m)*nscale(n1)*tmult
        IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
        sgn = SIGN(1, n)
        lk = 0
        DO j = 1, ntheta2
           DO k = 1, nzeta
              lk = lk + 1

              ! cos(mu - nv)
              tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) + sgn*sinmui(j,m)*sinnv(k,n1))

              ! sin(mu - nv)
              tsini = dmult*(sinmui(j,m)*cosnv(k,n1) - sgn*cosmui(j,m)*sinnv(k,n1))

              bmn(mn)     = bmn(mn)     + tcosi*bsq(js,lk)
              gmn(mn)     = gmn(mn)     + tcosi*gsqrt(js,lk)
              bsubumn(mn) = bsubumn(mn) + tcosi*bsubu(js,lk)
              bsubvmn(mn) = bsubvmn(mn) + tcosi*bsubv(js,lk)
              bsubsmn(mn) = bsubsmn(mn) + tsini*bsubs(js,lk)
              bsupumn(mn) = bsupumn(mn) + tcosi*bsupu(js,lk)
              bsupvmn(mn) = bsupvmn(mn) + tcosi*bsupv(js,lk)
           END DO
        END DO
     END DO MN2

     IF (js .eq. ns/2) bmodmn  = bmn(1:mnmax)
     IF (js .eq. ns  ) bmodmn1 = bmn(1:mnmax)
     gmnc(:,js) = gmn(:)
     bmnc(:,js) = bmn(:)
     bsubumnc(:,js) = bsubumn(:)
     bsubvmnc(:,js) = bsubvmn(:)
     bsubsmns(:,js) = bsubsmn(:)
     bsupumnc(:,js) = bsupumn(:)
     bsupvmnc(:,js) = bsupvmn(:)
  END DO RADIUS2

  ! endpoint values at magnetic axis
  gmnc(:,1) = 0
  bmnc(:,1) = 0
  bsubumnc(:,1) = 0
  bsubvmnc(:,1) = 0
  bsubsmns(:,1) = 2*bsubsmns(:,2) - bsubsmns(:,3) ! extrapolation on full grid
  bsupumnc(:,1) = 0
  bsupvmnc(:,1) = 0

  IF (lasym) then
     RADIUS3: DO js = 2, ns
        gmn = 0
        bmn = 0
        bsubumn = 0
        bsubvmn = 0
        bsubsmn = 0
        bsupumn = 0
        bsupvmn = 0
        MN3: DO mn = 1, mnmax_nyq
           n = NINT(xn_nyq(mn))/nfp
           m = NINT(xm_nyq(mn))
           n1 = ABS(n)
           dmult = mscale(m)*nscale(n1)*tmult
           IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
           sgn = SIGN(1, n)
           lk = 0
           jlk = js
           DO j = 1, ntheta2
              DO k = 1, nzeta
                 lk = lk + 1
                 tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) + sgn*sinmui(j,m)*sinnv(k,n1))
                 tsini = dmult*(sinmui(j,m)*cosnv(k,n1) - sgn*cosmui(j,m)*sinnv(k,n1))
                 bmn(mn)     = bmn(mn)     + tsini*bsqa(jlk)
                 gmn(mn)     = gmn(mn)     + tsini*gsqrta(jlk,0)
                 bsubumn(mn) = bsubumn(mn) + tsini*bsubua(jlk)
                 bsubvmn(mn) = bsubvmn(mn) + tsini*bsubva(jlk)
                 bsubsmn(mn) = bsubsmn(mn) + tcosi*bsubsa(jlk)
                 bsupumn(mn) = bsupumn(mn) + tsini*bsupua(jlk)
                 bsupvmn(mn) = bsupvmn(mn) + tsini*bsupva(jlk)
                 jlk = jlk+ns
              END DO
           END DO
        END DO MN3

        gmns(:,js)     = gmn(:)
        bmns(:,js)     = bmn(:)
        bsubumns(:,js) = bsubumn(:)
        bsubvmns(:,js) = bsubvmn(:)
        bsubsmnc(:,js) = bsubsmn(:)
        bsupumns(:,js) = bsupumn(:)
        bsupvmns(:,js) = bsupvmn(:)
     END DO RADIUS3

     gmns(:,1) = 0
     bmns(:,1) = 0
     bsubumns(:,1) = 0
     bsubvmns(:,1) = 0
     bsubsmnc(:,1) = 2*bsubsmnc(:,2) - bsubsmnc(:,3)
     bsupumns(:,1) = 0
     bsupvmns(:,1) = 0
  end if

  ! WRITE OUT ARRAYS
  CALL cdf_write(nwout, vn_racc, raxis_cc(0:ntor))
  CALL cdf_write(nwout, vn_zacs, zaxis_cs(0:ntor))
  CALL cdf_write(nwout, vn_rmnc, rmnc)
  CALL cdf_write(nwout, vn_zmns, zmns)
  CALL cdf_write(nwout, vn_lmns, lmns)
  CALL cdf_write(nwout, vn_gmnc, gmnc)              !Half mesh
  CALL cdf_write(nwout, vn_bmnc, bmnc)              !Half mesh
  CALL cdf_write(nwout, vn_bsubumnc, bsubumnc)      !Half mesh
  CALL cdf_write(nwout, vn_bsubvmnc, bsubvmnc)      !Half mesh
  CALL cdf_write(nwout, vn_bsubsmns, bsubsmns)      !Full mesh

  ! GET RID OF THESE EVENTUALLY: DON'T NEED THEM (can express in terms of lambdas)
  CALL cdf_write(nwout, vn_bsupumnc, bsupumnc)
  CALL cdf_write(nwout, vn_bsupvmnc, bsupvmnc)

  ! FULL-MESH quantities
  ! NOTE: jdotb is in units_of_A (1/mu0 incorporated in jxbforce...)
  ! prior to version 6.00, this was output in internal VMEC units...
  j = SIZE(am)-1
  CALL cdf_write(nwout, vn_am, am(0:j))
  j = SIZE(ac)-1
  CALL cdf_write(nwout, vn_ac, ac(0:j))
  j = SIZE(ai)-1
  CALL cdf_write(nwout, vn_ai, ai(0:j))

  j = SIZE(am_aux_s)
  CALL cdf_write(nwout, vn_am_aux_s, am_aux_s(1:j))
  j = SIZE(am_aux_f)
  CALL cdf_write(nwout, vn_am_aux_f, am_aux_f(1:j))
  j = SIZE(ac_aux_s)
  CALL cdf_write(nwout, vn_ac_aux_s, ac_aux_s(1:j))
  j = SIZE(ac_aux_f)
  CALL cdf_write(nwout, vn_ac_aux_f, ac_aux_f(1:j))
  j = SIZE(ai_aux_s)
  CALL cdf_write(nwout, vn_ai_aux_s, ai_aux_s(1:j))
  j = SIZE(ai_aux_f)
  CALL cdf_write(nwout, vn_ai_aux_f, ai_aux_f(1:j))

  CALL cdf_write(nwout, vn_iotaf, iotaf(1:ns))
  CALL cdf_write(nwout, vn_qfact, qfact(1:ns))
  CALL cdf_write(nwout, vn_presf, presf/mu0) ! NOTE: scaling !!!
  CALL cdf_write(nwout, vn_phi, phi)
  CALL cdf_write(nwout, vn_phipf, twopi*signgs*phipf) ! NOTE: scaling !!!
  CALL cdf_write(nwout, vn_chi, chi)
  CALL cdf_write(nwout, vn_chipf, twopi*signgs*chipf) ! NOTE: scaling !!!
  CALL cdf_write(nwout, vn_jcuru, jcuru/mu0) ! NOTE: scaling !!!
  CALL cdf_write(nwout, vn_jcurv, jcurv/mu0) ! NOTE: scaling !!!
  CALL cdf_write(nwout, vn_jdotb, jdotb)
  CALL cdf_write(nwout, vn_bgrv, bdotgradv)

  ! HALF-MESH quantities
  iotas(1) = 0
  mass(1) = 0
  pres(1) = 0
  phip(1) = 0
  buco(1) = 0
  bvco(1) = 0
  vp(1) = 0
  overr(1) = 0
  specw(1) = 1
  beta_vol(1) = 0
  CALL cdf_write(nwout, vn_iotah, iotas(1:ns))
  CALL cdf_write(nwout, vn_mass, mass/mu0) ! NOTE: scaling !!!
  CALL cdf_write(nwout, vn_presh, pres(1:ns)/mu0) ! NOTE: scaling !!!
  CALL cdf_write(nwout, vn_betah, beta_vol)
  CALL cdf_write(nwout, vn_buco, buco)
  CALL cdf_write(nwout, vn_bvco, bvco)
  CALL cdf_write(nwout, vn_vp, vp(1:ns))
  CALL cdf_write(nwout, vn_specw, specw)
  CALL cdf_write(nwout, vn_phip, phips(1:ns))
  CALL cdf_write(nwout, vn_overr, overr(1:ns))

  ! MERCIER_CRITERION
  CALL cdf_write(nwout, vn_merc, Dmerc(2:ns1))
  CALL cdf_write(nwout, vn_mshear, Dshear(2:ns1))
  CALL cdf_write(nwout, vn_mwell, Dwell(2:ns1))
  CALL cdf_write(nwout, vn_mcurr, Dcurr(2:ns1))
  CALL cdf_write(nwout, vn_mgeo, Dgeod(2:ns1))
  CALL cdf_write(nwout, vn_equif, equif(2:ns1))

  IF (lasym) THEN
     CALL cdf_write(nwout, vn_racs, raxis_cs(0:ntor))
     CALL cdf_write(nwout, vn_zacc, zaxis_cc(0:ntor))
     CALL cdf_write(nwout, vn_rmns, rmns)
     CALL cdf_write(nwout, vn_zmnc, zmnc)
     CALL cdf_write(nwout, vn_lmnc, lmnc)
     CALL cdf_write(nwout, vn_gmns, gmns)
     CALL cdf_write(nwout, vn_bmns, bmns)
     CALL cdf_write(nwout, vn_bsubumns, bsubumns)
     CALL cdf_write(nwout, vn_bsubvmns, bsubvmns)
     CALL cdf_write(nwout, vn_bsubsmnc, bsubsmnc)

     ! GET RID OF THESE EVENTUALLY: DON'T NEED THEM
     CALL cdf_write(nwout, vn_bsupumns, bsupumns)
     CALL cdf_write(nwout, vn_bsupvmns, bsupvmns)
  END IF

  CALL cdf_close(nwout)

  ! RESTORE nyq ENDPOINT VALUES
  IF (mnyq .ne. 0) cosmui(:,mnyq) = 2*cosmui(:,mnyq)
  IF (nnyq .ne. 0) cosnv (:,nnyq) = 2*cosnv (:,nnyq)

  ! DEALLOCATIONS
  IF (ALLOCATED(gmnc)) DEALLOCATE(gmnc, bmnc, bsubumnc, bsubvmnc, &
                                  bsubsmns, bsupumnc, bsupvmnc )
  IF (ALLOCATED(gmns)) DEALLOCATE(gmns, bmns, bsubumns, bsubvmns, &
                                  bsubsmnc, bsupumns, bsupvmns )
  IF (ALLOCATED(gmn))  DEALLOCATE (gmn, bmn, bsubumn, bsubvmn, &
                                   bsubsmn, bsupumn, bsupvmn, stat=istat)

  ! FREE BOUNDARY DATA
  CALL freeb_data(rmnc1, zmns1, rmns1, zmnc1, bmodmn, bmodmn1)

  rzl_array = 0

END SUBROUTINE wrout
