!> \file
!> \brief Program for computing local \f$\mathbf{K} \times \mathbf{B} = \nabla p\f$ force balance.

!> \brief Program for computing local \f$\mathbf{K} \times \mathbf{B} = \nabla p\f$ force balance.
!>
!> @param bsupu contravariant component of magnetic field \f$B^\theta\f$
!> @param bsupv contravariant component of magnetic field \f$B^\zeta\f$
!> @param bsubu covariant component of magnetic field \f$B_\theta\f$
!> @param bsubv covariant component of magnetic field \f$B_\zeta\f$
!> @param bsubsh covariant component of magnetic field \f$B_s\f$ on half grid from bss()
!> @param bsubsu tangential derivate of covariant component of magnetic field \f$\partial B_s / \partial \theta\f$ (?)
!> @param bsubsv tangential derivate of covariant component of magnetic field \f$\partial B_s / \partial \zeta\f$  (?)
!> @param gsqrt Jacobian \f$\sqrt{g}\f$
!> @param bsq modulus of magnetic field \f$|\mathbf{B}|^2\f$
!> @param itheta poloidal current
!> @param izeta toroidal current
!> @param brho radial component of magnetic field \f$B_\rho\f$ (?)
!> @param ier_flag error flag
SUBROUTINE jxbforce(bsupu, bsupv, bsubu, bsubv, bsubsh, &
                    bsubsu, bsubsv, &
                    gsqrt, bsq, itheta, izeta, brho, ier_flag)
  USE safe_open_mod
  USE vmec_main
  USE vmec_params, ONLY: signgs, mnyq, nnyq, successful_term_flag
  USE realspace,   ONLY: shalf, wint, guu, guv, gvv, r1, ru, rv, zu, zv, phip
  USE ezcdf

  use dbgout

  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,nznt), INTENT(in) :: bsupu
  REAL(rprec), DIMENSION(ns,nznt), INTENT(in) :: bsupv
  REAL(rprec), DIMENSION(ns,nznt,0:1), TARGET, INTENT(inout) :: bsubu
  REAL(rprec), DIMENSION(ns,nznt,0:1), TARGET, INTENT(inout) :: bsubv
  REAL(rprec), DIMENSION(ns,nznt), INTENT(in)  :: bsubsh
  REAL(rprec), DIMENSION(ns,nznt,0:1) :: bsubsu
  REAL(rprec), DIMENSION(ns,nznt,0:1) :: bsubsv
  REAL(rprec), DIMENSION(ns,nznt), INTENT(in) :: gsqrt
  REAL(rprec), DIMENSION(ns,nznt), INTENT(in) :: bsq
  REAL(rprec), DIMENSION(ns,nznt), INTENT(out) :: itheta
  REAL(rprec), DIMENSION(ns,nznt), INTENT(out) :: izeta
  REAL(rprec), DIMENSION(ns,nznt), INTENT(out) :: brho
  INTEGER, INTENT(in) :: ier_flag

  REAL(rprec), DIMENSION(ns,nznt), TARGET :: bsubs

  ! Prints out bsubs spectrum to fort.33
  LOGICAL, PARAMETER :: lprint = .false.

  INTEGER lk, lz, lt, k, m, js, j, n, injxbout, mparity, nznt1
  INTEGER :: njxbout = jxbout0, info

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bdotk, bsubuv, bsubvu
  REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: bsubsmn
  REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::   brhomn,           &
       bsubs3, bsubv3, bsubu3, jxb_gradp, jcrossb, sqrtg3,          &
       bsupv3, bsupu3, jsups3, jsupv3, jsupu3, jdotb_sqrtg

  REAL(rprec), POINTER :: bs1(:), bu1(:,:), bv1(:,:)

  REAL(rprec), DIMENSION(:), ALLOCATABLE     :: kperpu, kperpv,     &
      sqgb2, sqrtg, kp2, jxb, jxb2, bsupu1, bsupv1, bsubu1, bsubv1, &
      avforce, aminfor, amaxfor, toroidal_angle, phin, pprim, pprime
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE   :: bsubua, bsubva
  REAL(rprec) ::                                                    &
      bsubsmn1, bsubsmn2, bsubvmn1, bsubvmn2, bsubumn1, bsubumn2,   &
      bsubsmn3, bsubsmn4, bsubvmn3, bsubvmn4, bsubumn3, bsubumn4,   &
      dnorm1, tcos1, tcos2, tsini1, tsini2, tcosi1, tcosi2,         &
      tcosm1, tcosm2, tcosn1, tcosn2, tsinm1, tsinm2, tsin1, tsin2, &
      tsinn1, tsinn2, tjnorm, ovp, pnorm, brho00(ns)
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bsubu_s, bsubu_a, bsubv_s, bsubv_a
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsubs_s, bsubs_a
  CHARACTER(LEN=100) :: jxbout_file
  CHARACTER(LEN=100) :: legend(13)
  LOGICAL :: lprint_flag

  CHARACTER(LEN=*), PARAMETER ::                                    &
    vn_legend = 'legend',                                           &
    vn_radial_surfaces = 'radial_surfaces',                         &
    vn_poloidal_grid_points = 'poloidal_grid_points',               &
    vn_toroidal_grid_points = 'toroidal_grid_points',               &
    vn_mpol = 'mpol',                                               &
    vn_ntor = 'ntor',                                               &
    vn_phin = 'phin',                                               &
    vn_toroidal_angle = 'toroidal_angle',                           &
    vn_avforce = 'avforce',                                         &
    vn_jdotb = 'surf_av_jdotb',                                     &
    vn_sqg_bdotk = 'sqrt(g)*bdotk',                                 &
    vn_sqrtg = 'sqrt(g)',                                           &
    vn_bdotgradv = 'bdotgradv',                                     &
    vn_amaxfor = 'amaxfor',                                         &
    vn_aminfor = 'aminfor',                                         &
    vn_pprime = 'pprime',                                           &
    vn_jsupu = 'jsupu',                                             &
    vn_jsupv = 'jsupv',                                             &
    vn_jsups = 'jsups',                                             &
    vn_bsupu = 'bsupu',                                             &
    vn_bsupv = 'bsupv',                                             &
    vn_jcrossb = 'jcrossb',                                         &
    vn_jxb_gradp = 'jxb_gradp',                                     &
    vn_bsubu = 'bsubu',                                             &
    vn_bsubv = 'bsubv',                                             &
    vn_bsubs = 'bsubs'



   ! PROGRAM FOR COMPUTING LOCAL KXB = grad-p FORCE BALANCE
   !
   ! K = CURL(B) is the "effective" current (=J for isotropic pressure)
   ! Compute u (=theta), v (=zeta) derivatives of B sub s

  ! fixup input coming from bss: zero terms past magnetic axis
  ! NOTE: FIXME: This is actually the target array for the interpolation to the half-grid !!!
  ! TODO: This should probably be `bsubsh(1,:) = 0` instead !!!
  bsubs(1,:) = 0

  lprint_flag = (ier_flag.eq.successful_term_flag)
  IF (lprint_flag) THEN
     jxbout_file = 'jxbout_'//TRIM(input_extension)//'.nc'

     CALL cdf_open(njxbout,jxbout_file,'w',injxbout)
     IF (injxbout .ne. 0) THEN
        PRINT *,' Error opening JXBOUT file in jxbforce'
        RETURN
     END IF

     legend(1) = " S = normalized toroidal flux (0 - 1)"
     IF (lasym) THEN
        legend(2) = " U = VMEC poloidal angle (0 - 2*pi, FULL period)"
     ELSE
        legend(2) = " U = VMEC poloidal angle (0 - pi, HALF a period)"
     END IF
     legend(3) = " V = VMEC (geometric) toroidal angle (0 - 2*pi)"
     legend(4) = " SQRT(g') = |SQRT(g-VMEC)| / VOL': Cylindrical-to-s,u,v Jacobian normed to volume derivative"
     legend(5) = " VOL = Int_s'=0,s Int_u Int_v |SQRT(g_VMEC)| : plasma volume  enclosed by surface s'=s"
     legend(6) = " VOL' = d(VOL)/ds: differential volume element"
     legend(7) = " Es = SQRT(g') [grad(U) X grad(V)]: covariant radial unit vector (based on volume radial coordinate)"
     legend(8) = " BSUP{U,V} = B DOT GRAD{U,V}:  contravariant components of B"
     legend(9) = " JSUP{U,V} = SQRT(g) J DOT GRAD{U,V}"
     legend(10)= " K X B = Es DOT [K X B]: covariant component of K X B force"
     legend(11)= " K * B = K DOT B * SQRT(g')"
     legend(12)= " p' = dp/d(VOL): pressure gradient (based on volume radial coordinate)"
     legend(13)= " <KSUP{U,V}> = Int_u Int_v [KSUP{U,V}]/dV/ds"
  ENDIF

  nznt1 = nzeta*ntheta2
  ALLOCATE (avforce(ns),aminfor(ns),amaxfor(ns))
  ALLOCATE (bdotk(ns,nznt), bsubuv(ns,nznt),                        &
            bsubvu(ns,nznt), kperpu(nznt), kperpv(nznt),            &
            sqgb2(nznt), brhomn(0:mnyq,-nnyq:nnyq,0:1),kp2(nznt),   &
            jxb(nznt), jxb2(nznt), bsupu1(nznt),                    &
            bsubua(nznt1,0:1), bsubva(nznt1,0:1),                   &
            bsupv1(nznt), bsubu1(nznt), bsubv1(nznt),               &
            bsubsmn(ns,0:mnyq,-nnyq:nnyq,0:1),                      &
            bsubs_s(nznt),      bsubs_a(nznt),      sqrtg(nznt),    &
            bsubu_s(nznt1,0:1), bsubu_a(nznt1,0:1),                 &
            bsubv_s(nznt1,0:1), bsubv_a(nznt1,0:1), stat=j)
  IF (j .ne. 0) STOP 'Allocation error in jxbforce'

  ! NOTE: bsubuv, bsubvu are used to compute the radial current (should be zero)
  bsubsu = 0
  bsubsv = 0
  bsubuv = 0
  bsubvu = 0
  bsubsmn = 0

  radial: DO js = 1, ns

     ! Put bsubs on full mesh
     IF (js.gt.1 .and. js.lt.ns) THEN
        bsubs(js,:) = cp5*(bsubsh(js,:) + bsubsh(js+1,:))
     END IF

     ! undo radial scaling with sqrt(s) (done at iequi=1 case in bcovar)
     bsubu(js,:,1) = bsubu(js,:,1)/shalf(js)
     bsubv(js,:,1) = bsubv(js,:,1)/shalf(js)
     bsubua = 0
     bsubva = 0

     ! _s:     symmetric in u,v
     ! _a: antisymmetric in u,v on half (ntheta2) interval
     IF (lasym)  THEN
        bs1=>bsubs(js,:)
        bu1=>bsubu(js,:,:)
        bv1=>bsubv(js,:,:)
        CALL fsym_fft (bs1,     bu1,     bv1,     &
                       bsubs_s, bsubu_s, bsubv_s, & !     stellarator-symmetric parts
                       bsubs_a, bsubu_a, bsubv_a)   ! non-stellarator-symmetric parts
     ELSE
        ! only stellarator-symmetric parts
        bsubs_s(:) = bsubs(js,:)
        bsubu_s    = bsubu(js,:,:)
        bsubv_s    = bsubv(js,:,:)
     END IF

     ! FOURIER LOW-PASS FILTER bsubX
     DO m = 0, mpol1 ! 0, 1, ..., (mpol-1)
        mparity = MOD(m, 2)
        DO n = 0, ntor

           ! FOURIER TRANSFORM
           dnorm1 = one/r0scale**2
           IF (m .eq. mnyq)              dnorm1 = cp5*dnorm1
           IF (n .eq. nnyq .and. n.ne.0) dnorm1 = cp5*dnorm1

           bsubsmn1 = 0
           bsubsmn2 = 0
           bsubumn1 = 0
           bsubumn2 = 0
           bsubvmn1 = 0
           bsubvmn2 = 0
           IF (lasym) THEN
              ! SPH012314 in FixAray
              dnorm1 = 2*dnorm1

              bsubsmn3 = 0
              bsubsmn4 = 0
              bsubumn3 = 0
              bsubumn4 = 0
              bsubvmn3 = 0
              bsubvmn4 = 0
           END IF

           DO k = 1, nzeta
              lk = k
              DO j = 1, ntheta2
                 tsini1 = sinmui(j,m)*cosnv(k,n)*dnorm1 ! sin-cos
                 tsini2 = cosmui(j,m)*sinnv(k,n)*dnorm1 ! cos-sin
                 tcosi1 = cosmui(j,m)*cosnv(k,n)*dnorm1 ! cos-cos
                 tcosi2 = sinmui(j,m)*sinnv(k,n)*dnorm1 ! sin-sin

                 bsubsmn1 = bsubsmn1 + tsini1*bsubs_s(lk)          ! sin-cos
                 bsubsmn2 = bsubsmn2 + tsini2*bsubs_s(lk)          ! cos-sin
                 bsubvmn1 = bsubvmn1 + tcosi1*bsubv_s(lk, mparity) ! cos-cos
                 bsubvmn2 = bsubvmn2 + tcosi2*bsubv_s(lk, mparity) ! sin-sin
                 bsubumn1 = bsubumn1 + tcosi1*bsubu_s(lk, mparity) ! cos-cos
                 bsubumn2 = bsubumn2 + tcosi2*bsubu_s(lk, mparity) ! sin-sin

                 IF (lasym) THEN
                 bsubsmn3 = bsubsmn3 + tcosi1*bsubs_a(lk)          ! cos-cos
                 bsubsmn4 = bsubsmn4 + tcosi2*bsubs_a(lk)          ! sin-sin
                 bsubvmn3 = bsubvmn3 + tsini1*bsubv_a(lk, mparity) ! sin-cos
                 bsubvmn4 = bsubvmn4 + tsini2*bsubv_a(lk, mparity) ! cos-sin
                 bsubumn3 = bsubumn3 + tsini1*bsubu_a(lk, mparity) ! sin-cos
                 bsubumn4 = bsubumn4 + tsini2*bsubu_a(lk, mparity) ! cos-sin
                 END IF

                 lk = lk + nzeta

              END DO
           END DO

           ! FOURIER INVERSE TRANSFORM
           ! Compute on u-v grid (must add symmetric, antisymmetric parts for lasym=T)

           DO k = 1, nzeta
              lk = k
              DO j = 1, ntheta2
                 tcos1 = cosmu(j,m)*cosnv(k,n)
                 tcos2 = sinmu(j,m)*sinnv(k,n)
                 bsubua(lk,0) = bsubua(lk,0) + tcos1*bsubumn1 + tcos2*bsubumn2
                 bsubva(lk,0) = bsubva(lk,0) + tcos1*bsubvmn1 + tcos2*bsubvmn2

                 tcosm1 = cosmum(j,m)*cosnv(k,n)
                 tcosm2 = sinmum(j,m)*sinnv(k,n)
                 bsubsu(js,lk,0) = bsubsu(js,lk,0) + tcosm1*bsubsmn1 + tcosm2*bsubsmn2

                 tcosn1 = sinmu(j,m)*sinnvn(k,n)
                 tcosn2 = cosmu(j,m)*cosnvn(k,n)
                 bsubsv(js,lk,0) = bsubsv(js,lk,0) + tcosn1*bsubsmn1 + tcosn2*bsubsmn2

                 tsinm1 = sinmum(j,m)*cosnv(k,n)
                 tsinm2 = cosmum(j,m)*sinnv(k,n)
                 bsubvu(js,lk) = bsubvu(js,lk) + tsinm1*bsubvmn1 + tsinm2*bsubvmn2

                 tsinn1 = cosmu(j,m)*sinnvn(k,n)
                 tsinn2 = sinmu(j,m)*cosnvn(k,n)
                 bsubuv(js,lk) = bsubuv(js,lk) + tsinn1*bsubumn1 + tsinn2*bsubumn2

                 IF (lasym) THEN
                 tsin1 = sinmu(j,m)*cosnv(k,n)
                 tsin2 = cosmu(j,m)*sinnv(k,n)
                 bsubua(lk,1) = bsubua(lk,1) + tsin1*bsubumn3 + tsin2*bsubumn4
                 bsubva(lk,1) = bsubva(lk,1) + tsin1*bsubvmn3 + tsin2*bsubvmn4

                 bsubsu(js,lk,1) = bsubsu(js,lk,1) + tsinm1*bsubsmn3 + tsinm2*bsubsmn4
                 bsubsv(js,lk,1) = bsubsv(js,lk,1) + tsinn1*bsubsmn3 + tsinn2*bsubsmn4

                 bsubvu(js,lk) = bsubvu(js,lk) + tcosm1*bsubvmn3 + tcosm2*bsubvmn4
                 bsubuv(js,lk) = bsubuv(js,lk) + tcosn1*bsubumn3 + tcosn2*bsubumn4
                 END IF

                 lk = lk + nzeta

              END DO ! ntheta2
           END DO ! nzeta

           ! bsubsmn: coefficients of sin(mu)cos(nv), n>=0, cos(mu)sin(nv), n<0 (type=0)
           !                          cos(mu)cos(nv), n>=0, sin(mu)sin(nv), n<0 (type=1, nonzero only for lasym=T)!

           ! Don't need these except for comparison
           IF (lprint) then

             bsubsmn(js,m,n,0) = bsubsmn1
             IF (n .gt. 0) bsubsmn(js,m,-n,0) = bsubsmn2

             IF (.not.lasym) CYCLE

             bsubsmn(js,m,n,1) = bsubsmn3
             IF (n .gt. 0) bsubsmn(js,m,-n,0) = bsubsmn4

           end if ! lprint
        END DO
     END DO

     IF (lasym) THEN
        ! EXTEND FILTERED bsubu, bsubv TO NTHETA3 MESH
        ! NOTE: INDEX 0 - COS(mu-nv) SYMMETRY
        !             1 - SIN(mu-nv) SYMMETRY
        CALL fext_fft (bsubu(js,:,0), bsubua(:,0), bsubua(:,1))
        CALL fext_fft (bsubv(js,:,0), bsubva(:,0), bsubva(:,1))
     ELSE
        bsubu(js,:,0) = bsubua(:,0)
        bsubv(js,:,0) = bsubva(:,0)
     END IF

  END DO radial

  DEALLOCATE (bsubua, bsubva)

  ! EXTEND bsubsu, bsubsv TO NTHETA3 MESH
  IF (lasym) CALL fsym_invfft (bsubsu, bsubsv)


  if (open_dbg_context("jxbforce_bsub_lowpass", id=0)) then

    ! overwrite in-place; half-grid
    call add_real_3d("bsubu_e", ns, nzeta, ntheta3, bsubu(:,:,0))
    call add_real_3d("bsubv_e", ns, nzeta, ntheta3, bsubv(:,:,0))
    call add_real_3d("bsubu_o", ns, nzeta, ntheta3, bsubu(:,:,1))
    call add_real_3d("bsubv_o", ns, nzeta, ntheta3, bsubv(:,:,1))

    ! newly computed; full-grid
    call add_real_3d("bsubsu_e", ns, nzeta, ntheta3, bsubsu(:,:,0)) ! order=(/1, 4, 2, 3/))
    call add_real_3d("bsubsu_o", ns, nzeta, ntheta3, bsubsu(:,:,1)) ! order=(/1, 4, 2, 3/))
    call add_real_3d("bsubsv_e", ns, nzeta, ntheta3, bsubsv(:,:,0)) ! order=(/1, 4, 2, 3/))
    call add_real_3d("bsubsv_o", ns, nzeta, ntheta3, bsubsv(:,:,1)) ! order=(/1, 4, 2, 3/))

    ! newly computed; full-grid
    call add_real_3d("bsubuv", ns, nzeta, ntheta3, bsubuv)
    call add_real_3d("bsubvu", ns, nzeta, ntheta3, bsubvu)

    if (lprint) then
      call add_real_4d("bsubsmn", 2, ns, mpol, 2*ntor+1, bsubsmn, &
        order=(/ 4, 1, 2, 3 /))
    end if

    call close_dbg_out()
  end if






  ! SKIPS Bsubs Correction - uses Bsubs from metric elements
  IF (lbsubs) then

  ! Compute corrected Bsubs coefficients (brhomn) (impacts currents)
  ! by solving es dot (KXB - gradp_parallel) = 0 equation for brhomn in REAL SPACE
  ! Can be written Bsupu d(bs)/du + Bsupv d(bs)/dv = RHS (jxb below), bs==bsubs

  ! ANIMEC:
  ! brho==sigma B_s, pp1 and pp2 are the Jacobian times the hot particle parallel
  ! pressure radial gradient Amplitudes on the full integer mesh
  correct_bsubs: DO js = 2, ns-1

     jxb(:) = cp5*(gsqrt(js,:) + gsqrt(js+1,:)) ! re-use jxb array for Jacobian on full grid
     bsupu1(:) = cp5*(bsupu(js,:)*gsqrt(js,:) + bsupu(js+1,:)*gsqrt(js+1,:))
     bsupv1(:) = cp5*(bsupv(js,:)*gsqrt(js,:) + bsupv(js+1,:)*gsqrt(js+1,:))
     brho(js,:) = ohs* (   bsupu1(:)*(bsubu(js+1,:,0) - bsubu(js,:,0))    &
                         + bsupv1(:)*(bsubv(js+1,:,0) - bsubv(js,:,0))) &
                       + (pres(js+1) - pres(js))*ohs*jxb(:)

     ! SUBTRACT FLUX-SURFACE AVERAGE FORCE BALANCE FROM brho, OTHERWISE
     ! LOCAL FORCE BALANCE EQUATION B dot grad(Bs) = brho CAN'T BE SOLVED
     brho00(js) = SUM(brho(js,:)*wint(js:nrzt:ns))
     brho(js,:) = brho(js,:) - signgs*jxb(:)*brho00(js)/(cp5*(vp(js) + vp(js+1)))

     jxb(:) = brho(js,:)
     !              bsubsmn, frho, bsupu,  bsupv,  mmax, nmax, info
     CALL getbsubs (brhomn,  jxb,  bsupu1, bsupv1, mnyq, nnyq, info)
     IF (info .ne. 0) THEN
        PRINT *, 'Error in GETBRHO: info= ',info, ' js= ',js
     ELSE IF (lprint) THEN
        WRITE (33, *) ' JS = ', js
        IF (lasym) THEN
          WRITE (33, '(a)') '  M    N        BSUBS(old)        BSUBS(new)        BSUBS(old)        BSUBS(new)'
        ELSE
          WRITE (33, '(a)') '  M    N        BSUBS(old)        BSUBS(new)'
        END IF
        DO m = 0, mpol1
           DO n = -ntor, ntor
              IF (lasym) THEN
                WRITE(33,1223) m, n, bsubsmn(js,m,n,0), brhomn(m,n,0), bsubsmn(js,m,n,1), brhomn(m,n,1)
              ELSE
                WRITE(33,1224) m, n, bsubsmn(js,m,n,0), brhomn(m,n,0)
              END IF
           END DO
        END DO
     END IF
 1223    FORMAT (i4,1x,i4,4(6x,1p,e12.3))
 1224    FORMAT (i4,1x,i4,2(6x,1p,e12.3))

     ! Recompute bsubsu,v now using corrected bsubs
     ! Store old values (itheta,izeta) for checking force balance later
     itheta(js,:) = bsubsu(js,:,0)
     izeta (js,:) = bsubsv(js,:,0)

     IF (info .ne. 0) CYCLE ! skip until next iteration of radial loop

     bsubsu(js,:,:) = 0
     bsubsv(js,:,:) = 0
     bsubs_s = 0
     IF (lasym) bsubs_a = 0

     DO m = 0, mnyq
        DO n = 0, nnyq
           IF (n .eq. 0) THEN
              bsubsmn1 = brhomn(m,0,0)
              bsubsmn2 = 0
           ELSE
              bsubsmn1 = brhomn(m, n,0)
              bsubsmn2 = brhomn(m,-n,0)
           END IF

           IF (lasym) THEN
              IF (n .eq. 0) THEN
                 bsubsmn3 = brhomn(m,0,1)
                 bsubsmn4 = 0
              ELSE
                 bsubsmn3 = brhomn(m, n,1)
                 bsubsmn4 = brhomn(m,-n,1)
              END IF
           END IF

           DO k = 1, nzeta
              lk = k
              DO j = 1, ntheta2

                 tsin1 = sinmu(j,m)*cosnv(k,n)
                 tsin2 = cosmu(j,m)*sinnv(k,n)
                 bsubs_s(lk) = bsubs_s(lk) + tsin1*bsubsmn1 + tsin2*bsubsmn2

                 tcosm1 = cosmum(j,m)*cosnv(k,n)
                 tcosm2 = sinmum(j,m)*sinnv(k,n)
                 bsubsu(js,lk,0) = bsubsu(js,lk,0) + tcosm1*bsubsmn1 + tcosm2*bsubsmn2

                 tcosn1 = sinmu(j,m)*sinnvn(k,n)
                 tcosn2 = cosmu(j,m)*cosnvn(k,n)
                 bsubsv(js,lk,0) = bsubsv(js,lk,0) + tcosn1*bsubsmn1 + tcosn2*bsubsmn2

                 IF (lasym) THEN
                 tcos1 = cosmu(j,m)*cosnv(k,n)
                 tcos2 = sinmu(j,m)*sinnv(k,n)
                 bsubs_a(lk) = bsubs_a(lk) + tcos1*bsubsmn3 + tcos2*bsubsmn4

                 tsinm1 = sinmum(j,m)*cosnv(k,n)
                 tsinm2 = cosmum(j,m)*sinnv(k,n)
                 bsubsu(js,lk,1) = bsubsu(js,lk,1) + tsinm1*bsubsmn3 + tsinm2*bsubsmn4

                 tsinn1 = cosmu(j,m)*sinnvn(k,n)
                 tsinn2 = sinmu(j,m)*cosnvn(k,n)
                 bsubsv(js,lk,1) = bsubsv(js,lk,1) + tsinn1*bsubsmn3 + tsinn2*bsubsmn4
                 END IF

                 lk = lk + nzeta

              END DO
           END DO
        END DO
     END DO

     IF (lasym) THEN
        ! EXTEND TO FULL (ntheta3) u-GRID
        bs1 => bsubs(js,:)
        CALL fext_fft (bs1, bsubs_a, bsubs_s) ! TODO: a and s swapped ???
     ELSE
        bsubs(js,:) = bsubs_s(:)
     END IF

  END DO correct_bsubs ! radial loop: correct bsubs on every surface individually

  ! EXTEND bsubsu, bsubsv TO NTHETA3 MESH
  IF (lasym) CALL fsym_invfft (bsubsu, bsubsv)

  ! CHECK FORCE BALANCE: SQRT(g)*(bsupu*bsubsu + bsupv*bsubsv) = brho
  IF (lprint) then
    WRITE (33, '(/,2a,/)') 'ANGLE INDEX       B*grad(Bs)      Frhs          Fold'
    check_fb: DO js = 2, ns-1
       bsupu1(:) = cp5*(bsupu(js,:)*gsqrt(js,:) + bsupu(js+1,:)*gsqrt(js+1,:))
       bsupv1(:) = cp5*(bsupv(js,:)*gsqrt(js,:) + bsupv(js+1,:)*gsqrt(js+1,:))
       kp2(:) = bsupu1(:)*bsubsu(js,:,0) + bsupv1(:)*bsubsv(js,:,0)
       jxb(:) = bsupu1(:)*itheta(js,:)   + bsupv1(:)*izeta(js,:)

       WRITE (33, '(/,a,i4)') 'JS = ',js
       DO lk = 1, nznt
          WRITE(33,1230) lk, brho(js,lk),  kp2(lk),  jxb(lk)
   1230 FORMAT (i9,5x, 1p,3e14.4)
       END DO
    END DO check_fb
  end if ! force balance printout

  end if ! lbsubs

  DEALLOCATE (bsubs_s, bsubs_a, bsubu_s, bsubu_a, bsubv_s, bsubv_a, stat=lk)

  ! Compute end point values for bsubs
  bsubs(1,:)  = 2*bsubs(2,:)  - bsubs(3,:)
  !bsubs(ns,:) = 2*bsubs(ns,:) - bsubs(ns-1,:) ! TODO: from ns, ns-1 to ns ???
  bsubs(ns,:) = 2*bsubs(ns-1,:) - bsubs(ns-2,:)












  ! Now compute currents on the FULL radial mesh
  ! Here:
  !
  ! Itheta = sqrt(g) * Ksupu
  ! Izeta  = sqrt(g) * Ksupv
  ! Ksupx  = K dot grad(x)                          x=(u,v)
  ! jxb    = (K X B) dot (grad-u X grad-v) sqrt(g)
  ! bdotk  = sigma*sqrt(g)*K dot B
  ! kperpx = (B X gradp) dot grad(x) / |B|**2       x=(u,v)
  ! sqgb2  = sigma*sqrt(g)*|B|**2
  ! sqrtg  = sqrt(g)
  ! pprime = d(p||)/dV
  !
  ! kp2   == |k-perp|**2 = kperpu**2 * guu + 2*kperpu*kperpv*guv + kperpv**2 * gvv
  ! This was compared to the alternative expression (agreed very well):
  ! |j-perp|**2 = |grad-s|**2 * (dp/ds)**2 / |B|**2
  !
  ! Note: Multiply currents, pressure by 1/mu0 to get in mks units!
  !       TWOPI*TWOPI factor incorporated in vp (thru ovp factor below), so V' = (2pi)**2*vp
  !
  ALLOCATE(                                                         &
       bsubs3(ns,nzeta,ntheta3), bsubv3(ns,nzeta,ntheta3),          &
       bsubu3(ns,nzeta,ntheta3), jxb_gradp(ns,nzeta,ntheta3),       &
       jcrossb(ns,nzeta,ntheta3), bsupv3(ns,nzeta,ntheta3),         &
       bsupu3(ns,nzeta,ntheta3), jsups3(ns,nzeta,ntheta3),          &
       jsupv3(ns,nzeta,ntheta3), jsupu3(ns,nzeta,ntheta3),          &
       jdotb_sqrtg(ns,nzeta,ntheta3), sqrtg3(ns,nzeta,ntheta3),     &
       phin(ns), toroidal_angle(nzeta), stat=j)

  ! delete incoming leftovers
  itheta = 0.0_dp
  izeta  = 0.0_dp

  bsubs3      = 0
  bsubv3      = 0
  bsubu3      = 0
  jxb_gradp   = 0
  jcrossb     = 0
  bsupv3      = 0
  bsupu3      = 0
  jsups3      = 0
  jsupv3      = 0
  jsupu3      = 0
  phin        = 0
  phin(ns)    = 1
  jdotb_sqrtg = 0
  sqrtg3      = 0

  bdotk  = 0

  ALLOCATE (pprime(nznt), pprim(ns),stat=j)
  pprim = 0

  avforce = 0
  aminfor = 0
  amaxfor = 0

  dnorm1 = twopi*twopi

  DO js = 2, ns1

     ovp = c2p0/(vp(js+1) + vp(js))/dnorm1

     tjnorm = ovp*signgs

     sqgb2(:nznt) =   gsqrt(js+1,:nznt) * (bsq(js+1,:nznt)-pres(js+1)) &
                    + gsqrt(js  ,:nznt) * (bsq(js  ,:nznt)-pres(js  ))

     ! TAKE THIS OUT: MAY BE POORLY CONVERGED AT THIS POINT....
     ! IF (ANY(sqgb2(:nznt)*signgs .le. zero)) &
     !   STOP ' SQGB2 <= 0 in JXBFORCE'

     ! dp/ds here
     pprime(:) = ohs*(pres(js+1)-pres(js))/mu0

     kperpu(:nznt) = cp5*(bsubv(js+1,:nznt,0) + bsubv(js,:nznt,0))*pprime(:)/sqgb2
     kperpv(:nznt) =-cp5*(bsubu(js+1,:nznt,0) + bsubu(js,:nznt,0))*pprime(:)/sqgb2

     kp2(:nznt)=cp5*(    kperpu**2        * (guu(js+1:nrzt:ns) + guu(js:nrzt:ns)) &
                     + 2*kperpu*kperpv    * (guv(js+1:nrzt:ns) + guv(js:nrzt:ns)) &
                     +          kperpv**2 * (gvv(js+1:nrzt:ns) + gvv(js:nrzt:ns)))

     itheta(js,:nznt) =  bsubsv(js,:nznt,0) - ohs*(bsubv(js+1,:nznt,0) - bsubv(js,:nznt,0))
     izeta(js,:nznt)  = -bsubsu(js,:nznt,0) + ohs*(bsubu(js+1,:nznt,0) - bsubu(js,:nznt,0))

     itheta(js,:nznt) = itheta(js,:nznt)/mu0
     izeta(js,:nznt)  = izeta(js,:nznt)/mu0

     ! can be computed above (before lbsubs, where this appears as well)
     sqrtg(:) = cp5*(gsqrt(js,:) + gsqrt(js+1,:))

     bsupu1(:nznt) = cp5*(bsupu(js+1,:nznt)*gsqrt(js+1,:) + bsupu(js,:nznt)*gsqrt(js,:)) / sqrtg(:)
     bsupv1(:nznt) = cp5*(bsupv(js+1,:nznt)*gsqrt(js+1,:) + bsupv(js,:nznt)*gsqrt(js,:)) / sqrtg(:)

     bsubu1(:nznt) = cp5*(bsubu(js+1,:nznt,0) + bsubu(js,:nznt,0))
     bsubv1(:nznt) = cp5*(bsubv(js+1,:nznt,0) + bsubv(js,:nznt,0))

     jxb(:nznt) = ovp*(itheta(js,:nznt) * bsupv1(:nznt) - izeta (js,:nznt) * bsupu1(:nznt))

     bdotk(js,:nznt) = itheta(js,:nznt) * bsubu1(:nznt) + izeta (js,:nznt) * bsubv1(:nznt)

     pprime(:nznt) = ovp*pprime(:nznt)
     pnorm = one/(ABS(pprime(1)) + EPSILON(pprime(1)))

     amaxfor(js) = MAXVAL(jxb(:nznt) - pprime(:))*pnorm
     aminfor(js) = MINVAL(jxb(:nznt) - pprime(:))*pnorm

     amaxfor(js) = 100*MIN(amaxfor(js), 9.999_dp)
     aminfor(js) = 100*MAX(aminfor(js),-9.999_dp)

     avforce(js) = SUM(wint(2:nrzt:ns)*(jxb(:nznt) - pprime(:)))

     pprim(js) = SUM(wint(js:nrzt:ns)*pprime(:))

     ! Compute <K dot B>, <B sup v> = signgs*phip
     ! jpar2 = <j||**2>, jperp2 = <j-perp**2>,  with <...> = flux surface average

     jdotb(js) = dnorm1*tjnorm*SUM(bdotk(js,:nznt)*wint(2:nrzt:ns))
     bdotb(js) = dnorm1*tjnorm*SUM(sqgb2(:nznt)   *wint(2:nrzt:ns))

     bdotgradv(js) = cp5*dnorm1*tjnorm*(phip(js) + phip(js+1)) ! TODO: could also use phips here?

     jpar2(js) = dnorm1*tjnorm*SUM(bdotk(js,:nznt)**2 * wint(2:nrzt:ns)/sqgb2(:nznt))
     jperp2(js)= dnorm1*tjnorm*SUM(kp2(:nznt)*wint(2:nrzt:ns)*sqrtg(:nznt))

     IF (lprint_flag) THEN
        phin(js) = phi(js)/phi(ns)
        DO lz = 1, nzeta
           toroidal_angle(lz)=REAL(360*(lz-1),rprec)/nzeta ! TODO: unused!
           DO lt = 1, ntheta3
              lk = lz + nzeta*(lt-1)
              ! lu (js,lz,lt ) =  lt

              jsupu3 (js,lz,lt) = ovp*itheta(js,lk)
              jsupv3 (js,lz,lt) = ovp*izeta(js,lk)
              jsups3 (js,lz,lt) = ovp*(bsubuv(js,lk) - bsubvu(js,lk))/mu0

              bsupu3 (js,lz,lt) = bsupu1(lk)
              bsupv3 (js,lz,lt) = bsupv1(lk)

              jcrossb (js,lz,lt) = jxb(lk)
              jxb_gradp (js,lz,lt) = (jxb(lk) - pprime(lk))
              jdotb_sqrtg (js,lz,lt) = ovp*bdotk(js,lk)

              sqrtg3(js,lz,lt) = sqrtg(lk)*ovp

              bsubu3(js,lz,lt) = bsubu(js,lk,0)
              bsubv3(js,lz,lt) = bsubv(js,lk,0)
              bsubs3(js,lz,lt) = bsubs(js,lk)
           END DO
        END DO
     ENDIF ! lprint_flag

  END DO ! radial

  ! Need in wrout
  izeta( 1,:nznt) = c2p0*izeta(   2,:nznt) - izeta(   3,:nznt)
  izeta(ns,:nznt) = c2p0*izeta(ns-1,:nznt) - izeta(ns-2,:nznt)

  jdotb(1)  = c2p0*jdotb(   2) - jdotb(   3)
  jdotb(ns) = c2p0*jdotb(ns-1) - jdotb(ns-2)

  !bdotb(1)  = c2p0*bdotb(   3) - bdotb(   2) ! TODO: 2 <--> 3 ???
  bdotb(1)  = c2p0*bdotb(   2) - bdotb(   3)
  bdotb(ns) = c2p0*bdotb(ns-1) - bdotb(ns-2)

  bdotgradv(1)  = c2p0*bdotgradv(   2) - bdotgradv(   3)
  bdotgradv(ns) = c2p0*bdotgradv(ns-1) - bdotgradv(ns-2)

  jpar2(1)   = 0
  jpar2(ns)  = 0

  jperp2(1)  = 0
  jperp2(ns) = 0

  !pprim( 1) = 2*pprim(ns-1) - pprim(ns-2) ! TODO: what is going on here ??
  pprim( 1) = 2*pprim(   2) - pprim(   3)
  pprim(ns) = 2*pprim(ns-1) - pprim(ns-2)


  if (open_dbg_context("jxbout", id=0)) then

    call add_real_3d("itheta",  ns, nzeta, ntheta3, itheta)
    call add_real_3d("izeta",   ns, nzeta, ntheta3, izeta)
    call add_real_3d("bdotk",   ns, nzeta, ntheta3, bdotk)

    call add_real_1d("amaxfor",   ns, amaxfor)
    call add_real_1d("aminfor",   ns, aminfor)
    call add_real_1d("avforce",   ns, avforce)
    call add_real_1d("pprim",     ns, pprim)
    call add_real_1d("jdotb",     ns, jdotb)
    call add_real_1d("bdotb",     ns, bdotb)
    call add_real_1d("bdotgradv", ns, bdotgradv)
    call add_real_1d("jpar2",     ns, jpar2)
    call add_real_1d("jperp2",    ns, jperp2)

    call add_real_3d("jsupu3",      ns, nzeta, ntheta3, jsupu3)
    call add_real_3d("jsupv3",      ns, nzeta, ntheta3, jsupv3)
    call add_real_3d("jsups3",      ns, nzeta, ntheta3, jsups3)
    call add_real_3d("bsupu3",      ns, nzeta, ntheta3, bsupu3)
    call add_real_3d("bsupv3",      ns, nzeta, ntheta3, bsupv3)
    call add_real_3d("jcrossb",     ns, nzeta, ntheta3, jcrossb)
    call add_real_3d("jxb_gradp",   ns, nzeta, ntheta3, jxb_gradp)
    call add_real_3d("jdotb_sqrtg", ns, nzeta, ntheta3, jdotb_sqrtg)
    call add_real_3d("sqrtg3",      ns, nzeta, ntheta3, sqrtg3)
    call add_real_3d("bsubu3",      ns, nzeta, ntheta3, bsubu3)
    call add_real_3d("bsubv3",      ns, nzeta, ntheta3, bsubv3)
    call add_real_3d("bsubs3",      ns, nzeta, ntheta3, bsubs3)

    call close_dbg_out()
  end if





  IF (lprint_flag) THEN
     ! declare variables in netCDF file
     CALL cdf_define(njxbout, vn_legend, legend)
     CALL cdf_define(njxbout, vn_mpol, mpol)
     CALL cdf_define(njxbout, vn_ntor, ntor)
     CALL cdf_define(njxbout, vn_phin, phin)
     CALL cdf_define(njxbout, vn_radial_surfaces, ns)
     CALL cdf_define(njxbout, vn_poloidal_grid_points, ntheta3)
     CALL cdf_define(njxbout, vn_toroidal_grid_points, nzeta)
     CALL cdf_define(njxbout, vn_avforce, avforce)
     CALL cdf_define(njxbout, vn_jdotb, jdotb)

     CALL cdf_define(njxbout, vn_sqg_bdotk, jdotb_sqrtg)
     CALL cdf_define(njxbout, vn_sqrtg, sqrtg3)

     CALL cdf_define(njxbout, vn_bdotgradv, bdotgradv)
     CALL cdf_define(njxbout, vn_pprime, pprim)
     CALL cdf_define(njxbout, vn_aminfor, aminfor)
     CALL cdf_define(njxbout, vn_amaxfor, amaxfor)
     CALL cdf_define(njxbout, vn_jsupu, jsupu3)
     CALL cdf_define(njxbout, vn_jsupv, jsupv3)
     CALL cdf_define(njxbout, vn_jsups, jsups3)
     CALL cdf_define(njxbout, vn_bsupu, bsupu3)
     CALL cdf_define(njxbout, vn_bsupv, bsupv3)
     CALL cdf_define(njxbout, vn_jcrossb, jcrossb)
     CALL cdf_define(njxbout, vn_jxb_gradp, jxb_gradp)
     CALL cdf_define(njxbout, vn_bsubu, bsubu3)
     CALL cdf_define(njxbout, vn_bsubv, bsubv3)
     CALL cdf_define(njxbout, vn_bsubs, bsubs3)

     ! actually write data
     CALL cdf_write(njxbout, vn_legend,               legend     )
     CALL cdf_write(njxbout, vn_mpol,                 mpol       )
     CALL cdf_write(njxbout, vn_ntor,                 ntor       )
     CALL cdf_write(njxbout, vn_phin,                 phin       )
     CALL cdf_write(njxbout, vn_radial_surfaces,      ns         )
     CALL cdf_write(njxbout, vn_poloidal_grid_points, ntheta3    )
     CALL cdf_write(njxbout, vn_toroidal_grid_points, nzeta      )
     CALL cdf_write(njxbout, vn_avforce,              avforce    )
     CALL cdf_write(njxbout, vn_jdotb,                jdotb      )

     CALL cdf_write(njxbout, vn_sqg_bdotk,            jdotb_sqrtg)
     CALL cdf_write(njxbout, vn_sqrtg,                sqrtg3     )

     CALL cdf_write(njxbout, vn_bdotgradv,            bdotgradv  )
     CALL cdf_write(njxbout, vn_pprime,               pprim      )
     CALL cdf_write(njxbout, vn_aminfor,              aminfor    )
     CALL cdf_write(njxbout, vn_amaxfor,              amaxfor    )
     CALL cdf_write(njxbout, vn_jsupu,                jsupu3     )
     CALL cdf_write(njxbout, vn_jsupv,                jsupv3     )
     CALL cdf_write(njxbout, vn_jsups,                jsups3     )
     CALL cdf_write(njxbout, vn_bsupu,                bsupu3     )
     CALL cdf_write(njxbout, vn_bsupv,                bsupv3     )
     CALL cdf_write(njxbout, vn_jcrossb,              jcrossb    )
     CALL cdf_write(njxbout, vn_jxb_gradp,            jxb_gradp  )
     CALL cdf_write(njxbout, vn_bsubu,                bsubu3     )
     CALL cdf_write(njxbout, vn_bsubv,                bsubv3     )
     CALL cdf_write(njxbout, vn_bsubs,                bsubs3     )

     CALL cdf_close(njxbout)

     DEALLOCATE(                                                       &
          bsubs3, bsubv3, bsubu3, jxb_gradp, jcrossb, bsupv3,          &
          bsupu3, jsups3, jsupv3, jsupu3, jdotb_sqrtg, phin,           &
          toroidal_angle, sqrtg3, stat=j)
  END IF ! lprint_flag

  DEALLOCATE (kperpu, kperpv, sqgb2, sqrtg, kp2, brhomn, bsubsmn,   &
      jxb, jxb2, bsupu1, bsupv1, bsubu1, bsubv1, avforce, aminfor,  &
      amaxfor, pprim, stat=j)

  ! COMPUTE MERCIER CRITERION
  bdotk = mu0*bdotk

  !            gsqrt, bsq, bdotj, iotas, wint
  CALL Mercier(gsqrt, bsq, bdotk, iotas, wint, &
  !            r1, rt, rz, zt, zz
               r1, ru, rv, zu, zv, &
  !            bsubu, vp, phips, pres, ns, nznt
               bsubu, vp, phips, pres, ns, nznt)

  DEALLOCATE (bdotk, bsubuv, bsubvu, pprime, stat=j)

END SUBROUTINE jxbforce
