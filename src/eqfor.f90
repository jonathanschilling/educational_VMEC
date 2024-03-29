!> \file
!> \brief Basis physics analysis and evaluaton of force balance.
!>        This is where most of the contents of the \c threed1 output file is computed.

!> \brief Basis physics analysis and evaluaton of force balance.
!>        This is where most of the contents of the \c threed1 output file is computed.
!>
!> @param br cylindrical component of magnetic field \f$B^R\f$
!> @param bz cylindrical component of magnetic field \f$B^Z\f$
!> @param bsubu covariant component of magnetic field \f$B_\theta\f$
!> @param bsubv covariant component of magnetic field \f$B_\zeta\f$
!> @param tau Jacobian \f$\sqrt{g} = R \tau\f$
!> @param rzl_array state vector (all Fourier coefficients) of VMEC
!> @param ier_flag error flag
SUBROUTINE eqfor(br, bz, bsubu, bsubv, tau, rzl_array, ier_flag)

  USE vmec_main
  USE vmec_params
  USE realspace
  USE vforces, r12 => armn_o, bsupu => crmn_e, bsupv => czmn_e,         &
               gsqrt => azmn_o, bsq => bzmn_o, izeta => azmn_e,         &
               brho => bzmn_e, bphi => czmn_o, curtheta => brmn_e
  USE vacmod, only: bphiv, bsqvac, bsubvvac
  USE vmec_io
  USE mgrid_mod
  USE stel_constants, ONLY: pi

  use dbgout

  IMPLICIT NONE

  INTEGER :: ier_flag
  REAL(rprec), DIMENSION(ns,nznt,0:1), INTENT(in) :: bsubu, bsubv
  REAL(rprec), DIMENSION(nrzt), INTENT(out) :: br, bz
  REAL(rprec), DIMENSION(nrzt), INTENT(out) :: tau
  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax), TARGET, INTENT(in) :: rzl_array

  INTEGER :: i, icount, itheta, js, l, loff,                       &
     lpi, lt, n, n1, noff,                          &
     iv, iu, lk, nplanes
  REAL(rprec), DIMENSION(:), POINTER :: rmags, zmags, rmaga, zmaga
  REAL(rprec), DIMENSION(:,:,:), POINTER :: rmncc,zmnsc
  REAL(rprec), DIMENSION(ns) :: phi1, chi1, jPS2
  REAL(rprec) :: modb(nznt)
  REAL(rprec), DIMENSION(:), ALLOCATABLE ::                             &
     btor_vac, btor1, dbtor, phat, t12u, guu_1u, surf_area,             &
     r3v, redge, rbps1u, bpol2vac, phipf_loc
  REAL(rprec) :: aminr1, aminr2in, anorm,                       &
     aspectratio, betai, betstr, scaling_ratio,                         &
     bminz2, bminz2in, musubi,                           &
      cur0,                                 &
     delphid_exact, delta1, delta2, delta3, lambda,             &
     er, es, fac, facnorm, factor, fgeo,        &
     flao, fpsi0, pavg,                      &
     rcen, rcenin, rgeo,                            &
     rlao, rshaf, rshaf1, rshaf2, s11, s12,       &
     s13, s2, s3, sigr0, sigr1, sigz1, smaleli,                         &
     sumbpol, sumbtot, sumbtor, sump,                            &
     sump2, sump20, t1, tz, jpar_perp=0, jparPs_perp=0,                 &
     toroidal_flux, vnorm, xmax,                    &
     xmida, xmidb, xmin, rzmax, rzmin, zxmax, zxmin,    &
     zmax, zmin, yr1u, yz1u, waist(2), height(2)
  REAL(rprec) d_of_kappa, loc_jpar_perp, loc_jparPs_perp
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rax_symm, zax_symm, rax_asym, zax_asym


  ! POINTER ASSOCIATIONS
  rmags => rzl_array(1,:,0,rcc)
  zmags => rzl_array(1,:,0,zcs+ntmax)
  rmncc => rzl_array(:,:,:,rcc)
  zmnsc => rzl_array(:,:,:,zsc+ntmax)
  IF (lasym) THEN
    rmaga => rzl_array(1,:,0,rcs)
    zmaga => rzl_array(1,:,0,zcc+ntmax)
  END IF

  ! crmn_o == bss on half grid
  !         r12  rs    zs    ru12  zu12  bsubs   bsupu  bsupv  br  bphi  bz
  CALL bss (r12, bzmn, brmn, azmn, armn, crmn_o, bsupu, bsupv, br, bphi, bz)

  ! STORE EDGE VALUES OF B-FIELD
  IF (lfreeb .and. ivac.gt.1) THEN
     IF (ALLOCATED(bredge)) DEALLOCATE (bredge, bpedge, bzedge)
     ALLOCATE (bredge(2*nznt), bpedge(2*nznt), bzedge(2*nznt), stat=i)
     IF (i .ne. 0) STOP 'Error in EQFOR allocating bredge'
     DO iv = 1,nzeta
        DO iu = 1,ntheta3
           lk = iv + nzeta*(iu-1)
           n1 = ns*lk ! indices radially on LCFS
           bredge(lk) = 1.5_dp*br(n1)   - cp5*br(n1-1)
           bpedge(lk) = 1.5_dp*bphi(n1) - cp5*bphi(n1-1)
           bzedge(lk) = 1.5_dp*bz(n1)   - cp5*bz(n1-1)
        END DO
     END DO
  END IF

  ! NOTE: JXBFORCE ROUTINE MUST BE CALLED TO COMPUTE IZETA, JDOTB
  !       ON OUTPUT, J, IZETA, JDOTB ARE IN MKS UNITS (1/MU0 FACTOR)
  !
  ! CAUTION: THIS CALL WILL WRITE OVER br, bz
  !              bsupu, bsupv, bsubu, bsubv, bsubsh, bsubsu, bsubsv
  CALL jxbforce (bsupu, bsupv, bsubu, bsubv, crmn_o, rcon,   zcon, &
                 gsqrt, bsq, curtheta, izeta, brho, ier_flag)
  !              gsqrt, bsq, itheta,   izeta, brho, ier_flag








  ! HALF-MESH VOLUME-AVERAGED BETA
  tau(1) = 0
  tau(2:nrzt) = signgs*wint(2:nrzt)*gsqrt(2:nrzt) ! re-use tau for ready-to-integrate Jacobian
  DO i = 2, ns
     s2 = SUM( bsq(i:nrzt:ns)*tau(i:nrzt:ns))/vp(i) - pres(i)
     beta_vol(i) = pres(i)/s2

     ! hijack loop to also compute overr ???
     overr(i) = SUM(tau(i:nrzt:ns)/r12(i:nrzt:ns)) / vp(i)
  END DO

  ! extrapolate surface-averaged beta profile to magnetic axis
  betaxis = c1p5*beta_vol(2) - cp5*beta_vol(3)

  WRITE (nthreed, 5)
5 FORMAT(/,' NOTE:  S=normalized toroidal flux (0 - 1)',/,              &
           '        U=poloidal angle (0 - 2*pi)',/,                     &
           '        V=geometric toroidal angle (0 - 2*pi)',/,           &
           '       <RADIAL FORCE> = d(Ipol)/dPHI',                      &
           ' - IOTA*d(Itor)/dPHI - dp/dPHI * d(VOL)/dPHI',/,            &
           '                      = d(VOL)/dPHI*[<JSUPU>',              &
           ' - IOTA*<JSUPV> - SIGN(JAC)*dp/dPHI]',/,                    &
           '       (NORMED TO SUM OF INDIVIDUAL TERMS)',//,             &
           '      S     <RADIAL    TOROIDAL      IOTA     ',            &
           ' <JSUPU>    <JSUPV>     d(VOL)/',                           &
           '   d(PRES)/    <M>     PRESF    <BSUBU>    <BSUBV>',        &
           '      <J.B>      <B.B>',/,                                  &
           '             FORCE>      FLUX                  ',           &
           '                       d(PHI) ',                            &
           '    d(PHI)                             ',/,148('-'),/)

  ALLOCATE (phipf_loc(ns))

  ! interpolate pressure and phip onto full grid
  presf(1)     =               c1p5*pres(2) - cp5*pres(3)
  phipf_loc(1) = twopi*signgs*(c1p5*phip(2) - cp5*phip(3))
  DO i = 2,ns1
     presf(i)     = cp5      *       (pres(i) + pres(i+1))
     phipf_loc(i) = cp5*twopi*signgs*(phip(i) + phip(i+1))
  END DO
  presf(ns)     =               c1p5*pres(ns) - cp5*pres(ns1)
  phipf_loc(ns) = twopi*signgs*(c1p5*phip(ns) - cp5*phip(ns1))

  ! integrate flux differentials to get flux profiles
  phi1(1) = zero
  chi1(1) = zero
  DO i = 2, ns
     phi1(i) = phi1(i-1) + hs*phip(i)
     chi1(i) = chi1(i-1) + hs*(phip(i)*iotas(i))
  END DO
  chi = twopi*chi1




  CALL calc_fbal(bsubu, bsubv)

  ! NOTE:  jcuru, jcurv on FULL radial mesh coming out of calc_fbal
  !        They are local (surface-averaged) current densities (NOT integrated in s)
  !        jcurX = (dV/ds)/twopi**2 <JsupX>   for X=u,v
  bucof(1) = 0
  bvcof(1) = c1p5*bvco(2) - cp5*bvco(3)
  DO i = 2,ns1
     equif(i) = equif(i)*vpphi(i) /         &
                (  ABS(jcurv(i)*chipf(i))   &
                 + ABS(jcuru(i)*phipf(i))   &
                 + ABS(presgrad(i)*vpphi(i)))
     bucof(i) = cp5*(buco(i) + buco(i+1))
     bvcof(i) = cp5*(bvco(i) + bvco(i+1))
  END DO
  bucof(ns) = c1p5*buco(ns) - cp5*buco(ns1)
  bvcof(ns) = c1p5*bvco(ns) - cp5*bvco(ns1)


  ! extrapolate full-grid quantites to axis and LCFS
  equif(1)  = c2p0*equif(2)   - equif(3)
  equif(ns) = c2p0*equif(ns1) - equif(ns1-1)

  jcurv(1)  = c2p0*jcurv(2)   - jcurv(3)
  jcurv(ns) = c2p0*jcurv(ns1) - jcurv(ns1-1)

  jcuru(1)  = c2p0*jcuru(2)   - jcuru(3)
  jcuru(ns) = c2p0*jcuru(ns1) - jcuru(ns1-1)

  presgrad(1)  = c2p0*presgrad(2)   - presgrad(3)
  presgrad(ns) = c2p0*presgrad(ns1) - presgrad(ns1-1)

  vpphi(1)  = c2p0*vpphi(2)   - vpphi(3)
  vpphi(ns) = c2p0*vpphi(ns1) - vpphi(ns1-1)




  ! NOTE: phipf = phipf_loc/(twopi), phipf_loc ACTUAL (twopi factor) Toroidal flux derivative
  ! SPH/JDH (060211): remove twopi factors from <JSUPU,V> (agree with output in JXBOUT file)
  fac = twopi*signgs
  DO js = 1, ns
     es = (js - 1)*hs           ! full-grid s value
     cur0 = fac*vpphi(js)*twopi ! == dV/ds = dV/dPHI * d(PHI/ds)  (V=actual volume)
     WRITE (nthreed, 30)               &
       es,                             & ! S
       equif(js),                      & ! <RADIAL FORCE>
       fac*phi1(js),                   & ! TOROIDAL FLUX
       iotaf(js),                      & ! IOTA
       jcuru(js)/vpphi(js)/mu0,        & ! <JSUPU>
       jcurv(js)/vpphi(js)/mu0,        & ! <JSUPV>
       cur0/phipf_loc(js),             & ! d(VOL)/d(PHI)
       presgrad(js)/phipf_loc(js)/mu0, & ! d(PRES)/d(PHI)
       specw(js),                      & ! <M>
       presf(js)/mu0,                  & ! PRESF
       bucof(js),                      & ! <BSUBU>
       bvcof(js),                      & ! <BSUBV>
       jdotb(js),                      & ! <J.B>
       bdotb(js)                         ! <B.B>
  END DO
30 FORMAT(1p,2e10.2,2e12.4,4e11.3,0p,f7.3,1p,5e11.3)

  if (open_dbg_context("threed1_firstTable", id=0)) then

    call add_real_1d("beta_vol", ns-1, beta_vol(2:ns))
    call add_real_1d("overr",    ns-1, overr(2:ns))

    call add_real("betaxis", betaxis)

    call add_real_1d("presf",     ns, presf)
    call add_real_1d("phipf_loc", ns, phipf_loc)
    call add_real_1d("phi1", ns, phi1)
    call add_real_1d("chi1", ns, chi1)
    call add_real_1d("chi",  ns, chi)

    call add_real_1d("iotaf",  ns, iotaf)

    call add_real_1d("specw",  ns, specw)

    call add_real_1d("equif",  ns, equif)

    call add_real_1d("bucof",  ns, bucof)
    call add_real_1d("bvcof",  ns, bvcof)

    call add_real_1d("jcurv",  ns, jcurv)
    call add_real_1d("jcuru",  ns, jcuru)

    call add_real_1d("presgrad", ns, presgrad)
    call add_real_1d("vpphi",    ns, vpphi)

    call add_real_1d("jdotb",    ns, jdotb)
    call add_real_1d("bdotb",    ns, bdotb)

    call close_dbg_out()
  end if

  DEALLOCATE (phipf_loc)








  ! MAKE SURE WOUT FILE DOES NOT REQUIRE ANY STUFF COMPUTED BELOW....

  ! Calculate mean (toroidally averaged) poloidal cross section area & toroidal flux
  anorm = twopi*hs
  vnorm = twopi*anorm
  toroidal_flux = anorm * SUM(bsupv(2:nrzt)*tau(2:nrzt))

  ! Calculate poloidal circumference and normal surface area and aspect ratio
  ! Normal is | dr/du X dr/dv | = SQRT [R**2 guu + (RuZv - RvZu)**2]
  ALLOCATE (guu_1u(nznt), surf_area(nznt))
  guu_1u(:nznt) = ru0(ns:nrzt:ns)*ru0(ns:nrzt:ns) + zu0(ns:nrzt:ns)*zu0(ns:nrzt:ns)

  surf_area(:nznt) = wint(ns:nrzt:ns)*SQRT(guu_1u(:nznt)) ! re-use surf_area to compute circumference
  circum_p = twopi*SUM(surf_area(:nznt))

  surf_area(:nznt) = wint(ns:nrzt:ns)*SQRT(                               &
       + (r1(ns:nrzt:ns,0) + r1(ns:nrzt:ns,1))**2 * guu_1u(:nznt)         &
       +((rv(ns:nrzt:ns,0) + rv(ns:nrzt:ns,1))*zu0(ns:nrzt:ns)            &
       - (zv(ns:nrzt:ns,0) + zv(ns:nrzt:ns,1))*ru0(ns:nrzt:ns))**2 )
  surf_area_p = twopi**2*SUM(surf_area(:nznt))
  DEALLOCATE (guu_1u)

  aspect = aspectratio()

  ! Also, estimate mean elongation of plasma from the following relations
  ! for an axisymmetric torus with elliptical cross section and semi-axes
  ! a and a * kappa (kappa >= 1)
  !
  ! surf_area _p = 2*pi*R * 2*pi*a ctwiddle(kappa_p)
  ! volume_p    = 2*pi*R * pi*a ** 2 * kappa_p
  ! cross_area _p =   pi*a ** 2 * kappa_p
  !
  ! The cirumference of an ellipse of semi-axes a and a * kappa_p is
  !    2 * pi * a ctwiddle(kappa_p)
  ! The exact form for ctwiddle is 4 E(1 - kappa_p^2) / (2 pi), where
  !  E is the complete elliptic integral of the second kind
  ! (with parameter argument m, not modulus argument k)
  !
  ! The coding below implements an approximate inverse of the function
  ! d(kappa) = ctwiddle(kappa) / sqrt(kappa)
  ! The approximate inverse is
  !    kappa = 1 + (pi^2/8) * (d^2+sqrt(d^4-1)-1)
  ! Note that the variable aminor_p, for an elliptic cross section,
  ! would be a * sqrt(kappa)
  d_of_kappa = surf_area_p * aminor_p / ( 2 * volume_p)
  kappa_p = 1 + (pi * pi / 8) * (d_of_kappa ** 2 + SQRT(ABS(d_of_kappa ** 4 - 1)) -1)

  ! TODO: is this used anywhere or directly overwritten below?
  aminr1 = 2*volume_p/surf_area_p

  ! OUTPUT BETAS, INDUCTANCES, SAFETY FACTORS, ETC.
  ! (EXTRACTED FROM FQ-CODE, 9-10-92)
  !
  ! b poloidals (cylindrical estimates)
  rcen = cp5*(router + rinner)               !geometric center
  n = 0
  n1 = n + 1
  rcenin = DOT_PRODUCT(rmncc(ns,n1,:mpol1+1:2), mscale(:mpol1:2)*nscale(n)) ! TODO: used anywhere?

  l = (mpol1+1)/2
  ALLOCATE (t12u(l))
  t12u(:l) = mscale(1:mpol1:2)*nscale(n)
  aminr2in = DOT_PRODUCT(rmncc(ns,n1,2:mpol1+1:2),t12u(:l)) ! TODO: used anywhere?
  bminz2in = DOT_PRODUCT(zmnsc(ns,n1,2:mpol1+1:2),t12u(:l)) ! TODO: used anywhere?
  bminz2   = DOT_PRODUCT(zmnsc(ns,n1,2:mpol1+1:2),t12u(:l)) ! TODO: used anywhere?
  DEALLOCATE (t12u)

  ! vol av minor radius
  aminr1 = SQRT(c2p0*volume_p/(twopi*twopi*r00))

  ! cylindrical estimates for beta poloidal
  sump = vnorm*SUM(vp(2:ns)*pres(2:ns))
  pavg = sump/volume_p
  factor = 2*pavg

  ! delphid_exact = Integral[ (Bvac - B) * dSphi ]
  ! rshaf [= RT in Eq.(12), Phys Fluids B 5 (1993) 3119]
  !
  ! Note: tau = |gsqrt|*wint
  ALLOCATE (btor_vac(nznt), btor1(nznt), dbtor(nznt), phat(nznt), redge(nznt))
  ! Eq. 20 in Shafranov
  delphid_exact = zero
  musubi = zero
  rshaf1 = zero
  rshaf2 = zero
  DO js = 2, ns
     btor_vac(:nznt) = rbtor/r12(js:nrzt:ns) ! TODO: assumes B_tor ~ 1/R ???
     delphid_exact = delphid_exact + SUM( (btor_vac(:nznt)/r12(js:nrzt:ns) - bsupv(js:nrzt:ns))*tau(js:nrzt:ns) )

     btor1(:nznt) = r12(js:nrzt:ns)*bsupv(js:nrzt:ns)
     dbtor(:nznt) = btor1(:nznt)**2 - btor_vac(:nznt)**2
     musubi = musubi - SUM(dbtor(:nznt)*tau(js:nrzt:ns))

     ! phat is temporarily re-used for something else...
     phat(:nznt) = bsq(js:nrzt:ns) - cp5*btor_vac(:nznt)**2
     phat(:nznt) = (phat(:nznt) - dbtor(:nznt))*tau(js:nrzt:ns)
     rshaf1 = rshaf1 + SUM(phat(:nznt))
     rshaf2 = rshaf2 + SUM(phat(:nznt)/r12(js:nrzt:ns))
  END DO
  delphid_exact = anorm*delphid_exact
  rshaf = rshaf1/rshaf2
  fpsi0 = c1p5*bvco(2) - cp5*bvco(3)
  b0 = fpsi0/r00

  redge(:nznt) = r1(ns:nrzt:ns,0) + r1(ns:nrzt:ns,1)
  IF (lfreeb .and. ivac.gt.1) THEN
     phat = bsqvac - cp5*(bsubvvac/redge)**2
  ELSE
     phat = c1p5*bsq(ns:nrzt:ns) - cp5*bsq(ns-1:nrzt:ns) &
            - cp5*(rbtor/redge(:))**2
  END IF

  DEALLOCATE (btor_vac, btor1, dbtor)

  rmax_surf = MAXVAL(r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
  rmin_surf = MINVAL(r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
  zmax_surf = MAXVAL(ABS(z1(ns:nrzt:ns,0)+z1(ns:nrzt:ns,1)))

  DO js = 2, ns
     modb(:nznt) = SQRT(c2p0*(bsq(js:nrzt:ns) - pres(js)))
     CALL bextrema (modb, bmin(1,js), bmax(1,js), nzeta, ntheta2)
  END DO

  WRITE (nthreed, 75) bmin(1,ns), bmax(1,ns), bmin(ntheta2,ns), bmax(ntheta2,ns)
75 FORMAT(/' Magnetic field modulation (averaged over toroidal angle)',/,        &
          1x,71('-')/,' BMIN(u=0)             = ',f14.6/                         &
          ' BMAX(u=0)             = ',f14.6/' BMIN(u=pi)            = ',         &
          f14.6/' BMAX(u=pi)            = ',f14.6/)

  ! output geometrical, |B| quantities
  CALL elongation (r1, z1, waist, height)

  sumbtot = 2*(vnorm*SUM(bsq(2:nrzt)*tau(2:nrzt)) - sump)
  sumbtor = vnorm*SUM(tau(2:nrzt)*(r12(2:nrzt)*bsupv(2:nrzt))**2)
  sumbpol = sumbtot - sumbtor
  betapol = 2*sump/sumbpol
  sump20 = 2*sump
  sump2 = SUM(pres(2:ns)*pres(2:ns)*vp(2:ns)*vnorm)
  betatot = sump20/sumbtot
  betator = sump20/sumbtor
  VolAvgB = SQRT(ABS(sumbtot/volume_p))
  IonLarmor = 0.0032_dp/VolAvgB ! TODO: which ion is assumed here ?

  jPS2(2:ns1) = jpar2(2:ns1) - jdotb(2:ns1)**2/bdotb(2:ns1)
  jpar_perp   = SUM( jpar2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
  jparPS_perp = SUM(  jPS2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
  s2          = SUM(jperp2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
  IF (s2 .ne. zero) THEN
     jpar_perp   = jpar_perp  /s2
     jparPS_perp = jparPS_perp/s2
  END IF

  IF (ntor .gt. 1) THEN ! lthreed ?? ntor>1 means ntor .ge. 2 --> ntor=1 is already non-axisymmetric!
  WRITE (nthreed, 80) aspect, kappa_p, volume_p, cross_area_p,          &
     surf_area_p, circum_p, Rmajor_p, Aminor_p, rmin_surf,              &
     rmax_surf, zmax_surf, waist(1), height(1), waist(2), height(2)
  ELSE ! axisymmetric
  WRITE (nthreed, 80) aspect, kappa_p, volume_p, cross_area_p,          &
     surf_area_p, circum_p, Rmajor_p, Aminor_p, rmin_surf,              &
     rmax_surf, zmax_surf, waist(1), height(1)
  END IF
80 FORMAT(/,' Geometric and Magnetic Quantities',/,1x,71('-')/,         &
            ' Aspect Ratio          = ',f14.6, /                        &
            ' Mean Elongation       = ',f14.6, /                        &
            ' Plasma Volume         = ',f14.6,' [M**3]',/               &
            ' Cross Sectional Area  = ',f14.6,' [M**2]',/               &
            ' Normal Surface Area   = ',f14.6,' [M**2]',/               &
            ' Poloidal Circumference= ',f14.6,' [M]',/                  &
            ' Major Radius          = ',f14.6,' [M]',                   &
            ' (from Volume and Cross Section)',/                        &
            ' Minor Radius          = ',f14.6,' [M]',                   &
            ' (from Cross Section)',/                                   &
            ' Minimum (inboard)  R  = ',f14.6,' [M]',/                  &
            ' Maximum (outboard) R  = ',f14.6,' [M]',/                  &
            ' Maximum height     Z  = ',f14.6,' [M]',/                  &
            ' Waist (v = 0)   in R  = ',f14.6,' [M]',/                  &
            ' Full Height(v = 0)    = ',f14.6,' [M]',:,/                &
            ' Waist (v = pi)  in R  = ',f14.6,' [M]',:,/                &
            ' Full Height(v = pi)   = ',f14.6,' [M]')
  ! TODO: : means optional in above format?

  WRITE (nthreed, 85) toroidal_flux, 1.e-6_dp/mu0 * ctor, rbtor,        &
         rbtor0, VolAvgB, IonLarmor, jpar_perp, jparPS_perp
85 FORMAT(' Toroidal Flux         = ',f14.6,' [Wb]',/                   &
          ' Toroidal Current      = ',f14.6,' [MA]',/                   &
          ' RBtor(s=1)            = ',f14.6,' [T-m]',/                  &
          ' RBtor(s=0)            = ',f14.6,' [T-m]',/                  &
          ' Volume Average B      = ',f14.6,' [T]',/                    &
          ' Ion Larmor Radius     = ',f14.6,' [M] X Ti(keV)**0.5',/     &
          ' <J||**2>/<J-perp**2>  = ',f14.6,' (Vol. Averaged)',/        &
          ' <JPS**2>/<J-perp**2>  = ',f14.6,' (Vol. Averaged)',/ )

  WRITE (nthreed, 90)
90 FORMAT(/,71('-'),/,&
            ' MORE GEOMETRIC AND PHYSICS QUANTITIES',/,71('-'),/,&
            ' Toroidal Plane: Phi = 0',/,                     &
            5x,'j',3x,'psi-psiaxis',9x,'a [M]',3x,'ellipticity',3x,'indentation',&
            7x,'d-shape',4x,'rel. shift',6x,'<J||**2>/',4x,'<JPS**2>/',/,   &
            95x,                          '<J-perp**2>',3x,'<J-perp**2>'/, &
            ' -----',8(2x,12('-')))

  fac = twopi*hs*signgs
  psi(1) = zero
  ALLOCATE (r3v(ns-1))
  r3v(:ns-1) = fac*phip(2:ns)*iotas(2:ns)
  DO i = 1, ns - 1
     psi(1+i) = psi(i) + r3v(i)
  END DO
  DEALLOCATE (r3v)

  ! initialize to known conents for dbgout of axisymmetric run
  ygeo   = 0.0_dp
  yinden = 0.0_dp
  yellip = 0.0_dp
  ytrian = 0.0_dp
  yshift = 0.0_dp

  ! nphi-plane, noff = 1,....,nzeta
  PLANES: DO nplanes = 1, 2

     IF (nplanes .eq. 1) THEN
        ! nphi=0
        noff = 1
     ELSE
        IF (nzeta .eq. 1) EXIT
        WRITE (nthreed, 95)
        ! nphi=180
        noff = 1+nzeta/2
     END IF

     ygeo(nplanes, 1) = zero
     DO js = 2, ns

        zmin =  HUGE(zmin)
        zmax = -HUGE(zmax)
        xmin =  HUGE(xmin)
        xmax = -HUGE(xmax)

        rzmax = zero

        ! Theta = 0 to pi in upper half of X-Z plane
        DO icount = 1,2 ! TODO: why second loop over toroidal offset ?
           ! nphi-plane, n1 = noff,...,nzeta
           n1 = noff
           IF (icount .eq. 2) then
              ! (twopi-v), reflected plane
              n1 = MOD(nzeta-(noff-1),nzeta)+1
           end if
           loff = js + ns*(n1-1)

           t1 = one
           IF (icount .eq. 2) t1 = -one

           DO itheta = 1,ntheta2
              yr1u = r1(loff,0) + sqrts(js)*r1(loff,1)
              yz1u = z1(loff,0) + sqrts(js)*z1(loff,1)
              yz1u = t1*yz1u

              IF (yz1u .ge. zmax) THEN
                 zmax = ABS(yz1u)
                 rzmax = yr1u
              ELSEIF (yz1u .le. zmin) THEN
                 zmin = yz1u
                 rzmin = yr1u
              END IF

              IF (yr1u .ge. xmax) THEN
                 xmax = yr1u
                 zxmax = yz1u
              ELSEIF (yr1u .le. xmin) THEN
                 xmin = yr1u
                 zxmin = yz1u
              END IF

              loff = loff + ns*nzeta
           END DO ! itheta
        END DO ! icount

        ! ----------

        ! theta=180
        lpi = ns*((noff-1) + nzeta*(ntheta2-1))

        ! theta=0
        lt  = ns*(noff-1)

        xmida = r1(js+lpi,0) + sqrts(js)*r1(js+lpi,1)
        xmidb = r1(js+lt,0)  + sqrts(js)*r1(js+lt,1)

        ! ----------

        ! Geometric major radius
        rgeo = cp5*(xmidb + xmida)

        ! Geometric minor radius
        ygeo(nplanes, js) = cp5*(xmidb - xmida)

        ! Geometric indentation
        yinden(nplanes, js) = (xmida - xmin)/(xmax - xmin)

        ! Geometric ellipticity
        yellip(nplanes, js) = ( zmax - zmin)/(xmax - xmin)

        ! Geometric triangularity
        ytrian(nplanes, js) = (rgeo - rzmax)/(xmax - xmin)

        ! Geometric shift measured from magnetic axis
        yshift(nplanes, js) = (r1(1+lt,0)-rgeo)/(xmax - xmin)

        ! --- start independent stuff
        IF (jperp2(js) .eq. zero) &
            jperp2(js) = EPSILON(jperp2(js))

        loc_jpar_perp = jpar2(js)/jperp2(js)
        IF (js .lt. ns) THEN
           loc_jparPS_perp = jPS2(js)/jperp2(js)
        ELSE
           loc_jparPS_perp = zero
        END IF
        ! --- end independent stuff

        IF (nplanes .eq. 1) THEN
           WRITE (nthreed, 120) js, psi(js),            &
              ygeo(nplanes, js), yellip(nplanes, js),      &
              yinden(nplanes, js), ytrian(nplanes, js), &
              yshift(nplanes, js), loc_jpar_perp,       &
              loc_jparPS_perp
        ELSE
           WRITE (nthreed, 120) js, psi(js), &
              ygeo(nplanes, js), yellip(nplanes, js),   &
              yinden(nplanes, js), ytrian(nplanes, js), &
              yshift(nplanes, js)
        END IF
     END DO ! js
  END DO PLANES

 95 FORMAT(/,71('-'),/,' Toroidal Plane: Phi = 180/Nfp',/,71('-'),/)
120 FORMAT(1x,i5,6f14.5,1p,3e14.2)

  if (open_dbg_context("threed1_geomag", id=0)) then

    call add_real("anorm", anorm)
    call add_real("vnorm", vnorm)
    call add_real("toroidal_flux", toroidal_flux)

    call add_real_2d("surf_area", nzeta, ntheta3, surf_area)

    call add_real("circum_p", circum_p)
    call add_real("surf_area_p", surf_area_p)

    call add_real("volume_p",     volume_p)
    call add_real("cross_area_p", cross_area_p)
    call add_real("Rmajor_p",     Rmajor_p)
    call add_real("Aminor_p",     Aminor_p)
    call add_real("aspect",       aspect)

    call add_real("kappa_p", kappa_p)
    call add_real("rcen", rcen)
    call add_real("rcenin", rcenin)
    call add_real("aminr2in", aminr2in)
    call add_real("bminz2in", bminz2in)
    call add_real("bminz2", bminz2)
    call add_real("aminr1", aminr1)
    call add_real("sump", sump)
    call add_real("pavg", pavg)

    call add_real("delphid_exact", delphid_exact)
    call add_real("musubi", musubi)
    call add_real("rshaf1", rshaf1)
    call add_real("rshaf2", rshaf2)
    call add_real("rshaf", rshaf)
    call add_real("fpsi0", fpsi0)
    call add_real("b0", b0)

    call add_real_2d("redge", nzeta, ntheta3, redge)
    call add_real_2d("phat",  nzeta, ntheta3, phat)

    call add_real("rmax_surf", rmax_surf)
    call add_real("rmin_surf", rmin_surf)
    call add_real("zmax_surf", zmax_surf)

    call add_real("bmin_1_ns", bmin(1,ns))
    call add_real("bmax_1_ns", bmax(1,ns))
    call add_real("bmin_ntheta2_ns", bmin(ntheta2,ns))
    call add_real("bmax_ntheta2_ns", bmax(ntheta2,ns))

    IF (ntor .gt. 1) THEN
        call add_real_1d("waist",  2, waist)
        call add_real_1d("height", 2, height)
    else
        call add_real_1d("waist",  1, waist(1))
        call add_real_1d("height", 1, height(1))
    end if

    call add_real("sumbtot", sumbtot)
    call add_real("sumbtor", sumbtor)
    call add_real("sumbpol", sumbpol)
    call add_real("betapol", betapol)
    call add_real("sump20", sump20)
    call add_real("sump2", sump2)
    call add_real("betatot", betatot)
    call add_real("betator", betator)
    call add_real("VolAvgB", VolAvgB)
    call add_real("IonLarmor", IonLarmor)

    call add_real_1d("jPS2", ns-2, jPS2(2:ns1))

    call add_real("s2", s2)
    call add_real("jpar_perp", jpar_perp)
    call add_real("jparPS_perp", jparPS_perp)

    call add_real_1d("psi", ns, psi)

    call add_real_2d("ygeo",   2, ns,   ygeo)
    call add_real_2d("yinden", 2, ns-1, yinden(:,2:ns))
    call add_real_2d("yellip", 2, ns-1, yellip(:,2:ns))
    call add_real_2d("ytrian", 2, ns-1, ytrian(:,2:ns))
    call add_real_2d("yshift", 2, ns-1, yshift(:,2:ns))

    call close_dbg_out()
  end if

  WRITE (nthreed, 130)
130 FORMAT(//,' Magnetic Fields and Pressure',/,1x,71('-'))
  fac = cp5/mu0
  WRITE (nthreed, 140)                    &
     sump/mu0,      pavg/mu0,             &
     fac*sumbpol,   fac*sumbpol/volume_p, &
     fac*sumbtor,   fac*sumbtor/volume_p, &
     fac*sumbtot,   fac*sumbtot/volume_p, &
     c1p5*sump/mu0, c1p5*pavg/mu0
140 FORMAT(' Volume Integrals (Joules) and Volume Averages (Pascals)',/,&
           24x,'Integral',6x,'Average',/,        &
           ' pressure         = ',1p,2e14.6,/,&
           ' bpol**2 /(2 mu0) = ',2e14.6,/, &
           ' btor**2/(2 mu0)  = ',2e14.6,/, &
           ' b**2/(2 mu0)     = ',2e14.6,/,&
           ' EKIN (3/2p)      = ',2e14.6,/)

  if (open_dbg_context("threed1_volquant", id=0)) then

    call add_real("int_p", sump/mu0)
    call add_real("avg_p", pavg/mu0)

    call add_real("int_bpol", fac*sumbpol)
    call add_real("avg_bpol", fac*sumbpol/volume_p)

    call add_real("int_btor", fac*sumbtor)
    call add_real("avg_btor", fac*sumbtor/volume_p)

    call add_real("int_modb", fac*sumbtot)
    call add_real("avg_modb", fac*sumbtot/volume_p)

    call add_real("int_ekin", c1p5*sump/mu0)
    call add_real("avg_ekin", c1p5*pavg/mu0)

    call close_dbg_out()
  end if

  WRITE (nthreed, 800)
800 FORMAT(/,' MAGNETIC AXIS COEFFICIENTS'/,                            &
             '    n     rmag       zmag        rmag        zmag',/,     &
             '        (cos nv)   (sin nv)    (sin nv)    (cos nv)',/)
  loff = LBOUND(rmags,1)

  allocate(rax_symm(0:ntor), zax_symm(0:ntor), rax_asym(0:ntor), zax_asym(0:ntor))

  DO n = 0, ntor
     n1 = n + loff

     t1 = mscale(0)*nscale(n)

     tz = t1
     IF (.not.lthreed) tz = 0

     rax_symm(n) =  t1*rmags(n1)
     zax_symm(n) = -tz*zmags(n1)
     IF (lasym) THEN
       rax_asym(n) = -tz*rmaga(n1)
       zax_asym(n) =  t1*zmaga(n1)
     else
       rax_asym(n) = 0.0_dp
       zax_asym(n) = 0.0_dp
     end if

     IF (lasym) THEN
        WRITE (nthreed, 820) n, t1*rmags(n1), (-tz*zmags(n1)), -tz*rmaga(n1),   t1*zmaga(n1)
     ELSE
        WRITE (nthreed, 820) n, t1*rmags(n1), (-tz*zmags(n1))
     END IF
  END DO
820 FORMAT(i5,1p,4e12.4)

  if (open_dbg_context("threed1_axis", id=0)) then

    call add_real_1d("rax_symm", ntor+1, rax_symm)
    call add_real_1d("zax_symm", ntor+1, zax_symm)
    call add_real_1d("rax_asym", ntor+1, rax_asym)
    call add_real_1d("zax_asym", ntor+1, zax_asym)

    call close_dbg_out()
  end if

  betstr = c2p0*SQRT(sump2/volume_p)/(sumbtot/volume_p)

  WRITE (nthreed, 150) betatot, betapol, betator
150 FORMAT(/,' From volume averages over plasma, betas are',/,          &
             ' beta total    = ',f14.6,/, &
             ' beta poloidal = ',f14.6,/, &
             ' beta toroidal = ',f14.6,/  )

  WRITE (nthreed, 160) rbtor, betaxis, betstr ! TODO: should this maybe be bsubvvac ?
160 FORMAT(' R * Btor-vac         = ',f14.6,' [Wb/M]',/,                &
           ' Peak Beta            = ',f14.6,/,                          &
           ' Beta-star            = ',f14.6,/)

  if (open_dbg_context("threed1_beta", id=0)) then

    call add_real("betatot", betatot)
    call add_real("betapol", betapol)
    call add_real("betator", betator)

    call add_real("rbtor",   rbtor)
    call add_real("betaxis", betaxis)
    call add_real("betstr",  betstr)

    call close_dbg_out()
  end if

  ! Shafranov surface integrals s1,s2
  ! Plasma Physics vol 13, pp 757-762 (1971)
  ! Also, s3 = .5*S3, defined in Lao, Nucl. Fusion 25, p.1421 (1985)
  ! Note: if ctor = 0, use Int(Bsupu*Bsubu dV) for ctor*ctor/R
  ! Phys. Fluids B, Vol 5 (1993) p 3121, Eq. 9a-9d
  ALLOCATE (rbps1u(nznt), bpol2vac(nznt))
  IF (lfreeb .and. ivac.gt.1) THEN
     bpol2vac = 2*bsqvac - bphiv*bphiv
  ELSE
     bpol2vac = 2*(c1p5*bsq(ns:nrzt:ns)   - cp5*bsq(ns-1:nrzt:ns))       &
              -  ((c1p5*bsupv(ns:nrzt:ns) - cp5*bsupv(ns-1:nrzt:ns))     &
              * redge)**2
  END IF

  ! Compute current-like norm (factor) in Eq.(8), <a> * int(Bpol**2 * dA)
  ! where <a> == 2*pi*Rs in Eq. 8 is the effective minor radius = Vol/Asurf
  ! (corrects wrong description of Rs in paper, which is NOT the major radius)
  ! This aminr1 = 1/2 the "correct" aminr1
  aminr1 = volume_p/surf_area_p
  factor = twopi**2 * aminr1 * SUM(bpol2vac*surf_area)
  factor = one/factor
  facnorm = factor * twopi**2

  ! Lao's definition of normalization factor
  scaling_ratio = (mu0*curtor/circum_p)**2 * volume_p
  scaling_ratio = scaling_ratio*factor

  rbps1u(:nznt) = facnorm*redge(:nznt)*phat(:nznt)*wint(ns:nznt*ns:ns) ! TODO: why not use nrzt instead of nznt*ns ?
  sigr0 = SUM(rbps1u(:nznt)*zu0(ns:nrzt:ns))
  sigr1 = SUM(rbps1u(:nznt)*zu0(ns:nrzt:ns)*redge(:nznt))
  sigz1 =-SUM(rbps1u(:nznt)*ru0(ns:nrzt:ns)*(z1(ns:nrzt:ns,0) + z1(ns:nrzt:ns,1)))
  DEALLOCATE (redge, phat, rbps1u, bpol2vac, surf_area)

  er = sigr1 + sigz1

  ! LAO, NUCL.FUS. 25 (1985) 1421
  rlao = volume_p/(twopi*cross_area_p)
  flao = rshaf/rlao
  fgeo = rshaf/rcen

  smaleli = factor*sumbpol
  betai   = 2*factor*sump
  musubi  = vnorm*factor*musubi
  lambda = cp5*smaleli + betai

  ! Shafranov def. based on RT, Eq.(12)
  s11 = er - rshaf*sigr0

  ! R = Rgeometric
  s12 = er - rcen*sigr0

  ! R = RLao
  s13 = er - rlao*sigr0
  s2  = sigr0*rshaf

  ! 1/2 S3 in Eq.(14c)
  s3  = sigz1

  delta1 = zero
  delta2 = one - fgeo
  delta3 = one - flao

  WRITE (nthreed, 168)
168 FORMAT(' Shafranov Surface Integrals',/                             &
           ' Ref: S. P. Hirshman, Phys. Fluids B, 5, (1993) 3119',/,    &
           ' Note: s1 = S1/2, s2 = S2/2, where ',                       &
           ' s1,s2 are the Shafranov definitions,',/,                   &
           ' and s3 = S3/2, where S3 is Lao''s definition.',/,          &
           ' The quantity lsubi gives the ratio of volume poloidal',    &
         /,' field energy to the field energy estimated from the',      &
         /,' surface integral in Eq.8.',/,1x,22('-'),/)

  WRITE (nthreed, 170) rshaf, rcen, rlao, scaling_ratio,                &
     s3, smaleli, musubi, betai, lambda
170 FORMAT(' RT (Pressure-weighted)  = ',f14.6,' [M]',/,                &
           ' RG (Geometric)          = ',f14.6,' [M]',/,                &
           ' RL (Vol/2*pi*Area-Lao)  = ',f14.6,' [M]',/,                &
           ' Poloidal Field Energy',/,                                  &
           ' Normalization Ratio     = ',f14.6,' (Lao/Hirshman)',//,    &
           ' s3                      = ',f14.6,/,                       &
           ' lsubi                   = ',f14.6,/,                       &
           ' musubi                  = ',f14.6,/,                       &
           ' betai                   = ',f14.6,/,                       &
           ' lambda                  = ',f14.6,/)

  WRITE (nthreed, 174) &
     delta1, delta2, delta3,                             & ! delta = 1 - RT/R
     s11, s12, s13,                                      & ! s1
     s2, s2/fgeo, s2/flao,                               & ! s2
     musubi + s11, musubi + s12, musubi + s13,           & ! betai (Mui + s1)
     cp5*s11 + s2, cp5*s12 + s2/fgeo, cp5*s13 + s2/flao, & ! lambda (s1/2 + s2)
     cp5*(3*betai+smaleli-musubi)/(s11+s2     ) - one,   & ! 1st Shafr''v relation
     cp5*(3*betai+smaleli-musubi)/(s12+s2/fgeo) - one,   & ! (3*Betai + Li - Mui)/[2*(s1+s2)] - 1
     cp5*(3*betai+smaleli-musubi)/(s13+s2/flao) - one,   &
     cp5     *(betai+smaleli+musubi)/s2 - one,           & ! Radial force balance
     cp5*fgeo*(betai+smaleli+musubi)/s2 - one,           & !(Betai + Li + Mui)/(2*s2) - 1
     cp5*flao*(betai+smaleli+musubi)/s2 - one
174 FORMAT(/,32x,'R = RT',12x,'R = RG',12x,'R = RL',/,                  &
           20x,3(10x,8('-')),/,                                         &
           ' delta = 1 - RT/R     = ',3(f14.6,4x),/,                    &
           ' s1                   = ',3(f14.6,4x),/,                    &
           ' s2                   = ',3(f14.6,4x),/,                    &
           ' betai (Mui + s1)     = ',3(f14.6,4x),/,                    &
           ' lambda (s1/2 + s2)   = ',3(f14.6,4x),/,                    &
           ' 1st Shafr''v relation = ',3(f14.6,4x),/,                   &
           ' (3*Betai + Li - Mui)/[2*(s1+s2)] - 1',/,                   &
           ' Radial force balance = ',3(f14.6,4x),/,                    &
           ' (Betai + Li + Mui)/(2*s2) - 1',/)


  if (open_dbg_context("threed1_shafrint", id=0)) then

    call add_real("scaling_ratio", scaling_ratio)

    call add_real("rlao", rlao)
    call add_real("flao", flao)
    call add_real("fgeo", fgeo)

    call add_real("smaleli", smaleli)
    call add_real("betai",   betai)
    call add_real("musubi",  musubi)
    call add_real("lambda",  lambda)

    call add_real("s11", s11)
    call add_real("s12", s12)
    call add_real("s13", s13)
    call add_real("s2", s2)
    call add_real("s3", s3)

    call add_real("delta1", delta1)
    call add_real("delta2", delta2)
    call add_real("delta3", delta3)

    call close_dbg_out()
  end if

END SUBROUTINE eqfor
