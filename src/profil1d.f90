!> \file
!> \brief Compute phip and iota profiles on full grid.

!> \brief Compute phip and iota profiles on full grid.
!>
!> @param xc state vector of VMEC, i.e., all Fourier coefficients of \f$R\f$, \f$Z\f$ and \f$\lambda\f$
!> @param xcdot velocity vector in Fourier space
!> @param lreset xc will be zeroes if this is true
SUBROUTINE profil1d()
  USE vmec_main
  USE vmec_params, ONLY: signgs, lamscale, rcc, pdamp
  USE realspace, ONLY: shalf, sqrts

  use dbgout

  IMPLICIT NONE

  INTEGER :: i
  REAL(rprec) :: Itor, si, tf, pedge, vpnorm, torflux_edge, polflux_edge

  REAL(rprec), EXTERNAL :: pcurr, pmass, piota, &
                           torflux, torflux_deriv, polflux, polflux_deriv

  ! COMPUTE PHIP, IOTA PROFILES ON FULL-GRID
  ! COMPUTE MASS PROFILE ON HALF-GRID BY READING INPUT COEFFICIENTS.
  ! PRESSURE CONVERTED TO INTERNAL UNITS BY MULTIPLICATION BY mu0 = 4*pi*10**-7

  torflux_edge = signgs * phiedge / twopi
  si = torflux(one)
  IF (si .ne. zero) torflux_edge = torflux_edge/si

  polflux_edge = torflux_edge
  si = polflux(one)
  IF (si .ne. zero) polflux_edge = polflux_edge/si

  ! z00 gets assiged in reset_params and in funct3d
  r00 = rmn_bdy(0,0,rcc)

  ! zero fluxes and current at magnetic axis
  phips(1) = 0
  chips(1) = 0
  icurv(1) = 0

  ! half-grid quantities: s_i are shifted inwards by 0.5 grid points
  DO i = 2,ns
     si = hs*(i-c1p5)
     tf = MIN(one, torflux(si))
     phips(i) = torflux_edge * torflux_deriv(si)
     chips(i) = torflux_edge * polflux_deriv(si)
     iotas(i) = piota(tf) ! evaluate iota profile
     icurv(i) = pcurr(tf) ! evaluate current profile
  END DO
  phips(ns+1) = 2*phips(ns)-phips(ns-1) ! virtual point outside the LCFS

  ! Compute lamscale factor for "normalizing" lambda (needed for scaling hessian)
  lamscale = SQRT(hs*SUM(phips(2:ns)**2))
  IF (lamscale .EQ. 0) STOP 'PHIP == 0: ERROR!'

  IF (lflip) THEN
     iotas = -iotas
     chips = -chips
  END IF

  ! full-grid quantities
  DO i = 1,ns
     si = hs*(i-1)
     tf = MIN(one, torflux(si))
     phipf(i) = torflux_edge * torflux_deriv(si)
     chipf(i) = torflux_edge * polflux_deriv(si)
     iotaf(i) = piota(tf)
  ENDDO

  ! SCALE CURRENT TO MATCH INPUT EDGE VALUE, CURTOR
  ! FACTOR OF SIGNGS NEEDED HERE, SINCE MATCH IS MADE TO LINE
  ! INTEGRAL OF BSUBU (IN GETIOTA) ~ SIGNGS * CURTOR
  pedge = pcurr(one) ! pcurr gives current enclosed between axis and given parameter
  ! --> pedge is temporarily re-used to store the radial integral over current density profile (still in Amperes)
  Itor = 0
  ! TODO: what is this? might have something to do with preventing division-by-zero
  ! when computing scaled toroidal current profile from pcurr and curtor...
  IF (ABS(pedge) .gt. ABS(EPSILON(pedge)*curtor)) then
     Itor = signgs*currv/(twopi*pedge)
  end if
  icurv(2:ns) = Itor*icurv(2:ns)

  ! POSSIBLE PRESSURE PEDESTAL FOR S >= SPRES_PED
  spres_ped = ABS(spres_ped)
  DO i = 2,ns
    si = hs*(i - c1p5)
    ! NORMALIZE mass so dV/dPHI (or dV/dPSI) in pressure to mass relation
    ! See line 195 of bcovar: pres(2:ns) = mass(2:ns)/vp(2:ns)**gamma
    tf = MIN(one, torflux(si))
    vpnorm = torflux_edge * torflux_deriv(si)

    ! should either check if tf .gt. tf_pres_ped
    !          or   eval pmass at min(1, torflux(spres_ped))
    ! for consistency ...
    ! Also, pedge is re-used here; this actually is the pressure in the whole volume...
    IF (si .gt. spres_ped) THEN
       pedge = pmass(spres_ped)
    ELSE
       pedge = pmass(tf)
    END IF
    mass(i) = pedge*(ABS(vpnorm)*r00)**gamma
  END DO

  pres(:ns+1) = 0

  DO i = 1, ns

     ! sqrt(s_i) for half-grid s values (shifted inwards by -0.5 grid points)
     si = hs*ABS(i-1.5_dp)
     shalf(i:nrzt:ns) = SQRT(si)

     ! sqrt(s_i) for full-grid s values
     si = hs*(i-1)
     sqrts(i:nrzt:ns) = SQRT(si)

     ! TODO: what is this _exactly_?
     ! just observing, this is a linear profile on the half-grid
     ! which is at 2*pdamp == 0.1 at the axis
     ! and at              ca. 0  at the LCFS
     bdamp(i) = 2*pdamp*(1-si)
  END DO

  ! Avoid round-off
  sqrts(ns:nrzt:ns) = 1 ! boundary  value
  shalf(nrzt+1) = 1 ! no scaling for the weird hidden value at the end of shalf (see ndim in allocate_ns)
  sqrts(nrzt+1) = 1 ! no scaling for the weird hidden value at the end of sqrts (see ndim in allocate_ns)

  ! sm, sp used mainly in 1d preconditioner?
  DO i = 2,ns
     sm(i) = shalf(i)/sqrts(i)
     sp(i) = shalf(i+1)/sqrts(i)
  ENDDO

  sm(1) = 0
  sp(0) = 0
  sp(1) = sm(2)

  if (open_dbg_context("profil1d")) then

    call add_real("torflux_edge", torflux_edge)
    call add_real("polflux_edge", polflux_edge)
    call add_real("r00", r00)
    call add_real("lamscale", lamscale)
    call add_real("currv", currv)
    call add_real("Itor", Itor)

    call add_real_1d("shalf", ns+1, shalf(:ns+1))
    call add_real_1d("phips", ns,   phips(2:ns+1))
    call add_real_1d("chips", ns-1, chips(2:ns))
    call add_real_1d("iotas", ns-1, iotas(2:ns))
    call add_real_1d("icurv", ns-1, icurv(2:ns))
    call add_real_1d("mass",  ns-1, mass(2:ns))

    call add_real_1d("sqrts", ns,   sqrts(:ns))
    call add_real_1d("phipf", ns,   phipf)
    call add_real_1d("chipf", ns,   chipf)
    call add_real_1d("iotaf", ns,   iotaf)
    call add_real_1d("bdamp", ns,   bdamp)

    call add_real_1d("sm",    ns,   sm)
    call add_real_1d("sp",    ns+1, sp)

    call close_dbg_out()
  end if

END SUBROUTINE profil1d
