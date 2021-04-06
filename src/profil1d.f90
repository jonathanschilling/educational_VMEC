!> \file
SUBROUTINE profil1d(xc, xcdot, lreset)
  USE vmec_main
  USE vmec_params, ONLY: signgs, lamscale, rcc, pdamp
  USE realspace, ONLY: shalf, sqrts
  USE init_geometry, ONLY: lflip
  IMPLICIT NONE

  REAL(rprec), DIMENSION(neqs2), INTENT(out) :: xc, xcdot
  LOGICAL, INTENT(in) :: lreset

  REAL(rprec), PARAMETER :: c1p5 = 1.5_dp

  INTEGER :: i
  REAL(rprec) :: Itor, si, tf, pedge, vpnorm, torflux_edge, polflux_edge

  REAL(rprec), EXTERNAL :: pcurr, pmass, piota, &
                           torflux, torflux_deriv, polflux, polflux_deriv

  ! COMPUTE PHIP, IOTA PROFILES ON FULL-GRID
  ! COMPUTE MASS PROFILE ON HALF-GRID
  ! BY READING INPUT COEFFICIENTS.
  ! PRESSURE CONVERTED TO INTERNAL UNITS BY MULTIPLICATION BY mu0 = 4*pi*10**-7

  torflux_edge = signgs * phiedge / twopi
  si = torflux(one)
  IF (si .ne. zero) torflux_edge = torflux_edge/si
  polflux_edge = torflux_edge
  si = polflux(one)
  IF (si .ne. zero) polflux_edge = polflux_edge/si
  r00 = rmn_bdy(0,0,rcc)

  phips(1) = 0
  chips(1) = 0
  icurv(1) = 0

  DO i = 2,ns
     si = hs*(i-c1p5)
     tf = MIN(one, torflux(si))
     phips(i) = torflux_edge * torflux_deriv(si)
     chips(i) = torflux_edge * polflux_deriv(si)
     iotas(i) = piota(tf)
     icurv(i) = pcurr(tf)
  END DO

  ! Compute lamscale factor for "normalizing" lambda (needed for scaling hessian)
  lamscale = SQRT(hs*SUM(phips(2:ns)**2))
  phips(ns+1) = 2*phips(ns)-phips(ns-1)
  IF (lamscale .EQ. 0) STOP 'PHIP == 0: ERROR!'

  IF (lflip) THEN
     iotas = -iotas
     chips = -chips
  END IF

  DO i = 1,ns
     si = hs*(i-1)
     tf = MIN(one, torflux(si))
     iotaf(i) = piota(tf)
     phipf(i) = torflux_edge * torflux_deriv(si)
     chipf(i) = torflux_edge * polflux_deriv(si)
  ENDDO

  ! SCALE CURRENT TO MATCH INPUT EDGE VALUE, CURTOR
  ! FACTOR OF SIGNGS NEEDED HERE, SINCE MATCH IS MADE TO LINE
  ! INTEGRAL OF BSUBU (IN GETIOTA) ~ SIGNGS * CURTOR
  pedge = pcurr(one)
  Itor = 0
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
    IF (si .gt. spres_ped) THEN
       pedge = pmass(spres_ped)
    ELSE
       pedge = pmass(tf)
    END IF
    mass(i) = pedge*(ABS(vpnorm)*r00)**gamma
  END DO

  pres(:ns+1) = 0
  xcdot(:neqs2) = 0

  DO i = 1, ns
     si = hs*ABS(i-1.5_dp)
     shalf(i:nrzt:ns) = SQRT(si)
     si = hs*(i-1)
     sqrts(i:nrzt:ns) = SQRT(si)
     bdamp(i) = 2*pdamp*(1-si)
  END DO

  ! Avoid round-off
  sqrts(ns:nrzt:ns) = 1
  shalf(nrzt+1) = 1
  sqrts(nrzt+1) = 1

  DO i = 2,ns
     sm(i) = shalf(i)/sqrts(i)
     sp(i) = shalf(i+1)/sqrts(i)
  ENDDO

  sm(1) = 0
  sp(0) = 0
  sp(1) = sm(2)

  IF (lreset) THEN
    xc(:neqs1) = 0
  END IF

END SUBROUTINE profil1d
