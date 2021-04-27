!> \file
SUBROUTINE add_fluxes(overg, bsupu, bsupv)
  USE vmec_main
  USE realspace, ONLY: wint, guu, guv, chip

  IMPLICIT NONE

  REAL(rprec), DIMENSION(nrzt), INTENT(in)    :: overg
  REAL(rprec), DIMENSION(nrzt), INTENT(inout) :: bsupu, bsupv

  REAL(rprec), PARAMETER :: p5=0.5_dp, c1p5=1.5_dp
  REAL(rprec), PARAMETER :: iotaped = 0.10

  INTEGER :: js, l
  REAL(rprec) :: top, bot

  ! I think this is the "zero-current algorithm" published in section 2.3 of
  ! Hirshman, Hogan, "ORMEC: A Three-Dimensional MHD Spectral Inverse Equilibrium Code" (1986), J. Comp. Phys. 63 (2), 329-352
  ! https://doi.org/10.1016/0021-9991(86)90197-X

  ! ADD MAGNETIC FLUX (CHIP, PHIP) TERMS TO BSUPU=-OVERG*LAM_V, BSUPV=OVERG*LAM_U
  ! COMPUTE FLUX FROM ITOR = <B_u>, ITOR(s) = integrated toroidal current (icurv)
  IF (ncurr .eq. 1) then
     ! given current profile and lcurrent set --> compute fluxes, iota consistent with given current profile
     DO js = 2, ns
        ! solve Eqn. (11) of the ORMEC paper for each flux surface
        top = icurv(js) ! offset: this makes the zero-current algorithm a constrained-current algorithm
        bot = 0
        DO l = js, nrzt, ns
           ! bsupu contains -d(lambda)/d(zeta)*lamscale on entry (?)
           ! bsupv contains  d(lambda)/d(theta)*lamscale on entry (?)
           top = top - wint(l)*(guu(l)*bsupu(l) + guv(l)*bsupv(l))
           bot = bot + wint(l)* guu(l)*overg(l)
        END DO
        IF (bot .ne. zero) then
           chips(js) = top/bot
        end if
        IF (phips(js) .ne. zero) then
           iotas(js) = chips(js)/phips(js)
        end if
     END DO
  else ! ncurr .eq. 0
     ! given iota profile: compute chips from iotas, phips
     chips = iotas*phips
  END IF

  ! distribute chips (ns-sized array) into larger chip array (over full surface) (?)
  DO js = 2, ns
     chip(js:nrzt:ns) = chips(js)
  END DO

  ! half-grid to full-grid for chi-prime and iota below

  chipf(2:ns1) = (chips(2:ns1) + chips(3:ns1+1))/2
  chipf(ns)    = 2*chips(ns)-chips(ns1)

  ! Do not compute iota too near origin
  iotaf(1)  = c1p5*iotas(2) - p5*iotas(3)     !zero gradient near axis
  iotaf(ns) = c1p5*iotas(ns) - p5*iotas(ns-1)
  DO js = 2, ns-1
     iotaf(js) = p5*(iotas(js) + iotas(js+1))
  END DO

  ! bsupu contains -dLambda/dZeta*lamscale and now needs to get chip/sqrt(g) added, as outlined in bcovar above the call to this routine.
  bsupu(:nrzt) = bsupu(:nrzt) + chip(:nrzt)*overg(:nrzt)

END SUBROUTINE add_fluxes
