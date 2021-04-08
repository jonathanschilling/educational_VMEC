!> \file
SUBROUTINE add_fluxes(overg, bsupu, bsupv, lcurrent)
  USE vmec_main
  USE realspace, ONLY: wint, guu, guv, chip

  IMPLICIT NONE

  REAL(rprec), DIMENSION(nrzt), INTENT(in)    :: overg
  REAL(rprec), DIMENSION(nrzt), INTENT(inout) :: bsupu, bsupv
  LOGICAL, INTENT(in) :: lcurrent

  REAL(rprec), PARAMETER :: p5=0.5_dp, c1p5=1.5_dp
  REAL(rprec), PARAMETER :: iotaped = 0.10

  INTEGER :: js, l
  REAL(rprec) :: top, bot

  ! ADD MAGNETIC FLUX (CHIP, PHIP) TERMS TO BSUPU=-OVERG*LAM_V, BSUPV=OVERG*LAM_U
  ! COMPUTE FLUX FROM ITOR = <B_u>, ITOR(s) = integrated toroidal current (icurv)
  ! IF ncurr == 1
  !IF (.not.lcurrent .or. ncurr.eq.0) GOTO 100
  IF (lcurrent .and. ncurr.ne.0) then
     ! given current profile and lcurrent set --> compute fluxes, iota consistent with given current profile
     DO js = 2, ns
        top = icurv(js)
        bot = 0
        DO l = js, nrzt, ns
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
  end if

  IF (ncurr .eq. 0) THEN
     ! given iota profile: compute chips from iotas, phips
     chips = iotas*phips
  ELSE IF (.not.lcurrent) THEN
     ! given current profile, but .not. lcurrent (???): compute iotas from chips, phips
     WHERE (phips .ne. zero) iotas = chips/phips
  END IF

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

  ! what is this?
  bsupu(:nrzt) = bsupu(:nrzt)+chip(:nrzt)*overg(:nrzt)

END SUBROUTINE add_fluxes
