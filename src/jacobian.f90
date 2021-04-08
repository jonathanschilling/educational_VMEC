!> \file
SUBROUTINE jacobian
  USE vmec_main, ONLY: ohs, nrzt, irst
  USE vmec_params, ONLY: meven, modd
  USE realspace
  USE vmec_dim, ONLY: ns
  USE vforces, r12 => armn_o, ru12 => azmn_e, zu12 => armn_e, &
               rs => bzmn_e, zs => brmn_e, tau => azmn_o
  IMPLICIT NONE

  REAL(rprec), PARAMETER :: zero=0, p5=0.5_dp, p25=p5*p5

  INTEGER :: l
  REAL(rprec) :: taumax, taumin, dshalfds=p25, temp(nrzt/ns)

  ! (RS, ZS)=(R, Z) SUB S, (RU12, ZU12)=(R, Z) SUB THETA(=U)
  ! AND TAU=SQRT(G)/R ARE DIFFERENCED ON HALF MESH
  !
  ! SQRT(G) = R*TAU IS COMPUTED IN BCOVAR
  !
  ! HERE, TAU = (Ru * Zs - Rs * Zu).
  ! THE DERIVATIVES OF SHALF = SQRT(s) WERE COMPUTED EXPLICITLY AS: d(shalf)/ds = .5/shalf

  ! initially, all good
  irst = 1

  DO l = 2,nrzt
    zs(l)   = ohs*( z1(l,meven) - z1(l-1,meven) + shalf(l)*(z1(l,modd)  - z1(l-1,modd)) )
    rs(l)   = ohs*( r1(l,meven) - r1(l-1,meven) + shalf(l)*(r1(l,modd)  - r1(l-1,modd)) )
    r12(l)  =  p5*( r1(l,meven) + r1(l-1,meven) + shalf(l)*(r1(l,modd)  + r1(l-1,modd)) )
    ru12(l) =  p5*( ru(l,meven) + ru(l-1,meven) + shalf(l)*(ru(l,modd)  + ru(l-1,modd)) )
    zu12(l) =  p5*( zu(l,meven) + zu(l-1,meven) + shalf(l)*(zu(l,modd)  + zu(l-1,modd)) )

    tau(l)  = ru12(l)*zs(l) - rs(l)*zu12(l) + dshalfds*                               &
              (      ru(l,modd) *z1(l,modd) + ru(l-1,modd) *z1(l-1,modd)              &
                   - zu(l,modd) *r1(l,modd) - zu(l-1,modd) *r1(l-1,modd)              &
               + (   ru(l,meven)*z1(l,modd) + ru(l-1,meven)*z1(l-1,modd)              &
                   - zu(l,meven)*r1(l,modd) - zu(l-1,meven)*r1(l-1,modd) )/shalf(l) )
  END DO

  ! TEST FOR SIGN CHANGE IN JACOBIAN
  temp(:) = tau(2:nrzt:ns)
  tau(1:nrzt:ns) = temp(:)
  taumax = MAXVAL(tau(2:nrzt))
  taumin = MINVAL(tau(2:nrzt))
  IF (taumax*taumin .lt. zero) then
     ! bad jacobian !
     irst = 2
  end if

END SUBROUTINE jacobian
