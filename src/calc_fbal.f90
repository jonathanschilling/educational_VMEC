!> \file
SUBROUTINE calc_fbal(bsubu, bsubv)
  USE vmec_main, ONLY: buco, bvco, equif,             &
                       jcurv, jcuru, chipf, vp, pres, &
                       phipf, vpphi, presgrad, ohs
  USE vmec_params, ONLY: signgs
  USE vmec_dim, ONLY: ns, nrzt, ns1
  USE realspace, ONLY: wint
  USE stel_kinds, ONLY: dp

  implicit none

  REAL(dp), INTENT(in) :: bsubu(1:nrzt)
  REAL(dp), INTENT(in) :: bsubv(1:nrzt)

  INTEGER  :: js

  ! compute flux-surface averages of covariant magnetic field components
  DO js = 2, ns
     buco(js) = SUM(bsubu(js:nrzt:ns)*wint(js:nrzt:ns)) ! toroidal current I (?)
     bvco(js) = SUM(bsubv(js:nrzt:ns)*wint(js:nrzt:ns)) ! poloidal current G (?)
  END DO

  ! FROM AMPERE'S LAW, JcurX are angle averages of jac*JsupX, so
  !                    JcurX = (dV/ds)/twopi**2 <JsupX> where <...> is flux surface average
  DO js = 2, ns1
     ! radial derivatives of covariant magnetic field components give contravariant currents
     ! --> poloidal, toroidal derivative contributions to curl(B) are averaged out in flux-surface average ?
     jcurv(js) = (signgs*ohs)*(buco(js+1) - buco(js))
     jcuru(js) =-(signgs*ohs)*(bvco(js+1) - bvco(js))

     vpphi(js) = (vp(js+1) + vp(js))/2

     ! prescribed pressure gradient from user input
     presgrad(js) = (pres(js+1) - pres(js))*ohs

     ! total resulting force-imbalance: magnetic pressure + (prescribed) kinetic pressure
     ! F = grad(p) - <j x B> (or something like this...)
     equif(js) = (-phipf(js)*jcuru(js) + chipf(js)*jcurv(js))/vpphi(js) + presgrad(js)
  END DO

  equif(1) = 0
  equif(ns) = 0

END SUBROUTINE calc_fbal
