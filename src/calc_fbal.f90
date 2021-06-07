!> \file
!> \brief Compute flux-surface averaged radial force balance \f$\nabla p\, - <\mathbf{j} \times \mathbf{B}>\f$.

!> \brief Compute flux-surface averaged radial force balance \f$\nabla p\, - <\mathbf{j} \times \mathbf{B}>\f$.
!>
!> @param bsubu covariant component of magnetic field \f$B_\theta\f$ on half mesh
!> @param bsubv covariant component of magnetic field \f$B_\zeta\f$  on half mesh
SUBROUTINE calc_fbal(bsubu, bsubv)

  USE vmec_main, ONLY: buco, bvco, equif,             &
                       jcurv, jcuru, chipf, vp, pres, &
                       phipf, vpphi, presgrad, ohs,   &
                       input_extension, iter2
  USE vmec_params, ONLY: signgs
  USE vmec_dim, ONLY: ns, nrzt, ns1
  USE realspace, ONLY: wint
  USE stel_kinds, ONLY: dp

  implicit none

  REAL(dp), INTENT(in) :: bsubu(1:nrzt)
  REAL(dp), INTENT(in) :: bsubv(1:nrzt)

  INTEGER  :: js

  character(len=255) :: dump_filename
  logical            :: dump_calc_fbal = .false.


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


  ! check calc_fbal output
  if (dump_calc_fbal .and. iter2.le.2) then
    write(dump_filename, 998) ns, iter2, trim(input_extension)
998 format('calc_fbal_',i5.5,'_',i6.6,'.',a)

    open(unit=42, file=trim(dump_filename), status="unknown")

    write(42, *) "# ns"
    write(42, *) ns

    write(42, *) "# js buco bvco jcurv jcuru vpphi presgrad equif"
    DO js = 1, ns
      write(42,*) js, buco(js), bvco(js), jcurv(js), jcuru(js), &
                      vpphi(js), presgrad(js), equif(js)
    end do

    close(42)

    print *, "dumped calc_fbal output to '"//trim(dump_filename)//"'"
!     stop
  end if

END SUBROUTINE calc_fbal
