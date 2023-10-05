!> \file
!> \brief Computes br, bphi, bz, bsubs on half-radial mesh.

!> \brief Computes br, bphi, bz, bsubs on half-radial mesh.
!>
!> @param r12 \f$R^2\f$
!> @param rs \f$\partial R / \partial s\f$
!> @param zs \f$\partial Z / \partial s\f$
!> @param ru12 \f$\partial R / \partial \theta\f$ on half-grid from jacobian()
!> @param zu12 \f$\partial Z / \partial \theta\f$ on half-grid from jacobian()
!> @param bsubs covariant component of magnetic field \f$B_s\f$
!> @param bsupu contravariant component of magnetic field \f$B^\theta\f$
!> @param bsupv contravariant component of magnetic field \f$B^\zeta\f$
!> @param br   cylindrical component of magnetic field \f$B^R\f$
!> @param bphi cylindrical component of magnetic field \f$B^\varphi\f$
!> @param bz   cylindrical component of magnetic field \f$B^Z\f$
SUBROUTINE bss(r12, rs, zs, ru12, zu12, bsubs, bsupu, bsupv, br, bphi, bz)
  USE vmec_main
  USE realspace

  use dbgout

  IMPLICIT NONE

  REAL(rprec), DIMENSION(nrzt), INTENT(in)  :: r12, rs, zs, ru12, zu12, bsupu, bsupv
  REAL(rprec), DIMENSION(nrzt), INTENT(out) :: br, bphi, bz, bsubs

  REAL(rprec), PARAMETER :: p5 = 0.5_dp
  REAL(rprec), PARAMETER :: p25 = p5*p5, dshalfds=p25

  INTEGER :: l
  REAL(rprec), DIMENSION(:), allocatable :: rv12, zv12, rs12, zs12, gsu, gsv

  ! Computes br, bphi, bz, bsubs on HALF-RADIAL mesh
  ! bsubs will be averaged onto the FULL-RADIAL mesh in jxbforce before output to WOUT file

  allocate(rv12(nrzt), zv12(nrzt), &
           rs12(nrzt), zs12(nrzt), &
           gsu(nrzt), gsv(nrzt))

  ! initialize first entry, as this is otherwise never set and floats around -> noise in Git!
  rv12(1) = zero
  zv12(1) = zero

  rs12(1) = zero
  zs12(1) = zero

  gsu(1) = zero
  gsv(1) = zero

  bsubs(1) = zero

  br(1)   = zero
  bphi(1) = zero
  bz(1)   = zero

  DO l = 2, nrzt

     rv12(l)  = p5*(rv(l,0)+rv(l-1,0) + shalf(l)*(rv(l,1) + rv(l-1,1)))
     zv12(l)  = p5*(zv(l,0)+zv(l-1,0) + shalf(l)*(zv(l,1) + zv(l-1,1)))

     ! -------------

     rs12(l)  = rs(l) + dshalfds*(r1(l,1) + r1(l-1,1))/shalf(l)
     zs12(l)  = zs(l) + dshalfds*(z1(l,1) + z1(l-1,1))/shalf(l)

     gsu(l)   = rs12(l)*ru12(l) + zs12(l)*zu12(l)
     gsv(l)   = rs12(l)*rv12(l) + zs12(l)*zv12(l)

     bsubs(l) = bsupu(l)*gsu(l) + bsupv(l)*gsv(l)

     ! -------------

     br(l)    = bsupu(l)*ru12(l) + bsupv(l)*rv12(l)
     bphi(l)  = bsupv(l)*r12(l)
     bz(l)    = bsupu(l)*zu12(l) + bsupv(l)*zv12(l)

  END DO

  if (open_dbg_context("bss", id=0)) then

    call add_real_3d("rv12",  ns, nzeta, ntheta3, rv12)
    call add_real_3d("zv12",  ns, nzeta, ntheta3, zv12)

    call add_real_3d("rs12",  ns, nzeta, ntheta3, rs12)
    call add_real_3d("zs12",  ns, nzeta, ntheta3, zs12)

    call add_real_3d("gsu",   ns, nzeta, ntheta3, gsu)
    call add_real_3d("gsv",   ns, nzeta, ntheta3, gsv)

    call add_real_3d("br",    ns, nzeta, ntheta3, br  )
    call add_real_3d("bphi",  ns, nzeta, ntheta3, bphi)
    call add_real_3d("bz",    ns, nzeta, ntheta3, bz  )

    call add_real_3d("bsubs", ns, nzeta, ntheta3, bsubs)

    call close_dbg_out()
  end if

  deallocate(rv12, zv12, rs12, zs12, gsu, gsv)

END SUBROUTINE bss
