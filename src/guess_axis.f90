!> \file
!> \brief Computes guess for magnetic axis if user guess leads to initial sign change of Jacobian.

!> \brief Computes guess for magnetic axis if user guess leads to initial sign change of Jacobian.
!>
!> @param r1 \f$R\f$
!> @param z1 \f$Z\f$
!> @param ru0 \f$\partial R / \partial \theta\f$
!> @param zu0 \f$\partial Z / \partial \theta\f$
SUBROUTINE guess_axis(r1, z1, ru0, zu0)
  USE vmec_main
  USE vmec_params, ONLY: nscale, signgs
  USE realspace, ONLY: sqrts
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,nzeta,ntheta3,0:1), INTENT(in) :: r1, z1
  REAL(rprec), DIMENSION(ns,nzeta,ntheta3),     INTENT(in) :: ru0, zu0

  INTEGER, PARAMETER :: limpts = 61
  REAL(rprec), PARAMETER :: p5 = 0.5_dp, two = 2

  INTEGER :: iv, iu, iu_r, ivminus, nlim, ns12, klim, n
  REAL(rprec), DIMENSION(nzeta) :: rcom, zcom
  REAL(rprec), DIMENSION(ntheta1) :: r1b, z1b, rub, zub
  REAL(rprec), DIMENSION(ntheta1) :: r12, z12
  REAL(rprec), DIMENSION(ntheta1) :: rs, zs, tau, ru12, zu12, tau0
  REAL(rprec) :: rlim, zlim
  REAL(rprec) :: rmax, rmin, zmax, zmin, dzeta
  REAL(rprec) :: ds, mintau, mintemp

  character(len=255) :: dump_filename
  logical            :: dump_guess_axis = .false.

  ! COMPUTES GUESS FOR MAGNETIC AXIS IF USER GUESS
  ! LEADS TO INITIAL SIGN CHANGE OF JACOBIAN. DOES A GRID
  ! SEARCH (irgrid, izgrid) IN EACH PHI-PLANE FOR POINTS WHICH
  ! YIELD A VALUE FOR THE JACOBIAN WITH THE CORRECT SIGN (SIGNGS)
  ! CHOOSES THE AXIS POSITION SO THE MIN VALUE OF THE JACOBIAN IS MAXIMIZED

  ns12 = (ns+1)/2

  ! debugging output from guess_axis
  if (dump_guess_axis) then
    write(dump_filename, 998) trim(input_extension)
    open(unit=42, file=trim(dump_filename), status="unknown")

    write(42, *) "# nzeta ntheta1, ns12"
    write(42, *) nzeta, ntheta1, ns12

    write(42, *) "# iv iu r1b z1b rub zub r12 z12 rs zs tau ru12 zu12 tau0"
  end if

  planes: DO iv = 1, nzeta

     IF (.not.lasym .and. iv.gt.nzeta/2+1) THEN
        rcom(iv) = rcom(nzeta+2-iv)
        zcom(iv) =-zcom(nzeta+2-iv)
        CYCLE
     END IF

     r1b(:ntheta3) = r1(ns,iv,:,0) + r1(ns,iv,:,1)
     z1b(:ntheta3) = z1(ns,iv,:,0) + z1(ns,iv,:,1)
     r12(:ntheta3) = r1(ns12,iv,:,0) + r1(ns12,iv,:,1)*sqrts(ns12)
     z12(:ntheta3) = z1(ns12,iv,:,0) + z1(ns12,iv,:,1)*sqrts(ns12)
     rub(:ntheta3) = ru0(ns,iv,:)
     zub(:ntheta3) = zu0(ns,iv,:)
     ru12(:ntheta3) =  p5*(ru0(ns,iv,:) + ru0(ns12,iv,:))
     zu12(:ntheta3) =  p5*(zu0(ns,iv,:) + zu0(ns12,iv,:))

     IF (.not.lasym) THEN

        ! USE Z(v,-u) = -Z(twopi-v,u), R(v,-u) = R(twopi-v,u)
        ! TO DO EXTEND R,Z, etc. OVER ALL THETA (NOT JUST 0,PI)

        ! (twopi-v)
        ivminus = MOD(nzeta + 1 - iv, nzeta) + 1
        DO iu = 1+ntheta2, ntheta1
           iu_r = ntheta1 + 2 - iu ! (ntheta1 + 1) - (iu - 1)
           r1b(iu) = r1(ns,ivminus,iu_r,0) + r1(ns,ivminus,iu_r,1)
           z1b(iu) =-(z1(ns,ivminus,iu_r,0) + z1(ns,ivminus,iu_r,1))
           r12(iu) =  r1(ns12,ivminus,iu_r,0) + r1(ns12,ivminus,iu_r,1)*sqrts(ns12)
           z12(iu) =-(z1(ns12,ivminus,iu_r,0) + z1(ns12,ivminus,iu_r,1)*sqrts(ns12))
           rub(iu) =-ru0(ns,ivminus,iu_r)
           zub(iu) = zu0(ns,ivminus,iu_r)
           ru12(iu)=-p5*(ru0(ns,ivminus,iu_r) + ru0(ns12,ivminus,iu_r))
           zu12(iu)= p5*(zu0(ns,ivminus,iu_r) + zu0(ns12,ivminus,iu_r))
        END DO

     END IF

     ! Scan over r-z grid for interior point
     rmin = MINVAL(r1b)
     rmax = MAXVAL(r1b)
     zmin = MINVAL(z1b)
     zmax = MAXVAL(z1b)
     rcom(iv) = (rmax + rmin)/2
     zcom(iv) = (zmax + zmin)/2

     ! Estimate jacobian based on boundary and 1/2 surface
     ds = (ns - ns12)*hs
     DO iu = 1, ntheta1
        rs(iu) = (r1b(iu) - r12(iu))/ds + r1(1,iv,1,0)
        zs(iu) = (z1b(iu) - z12(iu))/ds + z1(1,iv,1,0)
        tau0(iu) = ru12(iu)*zs(iu) - zu12(iu)*rs(iu)
     END DO

     mintau = 0

     DO nlim = 1, limpts
        zlim = zmin + ((zmax - zmin)*(nlim-1))/(limpts-1)
        IF (.not.lasym .and. (iv.eq.1 .or. iv.eq.nzeta/2+1)) THEN
           zlim = 0
           IF (nlim .gt. 1) then
              EXIT
           end if
        END IF

        ! Find value of magnetic axis that maximizes the minimum jacobian value
        DO klim = 1, limpts
           rlim = rmin + ((rmax - rmin)*(klim-1))/(limpts-1)
           tau = signgs*(tau0 - ru12(:)*zlim + zu12(:)*rlim)

           mintemp = MINVAL(tau)
           IF (mintemp .gt. mintau) THEN
              mintau = mintemp
              rcom(iv) = rlim
              zcom(iv) = zlim
           ELSE IF (mintemp .eq. mintau) THEN
              ! If up-down symmetric and lasym=T, need this to pick z = 0
              IF (ABS(zcom(iv)).gt.ABS(zlim)) then
                 zcom(iv) = zlim
              end if
           END IF
        END DO
     END DO

     if (dump_guess_axis) then
       do iu = 1, ntheta1
         write(42, *) iv, iu,                                    &
           r1b(iu), z1b(iu), rub(iu), zub(iu),                   &
           r12(iu), z12(iu),                                     &
           rs(iu), zs(iu), tau(iu), ru12(iu), zu12(iu), tau0(iu)
       end do
     end if

  END DO planes

  ! FOURIER TRANSFORM RCOM, ZCOM
  dzeta = two/nzeta
  DO n = 0, ntor
     raxis_cc(n) = dzeta*SUM(cosnv(:,n)*rcom(:))/nscale(n)
     zaxis_cs(n) =-dzeta*SUM(sinnv(:,n)*zcom(:))/nscale(n)
     raxis_cs(n) =-dzeta*SUM(sinnv(:,n)*rcom(:))/nscale(n)
     zaxis_cc(n) = dzeta*SUM(cosnv(:,n)*zcom(:))/nscale(n)
     IF (n.eq.0 .or. n.eq.nzeta/2) THEN
        raxis_cc(n) = p5*raxis_cc(n)
        zaxis_cc(n) = p5*zaxis_cc(n)
     END IF
  END DO

  if (dump_guess_axis) then
    write(42, *) "# iv rcom zcom"
    do iv=1, nzeta
      write(42, *) iv, rcom(iv), zcom(iv)
    end do

    write(42, *) "# n raxis_cc zaxis_cs raxis_cs zaxis_cc"
    DO n = 0, ntor
      write(42, *) n, raxis_cc(n), zaxis_cs(n), raxis_cs(n), zaxis_cc(n)
    end do

    close(42)

    print *, "dumped guess_axis output to '"//trim(dump_filename)//"'"
    stop
  end if
998 format('guess_axis.',a)

END SUBROUTINE guess_axis