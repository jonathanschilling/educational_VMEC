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

  use dbgout

  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,nzeta,ntheta3,0:1), INTENT(in) :: r1, z1
  REAL(rprec), DIMENSION(ns,nzeta,ntheta3),     INTENT(in) :: ru0, zu0

  INTEGER, PARAMETER :: limpts = 61
  REAL(rprec), PARAMETER :: p5 = 0.5_dp, two = 2.0_dp

  INTEGER :: iv, iu, iu_r, ivminus, nlim, ns12, klim, n
!  REAL(rprec), DIMENSION(nzeta) :: rcom, zcom
!  REAL(rprec), DIMENSION(nzeta,ntheta1) :: r1b, z1b, rub, zub
!  REAL(rprec), DIMENSION(nzeta,ntheta1) :: r12, z12
!  REAL(rprec), DIMENSION(nzeta,ntheta1) :: rs, zs, tau, ru12, zu12, tau0

  REAL(rprec), DIMENSION(:),   allocatable :: rcom, zcom
  REAL(rprec), DIMENSION(:,:), allocatable :: r1b, z1b, rub, zub
  REAL(rprec), DIMENSION(:,:), allocatable :: r12, z12
  REAL(rprec), DIMENSION(:,:), allocatable :: rs, zs, tau, ru12, zu12, tau0

  REAL(rprec) :: rlim, zlim
  REAL(rprec) :: rmax, rmin, zmax, zmin, dzeta
  REAL(rprec) :: ds, mintau, mintemp

  integer :: ivmax

  if (.not.lasym) then
    ivmax = nzeta/2+1
  else
    ivmax = nzeta
  end if

  allocate(rcom(nzeta))
  allocate(zcom(nzeta))

  allocate(r1b (nzeta,ntheta1)); r1b  = cbig
  allocate(z1b (nzeta,ntheta1)); z1b  = cbig
  allocate(rub (nzeta,ntheta1)); rub  = cbig
  allocate(zub (nzeta,ntheta1)); zub  = cbig
  allocate(r12 (nzeta,ntheta1)); r12  = cbig
  allocate(z12 (nzeta,ntheta1)); z12  = cbig
  allocate(rs  (nzeta,ntheta1)); rs   = cbig
  allocate(zs  (nzeta,ntheta1)); zs   = cbig
  allocate(tau (nzeta,ntheta1)); tau  = cbig
  allocate(ru12(nzeta,ntheta1)); ru12 = cbig
  allocate(zu12(nzeta,ntheta1)); zu12 = cbig
  allocate(tau0(nzeta,ntheta1)); tau0 = cbig

  ! COMPUTES GUESS FOR MAGNETIC AXIS IF USER GUESS
  ! LEADS TO INITIAL SIGN CHANGE OF JACOBIAN. DOES A GRID
  ! SEARCH (irgrid, izgrid) IN EACH PHI-PLANE FOR POINTS WHICH
  ! YIELD A VALUE FOR THE JACOBIAN WITH THE CORRECT SIGN (SIGNGS)
  ! CHOOSES THE AXIS POSITION SO THE MIN VALUE OF THE JACOBIAN IS MAXIMIZED

  ns12 = (ns+1)/2

  planes: DO iv = 1, nzeta/2+1

     r1b(iv,:ntheta3) = r1(ns,iv,:,0) + r1(ns,iv,:,1)
     z1b(iv,:ntheta3) = z1(ns,iv,:,0) + z1(ns,iv,:,1)

     r12(iv,:ntheta3) = r1(ns12,iv,:,0) + r1(ns12,iv,:,1)*sqrts(ns12)
     z12(iv,:ntheta3) = z1(ns12,iv,:,0) + z1(ns12,iv,:,1)*sqrts(ns12)

     rub(iv,:ntheta3) = ru0(ns,iv,:)
     zub(iv,:ntheta3) = zu0(ns,iv,:)

     ru12(iv,:ntheta3) =  p5*(ru0(ns,iv,:) + ru0(ns12,iv,:))
     zu12(iv,:ntheta3) =  p5*(zu0(ns,iv,:) + zu0(ns12,iv,:))

     IF (.not.lasym) THEN

        ! USE Z(v,-u) = -Z(twopi-v,u), R(v,-u) = R(twopi-v,u)
        ! TO DO EXTEND R,Z, etc. OVER ALL THETA (NOT JUST 0,PI)

        ! (twopi-v)
        ivminus = MOD(nzeta - (iv - 1), nzeta) + 1
        DO iu = 1+ntheta2, ntheta1
           iu_r = ntheta1 + 2 - iu ! (ntheta1 + 1) - (iu - 1)

           r1b(iv,iu) =   r1(ns,  ivminus,iu_r,0) + r1( ns,  ivminus,iu_r,1)
           z1b(iv,iu) =-( z1(ns,  ivminus,iu_r,0) + z1( ns,  ivminus,iu_r,1))
           rub(iv,iu) =- ru0(ns,  ivminus,iu_r)
           zub(iv,iu) =  zu0(ns,  ivminus,iu_r)

           r12(iv,iu) =   r1(ns12,ivminus,iu_r,0) + r1( ns12,ivminus,iu_r,1)*sqrts(ns12)
           z12(iv,iu) =-( z1(ns12,ivminus,iu_r,0) + z1( ns12,ivminus,iu_r,1)*sqrts(ns12))
           ru12(iv,iu)=-(ru0(ns,  ivminus,iu_r)   + ru0(ns12,ivminus,iu_r))*p5
           zu12(iv,iu)= (zu0(ns,  ivminus,iu_r)   + zu0(ns12,ivminus,iu_r))*p5
        END DO
     END IF

     ! Scan over r-z grid for interior point
     rmin = MINVAL(r1b(iv,:))
     rmax = MAXVAL(r1b(iv,:))
     zmin = MINVAL(z1b(iv,:))
     zmax = MAXVAL(z1b(iv,:))

     ! initial guess for new axis: center of grid
     rcom(iv) = (rmax + rmin)/2.0_dp
     zcom(iv) = (zmax + zmin)/2.0_dp

     ! Estimate jacobian based on boundary and 1/2 surface
     ds = (ns - ns12)*hs
     DO iu = 1, ntheta1
        ! (1,iv,1,0) == axis, iv, theta=0, even-m
        rs(iv,iu) = (r1b(iv,iu) - r12(iv,iu))/ds + r1(1,iv,1,0)
        zs(iv,iu) = (z1b(iv,iu) - z12(iv,iu))/ds + z1(1,iv,1,0)
        tau0(iv,iu) = ru12(iv,iu)*zs(iv,iu) - zu12(iv,iu)*rs(iv,iu)
     END DO

     mintau = 0.0_dp

     DO nlim = 1, limpts
        zlim = zmin + ((zmax - zmin)*(nlim-1))/(limpts-1)

        IF (.not.lasym .and. (iv.eq.1 .or. iv.eq.nzeta/2+1)) THEN
           zlim = 0.0_dp
           IF (nlim .gt. 1) then
              EXIT
           end if
        END IF

        ! Find value of magnetic axis that maximizes the minimum jacobian value
        DO klim = 1, limpts
           rlim = rmin + ((rmax - rmin)*(klim-1))/(limpts-1)

           tau(iv,:) = signgs*(tau0(iv,:) - ru12(iv,:)*zlim + zu12(iv,:)*rlim)

           mintemp = MINVAL(tau(iv,:))
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

  END DO planes

  if (.not. lasym) then
     ! For a stellarator-symmetric case,
     ! mirror the rest of the toroidal range based on the first ~half.
     do iv = nzeta/2+2, nzeta
         rcom(iv) = rcom(nzeta+1-(iv-1))
         zcom(iv) =-zcom(nzeta+1-(iv-1))
     end do ! iv
  end if ! .not. lasym

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

  ! debugging output from guess_axis
  if (open_dbg_context("guess_axis")) then

    ! axis geometry on entry into this routine
    call add_real_1d("raxis_in" , ivmax, r1(1,1:ivmax,1,0) )
    call add_real_1d("zaxis_in" , ivmax, z1(1,1:ivmax,1,0) )

    call add_real_2d("r1b" , ivmax, ntheta1, r1b(1:ivmax,:) )
    call add_real_2d("z1b" , ivmax, ntheta1, z1b(1:ivmax,:) )
    call add_real_2d("rub" , ivmax, ntheta1, rub(1:ivmax,:) )
    call add_real_2d("zub" , ivmax, ntheta1, zub(1:ivmax,:) )

    call add_real_2d("r12" , ivmax, ntheta1, r12(1:ivmax,:) )
    call add_real_2d("z12" , ivmax, ntheta1, z12(1:ivmax,:) )
    call add_real_2d("ru12", ivmax, ntheta1, ru12(1:ivmax,:))
    call add_real_2d("zu12", ivmax, ntheta1, zu12(1:ivmax,:))

    call add_real_2d("rs"  , ivmax, ntheta1, rs(1:ivmax,:)  )
    call add_real_2d("zs"  , ivmax, ntheta1, zs(1:ivmax,:)  )

    call add_real_2d("tau0", ivmax, ntheta1, tau0(1:ivmax,:))

    ! Actually tau would need to be ivmax * limpts * limpts * ntheta1,
    ! since it gets set for every grid point.
    ! As of now, it contains the entries for the (rmax, zmax) grid point.
    call add_real_2d("tau" , ivmax, ntheta1, tau(1:ivmax,:) )

    call add_real_1d("rcom", nzeta, rcom)
    call add_real_1d("zcom", nzeta, zcom)

    call add_real_1d("raxis_cc", ntor+1, raxis_cc)
    call add_real_1d("zaxis_cs", ntor+1, zaxis_cs)
    call add_real_1d("raxis_cs", ntor+1, raxis_cs)
    call add_real_1d("zaxis_cc", ntor+1, zaxis_cc)

    call close_dbg_out()
  end if

  deallocate(rcom, zcom)
  deallocate(r1b, z1b, rub, zub, &
             r12, z12, &
             rs, zs, tau, ru12, zu12, tau0)

END SUBROUTINE guess_axis
