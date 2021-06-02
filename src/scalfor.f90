!> \file
!> \brief Build forces from different contributions

!> \brief Build forces from different contributions
!>
!> In below parameter names, x=R or Z.
!>
!> @param gcx force output
!> @param axm force contribution input
!> @param bxm force contribution input
!> @param axd force contribution input
!> @param bxd force contribution input
!> @param cx force contribution input
!> @param iflag subtract edge instability from preconditioner
SUBROUTINE scalfor(gcx, axm, bxm, axd, bxd, cx, iflag)

  USE vmec_main
  USE vmec_params
  USE vmec_dim, ONLY: ns

  IMPLICIT NONE

  INTEGER, INTENT(in) :: iflag
  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: gcx
  REAL(rprec), DIMENSION(ns+1,2), INTENT(in) :: axm, bxm, axd, bxd
  REAL(rprec), DIMENSION(ns), INTENT(in) :: cx

  REAL(rprec), PARAMETER :: ftol_edge = 1.e-9_dp
  REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
  REAL(rprec), PARAMETER :: fac=0.25_dp
  REAL(rprec), PARAMETER :: edge_pedestal= 0.05_dp
  INTEGER :: m , mp, n, js, jmax
  REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: ax, bx, dx
  REAL(rprec) :: mult_fac
  ! LOGICAL :: ledge ! improved convergence for free-boundary, see below

  character(len=255) :: dump_filename
  logical            :: dump_scalfor = .false.

  ALLOCATE (ax(ns,0:ntor,0:mpol1), bx(ns,0:ntor,0:mpol1), dx(ns,0:ntor,0:mpol1))

  ! clean state since not all entries are assigned
  ax = 0
  dx = 0
  bx = 0

  jmax = ns
  IF (ivac .lt. 1) jmax = ns1

  ! FOR SOME 3D PLASMAS, THIS SOMETIME HELPS (CHOOSE mult_fac =1 otherwise)
  ! TO AVOID JACOBIAN RESETS BY GIVING A SMOOTH TRANSITION FROM FIXED TO FREE ITERATIONS
  !  mult_fac = 1._dp/(1._dp + 10*(fsqr+fsqz))
  !  gcx(ns,:,:,:) = mult_fac*gcx(ns,:,:,:)

  DO m = 0, mpol1
     mp = MOD(m,2) + 1
     DO n = 0, ntor
        DO js = jmin2(m), jmax
           ax(js,n,m) = -(axm(js+1,mp) + bxm(js+1,mp)*m**2)
           bx(js,n,m) = -(axm(js,mp) + bxm(js,mp)*m**2)
           dx(js,n,m) = -(axd(js,mp) + bxd(js,mp)*m**2 + cx(js)*(n*nfp)**2)
        END DO

        IF (m .eq. 1) THEN
           dx(2,n,m) = dx(2,n,m) + bx(2,n,m)
           ! OFF 050311
           ! DO js = jmin2(m), jmax
           !    ax(js,n,m) = c1p5*ax(js,n,m)
           !    bx(js,n,m) = c1p5*bx(js,n,m)
           !    dx(js,n,m) = c1p5*dx(js,n,m)
           ! END DO
        END IF
     END DO
  END DO

  IF (jmax .ge. ns) THEN

    ! SMALL EDGE PEDESTAL NEEDED TO IMPROVE CONVERGENCE
    ! IN PARTICULAR, NEEDED TO ACCOUNT FOR POTENTIAL ZERO
    ! EIGENVALUE DUE TO NEUMANN (GRADIENT) CONDITION AT EDGE
     dx(ns,:,0:1)     = (1+  edge_pedestal)*dx(ns,:,0:1)
     dx(ns,:,2:mpol1) = (1+2*edge_pedestal)*dx(ns,:,2:mpol1)

     ! STABILIZATION ALGORITHM FOR ZC_00(NS)
     ! FOR UNSTABLE CASE, HAVE TO FLIP SIGN OF -FAC -> +FAC FOR CONVERGENCE
     ! COEFFICIENT OF < Ru (R Pvac)> ~ -fac*(z-zeq) WHERE fac (EIGENVALUE, OR
     ! FIELD INDEX) DEPENDS ON THE EQUILIBRIUM MAGNETIC FIELD AND CURRENT,
     ! AND zeq IS THE EQUILIBRIUM EDGE VALUE OF Z00
      mult_fac = MIN(fac, fac*hs*15)

      IF (iflag .eq. 1) THEN
         ! METHOD 1: SUBTRACT (INSTABILITY) Pedge ~ fac*z/hs FROM PRECONDITIONER AT EDGE
         dx(ns,0,0) = dx(ns,0,0)*(1-mult_fac)/(1+edge_pedestal)
      END IF
  ENDIF

  ! ACCELERATE (IMPROVE) CONVERGENCE OF FREE BOUNDARY.
  ! THIS WAS ADDED TO DEAL WITH CASES WHICH MIGHT OTHERWISE DIVERGE.
  ! BY DECREASING THE FSQ TOLERANCE LEVEL WHERE THIS KICKS IN (FTOL_EDGE),
  ! THE USER CAN TURN-OFF THIS FEATURE
  !
  ! DIAGONALIZE (DX DOMINANT) AND REDUCE FORCE (DX ENHANCED) AT EDGE
  ! TO IMPROVE CONVERGENCE FOR N != 0 TERMS

  ! ledge = .false.
  ! IF ((fsqr+fsqz) .lt. ftol_edge) ledge = .true.
  ! IF ((iter2-iter1).lt.400 .or. ivac.lt.1) ledge = .false.

  ! IF (ledge) THEN
  !    dx(ns,1:,1:) = 3*dx(ns,1:,1:)
  ! END IF

  ! check scalfor state == inputs to tridslv
  if (dump_scalfor) then

    ! prior knowledge about how this is called:
    ! iflag=0 --> R
    ! iflag=1 --> Z

    if (iflag.eq.0) then
      write(dump_filename, 995) trim(input_extension)
    elseif (iflag.eq.1) then
      write(dump_filename, 996) trim(input_extension)
    end if
995 format('scalfor_R.',a)
996 format('scalfor_Z.',a)

    open(unit=42, file=trim(dump_filename), status="unknown")

    write(42, *) "# ns ntor mpol1"
    write(42, *) ns, ntor, mpol1

    write(42, *) "# j n m ax dx bx"
    do js=1, ns
      do n=0, ntor
        do m=0, mpol1
          write(42, *) js, n, m, &
            ax(js, n, m), &
            dx(js, n, m), &
            bx(js, n, m)
        end do
      end do
    end do

    print *, "dumped scalfor state to '"//trim(dump_filename)//"'"
    close(42)

    ! should only stop after also Z has been dumped
    if (iflag.eq.1) then
      stop
    end if
  end if

  ! SOLVES AX(I)*X(I+1) + DX(I)*X(I) + BX(I)*X(I-1) = GCX(I), I=JMIN3,JMAX AND RETURNS ANSWER IN GCX(I)
  CALL tridslv (ax, dx, bx, gcx, jmin3, jmax, mnsize-1, ns, ntmax)

  DEALLOCATE (ax, bx, dx)

END SUBROUTINE scalfor
