!> \file
!> \brief Pre-compute constant terms required for NESTOR.

!> \brief Pre-compute constant terms required for NESTOR.
!>
SUBROUTINE precal
  USE vacmod

  use dbgout

  IMPLICIT NONE

  REAL(rprec), PARAMETER :: p25 = p5*p5
  REAL(rprec), PARAMETER :: bigno = 1.e50_dp
  REAL(rprec), PARAMETER :: epstan = EPSILON(one)

  INTEGER :: kp, ku, kuminus, kv, kvminus, i, m, n, mn,            &
             imn, jmn, kmn, l, istat1, smn
  REAL(rprec), DIMENSION(0:mf + nf,0:mf,0:nf) :: cmn
  REAL(rprec) :: argu, argv, argp, dn1, f1, f2, f3, alp_per


  ! THIS ROUTINE COMPUTES INITIAL CONSTANTS AND ARRAYS

  pi2 = 8*ATAN(one) ! 2 pi
  pi3 = p5*pi2**3   ! 4 pi^3 = pi * (2 pi)^2
  pi4 = 2*pi2       ! 4 pi
  onp = one/nfper   ! 1/nfp
  onp2 = onp*onp    ! 1/nfp^2
  alu = pi2/nu      ! (2 pi)/nu
  alv = pi2/nv      ! (2 pi)/nv
  alp = pi2*onp     ! (2 pi)/nfp
  alvp = onp*alv    ! (2 pi)/(nv * nfp)

  alp_per = pi2/nvper ! (2 pi)/nvper

  DO kp = 1, nvper
     cosper(kp) = COS(alp_per*(kp - 1))
     sinper(kp) = SIN(alp_per*(kp - 1))
  END DO

  ! IMIRR(I) GIVES THE INDEX OF THE POINT (TWOPI-THETA(I), TWOPI-ZETA(I))
  ! k(u,v)minus count from the other side of the grid
  ! imirr is then the linear index correponding to (kuminus, kvminus)
  DO ku = 1, nu
     kuminus = MOD(nu - (ku-1), nu) + 1
     DO kv = 1, nv
        kvminus = MOD(nv - (kv-1), nv) + 1

        i        = kv      + nv*(ku      - 1)
        imirr(i) = kvminus + nv*(kuminus - 1)

        cosuv(i) = COS(alvp*(kv - 1))
        sinuv(i) = SIN(alvp*(kv - 1))
     END DO
  END DO

  ! NOTE: ANGLE DIFFERENCE IS PI2*{[NUV + (KUP-1)] - (KU-1)}
  !       THIS DIFFERENCE IS ACCOUNTED FOR BY THE OFFSET IUOFF IN GREENF ROUTINE
  !
  !       THE KP SUM BELOW IS USED ONLY FOR NV == 1. IT PERFORMS THE V-INTEGRAL
  !       IN AN AXISYMMETRIC PLASMA

  if (nv .eq. 1) then
    ! Tokamak; nvper = 64
    DO kp = 1, nvper
      argp = p5*alp_per*(kp-1)
      IF (ABS(argp - p25*pi2) < epstan) THEN
         tanv_1d(kp) = bigno
      ELSE
         tanv_1d(kp) = 2*TAN(argp)
      ENDIF
    end do
  else
    ! Stellarator
    DO kv = 1, nv
       argv = p5*alv*(kv - 1)
       IF (ABS(argv - p25*pi2) < epstan) THEN
          tanv_1d(kv) = bigno
       ELSE
          tanv_1d(kv) = 2*TAN(argv)
       ENDIF
    end do
  end if

  DO ku = 1, 2*nu
    argu = p5*alu*(ku - 1)
    IF (ABS(argu - p25*pi2)<epstan .or. ABS(argu - 0.75_dp*pi2) < epstan) THEN
       tanu_1d(ku) = bigno
    ELSE
       tanu_1d(ku) = 2*TAN(argu)
    ENDIF
  end do

  i = 0
  DO kp = 1, nvper
     IF (kp.gt.1 .and. nv.ne.1) EXIT
     argp = p5*alp_per*(kp-1)
     DO ku = 1, 2*nu
        argu = p5*alu*(ku - 1)
        DO kv = 1, nv
           i = i + 1
           argv = p5*alv*(kv - 1) + argp

           IF (ABS(argu - p25*pi2)<epstan .or. ABS(argu - 0.75_dp*pi2) < epstan) THEN
              tanu(i) = bigno
           ELSE
              tanu(i) = 2*TAN(argu)
           ENDIF

           IF (ABS(argv - p25*pi2) < epstan) THEN
              tanv(i) = bigno
           ELSE
              tanv(i) = 2*TAN(argv)
           ENDIF

           ! print *, i, kp, argp, ku, argu, kv, argv, tanu(i), tanv(i)
        END DO
     END DO
  END DO

  ! poloidal Fourier transform basis functions
  DO m = 0, mf
     l40: DO ku = 1, nu
        cosu(m,ku) = COS(alu*(m*(ku - 1)))
        sinu(m,ku) = SIN(alu*(m*(ku - 1)))
        DO kv = 1, nv
           i = kv + nv*(ku - 1)
           IF (i > nuv2) then
              ! This crops the linear index i to the stellarator-symmetric half of a toroidal module.
              ! ku == poloidal is the slow direction, so there the index ku will only run up to nu2.
              ! Above ku loop was probably not directly ended at nu2 to have the full set of Fourier basis functions in cosu, sinu.
              CYCLE l40
           end if
           cosu1(i,m) = cosu(m,ku)
           sinu1(i,m) = sinu(m,ku)
        END DO
     END DO l40
     DO ku = 1, nu2
        cosui(m,ku) = cosu(m,ku)*alu*alv*2
        sinui(m,ku) = sinu(m,ku)*alu*alv*2
        IF (ku.eq.1 .or. ku.eq.nu2) then
           ! wint-like functionality to take stellarator-symmetric
           ! "folding-over" into first half-period into account
           cosui(m,ku) = p5*cosui(m,ku)
        end if
     END DO
  END DO

  ! toroidal Fourier transform basis functions
  DO n = -nf, nf
     dn1 = alvp*(n*nfper)
     csign(n) = SIGN(one,dn1)
     l50: DO ku = 1, nu
        DO kv = 1, nv
           i = kv + nv*(ku - 1)
           cosv(n,kv) = COS(dn1*(kv - 1))
           sinv(n,kv) = SIN(dn1*(kv - 1))
           IF (i.gt.nuv2 .or. n.lt.0) then
              ! This crops the linear index i to the stellarator-symmetric half of a toroidal module.
              ! ku == poloidal is the slow direction, so there the index ku will only run up to nu2.
              ! Above ku loop was probably not directly ended at nu2 to have the full set of Fourier basis functions in cosv, sinv.
              CYCLE  l50
           end if
           cosv1(i,n) = cosv(n,kv)
           sinv1(i,n) = sinv(n,kv)
        END DO
     END DO l50
  END DO

  ! Fourier mode number arrays for potsin, potcos
  mn = 0
  DO n = -nf, nf
     DO m = 0, mf
        mn = mn + 1
        xmpot(mn) = m
        xnpot(mn) = n*nfper
     END DO
  END DO

  ! COMPUTE CMNS AND THE COEFFICIENTS OF T+- IN EQ (A14 AND A13) IN J.COMP.PHYS PAPER (PKM)
  ! NOTE: HERE, THE INDEX L IN THE LOOP BELOW IS THE SUBSCRIPT OF T+-. THEREFORE,
  ! L = 2L' + Kmn (L' = INDEX IN EQ. A14, Kmn = |m-n|), WITH LMIN = K AND LMAX = Jmn == m+n.
  !
  ! THE FOLLOWING DEFINITIONS PERTAIN (NOTE: kmn <= L <= jmn):
  !
  ! F1 = [(L + jmn)/2]! / [(jmn - L)/2]! == [(jmn + kmn)/2 + L']!/[(jmn - kmn)/2 + L']!
  !
  ! F2 = [(L + kmn)/2]!  == (L' + kmn)!
  !
  ! F3 = [(L - kmn)/2]!  == (L')!

  DO m = 0, mf
     DO n = 0, nf
        jmn = m + n
        imn = m - n
        kmn = ABS(imn)

        ! Integer: J+K always even
        smn = (jmn + kmn)/2
        f1 = 1
        f2 = 1
        f3 = 1
        DO i = 1, kmn
           f1 = f1*(smn + 1 - i)
           f2 = f2*i
        END DO
        cmn(0:mf+nf,m,n) = 0
        DO l = kmn, jmn, 2
           cmn(l,m,n) = f1/(f2*f3)*((-1)**((l - imn)/2))
           f1 = f1*p25*((jmn + l + 2)*(jmn - l))
           f2 = f2*p5*(l + 2 + kmn)
           f3 = f3*p5*(l + 2 - kmn)
        END DO
     END DO
  END DO

  ! Now combine these into a single coefficient (cmns), Eq. A13).
  ! NOTE:  The ALP=2*pi/nfper factor is needed to normalize integral over field periods

  DO m = 1,mf
     DO n = 1,nf
        cmns(0:mf+nf,m,n) = p5*alp*(  cmn(0:mf+nf,m,n  ) + cmn(0:mf+nf,m-1,n  ) &
                                    + cmn(0:mf+nf,m,n-1) + cmn(0:mf+nf,m-1,n-1))
     END DO
  END DO
  cmns(0:mf+nf,1:mf,0) = (p5*alp)*(cmn(0:mf+nf,1:mf,0) + cmn(0:mf+nf,:mf-1,0))
  cmns(0:mf+nf,0,1:nf) = (p5*alp)*(cmn(0:mf+nf,0,1:nf) + cmn(0:mf+nf,0,:nf-1))
  cmns(0:mf+nf,0,0)    = (p5*alp)*(cmn(0:mf+nf,0,0)    + cmn(0:mf+nf,0,0))

  precal_done = .true.

  if (open_dbg_context("vac1n_precal", id=icall)) then

    call add_int("nvper", nvper)
    call add_int("nuv_tan", nuv_tan)

    call add_real_1d("cosper", nvper, cosper)
    call add_real_1d("sinper", nvper, sinper)

    call add_real_1d("tanu_1d", 2*nu,  tanu_1d)
    call add_real_1d("tanv_1d", nuv_tan/(2*nu), tanv_1d)

    call add_real_2d("tanu", 2*nu, nuv_tan/(2*nu), tanu)
    call add_real_2d("tanv", 2*nu, nuv_tan/(2*nu), tanv)

    call add_real_3d("cmns", mf+nf+1, mf+1, nf+1, cmns)

    call close_dbg_out()
  end if

END SUBROUTINE precal
