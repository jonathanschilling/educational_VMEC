!> \file
!> \brief allocate and fill some fixed-size arrays (only depending on Fourier resolution).

!> \brief allocate and fill some fixed-size arrays (only depending on Fourier resolution).
!>
SUBROUTINE fixaray

  USE vmec_main, p5 => cp5
  USE vmec_params, ONLY: jmin2, mscale, nscale, mnyq, nnyq, signgs

  IMPLICIT NONE

  REAL(rprec), PARAMETER :: two=2

  REAL(rprec), PARAMETER :: pexp=4 ! for <M> spectral width screen diagnostic

  INTEGER :: i, m, j, n, mn, mn1, nmin0, istat1, istat2
  INTEGER :: mnyq0, nnyq0
  REAL(rprec):: argi, arg, argj, dnorm

  logical, parameter :: dump_fixaray = .false.
  logical, parameter :: dump_spectral_constraint = .false.


 ! COMPUTE TRIGONOMETRIC FUNCTION ARRAYS
 ! NOTE: ARRAYS ALLOCATED HERE ARE GLOBAL AND ARE DEALLOCATED IN FILEOUT
 ! NOTE: NEED 2 X NYQUIST FOR FAST HESSIAN CALCULATIONS
  mnyq0  = ntheta1/2 ! maximum mode numbers supported by grid
  nnyq0  = nzeta/2

  ! make sure that mnyq, nnyq are at least twice mpol-1, ntor
  ! or large enough to fully represent the information held in realspace (mnyq0, nnyq0)
  mnyq = MAX(0, 2*mnyq0, 2*mpol1)
  nnyq = MAX(0, 2*nnyq0, 2*ntor)

  mnmax_nyq = nnyq/2+1+mnyq*(nnyq+1)/2

  ALLOCATE(cosmu(ntheta3,0:mnyq),  sinmu(ntheta3,0:mnyq),           &
           cosmum(ntheta3,0:mnyq), sinmum(ntheta3,0:mnyq),          &
           cosmui(ntheta3,0:mnyq), cosmumi(ntheta3,0:mnyq),         &
           cosmui3(ntheta3,0:mnyq),cosmumi3(ntheta3,0:mnyq),        &
           sinmui(ntheta3,0:mnyq), sinmumi(ntheta3,0:mnyq),         &
           cosnv(nzeta,0:nnyq),    sinnv(nzeta,0:nnyq),             &
           cosnvn(nzeta,0:nnyq),   sinnvn(nzeta,0:nnyq),            &
           cos01(nznt),            sin01(nznt),            stat=istat1 )
  ALLOCATE(xm(mnmax), xn(mnmax), ixm(mnsize), jmin3(0:mnsize-1),    &
           xm_nyq(mnmax_nyq), xn_nyq(mnmax_nyq),                    &
           mscale(0:mnyq), nscale(0:nnyq), stat=istat2)

  IF (istat1.ne.0) STOP 'allocation error in fixaray: istat1'
  IF (istat2.ne.0) STOP 'allocation error in fixaray: istat2'

  dnorm = one/(nzeta*(ntheta2-1))
  IF (lasym) dnorm = one/(nzeta*ntheta3)     !Fix, SPH012314

  ! (from vmec_params)
  ! array for norming theta-trig functions (internal use only)
  ! so that the discrete SUM[cos(mu)*cos(m'u)] = .5 delta(m,m')
  ! and analogously for zeta/v
  mscale(0) = 1
  nscale(0) = 1
  mscale(1:mnyq) = mscale(0)/osqrt2
  nscale(1:nnyq) = nscale(0)/osqrt2

  r0scale = mscale(0)*nscale(0)

  ! GENERALLY, ONLY NEED THIS FROM 1, ntheta2 EXCEPT IN GETBRHO ROUTINE
  DO i = 1, ntheta3
     argi = twopi*(i-1)/ntheta1
     DO m = 0, mnyq
        arg = argi*m
        cosmu(i,m) = COS(arg)*mscale(m)
        sinmu(i,m) = SIN(arg)*mscale(m)
        cosmui(i,m) = dnorm*cosmu(i,m)
        cosmui3(i,m) = cosmui(i,m)          ! Use this if integration over FULL 1,ntheta3 interval
        sinmui(i,m) = dnorm*sinmu(i,m)
        IF (i.EQ.1 .OR. i.EQ.ntheta2) then
           cosmui(i,m)=cosmui(i,m)/2
        end if
        IF (ntheta2 .EQ. ntheta3) then
           ! cosmui3 was preset from cosmui above,
           ! but in previous check cosmui could have changed,
           ! so update cosmui3 again in case it matters
           ! this is for stellarator symmetry, so lasym==.false.
           cosmui3(i,m) = cosmui(i,m)
        end if
        cosmum(i,m) = cosmu(i,m)*(m)
        sinmum(i,m) =-sinmu(i,m)*(m)
        cosmumi(i,m)= cosmui(i,m)*(m)
        cosmumi3(i,m) = cosmui3(i,m)*m
        sinmumi(i,m)=-sinmui(i,m)*(m)
     END DO
  END DO

  DO j = 1, nzeta
     argj = twopi*(j-1)/nzeta
     DO n = 0, nnyq
        arg = argj*(n)
        cosnv(j,n) = COS(arg)*nscale(n)
        sinnv(j,n) = SIN(arg)*nscale(n)
        cosnvn(j,n) =  cosnv(j,n)*(n*nfp)
        sinnvn(j,n) = -sinnv(j,n)*(n*nfp)
     END DO
  END DO

  ! R,Z,L / s**(m/2) ARE LINEAR NEAR ORIGIN
  mn = 0
  mn1 = 0
  DO m = 0, mpol1
     xmpq(m,1) = m*(m - 1)   ! used for spectral constraint force --> m^2-m

     ! xmpq(m,2:3) are ONLY used for screen diagnostic <M> !!!
     xmpq(m,2) = m**pexp     ! m^p     with p = pexp = 4
     xmpq(m,3) = m**(pexp+1) ! m^(p+q) with q = 1

     ! compute ixm == _i_nteger version of xm
     DO n = 0, ntor
        jmin3(mn) = jmin2(m) ! vmec_params: starting js(m) values for which R,Z are evolved
        mn = mn + 1
        ixm(mn) = m
     END DO

     ! use this loop also to compute xm, xn array contents
     nmin0 = -ntor
     IF (m .eq. 0) nmin0 = 0
     DO n = nmin0, ntor
        mn1 = mn1 + 1
        xm(mn1) = m
        xn(mn1) = n*nfp
     END DO
  END DO

  IF (mn1 .ne. mnmax) STOP 'mn1 != mnmax'

  ! COMPUTE NYQUIST-SIZED ARRAYS FOR OUTPUT.
  ! RESTORE m,n Nyquist TO 1 X ... (USED IN WROUT, JXBFORCE)
  !  mnyq = mnyq0;  nnyq = nnyq0
  mnyq = mnyq/2
  nnyq = nnyq/2

  mn1 = 0
  DO m = 0, mnyq
     nmin0 = -nnyq
     IF (m .eq. 0) nmin0 = 0
     DO n = nmin0, nnyq
        mn1 = mn1 + 1
        xm_nyq(mn1) = m
        xn_nyq(mn1) = n*nfp
     END DO
  END DO

  IF (mn1 .ne. mnmax_nyq) STOP 'mn1 != mnmax_nyq'

  ! cos01 and sin01 are used as trigmult for precondn in bcovar
  ! cos01:  m*cos(m u - n v) for m=1, n=0 --> d(sin(mu))/du
  ! sin01: -m*sin(m u - n v) for m=1, n=0 --> d(cos(mu))/du
  ! evaluated on (ntheta3 x nzeta) grid
  mn = 0 ! mn is re-used as linear grid index here
  m = 1
  DO i = 1, ntheta3
     argi = twopi*(i - 1)/ntheta1
     DO j = 1, nzeta
        mn = mn + 1
        cos01(mn) = m*COS(m*argi)*mscale(m)
        sin01(mn) =-m*SIN(m*argi)*mscale(m)
     END DO
  END DO

  ! _fac_tor for _con_straint
  faccon(0) = zero
  faccon(1:mpol1-1) = -0.25_dp*signgs/xmpq(2:mpol1,1)**2
  faccon(mpol1) = zero

  if (dump_fixaray) then
    open(unit=42, file="fixaray."//trim(input_extension), status="unknown")

    write(42, *) "# ntheta3 mnyq nzeta nnyq nznt mnmax mnsize mnmax_nyq"
    write(42, *) ntheta3, mnyq, nzeta, nnyq, nznt, mnmax, mnsize, mnmax_nyq

    write(42, *) "# i m cosmu sinmu cosmum sinmum cosmui sinmui cosmui3 cosmumi sinmumi cosmumi3"
    DO i = 1, ntheta3
      DO m = 0, mnyq
        write (42, *) i, m, cosmu(i,m), sinmu(i,m), cosmum(i,m), sinmum(i,m), &
                      cosmui(i,m), sinmui(i,m), cosmui3(i,m), &
                      cosmumi(i,m), sinmumi(i,m), cosmumi3(i,m)
      end do
    end do

    write(42, *) "# j n cosnv sinnv cosnvn sinnvn"
    DO j = 1, nzeta
      DO n = 0, nnyq
        write(42, *) j, n, cosnv(j,n), sinnv(j,n), cosnvn(j,n), sinnvn(j,n)
      END DO
    END DO

    write(42, *) "# i j mn cos01 sin01"
    mn = 0 ! mn is re-used as linear grid index here
    DO i = 1, ntheta3
      DO j = 1, nzeta
        mn = mn + 1
        write(42, *) i, j, mn, cos01(mn), sin01(mn)
      END DO
    END DO

    write(42, *) "# mn1 xm xn"
    DO mn1 = 1, mnmax
      write(42, *) mn1, xm(mn1), xn(mn1)
    END DO

    write(42, *) "# m n mn ixm"
    mn = 0
    DO m = 0, mpol1
      DO n = 0, ntor
        mn = mn + 1
        write(42, *) m, n, mn, ixm(mn)
      END DO
    end do

    write(42, *) "# mn1 xm_nyq xn_nyq"
    DO mn1 = 1, mnmax_nyq
        write(42, *) mn1, xm_nyq(mn1), xn_nyq(mn1)
     END DO

    write(42, *) "# m mscale"
    DO m = 0, mnyq
      write(42, *) m, mscale(m)
    end do

    write(42, *) "# n nscale"
    DO n = 0, nnyq
      write(42, *) n, nscale(n)
    end do

    close(42)
    stop "fixaray output dumped to fixaray.<ext>"
  end if

  if (dump_spectral_constraint) then
    open(unit=42, file="spectral_constraint."//trim(input_extension), status="unknown")

    write(42, *) "# mpol1"
    write(42, *) mpol1

    write(42, *) "# m xmpq(m,1) xmpq(m,2) xmpq(m,3) faccon"
    do m = 0, mpol1
      write(42, *) m, xmpq(m,1), xmpq(m,2), xmpq(m,3), faccon(m)
    end do

    close(42)
    stop "spectral constraint output dumped to spectral_constraint.<ext>"
  end if

END SUBROUTINE fixaray
