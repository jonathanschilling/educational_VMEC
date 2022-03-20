!> \file
!> \brief allocate and fill some fixed-size arrays (only depending on Fourier resolution).

!> \brief allocate and fill some fixed-size arrays (only depending on Fourier resolution).
!>
SUBROUTINE fixaray

  USE vmec_main, p5 => cp5
  USE vmec_params, ONLY: jmin2, mscale, nscale, mnyq, nnyq, signgs

  use dbgout

  IMPLICIT NONE

  REAL(rprec), PARAMETER :: two=2

  REAL(rprec), PARAMETER :: pexp=4 ! for <M> spectral width screen diagnostic

  INTEGER :: i, m, j, n, mn, mn1, nmin0, istat1, istat2
  INTEGER :: mnyq0, nnyq0
  REAL(rprec):: argi, arg, argj, dnorm, dnorm3
  logical :: dbg_fixaray

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

  ! Normalization factor for forward Fourier transforms (e.g. tomnsp*)
  ! In any case (symmetric or asymmetric), the Fourier integrals in tomnsp* and alias
  ! (the only place where dnorm is used through cosmui, sinmui, etc.)
  ! are ever only taken over [0, pi] corresponding to 1, ..., ntheta2.
  dnorm = one/(nzeta*(ntheta2-1))

  ! Normalization factor for surface integrals/averages (wint is based on cosmui3).
  ! For the asymmetric case, the norm in wint thus has to be 1/(nzeta*ntheta3).
  IF (lasym) then
    dnorm3 = one/(nzeta*ntheta1)
  else
    dnorm3 = one/(nzeta*(ntheta2-1))
  end if

  ! (from vmec_params)
  ! array for norming theta-trig functions (internal use only)
  ! so that the discrete SUM[cos(mu)*cos(m'u)] = .5 delta(m,m')
  ! and analogously for zeta/v
  mscale(0) = 1
  nscale(0) = 1
  mscale(1:mnyq) = mscale(0)/osqrt2 ! == sqrt(2)
  nscale(1:nnyq) = nscale(0)/osqrt2 ! == sqrt(2)

  r0scale = mscale(0)*nscale(0) ! == 1

  ! GENERALLY, ONLY NEED THIS FROM 1, ntheta2 EXCEPT IN GETBRHO ROUTINE
  DO i = 1, ntheta3
     argi = twopi*(i-1)/ntheta1
     DO m = 0, mnyq
        arg = argi*m
        cosmu(i,m) = COS(arg)*mscale(m)
        sinmu(i,m) = SIN(arg)*mscale(m)

        cosmui(i,m) = dnorm*cosmu(i,m)
        sinmui(i,m) = dnorm*sinmu(i,m)
        IF (i.EQ.1 .OR. i.EQ.ntheta2) then
           ! Trapezoidal integration requires a factor of 1/2 for the first and the last point.
           ! Note that this is done also in the case of an asymmetric run!
           ! There, cosmui is only used up to ntheta2 anyway...
           cosmui(i,m)=cosmui(i,m)/2.0_dp
        end if

        ! Use this if integration over FULL 1,ntheta3 interval
        cosmui3(i,m) = dnorm3*cosmu(i,m)
        IF (.not.lasym .and. (i.eq.1 .or. i.eq.ntheta2)) then
           ! Note that cosmui3 is only ever used to construct the wint array.
           ! where only the m=0 component of cosmui3 enters.
           cosmui3(i,m) = cosmui3(i,m)/2.0_dp
        end if

        cosmum(i,m)   = cosmu(i,m)*(m)
        sinmum(i,m)   =-sinmu(i,m)*(m)
        cosmumi(i,m)  = cosmui(i,m)*(m)
        sinmumi(i,m)  =-sinmui(i,m)*(m)

        cosmumi3(i,m) = cosmui3(i,m)*m ! not used... --> remove?
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

  dbg_fixaray = open_dbg_context("fixaray")
  if (dbg_fixaray) then

    call add_real_2d("cosmu",    ntheta2, mnyq+1, cosmu(1:ntheta2,:))
    call add_real_2d("sinmu",    ntheta2, mnyq+1, sinmu(1:ntheta2,:))
    call add_real_2d("cosmum",   ntheta2, mnyq+1, cosmum(1:ntheta2,:))
    call add_real_2d("sinmum",   ntheta2, mnyq+1, sinmum(1:ntheta2,:))
    call add_real_2d("cosmui",   ntheta2, mnyq+1, cosmui(1:ntheta2,:))
    call add_real_2d("sinmui",   ntheta2, mnyq+1, sinmui(1:ntheta2,:))
    call add_real_2d("cosmui3",  ntheta2, mnyq+1, cosmui3(1:ntheta2,:))
    call add_real_2d("cosmumi",  ntheta2, mnyq+1, cosmumi(1:ntheta2,:))
    call add_real_2d("sinmumi",  ntheta2, mnyq+1, sinmumi(1:ntheta2,:))
    call add_real_2d("cosmumi3", ntheta2, mnyq+1, cosmumi3(1:ntheta2,:))

    call add_real_2d("cosnv",  nzeta, nnyq+1, cosnv)
    call add_real_2d("sinnv",  nzeta, nnyq+1, sinnv)
    call add_real_2d("cosnvn", nzeta, nnyq+1, cosnvn)
    call add_real_2d("sinnvn", nzeta, nnyq+1, sinnvn)

    call add_real_1d("mscale", mnyq+1, mscale)
    call add_real_1d("nscale", nnyq+1, nscale)

    ! TODO: add r0scale output (currently in residue...)

  end if

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

  if (dbg_fixaray) then

    call add_int("ntheta3", ntheta3)
    call add_int("mnyq", mnyq)
    call add_int("nzeta", nzeta)
    call add_int("nnyq", nnyq)
    call add_int("nznt", nznt)
    call add_int("mnmax", mnmax)
    call add_int("mnsize", mnsize)
    call add_int("mnmax_nyq", mnmax_nyq)

    call add_real_1d("cos01", nznt, cos01)
    call add_real_1d("sin01", nznt, sin01)

    call add_real_1d("xm", mnmax, xm)
    call add_real_1d("xn", mnmax, xn)

    call add_real_1d("xm_nyq", mnmax_nyq, xm_nyq)
    call add_real_1d("xn_nyq", mnmax_nyq, xn_nyq)

    call add_int_1d("ixm", mnsize, ixm)

    call close_dbg_out()
  end if

  if (open_dbg_context("spectral_constraint")) then

    ! xmpq is allocated statically, so need size here explicitly!
    call add_real_2d("xmpq", mpol, 3, xmpq(0:mpol1,1:3))

    call add_real_1d("faccon", mpol, faccon)

    call close_dbg_out()
  end if

END SUBROUTINE fixaray
