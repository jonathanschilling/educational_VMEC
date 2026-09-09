!> \file
!> \brief Compute the analytical-and-numerical 4D Fourier integrals over the equivalently-singular functions.

!> \brief Compute the analytical-and-numerical 4D Fourier integrals over the equivalently-singular functions.
!>
!> @param grpmn
!> @param bvec
!> @param ivacskip
!> @param lasym
!> @param m_map
!> @param n_map
!> @param grpmn_m_map
!> @param grpmn_n_map
SUBROUTINE analyt(grpmn, bvec, ivacskip, lasym, m_map, n_map, grpmn_m_map, grpmn_n_map)
  USE vacmod, vm_grpmn => grpmn
  use dbgout
  use vmec_main, only: num_eqsolve_retries
  IMPLICIT NONE

  INTEGER, INTENT(in) :: ivacskip
  REAL(rprec), INTENT(out) :: grpmn(nuv2*mnpd2)
  REAL(rprec), INTENT(out) :: bvec(mnpd2)
  real(rprec), intent(out) :: m_map(mnpd2)
  real(rprec), intent(out) :: n_map(mnpd2)
  real(rprec), intent(out) :: grpmn_m_map(nuv2*mnpd2)
  real(rprec), intent(out) :: grpmn_n_map(nuv2*mnpd2)
  integer, intent(in) :: lasym

  INTEGER :: l, n, m, i
  REAL(rprec) :: fl, fl1, sign1

  ! ALL EQUATIONS REFER TO THE PAPER BY P. MERKEL (PKM)
  ! IN J. COMPUT. PHYSICS 66, p83 (1986)
  !
  ! IN GOING BETWEEN THE COMPLEX NOTATION OF (PKM) AND OUR REAL FORM,
  ! NOTE THAT THE INTEGRALS (APPENDIX, PKM) Imn AND Kmn ARE BOTH REAL.
  ! THUS, THE SIN(mu-nv) INTEGRALS OF THE SINGULAR PIECE (ANALYTIC CONTRIBUTION)
  ! VANISHES.
  !
  ! THE REQUIRED SOURCE-TERM INTEGRALS ARE (Eq.2.16-2.17):
  !
  ! BVECS(m,n) = Int< SIN(mu' - nv') han(u',v') >
  ! BVECC(m,n) = Int< COS(mu' - nv') han(u',v') >
  !
  ! Where Int<...> means integration over u (theta) and v (zeta) and
  ! summation over field periods. These can be written in terms of PKM integrals
  ! Imn(a,b,c), where a(u,v) = guu (g theta-theta), etc.:
  !
  ! BVECS(m,n) = ALP * Int<SIN(mu' - nv') * F * Im,-n(a,b,c)>
  ! BVECC(m,n) = ALP * Int<COS(mu' - nv') * F * Im,-n(a,b,c)>
  !
  ! Here, F = - BNORM(u',v') is defined in Eq.(2.13), and ALP = (2*pi/nfp).
  !
  ! Similarly, the analytic part of the matrix A(m,n;m',n') can be written:
  !
  ! A(m,n;m',n') = (2*pi/nfp) * Int<SIN(mu' - nv')*SIN(m'u' - n'v')
  !                          [Km,-n](a',b',c';A',B',C')>
  !
  ! On EXIT, GRPMN(ip,m,n) = ALP * SIN(ip,m,n) * K[m,-n](ip)
  !
  ! COMPUTE ALL QUANTITIES INDEPENDENT OF THE MODE INDICES L,M,N
  ! NOTE: 2b = guv_b HAS FACTOR OF 2 BUILT IN (see SUBROUTINE SURFACE)
  !
  ! ADP(M): a +(-)2b + c
  ! CMA:    c - a
  ! DELTA:  4*(ac - b**2)
  ! AZP(M): A +(-)2*B + C
  ! CMA1:   C - A
  ! R1P(M): Coefficient of l*Tl+(-) in eq (A17)
  ! R0P(M): Coefficient of l*T(l-1)+(-) in eq (A17)
  ! RA1P(M):Coefficient of Tl+(-) in eq (A17)

  adp  = guu_b  + guv_b  + gvv_b
  adm  = guu_b  - guv_b  + gvv_b
  cma  = gvv_b  - guu_b
  sqrtc  = two*SQRT(gvv_b)
  sqrta  = two*SQRT(guu_b)
  sqad1u = SQRT(adp)
  sqad2u = SQRT(adm)

  ! INITIALIZE VECTORS
  bvec = 0
  IF (ivacskip .eq. 0) THEN
     grpmn = 0

     delt1u  = adp*adm  - cma*cma
     azp1u  = auu  + auv  + avv
     azm1u  = auu  - auv  + avv
     cma11u  = avv  - auu
     r1p  = (azp1u*(delt1u - cma*cma)/adp - azm1u*adp + two*cma11u*cma)/delt1u
     r1m  = (azm1u*(delt1u - cma*cma)/adm - azp1u*adm + two*cma11u*cma)/delt1u
     r0p  = (-azp1u*adm*cma/adp - azm1u*cma + two*cma11u*adm)/delt1u
     r0m  = (-azm1u*adp*cma/adm - azp1u*cma + two*cma11u*adp)/delt1u
     ra1p = azp1u/adp
     ra1m = azm1u/adm
  ENDIF

  ! INITIALIZE T0+ and T0-
  !
  ! TLP(M): TL+(-)
  ! TLP(M)1:T(L-1)+(-)
  ! TLP(M)2:T(L-2)+(-)
  tlp1 = 0
  tlm1 = 0
  tlp  = one/sqad1u*log((sqad1u*sqrtc + adp + cma)/(sqad1u*sqrta - adp + cma))
  tlm  = one/sqad2u*log((sqad2u*sqrtc + adm + cma)/(sqad2u*sqrta - adm + cma))
  tlpm = tlp  + tlm

  ! TABULATE T+ AND T- FOR EVERY L BEFORE THE L-SUM
  !
  ! The two homogeneous solutions of Eq (A15) form a pair of modulus
  ! sqrt(adm/adp) for T+, and its reciprocal for T-, so stepping the
  ! recurrence upward over mf + nf steps multiplies the rounding of its
  ! inputs by that factor. Where the growth stays small the upward pass is
  ! used; where it does not, the same recurrence is run downward from a seed
  ! far enough above mf + nf that the seed has decayed out of the result by
  ! the time the pass reaches it.
  do i = 1, nuv2
     call recurrence_sum(adm(i), adp(i), sqrtc(i), sqrta(i), cma(i), tlp(i), all_tlp(:,i))
     call recurrence_sum(adp(i), adm(i), sqrtc(i), sqrta(i), cma(i), tlm(i), all_tlm(:,i))
  end do

  ! BEGIN L-SUM IN EQ (A14) TO COMPUTE Imn (and Kmn) INTEGRALS
  ! NOTE THAT IN THE LOOP OVER L BELOW: L == |m - n| + 2L_A14
  ! THUS, L BELOW IS THE INDEX OF THE T+- (S+-)
  sign1 = 1
  fl1 = 0

  LLOOP: DO l = 0, mf + nf
     fl = fl1

     ! T+- for this l, and for the one below it
     tlp = all_tlp(l,:)
     tlm = all_tlm(l,:)
     if (l .eq. 0) then
        tlp1 = 0
        tlm1 = 0
     else
        tlp1 = all_tlp(l-1,:)
        tlm1 = all_tlm(l-1,:)
     end if
     tlpm = tlp + tlm

     ! COMPUTE SL+ and SL- , Eq (A17)
     ! SLP(M): SL+(-)
     IF (ivacskip .eq. 0) THEN
        slp = (r1p*fl + ra1p)*tlp + r0p*fl*tlp1 - (r1p + r0p)/sqrtc + sign1*(r0p - r1p)/sqrta
        slm = (r1m*fl + ra1m)*tlm + r0m*fl*tlm1 - (r1m + r0m)/sqrtc + sign1*(r0m - r1m)/sqrta
        slpm = slp + slm

        ! here, slp/m are available for the current value of l
        ! --> save into matrix for debugging
        all_slp(l,:) = slp
        all_slm(l,:) = slm
     ENDIF

     ! BEGIN MODE NUMBER (m,n) LOOP
     DO n = 0, nf
        DO m = 0, mf
           ! always call analysum, analysum2 to get bvec fill pattern (m_map, n_map)
           ! IF (cmns(l,m,n) .eq. zero) CYCLE

           IF (n.eq.0 .or. m.eq.0) THEN
             ! 1. n = 0 and  m >= 0  OR n > 0 and m = 0
             CALL analysum (grpmn, bvec, slpm, tlpm, m, n, l, ivacskip, lasym, m_map, n_map, &
                            grpmn_m_map, grpmn_n_map)
           ELSE
             ! 2. n>=1  and  m>=1
             CALL analysum2 (grpmn, bvec, m, n, l, ivacskip, lasym, m_map, n_map, &
                             grpmn_m_map, grpmn_n_map)
           ENDIF
        END DO
     END DO

     ! next l
     fl1  = fl1 + 1

     ! (-1)**l (next l now)
     sign1 = -sign1
  END DO LLOOP

  if (open_dbg_context("vac1n_analyt", num_eqsolve_retries)) then

    call add_real_3d("all_tlp", mf+nf+1, nv, nu3, all_tlp)
    call add_real_3d("all_tlm", mf+nf+1, nv, nu3, all_tlm)
    call add_real_2d("bvec", mf1, nf1, bvec) ! (mf+1)x(2*nf+1)x(ndim: 1 or 2)

    if (ivacskip .eq. 0) then
      call add_real_3d("all_slp", mf+nf+1, nv, nu3, all_slp)
      call add_real_3d("all_slm", mf+nf+1, nv, nu3, all_slm)
      call add_real_4d("grpmn", mf1, nf1, nv, nu3, grpmn) ! missing dim: (ndim: 1 or 2)
    else
      call add_null("all_slp")
      call add_null("all_slm")
      call add_null("grpmn")
    end if

    call close_dbg_out()
  end if

CONTAINS

!> \brief One step of the T+- recurrence, Eq (A15).
!>
!> @param a     ad of the opposite parity to the one being stepped
!> @param b     ad of the parity being stepped
!> @param sqrtc 2*sqrt(gvv_b)
!> @param sqrta 2*sqrt(guu_b)
!> @param sgn   (-1)**l of the step
!> @param fl    l
!> @param fl1   l + 1
!> @param fl2   2*l + 1
!> @param cma   gvv_b - guu_b
!> @param below T(l-1)
!> @param here  T(l)
function recurrence(a, b, sqrtc, sqrta, sgn, fl, fl1, fl2, cma, below, here)
  implicit none
  real(rprec) :: recurrence
  real(rprec), intent(in) :: a, b, sqrtc, sqrta, sgn, cma, below, here
  integer, intent(in) :: fl, fl1, fl2

  recurrence = (sqrtc + sgn*sqrta - fl2*cma*here - fl*a*below)/(b*fl1)
end function recurrence

!> \brief Tabulate T+ or T- for l = 0 .. mf + nf at one surface point.
!>
!> @param a  ad of the opposite parity to the one being tabulated
!> @param b  ad of the parity being tabulated
!> @param sqrtc 2*sqrt(gvv_b)
!> @param sqrta 2*sqrt(guu_b)
!> @param cma   gvv_b - guu_b
!> @param t0    T(0), known in closed form
!> @param tl    the tabulated sequence
subroutine recurrence_sum(a, b, sqrtc, sqrta, cma, t0, tl)
  implicit none
  real(rprec), intent(in) :: a, b, sqrtc, sqrta, cma, t0
  real(rprec), intent(inout) :: tl(0:mf+nf)

  ! Upward stepping multiplies the rounding of its inputs by
  ! sqrt(a/b) per step; ten is the growth this tolerates over mf + nf steps.
  real(rprec), parameter :: max_forward_log_growth = 2.302585092994046_rprec
  ! A downward pass seeded with zero leaves (b/a)**((top - l)/2) of T(top) in
  ! T(l); the seed sits far enough above mf + nf for that to fall below 1e-17.
  real(rprec), parameter :: min_seed_log_decay = 39.14394658089878_rprec

  real(rprec) :: high, current, low, sgn
  integer :: l, ntail

  if (a .gt. b .and. b .gt. zero .and. (mf + nf)*log(a/b) .gt. two*max_forward_log_growth) then
     ntail = ceiling(two*min_seed_log_decay/log(a/b))
     high = 0
     current = 0
     if (mod(mf + nf + ntail, 2) .eq. 0) then
        sgn = -1
     else
        sgn = 1
     end if
     do l = mf + nf + ntail, mf + nf + 2, -1
        low = recurrence(b, a, sqrtc, sqrta, sgn, l + 1, l, 2*l + 1, cma, high, current)
        high = current
        current = low
        sgn = -sgn
     end do
     do l = mf + nf + 1, 1, -1
        tl(l - 1) = recurrence(b, a, sqrtc, sqrta, sgn, l + 1, l, 2*l + 1, cma, high, current)
        high = current
        current = tl(l - 1)
        sgn = -sgn
     end do
     ! the downward pass fixes the sequence only up to a factor
     tl(0) = t0
  else
     tl(0) = t0
     sgn = 1
     do l = 0, mf + nf - 1
        sgn = -sgn
        if (l .eq. 0) then
           tl(l + 1) = recurrence(a, b, sqrtc, sqrta, sgn, l, l + 1, 2*l + 1, cma, zero, tl(l))
        else
           tl(l + 1) = recurrence(a, b, sqrtc, sqrta, sgn, l, l + 1, 2*l + 1, cma, tl(l - 1), tl(l))
        end if
     end do
  end if
end subroutine recurrence_sum

END SUBROUTINE analyt
