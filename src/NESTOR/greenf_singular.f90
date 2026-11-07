!> \file
!> \brief Correction of the punctured trapezoidal sums of greenf at the source point.

!> \brief Correction of the punctured trapezoidal sums of greenf at the source point.
!>
!> About the source point the remainders G - ga1 and dG/dn - ga1*ga2 expand as
!> f0 + f1 + ... in the steps (s, t) of u and v, with f0 odd and homogeneous of
!> degree 0 and f1 even and homogeneous of degree 1. The punctured trapezoidal
!> sum of a remainder times a smooth function b differs from the integral by
!> L[f1] b + L[f0 s] db/du + L[f0 t] db/dv at third order in the grid step,
!> where L is the regularized lattice sum. In coordinates in which the metric
!> is the identity the six functions are |x| times even harmonics of the polar
!> angle up to order 8, so L follows from Epstein zeta functions with harmonic
!> weights, which are summed by Ewald's method.
module greenf_singular
  use stel_kinds, only: dp
  implicit none
  private
  public :: singular_defects
  real(dp), parameter :: pi = 3.14159265358979323846264338_dp
  ! the six functions are trigonometric polynomials of degree 8 with period
  ! pi in the polar angle
  integer, parameter :: nang = 10

contains

!> \brief Product of two polynomials in (s, t) truncated at degree 4; c(i,j) multiplies s**i*t**j.
pure function tmul(a, b) result(c)
  real(dp), intent(in) :: a(0:4,0:4), b(0:4,0:4)
  real(dp) :: c(0:4,0:4)
  integer :: i, j, k, l
  c = 0
  do i = 0, 4
    do j = 0, 4 - i
      if (a(i,j) .eq. 0) cycle
      do k = 0, 4 - i - j
        do l = 0, 4 - i - j - k
          c(i+k,j+l) = c(i+k,j+l) + a(i,j)*b(k,l)
        end do
      end do
    end do
  end do
end function tmul

!> \brief Homogeneous part of degree d, from the powers sp of s and tq of t.
pure real(dp) function tpart(a, d, sp, tq)
  real(dp), intent(in) :: a(0:4,0:4), sp(0:4), tq(0:4)
  integer, intent(in) :: d
  integer :: i
  tpart = 0
  do i = 0, d
    tpart = tpart + a(i,d-i)*sp(i)*tq(d-i)
  end do
end function tpart

!> \brief Regularized lattice sums with harmonic weights.
!>
!> esum(q) = A*Z_k, k = 2*q, q = 0, ..., 4: the regularized sum of
!> A*|x|**(1-k)*(x1 + i*x2)**k over the lattice spanned by e1 and e2, of cell
!> area A, less its integral. Z_k is the Epstein zeta function with the
!> harmonic weight at s = (k - 1)/2, split by Ewald's method between the
!> lattice and its dual; the terms fall off as exp(-pi*|x|**2) in the lattice
!> of unit cell area, and x and -x contribute alike.
subroutine harmonic_sums(e1, e2, esum)
  real(dp), intent(in) :: e1(2), e2(2)
  complex(dp), intent(out) :: esum(0:4)
  real(dp), parameter :: reach2 = 7
  real(dp) :: scale, f1(2), f2(2), d1(2), d2(2), b1(2), b2(2), w(2)
  real(dp) :: det, area, mu, g11, g22, g12, gdet, r2, x, xinv
  real(dp) :: ex, sx, a, xa, xp, g(0:6), pgam(0:4)
  complex(dp) :: direct(0:4), dual(0:4), z2, zk, total
  integer :: it, side, i, j, ni, nj, q, m

  ! unit cell area, then the Lagrange-Gauss reduced basis
  scale = sqrt(abs(e1(1)*e2(2) - e1(2)*e2(1)))
  f1 = e1/scale
  f2 = e2/scale
  do it = 1, 100
    if (sum(f1*f1) .gt. sum(f2*f2)) then
      w = f1
      f1 = f2
      f2 = w
      cycle
    end if
    mu = anint(sum(f1*f2)/sum(f1*f1))
    if (mu .eq. 0) exit
    f2 = f2 - mu*f1
  end do
  det = f1(1)*f2(2) - f1(2)*f2(1)
  d1 = [f2(2), -f2(1)]/det
  d2 = [-f1(2), f1(1)]/det
  area = abs(det)

  direct = 0
  dual = 0
  do side = 1, 2
    if (side .eq. 1) then
      b1 = f1
      b2 = f2
    else
      b1 = d1
      b2 = d2
    end if
    g11 = sum(b1*b1)
    g22 = sum(b2*b2)
    g12 = sum(b1*b2)
    gdet = g11*g22 - g12*g12
    ni = ceiling(sqrt(reach2*g22/gdet))
    nj = ceiling(sqrt(reach2*g11/gdet))
    do i = 0, ni
      do j = -nj, nj
        if (i .eq. 0 .and. j .le. 0) cycle
        w = i*b1 + j*b2
        r2 = sum(w*w)
        if (r2 .gt. reach2) cycle
        x = pi*r2
        xinv = 1/x
        ! g(m) = GAMMA(m - 1/2, x) by erfc and the recurrence
        ex = exp(-x)
        sx = sqrt(x)
        g(1) = sqrt(pi)*erfc(sx)
        g(0) = 2*(ex/sx - g(1))
        a = 0.5_dp
        xa = sx
        do m = 2, 6
          g(m) = a*g(m-1) + xa*ex
          a = a + 1
          xa = xa*x
        end do
        z2 = cmplx(w(1), w(2), dp)**2
        zk = 1
        if (side .eq. 1) then
          ! P_k x**(-s) GAMMA(s, x)
          xp = sx
          do q = 0, 4
            direct(q) = direct(q) + zk*xp*g(q)
            zk = zk*z2
            xp = xp*xinv
          end do
        else
          ! P_k x**(-(k + 3)/2) GAMMA((k + 3)/2, x)
          xp = sx*xinv*xinv
          do q = 0, 4
            dual(q) = dual(q) + zk*xp*g(q+2)
            zk = zk*z2
            xp = xp*xinv
          end do
        end if
      end do
    end do
  end do

  ! pi**s/GAMMA(s) at s = -1/2, 1/2, ..., 7/2
  pgam(0) = -1/(2*pi)
  do q = 1, 4
    pgam(q) = pi/(q - 1.5_dp)*pgam(q-1)
  end do
  do q = 0, 4
    ! (-i)**k = (-1)**q; the sums ran over half of each lattice
    total = 2*(direct(q) + (-1)**q*dual(q)/area)
    if (q .eq. 0) total = total - 2/(3*area) + 2
    esum(q) = area*pgam(q)*total*scale**3
  end do
end subroutine harmonic_sums

!> \brief Defects of the punctured trapezoidal sums at the source point.
!>
!> The defects are minus the regularized lattice sums, in the units of the
!> integral over (s, t): wg(1:3) for G - ga1 with the weights 1, s and t,
!> wk(1:3) the same for dG/dn - ga1*ga2.
!>
!> @param rt Taylor coefficients of R in the steps (s, t) of u and v
!> @param zt Taylor coefficients of Z in the steps (s, t) of u and v
!> @param x0 Cartesian x of the source
!> @param y0 Cartesian y of the source
!> @param z0 Cartesian z of the source
!> @param phi0 toroidal angle of the source
!> @param dphi rate of the toroidal angle in v
!> @param xn Cartesian x of the normal at the source
!> @param yn Cartesian y of the normal at the source
!> @param zn Cartesian z of the normal at the source
!> @param guu metric element guu_b at the source
!> @param guv metric element guv_b at the source
!> @param gvv metric element gvv_b at the source
!> @param auu second-form element auu at the source
!> @param auv second-form element auv at the source
!> @param avv second-form element avv at the source
!> @param hu grid step in u
!> @param hv grid step in v
!> @param wg defects of the single layer
!> @param wk defects of the double layer
subroutine singular_defects(rt, zt, x0, y0, z0, phi0, dphi, xn, yn, zn, &
                            guu, guv, gvv, auu, auv, avv, hu, hv, wg, wk)
  real(dp), intent(in) :: rt(0:4,0:4), zt(0:4,0:4), x0, y0, z0, phi0, dphi, &
                          xn, yn, zn, guu, guv, gvv, auu, auv, avv, hu, hv
  real(dp), intent(out) :: wg(3), wk(3)
  real(dp), parameter :: fact(0:4) = [1, 1, 2, 6, 24]
  real(dp), dimension(0:4,0:4) :: cphi, sphi, xc, yc, wc, r2, nr
  real(dp) :: rate, ca, cb, cc, th, s, t, sp(0:4), tq(0:4), ri, c3, d4
  real(dp) :: a2, a3, a4, s3t, t4, s4, i32, i52, i72, f0g, f1g, f0k
  real(dp) :: f1k, phi(6), e1(2), e2(2), el, cq(0:3)
  complex(dp) :: esum(0:4), coef(6,0:4), e2th, ek
  integer :: j, n, k, m

  ! Cartesian position relative to the source and the normal component, as
  ! polynomials in (s, t); cq holds cos(phi0 + j*pi/2)
  cq = [cos(phi0), -sin(phi0), -cos(phi0), sin(phi0)]
  cphi = 0
  sphi = 0
  rate = 1
  do j = 0, 4
    cphi(0,j) = rate/fact(j)*cq(mod(j, 4))
    sphi(0,j) = rate/fact(j)*cq(mod(j + 3, 4))
    rate = rate*dphi
  end do
  xc = tmul(rt, cphi)
  yc = tmul(rt, sphi)
  wc = zt
  xc(0,0) = xc(0,0) - x0
  yc(0,0) = yc(0,0) - y0
  wc(0,0) = wc(0,0) - z0
  r2 = tmul(xc, xc) + tmul(yc, yc) + tmul(wc, wc)
  nr = xn*xc + yn*yc + zn*wc

  ! (ca*s + cb*t, cc*t) are coordinates in which the metric is the identity
  ca = sqrt(guu)
  cb = guv/(2*ca)
  cc = sqrt(gvv - cb*cb)

  ! harmonics of the six functions on the unit circle of those coordinates
  coef = 0
  do n = 0, nang - 1
    th = pi*n/nang
    t = sin(th)/cc
    s = (cos(th) - cb*t)/ca
    sp(0) = 1
    tq(0) = 1
    do m = 1, 4
      sp(m) = sp(m-1)*s
      tq(m) = tq(m-1)*t
    end do
    ri = tpart(r2, 2, sp, tq)
    c3 = tpart(r2, 3, sp, tq)
    d4 = tpart(r2, 4, sp, tq)
    a2 = tpart(nr, 2, sp, tq)
    a3 = tpart(nr, 3, sp, tq)
    a4 = tpart(nr, 4, sp, tq)
    ! quartic parts of the quadratic forms of ga1 and ga2 in 2*tan(s/2) and
    ! 2*tan(t/2)
    s3t = sp(3)*t + s*tq(3)
    t4 = (guu*sp(4) + gvv*tq(4))/6 + guv*s3t/12
    s4 = (auu*sp(4) + avv*tq(4))/6 + auv*s3t/12
    i32 = 1/(ri*sqrt(ri))
    i52 = i32/ri
    i72 = i52/ri
    f0g = -c3*i32/2
    f1g = -(d4 - t4)*i32/2 + 3*c3*c3*i52/8
    f0k = a3*i32 - 3*a2*c3*i52/2
    f1k = (a4 - s4)*i32 - 3*(a3*c3 + a2*(d4 - t4))*i52/2 + 15*a2*c3*c3*i72/8
    phi = [f1g, f0g*s, f0g*t, f1k, f0k*s, f0k*t]/nang
    e2th = cmplx(cos(2*th), -sin(2*th), dp)
    ek = 1
    do k = 0, 4
      coef(:,k) = coef(:,k) + phi*ek
      ek = ek*e2th
    end do
  end do

  e1 = [ca*hu, 0.0_dp]
  e2 = [cb*hv, cc*hv]
  call harmonic_sums(e1, e2, esum)
  do m = 1, 6
    el = real(coef(m,0)*esum(0), dp)
    do k = 1, 4
      el = el + 2*real(coef(m,k)*esum(k), dp)
    end do
    el = -el/(ca*cc)
    if (m .le. 3) then
      wg(m) = el
    else
      wk(m-3) = el
    end if
  end do
end subroutine singular_defects

end module greenf_singular
