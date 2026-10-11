!> \file
!> \brief Extrema of the cross-section of a flux surface in one toroidal plane.

!> \brief Extrema of the cross-section of a flux surface in one toroidal plane,
!>        used for the geometry table in \c eqfor.
MODULE cross_section

  USE stel_kinds, ONLY: rprec

  IMPLICIT NONE

CONTAINS

!> \brief R and Z extrema of the cross-section of surface \c js in toroidal
!>        plane \c noff, and R at the point of largest Z.
!>
!> A lasym run stores the complete contour of the plane, which is scanned in
!> one pass. Under stellarator symmetry only theta in [0, pi] is stored, and the
!> reflected plane 2 pi - v with Z -> -Z supplies the other half.
!>
!> @param r1 R in real space, even (0) and odd (1) parity in m
!> @param z1 Z in real space, even (0) and odd (1) parity in m
!> @param sqrts_js \f$\sqrt{s}\f$ on surface \c js
!> @param js index of the surface
!> @param noff index of the toroidal plane
!> @param ns number of surfaces
!> @param nzeta number of toroidal grid points
!> @param ntheta2 number of poloidal grid points in [0, pi]
!> @param ntheta3 number of poloidal grid points stored per surface and plane
!> @param lasym whether the run is non-stellarator-symmetric
!> @param xmin smallest R
!> @param xmax largest R
!> @param zmin smallest Z
!> @param zmax largest Z
!> @param rzmax R at the largest Z
SUBROUTINE cross_section_extrema(r1, z1, sqrts_js, js, noff, ns, nzeta,    &
                                 ntheta2, ntheta3, lasym,                  &
                                 xmin, xmax, zmin, zmax, rzmax)
  REAL(rprec), DIMENSION(:,0:), INTENT(in) :: r1, z1
  REAL(rprec), INTENT(in) :: sqrts_js
  INTEGER, INTENT(in) :: js, noff, ns, nzeta, ntheta2, ntheta3
  LOGICAL, INTENT(in) :: lasym
  REAL(rprec), INTENT(out) :: xmin, xmax, zmin, zmax, rzmax

  INTEGER :: icount, itheta, loff, n1, npasses, ntheta
  REAL(rprec) :: t1, yr1u, yz1u

  zmin =  HUGE(zmin)
  zmax = -HUGE(zmax)
  xmin =  HUGE(xmin)
  xmax = -HUGE(xmax)
  rzmax = 0.0_rprec

  IF (lasym) THEN
     npasses = 1
     ntheta = ntheta3
  ELSE
     npasses = 2
     ntheta = ntheta2
  END IF

  DO icount = 1, npasses
     ! nphi-plane, n1 = noff,...,nzeta
     n1 = noff
     t1 = 1.0_rprec
     IF (icount .eq. 2) THEN
        ! (twopi-v), reflected plane
        n1 = MOD(nzeta-(noff-1),nzeta)+1
        t1 = -1.0_rprec
     END IF
     loff = js + ns*(n1-1)

     DO itheta = 1, ntheta
        yr1u = r1(loff,0) + sqrts_js*r1(loff,1)
        yz1u = t1*(z1(loff,0) + sqrts_js*z1(loff,1))

        IF (yz1u .ge. zmax) THEN
           zmax = yz1u
           rzmax = yr1u
        END IF
        zmin = MIN(zmin, yz1u)
        xmax = MAX(xmax, yr1u)
        xmin = MIN(xmin, yr1u)

        loff = loff + ns*nzeta
     END DO ! itheta
  END DO ! icount

END SUBROUTINE cross_section_extrema

END MODULE cross_section
