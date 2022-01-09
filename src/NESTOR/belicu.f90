!> \file
!> \brief Magnetic field due to net toroidal current modeled by a filament along the magnetic axis.

!> \brief Magnetic field due to net toroidal current modeled by a filament along the magnetic axis.
!>
!> @param torcur
!> @param bx
!> @param by
!> @param bz
!> @param cos1
!> @param sin1
!> @param rp
!> @param zp
SUBROUTINE belicu(torcur, bx, by, bz, cos1, sin1, rp, zp)

  USE vacmod, vm_bz => bz

  IMPLICIT NONE

  REAL(rprec), INTENT(in) :: torcur
  REAL(rprec), DIMENSION(nuv2), INTENT(in)  :: cos1, sin1, rp, zp
  REAL(rprec), DIMENSION(nuv2), INTENT(out) :: bx, by, bz

  REAL(rprec) :: current
  INTEGER :: i, j, kper, kv
  REAL(rprec), DIMENSION(3) :: xpt, bvec, dvec, Ri_vec

  ! quantities from Hanson & Hirshman, "Compact expressions for the Biot-Savart fields of a filamentary segment" (2002)
  REAL(rprec) :: L, Ri, Rf, Ri_p_Rf, Bmag

  ! B_External due to LIne CUrrent

  ! net toroidal plasma current in A
  current = torcur/mu0
  
  ! initialize target array
  bx = 0
  by = 0
  bz = 0

  ! first point (at index 0) is equal to last point --> closed curve
  xpts(1, 0) = raxis_nestor(nv)*(cosper(nvper)*cosuv(nv) - sinper(nvper)*sinuv(nv))
  xpts(2, 0) = raxis_nestor(nv)*(sinper(nvper)*cosuv(nv) + cosper(nvper)*sinuv(nv))
  xpts(3, 0) = zaxis_nestor(nv)

  ! loops over source geometry
  i = 1
  DO kper = 1, nvper
     DO kv = 1, nv

        ! xpts == xpt of _s_ource (current filament)
        xpts(1, i) = raxis_nestor(kv)*(cosper(kper)*cosuv(kv) - sinper(kper)*sinuv(kv))
        xpts(2, i) = raxis_nestor(kv)*(sinper(kper)*cosuv(kv) + cosper(kper)*sinuv(kv))
        xpts(3, i) = zaxis_nestor(kv)

        ! filament geometry: from current point (R_i == xpts(:,i)) to previous point (R_f == xpts(:,i-1))
        dvec = xpts(:,i)-xpts(:,i-1)
        L = norm2(dvec)

        ! loop over evaluation points
        DO j = 1,nuv2
           ! xpt == evaluation position
           xpt(1) = rp(j) * cos1(j)
           xpt(2) = rp(j) * sin1(j)
           xpt(3) = zp(j)

           Ri_vec = xpt - xpts(:,i-1)
           Ri = norm2(Ri_vec)
           Rf = norm2(xpt - xpts(:,i))
           Ri_p_Rf = Ri + Rf

           ! 1.0e-7 == mu0/4 pi
           Bmag = 1.0E-7_dp * current * 2.0_dp * Ri_p_Rf / ( Ri * Rf * (Ri_p_Rf*Ri_p_Rf - L*L) )

           ! cross product of L*hat(eps)==dvec with Ri_vec, scaled by Bmag
           bx(j) = bx(j) + Bmag * (dvec(2)*Ri_vec(3) - dvec(3)*Ri_vec(2))
           by(j) = by(j) + Bmag * (dvec(3)*Ri_vec(1) - dvec(1)*Ri_vec(3))
           bz(j) = bz(j) + Bmag * (dvec(1)*Ri_vec(2) - dvec(2)*Ri_vec(1))
        end do

        i = i + 1
     END DO
  END DO

END SUBROUTINE belicu
