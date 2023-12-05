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

  real(rprec), dimension(3, nuv2) :: eval_pos, magnetic_field

  ! If .true., use ABSCAB for computing the magnetic field contribution
  ! due to the net toridal current modeled as a filament along the magnetic axis.
  logical, parameter :: use_abscab_for_axis_current = .false.

  ! B_External due to LIne CUrrent

  ! net toroidal plasma current in A
  current = torcur/mu0

  ! loops over source geometry
  i = 0
  DO kper = 1, nvper
      DO kv = 1, nv
          i = i + 1

          ! xpts == xpt of _s_ource (current filament)
          xpts(1, i) = raxis_nestor(kv)*(cosper(kper)*cosuv(kv) - sinper(kper)*sinuv(kv))
          xpts(2, i) = raxis_nestor(kv)*(sinper(kper)*cosuv(kv) + cosper(kper)*sinuv(kv))
          xpts(3, i) = zaxis_nestor(kv)
      END DO
  END DO

  ! last point is equal to first point --> closed curve
  xpts(1, nvp + 1) = xpts(1, 1)
  xpts(2, nvp + 1) = xpts(2, 1)
  xpts(3, nvp + 1) = xpts(3, 1)

  if (use_abscab_for_axis_current) then

    DO j = 1, nuv2
        ! evaluation positions
        eval_pos(1, j) = rp(j) * cos1(j)
        eval_pos(2, j) = rp(j) * sin1(j)
        eval_pos(3, j) = zp(j)
    end do

    ! initialize target array
    magnetic_field = 0.0_dp

    ! use ABSCAB to compute the line-current-along-axis magnetic field contribution
    call magneticFieldPolygonFilament(nvper * nv + 1, xpts, current,  &
                                      nuv2, eval_pos, magnetic_field)

    bx(:) = magnetic_field(1,:)
    by(:) = magnetic_field(2,:)
    bz(:) = magnetic_field(3,:)

  else ! use_abscab_for_axis_current

    ! initialize target array
    bx = 0
    by = 0
    bz = 0

    ! iterate over all wire segments that make up the axis;
    ! the number of wire segments is one less than number of points of the closed loop
    DO i = 1, nvper * nv

        ! filament geometry: from current point (R_i == xpts(:,i)) to previous point (R_f == xpts(:,i-1))
        dvec = xpts(:,i+1)-xpts(:,i)
        L = norm2(dvec)

        ! loop over evaluation points
        DO j = 1,nuv2
            ! xpt == evaluation position
            xpt(1) = rp(j) * cos1(j)
            xpt(2) = rp(j) * sin1(j)
            xpt(3) = zp(j)

            Ri_vec = xpt - xpts(:,i)
            Ri = norm2(Ri_vec)
            Rf = norm2(xpt - xpts(:,i + 1))
            Ri_p_Rf = Ri + Rf

            ! 1.0e-7 == mu0/4 pi
            Bmag = 1.0E-7_dp * current * 2.0_dp * Ri_p_Rf / ( Ri * Rf * (Ri_p_Rf*Ri_p_Rf - L*L) )

            ! cross product of L*hat(eps)==dvec with Ri_vec, scaled by Bmag
            bx(j) = bx(j) + Bmag * (dvec(2)*Ri_vec(3) - dvec(3)*Ri_vec(2))
            by(j) = by(j) + Bmag * (dvec(3)*Ri_vec(1) - dvec(1)*Ri_vec(3))
            bz(j) = bz(j) + Bmag * (dvec(1)*Ri_vec(2) - dvec(2)*Ri_vec(1))
        end do

    END DO

  end if ! use_abscab_for_axis_current

END SUBROUTINE belicu
