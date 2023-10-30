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

  use abscab, only: magneticFieldPolygonFilament

  IMPLICIT NONE

  REAL(rprec), INTENT(in) :: torcur
  REAL(rprec), DIMENSION(nuv2), INTENT(in)  :: cos1, sin1, rp, zp
  REAL(rprec), DIMENSION(nuv2), INTENT(out) :: bx, by, bz

  REAL(rprec) :: current
  INTEGER :: i, j, kper, kv

  real(rprec), dimension(3, nuv2) :: eval_pos, magnetic_field

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
     end do
  end do

  ! last point is equal to first point --> closed curve
  xpts(1, nvp+1) = xpts(1, 1)
  xpts(2, nvp+1) = xpts(2, 1)
  xpts(3, nvp+1) = xpts(3, 1)

  DO j = 1, nuv2
    ! evaluation positions
    eval_pos(1, j) = rp(j) * cos1(j)
    eval_pos(2, j) = rp(j) * sin1(j)
    eval_pos(3, j) = zp(j)
  end do

  ! initialize target array
  magnetic_field = 0.0_dp

  ! use ABSCAB to compute the line-current-along-axis magnetic field contribution
  call magneticFieldPolygonFilament(nvper * nv + 1, xpts, current, &
                                    nuv2, eval_pos, magnetic_field)

  bx(:) = magnetic_field(1,:)
  by(:) = magnetic_field(2,:)
  bz(:) = magnetic_field(3,:)

END SUBROUTINE belicu
