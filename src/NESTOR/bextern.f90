!> \file
!> Compute the total magnetic field due to external coils and the net toroidal plasma current.

!> Compute the total magnetic field due to external coils and the net toroidal plasma current.
!>
!> @param plascur
!> @param wint
SUBROUTINE bextern(plascur, wint)
  USE vacmod
  USE mgrid_mod, ONLY: bvac

  use dbgout
  use vmec_main, only: num_eqsolve_retries

  IMPLICIT NONE

  REAL(rprec), INTENT(in) :: plascur
  REAL(rprec), DIMENSION(nuv2), INTENT(in) :: wint

  INTEGER :: i
  logical :: dbgout_active

  ! exterior Neumann problem

  IF (.not.ALLOCATED(bvac)) STOP 'BVAC unallocated in bextern'

  ! THIS ROUTINE COMPUTES THE B DOT DS ARISING FROM EXTERNAL COILS AND INTERNAL PLASMA CURRENT
  ! NOTE THAT BEXN = - BEX * DS IS THE EFFECTIVE SOURCE TERM
  !
  ! COMPUTE B FROM COILS ON THE PLASMA BOUNDARY
  ! This sets brad, bphi and bz to the interpolated field from the mgrid.
  CALL becoil(r1b, z1b, bvac(1,1), bvac(1,2), bvac(1,3))

  dbgout_active = open_dbg_context("vac1n_bextern", num_eqsolve_retries)
  if (dbgout_active) then

    ! these are only from the mgrid at this point
    call add_real_2d("mgrid_brad", nv, nu3, brad)
    call add_real_2d("mgrid_bphi", nv, nu3, bphi)
    call add_real_2d("mgrid_bz",   nv, nu3, bz)

    ! This should be in Amperes.
    call add_real("axis_current", plascur/mu0)

  end if ! dbgout_active

  ! COMPUTE B (ON PLASMA BOUNDARY) FROM NET TOROIDAL PLASMA CURRENT
  ! THE NET CURRENT IS MODELLED AS A WIRE AT THE MAGNETIC AXIS, AND THE
  ! BIOT-SAVART LAW IS USED TO COMPUTE THE FIELD AT THE PLASMA SURFACE
  !
  ! USE BEXU, BEXV, BEXN AS TEMPORARY STORAGE FOR BX, BY, BZ
  ! --> add to interpolated field from mgrid
  CALL belicu (plascur, bexu, bexv, bexn, cosuv, sinuv, r1b, z1b)
  DO i = 1, nuv2
     brad(i) = brad(i) + bexu(i)*cosuv(i) + bexv(i)*sinuv(i)
     bphi(i) = bphi(i) - bexu(i)*sinuv(i) + bexv(i)*cosuv(i)
     bz(i)   = bz(i) + bexn(i)
  END DO

  ! COMPUTE COVARIANT COMPONENTS OF EXTERNAL FIELD: BEXU = B0 dot dx/du,
  ! BEXV = B0 dot dx/dv. HERE, BEXN = -B0*SURF_NORM CORRESPONDS TO THE
  ! "exterior Neumann problem" convention of PKM (sign flipped as noted in PKM)
  ! THUS, THE UNIT NORMAL SHOULD POINT INTO THE PLASMA (OUTWARD FROM VACUUM),
  ! WHICH IT DOES FOR A NEGATIVE JACOBIAN (SIGNGS) SYSTEM
  DO i = 1, nuv2
    bexu(i) = rub(i)*brad(i) + zub(i)*bz(i)
    bexv(i) = rvb(i)*brad(i) + zvb(i)*bz(i) + r1b(i)*bphi(i)
    bexn(i) =-(brad(i)*snr(i) + bphi(i)*snv(i) + bz(i)*snz(i))
  END DO

  ! COMPUTE NORMALIZED [(2*pi)**2], READY-TO-INTEGRATE (WINT FACTOR) SOURCE TERM
  ! NOTE: BEXN == NP*F = -B0 dot [Xu cross Xv] NP        (see PKM, Eq. 2.13)
  bexni(:nuv2) = wint(:nuv2)*bexn(:nuv2)*pi2*pi2

  if (dbgout_active) then

    ! axis geometry used in belicu
    call add_real_2d("xpts_axis", 3, nvper * nv + 1, xpts)

    ! these are now mgrid + axis-current
    call add_real_2d("brad", nv, nu3, brad)
    call add_real_2d("bphi", nv, nu3, bphi)
    call add_real_2d("bz",   nv, nu3, bz)

    call add_real_2d("bexu", nv, nu3, bexu)
    call add_real_2d("bexv", nv, nu3, bexv)
    call add_real_2d("bexn", nv, nu3, bexn)

    call add_real_2d("bexni", nv, nu3, bexni(:nuv2))

    call close_dbg_out()
  end if

END SUBROUTINE bextern
