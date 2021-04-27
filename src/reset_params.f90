!> \file
SUBROUTINE reset_params
  USE vmec_main, ONLY: ivac, ftolv, fsqr, fsqz, fsq,      &
                       res0, delt0r, iter1, iter2, ijacob, irst, &
                       lconm1, z00, vacuum_calls
  USE vmec_input, ONLY: delt
  USE stel_constants, only: zero
  USE mgrid_mod, ONLY: nextcur, mgrid_mode
  IMPLICIT NONE

  ivac  = -1

  fsqr = 1
  fsqz = 1
  ftolv = fsqr

  fsq   = 1

  iter2 = 1
  iter1 = iter2

  ijacob = 0

  irst = 1

  res0 = -1

  delt0r = delt

! #ifndef _HBANGLE
  !> m=1 constraint (=t: apply correct, polar constraint; =f, apply approx. constraint)
  lconm1 = .TRUE.
! #end /* ndef _HBANGLE */

  z00 = zero
  ! r00 gets assiged in profil1d and in funct3d

  !> Assume scaled mode; read in from mgrid in free-bdy mode
  mgrid_mode = 'S'
  nextcur = 0

  vacuum_calls = 0

END SUBROUTINE reset_params
