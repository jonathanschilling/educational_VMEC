!> \file
!> \brief Reset some flow-control parameters to their default values.

!> \brief Reset some flow-control parameters to their default values.
!>
SUBROUTINE reset_params

  USE vmec_main, ONLY: ivac, ftolv, fsqr, fsqz, fsq,      &
                       res0, delt0r, iter1, iter2, ijacob, first, &
                       lconm1, z00, dp
  USE vmec_input, ONLY: delt
  USE stel_constants, only: zero
  USE mgrid_mod, ONLY: nextcur, mgrid_mode

  IMPLICIT NONE

  ivac  = -1

  fsqr = 1.0_dp
  fsqz = 1.0_dp
  ftolv = fsqr

  fsq   = 1.0_dp

  iter2 = 1
  iter1 = iter2

  ijacob = 0

  first = 1

  res0 = -1.0_dp

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

END SUBROUTINE reset_params
