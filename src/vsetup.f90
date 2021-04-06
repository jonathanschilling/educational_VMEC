!> \file
SUBROUTINE vsetup
  USE vmec_main
  USE vacmod
  USE realspace
  USE mgrid_mod, ONLY: nbcoil_max, nlim_max, nextcur, mgrid_mode
  IMPLICIT NONE

  ! Reset default initial values

  !> m=1 constraint (=t: apply correct, polar constraint; =f, apply approx. constraint)
  lconm1 = .TRUE.

  z00 = zero

  !> ! Assume scaled mode; read in from mgrid in free-bdy mode
  mgrid_mode = 'S'

  nextcur = 0

END SUBROUTINE vsetup
