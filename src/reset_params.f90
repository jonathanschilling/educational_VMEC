!> \file
SUBROUTINE reset_params
  USE vmec_main, ONLY: iequi, ivac, ftolv, fsqr, fsqz, fsq,     &
                       res0, delt0r, iter1, iter2, ijacob, irst, &
                       lconm1, z00
  USE vmec_input, ONLY: delt
  USE timer_sub, ONLY: timer
  USE vacmod
  USE realspace
  USE mgrid_mod, ONLY: nbcoil_max, nlim_max, nextcur, mgrid_mode
  IMPLICIT NONE

  iequi = 0
  ivac  = -1

  fsqr = 1
  fsqz = 1
  ftolv = fsqr

  fsq   = 1
  iter2 = 1
  iter1 = iter2
  ijacob = 0
  irst = 1
  res0 = -1       ! move       to vmec_main, remove from runvmec SAVE
  delt0r = delt   ! move delt0 to vmec_main, remove from runvmec SAVE (renamed delt0r)

  timer =  0

  ! below stuff is from vsetup

  !> m=1 constraint (=t: apply correct, polar constraint; =f, apply approx. constraint)
  lconm1 = .TRUE.

  z00 = zero

  !> ! Assume scaled mode; read in from mgrid in free-bdy mode
  mgrid_mode = 'S'

  nextcur = 0

END SUBROUTINE reset_params