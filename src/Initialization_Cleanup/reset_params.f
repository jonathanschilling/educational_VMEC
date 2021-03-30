      SUBROUTINE reset_params
      USE vmec_main, ONLY: iequi, ivac, ftolv, fsqr, fsqz, fsq,
     1                     res0, delt0r, iter1, iter2, ijacob, irst
      USE vmec_input, ONLY: delt
      USE timer_sub, ONLY: timer
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
      res0 = -1       !move to vmec_main, remove from runvmec SAVE
      delt0r = delt   !move delt0 to vmec_main from runvmec SAVE (renamed delt0r)

      timer =  0

      END SUBROUTINE reset_params
