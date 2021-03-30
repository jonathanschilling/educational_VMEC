      SUBROUTINE vsetup ()
      USE vmec_main
      USE vacmod
      USE realspace
      USE mgrid_mod, ONLY: nbcoil_max, nlim_max, nextcur, mgrid_mode
      USE gmres_mod, ONLY: nfcn
      IMPLICIT NONE
C-----------------------------------------------
!
!     Reset default initial values
!
!     m=1 constraint (=t: apply correct, polar constraint; =f, apply approx. constraint)
      lconm1 = .TRUE.

!     2d preconditioner
      nfcn = 0

      z00 = zero
      mgrid_mode = 'S'             !Assume scaled mode; read in from mgrid in free-bdy mode
      nextcur = 0


      END SUBROUTINE vsetup
