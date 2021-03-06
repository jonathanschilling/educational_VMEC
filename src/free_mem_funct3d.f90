!> \file
!> \brief Free memory required by funct3d()

!> \brief Free memory required by funct3d()
!>
SUBROUTINE free_mem_funct3d
  USE vmec_main
  USE realspace
  USE vforces
  IMPLICIT NONE

  INTEGER :: istat1 = 0

  IF (ALLOCATED(armn)) then
     DEALLOCATE (armn, azmn, brmn, bzmn, crmn, czmn, blmn, clmn,    &
        r1, ru, rv, z1, zu, zv, gcon, rcon, zcon, ru0, zu0,         &
        brmn_con, bzmn_con, &
        rcon0, zcon0, guu, guv, gvv, stat=istat1)
  IF (istat1 .ne. 0) STOP 'deallocation error#1 in funct3d'
  end if

  IF (ALLOCATED(extra1)) then
     DEALLOCATE (extra1, extra2, extra3, extra4, stat=istat1)
     IF (istat1 .ne. 0) STOP 'deallocation error#3 in funct3d'
  end if

END SUBROUTINE free_mem_funct3d
