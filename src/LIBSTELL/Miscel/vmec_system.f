      SUBROUTINE vmec_system(cmd, ierror)
      INTEGER, OPTIONAL :: ierror
      INTEGER :: ireturn
      CHARACTER(LEN=*), INTENT(in) :: cmd

      INTEGER :: system
      ireturn = system(TRIM(cmd))

      IF (PRESENT(ierror)) ierror = ireturn

      END SUBROUTINE vmec_system
