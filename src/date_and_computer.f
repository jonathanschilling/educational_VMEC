!> \file
      MODULE date_and_computer
      USE safe_open_mod
      IMPLICIT NONE

      CHARACTER(LEN=100) :: computer, os, os_release
      CONTAINS

      SUBROUTINE GetComputerInfo
      INTEGER :: ierror, ipid, iunit=10
      CHARACTER(LEN=200) :: fileId

      CALL GETENV('HOST',computer)
      CALL GETENV('OSTYPE',os)
      CALL GETENV('HOSTTYPE',os_release)

      END SUBROUTINE GetComputerInfo

      END MODULE date_and_computer
