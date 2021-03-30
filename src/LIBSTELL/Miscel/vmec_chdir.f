!> \file
      INTEGER FUNCTION vmec_chdir(new_path)
      IMPLICIT NONE
      CHARACTER*(*), INTENT(in) :: new_path
      INTEGER, EXTERNAL :: chdir

!      vmec_chdir = chdir(TRIM(new_path) // char(0))
      END FUNCTION vmec_chdir
