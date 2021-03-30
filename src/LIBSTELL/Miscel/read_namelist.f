!> \file
      SUBROUTINE read_namelist(iunit, io_stat, lc_name)
      USE vmec_input, ONLY: read_indata_namelist
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: iunit, io_stat
      CHARACTER(LEN=*) :: lc_name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ifind
      CHARACTER(LEN=1), PARAMETER :: lead = '&'
      CHARACTER(LEN=LEN_TRIM(lc_name)+1) :: namelc
!-----------------------------------------------

      io_stat = 0
      REWIND (iunit)
      namelc = lead // TRIM(ADJUSTL(lc_name))

      ifind = MIN(len_trim(namelc), 132)
      IF (namelc(1:ifind) .eq. '&indata') THEN
         CALL read_indata_namelist (iunit, io_stat)
      END IF

      END SUBROUTINE read_namelist
