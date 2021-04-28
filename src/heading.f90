!> \file
!> \brief Open output files and print banner message at the top.

!> \brief Open output files and print banner message at the top.
!>
!> @param extension input file "extension": part after \c 'input.' .
SUBROUTINE heading(extension)
  USE vmec_main, ONLY: rprec
  USE vparams, ONLY: nthreed
  USE vmec_params, ONLY: version_
  IMPLICIT NONE

  CHARACTER(LEN=50), PARAMETER :: banner     = ' THIS IS VMEC2000, A 3D EQUILIBRIUM CODE, VERSION '
  CHARACTER(LEN=*),  PARAMETER :: VersionID1 = ' Lambda: Full Radial Mesh. L-Force: hybrid full/half.'

  CHARACTER(LEN=*), intent(in) :: extension

  CHARACTER(LEN=50) :: Version
  LOGICAL :: lfirst

  CHARACTER(LEN=100) :: computer, os, os_release

  !     Open output files
  CALL open_output_files (extension, lfirst)

  IF (lfirst) then
     WRITE (*,'(2a)') '  PROCESSING INPUT.', TRIM(extension)

     Version = TRIM(ADJUSTL(version_))

     ! SAL some weird error about file not being ready
     CALL FLUSH(nthreed)

     WRITE(nthreed,'(a,1x,a,/,a,/)') TRIM(banner), TRIM(Version), TRIM(VersionID1)
     WRITE (*,'(1x,a,1x,a,/,1x,a,/)') TRIM(banner), TRIM(Version), TRIM(VersionID1)
  end if

END SUBROUTINE heading
