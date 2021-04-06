!> \file
SUBROUTINE heading(extension, lscreen)
  USE vmec_main, ONLY: rprec
  USE vparams, ONLY: nthreed
  USE vmec_params, ONLY: version_
  USE date_and_computer
  IMPLICIT NONE

  CHARACTER(LEN=*) :: extension
  LOGICAL :: lscreen

  CHARACTER(LEN=50), PARAMETER :: banner    = ' THIS IS VMEC2000, A 3D EQUILIBRIUM CODE, VERSION '
  CHARACTER(LEN=*), PARAMETER :: VersionID1 = ' Lambda: Full Radial Mesh. L-Force: hybrid full/half.'

  INTEGER :: imon, nout
  CHARACTER(LEN=10) :: date0, time0, zone0
  CHARACTER(LEN=50) :: dateloc, Version
  LOGICAL :: lfirst

  !     Open output files
  CALL open_output_files (extension, lscreen, lfirst)

  IF (.not.lfirst) RETURN

  ! FORTRAN-90 ROUTINE
  CALL DATE_AND_TIME(date0,time0,zone0)
  READ(date0(5:6),'(i2)')imon
  WRITE(dateloc,100)months(imon),date0(7:8),date0(1:4), time0(1:2),time0(3:4),time0(5:6)
100  FORMAT('DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)

  CALL GetComputerInfo

  IF (lscreen .and. lfirst) WRITE (*,'(2a)') '  PROCESSING INPUT.', TRIM(extension)

  Version = TRIM(ADJUSTL(version_))

  ! SAL some weird error about file not being ready
  CALL FLUSH(nthreed)

  WRITE(nthreed,'(a,1x,a,/,a,//,3(2a,2x),a)') TRIM(banner), &
       TRIM(Version), TRIM(VersionID1),                     &
       ' COMPUTER: ', TRIM(computer), ' OS: ', TRIM(os),    &
       ' RELEASE: ', TRIM(os_release), TRIM(dateloc)
  IF (lscreen .and. lfirst) then
     WRITE (*,'(1x,a,1x,a,/,1x,a,//,1x,3(2a,2x),a)') TRIM(banner), &
     TRIM(Version), TRIM(VersionID1),                              &
       ' COMPUTER: ', TRIM(computer), ' OS: ', TRIM(os),           &
       ' RELEASE: ', TRIM(os_release), TRIM(dateloc)
  end if

END SUBROUTINE heading
