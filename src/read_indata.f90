!> \file
SUBROUTINE read_indata(in_file, iunit, ier_flag)
  USE vmec_main
  USE vmec_input, ONLY: bloat, ncurr
  USE vmec_params
  USE vacmod
  USE safe_open_mod
  IMPLICIT NONE

  INTEGER ier_flag, iunit
  CHARACTER(LEN=*) :: in_file

  INTEGER :: ireadseq, iosnml = 0

  iunit = indata0
  CALL safe_open (iunit, ireadseq, in_file, 'old', 'formatted')
  IF (ireadseq .ne. 0) THEN
     WRITE (6, '(3a,i4)') ' In VMEC, error opening input file: ',TRIM(in_file),'. Iostat = ', ireadseq
     ier_flag = input_error_flag
     RETURN
  ENDIF

  CALL read_namelist (iunit, iosnml, 'indata')
  IF (iosnml .ne. 0) THEN
     WRITE (6, '(a,i4)') ' In VMEC, indata NAMELIST error: iostat = ', iosnml
     ier_flag = input_error_flag
     RETURN
  ENDIF

  IF (lfreeb .and. mgrid_file.eq.'NONE') then
     ! disable free-boundary mode if mgrid file is not specified
     lfreeb = .false.
  end if

  IF (bloat .eq. zero) bloat = one
  IF ((bloat.ne.one) .and. (ncurr.ne.1)) THEN
     ier_flag = 3 ! 'VMEC INDATA ERROR: NCURR.ne.1 but BLOAT.ne.1.'
     RETURN
  ENDIF

  ! fixup current profile
  IF (ncurr.eq.1 .and. ALL(ac.eq.cbig)) then
     ! previous version input of current profile: via ai (iota profile coeffs)
     ! Old FORMAT: may not be reading in ac
     ac = ai
  end if

  WHERE (ac .eq. cbig) ac = zero

  ! COMPUTE NTHETA, NZETA VALUES
  mpol = ABS(mpol)
  ntor = ABS(ntor)
  IF (mpol .gt. mpold) STOP 'mpol>mpold: lower mpol'
  IF (ntor .gt. ntord) STOP 'ntor>ntord: lower ntor'
  mpol1 = mpol - 1
  ntor1 = ntor + 1

  IF (ntheta .lt. 2*mpol+6 ) THEN
    ! number of theta grid points (>=2*mpol+6)
    ntheta = 2*mpol+6
  ENDIF

  ntheta1 = 2*(ntheta/2)
  ! u = pi
  ntheta2 = 1 + ntheta1/2
  IF (ntor .eq. 0) then
     lthreed = .false.
  else
     lthreed = .true.
  end if

  IF (ntor.eq.0 .and. nzeta.eq.0) then
     ! Tokamak (ntor=0) needs nzeta=1
     nzeta = 1
  end if

  IF (ntor.gt.0)then
    ! Stellarator case needs Nyquist criterion fulfilled for nzeta wrt. ntor
    IF (nzeta .lt. 2*ntor+4) THEN
      nzeta = 2*ntor+4      !number of zeta grid points (=1 IF ntor=0)
    ENDIF
  ENDIF

  ! SIZE of rmnc,  rmns,  ...
  mnmax = ntor1 + mpol1*(1 + 2*ntor)

  ! SIZE of rmncc, rmnss, ...
  mnsize = mpol*ntor1

  ! INDEXING FOR PACKED-ARRAY STRUCTURE OF XC, GC
  rcc = 1;  zsc = 1
  rss = 0;  rsc = 0;  rcs = 0
  zcc = 0;  zss = 0;  zcs = 0
  IF (.NOT.lasym) THEN
     ! can make use of Stellarator symmetry
     ntheta3 = ntheta2
     IF (lthreed) THEN
        ntmax = 2
        rss = 2;  zcs = 2
     ELSE
        ntmax = 1
     END IF
  ELSE
     ntheta3 = ntheta1
     IF (lthreed) THEN
         ntmax = 4
         rss = 2;  rsc = 3;  rcs = 4
         zcs = 2;  zcc = 3;  zss = 4
     ELSE
         ntmax = 2
         rsc = 2;  zcc = 2
     END IF
  END IF

  nznt = nzeta*ntheta3



  ! below stuff is required for NESTOR
  mf = mpol+1
  nf = ntor
  nu = ntheta1
  nv = nzeta
  mf1 = 1+mf
  nf1 = 2*nf+1
  mnpd = mf1*nf1

  IF (.NOT.lasym) THEN
    mnpd2 = mnpd
  ELSE
    mnpd2 = 2*mnpd
  END IF

  nuv = nu*nv
  nfper = nfp
  nu2 = nu/2 + 1
  nu3 = ntheta3
  nuv2 = nznt


END SUBROUTINE read_indata
