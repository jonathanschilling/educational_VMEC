!> \file
      SUBROUTINE read_indata(in_file, iunit, ier_flag)
      USE vmec_main
      USE vmec_input, ONLY: bloat, ncurr
      USE vmec_params
      USE vacmod
      USE safe_open_mod
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER ier_flag, iunit
      CHARACTER(LEN=*) :: in_file
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ireadseq, iosnml = 0
!-----------------------------------------------
      iunit = indata0
      CALL safe_open (iunit, ireadseq, in_file, 'old', 'formatted')
      IF (ireadseq .ne. 0) THEN
         WRITE (6, '(3a,i4)') ' In VMEC, error opening input file: ',
     1   TRIM(in_file), '. Iostat = ', ireadseq
         ier_flag = input_error_flag
         RETURN
      ENDIF

      CALL read_namelist (iunit, iosnml, 'indata')
      IF (iosnml .ne. 0) THEN
         WRITE (6, '(a,i4)')
     1   ' In VMEC, indata NAMELIST error: iostat = ', iosnml
         ier_flag = input_error_flag
         RETURN
      ENDIF

      IF (lfreeb .and. mgrid_file.eq.'NONE') lfreeb = .false.

      IF (bloat .eq. zero) bloat = one
      IF ((bloat.ne.one) .and. (ncurr.ne.1)) THEN
         ier_flag = 3
         RETURN
      ENDIF
!
!     COMPUTE NTHETA, NZETA VALUES
!
      mpol = ABS(mpol)
      ntor = ABS(ntor)
      IF (mpol .gt. mpold) STOP 'mpol>mpold: lower mpol'
      IF (ntor .gt. ntord) STOP 'ntor>ntord: lower ntor'
      mpol1 = mpol - 1
      ntor1 = ntor + 1
! 20130924 J.Geiger, assure minimum ntheta-value
!                    and add a notification
      IF (ntheta .lt. 2*mpol+6 ) THEN
!        WRITE(6,*)"Adjust NTHETA from ",ntheta,
!     1            " to new value: ",2*mpol+6
        ntheta = 2*mpol+6    !number of theta grid points (>=2*mpol+6)
      ENDIF
      ntheta1 = 2*(ntheta/2)
      ntheta2 = 1 + ntheta1/2                   !u = pi
      IF (ntor .eq. 0) lthreed = .false.
      IF (ntor .gt. 0) lthreed = .true.

      IF (ntor.eq.0 .and. nzeta.eq.0) nzeta = 1
! 20130924 J.Geiger, assure minimum nzeta-value for ntor>0
!                    and add a notification
      IF (ntor.gt.0)then
        IF (nzeta .lt. 2*ntor+4) THEN
!          WRITE(6,*)"Adjust NZETA from ",nzeta,
!     1              " to new value: ",2*ntor+4
          nzeta = 2*ntor+4      !number of zeta grid points (=1 IF ntor=0)
        ENDIF
      ENDIF
      mnmax = ntor1 + mpol1*(1 + 2*ntor)        !SIZE of rmnc,  rmns,  ...
      mnsize = mpol*ntor1                       !SIZE of rmncc, rmnss, ...

      mf = mpol+1
      nf = ntor
      nu = ntheta1
      nv = nzeta
      mf1 = 1+mf
      nf1 = 2*nf+1
      mnpd = mf1*nf1
!
!     INDEXING FOR PACKED-ARRAY STRUCTURE OF XC, GC
!
      rcc = 1;  zsc = 1
      rss = 0;  rsc = 0;  rcs = 0
      zcc = 0;  zss = 0;  zcs = 0
      IF (.NOT.lasym) THEN
         ntheta3 = ntheta2
         mnpd2 = mnpd
         IF (lthreed) THEN
            ntmax = 2
            rss = 2;  zcs = 2
         ELSE
            ntmax = 1
         END IF
      ELSE
         ntheta3 = ntheta1
         mnpd2 = 2*mnpd
         IF (lthreed) THEN
             ntmax = 4
             rss = 2;  rsc = 3;  rcs = 4
             zcs = 2;  zcc = 3;  zss = 4
         ELSE
             ntmax = 2
             rsc = 2;  zcc = 2
         END IF
      END IF

      nuv = nu*nv
      nznt = nzeta*ntheta3
      nfper = nfp
      nu2 = nu/2 + 1
      nu3 = ntheta3
      nuv2 = nznt

      IF (ncurr.eq.1 .and. ALL(ac.eq.cbig)) ac = ai            !!Old FORMAT: may not be reading in ac
      WHERE (ac .eq. cbig) ac = zero

      END SUBROUTINE read_indata