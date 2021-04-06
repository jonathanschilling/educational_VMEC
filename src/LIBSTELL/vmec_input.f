!> \file
      MODULE vmec_input
      USE vparams, ONLY: rprec, dp, mpol1d, ntord, ndatafmax
      USE vsvd0
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!   For variable descriptions, see VMEC "readin.f" routine
!-----------------------------------------------
      INTEGER, PARAMETER :: mpol_default = 6
      INTEGER, PARAMETER :: ntor_default = 0
      INTEGER, PARAMETER :: ns_default   = 31

      INTEGER :: nfp, ncurr, nsin, niter, nstep, nvacskip, mpol, ntor,
     1           ntheta, nzeta, mfilter_fbdy, nfilter_fbdy
      INTEGER, DIMENSION(100) :: ns_array, niter_array
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) ::
     1   rbs, zbc, rbc, zbs
      REAL(rprec) :: curtor, delt, ftol, tcon0,
     1   gamma, bloat, pres_scale,
     4   prec2d_threshold
      REAL(rprec) :: spres_ped !< value of s beyond which pressure profile is flat (pedestal)
      REAL(rprec) :: phiedge   !< value of real toroidal flux at plasma edge (s=1)
      REAL(rprec), DIMENSION(0:20) :: am !< array of coefficients in phi-series for mass (NWT/m**2)
      REAL(rprec), DIMENSION(0:20) :: ai !< array of coefficients in phi-series for iota (ncurr=0)
      REAL(rprec), DIMENSION(0:20) :: ac !< array of coefficients in phi-series for the quantity d(Icurv)/ds = toroidal
                                         !< current density * Vprime, so Icurv(s) = Itor(s) (used for ncurr=1)
      REAL(rprec), DIMENSION(1:20) :: aphi
      CHARACTER(len=20) :: pcurr_type  !  len=12 -> len=20 J Hanson 2010-03-16
      CHARACTER(len=20) :: piota_type
      CHARACTER(len=20) :: pmass_type
      REAL(rprec), DIMENSION(ndatafmax) :: am_aux_s, am_aux_f,                 &
     &   ai_aux_s, ai_aux_f, ac_aux_s, ac_aux_f

      REAL(rprec), DIMENSION(0:ntord) :: raxis_cc, raxis_cs,
     1                                   zaxis_cc, zaxis_cs
      REAL(rprec), DIMENSION(100) :: ftol_array
      REAL(rprec), DIMENSION(nigroup), TARGET :: extcur ! V3FIT needs a pointer to this.
      LOGICAL :: lfreeb, lasym
      LOGICAL :: lbsubs                   ! J Hanson See jxbforce coding

      CHARACTER(len=200) :: mgrid_file
      CHARACTER(len=10)  :: precon_type
      CHARACTER(len=120) :: arg1
      CHARACTER(len=100) :: input_extension

      NAMELIST /indata/ mgrid_file, nfp, ncurr, nsin,
     1   niter, nstep, nvacskip, delt, ftol, gamma, am, ai, ac, aphi,
     1   pcurr_type, pmass_type, piota_type, bloat,
     1   am_aux_s, am_aux_f, ai_aux_s, ai_aux_f, ac_aux_s, ac_aux_f,  ! J Hanson 2010-03-16
     2   rbc, zbs, rbs, zbc, spres_ped, pres_scale, raxis_cc, zaxis_cs,
     3   raxis_cs, zaxis_cc, mpol, ntor, ntheta, nzeta, mfilter_fbdy,
     3   nfilter_fbdy, niter_array,
     4   ns_array, ftol_array, tcon0, precon_type, prec2d_threshold,
     4   curtor, extcur,
     5   phiedge,
     A   lfreeb,                                                           ! S Lazerson 2010
     D   lasym,
     E   lbsubs                                                            ! 2014-01-12 See jxbforce

      CONTAINS

      SUBROUTINE read_indata_namelist (iunit, istat)
      INTEGER :: iunit, istat

!
!     INITIALIZATIONS
!
      gamma = 0
      spres_ped = 1
      mpol = mpol_default
      ntor = ntor_default
      ntheta = 0;  nzeta = 0
      ns_array = 0;  ns_array(1) = ns_default
      niter_array = -1;
      bloat = 1
      rbc = 0;  rbs = 0; zbs = 0; zbc = 0
      nfp = 1
      ncurr = 0
      nsin = ns_default
      niter = 100
      nstep = 10
      nvacskip = 1
      delt = 1
      ftol = 1.E-10_dp
      ftol_array = 0;  ftol_array(1) = ftol
      am = 0; ai = 0; ac = 0; aphi = 0; aphi(1) = 1
      pres_scale = 1
      raxis_cc = 0; zaxis_cs = 0; raxis_cs = 0; zaxis_cc = 0;
      mfilter_fbdy = -1; nfilter_fbdy = -1
      tcon0 = 1
      precon_type = 'NONE'; prec2d_threshold = 1.E-30_dp
      curtor = 0;
      extcur = 0;  phiedge = 1;
      mgrid_file = 'NONE'
      lfreeb = .true.
      lasym = .false.
      lbsubs = .false.         ! J Hanson. See jxbforce coding

      pcurr_type = 'power_series'
      piota_type = 'power_series'
      pmass_type = 'power_series'

      am_aux_s(:) = -1
      ac_aux_s(:) = -1
      ai_aux_s(:) = -1

      READ (iunit, nml=indata, iostat=istat)

      IF (ALL(niter_array == -1)) niter_array = niter

      END SUBROUTINE read_indata_namelist

      SUBROUTINE write_indata_namelist (iunit, istat)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: iunit
      INTEGER, INTENT(inout) :: istat
      INTEGER :: iftol,i,n,m
      INTEGER, DIMENSION(1) :: ins
      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outint1 = "(2X,A,1X,'=',1X,I1.1)"
      CHARACTER(LEN=*), PARAMETER :: outint2 = "(2X,A,1X,'=',1X,I2.2)"
      CHARACTER(LEN=*), PARAMETER :: outint3 = "(2X,A,1X,'=',1X,I3.3)"
      CHARACTER(LEN=*), PARAMETER :: outint4 = "(2X,A,1X,'=',1X,I4.4)"
      CHARACTER(LEN=*), PARAMETER :: outint5 = "(2X,A,1X,'=',1X,I5.5)"
      CHARACTER(LEN=*), PARAMETER :: outint6 = "(2X,A,1X,'=',1X,I6.6)"
      CHARACTER(LEN=*), PARAMETER :: outflt="(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outexp="(2X,A,1X,'=',1X,ES22.12E3)"
      IF (istat < 0) RETURN
      WRITE(iunit,'(A)') '!----- Runtime Parameters -----'
      WRITE(iunit,'(A)') '&INDATA'
      WRITE(iunit,outflt) 'DELT',delt
      WRITE(iunit,outint) 'NITER',niter
      WRITE(iunit,outint) 'NSTEP',nstep
      WRITE(iunit,outflt) 'TCON0',tcon0
      ins = MAXLOC(ns_array)
      WRITE(iunit,'(a,(1p,4i14))')  '  NS_ARRAY =    ',
     1     (ns_array(i), i=1,ins(1))
      iftol = 1
      DO WHILE(ftol_array(iftol).ne.0 .and. iftol.lt.100)
         iftol = iftol + 1
      END DO
      WRITE(iunit,'(a,(1p,4e14.6))')'  FTOL_ARRAY =  ',
     1     (ftol_array(i), i=1,iftol - 1)
      ins = MINLOC(niter_array)
      IF (ins(1) > 1)
     1 WRITE(iunit,'(a,(1p,4i14))') '  NITER_ARRAY = ',
     2      (niter_array(i), i=1,ins(1)-1)
      WRITE (iunit,'(2x,3a)') "PRECON_TYPE = '", TRIM(precon_type),"'"
      WRITE (iunit,'(2x,a,1p,e14.6)') "PREC2D_THRESHOLD = ",
     1                                prec2d_threshold
      WRITE(iunit,'(A)') '!----- Grid Parameters -----'
      WRITE(iunit,outboo) 'LASYM',lasym
      WRITE(iunit,outint4) 'NFP',nfp
      WRITE(iunit,outint4) 'MPOL',mpol
      WRITE(iunit,outint4) 'NTOR',ntor
      WRITE(iunit,outflt) 'PHIEDGE',phiedge
      WRITE(iunit,'(A)') '!----- Free Boundary Parameters -----'
      WRITE(iunit,outboo) 'LFREEB',lfreeb
      IF (lfreeb) THEN
         WRITE (iunit, '(2x,3a)') "MGRID_FILE = '",TRIM(mgrid_file),"'"
         WRITE(iunit,outint4) 'NZETA',nzeta
         DO n=1,SIZE(extcur)
            IF (extcur(n) == 0) CYCLE
            WRITE(iunit,'(2X,A,I3.3,A,ES22.12E3)')
     1      'EXTCUR(',n,') = ',extcur(n)
         END DO
         WRITE(iunit,outint4) 'NVACSKIP',nvacskip
      END IF
      WRITE(iunit,'(A)') '!----- Pressure Parameters -----'
      WRITE(iunit,outflt) 'GAMMA',gamma
      WRITE(iunit,outflt) 'BLOAT',bloat
      WRITE(iunit,outflt) 'SPRES_PED',spres_ped
      WRITE(iunit,outflt) 'PRES_SCALE',pres_scale
      WRITE(iunit,'(2x,3a)') "PMASS_TYPE = '",TRIM(pmass_type),"'"
      WRITE(iunit,'(a,(1p,4e22.14))')'  AM = ', (am(i-1), i=1,SIZE(am))
      i = minloc(am_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AM_AUX_S = ',
     1         (am_aux_s(n), n=1,i)
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AM_AUX_F = ',
     1         (am_aux_f(n), n=1,i)
      END IF

      WRITE(iunit,'(A)') '!----- Current/Iota Parameters -----'
      WRITE(iunit,outexp) 'CURTOR',curtor
      WRITE(iunit,outint) 'NCURR',ncurr
      WRITE (iunit, '(2x,3a)') "PIOTA_TYPE = '",TRIM(piota_type),"'"
      WRITE (iunit,'(a,(1p,4e22.14))') '  AI = ',(ai(n-1), n=1,SIZE(ai))
      i	= minloc(ai_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AI_AUX_S = ',
     1         (ai_aux_s(n), n=1,i)
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AI_AUX_F = ',
     1         (ai_aux_f(n), n=1,i)
      END IF
      WRITE (iunit, '(2x,3a)') "PCURR_TYPE = '",TRIM(pcurr_type),"'"
      WRITE (iunit,'(a,(1p,4ES22.12E3))')
     1 '  AC = ',(ac(n-1), n=1,SIZE(ac))
      i	= minloc(ac_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AC_AUX_S = ',
     1         (ac_aux_s(n), n=1,i)
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AC_AUX_F = ',
     1         (ac_aux_f(n), n=1,i)
      END IF
      WRITE(iunit,'(A)') '!----- Axis Parameters ----- '
      WRITE (iunit,'(a,(1p,4e22.14))')
     1     '  RAXIS_CC = ',(raxis_cc(n), n=0,ntor)
      IF (lasym)
     1      WRITE (iunit,'(a,(1p,4ES22.12E3))')
     2     '  RAXIS_CS = ',(raxis_cs(n), n=0,ntor)
      IF (lasym)
     1      WRITE (iunit,'(a,(1p,4ES22.12E3))')
     2     '  ZAXIS_CC = ',(zaxis_cc(n), n=0,ntor)
      WRITE (iunit,'(a,(1p,4ES22.12E3))')
     1     '  ZAXIS_CS = ',(zaxis_cs(n), n=0,ntor)
      WRITE(iunit,'(A)') '!----- Boundary Parameters -----'
      DO m = 0, mpol - 1
         DO n = -ntor, ntor
            IF ((rbc(n,m).ne.0) .or. (zbs(n,m).ne.0)) THEN
               WRITE(iunit,'(2(A,I4.3,A,I3.3,A,ES22.12E3))')
     1         '  RBC(',n,',',m,') = ',rbc(n,m),
     2         '    ZBS(',n,',',m,') = ',zbs(n,m)
               IF (.not. lasym) CYCLE
               WRITE(iunit,'(2(A,I4.3,A,I3.3,A,ES22.12E3))')
     1         '  RBS(',n,',',m,') = ',rbs(n,m),
     2         '    ZBC(',n,',',m,') = ',zbc(n,m)
            END IF
         END DO
      END DO
      WRITE(iunit,'(A)') '/'
      RETURN
      END SUBROUTINE write_indata_namelist

      END MODULE vmec_input


