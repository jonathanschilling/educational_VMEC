!> \file
!> \brief Read the input file.

!> \brief Read the input file.
!>
!> @param input_file input file to read from
!> @param ier_flag error flag
  SUBROUTINE readin(input_file, ier_flag)
  USE vmec_main
  USE vmec_params
  USE mgrid_mod, ONLY: nextcur, curlabel, read_mgrid, &
                       nr0b, np0b, nz0b, rminb, zminb, rmaxb, zmaxb
  IMPLICIT NONE

  INTEGER, INTENT(inout) :: ier_flag
  CHARACTER(LEN=*) :: input_file

  INTEGER :: n, iunit, ier_flag_init, i, ni, m, nsmin, mj, isgn, ioff, joff
  REAL(rprec), DIMENSION(:,:), POINTER :: rbcc, rbss, rbcs, rbsc, zbcs, zbsc, zbcc, zbss
  REAL(rprec) :: rtest, ztest, tzc, trc, delta
  REAL(rprec), ALLOCATABLE :: temp(:)
  CHARACTER(LEN=100) :: line, line2

  ier_flag_init = ier_flag
  ier_flag      = norm_term_flag

  ! READ IN DATA FROM INDATA FILE
  CALL read_indata(input_file, iunit, ier_flag)
  IF (ier_flag .ne. norm_term_flag) then
     RETURN
  end if

  ! Open output files here, print out heading to threed1 file
  CALL heading(input_extension)

  IF (lfreeb) THEN
     ! READ IN AND STORE (FOR SEQUENTIAL RUNNING) MAGNETIC FIELD DATA FROM MGRID_FILE
     CALL read_mgrid (mgrid_file, extcur, nzeta, nfp, .true., ier_flag)

     ! check again for lfreeb
     ! --> might have been reset to .false. if mgrid file was not found
     IF (lfreeb) THEN
        IF (ier_flag .ne. norm_term_flag) then
           RETURN
        end if

        WRITE (nthreed,20) nr0b, nz0b, np0b, rminb, rmaxb, zminb, zmaxb, TRIM(mgrid_file)
20 FORMAT(//,' VACUUM FIELD PARAMETERS:',/,1x,24('-'),/,       &
             '  nr-grid  nz-grid  np-grid      rmin      rmax      zmin', &
             '      zmax     input-file',/,3i9,4f10.3,5x,a)
     END IF
  END IF

  ! PARSE NS_ARRAY
  multi_ns_grid = 1
  nsmin = 1
  DO WHILE (ns_array(multi_ns_grid) .ge. nsmin .and. & ! increasing or constant ns
            multi_ns_grid .lt. 100)                    ! max 100 multi-grid iterations
     nsmin = MAX(nsmin, ns_array(multi_ns_grid))
     IF (nsmin .le. nsd) THEN ! max ns is given by nsd
        multi_ns_grid = multi_ns_grid + 1
     ELSE
        ! Optimizer, Boozer code overflows otherwise
        ns_array(multi_ns_grid) = nsd
        nsmin = nsd
        PRINT *,' NS_ARRAY ELEMENTS CANNOT EXCEED ',nsd
        PRINT *,' CHANGING NS_ARRAY(',multi_ns_grid,') to ', nsd
     END IF
  END DO
  multi_ns_grid = multi_ns_grid - 1

   ! consistency check on requested number of flux surfaces
  IF (ALL(ns_array.eq.0)) THEN
     stop 'NS ARRAY MUST NOT BE ALL ZEROES' ! was ier_flag=ns_error_flag
  END IF

  ns_maxval = nsmin

  ! WRITE OUT DATA TO THREED1 FILE

  !SPH121912 - SCALING TO RENDER LAMSCALE=1
  !      delta = twopi/phiedge    !phiedge=>twopi
  !      phiedge = phiedge*delta
  !      bcrit = bcrit*delta
  !      curtor = curtor*delta
  !      extcur = extcur*delta
  !      am = am*delta**2

  WRITE (nthreed,100)                                               &
    ns_array(multi_ns_grid),ntheta1,nzeta,mpol,ntor,nfp,            &
    gamma,spres_ped,phiedge,curtor
100 FORMAT(/,' COMPUTATION PARAMETERS: (u = theta, v = zeta)'/,       &
             1x,45('-'),/,                                                   &
             '     ns     nu     nv     mu     mv',/,5i7//,                  &
             ' CONFIGURATION PARAMETERS:',/,1x,39('-'),/,                    &
             '    nfp      gamma      spres_ped    phiedge(wb)     curtor(A)',/, &
             i7,1p,e11.3,2e15.3,e14.3/)

  IF (nvacskip.le.0) then
     ! default value for nvacskip: number of field periods ... ?
     ! how does that work for e.g. LHD with nfp=10 ???
     nvacskip = nfp
  end if
  WRITE (nthreed,110) ncurr,niter_array(multi_ns_grid),ns_array(1), &
    nstep,nvacskip,                                                 &
    ftol_array(multi_ns_grid),tcon0,lasym,lconm1,                   &
    mfilter_fbdy,nfilter_fbdy
110 FORMAT(' RUN CONTROL PARAMETERS:',/,1x,23('-'),/,                   &
           '  ncurr  niter   nsin  nstep  nvacskip      ftol     tcon0',&
           '    lasym   lconm1',/,                                      &
           4i7,i10,1p,2e10.2,2L9,/,                                     &
           '  mfilter_fbdy nfilter_fbdy',/,                             &
           2(6x,i7),/)

  IF (nextcur .gt. 0) THEN
     WRITE(nthreed, "(' EXTERNAL CURRENTS',/,1x,17('-'))")
     ni = 0
     IF (ALLOCATED(curlabel)) ni = MAXVAL(LEN_TRIM(curlabel(1:nextcur)))
     ni = MAX(ni+4, 14)
     WRITE (line,  '(a,i2.2,a)') "(5a",ni,")"
     WRITE (line2, '(a,i2.2,a)') "(5(",ni-12,"x,1p,e12.4))"
     DO i = 1,nextcur,5
        ni = MIN(i+4, nextcur)
        IF (ALLOCATED(curlabel)) then
           WRITE (nthreed, line, iostat=mj) (TRIM(curlabel(n)),n=i,ni)
        end if
        WRITE (nthreed, line2,iostat=mj) (extcur(n), n=i,ni)
     ENDDO
     WRITE (nthreed, *)
  ENDIF

  IF (bloat .ne. one) THEN
      WRITE (nthreed,'(" Profile Bloat Factor: ",1pe11.4)') bloat
      phiedge = phiedge*bloat
  ENDIF

  IF (pres_scale .ne. one) THEN
      WRITE (nthreed,'(" Pressure profile factor: ",1pe11.4,'//        &
            '" (multiplier for pressure)")') pres_scale
  END IF

  ! Print out am array
  WRITE(nthreed,130)
  WRITE(nthreed,131) TRIM(pmass_type)
  WRITE(nthreed,132)
130 FORMAT(' MASS PROFILE COEFFICIENTS am - newton/m**2 (EXPANSION IN NORMALIZED RADIUS):')
131 FORMAT(' PMASS parameterization type is ''', a,'''')
132 FORMAT(1x,35('-'))
  WRITE(nthreed,135)(am(i-1),i=1, SIZE(am))
135 FORMAT(1p,6e12.3)
  SELECT CASE(TRIM(pmass_type))
  CASE ('Akima_spline','cubic_spline')
     WRITE(nthreed,"(' am_aux_s is' )")
     WRITE(nthreed,135)(am_aux_s(i),i=1, SIZE(am_aux_s))
     WRITE(nthreed,"(' am_aux_f is' )")
     WRITE(nthreed,135)(am_aux_f(i),i=1, SIZE(am_aux_f))
  END SELECT

  IF (ncurr.eq.0) THEN
      WRITE(nthreed,140)
      WRITE(nthreed,143) TRIM(piota_type)
      WRITE(nthreed,132)
      !  Print out ai array
      WRITE(nthreed,135)(ai(i-1),i=1, SIZE(ai))
      SELECT CASE(TRIM(piota_type))
      CASE ('Akima_spline','cubic_spline')
         WRITE(nthreed,"(' ai_aux_s is' )")
         WRITE(nthreed,135)(ai_aux_s(i),i=1, SIZE(ai_aux_s))
         WRITE(nthreed,"(' ai_aux_f is' )")
         WRITE(nthreed,135)(ai_aux_f(i),i=1, SIZE(ai_aux_f))
      END SELECT
  ELSE
      !  Print out ac array
      WRITE(nthreed,145)
      WRITE(nthreed,146) TRIM(pcurr_type)
      WRITE(nthreed,147)
      WRITE(nthreed,135)(ac(i-1),i=1, SIZE(ac))
      SELECT CASE(TRIM(pcurr_type))
      CASE ('Akima_spline_Ip','Akima_spline_I',                            &
            'cubic_spline_Ip','cubic_spline_I')
         WRITE(nthreed,"(' ac_aux_s is' )")
         WRITE(nthreed,135)(ac_aux_s(i),i=1, SIZE(ac_aux_s))
         WRITE(nthreed,"(' ac_aux_f is' )")
         WRITE(nthreed,135)(ac_aux_f(i),i=1, SIZE(ac_aux_f))
      END SELECT
  ENDIF

140 FORMAT(/' IOTA PROFILE COEFFICIENTS ai (EXPANSION IN NORMALIZED RADIUS):')
143 FORMAT(' PIOTA parameterization type is ''', a,'''')
145 FORMAT(/' TOROIDAL CURRENT DENSITY (*V'') COEFFICIENTS ac (EXPANSION IN NORMALIZED RADIUS):')
146 FORMAT(' PCURR parameterization type is ''', a,'''')
147 FORMAT(1x,38('-'))

 WRITE(nthreed,150)
 WRITE(nthreed,135)(aphi(i),i=1, SIZE(aphi))
150 FORMAT(/' NORMALIZED TOROIDAL FLUX COEFFICIENTS aphi (EXPANSION IN S):',/,1x,35('-'))

!  Fourier Boundary Coefficients (header written always, but data only if in free-boundary mode)
  WRITE(nthreed,180)
180 FORMAT(/,' R-Z FOURIER BOUNDARY COEFFICIENTS AND',                  &
             ' MAGNETIC AXIS INITIAL GUESS',/,                          &
             ' R = RBC*COS(m*u - n*v) + RBS*SIN(m*u - n*v),',           &
             ' Z = ZBC*COS(m*u - n*v) + ZBS*SIN(m*u-n*v)'/1x,86('-'),   &
           /,'   nb  mb     rbc         rbs         zbc         zbs   ',&
             '    raxis(c)    raxis(s)    zaxis(c)    zaxis(s)')

  IF (lasym) THEN
     ! CONVERT TO REPRESENTATION WITH RBS(m=1) = ZBC(m=1)
     delta = ATAN( (rbs(0,1) - zbc(0,1))/(ABS(rbc(0,1)) + ABS(zbs(0,1))) )
     IF (delta .ne. zero) THEN
       DO m = 0,mpol1
         DO n = -ntor,ntor
           trc = rbc(n,m)*COS(m*delta) + rbs(n,m)*SIN(m*delta)
           rbs(n,m) = rbs(n,m)*COS(m*delta) - rbc(n,m)*SIN(m*delta)
           rbc(n,m) = trc
           tzc = zbc(n,m)*COS(m*delta) + zbs(n,m)*SIN(m*delta)
           zbs(n,m) = zbs(n,m)*COS(m*delta) - zbc(n,m)*SIN(m*delta)
           zbc(n,m) = tzc
         ENDDO
       ENDDO
     ENDIF
  ENDIF

  ! ALLOCATE MEMORY FOR NU, NV, MPOL, NTOR SIZED ARRAYS
  CALL allocate_nunv

  ! CONVERT TO INTERNAL REPRESENTATION OF MODES
  !
  ! R =   RBCC*COS(M*U)*COS(N*V) + RBSS*SIN(M*U)*SIN(N*V)
  !     + RBCS*COS(M*U)*SIN(N*V) + RBSC*SIN(M*U)*COS(N*V)
  ! Z =   ZBCS*COS(M*U)*SIN(N*V) + ZBSC*SIN(M*U)*COS(N*V)
  !     + ZBCC*COS(M*U)*COS(N*V) + ZBSS*SIN(M*U)*SIN(N*V)
  !
  !
  ! POINTER ASSIGNMENTS (NOTE: INDICES START AT 1, NOT 0, FOR POINTERS, EVEN THOUGH
  !                      THEY START AT ZERO FOR RMN_BDY)
  ! ARRAY STACKING ORDER DETERMINED HERE
  rbcc => rmn_bdy(:,:,rcc)
  zbsc => zmn_bdy(:,:,zsc)
  IF (lthreed) THEN
     rbss => rmn_bdy(:,:,rss)
     zbcs => zmn_bdy(:,:,zcs)
  END IF

  IF (lasym) THEN
     rbsc => rmn_bdy(:,:,rsc)
     zbcc => zmn_bdy(:,:,zcc)
     IF (lthreed) THEN
        rbcs => rmn_bdy(:,:,rcs)
        zbss => zmn_bdy(:,:,zss)
     END IF
  ENDIF

  rmn_bdy = 0;  zmn_bdy = 0

  ioff = LBOUND(rbcc,1)
  joff = LBOUND(rbcc,2)

  ! go over all mode number combinations that could be specified in input file for given Fourier resolution
  DO m=0,mpol1

     mj = m+joff

     ! in free-boundary, skip initial guess with Fourier mode number m larger than mfilter_fbdy
     IF (lfreeb .and. (mfilter_fbdy.gt.1 .and. m.gt.mfilter_fbdy)) CYCLE

     DO n=-ntor,ntor

        ! in free-boundary, skip initial guess with Fourier mode number n larger than nfilter_fbdy
        IF (lfreeb .and. (nfilter_fbdy.gt.0 .and. ABS(n).gt.nfilter_fbdy)) CYCLE

        ni = ABS(n) + ioff

        ! rbc, zbs etc are defined for the full 2d Fourier kernel with possibly negative n
        ! this sign is used to convert it to positive-only Fourier numbers
        IF (n .eq. 0) THEN
           isgn = 0
        ELSE IF (n .gt. 0) THEN
           isgn = 1
        ELSE
           isgn = -1
        END IF

        rbcc(ni,mj) = rbcc(ni,mj) + rbc(n,m)
        IF (m .gt. 0) zbsc(ni,mj) = zbsc(ni,mj) + zbs(n,m)

        IF (lthreed) THEN
           IF (m .gt. 0) rbss(ni,mj) = rbss(ni,mj) + isgn*rbc(n,m)
           zbcs(ni,mj) = zbcs(ni,mj) - isgn*zbs(n,m)
        END IF

        IF (lasym) THEN
           IF (m .gt. 0) rbsc(ni,mj) = rbsc(ni,mj) + rbs(n,m)
           zbcc(ni,mj) = zbcc(ni,mj) + zbc(n,m)

           IF (lthreed) THEN
              rbcs(ni,mj) = rbcs(ni,mj) - isgn*rbs(n,m)
              IF (m .gt. 0) zbss(ni,mj) = zbss(ni,mj) + isgn*zbc(n,m)
           END IF
        END IF

        ! TODO: why exit only here or is this leftover from something else?
        IF (ier_flag_init .ne. norm_term_flag) then
           print *, "skip printing of initial guess for geometry, since ier_flag_init = ",ier_flag_init
           CYCLE
        end if

        ! compute sum of quantities to be printed below
        trc = ABS(rbc(n,m)) + ABS(rbs(n,m)) + ABS(zbc(n,m)) + ABS(zbs(n,m))

        IF (m .eq. 0) THEN
           ! m=0: also axis coefficients can be printed

           ! above n-loop is re-used for printing:
           ! negative n come before positive n, so any |n| will have been completely assigned also its negative counterparts at this point
           IF (n .lt. 0) CYCLE

           IF (trc.eq.zero .and. ABS(raxis_cc(n)).eq.zero .and. ABS(zaxis_cs(n)).eq.zero) then
              ! all values to be printed are zero, so skip printing them
              CYCLE
           end if

           WRITE (nthreed,195) n, m, rbc(n,m), rbs(n,m), zbc(n,m), zbs(n,m), &
                               raxis_cc(n), raxis_cs(n), zaxis_cc(n), zaxis_cs(n)
        ELSE
           ! m > 0, so only print boundary coeffs (axis does not have coeffs with m>1

           IF (trc .eq. zero) then
              ! all values to be printed are zero, so skip printing them
              CYCLE
           end if
           WRITE (nthreed,195) n, m, rbc(n,m), rbs(n,m), zbc(n,m), zbs(n,m)
        END IF

     END DO
  END DO
195 FORMAT(i5,i4,1p,8e12.4)

  ! CHECK SIGN OF JACOBIAN (SHOULD BE SAME AS SIGNGS)
  m = 1
  mj = m+joff
  rtest = SUM(rbcc(1:ntor1,mj))
  ztest = SUM(zbsc(1:ntor1,mj))
  lflip=(rtest*ztest .lt. zero)
  signgs = -1
  IF (lflip) CALL flip_theta(rmn_bdy, zmn_bdy)

  ! CONVERT TO INTERNAL FORM FOR (CONSTRAINED) m=1 MODES
  ! INTERNALLY, FOR m=1: XC(rss) = .5(RSS+ZCS), XC(zcs) = .5(RSS-ZCS)
  ! WITH XC(zcs) -> 0 FOR POLAR CONSTRAINT ! TODO: jons: maybe XC(zsc) -> 0 is meant here ?
  ! FOR ASYMMETRIC CASE, XC(rsc) = .5(RSC+ZCC), XC(zcc) = .5(RSC-ZCC)
  ! WITH XC(zss) -> 0 FOR POLAR CONSTRAINT
  ! (see convert_sym, convert_asym in totzsp.f90 file)
  IF (lconm1 .AND. (lthreed.OR.lasym)) THEN
     ALLOCATE (temp(SIZE(rbcc,1)))
     IF (lthreed) THEN
        mj = 1+joff ! index of n=0, m=1 modes
        temp = rbss(:,mj)
        rbss(:,mj) = cp5*(temp(:) + zbcs(:,mj))
        zbcs(:,mj) = cp5*(temp(:) - zbcs(:,mj))
     END IF
     IF (lasym) THEN
        mj = 1+joff ! index of n=0, m=1 modes
        temp = rbsc(:,mj)
        rbsc(:,mj) = cp5*(temp(:) + zbcc(:,mj))
        zbcc(:,mj) = cp5*(temp(:) - zbcc(:,mj))
     END IF
     DEALLOCATE (temp)
  END IF

  ! Convert to Internal units
  currv = mu0*curtor

END SUBROUTINE readin
