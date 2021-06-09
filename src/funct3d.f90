!> \file
!> \brief Evaluate the three-dimensional MHD energy functional.
!>        Think of this as the "forward model" that
!>        tells you the MHD forces in Fourier space
!>        given the Fourier coefficients of the flux surface geometry.

!> \brief Evaluate the three-dimensional MHD energy functional.
!>        Think of this as the "forward model" that
!>        tells you the MHD forces in Fourier space
!>        given the Fourier coefficients of the flux surface geometry.
!>
!> @param ier_flag error flag
SUBROUTINE funct3d (ier_flag)

  USE vmec_main
  USE vacmod, ONLY: bsqvac, amatsav, bvecsav, mnpd2, bsubvvac
  use nestor_io, only: write_nestor_outputs
  USE vmec_params, ONLY: bad_jacobian_flag, signgs
  USE realspace
  USE vforces
  USE xstuff
  USE vparams, ONLY: twopi

  IMPLICIT NONE

  INTEGER, INTENT(inout) :: ier_flag

  INTEGER :: l0pi, l, lk, ivacskip, js, ku, m, n
  INTEGER :: nvskip0 = 0
  REAL(dp), DIMENSION(mnmax) :: rmnc, zmns, lmns, rmns, zmnc, lmnc
  REAL(dp), DIMENSION(:), POINTER :: lu, lv
  REAL(dp) :: presf_ns, delt0
  REAL(dp), EXTERNAL :: pmass

  CHARACTER(LEN=255) :: vac_file
  integer :: nvac, istat_vac

  character(len=255) :: nestor_cmd

!  character(len=*), parameter :: nestor_executable = &
!    "/home/IPP-HGW/jons/work/code/educational_VMEC/build/bin/xnestor"

!  character(len=*), parameter :: nestor_executable = &
!     "/data2/jonathan/work/code/educational_VMEC/build/bin/xnestor"

!    character(len=*), parameter :: nestor_executable = &
!     "python /home/jonathan/work/code/NESTOR/src/main/python/NESTOR.py"

   character(len=*), parameter :: nestor_executable = &
    "python /home/jonathan/work/code/NESTOR/src/main/python/ooNESTOR.py"

!   character(len=*), parameter :: nestor_executable = &
!     "python3 /home/IPP-HGW/jons/work/code/NESTOR/src/main/python/NESTOR.py"

  !> use system call to stand-alone NESTOR for vacuum computation
  logical :: lexternal_nestor = .false.

  !> dump reference input for and output of NESTOR when using internal NESTOR
  logical :: ldump_vacuum_ref = .false.

  character(len=255) :: dump_filename
  logical            :: dump_geometry = .false.
  logical            :: dump_constraint_force = .false.

  ! POINTER ALIASES
  lu => czmn
  lv => crmn

  ! CONVERT ODD M TO 1/SQRT(S) INTERNAL REPRESENTATION
  ! temprary use of gc (force) for scaled xc (position)
  gc(:neqs) = xc(:neqs)*scalxc(:neqs)

  ! INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
  ! R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
  ! FIRST, DO SYMMETRIC [ F(u,v) = F(-u,-v) ] PIECES
  ! ON THE RANGE u = 0,pi  and v = 0,2*pi
  CALL totzsps (gc, r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon)

  IF (lasym) THEN
     ! ANTI-SYMMETRIC CONTRIBUTIONS TO INVERSE TRANSFORMS
     CALL totzspa (gc, armn,   brmn, extra3, &
                       azmn,   bzmn, extra4, &
                       blmn,   clmn,         &
                     extra1, extra2           )

     ! SUM SYMMETRIC, ANTISYMMETRIC PIECES APPROPRIATELY
     ! TO GET R, Z, L, (AND RCON, ZCON) ON FULL RANGE OF u (0 to 2*pi)
     CALL symrzl (  r1,   ru,     rv,   z1,   zu,     zv,   lu,   lv,   rcon,   zcon, &
                  armn, brmn, extra3, azmn, bzmn, extra4, blmn, clmn, extra1, extra2    )
  ENDIF

  if (dump_geometry .and. iter2 .le. 2) then
      write(dump_filename, 999) ns, iter2, trim(input_extension)
      open(unit=42, file=trim(dump_filename), status="unknown")

      write(42, *) "# ns ntheta3 nzeta"
      write(42, *) ns, ntheta3, nzeta

      write(42, *) "# js ku lv m%2 r1 ru rv z1 zu zv lu lv rcon zcon"
      l = 1
      DO ku = 1, ntheta3
        DO lk = 1, nzeta
          DO js = 1, ns
            DO m=0,1
              !l = ((ku-1)*nzeta+(lk-1))*ns+js
              write (42, *) js, ku, lk, m, &
                            r1(l,m), ru(l,m), rv(l,m),  &
                            z1(l,m), zu(l,m), zv(l,m),  &
                            lu(l+m*nrzt), lv(l+m*nrzt), &
                            rcon(l,m), zcon(l,m)
            end do
            l = l+1
          end do
        end do
      end do

      close(42)

      print *, "dumped geometry output to '"//trim(dump_filename)//"'"
      !stop
  end if
999 format('funct3d_geometry_',i5.5,'_',i6.6,'.',a)

  ! now that we have the current real-space geometry, use the opportunity so do some analysis / statistics
  ! --> at ns, sqrt(s)==1 --> no need to scale r1(ns or l0pi, 1) by sqrt(s) anymore

  router = r1(ns,0) + r1(ns,1) ! index ns corresponds to (u=0, v=0, js=ns)

  ! l0pi is the index corresponding to (u = pi, v = 0, js = ns) (?)
  l0pi = ns*(1 + nzeta*(ntheta2 - 1))
  rinner = r1(l0pi,0) + r1(l0pi,1)

  r00 = r1(1,0) ! contrib from only even m
  z00 = z1(1,0) ! contrib from only even m



  ! COMPUTE CONSTRAINT RCON, ZCON

  ! --> see Hirshman, Schwenn & Nührenberg: summation of even-m and odd-m*sqrt(s)

! #ifndef _HBANGLE
  rcon(:nrzt,0) = rcon(:nrzt,0) + rcon(:nrzt,1)*sqrts(:nrzt) ! odd-m entries need to be scaled appropriately
  zcon(:nrzt,0) = zcon(:nrzt,0) + zcon(:nrzt,1)*sqrts(:nrzt)
! #end /* ndef _HBANGLE */

  ! Assemble dR/du and dZ/du, since they are needed for arnorm, aznorm in bcovar() and for guess_axis().
  ! Must store them in separate arrays ru0, zu0 since separation into even-m and odd-m
  ! must be kept for e.g jacobian and metric elements (I guess...).
  ru0(:nrzt)    = ru(:nrzt,0)   + ru(:nrzt,1)*sqrts(:nrzt)
  zu0(:nrzt)    = zu(:nrzt,0)   + zu(:nrzt,1)*sqrts(:nrzt)

  IF (iter2.eq.iter1 .and. ivac.le.0) THEN

     print *, "rcon0 <-- (rcon * s) into volume"
     ! iter2 == iter1 is true at start of a new multi-grid iteration
     ! ivac .le. 0 is always true for fixed-boundary,
     ! but only true for first iteration in free-boundary (?)

     ! COMPUTE RCON0, ZCON0 FOR FIXED BOUNDARY BY SCALING EDGE VALUES
     ! SCALE BY POWER OF SQRTS, RATHER THAN USE rcon0 = rcon, etc.
     ! THIS PREVENTS A DISCONTINUITY WHEN RESTARTING FIXED BOUNDARY WITH NEW RCON0....
     !
     ! NOTE: IN ORDER TO MAKE INITIAL CONSTRAINT FORCES SAME FOR FREE/FIXED
     ! BOUNDARY, WE SET RCON0,ZCON0 THE SAME INITIALLY, BUT TURN THEM OFF
     ! SLOWLY IN FREE-BOUNDARY VACUUM LOOP (BELOW)
     DO l = 1, ns
        ! value of rcon(ns) is scaled into the volume proportional to s
        rcon0(l:nrzt:ns) = rcon(ns:nrzt:ns,0) * sqrts(l:nrzt:ns)**2
        zcon0(l:nrzt:ns) = zcon(ns:nrzt:ns,0) * sqrts(l:nrzt:ns)**2
     END DO
  ENDIF
! #end /* ndef _HBANGLE */

  ! COMPUTE S AND THETA DERIVATIVE OF R AND Z AND JACOBIAN ON HALF-GRID
  CALL jacobian
  IF (first.eq.2 .and. iequi.eq.0) then
     ! bad jacobian and not final iteration yet (would be indicated by iequi.eq.1) --> need to restart
     ! except when computing output file --> ignore bad jacobian
     return
  end if



  ! NOTE: up to here, only worked on geometry so far...



  ! COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC AND KINETIC PRESSURE,
  ! AND METRIC ELEMENTS ON HALF-GRID
  CALL bcovar (lu, lv)

  ! COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
  ! NOTE: FOR FREE BOUNDARY RUNS, THE VALUE OF RBTOR=R*BTOR AT THE PLASMA EDGE
  ! SHOULD BE ADJUSTED TO APPROXIMATELY EQUAL THE VACUUM VALUE.
  ! THIS CAN BE DONE BY CHANGING EITHER PHIEDGE OR THE INITIAL CROSS SECTION
  ! ACCORDING TO THE SCALING LAW  R*BTOR .EQ. PHIEDGE/(R1 * Z1).
  IF (lfreeb .and. iter2.gt.1 .and. iequi.eq.0) THEN

     IF ((fsqr + fsqz) .le. 1.e-3_dp) then
        !print *, "force residuals decreased sufficiently => increment ivac=",ivac
        ! this does ivac=-1 --> ivac=0 to enable NESTOR at all

        ! initially, ivac is initialized to -1 by reset_params
        ! when R,Z forces are <1e-3, enable vacuum contribution
        ivac = ivac+1   ! decreased from e-1 to e-3 - sph12/04
        ! I guess this is where the vacuum pressure suddenly gets turned on ?
     end if

     IF (nvskip0 .eq. 0) then
        nvskip0 = MAX(1, nvacskip)
     end if

     IF (ivac .ge. 0) THEN
        ! IF INITIALLY ON, MUST TURN OFF rcon0, zcon0 SLOWLY
        rcon0 = 0.9_dp*rcon0
        zcon0 = 0.9_dp*zcon0

        ivacskip = MOD(iter2-iter1, nvacskip)
        IF (ivac .le. 2) then
           ivacskip = 0
           ! vacuum pressure not turned on yet (?)
           ! and do full vacuum calc on every iteration
        end if

        ! EXTEND NVACSKIP AS EQUILIBRIUM CONVERGES
        IF (ivacskip .eq. 0) THEN
           nvacskip = one/MAX(1.e-1_dp, 1.e11_dp*(fsqr+fsqz))
           nvacskip = MAX(nvacskip, nvskip0)
        END IF

        ! NOTE: gc contains correct edge values of r,z,l arrays
        ! convert_sym, convert_asym have been applied to m=1 modes
        CALL convert (rmnc, zmns, lmns, rmns, zmnc, lmnc, gc, ns)

        ! raxis_nestor(1:nzeta) = r1(1:ns*nzeta:ns,0)
        ! zaxis_nestor(1:nzeta) = z1(1:ns*nzeta:ns,0)

        if (ldump_vacuum_ref .or. lexternal_nestor) then

           ! build filename for NESTOR inputs
           if (ldump_vacuum_ref) then
             write(vac_file, "(A,I6.6,A)") "vac_ref/vacin_"//TRIM(input_extension)//"_", &
                vacuum_calls, ".nc"
           else ! lexternal_nestor
             write(vac_file, "(A,I6.6,A)") "vac/vacin_"//TRIM(input_extension)//"_", &
                vacuum_calls, ".nc"
           end if

           ! write NESTOR inputs
           call write_nestor_inputs(trim(vac_file),                                &
                  vacuum_calls, ier_flag, trim(mgrid_file), trim(input_extension), &
                  ivacskip, ivac, nfp, ntor, mpol, nzeta, ntheta,                  &
                  mnmax, xm, xn, rmnc, zmns, rmns, zmnc,                           &
                  rbtor, ctor, lasym, signgs, extcur,                              &
                  r1(1:ns*nzeta:ns,0), z1(1:ns*nzeta:ns,0), wint(ns:nznt*ns:ns), nznt, &
                  amatsav, bvecsav, mnpd2, bsubvvac)

           ! print *, "dumped NESTOR inputs to '"//trim(vac_file)//"'"
        end if

        if (.not. lexternal_nestor) then
           ! use internal NESTOR
           CALL vacuum (rmnc, rmns, zmns, zmnc, xm, xn,                                    &
                        ctor, rbtor, wint(ns:nznt*ns:ns), ivacskip, ivac, mnmax, ier_flag, &
                        lasym, signgs, r1(1:ns*nzeta:ns,0), z1(1:ns*nzeta:ns,0))

        else
           ! construct command with argument for stand-alone external NESTOR
           write(nestor_cmd, "(A,X,A)") trim(nestor_executable), trim(vac_file)

           ! do system call to external NESTOR
           !print *, "NESTOR command: '",trim(nestor_cmd),"'"
           call system(nestor_cmd)

           !print *, "system call to NESTOR finished"
        end if


        if (ldump_vacuum_ref) then
           ! construct filename for reference NESTOR output
           write(vac_file, "(A,I6.6,A)") "vac_ref/vacout_ref_"//TRIM(input_extension)//"_", &
              vacuum_calls, ".nc"

           call write_nestor_outputs(vac_file, lasym, ivac, ier_flag)

           ! print *, "dumped NESTOR outputs to '"//trim(vac_file)//"'"

        else if (lexternal_nestor) then
          ! contruct filename from which to read output of stand-alone NESTOR
          write(vac_file, "(A,I6.6,A)") "vac/vacout_"//TRIM(input_extension)//"_", &
             vacuum_calls, ".nc"

          !print *, "read NESTOR output from '"//trim(vac_file)//"'"

          ! read output of external NESTOR
          call read_nestor_outputs(trim(vac_file), ier_flag, ivac)
          !print *, "read in NESTOR output: ivac=",ivac," ier_flag=",ier_flag

        end if

        ! update counter for calls to NESTOR (initialized to 0 in reset_params)
        vacuum_calls = vacuum_calls + 1


        IF (ier_flag .ne. 0) then
           ! some error occured within NESTOR, so cancel the iterations
           return
        end if

        ! RESET FIRST TIME FOR SOFT START
        IF (ivac .eq. 1) THEN
           first = 2

           ! delt0 is never used --> ignore change to time step by restart_iter
           ! TODO: Since delt0 is never used and also nothing else is modified in restart_iter,
           ! why bother which value it has on entry to restart_iter?
           delt0 = delt
           CALL restart_iter(delt0)
           first = 1 ! already done in restart_iter for first.eq.2
        END IF

        ! IN CASE PRESSURE IS NOT ZERO AT EXTRAPOLATED EDGE...
        ! UNCOMMENT ALL "RPRES" COMMENTS HERE AND IN BCOVAR, FORCES ROUTINES
        ! IF NON-VARIATIONAL FORCES ARE DESIRED
        !
        ! presf_ns = 1.5_dp*pres(ns) - 0.5_dp*pres(ns1)
        ! MUST NOT BREAK TRI-DIAGONAL RADIAL COUPLING: OFFENDS PRECONDITIONER!
        presf_ns = pmass(hs*(ns-1.5_dp))
        IF (presf_ns .ne. zero) then
           presf_ns = (pmass(1._dp)/presf_ns) * pres(ns)
        end if

        lk = 0
        DO l = ns, nrzt, ns
           lk = lk + 1

           ! current extrapolation to LCFS of plasma magnetic field
           bsqsav(lk,3) = 1.5_dp*bzmn_o(l) - 0.5_dp*bzmn_o(l-1)

           ! total pressure (?) at LCFS
           ! (gcon(l) is only used as a temporary variable here,
           !  since it immediately gets overwritten when entering alias())
           gcon(l)      = bsqvac(lk) + presf_ns

           ! edge force contribution (see forces())
           ! --> *HERE* is where the free-boundary computation enters VMEC !
           rbsq(lk) = gcon(l)*(r1(l,0) + r1(l,1))*ohs

           ! residual magnetic field discontinuity at LCFS
           ! --> used in last column of printout (as flux surface avg.; relative to <bsqsav>)
           dbsq(lk) = ABS(gcon(l)-bsqsav(lk,3))
        END DO

        IF (ivac .eq. 1) THEN
           ! first time: debugging of initial soft-start ????
           print *,"bsqsav(:,1:2) are filled now"
           bsqsav(:nznt,1) = bzmn_o(ns:nrzt:ns) ! initial magnetic field at boundary
           bsqsav(:nznt,2) = bsqvac(:nznt)      ! initial NESTOR |B|^2 at boundary
        ENDIF
     ENDIF
  ENDIF

  IF (iequi .NE. 1) THEN
     ! normal iterations, not final call from fileout (which sets iequi=1)

! #ifndef _HBANGLE
     ! COMPUTE CONSTRAINT FORCE
     extra1(:nrzt,0) = (rcon(:nrzt,0) - rcon0(:nrzt))*ru0(:nrzt) &
                     + (zcon(:nrzt,0) - zcon0(:nrzt))*zu0(:nrzt)
     ! Fourier-space filter: only retain m=1, ..., (mpol1-1)==mpol-2 in gcon
     CALL alias (gcon, extra1(:,0), gc, gc(1+mns), gc(1+2*mns), extra1(:,1)) ! temporary re-use of extra1(:,1) for g_ss
! #end /* ndef _HBANGLE */

     if (dump_constraint_force .and. iter2.le.2) then
       write(dump_filename, 998) ns, iter2, trim(input_extension)
998 format('constraint_force_',i5.5,'_',i6.6,'.',a)

       open(unit=42, file=trim(dump_filename), status="unknown")

       write(42, *) "# ns ntheta3 nzeta mpol1 ntor"
       write(42, *) ns, ntheta3, nzeta, mpol1, ntor

       write(42, *) "# js ku lv extra1(:,0) gcons"
       DO js = 1, ns
         DO lk = 1, nzeta
           DO ku = 1, ntheta3
             l = ((ku-1)*nzeta+(lk-1))*ns+js
             write (42, *) js, ku, lk, &
                           extra1(l,0), gcon(l)
           end do
         end do
       end do

       write(42, *) "# js n m gcs gsc gcc gss"
       l=0
       DO m = 0, mpol1
         DO n = 0, ntor
           DO js = 1, ns
             l = l+1
             write (42, *) js, n, m, &
                           gc(l), gc(mns+l), gc(2*mns+l), extra1(l,1)
           end do
         end do
       end do

       close(42)

       print *, "dumped constraint force to '"//trim(dump_filename)//"'"
!        stop
     end if

     ! COMPUTE MHD FORCES ON INTEGER-MESH
     CALL forces

     ! SYMMETRIZE FORCES (in u-v space): NOTE - gc IS SMALL BY FACTOR 2 IF lasym=T
     IF (lasym) THEN
        CALL symforce (armn, brmn, crmn, azmn, bzmn, czmn,   blmn,   clmn,   rcon,   zcon, &
                         r1,   ru,   rv,   z1,   zu,   zv, extra3, extra4, extra1, extra2   )

        ! NOT NECESSARY (EVEN THOUGH CORRECT) --> why?
        ! gc = 2*gc
     END IF

     ! FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
     ! Note that gc is immediately zeroed on entry,
     ! so above use of gc is only for temporary values inside alias() !
     CALL tomnsps (gc,               &
                   armn, brmn, crmn, &
                   azmn, bzmn, czmn, &
                   blmn, clmn, rcon, zcon)
     IF (lasym) then
        CALL tomnspa (gc,             &
                      r1, ru, rv,     &
                      z1, zu, zv,     &
                      extra3, extra4, extra1, extra2)
     end if

     ! COMPUTE FORCE RESIDUALS (RAW AND PRECONDITIONED)
     gc = gc * scalxc    !!IS THIS CORRECT: SPH010214?
     CALL residue (gc, gc(1+irzloff), gc(1+2*irzloff))

     IF (iter2.eq.1 .and. (fsqr+fsqz+fsql).gt.1.E2_dp) then
         ! first iteration and gigantic force residuals --> what is going one here?
         first = 4 ! fatal error
     end if

!  ELSE
     ! iequi == 1 --> skip above remainder of funct3d
  END IF

END SUBROUTINE funct3d
