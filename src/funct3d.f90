!> \file
SUBROUTINE funct3d (ier_flag)
  USE vmec_main
  USE vacmod, ONLY: bsqvac, raxis_nestor, zaxis_nestor
  USE vmec_params, ONLY: bad_jacobian_flag, signgs
  USE realspace
  USE vforces
  USE xstuff
  USE vparams, ONLY: twopi
  IMPLICIT NONE

  INTEGER, INTENT(inout) :: ier_flag

  INTEGER :: l0pi, l, lk, ivacskip
  INTEGER :: nvskip0 = 0
  REAL(dp), DIMENSION(mnmax) :: rmnc, zmns, lmns, rmns, zmnc, lmnc
  REAL(dp), DIMENSION(:), POINTER :: lu, lv
  REAL(dp) :: presf_ns, delt0
  REAL(dp), EXTERNAL :: pmass

  ! POINTER ALIASES
  lu => czmn
  lv => crmn

  ! CONVERT ODD M TO 1/SQRT(S) INTERNAL REPRESENTATION
  gc(:neqs) = xc(:neqs)*scalxc(:neqs)

  ! INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
  ! R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
  ! FIRST, DO SYMMETRIC [ F(u,v) = F(-u,-v) ] PIECES
  ! ON THE RANGE u = 0,pi  and v = 0,2*pi
  CALL totzsps (gc, r1, ru, rv, z1, zu, zv, lu, lv, rcon, zcon)

  IF (lasym) THEN
     ! ANTI-SYMMETRIC CONTRIBUTIONS TO INVERSE TRANSFORMS
     CALL totzspa (gc, armn, brmn, extra3, azmn, bzmn, extra4,      &
                   blmn, clmn, extra1, extra2)

     ! SUM SYMMETRIC, ANTISYMMETRIC PIECES APPROPRIATELY
     ! TO GET R, Z, L, (AND RCON, ZCON) ON FULL RANGE OF u (0 to 2*pi)
     CALL symrzl (  r1,   ru,     rv,   z1,   zu,     zv,   lu,   lv,   rcon,   zcon, &
                  armn, brmn, extra3, azmn, bzmn, extra4, blmn, clmn, extra1, extra2    )
  ENDIF

  ! u = pi, v = 0, js = ns
  l0pi = ns*(1 + nzeta*(ntheta2 - 1))
  router = r1(  ns,0) + r1(  ns,1)
  rinner = r1(l0pi,0) + r1(l0pi,1)
  r00 = r1(1,0)
  z00 = z1(1,0)

  ! COMPUTE CONSTRAINT RCON, ZCON
  rcon(:nrzt,0) = rcon(:nrzt,0) + rcon(:nrzt,1)*sqrts(:nrzt)
  zcon(:nrzt,0) = zcon(:nrzt,0) + zcon(:nrzt,1)*sqrts(:nrzt)
  ru0(:nrzt)    = ru(:nrzt,0)   + ru(:nrzt,1)*sqrts(:nrzt)
  zu0(:nrzt)    = zu(:nrzt,0)   + zu(:nrzt,1)*sqrts(:nrzt)

  ! COMPUTE RCON0, ZCON0 FOR FIXED BOUNDARY BY SCALING EDGE VALUES
  ! SCALE BY POWER OF SQRTS, RATHER THAN USE rcon0 = rcon, etc.
  ! THIS PREVENTS A DISCONTINUITY WHEN RESTARTING FIXED BOUNDARY WITH NEW RCON0....
  !
  ! NOTE: IN ORDER TO MAKE INITIAL CONSTRAINT FORCES SAME FOR FREE/FIXED
  ! BOUNDARY, WE SET RCON0,ZCON0 THE SAME INITIALLY, BUT TURN THEM OFF
  ! SLOWLY IN FREE-BOUNDARY VACUUM LOOP (BELOW)
  IF (iter2.eq.iter1 .and. ivac.le.0) THEN
     DO l = 1, ns
        rcon0(l:nrzt:ns) = rcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
        zcon0(l:nrzt:ns) = zcon(ns:nrzt:ns,0)*sqrts(l:nrzt:ns)**2
     END DO
  ENDIF

  ! COMPUTE S AND THETA DERIVATIVE OF R AND Z AND JACOBIAN ON HALF-GRID
  CALL jacobian
  IF (irst.eq.2 .and. iequi.eq.0) then
     ! bad jacobian --> need to restart
     return
  end if

  ! COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC AND KINETIC
  ! PRESSURE, AND METRIC ELEMENTS ON HALF-GRID
  CALL bcovar (lu, lv)

  ! COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
  ! NOTE: FOR FREE BOUNDARY RUNS, THE VALUE OF RBTOR=R*BTOR
  ! AT THE PLASMA EDGE SHOULD BE ADJUSTED TO APPROXIMATELY
  ! EQUAL THE VACUUM VALUE. THIS CAN BE DONE BY CHANGING
  ! EITHER PHIEDGE OR THE INITIAL CROSS SECTION ACCORDING
  ! TO THE SCALING LAW  R*BTOR .EQ. PHIEDGE/(R1 * Z1).
  IF (lfreeb .and. iter2.gt.1 .and. iequi.eq.0) THEN

     IF ((fsqr + fsqz) .le. 1.e-3_dp) then
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

        raxis_nestor(1:nzeta) = r1(1:ns*nzeta:ns,0)
        zaxis_nestor(1:nzeta) = z1(1:ns*nzeta:ns,0)

        CALL vacuum (rmnc, rmns, zmns, zmnc, xm, xn,                         &
                     ctor, rbtor, wint, ns, ivacskip, ivac, mnmax, ier_flag, &
                     lasym, signgs)

        IF (ier_flag .ne. 0) then
           ! some error occured within NESTOR, so cancel the iterations
           return
        end if

        ! RESET FIRST TIME FOR SOFT START
        IF (ivac .eq. 1) THEN
           irst = 2
           delt0 = delt
           CALL restart_iter(delt0)
           irst = 1
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
           bsqsav(lk,3) = 1.5_dp*bzmn_o(l) - 0.5_dp*bzmn_o(l-1)
           gcon(l)      = bsqvac(lk) + presf_ns

           rbsq(lk) = gcon(l)*(r1(l,0) + r1(l,1))*ohs
           dbsq(lk) = ABS(gcon(l)-bsqsav(lk,3))
        END DO

        IF (ivac .eq. 1) THEN
           bsqsav(:nznt,1) = bzmn_o(ns:nrzt:ns)
           bsqsav(:nznt,2) = bsqvac(:nznt)
        ENDIF
     ENDIF
  ENDIF

  IF (iequi .NE. 1) THEN
     ! normal iterations, not final call from fileout (which sets iequi=1)

     ! COMPUTE CONSTRAINT FORCE
     extra1(:nrzt,0) = (rcon(:nrzt,0) - rcon0(:nrzt))*ru0(:nrzt) &
                     + (zcon(:nrzt,0) - zcon0(:nrzt))*zu0(:nrzt)
     CALL alias (gcon, extra1(:,0), gc, gc(1+mns), gc(1+2*mns), extra1(:,1))

     ! COMPUTE MHD FORCES ON INTEGER-MESH
     CALL forces

     ! SYMMETRIZE FORCES (in u-v space): NOTE - gc IS SMALL BY FACTOR 2 IF lasym=T
     IF (lasym) THEN
        CALL symforce (armn, brmn, crmn, azmn, bzmn,                   &
          czmn, blmn, clmn, rcon, zcon, r1, ru, rv, z1, zu, zv,        &
          extra3, extra4, extra1, extra2)

        ! NOT NECESSARY (EVEN THOUGH CORRECT)
        ! gc = 2*gc
     END IF

     ! FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
     CALL tomnsps (gc,               &
                   armn, brmn, crmn, &
                   azmn, bzmn, czmn, &
                         blmn, clmn, &
                   rcon, zcon         )
     IF (lasym) then
        CALL tomnspa (gc, r1, ru, rv, &
                          z1, zu, zv, &
                      extra3, extra4, &
                      extra1, extra2)
     end if

     ! COMPUTE FORCE RESIDUALS (RAW AND PRECONDITIONED)
     gc = gc * scalxc    !!IS THIS CORRECT: SPH010214?
     CALL residue (gc, gc(1+irzloff), gc(1+2*irzloff))

     IF (iter2.eq.1 .and. (fsqr+fsqz+fsql).gt.1.E2_dp) then
         ! first iteration and gigantic force residuals --> what is going one here?
         irst = 4
     end if

  ELSE
     ! iequi == 1 --> skip above remainder of funct3d
  END IF

END SUBROUTINE funct3d
