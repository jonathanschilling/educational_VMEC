!> \file
!> \brief Compute the real-space MHD forces.

!> \brief Compute the real-space MHD forces.
!>
SUBROUTINE forces
  USE vmec_main, p5 => cp5
  USE realspace

  ! I guess these double aliases are intended to show at which point in the code
  ! the respective arrays hold what quantity...
  ! E.g., ru12 gets filled by jacobian() and is overwritten with the force component azmn_e here.
  USE vforces,   ru12 => azmn_e,   zu12 => armn_e, &
               azmn_e => azmn_e, armn_e => armn_e, &
                 lv_e => crmn_e,   lu_e => czmn_e,   lu_o => czmn_o, &
               crmn_e => crmn_e, czmn_e => czmn_e, czmn_o => czmn_o
  IMPLICIT NONE

  REAL(rprec), PARAMETER :: p25 = p5*p5
  REAL(rprec), PARAMETER :: dshalfds=p25

  INTEGER :: l, js, lk, ku
  INTEGER :: ndim
  REAL(rprec), DIMENSION(:), POINTER :: bsqr
  REAL(rprec), DIMENSION(:), POINTER :: gvvs
  REAL(rprec), DIMENSION(:), POINTER :: guvs
  REAL(rprec), DIMENSION(:), POINTER :: guus

  character(len=255) :: dump_filename
  logical            :: dump_forces = .false.

  ! ON ENTRY, ARMN=ZU, BRMN=ZS, AZMN=RU, BZMN=RS, LU=R*BSQ, LV = BSQ*SQRT(G)/R12
  ! HERE, XS (X=Z,R) DO NOT INCLUDE DERIVATIVE OF EXPLICIT SQRT(S)
  ! BSQ = |B|**2/2 + p
  ! GIJ = (BsupI * BsupJ) * SQRT(G)  (I,J = U,V)
  ! IT IS ESSENTIAL THAT LU,LV AT j=1 ARE ZERO INITIALLY
  !
  ! SOME OF THE BIGGER LOOPS WERE SPLIT TO FACILITATE CACHE HITS, PIPELINING ON RISCS
  !
  ! ORIGIN OF VARIOUS TERMS
  !
  ! LU :  VARIATION OF DOMINANT RADIAL-DERIVATIVE TERMS IN JACOBIAN
  !
  ! LV :  VARIATION OF R-TERM IN JACOBIAN
  !
  ! GVV:  VARIATION OF R**2-TERM AND Rv**2,Zv**2 IN gvv
  !
  ! GUU, GUV: VARIATION OF Ru, Rv, Zu, Zv IN guu, guv

  ! inputs:
  ! lu_e, lv_e
  ! guu, guv, gvv
  ! ru12, zu12
  ! brmn_e, bzmn_e
  ! r1, z1, ru, zu, ru0, zu0
  ! rcon, zcon, rcon0, zcon0, gcon
  ! if lthreed:
  !  rv, zv

  if (dump_forces .and. iter2.le.2) then
    write(dump_filename, 999) ns, iter2, trim(input_extension)
999 format('forces_',i5.5,'_',i6.6,'.',a)
    open(unit=42, file=trim(dump_filename), status="unknown")

    write(42, *) "# ns ntheta3 nzeta"
    write(42, *) ns, ntheta3, nzeta

    if (lthreed) then
      write(42, *) "# js lv ku" // &
        " lu_e lv_e guu guv gvv ru12 zu12 brmn_e bzmn_e" // &
        " r1 z1 ru zu ru0 zu0 rcon zcon rcon0 zcon0 gcon rv zv"
      DO js = 1, ns
        DO lk = 1, nzeta
          DO ku = 1, ntheta3
              l = ((ku-1)*nzeta+(lk-1))*ns+js
              write (42, *) js, ku, lk, &
                lu_e(l), lv_e(l), &
                guu(l), guv(l), gvv(l), &
                ru12(l), zu12(l), &
                brmn_e(l), bzmn_e(l), &
                r1(l,:), z1(l,:), ru(l,:), zu(l,:), ru0(l), zu0(l), &
                rcon(l,0), zcon(l,0), rcon0(l), zcon0(l), gcon(l), &
                rv(l,:), zv(l,:)
          end do
        end do
      end do
    else ! lthreed
      write(42, *) "# js lv ku" // &
        " lu_e lv_e guu guv gvv ru12 zu12 brmn_e bzmn_e" // &
        " r1 z1 ru zu ru0 zu0 rcon zcon rcon0 zcon0 gcon"
      DO js = 1, ns
        DO lk = 1, nzeta
          DO ku = 1, ntheta3
              l = ((ku-1)*nzeta+(lk-1))*ns+js
              write (42, *) js, lk, ku, &
                lu_e(l), lv_e(l), &
                guu(l), guv(l), gvv(l), &
                ru12(l), zu12(l), &
                brmn_e(l), bzmn_e(l), &
                r1(l,:), z1(l,:), ru(l,:), zu(l,:), ru0(l), zu0(l), &
                rcon(l,0), zcon(l,0), rcon0(l), zcon0(l), gcon(l)
          end do
        end do
      end do
    end if ! lthreed
  end if ! dump_forces

  ndim = 1+nrzt ! TODO: remove this; one extra element at the end of a large vector sound like reconstruction stuff...

  ! POINTER ALIASES
  bsqr => extra1(:,1) ! output or temp
  gvvs => extra2(:,1) ! output or temp
  guvs => extra3(:,1) ! output or temp
  guus => extra4(:,1) ! output or temp

  ! zero values at axis
  lu_e(1:ndim:ns) = 0 ! fixup input
  lv_e(1:ndim:ns) = 0 ! fixup input
  guu(1:ndim:ns) = 0 ! fixup input
  guv(1:ndim:ns) = 0 ! fixup input
  gvv(1:ndim:ns) = 0 ! fixup input

  guus = guu*shalf ! output or temp
  guvs = guv*shalf ! output or temp
  gvvs = gvv*shalf ! output or temp

  armn_e  = ohs*zu12 * lu_e ! output or temp
  azmn_e  =-ohs*ru12 * lu_e ! output or temp
  brmn_e  = brmn_e * lu_e ! output or temp
  bzmn_e  =-bzmn_e * lu_e ! output or temp
  bsqr    = dshalfds*lu_e/shalf ! output or temp

  armn_o(1:ndim)  = armn_e(1:ndim) *shalf ! output or temp
  azmn_o(1:ndim)  = azmn_e(1:ndim) *shalf ! output or temp
  brmn_o(1:ndim)  = brmn_e(1:ndim) *shalf ! output or temp
  bzmn_o(1:ndim)  = bzmn_e(1:ndim) *shalf ! output or temp

  ! CONSTRUCT CYLINDRICAL FORCE KERNELS
  ! NOTE: presg(ns+1) == 0, AND WILL BE "FILLED IN" AT EDGE FOR FREE-BOUNDARY BY RBSQ
  DO l = 1, nrzt
     guu(l) = p5*(guu(l) + guu(l+1))
     gvv(l) = p5*(gvv(l) + gvv(l+1))
     bsqr(l) = bsqr(l) + bsqr(l+1)
     guus(l) = p5*(guus(l) + guus(l+1))
     gvvs(l) = p5*(gvvs(l) + gvvs(l+1))
     armn_e(l) = armn_e(l+1) - armn_e(l) + p5*(lv_e(l) + lv_e(l+1))
     azmn_e(l) = azmn_e(l+1) - azmn_e(l)
     brmn_e(l) = p5*(brmn_e(l) + brmn_e(l+1))
     bzmn_e(l) = p5*(bzmn_e(l) + bzmn_e(l+1))
  END DO

  armn_e(:nrzt) = armn_e(:nrzt)                            - (gvvs(:nrzt)*r1(:nrzt,1) + gvv(:nrzt)*r1(:nrzt,0))
  brmn_e(:nrzt) = brmn_e(:nrzt) +  bsqr(:nrzt)*z1(:nrzt,1) - (guus(:nrzt)*ru(:nrzt,1) + guu(:nrzt)*ru(:nrzt,0))
  bzmn_e(:nrzt) = bzmn_e(:nrzt) - (bsqr(:nrzt)*r1(:nrzt,1) +  guus(:nrzt)*zu(:nrzt,1) + guu(:nrzt)*zu(:nrzt,0))
  lv_e(1:ndim) = lv_e(1:ndim)*shalf(1:ndim)
  lu_o(1:ndim) = dshalfds*lu_e(1:ndim)

  DO l = 1, nrzt
     armn_o(l) = armn_o(l+1) - armn_o(l) - zu(l,0)*bsqr(l) + p5*(lv_e(l) + lv_e(l+1))
     azmn_o(l) = azmn_o(l+1) - azmn_o(l) + ru(l,0)*bsqr(l)
     brmn_o(l) = p5*(brmn_o(l) + brmn_o(l+1))
     bzmn_o(l) = p5*(bzmn_o(l) + bzmn_o(l+1))
     lu_o(l)   = lu_o(l) + lu_o(l+1)
  END DO

  guu(1:nrzt)  = guu(1:nrzt) * sqrts(1:nrzt)**2
  bsqr(1:nrzt) = gvv(1:nrzt) * sqrts(1:nrzt)**2

  armn_o(:nrzt) = armn_o(:nrzt) - (zu(:nrzt,1)*lu_o(:nrzt) + bsqr(:nrzt)*r1(:nrzt,1) + gvvs(:nrzt)*r1(:nrzt,0))
  azmn_o(:nrzt) = azmn_o(:nrzt) +  ru(:nrzt,1)*lu_o(:nrzt)
  brmn_o(:nrzt) = brmn_o(:nrzt) +  z1(:nrzt,1)*lu_o(:nrzt) - (guu(:nrzt)*ru(:nrzt,1) + guus(:nrzt)*ru(:nrzt,0))
  bzmn_o(:nrzt) = bzmn_o(:nrzt) - (r1(:nrzt,1)*lu_o(:nrzt) + guu(:nrzt)*zu(:nrzt,1) + guus(:nrzt)*zu(:nrzt,0))

  IF (lthreed) THEN
     DO l = 1, nrzt
        guv(l)  = p5*(guv(l) + guv(l+1))
        guvs(l) = p5*(guvs(l) + guvs(l+1))
     END DO

     brmn_e(:nrzt) = brmn_e(:nrzt) - (guv(:nrzt)*rv(:nrzt,0) + guvs(:nrzt)*rv(:nrzt,1))
     bzmn_e(:nrzt) = bzmn_e(:nrzt) - (guv(:nrzt)*zv(:nrzt,0) + guvs(:nrzt)*zv(:nrzt,1))
     crmn_e(:nrzt) = guv(:nrzt) *ru(:nrzt,0) + gvv(:nrzt) *rv(:nrzt,0) + gvvs(:nrzt)*rv(:nrzt,1) + guvs(:nrzt)*ru(:nrzt,1)
     czmn_e(:nrzt) = guv(:nrzt) *zu(:nrzt,0) + gvv(:nrzt) *zv(:nrzt,0) + gvvs(:nrzt)*zv(:nrzt,1) + guvs(:nrzt)*zu(:nrzt,1)

     guv(:nrzt) = guv(:nrzt) * sqrts(:nrzt)*sqrts(:nrzt)

     brmn_o(:nrzt) = brmn_o(:nrzt) - (guvs(:nrzt)*rv(:nrzt,0) + guv(:nrzt)*rv(:nrzt,1))
     bzmn_o(:nrzt) = bzmn_o(:nrzt) - (guvs(:nrzt)*zv(:nrzt,0) + guv(:nrzt)*zv(:nrzt,1))
     crmn_o(:nrzt) = guvs(:nrzt)*ru(:nrzt,0) + gvvs(:nrzt)*rv(:nrzt,0) + bsqr(:nrzt)*rv(:nrzt,1) + guv(:nrzt) *ru(:nrzt,1)
     czmn_o(:nrzt) = guvs(:nrzt)*zu(:nrzt,0) + gvvs(:nrzt)*zv(:nrzt,0) + bsqr(:nrzt)*zv(:nrzt,1) + guv(:nrzt) *zu(:nrzt,1)
  ENDIF

  ! ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY CALCULATION
  IF (ivac .ge. 1) THEN
     armn_e(ns:nrzt:ns) = armn_e(ns:nrzt:ns) + zu0(ns:nrzt:ns)*rbsq(1:nznt)
     armn_o(ns:nrzt:ns) = armn_o(ns:nrzt:ns) + zu0(ns:nrzt:ns)*rbsq(1:nznt)
     azmn_e(ns:nrzt:ns) = azmn_e(ns:nrzt:ns) - ru0(ns:nrzt:ns)*rbsq(1:nznt)
     azmn_o(ns:nrzt:ns) = azmn_o(ns:nrzt:ns) - ru0(ns:nrzt:ns)*rbsq(1:nznt)
  ENDIF

! #ifndef _HBANGLE
  ! COMPUTE CONSTRAINT FORCE KERNELS
  rcon(:nrzt,0) = (rcon(:nrzt,0) - rcon0(:nrzt)) * gcon(:nrzt)
  zcon(:nrzt,0) = (zcon(:nrzt,0) - zcon0(:nrzt)) * gcon(:nrzt)

  brmn_e(:nrzt) = brmn_e(:nrzt) + rcon(:nrzt,0)
  bzmn_e(:nrzt) = bzmn_e(:nrzt) + zcon(:nrzt,0)
  brmn_o(:nrzt) = brmn_o(:nrzt) + rcon(:nrzt,0) * sqrts(:nrzt)
  bzmn_o(:nrzt) = bzmn_o(:nrzt) + zcon(:nrzt,0) * sqrts(:nrzt)

  ! real-space B-type forces due to constraint only
  brmn_e_con(:nrzt) = brmn_e_con(:nrzt) + rcon(:nrzt,0)
  bzmn_e_con(:nrzt) = bzmn_e_con(:nrzt) + zcon(:nrzt,0)
  brmn_o_con(:nrzt) = brmn_o_con(:nrzt) + rcon(:nrzt,0) * sqrts(:nrzt)
  bzmn_o_con(:nrzt) = bzmn_o_con(:nrzt) + zcon(:nrzt,0) * sqrts(:nrzt)

  rcon(:nrzt,0) =  ru0(:nrzt) * gcon(:nrzt)
  zcon(:nrzt,0) =  zu0(:nrzt) * gcon(:nrzt)
  rcon(:nrzt,1) = rcon(:nrzt,0) * sqrts(:nrzt)
  zcon(:nrzt,1) = zcon(:nrzt,0) * sqrts(:nrzt)
! #end /* ndef _HBANGLE */

  if (dump_forces .and. iter2.le.2) then

    if (lthreed) then
      write(42, *) "# js lv ku" // &
        " armn_e armn_o brmn_e brmn_o crmn_e crmn_o" // &
        " azmn_e azmn_o bzmn_e bzmn_o czmn_e czmn_o" // &
        " guu guus guv guvs gvv gvvs" // &
        " bsqr lu_o lv_e rcon zcon"
      l = 1
      DO js = 1, ns
        DO lk = 1, nzeta
          DO ku = 1, ntheta3
              l = ((ku-1)*nzeta+(lk-1))*ns+js
              write (42, *) js, lk, ku, &
                armn_e(l), armn_o(l), brmn_e(l), brmn_o(l), crmn_e(l), crmn_o(l), &
                azmn_e(l), azmn_o(l), bzmn_e(l), bzmn_o(l), czmn_e(l), czmn_o(l), &
                guu(l), guus(l), guv(l), guvs(l), gvv(l), gvvs(l), &
                bsqr(l), lu_o(l), lv_e(l), rcon(l,:), zcon(l,:)
          end do
        end do
      end do
    else ! lthreed
      write(42, *) "# js lv ku" // &
        " armn_e armn_o brmn_e brmn_o" // &
        " azmn_e azmn_o bzmn_e bzmn_o" // &
        " guu guus guv guvs gvv gvvs" // &
        " bsqr lu_o lv_e rcon zcon"
      l = 1
      DO js = 1, ns
        DO lk = 1, nzeta
          DO ku = 1, ntheta3
              l = ((ku-1)*nzeta+(lk-1))*ns+js
              write (42, *) js, lk, ku, &
                armn_e(l), armn_o(l), brmn_e(l), brmn_o(l), &
                azmn_e(l), azmn_o(l), bzmn_e(l), bzmn_o(l), &
                guu(l), guus(l), guv(l), guvs(l), gvv(l), gvvs(l), &
                bsqr(l), lu_o(l), lv_e(l), rcon(l,:), zcon(l,:)
          end do
        end do
      end do
    end if ! lthreed

    close(42)

    print *, "dumped forces to '"//trim(dump_filename)//"'"
!     stop
  end if

END SUBROUTINE forces
