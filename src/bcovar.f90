!> \file
!> \brief Compute the covariant components of the magnetic field \f$B_\theta\f$, \f$B_\zeta\f$.

!> \brief Compute the covariant components of the magnetic field \f$B_\theta\f$, \f$B_\zeta\f$.
!>
!> @param lu \f$\partial\lambda / \partial\theta\f$
!> @param lv \f$- \partial\lambda / \partial\zeta\f$
SUBROUTINE bcovar (lu, lv)
  USE vmec_main, fpsi => bvco, p5 => cp5
  USE vmec_params, ONLY: ns4, signgs, pdamp, lamscale, ntmax, &
                         meven, modd
  USE realspace
  USE vforces, r12 => armn_o, ru12 => azmn_e, gsqrt => azmn_o,      &
               rs => bzmn_e, zs => brmn_e, zu12 => armn_e,          &
               bsubu_e => clmn_e, bsubv_e => blmn_e,                &
               bsubu_o => clmn_o, bsubv_o => blmn_o,                &
               bsq => bzmn_o, phipog => brmn_o
  USE xstuff, ONLY: xc

  use dbgout

  IMPLICIT NONE

  REAL(rprec), DIMENSION(nrzt,0:1), INTENT(inout) :: lu
  REAL(rprec), DIMENSION(nrzt,0:1), INTENT(inout) :: lv

  ! GENERALLY, IF TEMPORAL CONVERGENCE IS POOR, TRY TO INCREASE PDAMP (< 1) (STORED IN VMEC_PARAMS)

  INTEGER :: l, js, ndim, lk, ku, m, n, rzl
  REAL(rprec) :: r2, volume, curpol_temp

! #ifndef _HBANGLE
  REAL(rprec) :: arnorm, aznorm, tcon_mul
! #end /* ndef _HBANGLE */

  REAL(rprec), POINTER, DIMENSION(:) :: luu, luv, lvv, tau
  REAL(rprec), DIMENSION(:), POINTER :: bsupu, bsubuh, bsupv, bsubvh, r12sq

  ndim = 1+nrzt ! what is hidden at the end of these vectors? probably leftover from reconstruction stuff...

  ! POINTER ALIAS ASSIGNMENTS
  tau => extra1(:,1)
  luu => extra2(:,1)
  luv => extra3(:,1)
  lvv => extra4(:,1)

  bsupu  => luu
  bsupv  => luv
  bsubuh => bsubu_o
  bsubvh => bsubv_o
  r12sq  => bsq

  guu(ndim) = 0 ! TODO: can probably go away if last element is not used anyway...
  guv = 0
  gvv = 0

  ! COMPUTE METRIC ELEMENTS GIJ ON HALF MESH
  ! FIRST, GIJ = EVEN PART (ON FULL MESH), LIJ = ODD PART (ON FULL MESH)
  ! THEN, GIJ(HALF) = < GIJ(even)> + SHALF < GIJ(odd) >
  r12sq(1:nrzt) = sqrts(1:nrzt)*sqrts(1:nrzt) ! s on full grid
  guu(1:nrzt)   =   ru(1:nrzt,meven)*ru(1:nrzt,meven)                         &
                  + zu(1:nrzt,meven)*zu(1:nrzt,meven) + r12sq(1:nrzt)*         &
                 (  ru(1:nrzt,modd )*ru(1:nrzt,modd )                          &
                  + zu(1:nrzt,modd )*zu(1:nrzt,modd ))

  luu(1:nrzt)   = (  ru(1:nrzt,meven)*ru(1:nrzt,modd)                        &
                   + zu(1:nrzt,meven)*zu(1:nrzt,modd))*2

  ! temporary re-use of phipog for contribution to g_vv^o from R^2
  phipog(1:nrzt)= 2* r1(1:nrzt,meven)*r1(1:nrzt,modd)
  !do l=1,2*ns
  !  print *, l, r1(l,meven), r1(l,modd), phipog(l)
  !end do

  IF (lthreed) THEN
     guv(1:nrzt)   = ru(1:nrzt,meven)*rv(1:nrzt,meven)                      &
                   + zu(1:nrzt,meven)*zv(1:nrzt,meven) + r12sq(1:nrzt)*     &
                   ( ru(1:nrzt,modd )*rv(1:nrzt,modd )                      &
                   + zu(1:nrzt,modd )*zv(1:nrzt,modd ) )
     luv(1:nrzt)   = ru(1:nrzt,meven)*rv(1:nrzt,modd )                      &
                   + ru(1:nrzt,modd )*rv(1:nrzt,meven)                      &
                   + zu(1:nrzt,meven)*zv(1:nrzt,modd )                      &
                   + zu(1:nrzt,modd )*zv(1:nrzt,meven)

     gvv(1:nrzt)   = rv(1:nrzt,meven)*rv(1:nrzt,meven)                      &
                   + zv(1:nrzt,meven)*zv(1:nrzt,meven) + r12sq(1:nrzt)*     &
                   ( rv(1:nrzt,modd )*rv(1:nrzt,modd )                      &
                   + zv(1:nrzt,modd )*zv(1:nrzt,modd ) )
     lvv(1:nrzt)   =(rv(1:nrzt,meven)*rv(1:nrzt,modd )                      &
                   + zv(1:nrzt,meven)*zv(1:nrzt,modd ))*2
  END IF

  ! contribution to g_vv^e from R^2
  r12sq(1:nrzt) = r1(1:nrzt,meven)*r1(1:nrzt,meven) + r12sq(1:nrzt)*        &
                  r1(1:nrzt,modd )*r1(1:nrzt,modd )

  ! need to do this in reverse order since guu and r12sq were re-used
  ! for their full-grid even-m contributions
  DO l = nrzt, 2, -1
     guu(l) = p5*(guu(l) + guu(l-1) + shalf(l)*(luu(l) + luu(l-1)))
     ! g_uu done here

!     if (mod(l-1,ns).eq.0) then
!       print *, (l-1)/ns+1, l, "evn_contrib=",p5*(  r12sq(l) +  r12sq(l-1) )," from ", r12sq(l)," and ", r12sq(l-1)
!       print *, (l-1)/ns+1, l, "odd_contrib=",p5*( phipog(l) + phipog(l-1) )," from ",phipog(l)," and ",phipog(l-1)
!     end if
     ! r12sq = R^2 for g_vv here
     r12sq(l) = p5*( r12sq(l) + r12sq(l-1) + shalf(l)*(phipog(l) + phipog(l-1)) )
  END DO

  IF (lthreed) THEN
     DO l = nrzt, 2, -1
        guv(l) = p5*(guv(l) + guv(l-1) + shalf(l)*(luv(l) + luv(l-1)))
        gvv(l) = p5*(gvv(l) + gvv(l-1) + shalf(l)*(lvv(l) + lvv(l-1)))
     END DO
  END IF

  ! save tau as output from jacobian()
  ! TODO: for what later use is tau needed separately ?
  tau(1:nrzt) = gsqrt(1:nrzt)

  ! actually build \sqrt{g} now with tau from jacobian()
  ! and R^2 on the half-mesh from jacobian() as well
  gsqrt(1:nrzt) = r12(1:nrzt)*tau(1:nrzt) ! r12 = R on half grid

  ! constant extrapolation of sqrt(g) towards axis
  ! TODO: or does this overwrite the weird terms at the ghost surface outside the LCFS ???
  gsqrt(1:nrzt:ns) = gsqrt(2:nrzt:ns)

  ! g_vv contains all up to R^2 here; so add R^2 term now
  gvv(2:nrzt) = gvv(2:nrzt) + r12sq(2:nrzt)

  if (open_dbg_context("metric")) then

      call add_real_3d("gsqrt", ns, nzeta, ntheta3, gsqrt)
      call add_real_3d("guu",   ns, nzeta, ntheta3, guu  )
      call add_real_3d("r12sq", ns, nzeta, ntheta3, r12sq)
      call add_real_3d("gvv",   ns, nzeta, ntheta3, gvv  )

      if (lthreed) then
        call add_real_3d("guv", ns, nzeta, ntheta3, guv  )
      else
        call add_null("guv")
      end if

      call close_dbg_out()
  end if

  ! CATCH THIS AFTER WHERE LINE BELOW   phipog = 0

  ! this is where phipog == 1/sqrt(g) is actually assigned
  ! note that the phip factor in phipog is gone... (see comment below)
  WHERE (gsqrt(2:ndim) .ne. zero) phipog(2:ndim) = one/gsqrt(2:ndim)

  ! 1/sqrt(g) is zero(since undefined) at the magnetic axis
  phipog(1:ndim:ns) = 0

  ! compute plasma volume profile (vp) and total volume (voli)
  vp(1) = 0
  vp(ns+1) = 0
  DO js = 2, ns
     vp(js) = signgs*SUM(gsqrt(js:nrzt:ns)*wint(js:nrzt:ns))
  END DO
  IF (iter2 .eq. 1) then
    voli = twopi*twopi*hs*SUM(vp(2:ns))
  end if

  ! check plasma volume computation
  if (open_dbg_context("volume")) then

    call add_real_1d("vp", ns+1, vp)
    call add_real("voli", voli)

    call close_dbg_out()
  end if

  ! COMPUTE CONTRA-VARIANT COMPONENTS OF B (Bsupu,v) ON RADIAL HALF-MESH TO ACCOMODATE LRFP=T CASES.
  ! THE OVERALL PHIP FACTOR (PRIOR TO v8.46) HAS BEEN REMOVED FROM PHIPOG, SO NOW PHIPOG == 1/GSQRT!
  !
  ! NOTE: LU = LAMU == d(LAM)/du, LV = -LAMV == -d(LAM)/dv COMING INTO THIS ROUTINE.
  ! CHIP will be added IN CALL TO ADD_FLUXES.
  ! THE NET BSUPU, BSUPV ARE (PHIPOG=1/GSQRT AS NOTED ABOVE):
  !
  !      BSUPU = PHIPOG*(chip - LAMV*LAMSCALE),
  !      BSUPV = PHIPOG*(phip + LAMU*LAMSCALE)
  lu = lu*lamscale
  lv = lv*lamscale

  DO js=1,ns
     lu(js:nrzt:ns,0) = lu(js:nrzt:ns,0) + phipf(js)
  END DO

  bsupu(2:nrzt) = phipog(2:nrzt) * p5*(                 lv(2:nrzt,0) + lv(1:nrzt-1,0)  &
                                       + shalf(2:nrzt)*(lv(2:nrzt,1) + lv(1:nrzt-1,1))  )
  bsupv(2:nrzt) = phipog(2:nrzt) * p5*(                 lu(2:nrzt,0) + lu(1:nrzt-1,0)  &
                                       + shalf(2:nrzt)*(lu(2:nrzt,1) + lu(1:nrzt-1,1))  )

  ! first point at (u,v)=(0,0) on axis
  bsupu(1)=0
  bsupv(1)=0

  ! v8.49: add ndim points --> TODO: likely not needed anymore, reconstruction-related?
  bsupu(ndim)=0
  bsupv(ndim)=0

  if (open_dbg_context("bcontrav")) then

    call add_real_3d("bsupu", ns, nzeta, ntheta3, bsupu)
    call add_real_3d("bsupv", ns, nzeta, ntheta3, bsupv)

    call close_dbg_out()
  end if

  ! UPDATE IOTA EITHER OF TWO WAYS:
  ! 1)  FOR ictrl_prec2d = 0, SOLVE THE LINEAR ALGEBRAIC EQUATION <Bsubu> = icurv FOR iotas
  ! 2)  FOR ictrl_prec2d > 0, EVOLVE IOTAS IN TIME, USING Force-iota  = <Bsubu> - icurv.
  !
  ! NEED TO DO IT WAY (#2) TO EASILY COMPUTE THE HESSIAN ELEMENTS DUE TO LAMBDA-VARIATIONS.
  ! IOTAS IS "STORED" AT LOCATION LAMBDA-SC(0,0) IN XC-ARRAY [USE THIS COMPONENT SO IT
  ! WILL WORK EVEN FOR 2D PLASMA], ALTHOUGH ITS VARIATION IS LIKE THAT OF LV-CS(0,0),
  ! WITH N -> 1 IN THE HESSIAN CALCULATION ROUTINES (Compute_Hessian_Flam_lam, etc.)

  ! COMPUTE (IF NEEDED) AND ADD CHIP TO BSUPU
  CALL add_fluxes(phipog, bsupu, bsupv)

  ! COMPUTE COVARIANT B COMPONENT bsubu,v (LAMBDA FORCE KERNELS) ON RADIAL HALF-MESH
  bsubuh(1:nrzt) = guu(1:nrzt)*bsupu(1:nrzt) + guv(1:nrzt)*bsupv(1:nrzt)
  bsubvh(1:nrzt) = guv(1:nrzt)*bsupu(1:nrzt) + gvv(1:nrzt)*bsupv(1:nrzt)

  ! v8.49 --> TODO: likely not needed anymore, reconstruction-related?
  bsubuh(ndim) = 0
  bsubvh(ndim) = 0

  ! COMPUTE MAGNETIC AND KINETIC PRESSURE ON RADIAL HALF-MESH
  ! bsq = |B|^2/2 = 0.5*(B^u*B_u + B^v*B_v)
  bsq(:nrzt) = p5*(bsupu(:nrzt)*bsubuh(:nrzt) + bsupv(:nrzt)*bsubvh(:nrzt))
  pres(2:ns) = mass(2:ns)/vp(2:ns)**gamma

  ! magnetic energy
  wb = hs*ABS(SUM(wint(:nrzt)*gsqrt(:nrzt)*bsq(:nrzt)))

  ! kinetic == thermal enery
  wp = hs*SUM(vp(2:ns)*pres(2:ns))

!   write(*,*) "magnetic energy: ", wb
!   write(*,*) "kinetic energy: ", wp

  ! ADD KINETIC PRESSURE TO MAGNETIC PRESSURE
  DO js=2,ns
     bsq(js:nrzt:ns) = bsq(js:nrzt:ns) + pres(js)
  END DO

  if (open_dbg_context("bcov")) then

    call add_real_3d("bsubuh", ns, nzeta, ntheta3, bsubuh)
    call add_real_3d("bsubvh", ns, nzeta, ntheta3, bsubvh)
    call add_real_3d("bsq",    ns, nzeta, ntheta3, bsq   )

    call add_real_1d("pres", ns-1, pres(2:ns))

    call add_real("wb", wb)
    call add_real("wp", wp)

    call close_dbg_out()
  end if

  ! COMPUTE LAMBDA FULL MESH FORCES
  ! NOTE: bsubu_e is used here ONLY as a temporary array (TODO: for what?)
  lvv = phipog(:ndim)*gvv
  bsubv_e(1:nrzt) = p5*(lvv(1:nrzt)+lvv(2:ndim))*lu(1:nrzt,0)

!   if (iter2.eq.2) then
!     write(*,*) "bsubv_e(1)  first step: ", bsubv_e(1)
!   end if

  lvv = lvv*shalf
  bsubu_e(:nrzt) = guv(:nrzt)*bsupu(:nrzt)
  bsubu_e(ndim) = 0
  bsubv_e(1:nrzt) = bsubv_e(1:nrzt)                                 &
              + p5*((lvv(1:nrzt) + lvv(2:ndim))*lu(1:nrzt,1)        &
              +      bsubu_e(1:nrzt) + bsubu_e(2:ndim))

!   if (iter2.eq.2) then
!     write(*,*) "bsubu_e(1:2) = ", bsubu_e(1:2), " from ", &
!       guv(1:2), bsupu(1:2)
!     write(*,*) "bsubv_e(1) second step: ", bsubv_e(1), lvv(1), lvv(2), &
!       lu(1,1), bsubu_e(1), bsubu_e(2)
!   end if

   if (open_dbg_context("lambda_forces")) then

    call add_real_3d("lvv",     ns,    nzeta, ntheta3, lvv    )
    call add_real_4d("lu",      ns, 2, nzeta, ntheta3, lu, order=(/1, 3, 4, 2 /))
    call add_real_3d("bsubu_e", ns,    nzeta, ntheta3, bsubu_e)
    call add_real_3d("bsubv_e", ns,    nzeta, ntheta3, bsubv_e)

    call close_dbg_out()
  end if

  ! COMPUTE AVERAGE FORCE BALANCE AND TOROIDAL/POLOIDAL CURRENTS
  CALL calc_fbal(bsubuh, bsubvh)

  ! fpsi is simply an alias to bvco (which is filled in calc_fbal)
  ! --> why not use bvco here ???
  ! This computes the poloidal current close to the axis (rbtor0)
  ! and the poloidal current at the boundary (rbtor).
  rbtor0= c1p5*fpsi(2)  - p5*fpsi(3)
  rbtor = c1p5*fpsi(ns) - p5*fpsi(ns-1)

  ! (SPH:08/19/04)
  ! MUST AVOID BREAKING TRI-DIAGONAL RADIAL COUPLING AT EDGE WHEN USING PRECONDITIONER
  ! CTOR IS PASSED TO VACUUM TO COMPUTE EDGE BSQVAC, SO IT CAN ONLY DEPEND ON NS, NS-1
  ! THUS, CTOR ~ buco(ns) WORKS, WITH REMAINDER A FIXED CONSTANT.
  !
  ! ALSO, IF USING FAST SWEEP IN COMPUTE_BLOCKS, MUST MAKE CTOR CONSTANT
  ! TO AVOID BREAKING SYMMETRY OF A+(ns-1) AND B-(ns) HESSIAN ELEMENTS
  !
  ! TO GET CORRECT HESSIAN, USE THE CTOR=ctor_prec2d +... ASSIGNMENT
  ! FOR ictrl_prec2d.ne.0 (replace ictrl_prec2d.gt.1 with ictrl_prec2d.ne.0 in IF test below)

  ! NEXT COMPUTE COVARIANT BSUBV COMPONENT ~ lvv ON FULL RADIAL MESH BY AVERAGING HALF-MESH METRICS
  ! NOTE: EDGE VALUES AT JS=NS DOWN BY 1/2
  ! THIS IS NEEDED FOR NUMERICAL STABILITY

  ! This computes the net toroidal current enclosed by the LCFS (ctor).
  ctor = signgs*twopi*(c1p5*buco(ns) - p5*buco(ns1))

  DO l=1,ns
     lvv(l:nrzt:ns) = bdamp(l) ! --> profil1d()
     ! blending parameter: 0.1 at the axis, ca. 0 at the LCFS
     ! TODO: is it intentional that this is linked to pdamp (time-step algorithm) ?
  END DO

  ! COMMENTED OUT BY SAL --> why does this check hurt ?
  ! IF (ANY(bsubuh(1:ndim:ns) .ne. zero)) STOP 'BSUBUH != 0 AT JS=1'
  ! IF (ANY(bsubvh(1:ndim:ns) .ne. zero)) STOP 'BSUBVH != 0 AT JS=1'

  ! AVERAGE LAMBDA FORCES ONTO FULL RADIAL MESH
  ! USE BLENDING FOR bsubv_e FOR NUMERICAL STABILITY NEAR AXIS
  bsubu_e(1:nrzt) =  p5*                (bsubuh(1:nrzt) + bsubuh(2:ndim))
  bsubv_e(1:nrzt) =        lvv(1:nrzt) *        bsubv_e(1:nrzt)           &
                   + p5*(1-lvv(1:nrzt))*(bsubvh(1:nrzt) + bsubvh(2:ndim))

  if (open_dbg_context("bcov_full")) then

    call add_real("rbtor0", rbtor0)
    call add_real("rbtor",  rbtor)
    call add_real("ctor",   ctor)

    call add_real_1d("bdamp", ns, bdamp)

    call add_real_3d("bsubu_e", ns, nzeta, ntheta3, bsubu_e)
    call add_real_3d("bsubv_e", ns, nzeta, ntheta3, bsubv_e)

    call close_dbg_out()
  end if

  if (iequi .eq. 0) then
    ! COMPUTE R,Z AND LAMBDA PRE-CONDITIONING MATRIX ELEMENTS AND FORCE NORMS:

    ! COMPUTE PRECONDITIONING (1D) AND SCALING PARAMETERS
    ! NO NEED TO RECOMPUTE WHEN 2D-PRECONDITIONER ON
    IF (MOD(iter2-iter1,ns4).eq.0) THEN
       ! only update preconditioner every ns4==25 iterations (?) (for ns4, see vmec_params)

!        write(*,*) "update 1d preconditioner"

       phipog(:nrzt) = phipog(:nrzt)*wint(:nrzt) ! remember that actually phipog == 1/sqrt(g)

       CALL lamcal(phipog, guu, guv, gvv)

       CALL precondn(bsupv, bsq, gsqrt, r12, &
                     zs, zu12, zu, zu(1,1), z1(1,1), &
                     arm, ard, brm, brd, crd, cos01)

       CALL precondn(bsupv, bsq, gsqrt, r12, &
                     rs, ru12, ru, ru(1,1), r1(1,1), &
                     azm, azd, bzm, bzd, crd, sin01)

       ! check preconditioner output
       if (open_dbg_context("precondn")) then

         call add_real_2d("arm", ns+1, 2, arm)
         call add_real_2d("ard", ns+1, 2, ard)
         call add_real_2d("brm", ns+1, 2, brm)
         call add_real_2d("brd", ns+1, 2, brd)
         call add_real_1d("crd", ns+1,    crd)
         call add_real_2d("azm", ns+1, 2, azm)
         call add_real_2d("azd", ns+1, 2, azd)
         call add_real_2d("bzm", ns+1, 2, bzm)
         call add_real_2d("bzd", ns+1, 2, bzd)

         call close_dbg_out()
       end if

       volume = hs*SUM(vp(2:ns))

       r2 = MAX(wb,wp)/volume ! energy density ???

       !> R12 from RP in force
       guu(:nrzt) = guu(:nrzt)*r12(:nrzt)**2

       !> Norm, unpreconditioned R,Z forces
       fnorm = one/(SUM(guu(1:nrzt)*wint(1:nrzt))*(r2*r2))

       !> Norm for preconditioned R,Z forces
       fnorm1 = one/SUM(xc(1+ns:2*irzloff)**2)

       !> Norm for unpreconditioned Lambda force
       fnormL = one/(SUM((bsubuh(1:nrzt)**2 + bsubvh(1:nrzt)**2)*wint(1:nrzt))*lamscale**2)

       ! r3 = one/(2*r0scale)**2
       ! > Norm for preconditioned Lambda force
       ! fnorm2 = one/MAX(SUM(xc(2*irzloff+1:3*irzloff)**2),r3/4)

       ! COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
       ! OVERRIDE USER INPUT VALUE HERE

! #ifndef _HBANGLE
       r2 = ns ! temporary re-use of variable for floating-point version of ns

       ! ignore large tcon0 from old-style files
       tcon0 = MIN(ABS(tcon0), one)

       ! some parabola in ns, but why these specific values of the parameters ?
       tcon_mul = tcon0*(1 + r2*(one/60 + r2/(200*120)))

       tcon_mul = tcon_mul/((4*r0scale**2)**2)           ! Scaling of ard, azd (2*r0scale**2);
                                                         ! Scaling of cos**2 in alias (4*r0scale**2)

!       print *, "tcon_mul = ", tcon_mul

       DO js = 2, ns-1
         arnorm = SUM(wint(js:nrzt:ns)*ru0(js:nrzt:ns)**2)
         aznorm = SUM(wint(js:nrzt:ns)*zu0(js:nrzt:ns)**2)
         IF (arnorm.eq.zero .or. aznorm.eq.zero) then
            STOP 'arnorm or aznorm=0 in bcovar'
         end if

         tcon(js) = MIN(ABS(ard(js,1)/arnorm), ABS(azd(js,1)/aznorm)) * tcon_mul*(32*hs)**2
       END DO

!       print *, "tcon profile: ", tcon

       tcon(ns) = p5*tcon(ns-1)

       IF (lasym) tcon = p5*tcon
! #end /* ndef _HBANGLE */

       if (open_dbg_context("forceNorms_tcon")) then

         call add_real("volume", volume)
         call add_real("r2",     MAX(wb,wp)/volume)
         call add_real("fnorm",  fnorm)
         call add_real("fnorm1",  fnorm1)
         call add_real("fnormL",  fnormL)
         call add_real("tcon0",  tcon0)
         call add_real("tcon_mul",  tcon_mul)
         call add_real_1d("tcon", ns-1, tcon(2:ns))

         call add_real_3d("guu", ns, nzeta, ntheta3,        guu)
         call add_real_5d("xc",  3, ntmax, ns, ntor1, mpol, xc, order=(/ 3, 4, 5, 2, 1 /) )

         call close_dbg_out()
       end if

     ENDIF ! MOD(iter2-iter1,ns4).eq.0

     ! MINUS SIGN => HESSIAN DIAGONALS ARE POSITIVE
     bsubu_e = -lamscale*bsubu_e
     bsubv_e = -lamscale*bsubv_e
     bsubu_o(:nrzt)  = sqrts(:nrzt)*bsubu_e(:nrzt)
     bsubv_o(:nrzt)  = sqrts(:nrzt)*bsubv_e(:nrzt)

     ! STORE LU * LV COMBINATIONS USED IN FORCES
     lvv(2:nrzt)  = gsqrt(2:nrzt)
     guu(2:nrzt)  = bsupu(2:nrzt)*bsupu(2:nrzt)*lvv(2:nrzt)
     guv(2:nrzt)  = bsupu(2:nrzt)*bsupv(2:nrzt)*lvv(2:nrzt)
     gvv(2:nrzt)  = bsupv(2:nrzt)*bsupv(2:nrzt)*lvv(2:nrzt)
     lv(2:nrzt,0) = bsq(2:nrzt)*tau(2:nrzt)
     lu(2:nrzt,0) = bsq(2:nrzt)*r12(2:nrzt)

     if (open_dbg_context("lulv_comb")) then

       call add_real_3d("bsubu_e", ns,    nzeta, ntheta3, bsubu_e )
       call add_real_3d("bsubv_e", ns,    nzeta, ntheta3, bsubv_e )
       call add_real_3d("bsubu_o", ns,    nzeta, ntheta3, bsubu_o )
       call add_real_3d("bsubv_o", ns,    nzeta, ntheta3, bsubv_o )
       call add_real_3d("lvv",     ns,    nzeta, ntheta3, lvv     )
       call add_real_3d("guu",     ns,    nzeta, ntheta3, guu     )
       call add_real_3d("guv",     ns,    nzeta, ntheta3, guv     )
       call add_real_3d("gvv",     ns,    nzeta, ntheta3, gvv     )
       call add_real_4d("lv",      ns, 2, nzeta, ntheta3, lv, order=(/1, 3, 4, 2 /))
       call add_real_4d("lu",      ns, 2, nzeta, ntheta3, lu, order=(/1, 3, 4, 2 /))

       call close_dbg_out()
     end if

  ELSE ! (iequi .eq. 0)

    ! iequi == 1 --> final iter for fileout()

    ! NOTE THAT lu=>czmn, lv=>crmn externally
    !   SO THIS STORES bsupv in czmn_e, bsupu in crmn_e
    ! only executed at end of equilibrium
    lu(:nrzt,0) = bsupv(:nrzt)
    lv(:nrzt,0) = bsupu(:nrzt)

    ! COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
    ! FOR FORCE BALANCE AND RETURN (IEQUI=1)

    ! final call from fileout --> compute additional stuff
    DO js = ns-1,2,-1
       DO l = js, nrzt, ns
          bsubvh(l) = 2*bsubv_e(l) - bsubvh(l+1)
       END DO
    END DO

    ! ADJUST <bsubvh> AFTER MESH-BLENDING
    DO js = 2, ns
       curpol_temp = fpsi(js) - SUM(bsubvh(js:nrzt:ns)*wint(js:nrzt:ns))

       DO l = js, nrzt, ns
          bsubvh(l) = bsubvh(l) + curpol_temp
       END DO
    END DO

    ! NOTE: below four lines modify bsubvh above !!!
    bsubu_e(:nrzt) = bsubuh(:nrzt)
    bsubv_e(:nrzt) = bsubvh(:nrzt)

    bsubu_o(:nrzt) = shalf(:nrzt)*bsubu_e(:nrzt)
    bsubv_o(:nrzt) = shalf(:nrzt)*bsubv_e(:nrzt)

    if (open_dbg_context("bcovar_fileout")) then

      call add_real_3d("lu_e", ns, nzeta, ntheta3, lu(:nrzt,0))
      call add_real_3d("lv_e", ns, nzeta, ntheta3, lv(:nrzt,0))

      call add_real_1d("fpsi", ns-1, fpsi(2:ns))

!       call add_real_3d("bsubvh", ns, nzeta, ntheta3, bsubvh)

      call add_real_3d("bsubu_e", ns, nzeta, ntheta3, bsubu_e)
      call add_real_3d("bsubv_e", ns, nzeta, ntheta3, bsubv_e)

      call add_real_3d("bsubu_o", ns, nzeta, ntheta3, bsubu_o)
      call add_real_3d("bsubv_o", ns, nzeta, ntheta3, bsubv_o)

      call close_dbg_out()
    end if

  END IF ! (iequi .eq. 0)

END SUBROUTINE bcovar
