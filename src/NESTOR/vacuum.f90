!> \file
!> \brief Compute the vacuum contribution to the free-boundary energy functional.

!> \brief Compute the vacuum contribution to the free-boundary energy functional.
!>
!> @param rmnc
!> @param rmns
!> @param zmns
!> @param zmnc
!> @param xm
!> @param xn
!> @param plascur
!> @param rbtor
!> @param wint
!> @param ivac_skip
!> @param ivac
!> @param mnmax
!> @param ier_flag
!> @param lasym
!> @param signgs
!> @param raxis
!> @param zaxis
SUBROUTINE vacuum(rmnc, rmns, zmns, zmnc, xm, xn,             &
                  plascur, rbtor, wint, ivac_skip, ivac,  &
                  mnmax, ier_flag, lasym, signgs,             &
                  raxis, zaxis)
  USE vacmod
  USE vmec_params, ONLY: norm_term_flag, phiedge_error_flag
  use vmec_main, only: input_extension

  use dbgout
  use vmec_main, only: num_eqsolve_retries

  IMPLICIT NONE

  INTEGER, intent(in) :: ivac_skip, mnmax
  integer, intent(inout) :: ivac, ier_flag
  REAL(rprec), intent(in) :: plascur, rbtor
  REAL(rprec), DIMENSION(mnmax), INTENT(in) :: rmnc, rmns, zmns, zmnc, xm, xn
  REAL(rprec), DIMENSION(nuv2), INTENT(in) :: wint
  logical, intent(in) :: lasym
  real(rprec), intent(in) :: signgs
  real(rprec), dimension(nv), intent(in) :: raxis, zaxis

  INTEGER :: mn, n, n1, m, i, info
  REAL(rprec), DIMENSION(:), POINTER :: potcos, potsin
  REAL(rprec):: dn2, dm2, cosmn, sinmn, huv, hvv, det, bsupu, bsupv, bsubuvac, fac
  logical :: vac1n_solver_active

  ! THIS ROUTINE COMPUTES .5 * B**2 ON THE VACUUM / PLASMA SURFACE
  ! BASED ON THE PROGRAM BY P. MERKEL [J. Comp. Phys. 66, 83 (1986)]
  ! AND MODIFIED BY W. I. VAN RIJ AND S. P. HIRSHMAN (1987)

  ! THE USER MUST SUPPLY THE FILE << MGRID >> WHICH INCLUDES THE MAGNETIC
  ! FIELD DATA TO BE READ BY THE SUBROUTINE BECOIL

  ier_flag = norm_term_flag

  raxis_nestor = raxis
  zaxis_nestor = zaxis

  IF (.not.ALLOCATED(potvac)) STOP 'POTVAC not ALLOCATED in VACCUM'

  ! INDEX OF LOCAL VARIABLES
  !
  ! rmnc,rmns,zmns,zmnc:     Surface Fourier coefficients (m,n) of R,Z
  ! xm,xn:     m, n values corresponding to rc,zs array
  ! bsqvac:    B**2/2 at the vacuum INTERFACE
  ! plascur:   net toroidal current
  ! rbtor  :   net (effective) poloidal current (loop integrated R*Btor)
  ! mnmax:     number of R, Z modes in Fourier series of R,Z
  ! ivac_skip: regulates whether full (=0) or incremental (>0)
  !            update of matrix elements is necessary

  ! compute and store mean magnetic fields (due to
  ! toroidal plasma current and EXTERNAL tf-coils)
  ! note: these are fixed for a constant current iteration
  !
  ! bfield = rbtor*grad(zeta) + plascur*grad("theta") - grad(potential)
  !
  ! where "theta" is computed using Biot-Savart law for filaments
  ! Here, the potential term is needed to satisfy B dot dS = 0 and has the form:
  !
  ! potential = SUM potsin*SIN(mu - nv) + potcos*COS(mu - nv)

  ! write inputs to NESTOR
  if (open_dbg_context("vac1n_vacuum", num_eqsolve_retries)) then

    call add_real_1d("rmnc", mnmax, rmnc)
    call add_real_1d("zmns", mnmax, zmns)
    if (lasym) then
      call add_real_1d("rmns", mnmax, rmns)
      call add_real_1d("zmnc", mnmax, zmnc)
    else
      call add_null("rmns")
      call add_null("zmnc")
    end if
    call add_real_1d("xm", mnmax, xm)
    call add_real_1d("xn", mnmax, xn)

    call add_real("plascur", plascur)
    call add_real("rbtor",   rbtor)

    call add_real_2d("wint", nv, nu3, wint)

    call add_int("ivac_skip", ivac_skip)
    call add_int("ivac",      ivac)
    call add_int("mnmax",     mnmax)
    call add_int("ier_flag",  ier_flag)
    call add_logical("lasym", lasym)
    call add_real("signgs",   signgs)

    call add_real_1d("raxis_nestor", nv, raxis_nestor)
    call add_real_1d("zaxis_nestor", nv, zaxis_nestor)

    call close_dbg_out()
  end if

   IF (.not. precal_done) then
      CALL precal
   end if
   CALL surface (rmnc, rmns, zmns, zmnc, xm, xn, mnmax, lasym, signgs)
   CALL bextern (plascur, wint)

   ! NOTE: all fine up to here against NESTOR.py

   ! Determine scalar magnetic potential POTVAC
   CALL scalpot (potvac, amatrix, wint, ivac_skip, lasym, m_map_wrt, n_map_wrt)

   ! stand-alone for debugging: working on scalpot at the moment
   ! return

   vac1n_solver_active = open_dbg_context("vac1n_solver", num_eqsolve_retries)
   if (vac1n_solver_active) then
     call add_real_2d("amatrix", mnpd2, mnpd2, amatrix)
     call add_real_1d("potvac_in", mnpd2, potvac)
   end if

   CALL solver (amatrix, potvac, mnpd2, 1, info)
   IF (info .ne. 0) STOP 'Error in solver in VACUUM'

   potsin => potvac(1:mnpd)
   potcos => potvac(1+mnpd:)

   if (vac1n_solver_active) then
     call add_real_1d("potvac_out", mnpd2, potvac)

     call close_dbg_out()
   end if

   ! compute tangential covariant (sub u,v) and contravariant
   ! (super u,v) magnetic field components on the plasma surface
   potu(:nuv2) = zero
   potv(:nuv2) = zero

   mn = 0
   DO n = -nf, nf
      dn2 = -(n*nfper)
      n1 = ABS(n)
      DO m = 0, mf
         mn = mn + 1
         dm2 = m
         DO i = 1, nuv2
            cosmn = cosu1(i,m)*cosv1(i,n1) + csign(n)*sinu1(i,m)*sinv1(i,n1)
            potu(i) = potu(i) + dm2*potsin(mn)*cosmn
            potv(i) = potv(i) + dn2*potsin(mn)*cosmn
         END DO
         IF (lasym) then
            DO i = 1, nuv2
               sinmn = sinu1(i,m)*cosv1(i,n1) - csign(n)*cosu1(i,m)*sinv1(i,n1)
               potu(i) = potu(i) - dm2*potcos(mn)*sinmn
               potv(i) = potv(i) - dn2*potcos(mn)*sinmn
            END DO
         end if
      END DO
   END DO

   DO i = 1, nuv2
      ! Covariant components
      bsubu(i) = potu(i) + bexu(i)
      bsubv(i) = potv(i) + bexv(i)

      huv = p5*guv_b(i)*(nfper)
      hvv = gvv_b(i)*(nfper*nfper)
      det = one/(guu_b(i)*hvv-huv*huv)

      ! Contravariant components
      bsupu = (hvv*bsubu(i)-huv*bsubv(i))*det
      bsupv = ((-huv*bsubu(i))+guu_b(i)*bsubv(i))*det

      ! .5*|Bvac|**2
      bsqvac(i) = p5*(bsubu(i)*bsupu + bsubv(i)*bsupv)

!       if (ivac.eq.0) then
!         print *, i, bsqvac(i)
!       end if

      ! cylindrical components of vacuum magnetic field
      brv(i)   = rub(i)*bsupu + rvb(i)*bsupv
      bphiv(i) =                r1b(i)*bsupv
      bzv(i)   = zub(i)*bsupu + zvb(i)*bsupv
   END DO

   if (open_dbg_context("vac1n_bsqvac", num_eqsolve_retries)) then

    call add_real_1d("potsin", mnpd, potsin)
    if (lasym) then
      call add_real_1d("potcos", mnpd, potcos)
    else
      call add_null("potcos")
    end if

    call add_real_2d("potu", nv, nu3, potu)
    call add_real_2d("potv", nv, nu3, potv)

    call add_real_2d("bsubu", nv, nu3, bsubu)
    call add_real_2d("bsubv", nv, nu3, bsubv)

    call add_real_2d("bsqvac", nv, nu3, bsqvac)

    call add_real_2d("brv",   nv, nu3, brv)
    call add_real_2d("bphiv", nv, nu3, bphiv)
    call add_real_2d("bzv",   nv, nu3, bzv)

    call close_dbg_out()
  end if

   ! PRINT OUT VACUUM PARAMETERS
   IF (ivac .eq. 0) THEN
      ivac = ivac + 1

      WRITE (*, 200) nfper, mf, nf, nu, nv
      ! WRITE (nthreed, 200) nfper, mf, nf, nu, nv
 200 FORMAT(/,2x,'In VACUUM, np =',i3,2x,'mf =',i3,2x,'nf =',i3,' nu =',i3,2x,'nv = ',i4)

      ! -plasma current/pi2
      bsubuvac = SUM(bsubu(:nuv2)*wint(:nuv2))*signgs*pi2
      bsubvvac = SUM(bsubv(:nuv2)*wint(:nuv2))

      fac = 1.e-6_dp/mu0 ! currents in MA
      WRITE (*,1000) bsubuvac*fac, plascur*fac, bsubvvac, rbtor
      ! WRITE (nthreed, 1000)       bsubuvac*fac, plascur*fac, bsubvvac, rbtor
 1000 FORMAT(2x,'2*pi * a * -BPOL(vac) = ',1p,e10.2,                 &
         ' TOROIDAL CURRENT = ',e10.2,/,2x,'R * BTOR(vac) = ',       &
         e10.2,' R * BTOR(plasma) = ',e10.2)

      IF (rbtor*bsubvvac .lt. zero) THEN
         ! rbtor and bsubvvac must have the same sign
         ier_flag = phiedge_error_flag
      ENDIF

      IF (ABS((plascur - bsubuvac)/rbtor) .gt. 1.e-2_dp) THEN
!             ier_flag = 10 ! 'VAC-VMEC I_TOR MISMATCH : BOUNDARY MAY ENCLOSE EXT. COIL'
         print *, 'VAC-VMEC I_TOR MISMATCH : BOUNDARY MAY ENCLOSE EXT. COIL'
      ENDIF

   ENDIF

   !icall = icall + 1

END SUBROUTINE vacuum
