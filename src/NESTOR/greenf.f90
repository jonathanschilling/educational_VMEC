!> \file
!> \brief Compute the regularized evaluation of the Green's function and the source term.

!> \brief Compute the regularized evaluation of the Green's function and the source term.
!>
!> @param delgr
!> @param delgrp
!> @param ip
SUBROUTINE greenf(delgr, delgrp, ip)
  USE vacmod
  USE greenf_singular, only: singular_defects

!  use dbgout
  
  IMPLICIT NONE

  INTEGER, INTENT(in) :: ip
  REAL(rprec), DIMENSION(nuv), INTENT(out) :: delgr, delgrp

  INTEGER, DIMENSION(2) :: ilow, ihigh
  INTEGER :: ivoff, iskip, iuoff, i, kp, nloop, ivoff0
  INTEGER :: kua, kva, kv0, kvm, iup, ium, ivp, ivm
  REAL(rprec):: xip, yip, ftemp, htemp ! xper, yper, sxsave, sysave
  REAL(rprec) :: hv, su, du, sv, dv
  REAL(rprec) :: wg(3), wk(3), dphi, area

  REAL(rprec), dimension(nvper) :: xper, yper, sxsave, sysave
  integer,     dimension(nvper,nuv) :: idx_tanu, idx_tanv
  REAL(rprec), dimension(nvper,nuv) :: all_tanu, all_tanv


  ! ON ENTRANCE, IP IS THE INDEX OF THE PRIMED MESH POINT (lies in 1st field period)
  !
  ! ON EXIT:
  ! DELGR  IS THE DIFFERENCE OF "GREEN'S FUNCTION" AND ANALYTIC APPROXIMATION, SUMMED OVER ALL FIELD PERIODS
  ! DELGRP IS THE DIFFERENCE OF DERIVATIVE OF "GREEN'S FUNCTION" AND ANALYTIC APPROXIMATION.
  !
  ! BOTH THESE QUANTITIES ARE COMPUTED FOR ALL UNPRIMED U,V POINTS IN ONE FIELD PERIOD,
  ! FOR THIS FIXED PRIMED POINT (IP).

  ! COMPUTE OFFSETS FOR U,V ANGLE DIFFERENCES AND CONSTANTS

  ! first round: up before singularity
  ilow(1) = 1
  ihigh(1) = ip - 1

  ! second round: from after singularity onwards
  ilow(2) = ip + 1
  ihigh(2) = nuv

  ivoff0  = nuv - (ip - 1)

  iskip  = (ip - 1)/nv ! implicit floor by conversion to integer
  iuoff  = nuv - nv*iskip

  ! x == r*COS(ip), in 1st field period
  xip = rcosuv(ip)

  ! y == r*SIN(ip), in 1st field period
  yip = rsinuv(ip)

  delgr  = 0
  delgrp = 0

  ! COMPUTE FIELD-PERIOD INVARIANT VECTORS
  ! NOTE: |x - x'|**2 = gsave - 2*[x*x' + y*y']
  DO i = 1, nuv
     gsave(i) = rzb2(ip) + rzb2(i) - 2*z1b(ip)*z1b(i)
     dsave(i) = drv(ip) + z1b(i)*snz(ip)
  END DO

  ! SUM OVER FIELD-PERIODS (NVPER=NFPER) OR INTEGRATE OVER NV (NVPER=64) IF NV == 1
  !
  ! NOTE THE SURFACE NORMAL SNORM == Xu cross Xv = NP*[SNR, SNV, SNZ]
  ! IS PERIODIC ON EACH FIELD PERIOD: NOTE THE LOOP OVER KP IS A REDUCTION ON delgr, delgrp
  DO kp = 1, nvper
     ! add in offset due to toroidal module
     ivoff = ivoff0 + 2*nu*(kp-1)

     ! x(ip) in field period kp
     xper(kp) = xip*cosper(kp) - yip*sinper(kp)

     ! y(ip) in field period kp
     yper(kp) = xip*sinper(kp) + yip*cosper(kp)

     sxsave(kp) = (snr(ip)*xper(kp) - snv(ip)*yper(kp))/r1b(ip)
     sysave(kp) = (snr(ip)*yper(kp) + snv(ip)*xper(kp))/r1b(ip)

     IF (kp.EQ.1 .OR. nv.EQ.1) THEN
         ! Tokamak:     always
         ! Stellarator: first toroidal module
!
        ! INITIALIZE ANALYTIC APPROXIMATIONS GA1, GA2
        DO i = 1, nuv

           idx_tanu(kp, i) = i+iuoff
           idx_tanv(kp, i) = i+ivoff
           all_tanu(kp, i) = tanu(i+iuoff)
           all_tanv(kp, i) = tanv(i+ivoff)

           ga1(kp,i) = tanu(i+iuoff)*(  guu_b(ip)*tanu(i+iuoff) + guv_b(ip)*tanv(i+ivoff)) &
                                      + gvv_b(ip)*tanv(i+ivoff)*tanv(i+ivoff)
           ga2(kp,i) = tanu(i+iuoff)*(  auu  (ip)*tanu(i+iuoff) + auv  (ip)*tanv(i+ivoff)) &
                                      + avv  (ip)*tanv(i+ivoff)*tanv(i+ivoff)
        END DO

        DO nloop = 1, 2
           IF (kp.GT.1 .AND. nloop.EQ.2) then
              ! Tokamak (kp>1): only need to skip exactly singular point if in same module
              ! --> first round already goes to nuv, since ihigh(1) was updated below
              CYCLE
           end if

           ! loop over grid points; ilow, ihigh used to skip point at which exact singularity occurs: ip==i
           DO i = ilow(nloop), ihigh(nloop)
             ga2(kp,i) = ga2(kp,i)/ga1(kp,i)
             ga1(kp,i) = one/SQRT(ga1(kp,i))
             ftemp = one/(gsave(i) - 2*(xper(kp)*rcosuv(i) + yper(kp)*rsinuv(i)))
             htemp = SQRT(ftemp)
             delgrp(i) = delgrp(i) + ftemp*htemp*(rcosuv(i)*sxsave(kp) + rsinuv(i)*sysave(kp) + dsave(i)) - ga2(kp,i)*ga1(kp,i)
             delgr(i)  = delgr(i)  +       htemp                                                          -           ga1(kp,i)
           END DO
        END DO

        ! ga1 and ga1*ga2 vanish with a kink where u - u' or v - v' is pi,
        ! where the Green's function is smooth. The trapezoidal sum across a
        ! kink of slope s is off by h**2*s/6, so the nodes on the kink take
        ! h*s/6 in place of zero, and the node where the u and v kinks cross
        ! takes the v value. For odd nv the v kink lies midway between two
        ! nodes, and each of them takes h*s/24 less.
        hv = alv
        if (nv .eq. 1) then
           hv = pi2/nvper
        end if
        su = alu/(24*sqrt(guu_b(ip)))
        du = su*auu(ip)/guu_b(ip)
        sv = hv/(24*sqrt(gvv_b(ip)))
        dv = sv*avv(ip)/gvv_b(ip)
        kua = mod(iskip + nu/2, nu)*nv
        if (nv .eq. 1) then
           if (2*(kp - 1) .eq. nvper) then
              delgr  = delgr  - sv
              delgrp = delgrp - dv
           else
              delgr (kua + 1) = delgr (kua + 1) - su
              delgrp(kua + 1) = delgrp(kua + 1) - du
           end if
        else
           delgr (kua + 1:kua + nv) = delgr (kua + 1:kua + nv) - su
           delgrp(kua + 1:kua + nv) = delgrp(kua + 1:kua + nv) - du
           kv0 = mod(ip - 1, nv)
           if (mod(nv, 2) .eq. 0) then
              kva = mod(kv0 + nv/2, nv) + 1
              delgr (kua + kva) = delgr (kua + kva) + su
              delgrp(kua + kva) = delgrp(kua + kva) + du
              delgr (kva:nuv:nv) = delgr (kva:nuv:nv) - sv
              delgrp(kva:nuv:nv) = delgrp(kva:nuv:nv) - dv
           else
              do kvm = kv0 + (nv - 1)/2, kv0 + (nv + 1)/2
                 kva = mod(kvm, nv) + 1
                 delgr (kva:nuv:nv) = delgr (kva:nuv:nv) + sv/4
                 delgrp(kva:nuv:nv) = delgrp(kva:nuv:nv) + dv/4
                 delgr (kua + kva) = delgr (kua + kva) - sv/4
                 delgrp(kua + kva) = delgrp(kua + kva) - dv/4
              end do
           end if
        end if

        ! The sum skips the source point, where the remainders tend to a
        ! direction-dependent limit and rise linearly; its defects for the
        ! weights 1, u - u' and v - v' go onto the source node and its
        ! neighbours. For nv = 1 the v step is the image spacing in the
        ! toroidal angle, along which the modes do not vary.
        if (kp .eq. 1) then
           kv0 = mod(ip - 1, nv)
           dphi = onp
           if (nv .eq. 1) then
              dphi = one
           end if
           call singular_defects(rtay(:,:,ip), ztay(:,:,ip), xip, yip, z1b(ip), alvp*kv0, dphi, &
                                 sxsave(1), sysave(1), snz(ip), guu_b(ip), guv_b(ip), gvv_b(ip), &
                                 auu(ip), auv(ip), avv(ip), alu, hv, wg, wk)
           area = alu*hv
           iup = ip + nv
           if (iup .gt. nuv) iup = iup - nuv
           ium = ip - nv
           if (ium .lt. 1) ium = ium + nuv
           delgr (ip)  = delgr (ip)  + wg(1)/area
           delgrp(ip)  = delgrp(ip)  + wk(1)/area
           delgr (iup) = delgr (iup) + wg(2)/(2*alu*area)
           delgr (ium) = delgr (ium) - wg(2)/(2*alu*area)
           delgrp(iup) = delgrp(iup) + wk(2)/(2*alu*area)
           delgrp(ium) = delgrp(ium) - wk(2)/(2*alu*area)
           if (nv .gt. 1) then
              ivp = ip + 1
              if (kv0 .eq. nv - 1) ivp = ivp - nv
              ivm = ip - 1
              if (kv0 .eq. 0) ivm = ivm + nv
              delgr (ivp) = delgr (ivp) + wg(3)/(2*hv*area)
              delgr (ivm) = delgr (ivm) - wg(3)/(2*hv*area)
              delgrp(ivp) = delgrp(ivp) + wk(3)/(2*hv*area)
              delgrp(ivm) = delgrp(ivm) - wk(3)/(2*hv*area)
           end if
        end if

        ! update the upper bound of the first loop
        ! --> in Tokamak case, skip exact singularity only if in first toroidal "module"
        ihigh(1) = nuv

     ELSE
        ! Tokamak:     never
        ! Stellarator: all toroidal modules after first one

        DO i = 1,nuv
          ftemp = one/(gsave(i) - 2*(xper(kp)*rcosuv(i) + yper(kp)*rsinuv(i)))
          htemp = SQRT(ftemp)
          delgrp(i) = delgrp(i) + ftemp*htemp*(rcosuv(i)*sxsave(kp) + rsinuv(i)*sysave(kp) + dsave(i))
          delgr(i)  = delgr(i)  +       htemp
        END DO
     ENDIF
  END DO

  IF (nv.EQ.1) THEN
    ! Tokamak: toroidal summation is actually an integration
    ! --> need to divide by step length!
    delgrp = delgrp/nvper
    delgr  = delgr /nvper
  END IF

!   if (open_dbg_context("vac1n_greenf", id=icall)) then
!
!     call add_real_1d("xper",  nvper, xper)
!     call add_real_1d("yper",  nvper, yper)
!     call add_real_1d("sxsave",  nvper, sxsave)
!     call add_real_1d("sysave",  nvper, sysave)
!
!     call add_int_2d("idx_tanu",  nvper, nuv, idx_tanu)
!     call add_int_2d("idx_tanv",  nvper, nuv, idx_tanv)
!     call add_real_2d("all_tanu", nvper, nuv, all_tanu)
!     call add_real_2d("all_tanv", nvper, nuv, all_tanv)
!
!     call close_dbg_out()
!   end if


END SUBROUTINE greenf
