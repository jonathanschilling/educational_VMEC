!> \file
!> \brief Compute the regularized evaluation of the Green's function and the source term.

!> \brief Compute the regularized evaluation of the Green's function and the source term.
!>
!> @param delgr
!> @param delgrp
!> @param ip
SUBROUTINE greenf(delgr, delgrp, ip)
  USE vacmod
  IMPLICIT NONE

  INTEGER, INTENT(in) :: ip
  REAL(rprec), DIMENSION(nuv), INTENT(out) :: delgr, delgrp

  INTEGER, DIMENSION(2) :: ilow, ihigh
  INTEGER :: ivoff, iskip, iuoff, i, kp, nloop, ivoff0
  integer :: ku_i, kv_i, ku_ip, kv_ip, delta_ku, delta_kv
  REAL(rprec):: xip, yip, xper, yper, sxsave, sysave, ftemp, htemp

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

  iskip  = (ip - 1)/nv
  iuoff  = nuv - nv*iskip

  if (nv .eq. 1) then
     ! Tokamak: nvper toroidal "modules"
     ku_ip = (ip-1)/nvper + 1
  else
     ! Stellarator: nv toroidal grid points
     ku_ip = (ip-1)/nv + 1
  end if

  ! x == r*COS(ip), in 1st field period
  xip = rcosuv(ip)

  ! y == r*SIN(ip), in 1st field period
  yip = rsinuv(ip)

!   print *, ip, xip, yip

  delgr  = 0
  delgrp = 0

  ! COMPUTE FIELD-PERIOD INVARIANT VECTORS
  ! NOTE: |x - x'|**2 = gsave - 2*[x*x' + y*y']
  DO i = 1, nuv
     gsave(i) = rzb2(ip) + rzb2(i) - 2*z1b(ip)*z1b(i)
     dsave(i) = drv(ip) + z1b(i)*snz(ip)

     ! log_gsave_dsave.dat
     ! print *, ip, i, gsave(i), dsave(i)
  END DO

  ! SUM OVER FIELD-PERIODS (NVPER=NFPER) OR INTEGRATE OVER NV (NVPER=64) IF NV == 1
  !
  ! NOTE THE SURFACE NORMAL SNORM == Xu cross Xv = NP*[SNR, SNV, SNZ]
  ! IS PERIODIC ON EACH FIELD PERIOD: NOTE THE LOOP OVER KP IS A REDUCTION ON delgr, delgrp
  DO kp = 1, nvper
     ! add in offset due to toroidal module
     ivoff = ivoff0 + 2*nu*(kp-1)

     ! x(ip) in field period kp
     xper = xip*cosper(kp) - yip*sinper(kp)

     ! y(ip) in field period kp
     yper = yip*cosper(kp) + xip*sinper(kp)

     sxsave = (snr(ip)*xper - snv(ip)*yper)/r1b(ip)
     sysave = (snr(ip)*yper + snv(ip)*xper)/r1b(ip)

     ! log_xper_yper_sxsave_sysave.dat
!      print *, ip, kp, cosper(kp), sinper(kp), xper, yper, sxsave, sysave

     IF (kp.EQ.1 .OR. nv.EQ.1) THEN

!          if (ip .eq. 1) then
!            print *, "period ", kp
!          end if

         if (nv .eq. 1) then
            ! Tokamak: nvper toroidal "modules"
            !kv_ip = mod(ip-1, nvper) + 1
            kv_ip = kp
         else
            ! Stellarator: nv toroidal grid points
            kv_ip = mod(ip-1, nv) + 1
         end if

        ! INITIALIZE ANALYTIC APPROXIMATIONS GA1, GA2

        !if (ip .eq. 1) write(*,*)"#         ip,         kp,     ivoff0,      ivoff,      iskip,      iuoff"
        !print *,  ip, kp, ivoff0, ivoff, iskip, iuoff

        DO i = 1, nuv

          if (nv .eq. 1) then
              ! Tokamak: nvper toroidal "modules"
              ku_i = (i-1)/nvper + 1
              kv_i = mod(i-1, nvper) + 1

              delta_kv = mod(kv_i-kv_ip+nvper, nvper) + 1
           else
              ! Stellarator: nv toroidal grid points
              ku_i =    (i-1)/nv  + 1
              kv_i = mod(i-1, nv) + 1

              delta_kv = mod(kv_i-kv_ip+nv, nv) + 1
           end if

           delta_ku = ku_i-ku_ip+nu + 1

           ! this is where log_greenf_2.txt comes from
           !print *, ip, i, iuoff, ivoff, i+iuoff, i+ivoff, tanu(i+iuoff), tanv(i+ivoff)

!            ga1(i) = tanu(i+iuoff)*(  guu_b(ip)*tanu(i+iuoff) + guv_b(ip)*tanv(i+ivoff)) &
!                                    + gvv_b(ip)*tanv(i+ivoff)*tanv(i+ivoff)
!            ga2(i) = tanu(i+iuoff)*(  auu  (ip)*tanu(i+iuoff) + auv  (ip)*tanv(i+ivoff)) &
!                                    + avv  (ip)*tanv(i+ivoff)*tanv(i+ivoff)

           if (abs(tanu(i+iuoff) - tanu_1d(delta_ku)) .gt. 1.0e-9_dp) then
              print *, "mismatch in tanu:"
              print *, "  tanu(i+iuoff)    = ", tanu(i+iuoff)
              print *, "  tanu_1d(delta_ku)= ", tanu_1d(delta_ku)
              print *, "        difference = ", tanu(i+iuoff) - tanu_1d(delta_ku)
              print *, "                ip = ", ip
              print *, "                i  = ", i
              print *, "            iskip  = ", iskip
              print *, "            iuoff  = ", iuoff
              print *, "             ku_i  = ", ku_i
              print *, "             kv_i  = ", kv_i
              print *, "            ku_ip  = ", ku_ip
              print *, "            kv_ip  = ", kv_ip
              print *, "         delta_ku  = ", delta_ku
              stop
           end if


           if (abs(tanv(i+ivoff) - tanv_1d(delta_kv)) .gt. 1.0e-9_dp) then
              print *, "mismatch in tanv:"
              print *, "  tanv(i+ivoff)    = ", tanv(i+ivoff)
              print *, "  tanv_1d(delta_kv)= ", tanv_1d(delta_kv)
              print *, "        difference = ", tanv(i+iuoff) - tanv_1d(delta_kv)
              print *, "                ip = ", ip
              print *, "                i  = ", i
              print *, "           ivoff0  = ", ivoff0
              print *, "           ivoff   = ", ivoff
              print *, "             ku_i  = ", ku_i
              print *, "             kv_i  = ", kv_i
              print *, "            ku_ip  = ", ku_ip
              print *, "            kv_ip  = ", kv_ip
              print *, "         delta_kv  = ", delta_kv
              stop
           end if



           ga1(i) = tanu_1d(delta_ku)*(  guu_b(ip)*tanu_1d(delta_ku) + guv_b(ip)*tanv_1d(delta_kv)) &
                                   + gvv_b(ip)*tanv_1d(delta_kv)*tanv_1d(delta_kv)
           ga2(i) = tanu_1d(delta_ku)*(  auu  (ip)*tanu_1d(delta_ku) + auv  (ip)*tanv_1d(delta_kv)) &
                                   + avv  (ip)*tanv_1d(delta_kv)*tanv_1d(delta_kv)

           ! This was used to generate log_ga1_ga2.dat
           ! print *, ip, i, ga1(i), ga2(i)
        END DO

        DO nloop = 1, 2
           IF (kp.GT.1 .AND. nloop.EQ.2) then
              ! Tokamak (kp>1): only need to skip exactly singular point if in same module
              ! --> first round already goes to nuv, since ihigh(1) was updated below
              CYCLE
           end if

           ! loop over grid points; ilow, ihigh used to skip point at which exact singularity occurs: ip==i
           DO i = ilow(nloop), ihigh(nloop)
             ga2(i) = ga2(i)/ga1(i)
             ga1(i) = one/SQRT(ga1(i))
             ftemp = one/(gsave(i) - 2*(xper*rcosuv(i) + yper*rsinuv(i)))
             htemp = SQRT(ftemp)
             delgrp(i) = delgrp(i) + ftemp*htemp*(rcosuv(i)*sxsave + rsinuv(i)*sysave + dsave(i)) - ga2(i)*ga1(i)
             delgr(i)  = delgr(i)  +       htemp                                                  -        ga1(i)

             ! log_ftemp_htemp_green_greenp.dat
             ! print *, ip, i, ga1(i), ga2(i), ftemp, htemp, delgr(i), delgrp(i)

           END DO

!            if (nloop.eq.1) then
!               print *, ip, ip, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
!            end if
        END DO

        IF (kp.EQ.nvper .AND. nv.EQ.1) THEN
            ! Tokamak: toroidal summation is actually an integration
            ! --> need to divide by step length!
            delgrp = delgrp/nvper
            delgr  = delgr /nvper
        END IF

        ! this is not needed anymore since ivoff is computed directly from ivoff0 and kp above
        ! ivoff = ivoff + 2*nu

        ! update the upper bound of the first loop
        ! --> in Tokamak case, skip exact singularity only if in first toroidal "module"
        ihigh(1) = nuv

     ELSE
!         if (ip.eq.1) then
!           print *, "period ",kp
!         end if

        DO i = 1,nuv

          ! log_otherPeriods_gsave_green_etc.dat
          ! print *, ip, kp, i, gsave(i), dsave(i), xper, yper, rcosuv(i), rsinuv(i), sxsave, sysave, delgr(i), delgrp(i)

          ftemp = one/(gsave(i) - 2*(xper*rcosuv(i) + yper*rsinuv(i)))
          htemp = SQRT(ftemp)
          delgrp(i) = delgrp(i) + ftemp*htemp*(rcosuv(i)*sxsave + rsinuv(i)*sysave + dsave(i))
          delgr(i)  = delgr(i)  +       htemp

          ! print *, ip, i, ftemp, htemp, delgr(i), delgrp(i)

        END DO
     ENDIF
  END DO

END SUBROUTINE greenf
