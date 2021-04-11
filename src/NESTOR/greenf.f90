!> \file
SUBROUTINE greenf(delgr, delgrp, ip)
  USE vacmod
  USE vparams, ONLY: one
  IMPLICIT NONE

  INTEGER, INTENT(in) :: ip
  REAL(rprec), DIMENSION(nuv), INTENT(out) :: delgr, delgrp

  INTEGER, DIMENSION(2) :: ilow, ihigh
  INTEGER :: ivoff, iskip, iuoff, i, kp, nloop, ivoff0
  REAL(rprec):: xip, yip, xper, yper, sxsave, sysave, ftemp, htemp

  ! ON ENTRANCE, IP IS THE INDEX OF THE PRIMED MESH POINT (lies in 1st field period)
  !
  ! ON EXIT, DELGR IS THE DIFFERENCE OF "GREEN'S FUNCTION"
  ! AND ANALYTIC APPROXIMATION, SUMMED OVER ALL FIELD PERIODS
  ! DELGRP IS DIFFERENCE OF DERIVATIVE OF "GREEN'S FUNCTION"
  ! AND ANALYTIC APPROXIMATION.
  !
  ! BOTH THESE QUANTITIES ARE COMPUTED FOR ALL UNPRIMED U,V POINTS IN ONE FIELD PERIOD,
  ! FOR THIS FIXED PRIMED POINT (IP).

  ! COMPUTE OFFSETS FOR U,V ANGLE DIFFERENCES AND CONSTANTS
  ilow(1) = 1
  ilow(2) = ip + 1
  ihigh(1) = ip - 1
  ihigh(2) = nuv
  ivoff0 = nuv + 1 - ip
  iskip = (ip - 1)/nv
  iuoff = nuv - nv*iskip

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
     ivoff = ivoff0 + 2*nu*(kp-1)

     ! x(ip) in field period kp
     xper = xip*cosper(kp) - yip*sinper(kp)

     ! y(ip) in field period kp
     yper = yip*cosper(kp) + xip*sinper(kp)

     sxsave = (snr(ip)*xper - snv(ip)*yper)/r1b(ip)
     sysave = (snr(ip)*yper + snv(ip)*xper)/r1b(ip)

     IF (kp.EQ.1 .OR. nv.EQ.1) THEN
        ! INITIALIZE ANALYTIC APPROXIMATIONS GA1, GA2
        DO i = 1, nuv
           ga1(i) = tanu(i+iuoff)*(  guu_b(ip)*tanu(i+iuoff) + guv_b(ip)*tanv(i+ivoff)) &
                                   + gvv_b(ip)*tanv(i+ivoff)*tanv(i+ivoff)
           ga2(i) = tanu(i+iuoff)*(  auu(ip)*tanu(i+iuoff) + auv(ip)*tanv(i+ivoff)) &
                                   + avv(ip)*tanv(i+ivoff)*tanv(i+ivoff)
        END DO

        DO nloop = 1, 2
           IF (kp.GT.1 .AND. nloop.EQ.2) CYCLE
           DO i = ilow(nloop), ihigh(nloop)
             ga2(i) = ga2(i)/ga1(i)
             ga1(i) = one/SQRT(ga1(i))
             ftemp = one/(gsave(i) - 2*(xper*rcosuv(i) + yper*rsinuv(i)))
             htemp = SQRT(ftemp)
             delgrp(i) = delgrp(i) - ga2(i)*ga1(i) &
                + ftemp*htemp*(rcosuv(i)*sxsave + rsinuv(i)*sysave + dsave(i))
             delgr(i) = delgr(i) + htemp - ga1(i)
           END DO
        END DO

        IF (kp.EQ.nvper .AND. nv.EQ.1) THEN
            delgrp = delgrp/nvper
            delgr  = delgr /nvper
        END IF

        ! ivoff = ivoff + 2*nu
        ihigh(1) = nuv

     ELSE
        DO i = 1,nuv
          ftemp = one/(gsave(i) - 2*(xper*rcosuv(i) + yper*rsinuv(i)))
          htemp = SQRT(ftemp)
          delgrp(i) = delgrp(i) &
             + ftemp*htemp*(rcosuv(i)*sxsave + rsinuv(i)*sysave + dsave(i))
          delgr(i) = delgr(i) + htemp
        END DO
     ENDIF
  END DO

END SUBROUTINE greenf
