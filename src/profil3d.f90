!> \file
SUBROUTINE profil3d(rmn, zmn, lreset, linterp)
  USE vmec_main
  USE vmec_params
  USE realspace
  USE xstuff
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) ::  rmn, zmn
  LOGICAL, INTENT(in) :: lreset, linterp

  INTEGER :: js, l, lk, lt, lz, ntype, m, n, mn
  REAL(rprec), DIMENSION(0:ntor,ntmax) :: rold, zold
  REAL(rprec) :: sm0, t1, facj, si, rax1, zax1
  INTEGER :: jcount, jk, k

  DO js = 1, ns
     phip(js:nrzt:ns) = phips(js)
     chip(js:nrzt:ns) = chips(js)
  END DO

  phip(nrzt+1) = 0
  faclam = 0
  wint(1:nrzt:ns) = 0

  lk = 0
  DO lt = 1, ntheta3
     DO lz = 1, nzeta
        lk = lk + 1
        DO js=2,ns
           wint(js+ns*(lk-1)) = cosmui3(lt,0)/mscale(0)
        END DO
     END DO
  END DO

  ! COMPUTE ARRAY FOR REFLECTING v = -v (ONLY needed for lasym)
  jcount = 0
  DO k = 1, nzeta
     jk = nzeta + 2 - k
     IF (k .eq. 1) jk = 1
     DO js = 1,ns
       jcount = jcount+1
       ! Index for 2pi-zeta[k]
       ireflect(jcount) = js+ns*(jk-1)
     ENDDO
  END DO

  ! INDEX FOR u = -u (need for lasym integration in wrout)
  lk = 0
  IF (.NOT.ALLOCATED(uminus)) ALLOCATE (uminus(nznt))
  DO lt = 1, ntheta2
     k = ntheta1-lt+2
     IF (lt .eq. 1) then
        ! u=-0 => u=0
        k = 1
     end if
     DO lz = 1, nzeta
        lk = lk + 1
        ! (-u), for u = 0,pi
        uminus(lk) = k
     END DO
  END DO

  ! COMPUTE INITIAL R AND Z FOURIER COEFFICIENTS FROM SCALED BOUNDARY VALUES
  ! AND
  ! SCALXC ARRAY (1/SQRTS FACTOR FOR ODD M VALUES)
  DO js = 1, ns
     si = sqrts(js)*sqrts(js)
     sm0 = one - si
     DO ntype = 1, ntmax
        DO m = 0, mpol1
           DO n = 0, ntor
              t1 = one/(mscale(m)*nscale(n))

              mn = n + ntor1*m
              l = js + ns*mn + (ntype - 1)*mns
              IF (MOD(m,2) .eq. 0) THEN
                 scalxc(l) = one
              ELSE
                 scalxc(l) = one/MAX(sqrts(js),sqrts(2))
              ENDIF

              ! Do not overwrite r,z if read in from wout file AND in free bdy mode
              ! For fixed boundary, edge values MAY have been perturbed, so must execute this loop
              IF (.not.lreset .and. lfreeb) CYCLE
              IF (m .eq. 0) THEN
                 IF (.not.lreset) CYCLE        !Freeze axis if read in from wout file
                 rmn(js,n,m,ntype) = rmn(js,n,m,ntype) + si*(rmn_bdy(n,m,ntype)*t1 - rmn(ns,n,m,ntype))
                 zmn(js,n,m,ntype) = zmn(js,n,m,ntype) + si*(zmn_bdy(n,m,ntype)*t1 - zmn(ns,n,m,ntype))
                 IF (js .eq. 1) THEN
                    rold(n,ntype) = rmn(1,n,0,ntype)
                    zold(n,ntype) = zmn(1,n,0,ntype)
                 ENDIF
                 IF (ntype .eq. rcc) rax1 = raxis_cc(n)
                 IF (ntype .eq. zcs) zax1 =-zaxis_cs(n)
                 IF (ntype .eq. rcs) rax1 =-raxis_cs(n)
                 IF (ntype .eq. zcc) zax1 = zaxis_cc(n)
                 IF (ntype.eq.rcc .or. ntype.eq.rcs) THEN
                    rmn(js,n,m,ntype) = rmn(js,n,m,ntype) + sm0*(rax1*t1 - rold(n,ntype))
                 END IF
                 IF (ntype.eq.zcs .or. ntype.eq.zcc) THEN
                    zmn(js,n,m,ntype) = zmn(js,n,m,ntype) + sm0*(zax1*t1 - zold(n,ntype))
                 END IF
              ELSE
                 ! TURN OFF NEXT 3 LINES IF THIS ONE ACTIVATED
                 facj = sqrts(js)**m

                 ! IF (MOD(m,2) .eq. 0) THEN
                 !     facj = sqrts(js)*sqrts(js)
                 ! ELSE IF (MOD(m,2) .eq. 1) THEN
                 !    facj = sqrts(js)**MIN(m,3)
                 ! END IF
                 rmn(js,n,m,ntype) = rmn(js,n,m,ntype) + (rmn_bdy(n,m,ntype)*t1 - rmn(ns,n,m,ntype))*facj
                 zmn(js,n,m,ntype) = zmn(js,n,m,ntype) + (zmn_bdy(n,m,ntype)*t1 - zmn(ns,n,m,ntype))*facj
              ENDIF
           END DO
        END DO
     END DO
  END DO

  ! Z-components
  scalxc(1+irzloff:2*irzloff)   = scalxc(:irzloff)

  ! Lamda-components
  scalxc(1+2*irzloff:3*irzloff) = scalxc(:irzloff)

  ! STORE PHIFAC IN XC(NEQS1) ARRAY ELEMENT
  ! STORE DEL-RMSE IN XC(NEQS2) ARRAY ELEMENT
  xc(neqs1) = 1
  xc(neqs2) = 0
  scalxc(neqs1) = 1
  scalxc(neqs2) = 1

END SUBROUTINE profil3d