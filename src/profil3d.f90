!> \file
!> \brief Compute three-dimensional profiles of flux-surface geometry etc.

!> \brief Compute three-dimensional profiles of flux-surface geometry etc.
!>
!> @param rmn Fourier coefficients of \f$R\f$
!> @param zmn Fourier coefficients of \f$R\f$
!> @param lreset flag to indicate the geometry of the LCFS (and axis ?) should be used to interpolate into the plasma volume
SUBROUTINE profil3d(rmn, zmn, lreset)
  USE vmec_main
  USE vmec_params
  USE realspace
  USE xstuff
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) ::  rmn, zmn
  LOGICAL, INTENT(in) :: lreset

  INTEGER :: js, l, lk, lt, lz, ntype, m, n, mn
  REAL(rprec), DIMENSION(0:ntor,ntmax) :: rold, zold
  REAL(rprec) :: sm0, t1, facj, si, rax1, zax1
  INTEGER :: jcount, jk, k

  character(len=255) :: dump_filename
  logical, parameter :: dump_profil3d = .true.

  ! expant to full surface grid
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
     jk = nzeta + 2 - k ! (nzeta+1) - (k-1)
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
     k = ntheta1-lt+2 ! (ntheta1+1) - (lt-1)
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

  print *, "initial guess for R_mn, Z_mn in profil3d"

!   write(*,*) "raxis_cc(0)=",raxis_cc(0)
!   DO m = 0, mpol1
!      DO n = 0, ntor
!        if (abs(rmn_bdy(n,m,rcc)) .gt. 1.0e-12_dp) then
!          write(*,*) n, m, rmn_bdy(n,m,rcc)
!        end if
!        if (abs(zmn_bdy(n,m,zcs)) .gt. 1.0e-12_dp) then
!          write(*,*) n, m, zmn_bdy(n,m,zcs)
!        end if
!      end do
!    end do



  DO js = 1, ns

     si = sqrts(js)*sqrts(js) !   s(js) --> ==1 at boundary
     sm0 = one - si           ! 1-s(js) --> ==1 at magn. axis

     DO ntype = 1, ntmax
        DO m = 0, mpol1
           DO n = 0, ntor


              !-----------------------!
              ! compute scalxc(js, m) !
              !-----------------------!
              mn = n + ntor1*m ! linear index over Fourier modes
              l = js + ns*mn + (ntype - 1)*mns ! linear index over all Fourier coefficents on all surfaces
              IF (MOD(m,2) .eq. 0) THEN
                 ! m is even
                 scalxc(l) = one
              ELSE
                 ! m is odd

                 ! make sure to use at least sqrts(2)==1/(ns-1) as normalization/scaling factor,
                 ! since sqrts(1)==0 (see profil1d)
                 scalxc(l) = one/MAX(sqrts(js), sqrts(2))
              ENDIF



              ! ----------------------------------------------------------!
              ! extrapolate R_mn, Z_mn from axis and boundary into volume !
              ! ----------------------------------------------------------!
              t1 = one/(mscale(m)*nscale(n)) ! --> divide out mscale, nscale from user input

              ! Do not overwrite r,z if read in from wout file AND in free bdy mode
              ! For fixed boundary, edge values MAY have been perturbed, so must execute this loop
              IF (.not.lreset .and. lfreeb) CYCLE

              ! below code segment does the extrapolation
              ! of the boundary Fourier coefficients into the plasma volume

              IF (m .eq. 0) THEN

                 IF (.not.lreset) CYCLE        !Freeze axis if read in from wout file
                 ! above instruction cycles all n for m=0 --> skip axis, as said above!

                 ! subtraction of the edge value of rmn, zmn is probably left-over
                 ! from the restart feature from a previous wout file
                 !if (abs(rmn(ns,n,m,ntype)) .ne. zero) then
                 !  print *, "spurious remnant in rmn: ", rmn(ns,n,m,ntype)
                 if (abs(rmn(js,n,m,ntype)) .ne. zero) then
                   print *, "spurious remnant in rmn: ", rmn(js,n,m,ntype)
                   stop
                 end if
                 !if (abs(zmn(ns,n,m,ntype)) .ne. zero) then
                 !  print *, "spurious remnant in zmn: ", zmn(ns,n,m,ntype)
                 if (abs(zmn(js,n,m,ntype)) .ne. zero) then
                   print *, "spurious remnant in zmn: ", zmn(js,n,m,ntype)
                   stop
                 end if
                 !rmn(js,n,m,ntype) = rmn(js,n,m,ntype) + si*(rmn_bdy(n,m,ntype)*t1 - rmn(ns,n,m,ntype))
                 !zmn(js,n,m,ntype) = zmn(js,n,m,ntype) + si*(zmn_bdy(n,m,ntype)*t1 - zmn(ns,n,m,ntype))

                 ! first contribution: boundary scaled into volume
                 rmn(js,n,m,ntype) = si * rmn_bdy(n,m,ntype)*t1
                 zmn(js,n,m,ntype) = si * zmn_bdy(n,m,ntype)*t1

                 !IF (js .eq. 1) THEN
                 !   ! rold, zold will be subtracted
                 !   ! --> this zeroes out any contribution as extrapolated from the boundary ?
                 !   ! since js=1 is handled first, this gets assigned once and is available then
                 !   ! on the other hand, there should be nothing in here (as checked in above spurious... tests)
                 !   ! --> ignore!
                 !   ! --> also for js=1, there is no contribution from above, since si=0 for js=1 !
                 !   rold(n,ntype) = rmn(1,n,0,ntype)
                 !   zold(n,ntype) = zmn(1,n,0,ntype)
                 !
                 !   if (abs(rold(n,ntype)) .ne. zero) then
                 !     print *, "remnant in rold: ", rold(n,ntype)
                 !     stop
                 !   end if
                 !
                 !   if (abs(zold(n,ntype)) .ne. zero) then
                 !     print *, "remnant in zold: ", zold(n,ntype)
                 !     stop
                 !   end if
                 !
                 !ENDIF

                 ! second contribution: axis scaled towards boundary
                 IF (ntype .eq. rcc) rax1 = raxis_cc(n)
                 IF (ntype .eq. rcs) rax1 =-raxis_cs(n)
                 IF (ntype.eq.rcc .or. ntype.eq.rcs) THEN
                    !rmn(js,n,m,ntype) = rmn(js,n,m,ntype) + sm0*(rax1*t1 - rold(n,ntype))
                    rmn(js,n,m,ntype) = rmn(js,n,m,ntype) + sm0 * rax1*t1
                 END IF

                 IF (ntype .eq. zcs) zax1 =-zaxis_cs(n)
                 IF (ntype .eq. zcc) zax1 = zaxis_cc(n)
                 IF (ntype.eq.zcs .or. ntype.eq.zcc) THEN
                    !zmn(js,n,m,ntype) = zmn(js,n,m,ntype) + sm0*(zax1*t1 - zold(n,ntype))
                    zmn(js,n,m,ntype) = zmn(js,n,m,ntype) + sm0 * zax1*t1
                 END IF

              ELSE ! m != 0
                 ! TURN OFF BELOW LINES IF THIS ONE ACTIVATED
                 facj = sqrts(js)**m

                 ! IF (MOD(m,2) .eq. 0) THEN
                 !    facj = sqrts(js)*sqrts(js)
                 ! ELSE IF (MOD(m,2) .eq. 1) THEN
                 !    facj = sqrts(js)**MIN(m,3)
                 ! END IF

                 ! subtraction of the edge value of rmn, zmn is probably left-over
                 ! from the restart feature from a previous wout file
                 !if (abs(rmn(ns,n,m,ntype)) .ne. zero) then
                 !  print *, "spurious remnant in rmn: ", rmn(ns,n,m,ntype)
                 if (abs(rmn(js,n,m,ntype)) .ne. zero) then
                   print *, "spurious remnant in rmn: ", rmn(js,n,m,ntype)
                   stop
                 end if
                 !if (abs(zmn(ns,n,m,ntype)) .ne. zero) then
                 !  print *, "spurious remnant in zmn: ", zmn(ns,n,m,ntype)
                 if (abs(zmn(js,n,m,ntype)) .ne. zero) then
                   print *, "spurious remnant in zmn: ", zmn(js,n,m,ntype)
                   stop
                 end if
                 !rmn(js,n,m,ntype) = rmn(js,n,m,ntype) + (rmn_bdy(n,m,ntype)*t1 - rmn(ns,n,m,ntype))*facj
                 !zmn(js,n,m,ntype) = zmn(js,n,m,ntype) + (zmn_bdy(n,m,ntype)*t1 - zmn(ns,n,m,ntype))*facj

                 rmn(js,n,m,ntype) = facj * rmn_bdy(n,m,ntype)*t1
                 zmn(js,n,m,ntype) = facj * zmn_bdy(n,m,ntype)*t1
              ENDIF



           END DO
        END DO
     END DO
  END DO

  ! distribute scalxc content from R to Z     components
  scalxc(1+  irzloff:2*irzloff) = scalxc(:irzloff)

  ! distribute scalxc content from R to Lamda components
  scalxc(1+2*irzloff:3*irzloff) = scalxc(:irzloff)

  ! dump all relevant output to a text file
  if (dump_profil3d) then
    write(dump_filename, 999) ns, profil3d_calls, trim(input_extension)
999 format('profil3d_',i5.5,'_',i2.2,'.',a)

    open(unit=42, file=trim(dump_filename), status="unknown")

    write(42, *) "# ns ntmax mpol ntor"
    write(42, *) ns, ntmax, mpol, ntor

    write(42, *) "# js ntype m n mn l scalxc"
    DO js = 1, ns
      DO ntype = 1, ntmax
        DO m = 0, mpol1
          DO n = 0, ntor
            mn = n + ntor1*m
            l = js + ns*mn + (ntype - 1)*mns

            write (42, *) js, ntype, m, n, mn, l, scalxc(l)
          end do
        end do
      end do
    end do

    write(42, *) "# js ntype m n rmn zmn"
    DO js = 1, ns
      DO ntype = 1, ntmax
        DO m = 0, mpol1
          DO n = 0, ntor
            write (42, *) js, ntype, m, n, &
                          rmn(js,n,m,ntype), zmn(js,n,m,ntype)
          end do
        end do
      end do
    end do

    close(42)

    print *, "dumped profil3d output to '"//trim(dump_filename)//"'"

  end if


  profil3d_calls = profil3d_calls + 1

END SUBROUTINE profil3d
