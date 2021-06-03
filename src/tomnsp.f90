!> \file
!> \brief Fourier-transform forces from real space to Fourier space

!> \brief Fourier-transform symmetric forces from real space to Fourier space
!>
!> @param frzl_array
!> @param armn
!> @param brmn
!> @param crmn
!> @param azmn
!> @param bzmn
!> @param czmn
!> @param blmn
!> @param clmn
!> @param arcon
!> @param azcon
SUBROUTINE tomnsps(frzl_array,       &
                   armn, brmn, crmn, &
                   azmn, bzmn, czmn, &
                   blmn, clmn,       &
                   arcon, azcon       )
  USE vmec_main
  USE vmec_params, ONLY: jlam, jmin2, ntmax, rcc, rss, zsc, zcs
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax), TARGET, INTENT(out) :: frzl_array
  REAL(rprec), DIMENSION(ns*nzeta*ntheta3,0:1), INTENT(in) ::       &
     armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon

  INTEGER :: jmax, m, mparity, i, n, k, l, nsz
  INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl, js
  REAL(rprec), DIMENSION(:,:,:), POINTER :: frcc, frss, fzcs, fzsc, flcs, flsc
  REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: work1
  REAL(rprec), DIMENSION(:), ALLOCATABLE   :: tempr, tempz

  character(len=255) :: dump_filename
  logical            :: dump_tomnsps = .true.

  frcc => frzl_array(:,:,:,rcc)               !!COS(mu) COS(nv)
  fzsc => frzl_array(:,:,:,zsc+ntmax)         !!SIN(mu) COS(nv)
  flsc => frzl_array(:,:,:,zsc+2*ntmax)       !!SIN(mu) COS(nv)
  IF (lthreed) THEN
     frss => frzl_array(:,:,:,rss)            !!SIN(mu) SIN(nv)
     fzcs => frzl_array(:,:,:,zcs+ntmax)      !!COS(mu) SIN(nv)
     flcs => frzl_array(:,:,:,zcs+2*ntmax)    !!COS(mu) SIN(nv)
  END IF

  nsz = ns*nzeta

  ALLOCATE (work1(nsz,12), tempr(nsz), tempz(nsz), stat=i)
  IF (i .ne. 0) STOP 'Allocation error in VMEC2000 tomnsps'

  ioff = LBOUND(frcc,2)
  joff = LBOUND(frcc,3)

  frzl_array = 0

  jmax = ns
  IF (ivac .lt. 1) jmax = ns1

  ! BEGIN FOURIER TRANSFORM
  !
  ! FRmn = ARmn - d(BRmn)/du + d(CRmn)/dv
  ! FZmn = AZmn - d(BZmn)/du + d(CZmn)/dv
  ! FLmn =      - d(BLmn)/du + d(CLmn)/dv
  !
  ! NOTE: sinmumi = -m sin(mu),  sinnvn = -n sin(nv)
  DO m = 0, mpol1

     mparity = MOD(m,2)
     mj = m+joff

     j2 = jmin2(m)
     jl = jlam(m)

     work1 = 0

     ! DO THETA (U) INTEGRATION FIRST ON HALF INTERVAL (0 < U < PI)
     l = 0
     DO i = 1, ntheta2
        jll = l+1   ! start of poloidal slice
        nsl = l+nsz ! end of poloidal slice
        l   = l+nsz ! jump to next poloidal slice

        tempr(:) = armn(jll:nsl,mparity) &
! #ifndef _HBANGLE
            + xmpq(m,1)*arcon(jll:nsl,mparity)
! #end /* ndef _HBANGLE */

        tempz(:) = azmn(jll:nsl,mparity) &
! #ifndef _HBANGLE
            + xmpq(m,1)*azcon(jll:nsl,mparity)
! #end /* ndef _HBANGLE */

        work1(:,1) = work1(:,1) + tempr(:)*cosmui(i,m) + brmn(jll:nsl,mparity)*sinmumi(i,m)
        work1(:,7) = work1(:,7) + tempz(:)*sinmui(i,m) + bzmn(jll:nsl,mparity)*cosmumi(i,m)
        work1(:,11)= work1(:,11)+ blmn(jll:nsl,mparity)*cosmumi(i,m)

        IF (.not.lthreed) CYCLE

        work1(:,2) = work1(:,2) - crmn(jll:nsl,mparity)*cosmui(i,m)
        work1(:,3) = work1(:,3) + tempr(:)*sinmui(i,m) + brmn(jll:nsl,mparity)*cosmumi(i,m)
        work1(:,4) = work1(:,4) - crmn(jll:nsl,mparity)*sinmui(i,m)
        work1(:,5) = work1(:,5) + tempz(:)*cosmui(i,m) + bzmn(jll:nsl,mparity)*sinmumi(i,m)
        work1(:,6) = work1(:,6) - czmn(jll:nsl,mparity)*cosmui(i,m)

        work1(:,8) = work1(:,8) - czmn(jll:nsl,mparity)*sinmui(i,m)
        work1(:,9) = work1(:,9) + blmn(jll:nsl,mparity)*sinmumi(i,m)
        work1(:,10) =work1(:,10)- clmn(jll:nsl,mparity)*cosmui(i,m)

        work1(:,12) =work1(:,12)- clmn(jll:nsl,mparity)*sinmui(i,m)
     END DO

     ! NEXT, DO ZETA (V) TRANSFORM
     DO n = 0, ntor
        ni = n+ioff
        l = 0
        DO k = 1, nzeta
           j2l   = l+j2   ! start of radial slice for R,Z
           jmaxl = l+jmax ! end of radial slice for R,Z

           jll   = l+jl ! start of radial slice for lambda
           nsl   = l+ns ! end of radial slice for lambda

           l = l+ns ! jump to next radial slice

           frcc(j2:jmax,ni,mj) = frcc(j2:jmax,ni,mj) + work1(j2l:jmaxl,1)*cosnv(k,n)
           fzsc(j2:jmax,ni,mj) = fzsc(j2:jmax,ni,mj) + work1(j2l:jmaxl,7)*cosnv(k,n)
           flsc(jl:ns,ni,mj) = flsc(jl:ns,ni,mj) + work1(jll:nsl,11)*cosnv(k,n)

           IF (.not.lthreed) CYCLE

           frcc(j2:jmax,ni,mj) = frcc(j2:jmax,ni,mj) + work1(j2l:jmaxl,2)*sinnvn(k,n)
           fzsc(j2:jmax,ni,mj) = fzsc(j2:jmax,ni,mj) + work1(j2l:jmaxl,8)*sinnvn(k,n)
           frss(j2:jmax,ni,mj) = frss(j2:jmax,ni,mj) + work1(j2l:jmaxl,3)*sinnv(k,n)  &
                                                     + work1(j2l:jmaxl,4)*cosnvn(k,n)
           fzcs(j2:jmax,ni,mj) = fzcs(j2:jmax,ni,mj) + work1(j2l:jmaxl,5)*sinnv(k,n)  &
                                                     + work1(j2l:jmaxl,6)*cosnvn(k,n)

           flsc(jl:ns,ni,mj) = flsc(jl:ns,ni,mj) + work1(jll:nsl,12)*sinnvn(k,n)
           flcs(jl:ns,ni,mj) = flcs(jl:ns,ni,mj) + work1(jll:nsl,9)*sinnv(k,n)   &
                                                 + work1(jll:nsl,10)*cosnvn(k,n)
        END DO
     END DO
  END DO

  DEALLOCATE (work1, tempr, tempz)

  if (dump_tomnsps .and. iter2.le.2) then
      write(dump_filename, 999) ns, iter2, trim(input_extension)
999 format('tomnsps_',i5.5,'_',i6.6,'.',a)
      open(unit=42, file=trim(dump_filename), status="unknown")

      write(42, *) "# ns ntor mpol1"
      write(42, *) ns, ntor, mpol1

      if (lthreed) then
        write(42, *) "# js n m frcc frss fzsc fzcs flsc flcs"
        DO js = 1, ns
          do n=0, ntor
            ni = n+ioff
            do m=0, mpol1
              mj = m+joff
              write(42, *) js, n, m, &
                frcc(js,ni,mj), frss(js,ni,mj), &
                fzsc(js,ni,mj), fzcs(js,ni,mj), &
                flsc(js,ni,mj), flcs(js,ni,mj)
            end do
          end do
        end do
      else ! lthreed
        write(42, *) "# js n m frcc fzsc flsc"
        DO js = 1, ns
          do n=0, ntor
            ni = n+ioff
            do m=0, mpol1
              mj = m+joff
              write(42, *) js, n, m, &
                frcc(js,ni,mj), fzsc(js,ni,mj), flsc(js,ni,mj)
            end do
          end do
        end do

      end if ! lthreed

      close(42)

      print *, "dumped tomnsps output to '"//trim(dump_filename)//"'"
!       stop
  end if

END SUBROUTINE tomnsps

!> \brief Fourier-transform anti-symmetric forces from real space to Fourier space
!>
!> @param frzl_array
!> @param armn
!> @param brmn
!> @param crmn
!> @param azmn
!> @param bzmn
!> @param czmn
!> @param blmn
!> @param clmn
!> @param arcon
!> @param azcon
SUBROUTINE tomnspa(frzl_array,       &
                   armn, brmn, crmn, &
                   azmn, bzmn, czmn, &
                         blmn, clmn, &
                   arcon, azcon       )
  USE vmec_main
  USE vmec_params, ONLY: jlam, jmin2, ntmax, rsc, rcs, zcc, zss
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax), TARGET, INTENT(inout) :: frzl_array
  REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(in) ::     &
     armn, brmn, crmn, azmn, bzmn, czmn, blmn, clmn, arcon, azcon

  INTEGER :: jmax, m, mparity, i, n, k, l, nsz
  INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl, js
  REAL(rprec), DIMENSION(:,:,:), POINTER :: frcs, frsc, fzcc, fzss, flcc, flss
  REAL(rprec), DIMENSION(:), ALLOCATABLE   :: temp1, temp3
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1

  character(len=255) :: dump_filename
  logical            :: dump_tomnspa = .true.

  frsc => frzl_array(:,:,:,rsc)               !!R-SIN(mu) COS(nv)
  fzcc => frzl_array(:,:,:,zcc+ntmax)         !!Z-COS(mu) COS(nv)
  flcc => frzl_array(:,:,:,zcc+2*ntmax)       !!L-COS(mu) COS(nv)
  IF (lthreed) THEN
     frcs => frzl_array(:,:,:,rcs)            !!R-COS(mu) SIN(nv)
     fzss => frzl_array(:,:,:,zss+ntmax)      !!Z-SIN(mu) SIN(nv)
     flss => frzl_array(:,:,:,zss+2*ntmax)    !!L-SIN(mu) SIN(nv)
  END IF

  nsz = ns*nzeta
  ALLOCATE (work1(nsz,12), temp1(nsz), temp3(nsz), stat=i)
  IF (i .ne. 0) STOP 'Allocation error in VMEC tomnspa'

  ioff = LBOUND(frsc,2)
  joff = LBOUND(frsc,3)

  jmax = ns
  IF (ivac .lt. 1) jmax = ns1

  ! BEGIN FOURIER TRANSFORM
  DO m = 0, mpol1
     mparity = MOD(m,2)

     mj = m+joff

     j2 = jmin2(m)
     jl = jlam(m)

     work1 = 0

     ! DO THETA (U) TRANSFORM FIRST
     DO i = 1, ntheta2
        temp1(:) = armn(:,i,mparity) &
! #ifndef _HBANGLE
            + xmpq(m,1)*arcon(:,i,mparity)
! #end /* ndef _HBANGLE */

        temp3(:) = azmn(:,i,mparity) &
! #ifndef _HBANGLE
            + xmpq(m,1)*azcon(:,i,mparity)
! #end /* ndef _HBANGLE */

        work1(:,3) = work1(:,3) + temp1(:)*sinmui(i,m) + brmn(:,i,mparity)*cosmumi(i,m)
        work1(:,5) = work1(:,5) + temp3(:)*cosmui(i,m) + bzmn(:,i,mparity)*sinmumi(i,m)
        work1(:,9) = work1(:,9) + blmn(:,i,mparity)*sinmumi(i,m)

        IF (.not.lthreed) CYCLE

        work1(:,1) = work1(:,1) + temp1(:)*cosmui(i,m) + brmn(:,i,mparity)*sinmumi(i,m)
        work1(:,2) = work1(:,2) - crmn(:,i,mparity)*cosmui(i,m)
        work1(:,4) = work1(:,4) - crmn(:,i,mparity)*sinmui(i,m)
        work1(:,6) = work1(:,6) - czmn(:,i,mparity)*cosmui(i,m)
        work1(:,7) = work1(:,7) + temp3(:)*sinmui(i,m) + bzmn(:,i,mparity)*cosmumi(i,m)
        work1(:,8) = work1(:,8) - czmn(:,i,mparity)*sinmui(i,m)
        work1(:,10) = work1(:,10) - clmn(:,i,mparity)*cosmui(i,m)
        work1(:,11) = work1(:,11) + blmn(:,i,mparity)*cosmumi(i,m)
        work1(:,12) = work1(:,12) - clmn(:,i,mparity)*sinmui(i,m)
     END DO

     ! NEXT, DO ZETA (V) TRANSFORM
     DO n = 0, ntor
        ni = n+ioff
        DO k = 1, nzeta
           l = ns*(k-1) ! current slice offset

           j2l   = j2+l   ! start of radial slice for R,Z
           jmaxl = jmax+l ! end of radial slice for R,Z

           jll = jl+l ! start of radial slice for lambda
           nsl = ns+l ! end of radial slice for lambda

           frsc(j2:jmax,ni,mj) = frsc(j2:jmax,ni,mj) + work1(j2l:jmaxl,3)*cosnv(k,n)
           fzcc(j2:jmax,ni,mj) = fzcc(j2:jmax,ni,mj) + work1(j2l:jmaxl,5)*cosnv(k,n)
           flcc(jl:ns,ni,mj) = flcc(jl:ns,ni,mj) + work1(jll:nsl,9)*cosnv(k,n)

           IF (.not.lthreed) CYCLE

           frsc(j2:jmax,ni,mj) = frsc(j2:jmax,ni,mj) + work1(j2l:jmaxl,4)*sinnvn(k,n)
           fzcc(j2:jmax,ni,mj) = fzcc(j2:jmax,ni,mj) + work1(j2l:jmaxl,6)*sinnvn(k,n)
           frcs(j2:jmax,ni,mj) = frcs(j2:jmax,ni,mj) + work1(j2l:jmaxl,1)*sinnv(k,n)  &
                                                     + work1(j2l:jmaxl,2)*cosnvn(k,n)
           fzss(j2:jmax,ni,mj) = fzss(j2:jmax,ni,mj) + work1(j2l:jmaxl,7)*sinnv(k,n)  &
                                                     + work1(j2l:jmaxl,8)*cosnvn(k,n)
           flcc(jl:ns,ni,mj) = flcc(jl:ns,ni,mj) + work1(jll:nsl,10)*sinnvn(k,n)
           flss(jl:ns,ni,mj) = flss(jl:ns,ni,mj) + work1(jll:nsl,11)*sinnv(k,n)  &
                                                 + work1(jll:nsl,12)*cosnvn(k,n)
        END DO
     END DO
  END DO

  DEALLOCATE (work1, temp1, temp3)

  if (dump_tomnspa .and. iter2.le.2) then
      write(dump_filename, 998) ns, iter2, trim(input_extension)
998 format('tomnspa_',i5.5,'_',i6.6,'.',a)
      open(unit=42, file=trim(dump_filename), status="unknown")

      write(42, *) "# ns ntor mpol1"
      write(42, *) ns, ntor, mpol1

      if (lthreed) then
        write(42, *) "# js n m frsc frcs fzcc fzss flcc flss"
        DO js = 1, ns
          do n=0, ntor
            ni = n+ioff
            do m=0, mpol1
              mj = m+joff
              write(42, *) js, n, m, &
                frsc(js,ni,mj), frcs(js,ni,mj), &
                fzcc(js,ni,mj), fzss(js,ni,mj), &
                flcc(js,ni,mj), flss(js,ni,mj)
            end do
          end do
        end do
      else ! lthreed
        write(42, *) "# js n m frsc fzcc flcc"
        DO js = 1, ns
          do n=0, ntor
            ni = n+ioff
            do m=0, mpol1
              mj = m+joff
              write(42, *) js, n, m, &
                frsc(js,ni,mj), fzcc(js,ni,mj), flcc(js,ni,mj)
            end do
          end do
        end do

      end if ! lthreed

      close(42)

      print *, "dumped tomnspa output to '"//trim(dump_filename)//"'"
!       stop
  end if

END SUBROUTINE tomnspa
