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
SUBROUTINE tomnsps_con(frzl_array,       &
                       brmn_con, bzmn_con, &
                       arcon, azcon       )
  USE vmec_main
  USE vmec_params, ONLY: jlam, jmin2, ntmax, rcc, rss, zsc, zcs
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax), TARGET, INTENT(out) :: frzl_array
  REAL(rprec), DIMENSION(ns*nzeta*ntheta3,0:1), INTENT(in) :: &
    brmn_con, bzmn_con, arcon, azcon

  INTEGER :: jmax, m, mparity, i, n, k, l, nsz
  INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl, js
  REAL(rprec), DIMENSION(:,:,:), POINTER :: frcc, frss, fzcs, fzsc, flcs, flsc
  REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: work1
  REAL(rprec), DIMENSION(:), ALLOCATABLE   :: tempr, tempz

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

  frzl_array = 0.0_dp

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

     work1 = 0.0_dp

     ! DO THETA (U) INTEGRATION FIRST ON HALF INTERVAL (0 < U < PI)
     l = 0
     DO i = 1, ntheta2
        jll = l+1   ! start of poloidal slice
        nsl = l+nsz ! end of poloidal slice
        l   = l+nsz ! jump to next poloidal slice

        tempr(:) = xmpq(m,1)*arcon(jll:nsl,mparity)
        tempz(:) = xmpq(m,1)*azcon(jll:nsl,mparity)

        work1(:,1) = work1(:,1) + tempr(:)*cosmui(i,m) + brmn_con(jll:nsl,mparity)*sinmumi(i,m)
        work1(:,7) = work1(:,7) + tempz(:)*sinmui(i,m) + bzmn_con(jll:nsl,mparity)*cosmumi(i,m)

        IF (.not.lthreed) CYCLE

        work1(:,3) = work1(:,3) + tempr(:)*sinmui(i,m) + brmn_con(jll:nsl,mparity)*cosmumi(i,m)
        work1(:,5) = work1(:,5) + tempz(:)*cosmui(i,m) + bzmn_con(jll:nsl,mparity)*sinmumi(i,m)

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

           IF (.not.lthreed) CYCLE

           frss(j2:jmax,ni,mj) = frss(j2:jmax,ni,mj) + work1(j2l:jmaxl,3)*sinnv(k,n)
           fzcs(j2:jmax,ni,mj) = fzcs(j2:jmax,ni,mj) + work1(j2l:jmaxl,5)*sinnv(k,n)

        END DO
     END DO
  END DO

  DEALLOCATE (work1, tempr, tempz)

END SUBROUTINE tomnsps_con

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
SUBROUTINE tomnspa_con(frzl_array,       &
                       brmn_con, bzmn_con, &
                       arcon, azcon       )
  USE vmec_main
  USE vmec_params, ONLY: jlam, jmin2, ntmax, rsc, rcs, zcc, zss
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax), TARGET, INTENT(inout) :: frzl_array
  REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(in) :: &
    brmn_con, bzmn_con, arcon, azcon

  INTEGER :: jmax, m, mparity, i, n, k, l, nsz
  INTEGER :: ioff, joff, mj, ni, nsl, j2, j2l, jl, jll, jmaxl, js
  REAL(rprec), DIMENSION(:,:,:), POINTER :: frcs, frsc, fzcc, fzss, flcc, flss
  REAL(rprec), DIMENSION(:), ALLOCATABLE   :: temp1, temp3
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1

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

     work1 = 0.0_dp

     ! DO THETA (U) TRANSFORM FIRST
     DO i = 1, ntheta2
        temp1(:) = xmpq(m,1)*arcon(:,i,mparity)
        temp3(:) = xmpq(m,1)*azcon(:,i,mparity)

        work1(:,3) = work1(:,3) + temp1(:)*sinmui(i,m) + brmn_con(:,i,mparity)*cosmumi(i,m)
        work1(:,5) = work1(:,5) + temp3(:)*cosmui(i,m) + bzmn_con(:,i,mparity)*sinmumi(i,m)

        IF (.not.lthreed) CYCLE

        work1(:,1) = work1(:,1) + temp1(:)*cosmui(i,m) + brmn_con(:,i,mparity)*sinmumi(i,m)
        work1(:,7) = work1(:,7) + temp3(:)*sinmui(i,m) + bzmn_con(:,i,mparity)*cosmumi(i,m)

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

           IF (.not.lthreed) CYCLE

           frcs(j2:jmax,ni,mj) = frcs(j2:jmax,ni,mj) + work1(j2l:jmaxl,1)*sinnv(k,n)
           fzss(j2:jmax,ni,mj) = fzss(j2:jmax,ni,mj) + work1(j2l:jmaxl,7)*sinnv(k,n)

        END DO
     END DO
  END DO

  DEALLOCATE (work1, temp1, temp3)

END SUBROUTINE tomnspa_con
