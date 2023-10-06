!> \file
!> \brief Fourier-space bandpass filter on constraint force for spectral condensation

!> \brief Fourier-space bandpass filter on constraint force for spectral condensation
!>
!> Eliminates contributions to gcon from m=1 and m=mpol1-1=mpol-2 poloidal mode numbers.
!>
!> @param gcons ! output: Fourier-filtered spectral condensation force gcon
!> @param ztemp ! input: full-m gcon in real-space constraint force
!> @param gcs ! temporary storage
!> @param gsc ! temporary storage
!> @param gcc ! temporary storage
!> @param gss ! temporary storage
SUBROUTINE alias(gcons, ztemp, gcs, gsc, gcc, gss)
  USE vmec_main
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns*nzeta,ntheta3), INTENT(out) :: gcons
  REAL(rprec), DIMENSION(ns*nzeta,ntheta3), INTENT(in)  :: ztemp
  REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1), intent(inout) :: gcs, gsc, gcc, gss

  REAL(rprec), PARAMETER :: p5 = 0.5_dp

  INTEGER :: m, i, ir, jk, jka, n, k, js, l
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work, gcona

  ! Fourier transform alias force from ztemp to gcons
  ! and also return intermediate output in g(c,s)(c,s)

  ALLOCATE (work(ns*nzeta,4), gcona(ns*nzeta,ntheta3))

  gcons = 0.0_dp
  gcona = 0.0_dp

  gcs = 0.0_dp;  gsc = 0.0_dp
  gcc = 0.0_dp;  gss = 0.0_dp

  ! The start of this loop at m=1 and its end at mpol1-1=mpol-2
  ! is what makes this routine a Fourier-space bandpass filter.
  DO m = 1, mpol1-1
     work = 0
     DO i = 1, ntheta2
        DO jk = 1, ns*nzeta
           work(jk,1) = work(jk,1) + ztemp(jk,i)*cosmui(i,m)
           work(jk,2) = work(jk,2) + ztemp(jk,i)*sinmui(i,m)
        END DO
        IF (lasym) then
           ir = ntheta1 + 2 - i
           IF (i .eq. 1) ir = 1
           DO jk = 1, ns*nzeta
              jka = ireflect(jk)
              work(jk,3) = work(jk,3) + ztemp(jka,ir)*cosmui(i,m)
              work(jk,4) = work(jk,4) + ztemp(jka,ir)*sinmui(i,m)
           END DO
        end if
     END DO

     DO n = 0, ntor ! retain full toroidal resolution
        DO k = 1, nzeta
           l = ns*(k-1)
           IF (.not.lasym) THEN
              DO js = 2,ns
                 gcs(js,n,m) = gcs(js,n,m) + tcon(js)*work(js+l,1)*sinnv(k,n)
                 gsc(js,n,m) = gsc(js,n,m) + tcon(js)*work(js+l,2)*cosnv(k,n)
              END DO
           ELSE
              DO js = 2,ns
                 gcs(js,n,m) = gcs(js,n,m) + p5*tcon(js)*sinnv(k,n)*(work(js+l,1)-work(js+l,3))
                 gsc(js,n,m) = gsc(js,n,m) + p5*tcon(js)*cosnv(k,n)*(work(js+l,2)-work(js+l,4))
                 gss(js,n,m) = gss(js,n,m) + p5*tcon(js)*sinnv(k,n)*(work(js+l,2)+work(js+l,4))
                 gcc(js,n,m) = gcc(js,n,m) + p5*tcon(js)*cosnv(k,n)*(work(js+l,1)+work(js+l,3))
              END DO
           END IF
        END DO
     END DO

     ! INVERSE FOURIER TRANSFORM DE-ALIASED GCON
     work = 0.0_dp

     DO n = 0, ntor
        DO k = 1, nzeta
           l = ns*(k-1)
           DO js = 2, ns
              work(js+l,3) = work(js+l,3) + gcs(js,n,m)*sinnv(k,n)
              work(js+l,4) = work(js+l,4) + gsc(js,n,m)*cosnv(k,n)
           END DO
           IF (lasym) then
              DO js = 2, ns
                 work(js+l,1) = work(js+l,1) + gcc(js,n,m)*cosnv(k,n)
                 work(js+l,2) = work(js+l,2) + gss(js,n,m)*sinnv(k,n)
              END DO
           end if
        END DO
     END DO

     DO i = 1, ntheta2
        DO jk = 1, ns*nzeta
           gcons(jk,i) = gcons(jk,i) + (work(jk,3)*cosmu(i,m) + work(jk,4)*sinmu(i,m))*faccon(m)
        END DO
        IF (lasym) then
           DO jk = 1, ns*nzeta
              gcona(jk,i) = gcona(jk,i) + (work(jk,1)*cosmu(i,m) + work(jk,2)*sinmu(i,m))*faccon(m)
           END DO
        end if
     END DO

  END DO

  IF (lasym) THEN

     ! EXTEND GCON INTO THETA = PI,2*PI DOMAIN
     DO i = 1 + ntheta2, ntheta1
        ir = ntheta1 + 2 - i
        DO jk = 1, ns*nzeta
           jka = ireflect(jk)
           gcons(jk,i) = -gcons(jka,ir) + gcona(jka,ir)
        END DO
     END DO

     ! ADD SYMMETRIC, ANTI-SYMMETRIC PIECES IN THETA = 0,PI DOMAIN
     gcons(:,:ntheta2) = gcons(:,:ntheta2) + gcona(:,:ntheta2)

  END IF

  DEALLOCATE (work, gcona)

END SUBROUTINE alias
