!> \file
SUBROUTINE fouri(grpmn, gsource, amatrix, amatsq, bvec, wint, ndim, ns)
  USE vacmod
  IMPLICIT NONE

  INTEGER, INTENT(in) :: ns, ndim
  REAL(rprec), DIMENSION(mnpd,nv,nu3,ndim), INTENT(in) :: grpmn
  REAL(rprec), DIMENSION(nuv), INTENT(in) :: gsource
  REAL(rprec), DIMENSION(mnpd,mnpd,ndim**2), INTENT(out) :: amatrix
  REAL(rprec), DIMENSION(mnpd2,mnpd2), INTENT(out) :: amatsq
  REAL(rprec), DIMENSION(0:mf,-nf:nf,ndim), INTENT(inout) :: bvec
  REAL(rprec), DIMENSION(*) :: wint

  !> interior  (int_ext=-1), exterior  (int_ext=+1)  neumann problem
  REAL(rprec), PARAMETER :: int_ext = 1

  INTEGER :: k, i, j, n, kvi, kui, mn, m, mn0, isym
  REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: bcos, bsin, source
  REAL(rprec), ALLOCATABLE :: actemp(:,:,:,:), astemp(:,:,:,:)
  REAL(rprec) :: cosn, sinn, cosm, sinm

  ! AMATRIX(,1) = A(Sin)(Sin');  AMATRIX(,2) = A(Sin)(Cos');
  ! AMATRIX(,3) = A(Cos)(Sin');  AMATRIX(,4) = A(Cos)(Cos')
  !
  ! ARG OF TRIG FUNCTIONS: m'u' - n'v' CORRESPONDING TO THE PRIMED MESH IN
  ! PKM EQ.(2.14), AND mu  - nv (UNPRIMED MESH) IN EQ. (2.16)
  !
  ! ON ENTRY, GRPMN(MN,...) HAS BEEN FOURIER-TRANSFORMED OVER THE UNPRIMED
  ! COORDINATES (IN FOURP), SO THE FIRST INDEX OF GRPMN CORRESPONDS TO THE FIRST
  ! INDEX OF AMATRIX. THUS, THE FOURIER TRANSFORMS OF GRPMN HERE ARE OVER THE PRIMED MESH.
  !
  ! IN CONTRAST, THE INTEGRAL OF THE SOURCE TERM OVER THE PRIMED MESH WAS ALREADY
  ! DONE (IN SCALPOT), SO HERE THE FT ARE OVER THE UNPRIMED MESH FOR THE SOURCE.

  ALLOCATE (bcos(nu2,-nf:nf,ndim), bsin(nu2,-nf:nf,ndim),           &
     actemp(mnpd,-nf:nf,nu3,ndim), astemp(mnpd,-nf:nf,nu3,ndim),    &
     source(nv,nu2,ndim), stat = i)
  IF (i .ne. 0) STOP 'allocation error in fouri'

  ! SYMMETRIZE SOURCE TERMS (with respect to u,v and -u,-v)
  ! INDEX (1) IS ANTI-SYMMETRIC, INDEX (2) IS SYMMETRIC, PART
  !
  ! GSOURCE = - (2pi)**2 B * dS (h - hsing) * WINT
  !
  ! WINT: needed to normalize integral over angles to unity
  k = 0
  DO i = 1, nu2
     DO j = 1, nv
        k = k + 1

        ! anti-symmetric part
        source(j,i,1) = gsource(k) - gsource(imirr(k))

        IF (lasym) then
           ! symmetric part
           source(j,i,2) = gsource(k) + gsource(imirr(k))
        end if
     END DO
  END DO

  source = p5*onp*source

  ! INITIALIZE RUNNING-SUM ARRAYS
  bcos = 0;      bsin = 0
  actemp = 0;    astemp = 0
  amatrix = 0

  ! PERFORM KV (TOROIDAL ANGLE) TRANSFORM
  DO n = 0, nf
     DO kvi = 1, nv
        cosn = cosv(n,kvi)
        sinn = sinv(n,kvi)
        DO isym = 1, ndim
           bcos(:,n,isym) = bcos(:,n,isym) + cosn*source(kvi,:,isym)
           bsin(:,n,isym) = bsin(:,n,isym) + sinn*source(kvi,:,isym)
           actemp(:,n,:,isym) = actemp(:,n,:,isym) + cosn*grpmn(:,kvi,:,isym)
           astemp(:,n,:,isym) = astemp(:,n,:,isym) + sinn*grpmn(:,kvi,:,isym)

           IF (n .ne. 0) THEN
              bcos(:,(-n),isym) =  bcos(:,n,isym)
              bsin(:,(-n),isym) = -bsin(:,n,isym)
              actemp(:,(-n),:,isym) =  actemp(:,n,:,isym)
              astemp(:,(-n),:,isym) = -astemp(:,n,:,isym)
           ENDIF
        END DO
     END DO
  END DO

  ! PERFORM KU (POLOIDAL ANGLE) TRANSFORM
  DO m = 0, mf
     DO kui = 1, nu2
        cosm = cosui(m,kui)
        sinm = sinui(m,kui)
        bvec(m,-nf:nf,1)    = bvec(m,-nf:nf,1) + bcos(kui,-nf:nf,1)*sinm - bsin(kui,-nf:nf,1)*cosm
        IF (lasym) THEN
           bvec(m,-nf:nf,2) = bvec(m,-nf:nf,2) + bcos(kui,-nf:nf,2)*cosm + bsin(kui,-nf:nf,2)*sinm
        END IF
     END DO

    ! RECALL, LAST INDEX OF AS,CTEMP
    !                    = 1 CORRESPONDS TO SIN (UNPRIMED) TRANSFORM (FIRST INDEX OF AMATRIX)
    !                    = 2 CORRESPONDS TO COS (UNPRIMED) TRANSFORM
     DO kui = 1, nu3
        cosm = cosu(m,kui)*wint(kui*ns*nv)
        sinm = sinu(m,kui)*wint(kui*ns*nv)

        ! SIN SIN'
        amatrix(:,m+1:mnpd:mf1,1) = amatrix(:,m+1:mnpd:mf1,1) + sinm*actemp(:,-nf:nf,kui,1) - cosm*astemp(:,-nf:nf,kui,1)
        IF (.not.lasym) CYCLE

        ! SIN COS'
        amatrix(:,m+1:mnpd:mf1,2) = amatrix(:,m+1:mnpd:mf1,2) + cosm*actemp(:,-nf:nf,kui,1) + sinm*astemp(:,-nf:nf,kui,1)

        ! COS SIN'
        amatrix(:,m+1:mnpd:mf1,3) = amatrix(:,m+1:mnpd:mf1,3) + sinm*actemp(:,-nf:nf,kui,2) - cosm*astemp(:,-nf:nf,kui,2)

        ! COS COS'
        amatrix(:,m+1:mnpd:mf1,4) = amatrix(:,m+1:mnpd:mf1,4) + cosm*actemp(:,-nf:nf,kui,2) + sinm*astemp(:,-nf:nf,kui,2)
     END DO
  END DO

  DEALLOCATE (bcos, bsin, actemp, astemp, source, stat=i)

  amatrix = (pi2*pi2)*amatrix

  ! ZERO BVEC(0,n) FOR n < 0
  ! Fixed SPH081515: had -nf:0 before
  bvec(0,-nf:-1,1:ndim) = 0

  !     ZERO AMATRIX(0,n,m',n') M = 0 MODES FOR n < 0
  ! Index of m=0,n=0
  mn0 = 1+mf1*nf
  ! SPH082415: mn0-mf1: (m=0,n=-1 index)
  amatrix(1:mn0-mf1:mf1,:,1:ndim*ndim) = 0

  ! ADD DIAGONAL TERMS TO AMATRIX [THE FIRST TERM IN EQ(3.2) OF PKM]
  DO mn = 1, mnpd
     amatrix(mn,mn,1) = amatrix(mn,mn,1) + pi3*int_ext
  END DO

  IF (lasym) THEN
     DO mn = 1, mnpd
        amatrix(mn,mn,4) = amatrix(mn,mn,4) + pi3*int_ext
     END DO

     ! COS(0-0) mode *2
     amatrix(mn0,mn0,4) = amatrix(mn0,mn0,4) + pi3*int_ext
  END IF

  ! PUT ELEMENTS INTO SQUARE MATRIX
  ! Sin-Sin'
  amatsq(:mnpd,:mnpd) = amatrix(:,:,1)

  IF (.not.lasym) RETURN

  ! Sin-Cos'
  amatsq(:mnpd,1+mnpd:mnpd2) = amatrix(:,:,2)

  ! Cos-Sin'
  amatsq(1+mnpd:mnpd2,:mnpd) = amatrix(:,:,3)

  ! Cos-Cos'
  amatsq(1+mnpd:mnpd2,1+mnpd:mnpd2) = amatrix(:,:,4)

END SUBROUTINE fouri