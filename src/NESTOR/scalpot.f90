!> \file
SUBROUTINE scalpot(bvec, amatrix, wint, ns, ivacskip, lasym)
  USE vacmod
  IMPLICIT NONE

  INTEGER, INTENT(in) :: ns, ivacskip
  REAL(rprec), INTENT(out) :: bvec(mnpd2), amatrix(mnpd2*mnpd2)
  REAL(rprec), INTENT(in) :: wint(*)
  logical, intent(in) :: lasym

  INTEGER :: ip, ndim, istat
  REAL(rprec), ALLOCATABLE :: grpmn(:), green(:), gstore(:)
  REAL(rprec), ALLOCATABLE :: greenp(:,:)

  IF (.not.ALLOCATED(amatsav)) then
     STOP 'AMATSAV not allocated in scalpot'
  end if

  ALLOCATE (grpmn(nuv2*mnpd2), stat=ip)
  IF (ip .ne. 0) STOP 'GRPMN: Allocation error in scalpot'

  ! COMPUTE TRANFORM OF ANALYTIC SOURCE AND KERNEL
  ! ON EXIT, BVEC CONTAINS THE TRANSFORM OF THE ANALYTIC SOURCE
  ! AND GRPMN CONTAINS TRANSFORM OF NORMAL DERIVATIVE
  ! OF THE GREENS FUNCTION [PKM, EQ.(2.15)]
  !
  ! FOR ivacskip != 0, USE PREVIOUSLY COMPUTED bvecsav FOR SPEED

  ndim = mnpd2/mnpd
  CALL analyt (grpmn, bvec, ivacskip, ndim, lasym)

  IF (ivacskip .ne. 0) THEN
     ! Here, bvecsav contains the previous non-singular contribution to the "final" bvec from fouri.
     ! For ivacskip != 0, this contribution is used from the cache in bvecsav.
     bvec = bvec + bvecsav
  ELSE
     ! Here, bvecsav stores the singular part of bvec from analyt
     ! to be subtracted from the "final" bvec computed below by fouri
     bvecsav = bvec


     ALLOCATE (green(nuv), gstore(nuv), greenp(nuv,nuv2), stat=istat)
     if (istat.ne.0) then
        ! Below loop over nuv2 was previously chunked.
        ! Therefore, some extra care shall be used here to make sure
        ! everything still fits into memory...
        stop 'green allocation error in scalpot.'
     end if


     ! COMPUTE SURFACE INTEGRALS OF SOURCE AND GREENS FUNCTION NEEDED
     ! FOR SPECTRAL DECOMPOSITION OF POTENTIAL INTEGRAL EQUATION
     ! NOTE: SOURCE IS THE RHS OF EQ.(3.2), KERNEL IS THE LHS OF EQ (3.2).
     ! IP IS THE INDEX OF THE PRIMED VARIABLE MESH.

     gstore  = 0
     DO ip = 1, nuv2

        ! COMPUTE DIFFERENCE BETWEEN THE EXACT AND ANALYTIC GREENS FUNCTION AND GRADIENT
        ! [FIRST TERMS IN EQ.(2.14, 2.16)].
        CALL greenf (green, greenp(1,ip), ip)

        ! PERFORM INTEGRAL (SUM) OVER PRIMED MESH OF NON-SINGULAR SOURCE TERM
        ! [(h-hsing)(u,v,u',v') == bexni(ip)*green(u,v; ip) in Eq. 2.16]
        ! AND STORE IT - FOR UNPRIMED MESH VALUES - IN GSTORE
        gstore = gstore + bexni(ip)*green

     END DO

     ! PERFORM FOURIER INTEGRAL OF GRADIENT KERNEL (GREENP) OVER THE UNPRIMED MESH
     ! AND STORE IN GRPMN (NOTE THAT GRPMN IS ADDED TO THE ANALYTIC PIECE IN EQ. 2.14,
     ! - COMPUTED IN ANALYT - WHICH HAS THE APPROPRIATE SIN, COS FACTORS ALREADY)
     CALL fourp (grpmn, greenp, ndim)

     ! COMPUTE FOURIER INTEGRAL OF GRADIENT (GRPMN) OVER PRIMED MESH IN EQ. 2.14
     ! AND SOURCE (GSTORE) OVER UNPRIMED MESH IN EQ. 2.16
     CALL fouri (grpmn, gstore, amatrix, amatsav, bvec, wint, ndim, ns, lasym)

     DEALLOCATE (green, greenp, gstore)

     ! SAVE NON-SINGULAR CONTRIBUTION TO BVEC (IN BVECSAV)
     bvecsav(:mnpd2) = bvec - bvecsav(:mnpd2)

  ENDIF

  DEALLOCATE (grpmn)

  amatrix = amatsav

END SUBROUTINE scalpot
