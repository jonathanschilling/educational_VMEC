!> \file
SUBROUTINE analysum2(grpmn, bvec, slp, tlp, slm, tlm, m, n, l, ivacskip, ndim)
  USE vacmod
  IMPLICIT NONE

  INTEGER, INTENT(in) :: m, n, l, ivacskip, ndim
  REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,nuv2,ndim)
  REAL(rprec), INTENT(inout) :: bvec(0:mf,-nf:nf,ndim)
  REAL(rprec), DIMENSION(nuv2), INTENT(in) :: slp, tlp, slm, tlm

  INTEGER :: i
  REAL(rprec) :: sinp, sinm, cosp, cosm, temp

  IF (n .LT. 0) STOP 'error calling analysum2!'

  DO i = 1,nuv2
     sinp =  sinu1(i,m)*cosv1(i,n) * cmns(l,m,n)
     temp = -cosu1(i,m)*sinv1(i,n) * cmns(l,m,n)

     ! SIN(mu + |n|v) * cmns (l,m,|n|)
     sinm = sinp - temp

     ! SIN(mu - |n|v) * cmns (l,m,|n|)
     sinp = sinp + temp

     bvec(m, n,1) = bvec(m, n,1) + tlp(i)*bexni(i)*sinp
     bvec(m,-n,1) = bvec(m,-n,1) + tlm(i)*bexni(i)*sinm

     IF (ivacskip .EQ. 0) THEN
        grpmn(m, n,i,1) = grpmn(m, n,i,1) + slp(i)*sinp
        grpmn(m,-n,i,1) = grpmn(m,-n,i,1) + slm(i)*sinm
     END IF

     IF (lasym) THEN
        cosp = cosu1(i,m)*cosv1(i,n) * cmns(l,m,n)
        temp = sinu1(i,m)*sinv1(i,n) * cmns(l,m,n)

        ! COS(mu + |n|v) * cmns (l,m,|n|)
        cosm = cosp - temp

        ! COS(mu - |n|v) * cmns (l,m,|n|)
        cosp = cosp + temp

        bvec(m, n,2) = bvec(m, n,2) + tlp(i)*bexni(i)*cosp
        bvec(m,-n,2) = bvec(m,-n,2) + tlm(i)*bexni(i)*cosm

        IF (ivacskip .EQ. 0) THEN
           grpmn(m, n,i,2) = grpmn(m, n,i,2) + slp(i)*cosp
           grpmn(m,-n,i,2) = grpmn(m,-n,i,2) + slm(i)*cosm
        END IF
     END IF
  END DO

END SUBROUTINE analysum2
