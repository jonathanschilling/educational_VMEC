!> \file
SUBROUTINE analysum(grpmn, bvec, sl, tl, m, n, l, ivacskip, ndim)
  USE vacmod
  IMPLICIT NONE

  INTEGER, INTENT(in) :: m, n, l, ivacskip, ndim
  REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,nuv2,ndim)
  REAL(rprec), INTENT(inout) :: bvec(0:mf,-nf:nf,ndim)
  REAL(rprec), DIMENSION(nuv2), INTENT(in) :: sl, tl

  INTEGER :: i
  REAL(rprec) :: sinp, cosp

  IF (n .LT. 0) STOP 'error calling analysum!'

  DO i = 1, nuv2
     ! SIN(mu - |n|v)*cmns
     sinp = (sinu1(i,m)*cosv1(i,n) - sinv1(i,n)*cosu1(i,m)) * cmns(l,m,n)

     IF (ivacskip .EQ. 0) then
        grpmn(m,n,i,1) = grpmn(m,n,i,1)  + sl(i)*sinp
     end if
     bvec(m,n,1) = bvec(m,n,1) + tl(i)*bexni(i)*sinp

     IF (lasym) THEN
        ! COS(mu - |n|v)*cmns
        cosp = (cosu1(i,m)*cosv1(i,n) + sinv1(i,n)*sinu1(i,m)) * cmns(l,m,n)

        IF (ivacskip .EQ. 0) then
           grpmn(m,n,i,2) = grpmn(m,n,i,2)  + sl(i)*cosp
        end if
        bvec(m,n,2) = bvec(m,n,2) + tl(i)*bexni(i)*cosp
     END IF
  END DO

END SUBROUTINE analysum
