!> \file
SUBROUTINE analysum2(grpmn, bvec, m, n, l, ivacskip, lasym, m_map, n_map)
  USE vacmod, vm_grpmn => grpmn
  IMPLICIT NONE

  INTEGER, INTENT(in) :: m, n, l, ivacskip
  REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,nuv2,ndim)
  REAL(rprec), INTENT(inout) :: bvec(0:mf,-nf:nf,ndim)
  real(rprec), intent(inout) :: m_map(0:mf,-nf:nf)
  real(rprec), intent(inout) :: n_map(0:mf,-nf:nf)
  logical, intent(in) :: lasym

  INTEGER :: i
  REAL(rprec) :: sinp, sinm, cosp, cosm, temp

  IF (n .LT. 0) STOP 'error calling analysum2!'

  m_map(m,  n) =  m
  m_map(m, -n) =  m

  n_map(m,  n) =  n
  n_map(m, -n) = -n

  if (cmns(l,m,n) .eq. zero) then
     ! no need to compute zeros...
     return
  end if

  DO i = 1,nuv2

     sinp =  sinu1(i,m)*cosv1(i,n) * cmns(l,m,n)
     temp = -cosu1(i,m)*sinv1(i,n) * cmns(l,m,n)

     ! SIN(mu + |n|v) * cmns (l,m,|n|)
     sinm = sinp - temp

     ! SIN(mu - |n|v) * cmns (l,m,|n|)
     sinp = sinp + temp

     IF (ivacskip .EQ. 0) THEN
     grpmn(m, n,i,1) = grpmn(m, n,i,1) + slp(i)         *sinp
     grpmn(m,-n,i,1) = grpmn(m,-n,i,1) + slm(i)         *sinm
     END IF

     bvec (m, n,  1) = bvec (m, n,  1) + tlp(i)*bexni(i)*sinp
     bvec (m,-n,  1) = bvec (m,-n,  1) + tlm(i)*bexni(i)*sinm


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
