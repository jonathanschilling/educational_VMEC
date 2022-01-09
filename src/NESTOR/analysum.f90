!> \file
!> \brief Compute the (m=0 or n=0) part of the DFT of the analytical Fourier transforms of the equivalently-singular integrals.

!> \brief Compute the (m=0 or n=0) part of the DFT of the analytical Fourier transforms of the equivalently-singular integrals.
!>
!> @param grpmn
!> @param bvec
!> @param sl
!> @param tl
!> @param m
!> @param n
!> @param l
!> @param ivacskip
!> @param lasym
!> @param m_map
!> @param n_map
!> @param grpmn_m_map
!> @param grpmn_n_map
SUBROUTINE analysum(grpmn, bvec, sl, tl, m, n, l, ivacskip, lasym, m_map, n_map, &
                    grpmn_m_map, grpmn_n_map)
  USE vacmod, vm_grpmn => grpmn
  IMPLICIT NONE

  INTEGER, INTENT(in) :: m, n, l, ivacskip
  REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,nuv2,ndim)
  REAL(rprec), INTENT(inout) :: bvec(0:mf,-nf:nf,ndim)
  real(rprec), intent(inout) :: m_map(0:mf,-nf:nf)
  real(rprec), intent(inout) :: n_map(0:mf,-nf:nf)
  real(rprec), intent(inout) :: grpmn_m_map(0:mf,-nf:nf,nuv2)
  real(rprec), intent(inout) :: grpmn_n_map(0:mf,-nf:nf,nuv2)
  REAL(rprec), DIMENSION(nuv2), INTENT(in) :: sl, tl
  logical, intent(in) :: lasym

  INTEGER :: i
  REAL(rprec) :: sinp, cosp

  IF (n .LT. 0) STOP 'error calling analysum!'

  m_map(m, n) =  m
  m_map(m,-n) =  m

  n_map(m, n) =  n
  n_map(m,-n) = -n

!   if (cmns(l,m,n) .eq. zero) then
!      ! no need to compute zeros....
!      return
!   end if

  DO i = 1, nuv2

     grpmn_m_map(m, n, i) =  m
     grpmn_m_map(m,-n, i) =  m

     grpmn_n_map(m, n, i) =  n
     grpmn_n_map(m,-n, i) = -n

     ! SIN(mu - |n|v)*cmns
     sinp = sinu1(i,m)*cosv1(i,n) - cosu1(i,m)*sinv1(i,n)

     ! Fourier-transform S_l or T_l*bexni and then
     ! add up Fourier coefficients, weighted by cmns
     bvec (m,n,1)   = bvec (m,n,1)    + cmns(l,m,n) * bexni(i)*tl(i)*sinp
     IF (ivacskip .EQ. 0) then
     grpmn(m,n,i,1) = grpmn(m,n,i,1)  + cmns(l,m,n) *          sl(i)*sinp
     end if

     IF (lasym) THEN
        ! COS(mu - |n|v)*cmns
        cosp = cosu1(i,m)*cosv1(i,n) + sinv1(i,n)*sinu1(i,m)

        bvec(m,n,2) = bvec(m,n,2) + tl(i)*cmns(l,m,n)*bexni(i)*cosp
        IF (ivacskip .EQ. 0) then
           grpmn(m,n,i,2) = grpmn(m,n,i,2)  + sl(i)*cmns(l,m,n)*cosp
        end if
     END IF
  END DO

END SUBROUTINE analysum
