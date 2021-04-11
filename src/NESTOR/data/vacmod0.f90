!> \file
MODULE vacmod0
  IMPLICIT NONE

  INTEGER :: mf
  INTEGER :: nf
  INTEGER :: nu
  INTEGER :: nv
  INTEGER :: mf1
  INTEGER :: nf1
  INTEGER :: mnpd
  INTEGER :: mnpd2
  INTEGER :: nuv
  INTEGER :: nu2
  INTEGER :: nu3
  INTEGER :: nuv2
  INTEGER :: nfper

  INTEGER :: nvper
  INTEGER :: nuv_tan
  INTEGER :: nvp

  INTEGER :: ndim

  contains

subroutine set_nestor_sizes(nfp, ntor, mpol, nzeta, ntheta, lasym)

  integer, intent(in) :: nfp
  integer, intent(in) :: ntor
  integer, intent(in) :: mpol
  integer, intent(in) :: nzeta
  integer, intent(in) :: ntheta
  logical, intent(in) :: lasym

  ! copied from vmec:read_indata
  integer :: ntheta1, ntheta2, ntheta3, nznt

  print *, "set_nestor_sizes"

  ! even (rounded down) ntheta
  ntheta1 = 2*(ntheta/2)

  ! odd stellarator-symmetric little-more-than-half of ntheta
  ntheta2 = 1 + ntheta1/2

  IF (.NOT.lasym) THEN
     ntheta3 = ntheta2
  ELSE
     ntheta3 = ntheta1
  END IF

  nznt = nzeta*ntheta3

  ! from here on, for NESTOR only
  mf = mpol+1
  nf = ntor
  nu = ntheta1
  nv = nzeta
  mf1 = 1+mf
  nf1 = 2*nf+1
  mnpd = mf1*nf1

  IF (.NOT.lasym) THEN
    mnpd2 = mnpd
    ndim = 1
  ELSE
    mnpd2 = 2*mnpd
    ndim = 2
  END IF
  ! ndim = mnpd2/mnpd

  nuv = nu*nv
  nfper = nfp
  nu2 = nu/2 + 1 ! nu2 is equal to ntheta2
  nu3 = ntheta3
  nuv2 = nznt

  IF (nv == 1) THEN
     ! AXISYMMETRIC CASE: DO FP SUM TO INTEGRATE IN V
     nvper = 64
     nuv_tan = 2*nu*nvper
  ELSE
     nvper = nfper
     nuv_tan = 2*nuv
  END IF

  nvp = nv*nvper



  ! in the end:
  ! nu  = ntheta1
  ! nu2 = ntheta2
  ! nu3 = ntheta3

end subroutine set_nestor_sizes

END MODULE vacmod0
