!> \file
MODULE vacmod
  USE vacmod0
  USE vac_persistent
  USE vparams, ONLY: zero, one, c2p0, cp5

  IMPLICIT NONE

  REAL(rprec), PARAMETER :: p5 = cp5
  REAL(rprec), PARAMETER :: two = c2p0

  REAL(rprec) :: bsubvvac
  REAL(rprec) :: pi2
  REAL(rprec) :: pi3
  REAL(rprec) :: pi4
  REAL(rprec) :: alp
  REAL(rprec) :: alu
  REAL(rprec) :: alv
  REAL(rprec) :: alvp
  REAL(rprec) :: onp
  REAL(rprec) :: onp2

  logical :: precal_done

  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: potvac

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bvecsav
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: amatsav

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bexni

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: brv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bphiv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bzv

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsqvac

  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: bsqsav

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: r1b
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rub
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rvb
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: z1b
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: zub
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: zvb

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bexu
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bexv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bexn

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: auu
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: auv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: avv

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: snr
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: snv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: snz

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: drv

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: guu_b
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: guv_b
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: gvv_b

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rzb2

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rcosuv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: rsinuv

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bredge
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bpedge
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bzedge

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: raxis_nestor
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: zaxis_nestor

  ! from vacuum
  REAL(rprec), ALLOCATABLE :: bsubu(:), bsubv(:), potu(:), potv(:)
  REAL(rprec), ALLOCATABLE :: amatrix(:)



  CONTAINS

subroutine allocate_nestor
  integer :: istat1, istat2, i, ndim

  precal_done = .false.

  ! nuv2 = nznt in read_indata

  ALLOCATE (amatsav(mnpd2*mnpd2), bvecsav(mnpd2), &
            bsqsav(nuv2,3), potvac(2*mnpd),          &
            raxis_nestor(nv), zaxis_nestor(nv), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #3 in allocate_nestor'

  ALLOCATE (brv(nuv2), bphiv(nuv2), bzv(nuv2), bsqvac(nuv2), stat=istat2)
  IF (istat2.ne.0) STOP 'allocation error #2 in allocate_nestor'

 brv=0
 bphiv=0
 bzv=0

 bsqvac=0

 ! from vacuum()
 ALLOCATE (amatrix(mnpd2*mnpd2), bsubu(nuv2), bsubv(nuv2), potu(nuv2), potv(nuv2), stat = i)
  IF (i .ne. 0) STOP 'Allocation error in vacuum'

  ALLOCATE (bexu(nuv2), bexv(nuv2), bexn(nuv2), bexni(nuv2),                    &
            r1b(nuv), rub(nuv2), rvb(nuv2), z1b(nuv), zub(nuv2), zvb(nuv2), &
            auu(nuv2), auv(nuv2), avv(nuv2), snr(nuv2), snv(nuv2), snz(nuv2), &
            drv(nuv2), guu_b(nuv2), guv_b(nuv2), gvv_b(nuv2), &
            rzb2(nuv), rcosuv(nuv), rsinuv(nuv), stat=i)
  IF (i .ne. 0) STOP 'Allocation error in vacuum'

  ! from precal
! ALLOCATE PERSISTENT ARRAYS
  IF (nv == 1) THEN
     ! AXISYMMETRIC CASE: DO FP SUM TO INTEGRATE IN V
     nvper = 64
     nuv_tan = 2*nu*nvper
  ELSE
     nvper = nfper
     nuv_tan = 2*nuv
  END IF

  ALLOCATE (tanu(nuv_tan), tanv(nuv_tan),                           &
       sinper(nvper), cosper(nvper), sinuv(nuv), cosuv(nuv),        &
       sinu(0:mf,nu), cosu(0:mf,nu), sinv(-nf:nf,nv),               &
       cosv(-nf:nf,nv), sinui(0:mf,nu2), cosui(0:mf,nu2),           &
       cmns(0:(mf+nf),0:mf,0:nf), csign(-nf:nf),                    &
       sinu1(nuv2,0:mf), cosu1(nuv2,0:mf),                          &
       sinv1(nuv2,0:nf), cosv1(nuv2,0:nf), imirr(nuv),              &
       xmpot(mnpd), xnpot(mnpd), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error in precal'

end subroutine allocate_nestor

subroutine free_mem_nestor

  integer :: istat1, istat3, istat4, i

  IF (ALLOCATED(amatsav)) then
      DEALLOCATE (amatsav, bvecsav, bsqsav, potvac,  &
                  raxis_nestor, zaxis_nestor, stat=istat3)
      IF (istat3.ne.0) STOP 'dealloc error #3 in free_mem_nestor'
  end if

  IF (ALLOCATED(brv)) then
     DEALLOCATE (brv, bphiv, bzv, bsqvac, stat=istat4)
     IF (istat4.ne.0) STOP 'dealloc error #4 in free_mem_nestor'
  end if

  ! from vacuum()
IF (ALLOCATED(bexu)) then
     DEALLOCATE (bexu, bexv, bexn, bexni, &
        r1b, rub, rvb, z1b, zub, zvb,   &
        auu, auv, avv, snr, snv, snz, &
        drv, guu_b, guv_b, gvv_b, &
        rzb2, rcosuv, rsinuv, stat=i)
     IF (i .ne. 0) STOP 'Deallocation error in vacuum'
  end if

  if (allocated(amatrix)) then
     DEALLOCATE (amatrix, bsubu, bsubv, potu, potv, stat = i)
     IF (i .ne. 0) STOP 'Deallocation error in vacuum'
  end if

  IF (ALLOCATED(tanu)) then
     DEALLOCATE(tanu, tanv, sinper, cosper, sinuv, cosuv, cmns,         &
                sinu, cosu, sinv, cosv, sinui, cosui, csign, sinu1,     &
                cosu1, sinv1, cosv1, imirr, xmpot, xnpot, stat=istat1)
     IF (istat1 .ne. 0) STOP 'Deallocation error in vacuum'
  end if
end subroutine free_mem_nestor


END MODULE vacmod
