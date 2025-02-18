!> \file
MODULE vacmod

  USE vacmod0
  USE vac_persistent
  USE vparams, ONLY: zero, one, c2p0, cp5, mu0
  use vmec_input, only: vac_1_2

  IMPLICIT NONE

  ! integer, save :: icall = 0

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

  real(rprec), dimension(:), allocatable :: m_map_wrt
  real(rprec), dimension(:), allocatable :: n_map_wrt

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bvecsav
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: amatsav

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bexni

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: brv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bphiv
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bzv

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: bsqvac

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

  REAL(rprec), DIMENSION(:), ALLOCATABLE :: raxis_nestor
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: zaxis_nestor

  ! from vacuum
  REAL(rprec), ALLOCATABLE :: bsubu(:)
  REAL(rprec), ALLOCATABLE :: bsubv(:)
  REAL(rprec), ALLOCATABLE :: potu(:)
  REAL(rprec), ALLOCATABLE :: potv(:)
  REAL(rprec), ALLOCATABLE :: amatrix(:)

  ! from surface
  REAL(rprec), ALLOCATABLE, DIMENSION(:) :: ruu
  REAL(rprec), ALLOCATABLE, DIMENSION(:) :: ruv
  REAL(rprec), ALLOCATABLE, DIMENSION(:) :: rvv
  REAL(rprec), ALLOCATABLE, DIMENSION(:) :: zuu
  REAL(rprec), ALLOCATABLE, DIMENSION(:) :: zuv
  REAL(rprec), ALLOCATABLE, DIMENSION(:) :: zvv

  ! from bextern
  REAL(rprec), ALLOCATABLE :: brad(:)
  REAL(rprec), ALLOCATABLE :: bphi(:)
  REAL(rprec), ALLOCATABLE :: bz(:)

  ! from tolicu
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: xpts

  ! from scalpot
  REAL(rprec), ALLOCATABLE :: grpmn(:)
  real(rprec), dimension(:), allocatable :: grpmn_m_map_wrt
  real(rprec), dimension(:), allocatable :: grpmn_n_map_wrt

  REAL(rprec), ALLOCATABLE :: gstore(:)
  REAL(rprec), ALLOCATABLE :: green(:,:)
  REAL(rprec), ALLOCATABLE :: greenp(:,:)

  ! from analyt
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: r0p
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: r1p
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: r0m
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: r1m
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: sqrtc
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: sqrta
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: tlp2
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: tlp1
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: tlp
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: tlm2
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: tlm1
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: tlm
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: adp
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: adm
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: cma
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: ra1p
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: ra1m
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: slm
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: slp
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: tlpm
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: slpm
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: delt1u
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: azp1u
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: azm1u
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: cma11u
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: sqad1u
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: sqad2u

  real(rprec), dimension(:,:), allocatable :: all_tlp
  real(rprec), dimension(:,:), allocatable :: all_tlm
  real(rprec), dimension(:,:), allocatable :: all_slp
  real(rprec), dimension(:,:), allocatable :: all_slm

  ! from greenf
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: gsave
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ga1
  REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: ga2
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: dsave

  ! from fourp
  REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: g1
  REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: g2

  ! from fouri
  REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: bcos
  REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: bsin
  REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: source

  REAL(rprec), ALLOCATABLE :: actemp(:,:,:,:)
  REAL(rprec), ALLOCATABLE :: astemp(:,:,:,:)

contains

subroutine allocate_nestor
  integer :: istat = 0
  integer :: istat1 = 0
  integer :: istat2 = 0
  integer :: i = 0
  integer :: l = 0
  integer :: m = 0
  integer :: ip = 0

  precal_done = .false.

  ALLOCATE (amatsav(mnpd2*mnpd2), bvecsav(mnpd2), &
            potvac(2*mnpd),          &
            raxis_nestor(nv), zaxis_nestor(nv), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error #3 in allocate_nestor'

  allocate(m_map_wrt(mnpd2), n_map_wrt(mnpd2), stat=istat1)
  if (istat1 .ne. 0) then
     stop 'could not allocate m_map, n_map'
  end if

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
  allocate(tanu_1d(2*nu))
  if (nv .eq. 1) then
    allocate(tanv_1d(nvper))
  else
    allocate(tanv_1d(nv))
  end if

  ALLOCATE (tanu(nuv_tan), tanv(nuv_tan),                           &
       sinper(nvper), cosper(nvper), sinuv(nuv), cosuv(nuv),        &
       sinu(0:mf,nu), cosu(0:mf,nu), sinv(-nf:nf,nv),               &
       cosv(-nf:nf,nv), sinui(0:mf,nu2), cosui(0:mf,nu2),           &
       cmns(0:(mf+nf),0:mf,0:nf), csign(-nf:nf),                    &
       sinu1(nuv2,0:mf), cosu1(nuv2,0:mf),                          &
       sinv1(nuv2,0:nf), cosv1(nuv2,0:nf), imirr(nuv),              &
       xmpot(mnpd), xnpot(mnpd), stat=istat1)
  IF (istat1.ne.0) STOP 'allocation error in precal'

  ! from surface()
  ALLOCATE (ruu(nuv2), ruv(nuv2), rvv(nuv2), zuu(nuv2), zuv(nuv2), zvv(nuv2), stat = i)
  IF (i .NE. 0) STOP 'Allocation error in SURFACE'

  ! from bextern
  ALLOCATE (brad(nuv2), bphi(nuv2), bz(nuv2), stat=i)
  IF (i .ne. 0) STOP 'allocation error in bextern'

  ! from tolicu
  ! need nvp+1 for "virtual" point at last index,
  ! which is equal to first point for a closed curve
  ALLOCATE (xpts(3,nvp+1), stat=i)
  IF (i .ne. 0) STOP ' allocation error in tolicu'

  ! from scalpot
  ALLOCATE (grpmn(nuv2*mnpd2), stat=ip)
  IF (ip .ne. 0) STOP 'GRPMN: Allocation error in scalpot'

  ALLOCATE (grpmn_m_map_wrt(nuv2*mnpd2), grpmn_n_map_wrt(nuv2*mnpd2), stat=ip)
  IF (ip .ne. 0) STOP 'GRPMN: Allocation error in scalpot'

  ALLOCATE (gstore(nuv), green(nuv,nuv2), greenp(nuv,nuv2), stat=istat)
  if (istat.ne.0) then
     ! Below loop over nuv2 was previously chunked.
     ! Therefore, some extra care shall be used here to make sure
     ! everything still fits into memory...
     stop 'green allocation error in scalpot.'
  end if

  ! from analyt
  ALLOCATE (r0p(nuv2), r1p(nuv2), r0m(nuv2), r1m(nuv2),             &
            sqrtc(nuv2), sqrta(nuv2), tlp2(nuv2), tlp1(nuv2),       &
            tlp(nuv2), tlm2(nuv2), tlm1(nuv2), tlm(nuv2), adp(nuv2),&
            adm(nuv2), cma(nuv2), ra1p(nuv2), ra1m(nuv2), slm(nuv2),&
            slp(nuv2), tlpm(nuv2), slpm(nuv2), delt1u(nuv2),        &
            azp1u(nuv2), azm1u(nuv2), cma11u(nuv2), sqad1u(nuv2),   &
            sqad2u(nuv2), stat = l)
  IF (l .ne. 0) STOP 'Allocation error in SUBROUTINE analyt'

  allocate(all_tlp(0:mf+nf, nuv2), all_tlm(0:mf+nf, nuv2), &
           all_slp(0:mf+nf, nuv2), all_slm(0:mf+nf, nuv2), stat=l)
  IF (l .ne. 0) STOP 'Allocation error for debugging analyt'

  ! from greenf
  ALLOCATE (gsave(nuv), ga1(nvper,nuv), ga2(nvper,nuv), dsave(nuv), stat=i)
  IF (i .ne. 0) STOP 'allocation error in greenf'

  ! print *, "allocations; ndim=",ndim

  ! from fourp
  ALLOCATE (g1(nuv2,0:nf,ndim), g2(nuv2,0:nf,ndim), stat=m)
  IF (m .NE. 0) STOP 'Allocation error in fourp'

  ! from fouri
  ALLOCATE (bcos(nu2,-nf:nf,ndim), bsin(nu2,-nf:nf,ndim),           &
     actemp(mnpd,-nf:nf,nu3,ndim), astemp(mnpd,-nf:nf,nu3,ndim),    &
     source(nv,nu2,ndim), stat = i)
  IF (i .ne. 0) STOP 'allocation error in fouri'

end subroutine allocate_nestor



subroutine free_mem_nestor

  integer :: istat1 = 0
  integer :: istat3 = 0
  integer :: istat4 = 0
  integer :: i = 0
  integer :: l = 0
  integer :: m = 0

  IF (ALLOCATED(amatsav)) then
      DEALLOCATE (amatsav, bvecsav, potvac,  &
                  raxis_nestor, zaxis_nestor, stat=istat3)
      IF (istat3.ne.0) STOP 'dealloc error #3 in free_mem_nestor'
  end if

  if (allocated(m_map_wrt)) then
     deallocate(m_map_wrt, n_map_wrt, stat=istat3)
     if (istat3 .ne. 0) stop 'could not deallocate m_map, n_map'
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

  if (allocated(tanu_1d)) then
    deallocate(tanu_1d)
    deallocate(tanv_1d)
  end if

  IF (ALLOCATED(tanu)) then
     DEALLOCATE(tanu, tanv, sinper, cosper, sinuv, cosuv, cmns,         &
                sinu, cosu, sinv, cosv, sinui, cosui, csign, sinu1,     &
                cosu1, sinv1, cosv1, imirr, xmpot, xnpot, stat=istat1)
     IF (istat1 .ne. 0) STOP 'Deallocation error in vacuum'
  end if

  ! from surface
  if (allocated(ruu)) then
     DEALLOCATE (ruu, ruv, rvv, zuu, zuv, zvv, stat=i)
  end if

  ! from bextern
  if (allocated(brad)) then
     DEALLOCATE (brad, bphi, bz)
  end if

  ! added for tolicu
  if (allocated(xpts)) then
     deallocate(xpts)
  end if

  ! from scalpot
  if (allocated(grpmn)) then
    DEALLOCATE (grpmn)
    DEALLOCATE (green, greenp, gstore)
  end if

  if (allocated(grpmn_m_map_wrt)) then
     deallocate(grpmn_m_map_wrt, grpmn_n_map_wrt)
  end if

  ! from analyt
  if (allocated(r0p)) then
    DEALLOCATE (r0p, r1p, r0m, r1m, sqrtc, sqrta, tlp2, tlp1,         &
              tlp, tlm2, tlm1, tlm, adp, adm, cma, ra1p, ra1m, slm,   &
              slp, tlpm, slpm, delt1u, azp1u, azm1u, cma11u, sqad1u,  &
              sqad2u, stat = l)
  end if

  ! from greenf
  if (allocated(gsave)) then
    DEALLOCATE (gsave, ga1, ga2, dsave, stat=i)
  end if

  ! from fourp
  if (allocated(g1)) then
    DEALLOCATE (g1, g2, stat=m)
  end if

  ! from fouri
  if (allocated(bcos)) then
    DEALLOCATE (bcos, bsin, actemp, astemp, source, stat=i)
  end if


end subroutine free_mem_nestor


END MODULE vacmod
