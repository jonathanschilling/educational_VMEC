!> \file
!> \brief Perform Fourier integrals of Green's function kernel.

!> \brief Perform Fourier integrals of Green's function kernel.
!>
!> @param grpmn
!> @param grp
SUBROUTINE fourp (grpmn, grp)
  USE vacmod, vm_grpmn => grpmn

  use dbgout

  IMPLICIT NONE

  REAL(rprec), INTENT(inout) :: grpmn(0:mf,-nf:nf,nuv2,ndim)
  REAL(rprec), INTENT(in)    :: grp(nuv,nuv2)

  INTEGER :: n, kv, ku, ip, iuv, m, ireflect, isym
  REAL(rprec) :: cosm, sinm, cosn, sinn, kernel, gcos, gsin

  IF (ndim .GT. 2) STOP 'NDIM > 2'

  DO ku = 1,nu2
     g1 = 0
     g2 = 0

     ! PERFORM KV (TOROIDAL ANGLE) TRANSFORM (OVER UNPRIMED MESH IN EQ. 2.14)
     ! THUS, THE (m,n) INDEX HERE CORRESPONDS TO THE FIRST INDEX OF AMATRIX
     !
     ! NOTE: THE .5 FACTOR (IN COSN,SINN) ACCOUNTS FOR THE SUM IN KERNEL.
     ! ON ENTRY THE FIRST TIME, GRPMN IS SIN,COS * Kmn(analytic)
     !
     ! THE 3rd INDEX OF GRPMN IS THE PRIMED U,V MESH COORDINATE
     DO kv = 1,nv
        iuv = kv+nv*(ku-1)
        ireflect = imirr(iuv)

        DO n = 0,nf
           cosn = p5*onp*cosv(n,kv)
           sinn = p5*onp*sinv(n,kv)

           DO isym = 1, ndim ! ndim == 1 for stellarator-symmetry, 2 for asymmetric (lasym=T)
              DO ip = 1,nuv2
                 IF (isym .eq. 1) THEN ! only contrib for stellarator-symmetry
                    ! anti-symmetric part (u,v -> -u,-v)

                    kernel = grp(iuv,ip) - grp(ireflect,ip)

!                     if (iuv .eq. ireflect) then
!                       print *, "reflecting on self at iuv=", iuv, " ip=", ip, " n=", n, " ku=", ku, " kv=", kv, &
!                                " => kernel=", kernel
!                     end if

!                      if (n.eq.0) then
!                         ! fourp_kernel.dat
!                         print *, iuv, ip, kernel
!                      end if
                 ELSEIF (isym .eq. 2) THEN
                    ! symmetric part
                    kernel = grp(iuv,ip) + grp(ireflect,ip)
                 END IF

                 g1(ip,n,isym)=g1(ip,n,isym) + cosn*kernel
                 g2(ip,n,isym)=g2(ip,n,isym) + sinn*kernel
              END DO ! ip
           END DO ! isym
        END DO ! n
     END DO ! kv

     ! PERFORM KU (POLOIDAL ANGLE) TRANFORM [COMPLETE SIN(mu-nv) / COS(mu-nv) TRANSFORM]
     DO m = 0,mf
        DO isym = 1, ndim
           IF (isym .EQ. 1) THEN
              cosm = -cosui(m,ku)
              sinm =  sinui(m,ku)
           ELSEIF (isym .EQ. 2) THEN
              sinm = cosui(m,ku)
              cosm = sinui(m,ku)
           END IF

           DO n= 0,nf
             DO ip = 1,nuv2
                 gcos = g1(ip,n,isym)*sinm
                 gsin = g2(ip,n,isym)*cosm ! has -1
                 grpmn(m, n,ip,isym)    = grpmn(m, n,ip,isym) + gcos + gsin
                 IF (n .NE. 0) THEN
                    grpmn(m,-n,ip,isym) = grpmn(m,-n,ip,isym) + gcos - gsin
                 ENDIF
             END DO ! ip
           END DO ! n
        END DO ! isym
     END DO ! m

  END DO ! ku

  if (open_dbg_context("vac1n_fourp", id=icall)) then

    ! NOTE: This is the sum of the results from analyt and fourp!
    call add_real_4d("grpmn", mf1, nf1, nv, nu3, grpmn)

    call close_dbg_out()
  end if

END SUBROUTINE fourp
