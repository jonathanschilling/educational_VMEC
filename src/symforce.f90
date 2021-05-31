!> \file
!> \brief Symmetrize forces on restricted \f$\theta\f$ interval \f$(0 \leq \theta \leq \pi \f$
!>        so cos, sin integrals can be performed.

!> \brief Symmetrize forces on restricted \f$\theta\f$ interval \f$(0 \leq \theta \leq \pi \f$
!>        so cos, sin integrals can be performed.
!>
!> @param ars
!> @param brs
!> @param crs
!> @param azs
!> @param bzs
!> @param czs
!> @param bls
!> @param cls
!> @param rcs
!> @param zcs
!> @param ara
!> @param bra
!> @param cra
!> @param aza
!> @param bza
!> @param cza
!> @param bla
!> @param cla
!> @param rca
!> @param zca
SUBROUTINE symforce(ars, brs, crs, azs, bzs, czs, bls, cls, rcs, zcs, &
                    ara, bra, cra, aza, bza, cza, bla, cla, rca, zca)
  USE vmec_main, p5 => cp5
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(inout) :: &
     ars, brs, crs, azs, bzs, czs, bls, cls, rcs, zcs
  REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(out) ::   &
     ara, bra, cra, aza, bza, cza, bla, cla, rca, zca

  INTEGER :: mpar, ir, i, jk, jka
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: ars_0, brs_0, azs_0, &
                bzs_0, bls_0, rcs_0, zcs_0, crs_0, czs_0, cls_0

  i = ns*nzeta
  ALLOCATE (ars_0(i), brs_0(i), azs_0(i), bzs_0(i), bls_0(i),         &
            rcs_0(i), zcs_0(i), crs_0(i), czs_0(i), cls_0(i), stat=ir)

  ! SYMMETRIZE FORCES ON RESTRICTED THETA INTERVAL (0 <= u <= pi)
  ! SO COS,SIN INTEGRALS CAN BE PERFORMED. FOR EXAMPLE,
  !
  ! ARS(v,u) = .5*( ARS(v,u) + ARS(-v,-u) )     ! * COS(mu - nv)
  ! ARA(v,u) = .5*( ARS(v,u) - ARS(-v,-u) )     ! * SIN(mu - nv)
  DO mpar = 0, 1
     DO i = 1, ntheta2
        ir = ntheta1 + 2 - i                 !-theta
        IF (i .eq. 1) ir = 1
        DO jk = 1, ns*nzeta
           jka = ireflect(jk)                !-zeta
           ara(jk,i,mpar) = p5*(ars(jk,i,mpar)-ars(jka,ir,mpar))
           ars_0(jk)      = p5*(ars(jk,i,mpar)+ars(jka,ir,mpar))
           bra(jk,i,mpar) = p5*(brs(jk,i,mpar)+brs(jka,ir,mpar))
           brs_0(jk)      = p5*(brs(jk,i,mpar)-brs(jka,ir,mpar))
           aza(jk,i,mpar) = p5*(azs(jk,i,mpar)+azs(jka,ir,mpar))
           azs_0(jk)      = p5*(azs(jk,i,mpar)-azs(jka,ir,mpar))
           bza(jk,i,mpar) = p5*(bzs(jk,i,mpar)-bzs(jka,ir,mpar))
           bzs_0(jk)      = p5*(bzs(jk,i,mpar)+bzs(jka,ir,mpar))
           bla(jk,i,mpar) = p5*(bls(jk,i,mpar)-bls(jka,ir,mpar))
           bls_0(jk)      = p5*(bls(jk,i,mpar)+bls(jka,ir,mpar))
           rca(jk,i,mpar) = p5*(rcs(jk,i,mpar)-rcs(jka,ir,mpar))
           rcs_0(jk)      = p5*(rcs(jk,i,mpar)+rcs(jka,ir,mpar))
           zca(jk,i,mpar) = p5*(zcs(jk,i,mpar)+zcs(jka,ir,mpar))
           zcs_0(jk)      = p5*(zcs(jk,i,mpar)-zcs(jka,ir,mpar))
        END DO

        ars(:,i,mpar) = ars_0(:)
        brs(:,i,mpar) = brs_0(:)
        azs(:,i,mpar) = azs_0(:)
        bzs(:,i,mpar) = bzs_0(:)
        bls(:,i,mpar) = bls_0(:)
        rcs(:,i,mpar) = rcs_0(:)
        zcs(:,i,mpar) = zcs_0(:)

        IF (lthreed) THEN
           DO jk = 1, ns*nzeta
              jka = ireflect(jk)
              cra(jk,i,mpar) = p5*(crs(jk,i,mpar)+crs(jka,ir,mpar))
              crs_0(jk)      = p5*(crs(jk,i,mpar)-crs(jka,ir,mpar))
              cza(jk,i,mpar) = p5*(czs(jk,i,mpar)-czs(jka,ir,mpar))
              czs_0(jk)      = p5*(czs(jk,i,mpar)+czs(jka,ir,mpar))
              cla(jk,i,mpar) = p5*(cls(jk,i,mpar)-cls(jka,ir,mpar))
              cls_0(jk)      = p5*(cls(jk,i,mpar)+cls(jka,ir,mpar))
           END DO

           crs(:,i,mpar) = crs_0(:)
           czs(:,i,mpar) = czs_0(:)
           cls(:,i,mpar) = cls_0(:)
        ENDIF

     END DO
  END DO

  DEALLOCATE (ars_0, brs_0, azs_0, bzs_0, bls_0,          &
              rcs_0, zcs_0, crs_0, czs_0, cls_0, stat=ir)

END SUBROUTINE symforce

!> \brief Symmetrize some quantities so that they can be output (?)
!>
!> @param bsq
!> @param gsqrt
!> @param bsubu
!> @param bsubv
!> @param bsupu
!> @param bsupv
!> @param bsubs
!> @param bsqa
!> @param gsqrta
!> @param bsubua
!> @param bsubva
!> @param bsupua
!> @param bsupva
!> @param bsubsa
SUBROUTINE symoutput (bsq , gsqrt , bsubu , bsubv ,bsupu,  bsupv , bsubs , &
                      bsqa, gsqrta, bsubua, bsubva,bsupua, bsupva, bsubsa   )

  USE vmec_main, p5 => cp5

  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns*nzeta,ntheta3), INTENT(inout) :: &
     bsq, gsqrt, bsubu, bsubv, bsupu, bsupv, bsubs
  REAL(rprec), DIMENSION(ns*nzeta,ntheta3), INTENT(out)   :: &
     bsqa,gsqrta,bsubua,bsubva,bsupua,bsupva,bsubsa

  INTEGER :: ir, i, jk, jka
  REAL(rprec), DIMENSION(ns*nzeta) :: bsq_0, gsqrt_0, bsubu_0, &
      bsubv_0, bsupu_0, bsupv_0, bsubs_0

  ! SYMMETRIZE FORCES ON RESTRICTED THETA INTERVAL (0 <= u <= pi)
  ! SO COS,SIN INTEGRALS CAN BE PERFORMED. FOR EXAMPLE,
  !
  ! BSQ-S(v,u) = .5*( BSQ(v,u) + BSQ(-v,-u) )     ! * COS(mu - nv)
  ! BSQ-A(v,u) = .5*( BSQ(v,u) - BSQ(-v,-u) )     ! * SIN(mu - nv)
  !
  ! FOR BSUBS, THIS IS REVERSED, S-PIECE ~ SIN, A-PIECE ~ COS
  DO i = 1, ntheta2
     ir = ntheta1 + 2 - i                 !-theta
     IF (i == 1) ir = 1
     DO jk = 1, ns*nzeta
        jka = ireflect(jk)                !-zeta
        bsqa(jk,i)      = p5*(bsq(jk,i)     -bsq(jka,ir))
        bsq_0(jk)       = p5*(bsq(jk,i)     +bsq(jka,ir))
        gsqrta(jk,i)    = p5*(gsqrt(jk,i)   -gsqrt(jka,ir))
        gsqrt_0(jk)     = p5*(gsqrt(jk,i)   +gsqrt(jka,ir))
        bsubua(jk,i)    = p5*(bsubu(jk,i)   -bsubu(jka,ir))
        bsubu_0(jk)     = p5*(bsubu(jk,i)   +bsubu(jka,ir))
        bsubva(jk,i)    = p5*(bsubv(jk,i)   -bsubv(jka,ir))
        bsubv_0(jk)     = p5*(bsubv(jk,i)   +bsubv(jka,ir))
        bsupua(jk,i)    = p5*(bsupu(jk,i)   -bsupu(jka,ir))
        bsupu_0(jk)     = p5*(bsupu(jk,i)   +bsupu(jka,ir))
        bsupva(jk,i)    = p5*(bsupv(jk,i)   -bsupv(jka,ir))
        bsupv_0(jk)     = p5*(bsupv(jk,i)   +bsupv(jka,ir))
        ! Dominant symmetry reversed
        bsubsa(jk,i)    = p5*(bsubs(jk,i)   +bsubs(jka,ir))
        bsubs_0(jk)     = p5*(bsubs(jk,i)   -bsubs(jka,ir))
     END DO

     bsq(:,i)      = bsq_0(:)
     gsqrt(:,i)    = gsqrt_0(:)
     bsubu(:,i)    = bsubu_0(:)
     bsubv(:,i)    = bsubv_0(:)
     bsupu(:,i)    = bsupu_0(:)
     bsupv(:,i)    = bsupv_0(:)
     bsubs(:,i)    = bsubs_0(:)

  END DO

END SUBROUTINE symoutput

