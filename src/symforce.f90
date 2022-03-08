!> \file
!> \brief Symmetrize forces on restricted \f$\theta\f$ interval \f$(0 \leq \theta \leq \pi \f$
!>        so cos, sin integrals can be performed.

!> \brief Symmetrize forces on restricted \f$\theta\f$ interval \f$(0 \leq \theta \leq \pi \f$
!>        so cos, sin integrals can be performed.
!>
!> @param ars contribution to A^R       with even parity (equal to parity of A^R       in symmetric case)
!> @param brs contribution to B^R       with  odd parity (equal to parity of B^R       in symmetric case)
!> @param crs contribution to C^R       with  odd parity (equal to parity of C^R       in symmetric case)
!> @param azs contribution to A^Z       with  odd parity (equal to parity of A^Z       in symmetric case)
!> @param bzs contribution to B^Z       with even parity (equal to parity of B^Z       in symmetric case)
!> @param czs contribution to C^Z       with even parity (equal to parity of C^Z       in symmetric case)
!> @param bls contribution to B^lambda  with even parity (equal to parity of B^lambda  in symmetric case)
!> @param cls contribution to C^lambda  with even parity (equal to parity of C^lambda  in symmetric case)
!> @param rcs contribution to F^{R_con} with even parity (equal to parity of F^{R_con} in symmetric case)
!> @param zcs contribution to F^{Z_con} with  odd parity (equal to parity of F^{Z_con} in symmetric case)
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

  use dbgout

  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(inout) :: &
     ars, brs, crs, azs, bzs, czs, bls, cls, rcs, zcs
  REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(out) ::   &
     ara, bra, cra, aza, bza, cza, bla, cla, rca, zca

  INTEGER :: mpar, ir, i, jk, jka, j, k
  REAL(rprec), DIMENSION(:), ALLOCATABLE :: ars_0, brs_0, azs_0, &
                bzs_0, bls_0, rcs_0, zcs_0, crs_0, czs_0, cls_0

  logical :: dbg_symforce

  i = ns*nzeta
  ALLOCATE (ars_0(i), brs_0(i), azs_0(i), bzs_0(i), bls_0(i),         &
            rcs_0(i), zcs_0(i), crs_0(i), czs_0(i), cls_0(i), stat=ir)

  dbg_symforce = open_dbg_context("symforce")
  if (dbg_symforce) then

    call add_real_4d("ars", ns, 2, nzeta, ntheta3, ars, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("brs", ns, 2, nzeta, ntheta3, brs, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("azs", ns, 2, nzeta, ntheta3, azs, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("bzs", ns, 2, nzeta, ntheta3, bzs, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("bls", ns, 2, nzeta, ntheta3, bls, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("rcs", ns, 2, nzeta, ntheta3, rcs, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("zcs", ns, 2, nzeta, ntheta3, zcs, order=(/ 1, 3, 4, 2 /) )

    if (lthreed) then
      call add_real_4d("crs", ns, 2, nzeta, ntheta3, crs, order=(/ 1, 3, 4, 2 /) )
      call add_real_4d("czs", ns, 2, nzeta, ntheta3, czs, order=(/ 1, 3, 4, 2 /) )
      call add_real_4d("cls", ns, 2, nzeta, ntheta3, cls, order=(/ 1, 3, 4, 2 /) )
    else
      call add_null("crs")
      call add_null("czs")
      call add_null("cls")
    end if

  end if ! dump_symforce

  ! SYMMETRIZE FORCES ON RESTRICTED THETA INTERVAL (0 <= u <= pi)
  ! SO COS,SIN INTEGRALS CAN BE PERFORMED. FOR EXAMPLE,
  !
  ! ARS(v,u) = .5*( ARS(v,u) + ARS(-v,-u) )     ! * COS(mu - nv)
  ! ARA(v,u) = .5*( ARS(v,u) - ARS(-v,-u) )     ! * SIN(mu - nv)
  !
  ! See also Theorem 26 (Parity Decomposition) in Boyd, "Chebychev and Fourier Spectral Methods".
  DO mpar = 0, 1
     DO i = 1, ntheta2

        ! (ntheta1 + 1) - (i-1)
        ir = ntheta1 + 2 - i                 !-theta
        IF (i .eq. 1) ir = 1


        DO jk = 1, ns*nzeta
           jka = ireflect(jk)                !-zeta

           ara(jk,i,mpar) = p5*(ars(jk,i,mpar)-ars(jka,ir,mpar)) ! std. parity
           ars_0(jk)      = p5*(ars(jk,i,mpar)+ars(jka,ir,mpar))
           bra(jk,i,mpar) = p5*(brs(jk,i,mpar)+brs(jka,ir,mpar)) ! rev. parity
           brs_0(jk)      = p5*(brs(jk,i,mpar)-brs(jka,ir,mpar))
           ! cr only in lthreed case (see below)

           aza(jk,i,mpar) = p5*(azs(jk,i,mpar)+azs(jka,ir,mpar)) ! rev. parity
           azs_0(jk)      = p5*(azs(jk,i,mpar)-azs(jka,ir,mpar))
           bza(jk,i,mpar) = p5*(bzs(jk,i,mpar)-bzs(jka,ir,mpar)) ! std. parity
           bzs_0(jk)      = p5*(bzs(jk,i,mpar)+bzs(jka,ir,mpar))
           ! cz only in lthreed case (see below)

           bla(jk,i,mpar) = p5*(bls(jk,i,mpar)-bls(jka,ir,mpar)) ! std. parity
           bls_0(jk)      = p5*(bls(jk,i,mpar)+bls(jka,ir,mpar))
           ! cl only in lthreed case (see below)

           rca(jk,i,mpar) = p5*(rcs(jk,i,mpar)-rcs(jka,ir,mpar)) ! std. parity
           rcs_0(jk)      = p5*(rcs(jk,i,mpar)+rcs(jka,ir,mpar))
           zca(jk,i,mpar) = p5*(zcs(jk,i,mpar)+zcs(jka,ir,mpar)) ! rev. parity
           zcs_0(jk)      = p5*(zcs(jk,i,mpar)-zcs(jka,ir,mpar))
        END DO

        ! need to store symmetric part in temp array,
        ! since reflected indices must not be overwritten above!
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

              cra(jk,i,mpar) = p5*(crs(jk,i,mpar)+crs(jka,ir,mpar)) ! rev. parity
              crs_0(jk)      = p5*(crs(jk,i,mpar)-crs(jka,ir,mpar))

              cza(jk,i,mpar) = p5*(czs(jk,i,mpar)-czs(jka,ir,mpar)) ! std. parity
              czs_0(jk)      = p5*(czs(jk,i,mpar)+czs(jka,ir,mpar))

              cla(jk,i,mpar) = p5*(cls(jk,i,mpar)-cls(jka,ir,mpar)) ! std. parity
              cls_0(jk)      = p5*(cls(jk,i,mpar)+cls(jka,ir,mpar))
           END DO

           crs(:,i,mpar) = crs_0(:)
           czs(:,i,mpar) = czs_0(:)
           cls(:,i,mpar) = cls_0(:)
        ENDIF

     END DO

     ! clear remainder of arrays for debug output
     DO i = ntheta2+1, ntheta3
       ara(:,i,mpar) = 0.0_dp
       bra(:,i,mpar) = 0.0_dp
       aza(:,i,mpar) = 0.0_dp
       bza(:,i,mpar) = 0.0_dp
       bla(:,i,mpar) = 0.0_dp
       rca(:,i,mpar) = 0.0_dp
       zca(:,i,mpar) = 0.0_dp
       if (lthreed) then
         cra(:,i,mpar) = 0.0_dp
         cza(:,i,mpar) = 0.0_dp
         cla(:,i,mpar) = 0.0_dp
       end if
     end do
  END DO

  DEALLOCATE (ars_0, brs_0, azs_0, bzs_0, bls_0,          &
              rcs_0, zcs_0, crs_0, czs_0, cls_0, stat=ir)

  if (dbg_symforce) then
    call add_real_4d("ars_out", ns, 2, nzeta, ntheta3, ars, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("ara_out", ns, 2, nzeta, ntheta3, ara, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("brs_out", ns, 2, nzeta, ntheta3, brs, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("bra_out", ns, 2, nzeta, ntheta3, bra, order=(/ 1, 3, 4, 2 /) )

    call add_real_4d("azs_out", ns, 2, nzeta, ntheta3, azs, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("aza_out", ns, 2, nzeta, ntheta3, aza, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("bzs_out", ns, 2, nzeta, ntheta3, bzs, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("bza_out", ns, 2, nzeta, ntheta3, bza, order=(/ 1, 3, 4, 2 /) )

    call add_real_4d("bls_out", ns, 2, nzeta, ntheta3, bls, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("bla_out", ns, 2, nzeta, ntheta3, bla, order=(/ 1, 3, 4, 2 /) )

    call add_real_4d("rcs_out", ns, 2, nzeta, ntheta3, rcs, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("rca_out", ns, 2, nzeta, ntheta3, rca, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("zcs_out", ns, 2, nzeta, ntheta3, zcs, order=(/ 1, 3, 4, 2 /) )
    call add_real_4d("zca_out", ns, 2, nzeta, ntheta3, zca, order=(/ 1, 3, 4, 2 /) )

    if (lthreed) then
      call add_real_4d("crs_out", ns, 2, nzeta, ntheta3, crs, order=(/ 1, 3, 4, 2 /) )
      call add_real_4d("cra_out", ns, 2, nzeta, ntheta3, cra, order=(/ 1, 3, 4, 2 /) )
      call add_real_4d("czs_out", ns, 2, nzeta, ntheta3, czs, order=(/ 1, 3, 4, 2 /) )
      call add_real_4d("cza_out", ns, 2, nzeta, ntheta3, cza, order=(/ 1, 3, 4, 2 /) )
      call add_real_4d("cls_out", ns, 2, nzeta, ntheta3, cls, order=(/ 1, 3, 4, 2 /) )
      call add_real_4d("cla_out", ns, 2, nzeta, ntheta3, cla, order=(/ 1, 3, 4, 2 /) )
    else
      call add_null("crs_out")
      call add_null("cra_out")
      call add_null("czs_out")
      call add_null("cza_out")
      call add_null("cls_out")
      call add_null("cla_out")
    end if

    call close_dbg_out()
  end if ! dump_symforce

END SUBROUTINE symforce
