!> \file
!> \brief Symmetrize some quantities so that they can be output (?)

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
