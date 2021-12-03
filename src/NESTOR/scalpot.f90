!> \file
!> \brief Compute all required terms for solving for the scalar magnetic potential

!> \brief Compute all required terms for solving for the scalar magnetic potential
!>
!> @param bvec
!> @param amatrix
!> @param wint
!> @param ivacskip
!> @param lasym
!> @param m_map
!> @param n_map
SUBROUTINE scalpot(bvec, amatrix, wint, ivacskip, lasym, m_map, n_map)
   USE vacmod, vm_amatrix => amatrix
   use dbgout
   IMPLICIT NONE

   INTEGER, INTENT(in) :: ivacskip
   REAL(rprec), INTENT(out) :: bvec(mnpd2), amatrix(mnpd2*mnpd2), m_map(mnpd2), n_map(mnpd2)
   REAL(rprec), dimension(nuv2), INTENT(in) :: wint
   logical, intent(in) :: lasym

   INTEGER :: ip, istat

   IF (.not.ALLOCATED(amatsav)) then
      STOP 'AMATSAV not allocated in scalpot'
   end if

   ! COMPUTE TRANFORM OF ANALYTIC SOURCE AND KERNEL.
   ! ON EXIT:
   ! BVEC  CONTAINS THE TRANSFORM OF THE ANALYTIC SOURCE AND
   ! GRPMN CONTAINS THE TRANSFORM OF THE NORMAL DERIVATIVE OF THE GREENS FUNCTION [PKM, EQ.(2.15)]
   ! GReen's function Primed (normal derivative...) and Fourier-transformed to MN mode numbers --> "GR P MN"
   CALL analyt (grpmn, bvec, ivacskip, lasym, m_map, n_map, grpmn_m_map_wrt, grpmn_n_map_wrt)

   IF (ivacskip .ne. 0) THEN
      ! FOR ivacskip != 0, USE PREVIOUSLY COMPUTED bvecsav FOR SPEED

      ! Here, bvecsav contains the previous non-singular contribution to the "final" bvec from fouri.
      ! For ivacskip != 0, this contribution is used from the cache in bvecsav.
      bvec = bvec + bvecsav
   ELSE
      ! Here, bvecsav stores the singular part of bvec from analyt
      ! to be subtracted from the "final" bvec computed below by fouri
      bvecsav = bvec

      ! COMPUTE SURFACE INTEGRALS OF SOURCE AND GREENS FUNCTION NEEDED
      ! FOR SPECTRAL DECOMPOSITION OF POTENTIAL INTEGRAL EQUATION
      ! NOTE: SOURCE IS THE RHS OF EQ.(3.2), KERNEL IS THE LHS OF EQ (3.2).
      ! IP IS THE INDEX OF THE PRIMED VARIABLE MESH.
      gstore = 0
      DO ip = 1, nuv2

         ! COMPUTE DIFFERENCE BETWEEN THE EXACT AND ANALYTIC GREENS FUNCTION AND GRADIENT
         ! [FIRST TERMS IN EQ.(2.14, 2.16)].
         CALL greenf (green(1,ip), greenp(1,ip), ip)

         ! PERFORM INTEGRAL (SUM) OVER PRIMED MESH OF NON-SINGULAR SOURCE TERM
         ! [(h-hsing)(u,v,u',v') == bexni(ip)*green(u,v; ip) in Eq. 2.16]
         ! AND STORE IT - FOR UNPRIMED MESH VALUES - IN GSTORE
         gstore = gstore + bexni(ip)*green(:,ip)

      END DO

      if (open_dbg_context("vac1n_greenf", id=icall)) then

        call add_real_4d("green",  nv, nu, nv, nu3, green)
        call add_real_4d("greenp", nv, nu, nv, nu3, greenp)

        call add_real_4d("gstore", nv, nu, nv, nu3, gstore)

        call close_dbg_out()
      end if

      ! PERFORM FOURIER INTEGRAL OF GRADIENT KERNEL (GREENP) OVER THE UNPRIMED MESH
      ! AND STORE IN GRPMN (NOTE THAT GRPMN IS ADDED TO THE ANALYTIC PIECE IN EQ. 2.14,
      ! - COMPUTED IN ANALYT - WHICH HAS THE APPROPRIATE SIN, COS FACTORS ALREADY)
      CALL fourp (grpmn, greenp)

      ! COMPUTE FOURIER INTEGRAL OF GRADIENT (GRPMN) OVER PRIMED MESH IN EQ. 2.14
      ! AND SOURCE (GSTORE) OVER UNPRIMED MESH IN EQ. 2.16
      CALL fouri (grpmn, gstore, amatrix, amatsav, bvec, wint, lasym)

      ! debugging: focus on Fourier transforms in fouri for now
      ! return

      ! SAVE NON-SINGULAR CONTRIBUTION TO BVEC (IN BVECSAV)
      bvecsav(:mnpd2) = bvec - bvecsav(:mnpd2)

   ENDIF

   amatrix = amatsav

END SUBROUTINE scalpot
