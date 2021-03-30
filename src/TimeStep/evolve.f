!> \file
      SUBROUTINE evolve(time_step, ier_flag, liter_flag, lscreen)
      USE vmec_main
      USE vmec_params, ONLY: bad_jacobian_flag, successful_term_flag,
     1                       norm_term_flag
      USE xstuff
      USE timer_sub
!       USE gmres_mod
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec) :: time_step          !, r0dot
      INTEGER, INTENT(inout) :: ier_flag
      LOGICAL, INTENT(inout) :: liter_flag
      LOGICAL, INTENT(in)  :: lscreen
C-----------------------------------------------``
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: fcn_message =
     1 "External calls to FUNCT3D: "
!      REAL(rprec), PARAMETER :: r0dot_threshold = 5.E-06_dp
      REAL(rprec) :: fsq1, dtau, b1, bprec, fac
      INTEGER :: lcount
      INTEGER, SAVE :: iter_on
C-----------------------------------------------

      IF (iter2 .lt. 10) THEN
         iter_on = -1
      end if

!
!     COMPUTE MHD FORCES
!     MUST CALL funct3d EVEN WHEN IN 2D PRECONDITIONING MODE, SINCE
!     INITIAL RESIDUALS MUST BE KNOWN WHEN CALLING gmres_fun, etc.
!
      CALL funct3d (lscreen, ier_flag)

!
!     COMPUTE ABSOLUTE STOPPING CRITERION
      IF (iter2.eq.1 .and. irst.eq.2) THEN
         ier_flag = bad_jacobian_flag

      ELSE IF (fsqr.le.ftolv .and. fsqz.le.ftolv .and.
     1         fsql.le.ftolv) THEN
         liter_flag = .false.
         ier_flag = successful_term_flag
      ENDIF

      IF (ier_flag.ne.norm_term_flag .or. .not.liter_flag) RETURN


!     COMPUTE DAMPING PARAMETER (DTAU) AND EVOLVE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE

      fsq1 = fsqr1 + fsqz1 + fsql1

      IF (iter2 .eq. iter1) otau(:ndamp) = cp15/time_step

      bprec = 1

      IF (iter2.gt.iter1 .and. fsq1.ne.zero)
     1    dtau = MIN(ABS(LOG(fsq1/fsq)), bprec*cp15)

      fsq = fsq1

      otau(1:ndamp-1) = otau(2:ndamp)

      IF (iter2 .gt. iter1) otau(ndamp) = dtau/time_step
!REMOVED 071505: OTHERWISE I=1 STATE REPEATED (SKIP THIS TO GET OUT OF ITER2=1 STATE)
!     IF (iter2 .le. 1) RETURN

      otav = SUM(otau(:ndamp))/ndamp
      dtau = time_step*otav/2

      b1  = one - dtau
      fac = one/(one + dtau)

!
!     THIS IS THE TIME-STEP ALGORITHM. IT IS ESSENTIALLY A CONJUGATE
!     GRADIENT METHOD, WITHOUT THE LINE SEARCHES (FLETCHER-REEVES),
!     BASED ON A METHOD GIVEN BY P. GARABEDIAN
!
      xcdot = fac*(b1*xcdot + time_step*gc)
      xc    = xc + time_step*xcdot

      END SUBROUTINE evolve
