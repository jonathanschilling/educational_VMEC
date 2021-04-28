!> \file
!> \brief Compute total force residual on flux surfaces.

!> \brief Compute total force residual on flux surfaces.
!>
!> @param gcr \f$R\f$-component of force
!> @param gcz \f$Z\f$-component of force
!> @param gnormr normalized total force residual in \f$R\f$
!> @param gnormz normalized total force residual in \f$Z\f$
!> @param gnorm normalization factor for forces
!> @param medge =0: exclude contribution from LCFS; =1: include LCFS contribution
SUBROUTINE getfsq(gcr, gcz, gnormr, gnormz, gnorm, medge)
  USE vmec_main, ONLY: rprec, ns, ns1, mnsize
  USE vmec_params, ONLY: ntmax
  IMPLICIT NONE

  INTEGER, INTENT(in) :: medge
  REAL(rprec), INTENT(out) :: gnormr, gnormz
  REAL(rprec), INTENT(in)  :: gnorm
  REAL(rprec), DIMENSION(ns,mnsize*ntmax), INTENT(in) :: gcr, gcz

  INTEGER :: jsmax

  ! ns1 is ns-1
  ! so  if medge==0, exclude contribution from LCFS
  ! and if medge==1, also take contribution from LCFS into account
  jsmax = ns1 + medge
  gnormr = gnorm * SUM(gcr(:jsmax,:)**2)
  gnormz = gnorm * SUM(gcz(:jsmax,:)**2)

END SUBROUTINE getfsq
