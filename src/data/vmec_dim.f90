!> \file
MODULE vmec_dim

  IMPLICIT NONE

  INTEGER :: mpol1 !< mpol-1
  INTEGER :: ntor1 !< ntor+1
  INTEGER :: mnmax
  INTEGER :: ntheta1
  INTEGER :: ntheta2
  INTEGER :: ntheta3 !< effective number of poloidal grid points
  INTEGER :: nznt
  INTEGER :: nrzt
  INTEGER :: mns
  INTEGER :: mnsize
  INTEGER :: mnmax_nyq
  INTEGER :: ns !< number of flux surfaces
  INTEGER :: ns1 !< ns-1
  INTEGER :: ns_maxval

END MODULE vmec_dim
