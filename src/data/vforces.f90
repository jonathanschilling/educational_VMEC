!> \file
MODULE vforces
  USE stel_kinds, ONLY: rprec
  IMPLICIT NONE

  REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: &
     armn, azmn, brmn, bzmn, crmn, czmn, blmn, clmn
  REAL(rprec), POINTER, DIMENSION(:) :: &
      armn_e, armn_o, azmn_e, azmn_o, &
      brmn_e, brmn_o, bzmn_e, bzmn_o, &
      crmn_e, crmn_o, czmn_e, czmn_o, &
      blmn_e, blmn_o, clmn_e, clmn_o

END MODULE vforces
