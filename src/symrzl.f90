!> \file
SUBROUTINE symrzl(r1s, rus, rvs, z1s, zus, zvs, lus, lvs, rcons, &
   zcons, r1a, rua, rva, z1a, zua, zva, lua, lva, rcona, zcona)
  USE vmec_main
  IMPLICIT NONE

  REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(inout) :: &
     r1s, rus, rvs, z1s, zus, zvs, lus, lvs, rcons, zcons
  REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(in)    :: &
     r1a, rua, rva, z1a, zua, zva, lua, lva, rcona, zcona

  INTEGER :: mpar, ir, i, jk, jka, n2

  ! FIRST SUM SYMMETRIC, ANTISYMMETRIC PIECES ON EXTENDED INTERVAL, THETA = [PI,2*PI]
  DO mpar = 0, 1
     DO i = 1 + ntheta2, ntheta1
        ir = ntheta1 + 2 - i                 !-theta
        DO jk = 1, ns*nzeta
           jka = ireflect(jk)                !-zeta
           r1s(jk,i,mpar) = r1s(jka,ir,mpar) - r1a(jka,ir,mpar)
           rus(jk,i,mpar) =-rus(jka,ir,mpar) + rua(jka,ir,mpar)
           z1s(jk,i,mpar) =-z1s(jka,ir,mpar) + z1a(jka,ir,mpar)
           zus(jk,i,mpar) = zus(jka,ir,mpar) - zua(jka,ir,mpar)
           lus(jk,i,mpar) = lus(jka,ir,mpar) - lua(jka,ir,mpar)
           rcons(jk,i,mpar)= rcons(jka,ir,mpar)-rcona(jka,ir,mpar)
           zcons(jk,i,mpar)=-zcons(jka,ir,mpar)+zcona(jka,ir,mpar)
        END DO
        IF (lthreed) THEN
           DO jk = 1, ns*nzeta
              jka = ireflect(jk)
              rvs(jk,i,mpar)=(-rvs(jka,ir,mpar))+rva(jka,ir,mpar)
              zvs(jk,i,mpar) = zvs(jka,ir,mpar) - zva(jka,ir,mpar)
              lvs(jk,i,mpar) = lvs(jka,ir,mpar) - lva(jka,ir,mpar)
           END DO
        ENDIF
     END DO

     ! NOW SUM SYMMETRIC, ANTISYMMETRIC PIECES FOR THETA = [0,PI]
     n2 = ntheta2
     r1s(:,:n2,mpar) = r1s(:,:n2,mpar) + r1a(:,:n2,mpar)
     rus(:,:n2,mpar) = rus(:,:n2,mpar) + rua(:,:n2,mpar)
     z1s(:,:n2,mpar) = z1s(:,:n2,mpar) + z1a(:,:n2,mpar)
     zus(:,:n2,mpar) = zus(:,:n2,mpar) + zua(:,:n2,mpar)
     lus(:,:n2,mpar) = lus(:,:n2,mpar) + lua(:,:n2,mpar)
     rcons(:,:n2,mpar) = rcons(:,:n2,mpar) + rcona(:,:n2,mpar)
     zcons(:,:n2,mpar) = zcons(:,:n2,mpar) + zcona(:,:n2,mpar)
     IF (lthreed) THEN
        rvs(:,:n2,mpar) = rvs(:,:n2,mpar) + rva(:,:n2,mpar)
        zvs(:,:n2,mpar) = zvs(:,:n2,mpar) + zva(:,:n2,mpar)
        lvs(:,:n2,mpar) = lvs(:,:n2,mpar) + lva(:,:n2,mpar)
     ENDIF
  END DO

END SUBROUTINE symrzl
