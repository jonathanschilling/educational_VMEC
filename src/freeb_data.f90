!> \file
!> \brief Write out edge values of fields.

!> \brief Write out edge values of fields.
!>
!> @param rmnc stellarator-symmetric Fourier coefficients of \f$R\f$
!> @param zmns stellarator-symmetric Fourier coefficients of \f$Z\f$
!> @param rmns non-stellarator-symmetric Fourier coefficients of \f$R\f$
!> @param zmnc non-stellarator-symmetric Fourier coefficients of \f$Z\f$
!> @param bmodmn stellarator-symmetric Fourier coefficients of \f$|\mathbf{B}|\f$
!> @param bmodmn1 non-stellarator-symmetric Fourier coefficients of \f$|\mathbf{B}|\f$
SUBROUTINE freeb_data (rmnc, zmns, rmns, zmnc, bmodmn, bmodmn1)
  USE vmec_main
  USE vacmod, only: brv, bphiv, bzv, bsqvac, potvac, mnpd, xmpot, xnpot
  USE realspace, ONLY: r1, z1
  IMPLICIT NONE

  REAL(rprec), DIMENSION(mnmax) :: rmnc, zmns, rmns, zmnc, bmodmn, bmodmn1

  INTEGER :: iprint, nzskip, l, k, lk, mn,                           &
             mn0, n, nedge, iu, iv, nl, lkr
  REAL(rprec) :: zeta, potsin, potcos
  REAL(rprec), ALLOCATABLE, DIMENSION(:) :: rb, phib, zb

  ! WRITE OUT EDGE VALUES OF FIELDS TO FORT.NEDGE0 (INCLUDE REFLECTED POINT)
  ! NOTE: BR, BPHI, BZ WERE COMPUTED IN BSS, CALLED FROM EQFOR
  IF (ivac.le.0 .or. .not.lfreeb) RETURN

  ALLOCATE (rb(2*nznt), phib(2*nznt), zb(2*nznt), stat=l)
  IF (l .ne. 0) STOP 'allocation error in freeb_data'

  nedge = 0
  lkr = nznt
  DO iv = 1,nzeta
     zeta = (twopi*(iv-1))/(nzeta*nfp)
     DO iu = 1,ntheta3
        lk = iv + nzeta*(iu-1)
        nl = ns*lk
        nedge = nedge+1

        rb(lk)   = r1(nl,0) + r1(nl,1)
        phib(lk) = zeta
        zb(lk)   = z1(nl,0) + z1(nl,1)

        ! INCLUDE -u,-v POINTS HERE BY STELLARATOR SYMMETRY
        IF (.not.lasym .and. (iu.ne.1 .and. iu.ne.ntheta2)) THEN
          lkr = lkr + 1
          nedge = nedge+1

          rb(lkr)   = rb(lk)
          phib(lkr) =-phib(lk)
          zb(lkr)   =-zb(lk)

          bredge(lkr) = -bredge(lk)
          bpedge(lkr) =  bpedge(lk)
          bzedge(lkr) =  bzedge(lk)
        ENDIF
     END DO
  END DO

  ! WRITE OUT (TO THREED1 FILE) VACUUM INFORMATION
  IF (.not.lfreeb) THEN ! TODO: should be handled by check for (... .or. .not.lfreeb) above already
     DEALLOCATE (rb, phib, zb, stat=l)
     RETURN
  END IF

  ! TODO: below outputs only up to ntheta2, so why compute full theta range above?

  nzskip = 1 + nzeta/6

  ! iprint == 1
  WRITE (nthreed, 750)
750 FORMAT(/,3x,'NF*PHI',7x,' Rb ',8x,' Zb ',&
             6x,'BSQMHDI',5x,'BSQVACI',      &
             5x,'BSQMHDF',5x,'BSQVACF',/)
  DO l = 1, nzeta, nzskip
    zeta = (360.0_dp*(l - 1))/nzeta

    DO k = 1, ntheta2
      lk = l + nzeta*(k - 1)
      WRITE (nthreed, 770) zeta, rb(lk), zb(lk), &
                           (bsqsav(lk,n),n=1,3), bsqvac(lk)
770 FORMAT(1p,e10.2,6e12.4)
    END DO
  end do

  ! iprint == 2
  WRITE (nthreed, 760)
760 FORMAT(/,3x,'NF*PHI',7x,' Rb ',8x,' Zb ',&
             6x,'BR', 8x,'BPHI', 6x,'BZ',  &
             8x,'BRv',7x,'BPHIv',5x,'BZv',/)
  DO l = 1, nzeta, nzskip
    zeta = (360.0_dp*(l - 1))/nzeta
    DO k = 1, ntheta2
      lk = l + nzeta*(k - 1)
      WRITE (nthreed, 780) zeta,       rb(lk),     zb(lk),      &
                           bredge(lk), bpedge(lk), bzedge(lk),  &
                           brv(lk),    bphiv(lk),  bzv(lk)
780 FORMAT(1p,e10.2,2e12.4,6e10.2)
    END DO
  end do

  ! allocated in eqfor
  DEALLOCATE (rb, phib, zb, bredge, bpedge, bzedge, stat=l)

  ! DIAGNO v1 input ???
  IF (lasym) THEN
     WRITE (nthreed, 900)
900 FORMAT(//,3x,'nb',2x,'mb',&
              6x,'rbc',9x,'zbs',9x,'rbs',9x,'zbc',      &
              6x,'vacpot_s',   4x,'vacpot_c', &
              2x,'|B|_c(s=.5)',1x,'|B|_c(s=1.)'/)
     DO mn = 1, mnmax
        potsin = 0
        potcos = 0
        DO mn0 = 1, mnpd
           IF ( (NINT(xnpot(mn0)).eq.NINT(xn(mn))) .and.                &
                (NINT(xmpot(mn0)).eq.NINT(xm(mn))) ) THEN
              potsin = potvac(mn0)
              potcos = potvac(mn0+mnpd)
              EXIT
           END IF
        END DO
     WRITE (nthreed, 910) NINT(xn(mn)/nfp), NINT(xm(mn)), &
        rmnc(mn), zmns(mn), rmns(mn), zmnc(mn),           &
        potsin,     potcos,                               &
        bmodmn(mn), bmodmn1(mn)
910 FORMAT(i5,i4,1p,10e12.4) ! TODO: only 8 real?
     END DO

  ELSE
     WRITE (nthreed, 800)
800 FORMAT(//,3x,'nb',2x,'mb',  &
              6x,'rbc',9x,'zbs',&
              6x,'vacpot_s',    &
              2x,'|B|_c(s=.5)',1x,'|B|_c(s=1.)'/)
     DO mn = 1, mnmax
        potsin = 0
        DO mn0 = 1, mnpd
           IF ( (NINT(xnpot(mn0)).eq.NINT(xn(mn))) .and.                &
                (NINT(xmpot(mn0)).eq.NINT(xm(mn))) ) THEN
              potsin = potvac(mn0)
              EXIT
           END IF
        END DO
        WRITE (nthreed, 810) NINT(xn(mn)/nfp), NINT(xm(mn)),            &
            rmnc(mn), zmns(mn), &
            potsin,             &
            bmodmn(mn), bmodmn1(mn)
810 FORMAT(i5,i4,1p,7e12.4) ! TODO: only 5 REAL ???
     END DO
  END IF

  WRITE (nthreed, *)

END SUBROUTINE freeb_data
