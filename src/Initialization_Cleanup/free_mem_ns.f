!> \file
      SUBROUTINE free_mem_ns(lreset)
      USE vmec_main
      USE realspace
      USE vforces
      USE xstuff
!       USE csplinx
      USE fbal
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      LOGICAL, INTENT(in) :: lreset
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat1 = 0, istat2 = 0, istat3 = 0, istat4 = 0,
     1           istat5 = 0, istat6 = 0, istat7 = 0, istat8 = 0,
     2           istat9 = 0, istat10 = 0
C-----------------------------------------------
      IF (ALLOCATED(phip))
     1  DEALLOCATE (phip, chip, shalf, sqrts, wint, stat=istat3)

      IF (ALLOCATED(ireflect))
     1  DEALLOCATE (ireflect, stat=istat4)

      IF (ALLOCATED(ard))
     1  DEALLOCATE (ard,arm,brd,brm,crd,azd,azm,bzd,bzm, sm,sp,
     2        bmin, bmax,stat=istat7)

      IF (ALLOCATED(iotaf))
     1  DEALLOCATE (iotaf,mass,phi,presf,jcuru,jcurv,jdotb,buco,bvco,
     1     bucof, bvcof, chi,
     2     bdotgradv,equif,specw,tcon,psi,yellip,yinden,
     3     ytrian,yshift,ygeo,overr,faclam,iotas,phips,chips,pres,vp,
     4     beta_vol, jperp2, jpar2, bdotb, clam, blam, dlam, phipf,
     5     chipf, rru_fac, rzu_fac, frcc_fac, fzsc_fac, icurv, vpphi,
     6     presgrad, r01, z01, bdamp, stat=istat8)

      IF (ALLOCATED(gc))
     1  DEALLOCATE (gc, xsave, xstore, xcdot, stat=istat10)
      IF (ALLOCATED(xc) .and. lreset) DEALLOCATE (xc, scalxc)

      IF (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0 .or.
     1      istat4.ne.0 .or. istat5.ne.0 .or. istat6.ne.0 .or.
     2      istat7.ne.0 .or. istat8.ne.0 .or. istat9.ne.0 .or.
     3      istat10.ne.0) THEN
          PRINT *,' deallocation problem in free_mem_ns'
          PRINT *,' istat1 = ',istat1,' istat2 = ',istat2
          PRINT *,' istat3 = ',istat3,' istat4 = ',istat4
          PRINT *,' istat5 = ',istat5,' istat6 = ',istat6
          PRINT *,' istat7 = ',istat7,' istat8 = ',istat8
          PRINT *,' istat9 = ',istat9,' istat10= ',istat10
       ENDIF

      END SUBROUTINE free_mem_ns
