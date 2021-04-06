!> \file
SUBROUTINE free_mem_ns(lreset)
  USE vmec_main
  USE realspace
  USE vforces
  USE xstuff
  USE fbal
  IMPLICIT NONE

  LOGICAL, INTENT(in) :: lreset

  INTEGER :: istat1 = 0, istat2 = 0, istat3 = 0, istat4 = 0, &
             istat5 = 0, istat6 = 0, istat7 = 0, istat8 = 0, &
             istat9 = 0, istat10 = 0

  IF (ALLOCATED(phip)) &
    DEALLOCATE (phip, chip, shalf, sqrts, wint, stat=istat3)

  IF (ALLOCATED(ireflect)) &
    DEALLOCATE (ireflect, stat=istat4)

  IF (ALLOCATED(ard)) &
    DEALLOCATE (ard,arm,brd,brm,crd,azd,azm,bzd,bzm, sm,sp, bmin, bmax,stat=istat7)

  IF (ALLOCATED(iotaf)) &
    DEALLOCATE (iotaf,mass,phi,presf,jcuru,jcurv,jdotb,buco,bvco,   &
       bucof, bvcof, chi,                                           &
       bdotgradv,equif,specw,tcon,psi,yellip,yinden,                &
       ytrian,yshift,ygeo,overr,faclam,iotas,phips,chips,pres,vp,   &
       beta_vol, jperp2, jpar2, bdotb, clam, blam, dlam, phipf,     &
       chipf, rru_fac, rzu_fac, frcc_fac, fzsc_fac, icurv, vpphi,   &
       presgrad, r01, z01, bdamp, stat=istat8)

  IF (ALLOCATED(gc)) then
    DEALLOCATE (gc, xsave, xstore, xcdot, stat=istat10)
  end if

  IF (ALLOCATED(xc) .and. lreset) then
    DEALLOCATE (xc, scalxc)
  end if

  IF (istat1.ne.0 .or. istat2.ne.0 .or. istat3.ne.0 .or. &
      istat4.ne.0 .or. istat5.ne.0 .or. istat6.ne.0 .or. &
      istat7.ne.0 .or. istat8.ne.0 .or. istat9.ne.0 .or. istat10.ne.0) THEN

      PRINT *,' deallocation problem in free_mem_ns'
      PRINT *,' istat1 = ',istat1,' istat2 = ',istat2
      PRINT *,' istat3 = ',istat3,' istat4 = ',istat4
      PRINT *,' istat5 = ',istat5,' istat6 = ',istat6
      PRINT *,' istat7 = ',istat7,' istat8 = ',istat8
      PRINT *,' istat9 = ',istat9,' istat10= ',istat10
  ENDIF

END SUBROUTINE free_mem_ns
