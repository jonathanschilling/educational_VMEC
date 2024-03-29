!> \file
!> \brief Print iteration progress to screen and \c threed1 output file.

!> \brief Print iteration progress to screen and \c threed1 output file.
!>
!> @param i0 current iteration number (iter2)
!> @param delt0 current time step
!> @param w0 current MHD energy
SUBROUTINE printout(i0, delt0, w0)
  USE vmec_main
  USE realspace
  USE xstuff
  use vmec_params, only: ntmax

  use dbgout

  IMPLICIT NONE

  INTEGER :: i0
  REAL(rprec) :: delt0, w0

! #ifndef _HBANGLE
  CHARACTER(LEN=*), PARAMETER :: iter_line  = "  ITER    FSQR      FSQZ      CONR      CONZ      MHDR      MHDZ      FSQL   "
  CHARACTER(LEN=*), PARAMETER :: fsq_line   = "   fsqr      fsqz      conr      conz      mhdr      mhdz      fsql      DELT    "
  CHARACTER(LEN=*), PARAMETER :: iter_lines = iter_line
  CHARACTER(LEN=*), PARAMETER :: fsq_lines  = fsq_line
  CHARACTER(LEN=*), PARAMETER :: raxis_line = "RAX(v=0) "
! #end /* ndef _HBANGLE */

  CHARACTER(LEN=*), PARAMETER :: delt_line  = "    DELT  "
  CHARACTER(LEN=*), PARAMETER :: zaxis_line = "  ZAX(v=0)      "

  REAL(rprec) :: betav, w, avm, den
  CHARACTER(len=LEN(iter_line) + LEN(fsq_line) + LEN(raxis_line) + LEN(zaxis_line)) :: print_line
  logical :: dbgout_printout

  betav = wp/wb
  w = w0*twopi*twopi

  den = zero ! TODO: why? will be set of sum(vp(2:ns)) below anyway...
  specw(1) = one
  gc = xstore ! TODO: why compute spectral width from backup and not current gc (== physical xc) --> <M> includes scalxc ???

  dbgout_printout = open_dbg_context("printout", num_eqsolve_retries)
  if (dbgout_printout) then
    ! dump gc before it gets modified by spectrum() below
    call add_real_5d("gc", 3, ntmax, ns, ntor1, mpol, gc, order=(/ 3, 4, 5, 2, 1 /) )
  end if ! dbgout_printout

  CALL spectrum (gc(:irzloff), gc(1+irzloff:2*irzloff))
  den = SUM(vp(2:ns))
  avm = DOT_PRODUCT(vp(2:ns), specw(2:ns)+specw(1:ns-1))
  avm = 0.5_dp*avm/den ! volume-averaged spectral width (_av_erage _M_)

  delbsq = 0.0_dp ! default output in case of fixed-boundary run
  IF (ivac .ge. 1 .and. iter2.gt.1) then
     delbsq = SUM( dbsq(:nznt)*wint(2:nrzt:ns) ) / SUM( bsqsav(:nznt,3)*wint(2:nrzt:ns) )
  end if

  if (dbgout_printout) then
    call add_real("betav", betav)
    call add_real("avm", avm)
    call add_real("delbsq", delbsq)

    call add_real_1d("specw", ns, specw)

    call close_dbg_out()
  end if ! printout

  IF (i0.eq.1 .and. lfreeb) THEN
     print_line = iter_lines // " " // raxis_line
     IF (lasym) print_line = TRIM(print_line) // " " // zaxis_line
     PRINT 20, TRIM(print_line)//delt_line
     print_line = iter_line // fsq_line // raxis_line
     IF (lasym) print_line = TRIM(print_line) // " " // zaxis_line
     WRITE (nthreed, 16) TRIM(print_line)
  ELSE IF (i0.eq.1 .and. .not.lfreeb) THEN
     print_line = raxis_line
     IF (lasym) print_line = raxis_line // zaxis_line
     PRINT 30, iter_lines, TRIM(print_line)//delt_line
     print_line = iter_line // fsq_line // raxis_line // "     "
     IF (lasym) then
        print_line = iter_line // fsq_line // raxis_line // zaxis_line
     end if
     WRITE (nthreed, 25) TRIM(print_line)
  ENDIF
16 FORMAT(/,a,6x,'WMHD      BETA     PHIEDGE  DEL-BSQ    FEDGE',/)
20 FORMAT(/,a,6x,'WMHD      DEL-BSQ',/)
25 FORMAT(/,a,6x,'WMHD      BETA      <M>        ',/)
30 FORMAT(/,a,1x,a,5x,'WMHD',/)

  IF (.not. lasym) THEN
     IF (.not.lfreeb) THEN
        PRINT 45,           i0, fsqr, fsqz, fsqr_con, fsqz_con, fsqr_mhd, fsqz_mhd, &
           fsql, r00, delt0, w
        WRITE (nthreed, 40) i0, fsqr, fsqz, fsqr_con, fsqz_con, fsqr_mhd, fsqz_mhd, &
           fsql, fsqr1, fsqz1, fsqr1_con, fsqz1_con, fsqr1_mhd, fsqz1_mhd, fsql1, &
           delt0, r00, w, betav, avm
        RETURN
     ENDIF
     PRINT 50,           i0, fsqr, fsqz, fsqr_con, fsqz_con, fsqr_mhd, fsqz_mhd, &
        fsql, r00, delt0, w, delbsq
     WRITE (nthreed, 42) i0, fsqr, fsqz, fsqr_con, fsqz_con, fsqr_mhd, fsqz_mhd, &
        fsql, fsqr1, fsqz1, fsqr1_con, fsqz1_con, fsqr1_mhd, fsqz1_mhd, &
        fsql1, delt0, r00, w, betav, ABS(phiedge), delbsq, fedge

  ELSE ! (.not. lasym)
     IF (.not.lfreeb) THEN
        PRINT 65,           i0, fsqr, fsqz, fsqr_con, fsqz_con, fsqr_mhd, fsqz_mhd, &
           fsql, r00, z00, delt0, w
        WRITE (nthreed, 60) i0, fsqr, fsqz, fsqr_con, fsqz_con, fsqr_mhd, fsqz_mhd, &
           fsql, fsqr1, fsqz1, fsqr1_con, fsqz1_con, fsqr1_mhd, fsqz1_mhd, &
           fsql1, delt0, r00, z00, w, betav, avm
        RETURN
     ENDIF
     PRINT 70,           i0, fsqr, fsqz, fsqr_con, fsqz_con, fsqr_mhd, fsqz_mhd, &
        fsql, r00, z00, delt0, w, delbsq
     WRITE (nthreed, 60) i0, fsqr, fsqz, fsqr_con, fsqz_con, fsqr_mhd, fsqz_mhd, &
        fsql, fsqr1, fsqz1, fsqr1_con, fsqz1_con, fsqr1_mhd, fsqz1_mhd, &
        fsql1, delt0, r00, z00, w, betav, ABS(phiedge), delbsq, fedge
  END IF

40 FORMAT(i6,1x,1p,15e10.2,e11.3,e12.4,e11.3,0p,f7.3,1p,2e9.2)
42 FORMAT(i5,1p,15e10.2,e11.3,e12.4,2e11.3,0p,f7.3,1p,e9.2)
45 FORMAT(i5,1p,7e10.2,e11.3,e10.2,e12.4)
50 FORMAT(i5,1p,7e10.2,e11.3,e10.2,e12.4,e11.3)
60 FORMAT(i6,1x,1p,15e10.2,2e11.3,e12.4,e11.3,0p,f7.3,1p,2e9.2)
65 FORMAT(i5,1p,7e10.2,2e11.3,e10.2,e12.4)
70 FORMAT(i5,1p,7e10.2,2e11.3,e10.2,e12.4,e11.3)

END SUBROUTINE printout
