!> \file
SUBROUTINE printout(i0, delt0, w0, lscreen)
  USE vmec_main
  USE realspace
  USE xstuff
  IMPLICIT NONE

  INTEGER :: i0
  REAL(rprec) :: delt0, w0
  LOGICAL :: lscreen

  CHARACTER(LEN=*), PARAMETER :: iter_line  = "  ITER    FSQR      FSQZ      FSQL   "
  CHARACTER(LEN=*), PARAMETER :: fsq_line   = "   fsqr      fsqz      fsql      DELT    "
  CHARACTER(LEN=*), PARAMETER :: iter_lines = iter_line
  CHARACTER(LEN=*), PARAMETER :: fsq_lines  = fsq_line
  CHARACTER(LEN=*), PARAMETER :: raxis_line = "RAX(v=0) "
  CHARACTER(LEN=*), PARAMETER :: delt_line  = "    DELT  "
  CHARACTER(LEN=*), PARAMETER :: zaxis_line = "  ZAX(v=0)      "

  REAL(rprec) :: betav, w, avm, den
  CHARACTER(len=LEN(iter_line) + LEN(fsq_line) + LEN(raxis_line) + LEN(zaxis_line)) :: print_line

  betav = wp/wb
  w = w0*twopi*twopi
  den = zero
  specw(1) = one

  gc = xstore
  CALL spectrum (gc(:irzloff), gc(1+irzloff:2*irzloff))

  den = SUM(vp(2:ns))
  avm = DOT_PRODUCT(vp(2:ns), specw(2:ns)+specw(1:ns-1))
  avm = 0.5_dp*avm/den
  IF (ivac .ge. 1 .and. iter2.gt.1) then
     delbsq = SUM(dbsq(:nznt)*wint(2:nrzt:ns))/SUM(bsqsav(:nznt,3)*wint(2:nrzt:ns))
  end if
  IF (i0.eq.1 .and. lfreeb) THEN
     print_line = iter_lines // " " // raxis_line
     IF (lasym) print_line = TRIM(print_line) // " " // zaxis_line
     IF (lscreen) PRINT 20, TRIM(print_line)//delt_line
     print_line = iter_line // fsq_line // raxis_line
     IF (lasym) print_line = TRIM(print_line) // " " // zaxis_line
     WRITE (nthreed, 16) TRIM(print_line)
  ELSE IF (i0.eq.1 .and. .not.lfreeb) THEN
     print_line = raxis_line
     IF (lasym) print_line = raxis_line // zaxis_line
     IF (lscreen) PRINT 30, iter_lines, TRIM(print_line)//delt_line
     print_line = iter_line // fsq_line // raxis_line // "     "
     IF (lasym) then
        print_line = iter_line // fsq_line // raxis_line // zaxis_line
     end if
     WRITE (nthreed, 25) TRIM(print_line)
  ENDIF
15 FORMAT(/,a,6x,'WMHD      BETA      <M>   DEL-BSQ   FEDGE',/)
16 FORMAT(/,a,6x,'WMHD      BETA     PHIEDGE  DEL-BSQ    FEDGE',/)
20 FORMAT(/,a,6x,'WMHD      DEL-BSQ',/)
25 FORMAT(/,a,6x,'WMHD      BETA      <M>        ',/)
30 FORMAT(/,a,1x,a,5x,'WMHD',/)

  IF (.not. lasym) THEN
     IF (.not.lfreeb) THEN
        IF (lscreen) PRINT 45, i0, fsqr, fsqz, fsql, r00, delt0, w
        WRITE (nthreed, 40) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, fsql1, &
           delt0, r00, w, betav, avm
     RETURN
     ENDIF
     IF (lscreen) PRINT 50, i0, fsqr, fsqz, fsql, r00, delt0, w, delbsq
     WRITE (nthreed, 42) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, &
        fsql1, delt0, r00, w, betav, ABS(phiedge), delbsq, fedge

  ELSE ! (.not. lasym)
     IF (.not.lfreeb) THEN
        IF (lscreen) PRINT 65, i0, fsqr, fsqz, fsql, r00, z00, delt0, w
        WRITE (nthreed, 60) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, &
           fsql1, delt0, r00, z00, w, betav, avm
     RETURN
     ENDIF
     IF (lscreen) PRINT 70, i0, fsqr, fsqz, fsql, r00, z00, delt0, w, delbsq
     WRITE (nthreed, 60) i0, fsqr, fsqz, fsql, fsqr1, fsqz1, &
        fsql1, delt0, r00, z00, w, betav, ABS(phiedge), delbsq, fedge
  END IF

40 FORMAT(i6,1x,1p,7e10.2,e11.3,e12.4,e11.3,0p,f7.3,1p,2e9.2)
42 FORMAT(i5,1p,7e10.2,e11.3,e12.4,2e11.3,0p,f7.3,1p,e9.2)
45 FORMAT(i5,1p,3e10.2,e11.3,e10.2,e12.4)
50 FORMAT(i5,1p,3e10.2,e11.3,e10.2,e12.4,e11.3)
60 FORMAT(i6,1x,1p,7e10.2,2e11.3,e12.4,e11.3,0p,f7.3,1p,2e9.2)
65 FORMAT(i5,1p,3e10.2,2e11.3,e10.2,e12.4)
70 FORMAT(i5,1p,3e10.2,2e11.3,e10.2,e12.4,e11.3)

END SUBROUTINE printout
