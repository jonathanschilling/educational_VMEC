!> \file
MODULE timer_sub
  USE stel_kinds, ONLY: rprec
  IMPLICIT NONE

  INTEGER, PARAMETER :: tsum = 0, tvac = 1, tread = 2, twout= 3,    &
                        teqf = 4, tfun = 5, trecon= 6, tfft = 7,    &
                        tffi = 8, tfor = 9, tbcov =10, tres = 11,   &
                        tprec2d = 12

  REAL(rprec) :: tvacon, tvacoff, tfunon, tfunoff,                  &
                 twouton, twoutoff, teqfon, teqfoff,                &
                 treadon, treadoff, timeon, timeoff,                &
                 treconon, treconoff, tffton, tfftoff,              &
                 tbcovon, tbcovoff, tforon, tforoff,                &
                 treson, tresoff, tprec2don, tprec2doff,            &
                 timer_tsum, timer_tfun
  REAL(rprec), DIMENSION(0:12) :: timer

CONTAINS

SUBROUTINE write_times (nthreed, lscreen, lfreeb, lprec2d)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: nthreed
  LOGICAL, INTENT(in) :: lscreen, lfreeb, lprec2d
  INTEGER             :: i, nform
  CHARACTER(LEN=*), DIMENSION(0:12), PARAMETER :: form =            &
      (/ 'TOTAL COMPUTATIONAL TIME       ',                         &
         'TIME IN VACUUM LOOP            ',                         &
         'TIME TO READ IN DATA           ',                         &
         'TIME TO WRITE DATA TO WOUT     ',                         &
         'TIME IN EQFORCE                ',                         &
         'TIME (REMAINDER) IN FUNCT3D    ',                         &
         'TIME IN PROFILE RECONSTRUCTION ',                         &
         'TIME IN FOURIER TRANSFORM      ',                         &
         'TIME IN INVERSE FOURIER XFORM  ',                         &
         'TIME IN FORCES + SYMFORCES     ',                         &
         'TIME IN BCOVAR                 ',                         &
         'TIME IN RESIDUE                ',                         &
         'TIME IN PRECON2D SETUP         '                          &
       /)

  timer_tsum = timer(tsum) + timer(twout) + timer(teqf)

  timer_tfun = timer(tfun) - timer(tvac) - timer(trecon)            &
             - timer(tfft) - timer(tffi) - timer(tfor)              &
             - timer(tbcov) - timer(tres)

  DO i = 1,2
     IF (i .eq. 1) nform = 6
     IF (i .eq. 2) nform = nthreed
     IF (.not.lscreen .and. i.eq.1) CYCLE
     IF (lfreeb) THEN
        WRITE (nform, 20)                                           &
              form(tsum),timer_tsum,form(tread),timer(tread),       &
              form(twout),timer(twout),form(teqf),timer(teqf),      &
              form(tvac),timer(tvac),                               &
              form(tfft), timer(tfft), form(tffi), timer(tffi),     &
              form(tfor), timer(tfor), form(tbcov),timer(tbcov),    &
              form(tres), timer(tres)
     ELSE
        WRITE (nform, 20)                                           &
              form(tsum),timer_tsum,form(tread),timer(tread),       &
              form(twout),timer(twout),form(teqf),timer(teqf),      &
              form(tfft), timer(tfft), form(tffi), timer(tffi),     &
              form(tfor), timer(tfor), form(tbcov),timer(tbcov),    &
              form(tres), timer(tres)
     END IF
     IF (lprec2d)  WRITE (nform, 20) form(tprec2d), timer(tprec2d)
     WRITE (nform, 20) form(tfun), timer_tfun
  END DO

20 FORMAT(a35,f12.2,' SECONDS')

END SUBROUTINE write_times

END MODULE timer_sub
