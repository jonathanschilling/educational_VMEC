!> \file
!> \brief compute aspect-ratio (independent of elongation): \f$A = <R>/\sqrt{<a b>}\f$

!> \brief compute aspect-ratio (independent of elongation): \f$A = <R>/\sqrt{<a b>}\f$
!>        where \f$\pi <a>^2 = \textrm{Area (toroidally averaged)}\f$
!>        and   \f$2 \pi <R> \textrm{Area} = \textrm{Volume}\f$
FUNCTION aspectratio()
  use stel_constants, only: pi
  USE vmec_main
  USE realspace
  USE vmec_io
  IMPLICIT NONE

  INTEGER :: lk, l
  REAL(rprec) :: rb, zub, t1, aspectratio ! pi,

 ! routine for computing aspect-ratio (independent of elongation):
 ! A = <R>/<ab>**.5
 !
 ! WHERE pi <a>**2 = Area (toroidally averaged)
 !       2*pi * <R> * Area = Volume
 ! Use integration by parts to compute as surface integral (Stoke''s theorem)

  !pi = 4*ATAN(one) ! now from stel_constants

  ! Compute Volume and Mean (toroidally averaged) Cross Section Area

  volume_p = 0
  cross_area_p = 0
  DO lk = 1, nznt
     l = ns*lk
     rb  = r1(l,0) + r1(l,1)
     zub = zu(l,0) + zu(l,1)
     t1  = rb*zub*wint(l)
     volume_p = volume_p + rb*t1
     cross_area_p = cross_area_p + t1
  END DO

  volume_p = 2*pi*pi*ABS(volume_p)
  cross_area_p = 2*pi*ABS(cross_area_p)

  Rmajor_p = volume_p/(2*pi*cross_area_p)
  Aminor_p = SQRT(cross_area_p/pi)

  aspectratio = Rmajor_p/Aminor_p

END FUNCTION aspectratio
