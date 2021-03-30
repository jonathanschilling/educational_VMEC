!---------------------------------------------------------------------!
! adas_bms_splines.f90:                                               !
!     by Michael Kraus (michael.kraus@ipp.mpg.de)                     !
!     created 2009/06/11                                              !
!                                                                     !
! EzSpline wrapper for adas beam stopping data                        !
! takes care of creating, evaluation and freeing of spline objects    !
!                                                                     !
! - spline objects stay in memory after creation                      !
! - saves beam and plasma species for checking if reloading data is   !
!   necessary                                                         !
! - at first adas_initSplines() must be called to allocate all        !
!   arrays with a size equivalent to the number of plasma species     !
!   and create the spline objects                                     !
! - then for each species adas_setSplineData() has to be called to    !
!   feed the splines with data                                        !
! - adas_evalSpline() evaluates the splines                           !
! - after finishing all data requests adas_closeSplines() should be   !
!   called to free splines and deallocate arrays                      !
! - for error codes see adas_bms.f90                                  !
!                                                                     !
!                                                                     !
! References:                                                         !
!                                                                     !
!    1. NTCC PSPLINE Module                                           !
!       http://w3.pppl.gov/ntcc/PSPLINE/                              !
!                                                                     !
!---------------------------------------------------------------------!

module adas_bms_splines

end module
