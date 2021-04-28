!> \file
!> This module containes functions used by the profiles.

      MODULE functions

      USE stel_kinds

      implicit none

      PUBLIC :: two_power, two_power_gs
      PRIVATE :: check

      CONTAINS

!> \brief Profile function for the \c two\_power profile.
!>        \f$b(0) * (1 - x^{b(1)})^{b(2)}\f$
!>
!> @param x evaluation location
!> @param b parameter vector
      REAL(rprec) FUNCTION two_power(x, b)

      REAL(rprec), INTENT(in)                  :: x
      REAL(rprec), DIMENSION(0:20), INTENT(in) :: b

      two_power = b(0)*((1 - x**b(1))**b(2))

      END FUNCTION

!> \brief Profile function for the \c two\_power\_gs profile.
!>        \f$\texttt{two\_power}(x)*(1 + \sum\left[b(i)*\exp(-(x - b(i+1))/b(i+2))^2\right])\f$
!>
!> @param x evaluation location
!> @param b parameter vector
      REAL(rprec) FUNCTION two_power_gs(x, b)

      REAL(rprec), INTENT(in)                  :: x
      REAL(rprec), DIMENSION(0:20), INTENT(in) :: b

      INTEGER :: i

      two_power_gs = 1.0
      DO i = 3, 18, 3
         two_power_gs = two_power_gs +                                         &
     &                  b(i)*exp(-((x - b(i+1))/b(i+2))**2)
      END DO
      two_power_gs = two_power_gs*two_power(x,b)

      END FUNCTION

!> \brief Main test function
      FUNCTION function_test()

      LOGICAL :: function_test
      REAL(rprec) :: result
      REAL(rprec), DIMENSION(0:20) :: b(0:20) = 0

!> Test \c two\_power function for x = 0, b = {1,10,2} is 1
      b(0:2) = (/ 1.0d+0, 10.0d+0, 2.0d+0 /)
      result = two_power(0.0d+0, b)
      function_test = check(1.0d+0,result,1,"two_power()")
      IF (.not.function_test) RETURN

!> Test \c two\_power function for x = 1, b = {1,10,2} is 0
      result = two_power(1.0d+0, b)
      function_test = check(0.0d+0,result,2,"two_power()")
      IF (.not.function_test) RETURN

!> Test \c two\_power function for x = 0.5, b = {1,1,1} is 0.5
      b(0:2) = (/ 1.0d+0, 1.0d+0, 1.0d+0 /)
      result = two_power(0.5d+0, b)
      function_test = check(0.5d+0,result,3,"two_power()")
      IF (.not.function_test) RETURN

!> Test \c two\_power function for x = 0.5, b = {1,1,2} is 0.25
      b(0:2) = (/ 1.0d+0, 1.0d+0, 2.0d+0 /)
      result = two_power(0.5d+0, b)
      function_test = check(0.25d+0,result,4,"two_power()")
      IF (.not.function_test) RETURN

!> Test \c two\_power_gs function for x = 0.4, b = {1,1,1,0,0,1} is two\_power(x,b)
      b(0:5) = (/ 1.0d+0, 1.0d+0, 1.0d+0, 0.0d+0, 0.0d+0, 1.0d+0 /)
      result = two_power_gs(0.4d+0, b)
      function_test = check(two_power(0.4d+0, b),                              &
     &                      result,1,"two_power_gs")
      IF (.not.function_test) RETURN

!> Test \c two\_power\_gs function for x = 0.8, b = {1,1,0,1,0.8,0.1} is 2
      b(0:5) = (/ 1.0d+0, 1.0d+0, 0.0d+0, 1.0d+0, 0.8d+0, 0.1d+0 /)
      result = two_power_gs(0.8d+0, b)
      function_test = check(2.0d+0,result,1,"two_power_gs")
      IF (.not.function_test) RETURN

      END FUNCTION

!> \brief Check Test result
      FUNCTION check(expected, recieved, testNum, name)

      LOGICAL :: check
      REAL(rprec), INTENT(in) :: expected, recieved
      INTEGER, INTENT(in) :: testNum
      CHARACTER (LEN=*), INTENT(in) :: name

      check = expected .eq. recieved
      IF (.not.check) THEN
         write(*,*) "functions.f: ", name, " test", testNum,                   &
     &              "failed."
         write(*,*) "Expected", expected, "Recieved", recieved
      END IF

      END FUNCTION

      END MODULE
