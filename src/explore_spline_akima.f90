program hello
  implicit none
  INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND(12,100)
  integer, parameter :: dp = rprec

  integer, parameter :: l = 91

  real(rprec), dimension(1:l) :: am_aux_s = (/ 0.0000000E+00_dp, &
     1.8181818E-02_dp,   3.6363636E-02_dp,   5.4545455E-02_dp,   7.2727273E-02_dp,   9.0909091E-02_dp, &
     1.0909091E-01_dp,   1.2727273E-01_dp,   1.4545455E-01_dp,   1.6363636E-01_dp,   1.8181818E-01_dp, &
     2.0000000E-01_dp,   2.1818182E-01_dp,   2.3636364E-01_dp,   2.5454545E-01_dp,   2.7272727E-01_dp, &
     2.9090909E-01_dp,   3.0909091E-01_dp,   3.2727273E-01_dp,   3.4545455E-01_dp,   3.6363636E-01_dp, &
     3.8181818E-01_dp,   4.0000000E-01_dp,   4.1818182E-01_dp,   4.3636364E-01_dp,   4.5454545E-01_dp, &
     4.7272727E-01_dp,   4.9090909E-01_dp,   5.0909091E-01_dp,   5.2727273E-01_dp,   5.4545455E-01_dp, &
     5.6363636E-01_dp,   5.8181818E-01_dp,   6.0000000E-01_dp,   6.1818182E-01_dp,   6.3636364E-01_dp, &
     6.5454545E-01_dp,   6.7272727E-01_dp,   6.9090909E-01_dp,   7.0909091E-01_dp,   7.2727273E-01_dp, &
     7.4545455E-01_dp,   7.6363636E-01_dp,   7.8181818E-01_dp,   8.0000000E-01_dp,   8.0444444E-01_dp, &
     8.0888889E-01_dp,   8.1333333E-01_dp,   8.1777778E-01_dp,   8.2222222E-01_dp,   8.2666667E-01_dp, &
     8.3111111E-01_dp,   8.3555556E-01_dp,   8.4000000E-01_dp,   8.4444444E-01_dp,   8.4888889E-01_dp, &
     8.5333333E-01_dp,   8.5777778E-01_dp,   8.6222222E-01_dp,   8.6666667E-01_dp,   8.7111111E-01_dp, &
     8.7555556E-01_dp,   8.8000000E-01_dp,   8.8444444E-01_dp,   8.8888889E-01_dp,   8.9333333E-01_dp, &
     8.9777778E-01_dp,   9.0222222E-01_dp,   9.0666667E-01_dp,   9.1111111E-01_dp,   9.1555556E-01_dp, &
     9.2000000E-01_dp,   9.2444444E-01_dp,   9.2888889E-01_dp,   9.3333333E-01_dp,   9.3777778E-01_dp, &
     9.4222222E-01_dp,   9.4666667E-01_dp,   9.5111111E-01_dp,   9.5555556E-01_dp,   9.6000000E-01_dp, &
     9.6444444E-01_dp,   9.6888889E-01_dp,   9.7333333E-01_dp,   9.7777778E-01_dp,   9.8222222E-01_dp, &
     9.8666667E-01_dp,   9.9111111E-01_dp,   9.9555556E-01_dp,   1.0000000E+00_dp, 0.0_dp /)

  real(rprec), dimension(1:l) :: am_aux_f = (/ 1.1343073E+05_dp, &
     1.0511264E+05_dp,   9.7719901E+04_dp,   9.1142143E+04_dp,   8.5238372E+04_dp,   7.9893902E+04_dp, &
     7.5021684E+04_dp,   7.0555497E+04_dp,   6.6444168E+04_dp,   6.2647239E+04_dp,   5.9132114E+04_dp, &
     5.5871944E+04_dp,   5.2844250E+04_dp,   5.0029862E+04_dp,   4.7412214E+04_dp,   4.4976739E+04_dp, &
     4.2710495E+04_dp,   4.0601792E+04_dp,   3.8639948E+04_dp,   3.6815020E+04_dp,   3.5117615E+04_dp, &
     3.3538679E+04_dp,   3.2069304E+04_dp,   3.0700518E+04_dp,   2.9423071E+04_dp,   2.8227198E+04_dp, &
     2.7102369E+04_dp,   2.6037002E+04_dp,   2.5018180E+04_dp,   2.4031951E+04_dp,   2.3067082E+04_dp, &
     2.2117542E+04_dp,   2.1181454E+04_dp,   2.0260068E+04_dp,   1.9356991E+04_dp,   1.8477613E+04_dp, &
     1.7628651E+04_dp,   1.6817841E+04_dp,   1.6053656E+04_dp,   1.5345090E+04_dp,   1.4699130E+04_dp, &
     1.4113009E+04_dp,   1.3572603E+04_dp,   1.3055378E+04_dp,   1.2532538E+04_dp,   1.2400160E+04_dp, &
     1.2264938E+04_dp,   1.2126344E+04_dp,   1.1983859E+04_dp,   1.1836908E+04_dp,   1.1684914E+04_dp, &
     1.1527343E+04_dp,   1.1363692E+04_dp,   1.1193491E+04_dp,   1.1016295E+04_dp,   1.0831470E+04_dp, &
     1.0638256E+04_dp,   1.0435957E+04_dp,   1.0224325E+04_dp,   1.0003778E+04_dp,   9.7747634E+03_dp, &
     9.5363620E+03_dp,   9.2816431E+03_dp,   9.0022662E+03_dp,   8.6903284E+03_dp,   8.3409886E+03_dp, &
     7.9527438E+03_dp,   7.5237757E+03_dp,   7.0538828E+03_dp,   6.5554498E+03_dp,   6.0461036E+03_dp, &
     5.5435531E+03_dp,   5.0630772E+03_dp,   4.6140337E+03_dp,   4.2046228E+03_dp,   3.8422577E+03_dp, &
     3.5274335E+03_dp,   3.2543658E+03_dp,   3.0169701E+03_dp,   2.8089885E+03_dp,   2.6245185E+03_dp, &
     2.4577789E+03_dp,   2.3032608E+03_dp,   2.1572214E+03_dp,   2.0204623E+03_dp,   1.8947411E+03_dp, &
     1.7810514E+03_dp,   1.6799853E+03_dp,   1.5917883E+03_dp,   1.5163600E+03_dp, 0.0_dp /)

  integer :: i, j, iflag
  integer, parameter :: ns = 25
  real(rprec) :: hs, si
  real(rprec) :: y

  j = minloc(am_aux_s(2:), dim=1)
  print *, "length : ", j

  hs = 1.0_rprec/(ns-1)

  DO i = 2, ns
    si = hs * (i-1.5_rprec)

    call spline_akima(si, y, am_aux_s, am_aux_f, j, iflag)

    if (iflag .eq. 0) then
      print *, "si = ", si, " => y = ", y
    else
      print *, "error: iflag = ", iflag
    end if

  end do

end program Hello

SUBROUTINE spline_akima(x, y, xx, yy, npts, iflag)
      implicit none

      INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND(12,100)
      INTEGER, PARAMETER :: dp = rprec

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,     intent(in)                  :: npts  !< number of active points
      real(rprec), intent(in)                  :: x     !< location at which to evaluate
      real(rprec), intent(out)                 :: y     !< evaluated function value
      real(rprec), dimension(npts), intent(in) :: xx    !< knots  of data to interpolate
      real(rprec), dimension(npts), intent(in) :: yy    !< values of data to interpolate
      integer,     intent(inout)               :: iflag !< error flag; 0 means all good
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(rprec), dimension(-1:size(xx)+2) :: a !<  constant coefficients
      real(rprec), dimension(-1:size(xx)+2) :: b !< 1st-order coefficients
      real(rprec), dimension(-1:size(xx)+2) :: c !< 2nd-order coefficients
      real(rprec), dimension(-1:size(xx)+2) :: d !< 3rd-order coefficients

      real(rprec), dimension(-1:size(xx)+2) :: xloc !< copy of xx, extended by 2 points at start and end
      real(rprec), dimension(-1:size(xx)+2) :: yloc !< copy of yy, quadratically-extrapolated beyond edge

      real(rprec), dimension(-1:size(xx)+2) :: m  !< linear derivatives
      real(rprec), dimension(-1:size(xx)+2) :: t  !< intermediate quantity
      real(rprec), dimension(-1:size(xx)+2) :: dm !< intermediate quantity
      real(rprec), dimension(-1:size(xx)+2) :: p  !< intermediate quantity
      real(rprec), dimension(-1:size(xx)+2) :: q  !< intermediate quantity

      integer :: i    !< iteration variable
      integer :: ix   !< size of xx, yy
      integer :: iv   !< number of points to actually use: 1:iv
      integer :: ivm1 !< iv-1

      real(rprec) :: cl !< intermediate quantity
      real(rprec) :: bl !< intermediate quantity
      real(rprec) :: cr !< intermediate quantity
      real(rprec) :: br !< intermediate quantity

      real(rprec) :: dx !< temporary variable

!-----------------------------------------------

      !WRITE(6,*)  x, y
      !WRITE(6,*)npts
      !WRITE(6,*) xx
      !WRITE(6,*) yy
      !stop 'Test!'

      iflag = 0  !initialization

      ix = size(xx)
!     print *, size(xx)
!     print *, size(yy)
      if(npts > ix) &
         stop'spline_akima: more active points requested than available'
      if(ix /= size(yy)) stop 'size mismatch of xx and yy!'

      iv=npts

      ivm1  = iv-1
! initialize local variables
      a = 0._DP ; b = 0._DP ; c = 0._DP ; d = 0._DP
      xloc = 0._DP ; yloc = 0._DP
      m = 0._DP ; t = 0._DP ; dm = 0._DP
      p = 0._DP ; q = 0._DP

      xloc(1:iv)=xx
      xloc(-1)= 2*xloc(1)-xloc(3)
      xloc( 0)= xloc(1)+xloc(2)-xloc(3)
      xloc(iv+2)= 2*xloc(iv)-xloc(iv-2)
      xloc(iv+1)= xloc(iv)+xloc(iv-1)-xloc(iv-2)

      !print *, xloc(-1), xloc(0), xloc(iv+1), xloc(iv+2)

      yloc(1:iv)=yy

! calculate linear derivatives as far as existent
      m(-1:iv+1) = (yloc(0:iv+2)-yloc(-1:iv+1))/ &
                   (xloc(0:iv+2)-xloc(-1:iv+1))
! values for i=0, -1 and iv, iv+1 by quadratic extrapolation:
      cl = (m(2)-m(1))/(xloc(3)-xloc(1))
      bl = m(1) - cl*(xloc(2)-xloc(1))
      cr = (m(iv-2)-m(iv-1))/(xloc(iv)-xloc(iv-2))
      br = m(iv-2) - cr*(xloc(iv-1)-xloc(iv-2))

      !print *, cl, bl, cr, br

      yloc( 0)=yloc(1)+bl*(xloc( 0)-xloc(1))+ &
                       cl*(xloc( 0)-xloc(1))**2
      yloc(-1)=yloc(1)+bl*(xloc(-1)-xloc(1))+ &
                       cl*(xloc(-1)-xloc(1))**2
      yloc(iv+1)=yloc(iv)+br*(xloc(iv+1)-xloc(iv))+ &
                          cl*(xloc(iv+1)-xloc(iv))**2
      yloc(iv+2)=yloc(iv)+br*(xloc(iv+2)-xloc(iv))+ &
                          cl*(xloc(iv+2)-xloc(iv))**2

      !print *, yloc(-1), yloc(0), yloc(iv+1), yloc(iv+2)

! rest of linear derivatives
      m(-1) = (yloc(0)-yloc(-1))/(xloc(0)-xloc(-1))
      m( 0) = (yloc(1)-yloc( 0))/(xloc(1)-xloc( 0))
      m(iv  ) = (yloc(iv+1)-yloc(iv  ))/(xloc(iv+1)-xloc(iv  ))
      m(iv+1) = (yloc(iv+2)-yloc(iv+1))/(xloc(iv+2)-xloc(iv+1))

      !print *, m(-1), m(0), m(iv), m(iv+1)

! calculate weights for derivatives
      dm(-1:iv)= abs(m(0:iv+1)-m(-1:iv))
      where (dm /= 0._DP) !exclude division by zero
        p(1:iv) = dm(1:iv)/(dm(1:iv)+dm(-1:iv-2))
      end where
      where (dm /= 0._DP) !exclude division by zero
        q(1:iv) = dm(-1:iv-2)/(dm(1:iv)+dm(-1:iv-2))
      end where

      print *, "    i           p(i)          q(i)"
      do i=1, iv
        print *, i, p(i), q(i)
      end do

      stop

      t(1:iv) = p(1:iv)*m(0:iv-1)+q(1:iv)*m(1:iv)
      where ( p(1:iv)+q(1:iv) < TINY(1._DP)) ! in case of two zeros give equal weight
        t(1:iv) = 0.5_DP*m(0:iv-1)+0.5_DP*m(1:iv)
      end where
! fix coefficients
      a = yloc
      b = t
      c(1:iv-1) = (3*m(1:iv-1)-t(2:iv)-2*t(1:iv-1))/ &
                     (xloc(2:iv)-xloc(1:iv-1))
      d(1:iv-1) = (t(2:iv)+t(1:iv-1)-2*m(1:iv-1))/ &
                     (xloc(2:iv)-xloc(1:iv-1))**2

! calculation
      if(x<xloc(1) .or. x>xloc(iv)) then
        y=0.0
        iflag=-1
        return
      endif

      if(x == xloc(iv)) then
        y = yy(iv)
        iflag = 0
        return
      endif

      do i=1,iv-1
        if(x >= xloc(i) .and. x < xloc(i+1))then
          dx=x-xloc(i)
          y=a(i)+dx*(b(i)+dx*(c(i)+d(i)*dx))
          iflag = 0
          return
        endif
      enddo

      END SUBROUTINE spline_akima
