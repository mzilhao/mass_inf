!> Polynomial interpolation routine
!! Copied from Numerical Recipes in C++, 2nd edition.
!! Performs Lagrange polynomial interpolation/extrapolation.

module polint_mod
  use precision
  implicit none
  private
  public :: polint

contains

!> Lagrange polynomial interpolation/extrapolation
!! @param x      Point to interpolate at
!! @param xa     Known x data points (array of size N)
!! @param ya     Known y data points (array of size N)
!! @return       Interpolated value at x
function polint(x, xa, ya)
  implicit none
  real(dp), intent(in) :: x
  real(dp), dimension(:), intent(in) :: xa, ya
  real(dp) :: polint

  integer  :: i, m, ns, N
  real(dp) :: ho, hp, den, diff, difftemp, w
  real(dp), dimension(size(xa)) :: C, D

  ns = 1
  N = size(ya)

  ! Find the closest point to x
  diff = abs(x - xa(1))
  do i = 1, N
    difftemp = abs(x - xa(i))
    if (difftemp < diff) then
      ns = i
      diff = difftemp
    end if

    C(i) = ya(i)   ! Initialize C and D arrays
    D(i) = ya(i)
  end do

  polint = ya(ns)  ! Initial estimate
  ns = ns - 1

  ! Build the interpolation table
  do m = 1, N - 1
    do i = 1, N - m
      ho = xa(i) - x
      hp = xa(i + m) - x
      w = C(i + 1) - D(i)

      den = ho - hp
      ! Error if two xa values are identical (within numerical precision)
      if (abs(den) < epsilon(den)*max(abs(ho), abs(hp))) then
        write(*,*) "Error in polint: denominator too small (coincident xa values)"
        stop
      end if

      den = w/den
      D(i) = hp*den   ! Update C and D
      C(i) = ho*den
    end do

    if (2*ns < N - m) then
      polint = polint + C(ns + 1)
    else
      polint = polint + D(ns)
      ns = ns - 1
    end if
  end do

end function polint

end module polint_mod
