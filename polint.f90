! Polynomial interpolation routine
! Copied from Numerical Recipes in C++, 2nd edition.
! Performs Lagrange polynomial interpolation/extrapolation.

module polint_mod
  implicit none
  private
  public :: polint

contains

!> Lagrange polynomial interpolation/extrapolation
!! Input:
!!   xa, ya : known data points (arrays of size N)
!!   x      : point to interpolate at
!! Output:
!!   polint : interpolated value at x
function polint(x, xa, ya)
  implicit none

  double precision :: polint
  double precision, dimension(:), intent(in) :: xa, ya
  double precision, intent(in) :: x

  integer :: i, m, ns, N
  double precision :: ho, hp, den, diff, difftemp, w
  double precision, allocatable, dimension(:) :: C, D

  ns = 1
  N = size(ya)

  allocate(C(N), D(N))

  ! Find the closest point to x
  diff = abs(x - xa(1))
  do i = 1, N
    difftemp = abs(x - xa(i))
    if (difftemp < diff) then
      ns = i
      diff = difftemp
    end if

    C(i) = ya(i)  ! Initialize C and D arrays
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

  deallocate(C, D)

end function polint

end module polint_mod
