module utils
  implicit none
  private
  public :: relative_difference

contains

!> Compute relative difference between two values
!! Used for AMR criterion based on the variation of the r variable.
!!
!! @param[in] val1
!! @param[in] val2
!! @return 2 * |val1 - val2| / (|val1| + |val2|)
pure function relative_difference(val1, val2) result(grad)
  double precision, intent(in) :: val1, val2
  double precision :: grad
  double precision, parameter :: EPSILON_GUARD = 1.0d-16

  ! The factor of 2 normalizes the gradient relative to the average value.
  grad = 2.0d0 * abs(val1 - val2) / (abs(val1) + abs(val2) + EPSILON_GUARD)
end function relative_difference

end module utils
