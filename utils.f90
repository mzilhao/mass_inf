module utils
  implicit none
  private
  public :: relative_difference, print_status

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


!> Report progress to stdout: iteration, v, elapsed time, rate, min/max(r)
subroutine print_status(iter, v_cur, v_initial, v_final, start_cpu, h_v0, next_idx, progress_stride, header_stride)
  integer, intent(in) :: iter, next_idx
  double precision, intent(in) :: v_cur, v_initial, v_final, start_cpu
  double precision, dimension(:,:), intent(in) :: h_v0
  integer, intent(in) :: progress_stride, header_stride

  double precision :: elapsed, rate, min_r, max_r
  logical :: should_print, should_header

  if (progress_stride <= 0) return

  should_print = (iter == 1) .or. (mod(iter, progress_stride) == 0)
  if (.not. should_print) return

  should_header = (iter == 1)
  if (header_stride > 0) then
    should_header = should_header .or. (mod(iter, header_stride) == 0)
  end if

  if (should_header) then
    write(*,'(a)') '-----------------------------------------------------------------------'
    write(*,'(a)') 'Iteration       v / vf     |         Δv/min |               r          '
    write(*,'(a)') '                           |                |     minimum       maximum'
    write(*,'(a)') '-----------------------------------------------------------------------'
  end if

  call cpu_time(elapsed)
  elapsed = (elapsed - start_cpu)/60.0d0  ! convert to minutes
  if (elapsed > 1.0d-9) then
    rate = (v_cur - v_initial) / elapsed
  else
    rate = 0.0d0
  end if

  ! this assumes that the r variable is stored in h_v0(:,1).
  ! needs to be adapted adapt if that ever changes.
  min_r = minval(h_v0(1:next_idx-1, 1))
  max_r = maxval(h_v0(1:next_idx-1, 1))

  write(*,'(i8,2x,f7.3,a,f6.2,a,f15.3,a,f12.4,2x,f12.4)') &
    iter, v_cur, ' /', v_final, '  |', rate, ' |', min_r, max_r
end subroutine print_status

end module utils
