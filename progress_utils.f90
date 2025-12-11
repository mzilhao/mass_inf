module progress_utils
  implicit none
  private
  public :: report_progress

contains

  !> Report progress to stdout: iteration, v, elapsed time, rate, min/max(r)
  subroutine report_progress(iter, v_cur, v_final, start_cpu, h_v0, next_idx, progress_stride, header_stride)
    integer, intent(in) :: iter, next_idx
    double precision, intent(in) :: v_cur, v_final, start_cpu
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
      write(*,'(a)') 'Iteration       v / vf     | rate (iters/s) |               r          '
      write(*,'(a)') '                           |                |     minimum       maximum'
      write(*,'(a)') '-----------------------------------------------------------------------'
    end if

    call cpu_time(elapsed)
    elapsed = max(elapsed - start_cpu, 0.0d0)
    if (elapsed > 1.0d-9) then
      rate = dble(iter) / elapsed
    else
      rate = 0.0d0
    end if

    min_r = minval(h_v0(1:next_idx-1, 1))
    max_r = maxval(h_v0(1:next_idx-1, 1))

    write(*,'(i8,2x,f7.3,a,f6.2,a,f15.3,a,f12.4,2x,f12.4)') &
      iter, v_cur, ' /', v_final, '  |', rate, ' |', min_r, max_r
  end subroutine report_progress

end module progress_utils
