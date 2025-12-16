module utils_mod
  implicit none
  private
  public :: relative_difference, print_status, startup, trim_filename

contains

!> Trim directory path and file extension from a filename
!!
!! @param[in] filename Full filename with path and extension
!! @return Trimmed filename without path and extension
function trim_filename(filename) result(trimmed)
  character(len=*), intent(in) :: filename
  character(len=len(filename)) :: trimmed
  integer :: last_slash, last_dot

  ! Strip any leading path using back=.true. to find the last '/'
  last_slash = index(filename, '/', back=.true.)
  if (last_slash > 0) then
    trimmed = filename(last_slash + 1:)
  else
    trimmed = filename
  end if

  ! Strip extension if present
  last_dot = index(trimmed, '.', back=.true.)
  if (last_dot > 0) then
    trimmed = trimmed(:last_dot - 1)
  end if

end function trim_filename

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


!> Print ASCII banner and general info to stdout at startup
subroutine startup(param_file, out_dir, force_overwrite, output_base_dir)
  character(len=256), intent(in) :: param_file
  character(len=256), intent(inout) :: out_dir
  logical, intent(in) :: force_overwrite
  character(len=256), intent(in) :: output_base_dir

  character(len=8) :: date_str
  character(len=10) :: time_str
  character(len=256) :: username, hostname
  integer :: stat
  logical :: out_dir_exists

  ! Get runtime info
  call date_and_time(DATE=date_str, TIME=time_str)
  call get_environment_variable('USER', username, status=stat)
  if (stat /= 0) username = 'unknown'
  call get_environment_variable('HOSTNAME', hostname, status=stat)
  if (stat /= 0) hostname = 'unknown'

  ! Setup output directory
  if (len_trim(output_base_dir) > 0) then
    out_dir = trim(output_base_dir) // '/' // trim(out_dir)
  end if

  if (len_trim(out_dir) > 0) then
    ! Check if output directory already exists
    inquire(file=trim(out_dir)//'/.', exist=out_dir_exists)

    if (out_dir_exists) then
      if (force_overwrite) then
        write(*, '(a,a,a)') 'Warning: removing existing directory "', trim(out_dir), '"'
        call execute_command_line('rm -rf "' // trim(out_dir) // '"')
      else
        write(*, '(a,a,a)') 'Error: output directory "', trim(out_dir), '" already exists.'
        write(*, '(a)') 'Use -f flag to force overwrite, or remove the directory manually.'
        call exit(1)
      end if
    end if

    call execute_command_line('mkdir -p ' // trim(out_dir))
    ! Copy input parameter file to output directory for reproducibility
    call execute_command_line('cp "' // trim(param_file) // '" "' // trim(out_dir) // '/"')
  end if

  ! Print banner
  write(*,'(a)') '======================================================================='
  write(*,'(a)') ''
  write(*,'(a)') '                    x\                /x'
  write(*,'(a)') '                    x \              / x'
  write(*,'(a)') '                    x  \            /  x'
  write(*,'(a)') '                    x   \          /   x'
  write(*,'(a)') '                    x    \        /    x'
  write(*,'(a)') '                    x     \      /     x'
  write(*,'(a)') '                    x      \    /      x'
  write(*,'(a)') '                    x       \  /       x'
  write(*,'(a)') '                    x        \/        x'
  write(*,'(a)') '                    x        /\        x'
  write(*,'(a)') '                    x       /  \       x'
  write(*,'(a)') '                    x      /    \      x'
  write(*,'(a)') '                    x     /  /\_/\     x'
  write(*,'(a)') '                    x    /  ( O_O )    x'
  write(*,'(a)') '                    x   /    > ^ < \   x'
  write(*,'(a)') '                    x  /            \  x'
  write(*,'(a)') '                    x /              \ x'
  write(*,'(a)') '                    x/                \x'
  write(*,'(a)') '                     \                /\'
  write(*,'(a)') '                      \              /  \'
  write(*,'(a)') '                       \            /    \'
  write(*,'(a)') '                        \          /      \'
  write(*,'(a)') '                         \        /        \'
  write(*,'(a)') '                          \      /          \'
  write(*,'(a)') '                           \    /            \'
  write(*,'(a)') '                            \  /              \'
  write(*,'(a)') '                             \/                \'
  write(*,'(a)') '                              \                /'
  write(*,'(a)') '                               \              /'
  write(*,'(a)') '                                \            /'
  write(*,'(a)') '                                 \          /'
  write(*,'(a)') '                                  \        /'
  write(*,'(a)') '                                   \      /'
  write(*,'(a)') '                                    \    /'
  write(*,'(a)') '                                     \  /'
  write(*,'(a)') '                                      \/'
  write(*,'(a)') ''
  write(*,'(a)') '======================================================================='
  write(*,'(a)') ''

  write(*,'(a,a4,a,a2,a,a2)') 'Run date: ', date_str(1:4), '-', date_str(5:6), '-', date_str(7:8)
  write(*,'(a,a2,a,a2,a,a2)') 'Run time: ', time_str(1:2), ':', time_str(3:4), ':', time_str(5:6)
  write(*,'(a,a)') 'User:     ', trim(username)
  write(*,'(a,a)') 'Host:     ', trim(hostname)
  write(*,'(a)') ''

end subroutine startup

end module utils_mod
