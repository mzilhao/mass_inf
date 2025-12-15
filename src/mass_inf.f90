program mass_inflation
  use precision
  use model_config_mod
  use grid_config_mod
  use model_mod
  use evolve_wrapper, only: step
  use polint_mod
  use utils
  use amr_mod
  implicit none

  ! Physics and simulation configuration
  type(model_config)    :: model_cfg
  type(grid_config)     :: grid_cfg

  real(dp), allocatable, dimension(:,:) :: h_u0, h_v0
  integer,  allocatable, dimension(:)   :: plus, minus
  real(dp), allocatable                 :: u(:)
  real(dp), allocatable, dimension(:,:) :: h_v1
  real(dp), dimension(NEQ)              :: h_S, h_E, h_W, h_N

  real(dp) :: u_cur, v_cur, reldiff_r = 0.0_dp
  real(dp) :: du, dv, u_min, v_min, u_max, v_max
  real(dp) :: start_time_cpu
  integer  :: Nv, Nu, Nu_max
  integer  :: i, j, k, jm1, next_idx

  integer, parameter :: N_PICARD_ITERATIONS = 4
  integer, parameter :: N_INTERP_POINTS = 5

  ! IO and file management
  character(len=256) :: out_dir, param_file, arg, prog_name
  logical :: param_file_exists, force_overwrite
  integer :: num_args, arg_idx


  ! Parse command line arguments
  call get_command_argument(0, prog_name)
  num_args = command_argument_count()
  force_overwrite = .false.

  if (num_args < 1 .or. num_args > 2) then
    write(*, '(a,a,a)') 'Usage: ', trim(prog_name), ' [-f] <parameter_file.nml>'
    write(*, '(a)') ''
    write(*, '(a)') 'Options:'
    write(*, '(a)') '  -f                 force overwrite of existing output directory'
    call exit(1)
  end if

  ! Check for -f flag
  arg_idx = 1
  if (num_args == 2) then
    call get_command_argument(1, arg)
    if (trim(arg) == '-f') then
      force_overwrite = .true.
      arg_idx = 2
    else
      write(*, '(a,a,a)') 'Error: unknown option "', trim(arg), '"'
      write(*, '(a,a,a)') 'Usage: ', trim(prog_name), ' [-f] <parameter_file.nml>'
      call exit(1)
    end if
  end if

  ! Get parameter file path
  call get_command_argument(arg_idx, param_file)
  inquire(file=param_file, exist=param_file_exists)
  if (.not. param_file_exists) then
    write(*, '(a,a,a)') 'Error: file "', trim(param_file), '" not found.'
    call exit(1)
  end if

  ! Output directory name (parameter file basename without extension)
  out_dir = trim(param_file)
  ! Strip any leading path
  i = index(out_dir, '/', back=.true.)
  if (i > 0) out_dir = out_dir(i+1:)
  ! Strip extension if present
  j = index(out_dir, '.', back=.true.)
  if (j > 0) out_dir = out_dir(:j-1)

  ! Read simulation and physics configurations
  call read_grid_config_from_file(grid_cfg, param_file)
  call read_model_config_from_file(model_cfg, param_file)

  ! Setup output directory and print startup banner
  call startup(param_file, out_dir, force_overwrite, grid_cfg%output_base_dir)

  ! Create local aliases for readability
  du     = grid_cfg%du
  dv     = grid_cfg%dv
  u_min  = grid_cfg%u_min
  v_min  = grid_cfg%v_min
  u_max  = grid_cfg%u_max
  v_max  = grid_cfg%v_max
  Nu     = grid_cfg%Nu
  Nv     = grid_cfg%Nv
  Nu_max = grid_cfg%Nu_max

  ! Open output files
  call open_output_files(out_dir)

  ! Allocate boundary condition arrays
  allocate(h_u0(Nv, NEQ))
  allocate(h_v0(Nu_max, NEQ))
  allocate(h_v1(Nu_max, NEQ))
  h_u0 = 0.0d0
  h_v0 = 0.0d0
  h_v1 = 0.0d0

  ! Array 'u' will hold all values of u on v slices. We start with
  ! a uniform grid, but this may change with AMR. Linked lists 'plus' and 'minus'
  ! help us navigate the non-uniform grid.
  allocate(u(Nu_max))
  allocate(plus(Nu_max), minus(Nu_max))
  u = 0.0d0
  plus  = 0
  minus = 0

  u_cur = u_min
  do j = 1, Nu
    u(j) = u_min + (j - 1) * du
  end do

  do j = 1, Nu_max
    minus(j) = j - 1
    plus(j)  = j + 1
  end do
  next_idx = Nu + 1


  ! Initialize boundary conditions
  call init_cond(h_u0, h_v0, grid_cfg, model_cfg)
  call cpu_time(start_time_cpu)

  v_cur = v_min

  ! Start the main integration loop. i is the step in 'v'; j the step in 'u'.
  ! At each step we assume we are at the point (u,v).
  do i = 1, Nv - 1
    ! Print progress to stdout (cadence controlled by simulation config)
    call print_status(i, v_cur, v_min, v_max, start_time_cpu, h_v0, next_idx, &
                      grid_cfg%progress_stride, grid_cfg%progress_header_stride)

    ! Reset u position each time we advance in v
    u_cur = u_min

    ! h_N <- h(u, v + dv)
    h_N(:) = h_u0(i + 1, :)

    ! h_v1 will store all values at points (u, v + dv)
    h_v1(1, :) = h_N(:)

    ! Advance from u to u + du
    j = plus(1)
    do while (u_cur < u_max)
      jm1 = minus(j)

      ! h_E holds the values that h_N had in the previous step: h(u, v + dv)
      h_E(:) = h_N(:)        ! h(u, v + dv)
      h_S(:) = h_v0(jm1, :)  ! h(u, v)
      h_W(:) = h_v0(j, :)    ! h(u + du, v)

      ! Adaptive Mesh Refinement
      if (grid_cfg%AMR) call refine_u_grid(u, h_v0, h_W, h_S, j, next_idx, u_max, grid_cfg%reldiff_max, plus, minus)

      du = u(j) - u(minus(j))  ! local du may change during AMR

      ! step returns h_N = h(u + du, v + dv)
      call step(h_N, h_S, h_E, h_W, du, dv, model_cfg, N_PICARD_ITERATIONS)
      h_v1(j, :) = h_N(:)

      call write_output(u_cur, v_cur, h_N, h_S, h_E, h_W, du, dv, grid_cfg, model_cfg)

      j = plus(j)
      u_cur = u_cur + du
    end do ! do while in the u direction

    ! Copy active grid points from h_v1 back to h_v0 for the next v-step.
    ! Only copy rows 1 to next_idx (active points after AMR refinement) to avoid
    ! unnecessary copying of unused array memory.
    h_v0(1:next_idx, :) = h_v1(1:next_idx, :)

    ! Advance in v
    v_cur = v_cur + dv
  end do

  ! Close output files
  call close_output_files()

  deallocate(h_u0, h_v0, u, minus, plus, h_v1)

end program mass_inflation
