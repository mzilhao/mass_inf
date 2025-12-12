program mass_inflation
  use physics_config_mod
  use simulation_config_mod
  use functions
  use evolve_wrapper, only: step
  use polint_mod
  use utils
  use progress_utils, only: report_progress, print_banner
  implicit none

  ! Physics and simulation configuration
  type(physics_config)    :: cfg
  type(simulation_config) :: sim_cfg
  integer                 :: neq, output_unit
  double precision        :: upos, v, reldiff_r = 0.0d0

  double precision, allocatable, dimension(:,:) :: h_u0, h_v0
  integer,          allocatable, dimension(:)   :: plus, minus
  double precision, allocatable                 :: u(:)
  double precision, allocatable, dimension(:,:) :: h_v1
  double precision, allocatable, dimension(:)   :: h_S, h_E, h_W, h_N

  double precision :: mass, drdv, ricci
  double precision :: du, dv, u0, v0, uf, vf
  double precision :: start_time_cpu
  integer          :: Nv, Nu, Nu_max
  double precision :: tempu, tempv

  integer :: i, j, k, jm1, jm2, jm3, jp1, next_idx
  character(len=20), parameter :: filename = 'data.dat'

  ! Named constants for clarity
  integer, parameter :: N_PICARD_ITERATIONS = 4
  integer, parameter :: N_INTERP_POINTS = 5
  double precision, parameter :: OUTPUT_TOLERANCE_V = 1.0d-6
  double precision, parameter :: OUTPUT_TOLERANCE_U = 1.0d-7
  double precision, dimension(N_INTERP_POINTS) :: interp_x, interp_y

  call init_physics_config(cfg)
  call init_simulation_config(sim_cfg)
  neq = cfg%neq

  call print_banner()

  ! Create local aliases for readability
  du = sim_cfg%du
  dv = sim_cfg%dv
  u0 = sim_cfg%u0
  v0 = sim_cfg%v0
  uf = sim_cfg%uf
  vf = sim_cfg%vf
  Nu = sim_cfg%Nu
  Nv = sim_cfg%Nv
  Nu_max = sim_cfg%Nu_max

  call cpu_time(start_time_cpu)

  open(newunit=output_unit, file=filename, status='replace')
  call write_output_header(output_unit, cfg)

  ! Initialize boundary conditions (returns allocated h_u0, h_v0)
  call init_cond(h_u0, h_v0, sim_cfg, cfg)

  allocate(plus(Nu_max), minus(Nu_max))
  allocate(u(Nu_max))
  allocate(h_v1(Nu_max, neq))
  allocate(h_S(neq), h_E(neq), h_W(neq), h_N(neq))

  ! Array 'u' will hold all values of u on v slices. We start with
  ! a uniform grid, but this may change with AMR. Linked lists 'plus' and 'minus'
  ! help us navigate the non-uniform grid.
  u = 0.0d0
  plus = 0
  minus = 0

  upos = u0
  do j = 1, Nu
    u(j) = u0 + (j - 1) * du
  end do

  do j = 1, Nu_max
    minus(j) = j - 1
    plus(j)  = j + 1
  end do

  next_idx = Nu + 1
  v = sim_cfg%v0


  ! Start the main integration loop. i is the step in 'v'; j the step in 'u'.
  ! At each step we assume we are at the point (u,v).
  do i = 1, Nv - 1

    ! FIXME
    ! Output separator line, for easier parsing of data files
    tempv = abs(v - v0) * sim_cfg%resv
    if (tempv - int(tempv + dv*0.1d0) < OUTPUT_TOLERANCE_V) then
      call write_output_separator(output_unit)
    end if

    ! Progress output to stdout (cadence controlled by simulation config)
    call report_progress(i, v, vf, start_time_cpu, h_v0, next_idx, sim_cfg%progress_stride, sim_cfg%progress_header_stride)

    ! Reset u position each time we advance in v
    upos = u0

    v = v + dv

    ! h_N <- h(u, v + dv)
    h_N(:) = h_u0(i + 1, :)

    ! h_v1 will store all values at points (u, v + dv)
    h_v1(1, :) = h_N(:)


    ! Advance from u to u + du at v + dv
    j = plus(1)
    do while (upos < uf)
      jm1 = minus(j)

      ! h_E holds the values that h_N had in the previous step: h(u, v + dv)
      h_E(:) = h_N(:)        ! h(u, v + dv)
      h_S(:) = h_v0(jm1, :)  ! h(u, v)
      h_W(:) = h_v0(j, :)    ! h(u + du, v)

      ! Adaptive Mesh Refinement (AMR) in the 'u' direction. If the relative
      ! variation in 'r' between points (u, v) and (u + du, v) exceeds the threshold,
      ! we add a new point halfway between them by interpolating all
      ! field values using polynomial interpolation.
      if (sim_cfg%AMR) then
        reldiff_r = relative_difference(h_W(1), h_S(1))

        ! We keep adding points in u until the gradient in r is small enough
        do while (reldiff_r > sim_cfg%gradmax .and. j >= 4)
          jm1 = minus(j)
          jm2 = minus(jm1)
          jm3 = minus(jm2)
          jp1 = plus(j)

          do k = 1, neq
            if (u(jp1) < uf) then
              interp_x = (/ u(1), u(jm2), u(jm1), u(j), u(jp1) /)
              interp_y = (/ h_v0(1,k), h_v0(jm2,k), h_v0(jm1,k), h_v0(j,k), &
                            h_v0(jp1,k) /)
            else
              interp_x = (/ u(1), u(jm3), u(jm2), u(jm1), u(j) /)
              interp_y = (/ h_v0(1,k), h_v0(jm3,k), h_v0(jm2,k), h_v0(jm1,k), &
                            h_v0(j,k) /)
            end if
            h_W(k) = polint((u(j) + u(jm1)) * 0.5d0, interp_x, interp_y)
          end do

          u(next_idx)       = (u(j) + u(jm1)) * 0.5d0
          h_v0(next_idx, :) = h_W(:)

          minus(next_idx) = jm1
          plus(next_idx)  = j
          plus(jm1)       = next_idx
          minus(j)        = next_idx
          j               = next_idx
          next_idx        = next_idx + 1

          reldiff_r = relative_difference(h_W(1), h_S(1))
        end do
      end if ! end AMR

      du = u(j) - u(minus(j))  ! local du may change during AMR

      ! step returns h_N = h(u + du, v + dv)
      call step(h_N, h_S, h_E, h_W, du, dv, cfg, N_PICARD_ITERATIONS)

      h_v1(j, :) = h_N(:)
      upos = upos + du

      call write_output_if_needed(output_unit, u(j), v, h_N, h_S, h_E, h_W, du, dv, &
                                   u0, v0, sim_cfg, cfg, OUTPUT_TOLERANCE_U, OUTPUT_TOLERANCE_V)

      j = plus(j)

    end do ! do while in the u direction

    ! Copy active grid points from h_v1 back to h_v0 for the next v-step.
    ! Only copy rows 1 to next_idx (active points after AMR refinement) to avoid
    ! unnecessary copying of unused array memory.
    h_v0(1:next_idx, :) = h_v1(1:next_idx, :)

  end do

  close(output_unit)

  deallocate(h_u0, h_v0, u, minus, plus, h_v1)

end program mass_inflation
