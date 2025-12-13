!> Simulation Configuration Module
!! 
!! Encapsulates all numerical/integration parameters independent of physics.
!! This allows swapping physics implementations (functions.f90) while keeping
!! the numerical solver infrastructure constant.

module simulation_config_mod
  implicit none

  !> Numerical simulation configuration type
  !! Stores all grid, integration, and AMR parameters
  type :: simulation_config
    ! Integration domain bounds
    double precision :: u_min = 0.0d0, v_min = 5.0d0   ! Domain start
    double precision :: u_max = 30.0d0, v_max = 10.0d0 ! Domain end

    ! Integration step sizes
    double precision :: du = 0.01d0, dv = 0.0005d0     ! Step sizes in u, v directions

    ! AMR (Adaptive Mesh Refinement) parameters
    logical :: AMR = .true.                            ! Enable AMR
    double precision :: reldiff_max = 0.0001d0         ! Thresholds for refinement

    ! Output sampling spacing (absolute values in domain units)
    double precision :: output_du = 0.05d0             ! Write every this Δu
    double precision :: output_dv = 0.05d0             ! Write every this Δv

    ! Progress reporting cadence (stdout)
    integer :: progress_stride = 100
    integer :: progress_header_stride = 5000

    ! FIXME: move to physics config?
    ! Initial conditions
    double precision :: m0 = 1.0d0                     ! Initial mass parameter

    ! Derived values (dummy initializations, to be computed later)
    integer :: Nu = 0, Nv = 0                          ! Number of grid points
    integer :: Nu_max = 0                              ! Maximum u-arrays allocation size
  end type simulation_config

contains

!> Read simulation configuration from a namelist file
!! Reads &simulation namelist and computes derived grid dimensions.
!! Falls back to defaults if file doesn't exist or read fails.
subroutine read_simulation_config_from_file(sim_cfg, filename)
  type(simulation_config), intent(out) :: sim_cfg
  character(len=*), intent(in)         :: filename

  ! Local variables for namelist reading
  double precision :: u_min, v_min, u_max, v_max, du, dv, reldiff_max
  double precision :: output_du, output_dv, m0
  logical :: AMR
  integer :: progress_stride, progress_header_stride
  namelist /simulation/ u_min, v_min, u_max, v_max, du, dv, AMR, reldiff_max, &
                        output_du, output_dv, progress_stride, progress_header_stride, m0

  integer :: unit, ierr

  ! Initialize with type defaults
  sim_cfg = simulation_config()
  u_min = sim_cfg%u_min
  v_min = sim_cfg%v_min
  u_max = sim_cfg%u_max
  v_max = sim_cfg%v_max
  du    = sim_cfg%du
  dv    = sim_cfg%dv
  AMR   = sim_cfg%AMR
  reldiff_max = sim_cfg%reldiff_max
  output_du = sim_cfg%output_du
  output_dv = sim_cfg%output_dv
  progress_stride        = sim_cfg%progress_stride
  progress_header_stride = sim_cfg%progress_header_stride

  m0 = sim_cfg%m0

  ! Read namelist
  open(newunit=unit, file=filename, status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
    write(*, '(a,a,a)') 'Error: could not open parameter file "', trim(filename), '"'
  else
    read(unit, nml=simulation, iostat=ierr)
    if (ierr /= 0) then
      write(*, '(a,a,a)') 'Error: could not read &simulation namelist from "', trim(filename), '"'
    end if
    close(unit)
  end if

  ! Update sim_cfg with (possibly modified) namelist values
  sim_cfg%u_min = u_min
  sim_cfg%v_min = v_min
  sim_cfg%u_max = u_max
  sim_cfg%v_max = v_max
  sim_cfg%du = du
  sim_cfg%dv = dv
  sim_cfg%AMR = AMR
  sim_cfg%reldiff_max = reldiff_max
  sim_cfg%output_du = output_du
  sim_cfg%output_dv = output_dv
  sim_cfg%progress_stride = progress_stride
  sim_cfg%progress_header_stride = progress_header_stride

  sim_cfg%m0 = m0

  ! Compute derived grid dimensions
  call compute_grid_dimensions(sim_cfg)

end subroutine read_simulation_config_from_file

!> Compute derived grid dimensions from domain and step sizes
subroutine compute_grid_dimensions(sim_cfg)
  type(simulation_config), intent(inout) :: sim_cfg
  sim_cfg%Nu = int((sim_cfg%u_max - sim_cfg%u_min) / sim_cfg%du + 1.001d0)
  sim_cfg%Nv = int((sim_cfg%v_max - sim_cfg%v_min) / sim_cfg%dv + 1.001d0)
  sim_cfg%Nu_max = int(2.0d0 * (sim_cfg%u_max - sim_cfg%u_min) / sim_cfg%reldiff_max)
end subroutine compute_grid_dimensions

end module simulation_config_mod
