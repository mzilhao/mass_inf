!> Grid Configuration Module
!! 
!! Encapsulates all numerical/integration parameters independent of the model.
!! This allows swapping model implementations (models/) while keeping
!! the numerical solver infrastructure.

module grid_config_mod
  use precision
  implicit none
  private
  public :: grid_config, load

  !> Numerical grid configuration type
  !! Stores all grid, integration, and AMR parameters
  type :: grid_config
    ! Integration domain bounds
    real(dp) :: u_min = 0.0_dp, v_min = 5.0_dp   ! Domain start
    real(dp) :: u_max = 30.0_dp, v_max = 10.0_dp ! Domain end

    ! Integration step sizes
    real(dp) :: du = 0.01_dp, dv = 0.0005_dp     ! Step sizes in u, v directions

    ! AMR (Adaptive Mesh Refinement) parameters
    logical :: AMR = .true.                      ! Enable AMR
    real(dp) :: reldiff_max = 0.0001_dp          ! Thresholds for refinement

    ! Output sampling spacing (absolute values in domain units)
    real(dp) :: output_du = 0.05_dp              ! Write every this Δu
    real(dp) :: output_dv = 0.05_dp              ! Write every this Δv
    character(len=256) :: output_base_dir = ''   ! Base dir under which run folder is created ('' = CWD)

    ! Progress reporting cadence (stdout)
    integer :: progress_stride = 100
    integer :: progress_header_stride = 5000

    ! Derived values (dummy initializations, to be computed later)
    integer :: Nu = 0, Nv = 0                    ! Number of grid points
    integer :: Nu_max = 0                        ! Maximum u-arrays allocation size
  end type grid_config

contains

!> Read grid configuration from a namelist file
!! Reads &grid namelist and computes derived grid dimensions.
!! Falls back to defaults for parameters not specified in the file.
subroutine load(grid_cfg, filename)
  type(grid_config), intent(out) :: grid_cfg
  character(len=*), intent(in)   :: filename

  ! Local variables for namelist reading
  real(dp) :: u_min, v_min, u_max, v_max, du, dv, reldiff_max
  real(dp) :: output_du, output_dv
  character(len=256) :: output_base_dir
  logical :: AMR
  integer :: progress_stride, progress_header_stride
  namelist /grid/ u_min, v_min, u_max, v_max, du, dv, AMR, reldiff_max, &
                        output_du, output_dv, output_base_dir, progress_stride, progress_header_stride

  integer :: unit, ierr

  ! Initialize with type defaults
  grid_cfg = grid_config()

  u_min = grid_cfg%u_min
  v_min = grid_cfg%v_min
  u_max = grid_cfg%u_max
  v_max = grid_cfg%v_max
  du    = grid_cfg%du
  dv    = grid_cfg%dv
  AMR   = grid_cfg%AMR
  reldiff_max = grid_cfg%reldiff_max
  output_du = grid_cfg%output_du
  output_dv = grid_cfg%output_dv
  output_base_dir = grid_cfg%output_base_dir
  progress_stride        = grid_cfg%progress_stride
  progress_header_stride = grid_cfg%progress_header_stride

  ! Read namelist
  open(newunit=unit, file=filename, status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
    write(*, '(a,a,a)') 'Error: could not open parameter file "', trim(filename), '"'
    call exit(1)
  end if

  read(unit, nml=grid, iostat=ierr)
  if (ierr /= 0) then
    write(*, '(a,a,a)') 'Error: could not read &grid namelist from "', trim(filename), '"'
    close(unit)
    call exit(1)
  end if
  close(unit)

  ! Update grid_cfg with (possibly modified) namelist values
  grid_cfg%u_min = u_min
  grid_cfg%v_min = v_min
  grid_cfg%u_max = u_max
  grid_cfg%v_max = v_max
  grid_cfg%du = du
  grid_cfg%dv = dv
  grid_cfg%AMR = AMR
  grid_cfg%reldiff_max = reldiff_max
  grid_cfg%output_du = output_du
  grid_cfg%output_dv = output_dv
  grid_cfg%output_base_dir = output_base_dir
  grid_cfg%progress_stride = progress_stride
  grid_cfg%progress_header_stride = progress_header_stride

  ! Compute derived grid dimensions
  call compute_grid_dimensions(grid_cfg)

end subroutine load

!> Compute derived grid dimensions from domain and step sizes
subroutine compute_grid_dimensions(grid_cfg)
  type(grid_config), intent(inout) :: grid_cfg

  grid_cfg%Nu = int((grid_cfg%u_max - grid_cfg%u_min) / grid_cfg%du + 1.001_dp)
  grid_cfg%Nv = int((grid_cfg%v_max - grid_cfg%v_min) / grid_cfg%dv + 1.001_dp)
  grid_cfg%Nu_max = int(2.0_dp * (grid_cfg%u_max - grid_cfg%u_min) / grid_cfg%reldiff_max)

end subroutine compute_grid_dimensions

end module grid_config_mod
