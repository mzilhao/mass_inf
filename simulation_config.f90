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
    double precision :: u_min, v_min            ! Domain start
    double precision :: u_max, v_max            ! Domain end

    ! Integration step sizes
    double precision :: du, dv                  ! Step sizes in u, v directions

    ! AMR (Adaptive Mesh Refinement) parameters
    logical :: AMR                              ! Enable AMR
    double precision :: reldiff_max             ! Thresholds for refinement

    ! Output sampling spacing (absolute values in domain units)
    double precision :: output_du               ! Write every this Δu
    double precision :: output_dv               ! Write every this Δv

    ! Progress reporting cadence (stdout)
    integer :: progress_stride
    integer :: progress_header_stride

    ! FIXME: move to physics config?
    ! Initial conditions
    double precision :: m0                      ! Initial mass parameter

    ! Derived values (computed)
    integer :: Nu, Nv                           ! Number of grid points
    integer :: Nu_max                           ! Maximum u-arrays allocation size
  end type simulation_config

contains

!> Initialize simulation configuration with default values
!! Sets reasonable defaults for a mass inflation simulation.
!! Users can override individual fields after calling this.
subroutine init_simulation_config(sim_cfg)
  type(simulation_config), intent(out) :: sim_cfg

  ! Integration domain
  sim_cfg%u_min = 0.0d0
  sim_cfg%v_min = 5.0d0
  sim_cfg%u_max = 30.0d0
  sim_cfg%v_max = 10.0d0

  ! Integration step sizes
  sim_cfg%du = 0.01d0
  sim_cfg%dv = 0.0005d0

  ! AMR parameters
  sim_cfg%AMR = .true.
  sim_cfg%reldiff_max = 0.0001d0

  ! Output sampling spacings
  sim_cfg%output_du = 0.05d0
  sim_cfg%output_dv = 0.05d0

  ! Progress reporting cadence (stdout)
  sim_cfg%progress_stride = 100
  sim_cfg%progress_header_stride = 5000

  ! Initial condition scale
  sim_cfg%m0 = 1.0d0

  ! Compute grid dimensions
  sim_cfg%Nu = int((sim_cfg%u_max - sim_cfg%u_min) / sim_cfg%du + 1.001d0)
  sim_cfg%Nv = int((sim_cfg%v_max - sim_cfg%v_min) / sim_cfg%dv + 1.001d0)

  ! Estimate maximum array allocation size (used for AMR refinement)
  sim_cfg%Nu_max = int(2.0d0 * (sim_cfg%u_max - sim_cfg%u_min) / sim_cfg%reldiff_max)
end subroutine init_simulation_config

end module simulation_config_mod
