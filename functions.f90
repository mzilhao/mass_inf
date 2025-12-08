module physics_config_mod
  implicit none

  !> Physics configuration type - encapsulates all physics parameters
  type :: physics_config
    integer :: D = 4                            ! Spacetime dimension
    integer :: neq = 3                          ! Number of equations
    double precision :: lambda = 0.0d0          ! Cosmological constant
    double precision :: q = 0.95d0              ! Electric charge
    double precision :: q2, qq2, qq            ! Derived constants
    double precision :: Pi = 3.1415926535897932384626433d0
  end type physics_config

  !> Simulation configuration type - grid and integration parameters
  type :: simulation_config
    double precision :: u0, v0                  ! Integration domain start
    double precision :: uf, vf                  ! Integration domain end
    double precision :: Du, Dv                  ! Integration step sizes
    double precision :: gradmax, gradmin        ! AMR gradient thresholds
    logical :: AMR                              ! Adaptive mesh refinement enabled
    integer :: resu, resv                       ! Output resolution
    double precision :: m0                      ! Initial mass
    integer :: Nu, Nv                           ! Grid points computed
    integer :: big_dim                          ! Array allocation size
  end type simulation_config

end module physics_config_mod

module functions
  use physics_config_mod
  implicit none
  private
  public :: F, init_physics_config, init_simulation_config, init_cond, print_simulation_header

contains

  !> Initialize physics configuration with derived constants
  subroutine init_physics_config(cfg)
    type(physics_config), intent(inout) :: cfg
    cfg%q2 = cfg%q * cfg%q
    cfg%qq2 = 0.5d0 * cfg%q2 * (cfg%D - 3) * (cfg%D - 2)
    cfg%qq = sqrt(cfg%qq2)
  end subroutine init_physics_config

  !> Initialize simulation configuration with default values
  subroutine init_simulation_config(sim_cfg, physics_cfg)
    type(simulation_config), intent(out) :: sim_cfg
    type(physics_config), intent(in) :: physics_cfg

    ! Integration domain
    sim_cfg%u0 = 0.0d0
    sim_cfg%v0 = 5.0d0
    sim_cfg%uf = 30.0d0
    sim_cfg%vf = 15.0d0

    ! Integration step sizes
    sim_cfg%Du = 0.01d0
    sim_cfg%Dv = 0.0005d0

    ! AMR parameters
    sim_cfg%AMR = .true.
    sim_cfg%gradmax = 0.0001d0
    sim_cfg%gradmin = 0.1d0  ! Currently unused

    ! Output resolution
    sim_cfg%resu = 20
    sim_cfg%resv = 20

    ! Initial conditions
    sim_cfg%m0 = 1.0d0

    ! Compute grid points
    sim_cfg%Nu = int((sim_cfg%uf - sim_cfg%u0) / sim_cfg%Du + 1.001d0)
    sim_cfg%Nv = int((sim_cfg%vf - sim_cfg%v0) / sim_cfg%Dv + 1.001d0)

    ! Allocate array size
    sim_cfg%big_dim = int(2.0d0 * (sim_cfg%uf - sim_cfg%u0) / sim_cfg%gradmax)
  end subroutine init_simulation_config

  !> Print simulation header to output file
  subroutine print_simulation_header(id, physics_cfg, A)
    integer, intent(in) :: id
    type(physics_config), intent(in) :: physics_cfg
    double precision, intent(in) :: A

    write(id, '(a,i2)')    '# D         = ', physics_cfg%D
    write(id, '(a,g10.4)') '# lambda    = ', physics_cfg%lambda
    write(id, '(a,g10.4)') '# q         = ', sqrt(physics_cfg%q2)
    write(id, '(a,g10.4)') '# A         = ', A
    write(id, '(a,g10.4)') '# sigma_0   = ', -0.5d0 * log(2.0d0)
    write(id, '(a,g10.4)') '# m0        = ', 1.0d0
    write(id, '(a)') '#'
  end subroutine print_simulation_header

  !====================================================================================
  !> Right-hand side of the PDE system
  !! Inputs: h (field values), dhdu, dhdv (derivatives)
  !! Output: dhduv (mixed derivatives)
  subroutine F(dhduv, h, dhdu, dhdv, neq, cfg)
    implicit none

    integer, intent(in) :: neq
    double precision, dimension(neq), intent(out) :: dhduv
    double precision, dimension(neq), intent(in)  :: h, dhdu, dhdv
    type(physics_config), intent(in) :: cfg

    ! Local copies for readability
    integer :: D
    double precision :: lambda, qq2

    if (neq /= cfg%neq) error stop "F: neq mismatch with config"

    ! Extract config values for cleaner code
    D = cfg%D
    lambda = cfg%lambda
    qq2 = cfg%qq2

    ! h(1) = r, h(2) = phi, h(3) = sigma
    dhduv(1) = qq2/(D-2) * exp(2*h(3)) / (h(1)**(2*D-5)) &
         + (D-1)/6.0d0 * lambda * h(1) * exp(2*h(3)) &
         - (D-3)/2.0d0 * exp(2*h(3)) / h(1) &
         - (D-3) * dhdu(1) * dhdv(1) / h(1)

    dhduv(2) = -(D-2) / (2.0d0*h(1)) * (dhdv(1)*dhdu(2) + dhdu(1)*dhdv(2))

    dhduv(3) = -dhdu(2) * dhdv(2) &
         - (3*D-8) / (2.0d0*D-4) * qq2 * exp(2*h(3)) / (h(1)**(2*(D-2))) &
         - (D-4)*(D-1) * lambda / 12.0d0 * exp(2*h(3)) &
         + (D-3)*(D-2) / 4.0d0 * exp(2*h(3)) / (h(1)*h(1)) &
         + (D-3)*(D-2) / 2.0d0 * dhdu(1) * dhdv(1) / (h(1)*h(1))
  end subroutine F

  !====================================================================================
  !> Initialize boundary conditions at u=u0 and v=v0
  !! Returns: h_u0 (IC along u0), h_v0 (IC along v0)
  subroutine init_cond(h_u0, h_v0, sim_cfg, physics_cfg)
    implicit none

    double precision, dimension(:,:), allocatable, intent(out) :: h_u0, h_v0
    type(simulation_config), intent(inout) :: sim_cfg
    type(physics_config), intent(in) :: physics_cfg

    ! Local variables
    integer :: i, D
    double precision :: u, v, Pi, lambda, q2
    double precision :: r00, sigma_0, m0, ru0, v1
    double precision :: A, Delta
    logical :: scalarfield

    ! Extract config values
    D = physics_cfg%D
    Pi = physics_cfg%Pi
    lambda = physics_cfg%lambda
    q2 = physics_cfg%q2

    ! Allocate boundary condition arrays
    allocate(h_u0(physics_cfg%neq, sim_cfg%big_dim))
    allocate(h_v0(physics_cfg%neq, sim_cfg%big_dim))
    h_u0 = 0.0d0
    h_v0 = 0.0d0

    ! Initial condition parameters
    r00 = sim_cfg%v0
    sigma_0 = -0.5d0 * log(2.0d0)
    m0 = sim_cfg%m0
    ru0 = 0.25d0 * (2.0d0 / (r00**(D-3)) * (m0 - q2/(2.0d0*r00**(D-3))) &
         - 1.0d0 + lambda * r00 * r00 / 3.0d0)

    ! Perturbation parameters
    A = 0.0d0
    Delta = 1.0d0
    scalarfield = .true.
    v1 = sim_cfg%v0 + Delta

    ! Boundary conditions at u = u0
    do i = 1, sim_cfg%Nv
      v = sim_cfg%v0 + (i-1) * sim_cfg%Dv
      h_u0(1, i) = v  ! r(u0,v)

      if (scalarfield) then
        if (v <= sim_cfg%v0 + Delta) then
          ! Perturbation phase
          h_u0(2, i) = A / (4*Pi) * (2*Pi*(v - sim_cfg%v0) - Delta * sin(2.0d0*Pi*(v - sim_cfg%v0)/Delta))
          h_u0(3, i) = sigma_0 + 2.0d0/(D-2.0d0) * A*A/(256*Pi*Pi) &
               * (15*sim_cfg%v0**2 - 24*Pi*Pi*sim_cfg%v0**2 - 30*sim_cfg%v0*v1 &
               + 15*v1**2 + 24*Pi*Pi*v**2 &
               - 16*Delta**2*cos(2*Pi*(sim_cfg%v0-v)/Delta) &
               + Delta**2*cos(4*Pi*(sim_cfg%v0-v)/Delta) &
               - 32*Pi*sim_cfg%v0*v*sin(2*Pi*(sim_cfg%v0-v)/Delta) &
               + 32*Pi*v1*v*sin(2*Pi*(sim_cfg%v0-v)/Delta) &
               + 4*Pi*sim_cfg%v0*v*sin(4*Pi*(sim_cfg%v0-v)/Delta) &
               - 4*Pi*v1*v*sin(4*Pi*(sim_cfg%v0-v)/Delta))
        else
          ! No perturbation phase
          h_u0(2, i) = 0.5d0 * A * Delta
          h_u0(3, i) = sigma_0 + 2.0d0/(D-2.0d0) * 3*A*A*Delta*(sim_cfg%v0+v1)/32.0d0
        end if
      else
        ! No scalar field
        h_u0(2, i) = 0.0d0
        h_u0(3, i) = sigma_0
      end if
    end do

    ! Boundary conditions at v = v0
    do i = 1, sim_cfg%Nu
      u = sim_cfg%u0 + (i-1) * sim_cfg%Du
      h_v0(1, i) = r00 + u * ru0  ! r(u,v0)
      h_v0(2, i) = 0.0d0          ! phi(u,v0)
      h_v0(3, i) = sigma_0        ! sigma(u,v0)
    end do
  end subroutine init_cond

  !====================================================================================
end module functions
