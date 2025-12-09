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

end module physics_config_mod

module functions
  use physics_config_mod
  use simulation_config_mod
  implicit none
  private
  public :: F, init_physics_config, init_cond, print_simulation_header
  public :: compute_diagnostics, write_output_if_needed, write_output_header, write_output_separator

contains

  !> Initialize physics configuration with derived constants
  subroutine init_physics_config(cfg)
    type(physics_config), intent(inout) :: cfg
    cfg%q2 = cfg%q * cfg%q
    cfg%qq2 = 0.5d0 * cfg%q2 * (cfg%D - 3) * (cfg%D - 2)
    cfg%qq = sqrt(cfg%qq2)
  end subroutine init_physics_config

  !> Compute model-dependent diagnostics (mass, drdv, Ricci)
  subroutine compute_diagnostics(h, dhdu, dhdv, dhduv, cfg, mass, drdv, ricci)
    implicit none
    double precision, dimension(:), intent(in)  :: h, dhdu, dhdv, dhduv
    type(physics_config), intent(in)            :: cfg
    double precision,           intent(out)     :: mass, drdv, ricci

    integer :: D
    double precision :: lambda, q2

    D = cfg%D
    lambda = cfg%lambda
    q2 = cfg%q2

    drdv = dhdv(1)

    mass = 0.5d0 * h(1)**(D-3) * ( 1.d0 - lambda/3.d0 * h(1)*h(1)           &
          + q2 / ( (h(1)*h(1))**(D-3) )                                     &
          + 2.d0 * dhdu(1) * dhdv(1) / exp(2.d0 * h(3)) )

    ricci = (D-3)*(D-2)/( h(1)*h(1) ) * (                                   &
          1 + 2*exp(-2*h(3)) * dhdu(1)*dhdv(1)                             &
          )                                                                &
          + 4*exp(-2*h(3))*dhduv(3) + 4*(D-2)*dhduv(1)*exp(-2*h(3)) / h(1)
  end subroutine compute_diagnostics

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
    allocate(h_u0(sim_cfg%Nv, physics_cfg%neq))
    allocate(h_v0(sim_cfg%Nu_max, physics_cfg%neq))
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
      v = sim_cfg%v0 + (i-1) * sim_cfg%dv
      h_u0(i, 1) = v  ! r(u0,v)

      if (scalarfield) then
        if (v <= sim_cfg%v0 + Delta) then
          ! Perturbation phase
          h_u0(i, 2) = A / (4*Pi) * (2*Pi*(v - sim_cfg%v0) - Delta * sin(2.0d0*Pi*(v - sim_cfg%v0)/Delta))
          h_u0(i, 3) = sigma_0 + 2.0d0/(D-2.0d0) * A*A/(256*Pi*Pi) &
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
          h_u0(i, 2) = 0.5d0 * A * Delta
          h_u0(i, 3) = sigma_0 + 2.0d0/(D-2.0d0) * 3*A*A*Delta*(sim_cfg%v0+v1)/32.0d0
        end if
      else
        ! No scalar field
        h_u0(i, 2) = 0.0d0
        h_u0(i, 3) = sigma_0
      end if
    end do

    ! Boundary conditions at v = v0
    do i = 1, sim_cfg%Nu
      u = sim_cfg%u0 + (i-1) * sim_cfg%du
      h_v0(i, 1) = r00 + u * ru0  ! r(u,v0)
      h_v0(i, 2) = 0.0d0          ! phi(u,v0)
      h_v0(i, 3) = sigma_0        ! sigma(u,v0)
    end do
  end subroutine init_cond

  !====================================================================================
  !> Initialize output file with physics model header
  subroutine write_output_header(output_unit, cfg)
    integer, intent(in) :: output_unit
    type(physics_config), intent(in) :: cfg

    call print_simulation_header(output_unit, cfg, 0.0d0)
    write(output_unit,'(a)') '# | u | v | r | phi | sigma | mass | drdv | Ricci'
  end subroutine write_output_header

  !====================================================================================
  !> Compute diagnostics and write them if output cadence is met
  !!
  !! This routine handles both diagnostic computation and output filtering.
  !! It encapsulates which quantities are physically meaningful to output,
  !! making it easy to swap for a different physics model.
  subroutine write_output_if_needed(output_unit, u_val, v_val, h_N, h_S, h_E, h_W, du, dv, &
                                     u0, v0, sim_cfg, cfg, OUTPUT_TOL_U, OUTPUT_TOL_V)
    integer, intent(in)                 :: output_unit
    double precision, intent(in)        :: u_val, v_val, du, dv, u0, v0
    double precision, dimension(:), intent(in) :: h_N, h_S, h_E, h_W
    type(simulation_config), intent(in) :: sim_cfg
    type(physics_config),    intent(in) :: cfg
    double precision, intent(in)        :: OUTPUT_TOL_U, OUTPUT_TOL_V

    integer :: neq_local
    double precision :: tempu, tempv
    double precision, dimension(:), allocatable :: h_P, dhdu_P, dhdv_P, dhduv_P
    double precision :: mass, drdv, ricci

    neq_local = size(h_N)
    allocate(h_P(neq_local), dhdu_P(neq_local), dhdv_P(neq_local), dhduv_P(neq_local))

    h_P     = 0.25d0 * (h_N + h_S + h_E - h_W)
    dhdu_P  = (h_W - h_S + h_N - h_E) * 0.5d0 / du
    dhdv_P  = (h_E - h_S + h_N - h_W) * 0.5d0 / dv

    call F(dhduv_P, h_P, dhdu_P, dhdv_P, neq_local, cfg)
    call compute_diagnostics(h_P, dhdu_P, dhdv_P, dhduv_P, cfg, mass, drdv, ricci)

    tempv = abs(v_val - v0) * sim_cfg%resv
    tempu = abs(u_val - u0) * sim_cfg%resu
    if (abs(tempu - int(tempu + du*0.5d0)) < OUTPUT_TOL_U .and. &
        abs(tempv - int(tempv + dv*0.5d0)) < OUTPUT_TOL_V) then
      write(output_unit,*) (/ u_val, v_val, h_N, (/ mass, drdv, ricci /) /)
    end if

    deallocate(h_P, dhdu_P, dhdv_P, dhduv_P)
  end subroutine write_output_if_needed

  !====================================================================================
  !> Write output separator (blank line) for gnuplot block separation
  !! This is a physics model choice: some output formats may need separators
  !! between data blocks for proper parsing.
  subroutine write_output_separator(output_unit)
    integer, intent(in) :: output_unit
    write(output_unit,'(a)') ''
  end subroutine write_output_separator

  !====================================================================================
end module functions
