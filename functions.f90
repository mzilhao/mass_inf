module physics_config_mod
  implicit none

  !> Physics configuration type - encapsulates all physics parameters
  type :: physics_config
    integer :: neq = 3                          ! Number of equations
    integer :: D = 4                            ! Spacetime dimension
    double precision :: lambda = 0.0d0          ! Cosmological constant
    double precision :: A = 0.0d0               ! Scalar field amplitude
    double precision :: Delta = 1.0d0           ! Scalar field width
    double precision :: q = 0.95d0              ! Electric charge
    double precision :: m0 = 1.0d0              ! Initial mass parameter

    double precision :: q2 = 0.0d0, qq2 = 0.0d0, qq = 0.0d0  ! Derived constants, dummy values
  end type physics_config

end module physics_config_mod

module functions
  use physics_config_mod
  use simulation_config_mod
  implicit none
  double precision, parameter :: PI = 4.0d0 * atan(1.0d0)
  private
  public :: F, init_physics_config, read_physics_config_from_file, init_cond
  public :: compute_diagnostics, write_output, write_output_header

contains

!> Initialize physics configuration with derived constants
subroutine init_physics_config(cfg)
  type(physics_config), intent(inout) :: cfg
  cfg%q2 = cfg%q * cfg%q
  cfg%qq2 = 0.5d0 * cfg%q2 * (cfg%D - 3) * (cfg%D - 2)
  cfg%qq = sqrt(cfg%qq2)
end subroutine init_physics_config

!> Read physics configuration from a namelist file
!! Reads &physics namelist and computes derived constants.
!! Reads from an existing namelist file.
subroutine read_physics_config_from_file(cfg, filename)
  type(physics_config), intent(out) :: cfg
  character(len=*), intent(in)      :: filename

  ! Local variables for namelist reading
  integer :: D
  double precision :: lambda, A, Delta, q, m0
  namelist /physics/ D, lambda, A, Delta, q, m0

  integer :: unit, ierr

  ! Initialize with type defaults
  cfg = physics_config()
  D      = cfg%D
  lambda = cfg%lambda
  A      = cfg%A
  Delta  = cfg%Delta
  q      = cfg%q
  m0     = cfg%m0

  ! Read namelist
  open(newunit=unit, file=filename, status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
    write(*, '(a,a,a)') 'Error: could not open parameter file "', trim(filename), '"'
  else
    read(unit, nml=physics, iostat=ierr)
    if (ierr /= 0) then
      write(*, '(a,a,a)') 'Error: could not read &physics namelist from "', trim(filename), '"'
    end if
    close(unit)
  end if

  ! Update cfg with (possibly modified) namelist values
  cfg%D      = D
  cfg%lambda = lambda
  cfg%A      = A
  cfg%Delta  = Delta
  cfg%q      = q
  cfg%m0     = m0

  ! Compute derived constants
  call init_physics_config(cfg)

end subroutine read_physics_config_from_file

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
!> Initialize boundary conditions at u=u_min and v=v_min
!! Returns: h_u0 (IC along u_min), h_v0 (IC along v_min)
subroutine init_cond(h_u0, h_v0, sim_cfg, cfg)
  implicit none

  double precision, dimension(:,:), allocatable, intent(out) :: h_u0, h_v0
  type(simulation_config), intent(inout) :: sim_cfg
  type(physics_config), intent(in) :: cfg

  ! Local variables
  integer :: i, D
  double precision :: u, v, lambda, q2
  double precision :: r00, sigma_0, m0, ru0, v1
  double precision :: A, Delta
  logical :: scalarfield

  ! Extract config values
  D      = cfg%D
  lambda = cfg%lambda
  q2     = cfg%q2
  A      = cfg%A
  Delta  = cfg%Delta
  m0     = cfg%m0

  ! Allocate boundary condition arrays
  allocate(h_u0(sim_cfg%Nv, cfg%neq))
  allocate(h_v0(sim_cfg%Nu_max, cfg%neq))
  h_u0 = 0.0d0
  h_v0 = 0.0d0

  ! Initial condition parameters
  r00 = sim_cfg%v_min
  sigma_0 = -0.5d0 * log(2.0d0)
  ru0 = 0.25d0 * (2.0d0 / (r00**(D-3)) * (m0 - q2/(2.0d0*r00**(D-3))) &
       - 1.0d0 + lambda * r00 * r00 / 3.0d0)

  ! FIXME
  scalarfield = .true.
  v1 = sim_cfg%v_min + Delta

  ! Boundary conditions at u = u_min
  do i = 1, sim_cfg%Nv
    v = sim_cfg%v_min + (i-1) * sim_cfg%dv
    h_u0(i, 1) = v  ! r(u_min,v)

    if (scalarfield) then
      if (v <= sim_cfg%v_min + Delta) then
        ! Perturbation phase
        h_u0(i, 2) = A / (4*Pi) * (2*Pi*(v - sim_cfg%v_min) - Delta * sin(2.0d0*Pi*(v - sim_cfg%v_min)/Delta))
        h_u0(i, 3) = sigma_0 + 2.0d0/(D-2.0d0) * A*A/(256*Pi*Pi) &
             * (15*sim_cfg%v_min**2 - 24*Pi*Pi*sim_cfg%v_min**2 - 30*sim_cfg%v_min*v1 &
             + 15*v1**2 + 24*Pi*Pi*v**2 &
             - 16*Delta**2*cos(2*Pi*(sim_cfg%v_min-v)/Delta) &
             + Delta**2*cos(4*Pi*(sim_cfg%v_min-v)/Delta) &
             - 32*Pi*sim_cfg%v_min*v*sin(2*Pi*(sim_cfg%v_min-v)/Delta) &
             + 32*Pi*v1*v*sin(2*Pi*(sim_cfg%v_min-v)/Delta) &
             + 4*Pi*sim_cfg%v_min*v*sin(4*Pi*(sim_cfg%v_min-v)/Delta) &
             - 4*Pi*v1*v*sin(4*Pi*(sim_cfg%v_min-v)/Delta))
      else
        ! No perturbation phase
        h_u0(i, 2) = 0.5d0 * A * Delta
        h_u0(i, 3) = sigma_0 + 2.0d0/(D-2.0d0) * 3*A*A*Delta*(sim_cfg%v_min+v1)/32.0d0
      end if
    else
      ! No scalar field
      h_u0(i, 2) = 0.0d0
      h_u0(i, 3) = sigma_0
    end if
  end do

  ! Boundary conditions at v = v_min
  do i = 1, sim_cfg%Nu
    u = sim_cfg%u_min + (i-1) * sim_cfg%du
    h_v0(i, 1) = r00 + u * ru0  ! r(u,v_min)
    h_v0(i, 2) = 0.0d0          ! phi(u,v_min)
    h_v0(i, 3) = sigma_0        ! sigma(u,v_min)
  end do
end subroutine init_cond

!====================================================================================
!> Initialize output file with physics model header
subroutine write_output_header(output_unit, cfg, sim_cfg)
  integer, intent(in) :: output_unit
  type(physics_config), intent(in) :: cfg
  type(simulation_config), intent(in) :: sim_cfg

  ! Physics parameters
  write(output_unit, '(a,i2)')    '# D         = ', cfg%D
  write(output_unit, '(a,g10.4)') '# lambda    = ', cfg%lambda
  write(output_unit, '(a,g10.4)') '# q         = ', sqrt(cfg%q2)
  write(output_unit, '(a,g10.4)') '# A         = ', cfg%A
  write(output_unit, '(a,g10.4)') '# Delta     = ', cfg%Delta
  ! FIXME: hardcoded initial conditions info
  write(output_unit, '(a,g10.4)') '# sigma_0   = ', -0.5d0 * log(2.0d0)
  write(output_unit, '(a,g10.4)') '# m0        = ', 1.0d0

  ! Grid parameters
  write(output_unit, '(a)') '#'
  write(output_unit, '(a,g14.6)') '# u_min     = ', sim_cfg%u_min
  write(output_unit, '(a,g14.6)') '# u_max     = ', sim_cfg%u_max
  write(output_unit, '(a,g14.6)') '# v_min     = ', sim_cfg%v_min
  write(output_unit, '(a,g14.6)') '# v_max     = ', sim_cfg%v_max
  write(output_unit, '(a,g14.6)') '# du        = ', sim_cfg%du
  write(output_unit, '(a,g14.6)') '# dv        = ', sim_cfg%dv
  write(output_unit, '(a,g14.6)') '# output_du = ', sim_cfg%output_du
  write(output_unit, '(a,g14.6)') '# output_dv = ', sim_cfg%output_dv
  write(output_unit, '(a)') '#'

  ! Column headers
  write(output_unit, '(a)') '# Columns: u, r, phi, sigma, mass, drdv, Ricci'
end subroutine write_output_header

!====================================================================================
!> Compute diagnostics and write them if output cadence is met
!!
!! Writes columnar ASCII output: u, r, phi, sigma, mass, drdv, Ricci
!! v-slices are marked with comment lines: # v = X.XXXXX
subroutine write_output(output_unit, u_val, v_val, h_N, h_S, h_E, h_W, du, dv, sim_cfg, cfg)
  integer, intent(in)                 :: output_unit
  double precision, intent(in)        :: u_val, v_val, du, dv
  double precision, dimension(:), intent(in) :: h_N, h_S, h_E, h_W
  type(simulation_config), intent(in) :: sim_cfg
  type(physics_config),    intent(in) :: cfg

  integer :: neq_local
  double precision :: temp
  logical :: u_ok, v_ok
  double precision, dimension(:), allocatable :: h_P, dhdu_P, dhdv_P, dhduv_P
  double precision :: u_P, v_P, mass, drdv, ricci

  double precision, save :: last_v_marked_val = -1.0d99

  ! Check output condition first, before doing any work
  ! Output condition: write when current (u,v) aligns with sampling spacings.
  ! We check if (u - u_min)/output_du and (v - v_min)/output_dv are near integers.
  ! This is robust to floating-point drift and local AMR changes.
  if (sim_cfg%output_du > 0.0d0) then
    temp = abs((u_val - sim_cfg%u_min) / sim_cfg%output_du)
    u_ok = abs(temp - dnint(temp)) < 1.0d-7
  else
    u_ok = .true.
  end if

  if (sim_cfg%output_dv > 0.0d0) then
    temp = abs((v_val - sim_cfg%v_min) / sim_cfg%output_dv)
    v_ok = abs(temp - dnint(temp)) < 1.0d-6
  else
    v_ok = .true.
  end if

  if (.not. (u_ok .and. v_ok)) return

  ! The diagnostics are computed at the midpoint P = (u+du/2, v+dv/2)
  u_P = u_val + 0.5d0 * du
  v_P = v_val + 0.5d0 * dv

  ! Write a new v-slice block if v_val is more than one half output_dv away from last marked
  if (abs(v_val - last_v_marked_val) > 0.5d0 * sim_cfg%output_dv) then
    write(output_unit, '(a)')
    write(output_unit, '(a,f10.6)') '# v = ', v_P
    last_v_marked_val = v_val
  end if

  neq_local = size(h_N)
  allocate(h_P(neq_local), dhdu_P(neq_local), dhdv_P(neq_local), dhduv_P(neq_local))

  h_P     = 0.25d0 * (h_N + h_S + h_E - h_W)
  dhdu_P  = (h_W - h_S + h_N - h_E) * 0.5d0 / du
  dhdv_P  = (h_E - h_S + h_N - h_W) * 0.5d0 / dv

  call F(dhduv_P, h_P, dhdu_P, dhdv_P, neq_local, cfg)
  call compute_diagnostics(h_P, dhdu_P, dhdv_P, dhduv_P, cfg, mass, drdv, ricci)

  ! Write columnar output: u, r, phi, sigma, mass, drdv, Ricci
  write(output_unit, '(7e16.8)') u_P, h_P(1), h_P(2), h_P(3), mass, drdv, ricci

  deallocate(h_P, dhdu_P, dhdv_P, dhduv_P)
end subroutine write_output


!====================================================================================
end module functions
