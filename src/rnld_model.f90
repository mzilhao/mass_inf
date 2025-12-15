module model_config_mod
  use precision
  implicit none

  !> Model configuration type - encapsulates all physics parameters
  type :: model_config
    integer  :: D = 4                    ! Spacetime dimension
    real(dp) :: lambda = 0.0_dp          ! Cosmological constant
    real(dp) :: A = 0.0_dp               ! Scalar field amplitude
    real(dp) :: Delta = 1.0_dp           ! Scalar field width
    real(dp) :: q = 0.95_dp              ! Electric charge
    real(dp) :: m0 = 1.0_dp              ! Initial mass parameter

    real(dp) :: q2 = 0.0_dp, qq2 = 0.0_dp, qq = 0.0_dp  ! Derived constants, dummy values
  end type model_config

end module model_config_mod

module rnld_model
  use precision
  use model_config_mod
  use grid_config_mod
  implicit none
  real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)
  integer, parameter :: NEQ = 3      ! Number of equations/fields: r, phi, sigma
  integer, save :: diag_unit = -1, fields_unit = -1, derivs_unit = -1
  logical, save :: diag_open = .false., fields_open = .false., derivs_open = .false.

  private
  public :: F, init_model_config, read_model_config_from_file, init_cond
  public :: compute_diagnostics, write_output
  public :: open_output_files, close_output_files
  public :: NEQ

contains

!> Initialize physics configuration with derived constants
subroutine init_model_config(model_cfg)
  type(model_config), intent(inout) :: model_cfg
  model_cfg%q2 = model_cfg%q * model_cfg%q
  model_cfg%qq2 = 0.5d0 * model_cfg%q2 * (model_cfg%D - 3) * (model_cfg%D - 2)
  model_cfg%qq = sqrt(model_cfg%qq2)
end subroutine init_model_config

!> Read physics configuration from a namelist file
!! Reads &physics namelist and computes derived constants.
!! Reads from an existing namelist file.
subroutine read_model_config_from_file(model_cfg, filename)
  type(model_config), intent(out) :: model_cfg
  character(len=*), intent(in)      :: filename

  ! Local variables for namelist reading
  integer :: D
  real(dp) :: lambda, A, Delta, q, m0
  namelist /physics/ D, lambda, A, Delta, q, m0

  integer :: unit, ierr

  ! Initialize with type defaults
  model_cfg = model_config()
  D      = model_cfg%D
  lambda = model_cfg%lambda
  A      = model_cfg%A
  Delta  = model_cfg%Delta
  q      = model_cfg%q
  m0     = model_cfg%m0

  ! Read namelist
  open(newunit=unit, file=filename, status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
    write(*, '(a,a,a)') 'Error: could not open parameter file "', trim(filename), '"'
    call exit(1)
  end if

  read(unit, nml=physics, iostat=ierr)
  if (ierr /= 0) then
    write(*, '(a,a,a)') 'Error: could not read &physics namelist from "', trim(filename), '"'
    close(unit)
    call exit(1)
  end if
  close(unit)

  ! Update cfg with (possibly modified) namelist values
  model_cfg%D      = D
  model_cfg%lambda = lambda
  model_cfg%A      = A
  model_cfg%Delta  = Delta
  model_cfg%q      = q
  model_cfg%m0     = m0

  ! Compute derived constants
  call init_model_config(model_cfg)
end subroutine read_model_config_from_file

!> Compute model-dependent diagnostics (mass, Ricci)
subroutine compute_diagnostics(mass, ricci, h, dhdu, dhdv, dhduv, model_cfg)
  implicit none
  real(dp), intent(out)               :: mass, ricci
  real(dp), dimension(:), intent(in)  :: h, dhdu, dhdv, dhduv
  type(model_config), intent(in)    :: model_cfg

  integer :: D
  real(dp) :: lambda, q2

  D      = model_cfg%D
  lambda = model_cfg%lambda
  q2     = model_cfg%q2

  mass = 0.5d0 * h(1)**(D-3) * ( 1.d0 - lambda/3.d0 * h(1)*h(1)           &
        + q2 / ( (h(1)*h(1))**(D-3) )                                     &
        + 2.d0 * dhdu(1) * dhdv(1) / exp(2.d0 * h(3)) )

  ricci = (D-3)*(D-2)/( h(1)*h(1) ) * (                                   &
        1 + 2*exp(-2*h(3)) * dhdu(1)*dhdv(1)                              &
        )                                                                 &
        + 4*exp(-2*h(3))*dhduv(3) + 4*(D-2)*dhduv(1)*exp(-2*h(3)) / h(1)
end subroutine compute_diagnostics

!====================================================================================
!> Right-hand side of the PDE system
!! Inputs: h (field values), dhdu, dhdv (derivatives)
!! Output: dhduv (mixed derivatives)
subroutine F(dhduv, h, dhdu, dhdv, model_cfg)
  implicit none

  real(dp), dimension(NEQ), intent(out) :: dhduv
  real(dp), dimension(NEQ), intent(in)  :: h, dhdu, dhdv
  type(model_config), intent(in)      :: model_cfg

  ! Local copies for readability
  integer :: D
  real(dp) :: lambda, qq2

  ! Extract config values for cleaner code
  D      = model_cfg%D
  lambda = model_cfg%lambda
  qq2    = model_cfg%qq2

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
subroutine init_cond(h_u0, h_v0, grid_cfg, model_cfg)
  implicit none

  real(dp), dimension(:,:), intent(inout) :: h_u0, h_v0
  type(grid_config),  intent(in)          :: grid_cfg
  type(model_config), intent(in)          :: model_cfg

  ! Local variables
  integer :: i, D
  real(dp) :: u, v, lambda, q2
  real(dp) :: r00, sigma_0, m0, ru0, v0, v1
  real(dp) :: A, Delta
  logical :: scalarfield

  ! Extract config values
  D      = model_cfg%D
  lambda = model_cfg%lambda
  q2     = model_cfg%q2
  A      = model_cfg%A
  Delta  = model_cfg%Delta
  m0     = model_cfg%m0
  v0     = grid_cfg%v_min

  ! Initial condition parameters
  r00 = v0
  sigma_0 = -0.5d0 * log(2.0d0)
  ru0 = 0.25d0 * (2.0d0 / (r00**(D-3)) * (m0 - q2/(2.0d0*r00**(D-3))) &
       - 1.0d0 + lambda * r00 * r00 / 3.0d0)

  ! FIXME
  scalarfield = .true.
  v1 = v0 + Delta

  ! Boundary conditions at u = u_min
  do i = 1, grid_cfg%Nv
    v = v0 + (i-1) * grid_cfg%dv
    h_u0(i, 1) = v  ! r(u_min,v)

    if (scalarfield) then
      if (v <= v0 + Delta) then
        ! Perturbation phase
        h_u0(i, 2) = A / (4*Pi) * (2*Pi*(v - v0) - Delta * sin(2.0d0*Pi*(v - v0)/Delta))
        h_u0(i, 3) = sigma_0 + 2.0d0/(D-2.0d0) * A*A/(256*Pi*Pi) &
             * (15*v0**2 - 24*Pi*Pi*v0**2 - 30*v0*v1 &
             + 15*v1**2 + 24*Pi*Pi*v**2 &
             - 16*Delta**2*cos(2*Pi*(v0-v)/Delta) &
             + Delta**2*cos(4*Pi*(v0-v)/Delta) &
             - 32*Pi*v0*v*sin(2*Pi*(v0-v)/Delta) &
             + 32*Pi*v1*v*sin(2*Pi*(v0-v)/Delta) &
             + 4*Pi*v0*v*sin(4*Pi*(v0-v)/Delta) &
             - 4*Pi*v1*v*sin(4*Pi*(v0-v)/Delta))
      else
        ! No perturbation phase
        h_u0(i, 2) = 0.5d0 * A * Delta
        h_u0(i, 3) = sigma_0 + 2.0d0/(D-2.0d0) * 3*A*A*Delta*(v0+v1)/32.0d0
      end if
    else
      ! No scalar field
      h_u0(i, 2) = 0.0d0
      h_u0(i, 3) = sigma_0
    end if
  end do

  ! Boundary conditions at v = v_min
  do i = 1, grid_cfg%Nu
    u = grid_cfg%u_min + (i-1) * grid_cfg%du
    h_v0(i, 1) = r00 + u * ru0  ! r(u,v_min)
    h_v0(i, 2) = 0.0d0          ! phi(u,v_min)
    h_v0(i, 3) = sigma_0        ! sigma(u,v_min)
  end do
end subroutine init_cond

!====================================================================================
!> Open all output files
subroutine open_output_files(out_dir)
  character(len=*), intent(in)      :: out_dir

  if (diag_open .or. fields_open .or. derivs_open) call close_output_files()

  open(newunit=fields_unit, file=trim(out_dir)//'/fields.dat', status='replace')
  open(newunit=diag_unit, file=trim(out_dir)//'/diagnostics.dat', status='replace')
  open(newunit=derivs_unit, file=trim(out_dir)//'/derivatives.dat', status='replace')

  write(fields_unit, '(a)') '# Columns: u, r, phi, sigma'
  write(diag_unit,   '(a)') '# Columns: u, mass, Ricci'
  write(derivs_unit, '(a)') '# Columns: u, drdu, drdv'

  fields_open = .true.
  diag_open   = .true.
  derivs_open = .true.
end subroutine open_output_files

!> Close all output files if open
subroutine close_output_files()
  if (fields_open .and. fields_unit > 0) then
    close(fields_unit)
  end if
  if (diag_open .and. diag_unit > 0) then
    close(diag_unit)
  end if
  if (derivs_open .and. derivs_unit > 0) then
    close(derivs_unit)
  end if
  fields_unit = -1
  fields_open = .false.
  diag_unit   = -1
  diag_open   = .false.
  derivs_unit = -1
  derivs_open = .false.
end subroutine close_output_files

!====================================================================================
!> Compute diagnostics and write them if output cadence is met
!!
!! Writes ASCII output in columns.
!! v-slices are marked with comment lines: # v = X.XXXXX
subroutine write_output(u_val, v_val, h_N, h_S, h_E, h_W, du, dv, grid_cfg, model_cfg)
  real(dp), intent(in)                :: u_val, v_val, du, dv
  real(dp), dimension(:), intent(in)  :: h_N, h_S, h_E, h_W
  type(grid_config), intent(in)       :: grid_cfg
  type(model_config),    intent(in)   :: model_cfg

  real(dp), dimension(NEQ) :: h_P, dhdu_P, dhdv_P, dhduv_P
  logical  :: u_ok, v_ok
  real(dp) :: temp, u_P, v_P, mass, ricci

  real(dp), save :: last_v_marked_val = -1.0e99_dp

  if (.not. (fields_open .and. diag_open .and. derivs_open)) error stop 'write_output: output files not open'

  ! Check output condition first, before doing any work
  ! Output condition: write when current (u,v) aligns with sampling spacings.
  ! We check if (u - u_min)/output_du and (v - v_min)/output_dv are near integers.
  ! This is robust to floating-point drift and local AMR changes.
  if (grid_cfg%output_du > 0.0_dp) then
    temp = abs((u_val - grid_cfg%u_min) / grid_cfg%output_du)
    u_ok = abs(temp - nint(temp)) < 1.0e-7_dp
  else
    u_ok = .true.
  end if

  if (grid_cfg%output_dv > 0.0_dp) then
    temp = abs((v_val - grid_cfg%v_min) / grid_cfg%output_dv)
    v_ok = abs(temp - nint(temp)) < 1.0e-6_dp
  else
    v_ok = .true.
  end if

  if (.not. (u_ok .and. v_ok)) return

  ! The diagnostics are computed at the midpoint P = (u+du/2, v+dv/2)
  u_P = u_val + 0.5_dp * du
  v_P = v_val + 0.5_dp * dv

  ! Write a new v-slice block if v_val is more than one half output_dv away from last marked
  if (abs(v_val - last_v_marked_val) > 0.5_dp * grid_cfg%output_dv) then
    write(fields_unit, '(a)')
    write(fields_unit, '(a,f10.6)') '# v = ', v_P
    write(diag_unit,   '(a)')
    write(diag_unit,   '(a,f10.6)') '# v = ', v_P
    write(derivs_unit, '(a)')
    write(derivs_unit, '(a,f10.6)') '# v = ', v_P
    last_v_marked_val = v_val
  end if

  h_P     = 0.25_dp * (h_N + h_S + h_E - h_W)
  dhdu_P  = (h_W - h_S + h_N - h_E) * 0.5_dp / du
  dhdv_P  = (h_E - h_S + h_N - h_W) * 0.5_dp / dv

  call F(dhduv_P, h_P, dhdu_P, dhdv_P, model_cfg)
  call compute_diagnostics(mass, ricci, h_P, dhdu_P, dhdv_P, dhduv_P, model_cfg)

  ! Write columnar output
  write(fields_unit, '(7e16.8)') u_P, h_P(1), h_P(2), h_P(3) ! u, r, phi, sigma
  write(diag_unit,   '(7e16.8)') u_P, mass, ricci            ! u, mass, Ricci
  write(derivs_unit, '(7e16.8)') u_P, dhdu_P(1), dhdv_P(1)   ! u, drdu, drdv

end subroutine write_output

end module rnld_model
