! Template model implementation
! Copy this file to create a new model and replace TODOs.
! Implement the stable API expected by callers:
! - module model_config_mod: defines type(model_config), init_model_config, load
! - module model_mod: provides NEQ and routines:
!     F, init_cond, compute_diagnostics,
!     open_output_files, write_output, close_output_files

module model_config_mod
  use precision_mod
  implicit none

  ! TODO: customize configuration fields for your physics
  type :: model_config
    integer  :: D = 4
    real(dp) :: lambda = 0.0_dp
    real(dp) :: A = 0.0_dp
    real(dp) :: Delta = 0.0_dp
    real(dp) :: q = 0.0_dp
    real(dp) :: m0 = 0.0_dp
    ! Derived constants (example placeholders)
    real(dp) :: q2 = 0.0_dp, qq2 = 0.0_dp, qq = 0.0_dp
  end type model_config

contains

! Derive any constants from input config
subroutine init_model_config(model_cfg)
  type(model_config), intent(inout) :: model_cfg
  ! TODO: derive constants (examples below)
  model_cfg%q2  = model_cfg%q * model_cfg%q
  model_cfg%qq2 = 0.0_dp
  model_cfg%qq  = 0.0_dp
end subroutine init_model_config

! Read &physics namelist and initialize model_cfg
subroutine load(model_cfg, filename)
  type(model_config), intent(out) :: model_cfg
  character(len=*), intent(in)    :: filename
  integer :: unit, ierr
  integer :: D
  real(dp) :: lambda, A, Delta, q, m0
  namelist /physics/ D, lambda, A, Delta, q, m0

  model_cfg = model_config()
  D      = model_cfg%D
  lambda = model_cfg%lambda
  A      = model_cfg%A
  Delta  = model_cfg%Delta
  q      = model_cfg%q
  m0     = model_cfg%m0

  open(newunit=unit, file=filename, status='old', action='read', iostat=ierr)
  if (ierr == 0) then
    read(unit, nml=physics, iostat=ierr)
    close(unit)
    if (ierr == 0) then
      model_cfg%D      = D
      model_cfg%lambda = lambda
      model_cfg%A      = A
      model_cfg%Delta  = Delta
      model_cfg%q      = q
      model_cfg%m0     = m0
    end if
  end if

  call init_model_config(model_cfg)
end subroutine load

end module model_config_mod

module model_mod
  use precision
  use model_config_mod
  use grid_config_mod
  implicit none

  ! TODO: set number of equations/fields for your model
  integer, parameter :: NEQ = 3

  integer, save :: diag_unit = -1, fields_unit = -1, derivs_unit = -1
  logical, save :: diag_open = .false., fields_open = .false., derivs_open = .false.

  private
  public :: F, init_cond
  public :: compute_diagnostics, write_output
  public :: open_output_files, close_output_files
  public :: NEQ

contains

! Compute diagnostics (e.g., mass, Ricci). Placeholder returns zeros.
subroutine compute_diagnostics(mass, ricci, h, dhdu, dhdv, dhduv, model_cfg)
  real(dp), intent(out)               :: mass, ricci
  real(dp), dimension(:), intent(in)  :: h, dhdu, dhdv, dhduv
  type(model_config), intent(in)      :: model_cfg
  ! TODO: implement diagnostics for your model
  mass = 0.0_dp
  ricci = 0.0_dp
end subroutine compute_diagnostics

! PDE right-hand side: dhduv = F(h, dhdu, dhdv)
subroutine F(dhduv, h, dhdu, dhdv, model_cfg)
  real(dp), dimension(NEQ), intent(out) :: dhduv
  real(dp), dimension(NEQ), intent(in)  :: h, dhdu, dhdv
  type(model_config), intent(in)        :: model_cfg
  ! TODO: implement your RHS
  dhduv = 0.0_dp
end subroutine F

! Initialize boundary conditions on u_min and v_min
subroutine init_cond(h_u0, h_v0, grid_cfg, model_cfg)
  real(dp), dimension(:,:), intent(inout) :: h_u0, h_v0
  type(grid_config),  intent(in)          :: grid_cfg
  type(model_config), intent(in)          :: model_cfg
  integer :: i
  ! TODO: implement meaningful initial data for your physics
  do i = 1, grid_cfg%Nv
    h_u0(i, 1) = grid_cfg%v_min + (i-1) * grid_cfg%dv
    if (NEQ >= 2) h_u0(i, 2) = 0.0_dp
    if (NEQ >= 3) h_u0(i, 3) = 0.0_dp
  end do
  do i = 1, grid_cfg%Nu
    h_v0(i, 1) = h_u0(1, 1)
    if (NEQ >= 2) h_v0(i, 2) = 0.0_dp
    if (NEQ >= 3) h_v0(i, 3) = 0.0_dp
  end do
end subroutine init_cond

! Open output files and write headers
subroutine open_output_files(out_dir)
  character(len=*), intent(in) :: out_dir
  if (diag_open .or. fields_open .or. derivs_open) call close_output_files()
  open(newunit=fields_unit, file=trim(out_dir)//'/fields.dat', status='replace')
  open(newunit=diag_unit,   file=trim(out_dir)//'/diagnostics.dat', status='replace')
  open(newunit=derivs_unit, file=trim(out_dir)//'/derivatives.dat', status='replace')
  write(fields_unit, '(a)') '# Columns: u, r, phi, sigma'
  write(diag_unit,   '(a)') '# Columns: u, mass, Ricci'
  write(derivs_unit, '(a)') '# Columns: u, drdu, drdv'
  fields_open = .true.
  diag_open   = .true.
  derivs_open = .true.
end subroutine open_output_files

! Close output files
subroutine close_output_files()
  if (fields_open .and. fields_unit > 0) close(fields_unit)
  if (diag_open   .and. diag_unit   > 0) close(diag_unit)
  if (derivs_open .and. derivs_unit > 0) close(derivs_unit)
  fields_unit = -1; fields_open = .false.
  diag_unit   = -1; diag_open   = .false.
  derivs_unit = -1; derivs_open = .false.
end subroutine close_output_files

! Write outputs at midpoint P
subroutine write_output(u_val, v_val, h_N, h_S, h_E, h_W, du, dv, grid_cfg, model_cfg)
  real(dp), intent(in)               :: u_val, v_val, du, dv
  real(dp), dimension(:), intent(in) :: h_N, h_S, h_E, h_W
  type(grid_config), intent(in)      :: grid_cfg
  type(model_config), intent(in)     :: model_cfg

  real(dp), dimension(NEQ) :: h_P, dhdu_P, dhdv_P, dhduv_P
  real(dp) :: u_P, v_P, mass, ricci
  real(dp), save :: last_v_marked_val = -1.0e99_dp

  if (.not. (fields_open .and. diag_open .and. derivs_open)) error stop 'write_output: output files not open'

  u_P = u_val + 0.5_dp * du
  v_P = v_val + 0.5_dp * dv

  if (abs(v_val - last_v_marked_val) > max(1.0e-12_dp, 0.5_dp * grid_cfg%output_dv)) then
    write(fields_unit, '(a)')
    write(fields_unit, '(a,f10.6)') '# v = ', v_P
    write(diag_unit,   '(a)')
    write(diag_unit,   '(a,f10.6)') '# v = ', v_P
    write(derivs_unit, '(a)')
    write(derivs_unit, '(a,f10.6)') '# v = ', v_P
    last_v_marked_val = v_val
  end if

  h_P    = 0.25_dp * (h_N + h_S + h_E - h_W)
  dhdu_P = (h_W - h_S + h_N - h_E) * 0.5_dp / du
  dhdv_P = (h_E - h_S + h_N - h_W) * 0.5_dp / dv

  call F(dhduv_P, h_P, dhdu_P, dhdv_P, model_cfg)
  call compute_diagnostics(mass, ricci, h_P, dhdu_P, dhdv_P, dhduv_P, model_cfg)

  write(fields_unit, '(7e16.8)') u_P, h_P(1), merge(h_P(2), 0.0_dp, NEQ>=2), merge(h_P(3), 0.0_dp, NEQ>=3)
  write(diag_unit,   '(7e16.8)') u_P, mass, ricci
  write(derivs_unit, '(7e16.8)') u_P, dhdu_P(1), dhdv_P(1)
end subroutine write_output

end module model_mod
