! Flat model implementation
! Required public API (to be implemented by any model):
! - module model_config_mod
!     type(model_config)
!     subroutine load(model_cfg, filename)
! - module model_mod
!     integer, parameter :: NEQ
!     subroutine F(dhduv, h, dhdu, dhdv, model_cfg)
!     subroutine init_cond(h_u0, h_v0, grid_cfg, model_cfg)
!     subroutine open_output_files(out_dir)
!     subroutine write_output(u_val, v_val, h_N, h_S, h_E, h_W, du, dv, grid_cfg, model_cfg)
!     subroutine close_output_files()
!
! This flat-space implementation compiles and runs, returning trivial results.
! Use it as a reference for the minimal set of functions a model must provide.
module model_config_mod
  use precision
  implicit none
  private
  public :: model_config, load

  !> Model configuration type - encapsulates all model parameters
  type :: model_config
    real(dp) :: A = 0.0_dp               ! Flat model parameter
    real(dp) :: q = 0.95_dp              ! Electric charge
    real(dp) :: m0 = 1.0_dp              ! Initial mass parameter
  end type model_config

contains

! Read &model namelist (optional) and compute derived constants
subroutine load(model_cfg, filename)
  type(model_config), intent(out) :: model_cfg
  character(len=*), intent(in)    :: filename

  ! Local variables for namelist reading
  real(dp) :: A, q, m0
  namelist /model/ A, q, m0

  integer :: unit, ierr

  ! Initialize with type defaults
  model_cfg = model_config()
  A         = model_cfg%A
  q         = model_cfg%q
  m0        = model_cfg%m0

  ! Read namelist
  open(newunit=unit, file=filename, status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
    write(*, '(a,a,a)') 'Error: could not open parameter file "', trim(filename), '"'
    call exit(1)
  end if

  read(unit, nml=model, iostat=ierr)
  if (ierr /= 0) then
    write(*, '(a,a,a)') 'Error: could not read &model namelist from "', trim(filename), '"'
    close(unit)
    call exit(1)
  end if
  close(unit)

  ! Update cfg with (possibly modified) namelist values
  model_cfg%A      = A
  model_cfg%q      = q
  model_cfg%m0     = m0

end subroutine load

end module model_config_mod

!====================================================================================
!====================================================================================
module model_mod
  use precision
  use model_config_mod
  use grid_config_mod
  implicit none

  integer, parameter  :: NEQ = 3      ! Number of equations/fields: r, phi, sigma

  integer, save :: fields_unit = -1, derivs_unit = -1
  logical, save :: fields_open = .false., derivs_open = .false.

  private
  public :: NEQ
  public :: F, init_cond
  public :: open_output_files, write_output, close_output_files

contains

!====================================================================================
!> Right-hand side of the PDE system
!! Inputs: h (field values), dhdu, dhdv (derivatives)
!! Output: dhduv (mixed derivatives)
subroutine F(dhduv, h, dhdu, dhdv, model_cfg)
  implicit none

  real(dp), dimension(NEQ), intent(out) :: dhduv
  real(dp), dimension(NEQ), intent(in)  :: h, dhdu, dhdv
  type(model_config), intent(in)        :: model_cfg

  ! h(1) = r, h(2) = phi, h(3) = sigma
  dhduv(1) = 0.0_dp

  dhduv(2) = 1.0_dp/ h(1) * (dhdv(1)*dhdu(2) + dhdu(1)*dhdv(2))

  dhduv(3) = 0.0_dp

end subroutine F

!====================================================================================
!> Initialize boundary conditions at u=u_min and v=v_min
!! Returns: h_u0 (IC along u_min), h_v0 (IC along v_min)
subroutine init_cond(h_u0, h_v0, grid_cfg, model_cfg)
  implicit none

  real(dp), dimension(:,:), intent(inout) :: h_u0, h_v0
  type(grid_config),  intent(in)          :: grid_cfg
  type(model_config), intent(in)          :: model_cfg

  integer  :: i
  real(dp) :: u, r00, ru0

  r00 = grid_cfg%v_min
  ru0 = 0.5_dp / r00 * (model_cfg%m0 - 0.5_dp * model_cfg%q**2 / r00 - 1.0_dp)

  ! r = v along u_min; phi = 0; sigma = 0
  do i = 1, grid_cfg%Nv
    h_u0(i, 1) = grid_cfg%v_min + (i-1) * grid_cfg%dv
    h_u0(i, 2) = 0.0_dp
    h_u0(i, 3) = 0.0_dp
  end do

  ! Boundary conditions at v = v_min
  do i = 1, grid_cfg%Nu
    u = grid_cfg%u_min + (i-1) * grid_cfg%du
    h_v0(i, 1) = r00 + u * ru0  ! r(u,v_min)
    h_v0(i, 2) = 0.0_dp
    h_v0(i, 3) = 0.0_dp
  end do
end subroutine init_cond

!====================================================================================
!> Open all output files
subroutine open_output_files(out_dir)
  character(len=*), intent(in)      :: out_dir

  if (fields_open .or. derivs_open) call close_output_files()

  open(newunit=fields_unit, file=trim(out_dir)//'/fields.dat', status='replace')
  open(newunit=derivs_unit, file=trim(out_dir)//'/derivatives.dat', status='replace')

  write(fields_unit, '(a)') '# Columns: u, r, phi, sigma'
  write(derivs_unit, '(a)') '# Columns: u, dphidu, dphidv'

  fields_open = .true.
  derivs_open = .true.
end subroutine open_output_files

!====================================================================================
!> Close all output files if open
subroutine close_output_files()
  if (fields_open .and. fields_unit > 0) close(fields_unit)
  if (derivs_open .and. derivs_unit > 0) close(derivs_unit)
  fields_unit = -1; fields_open = .false.
  derivs_unit = -1; derivs_open = .false.
end subroutine close_output_files

!====================================================================================
!> Compute diagnostics and write them if output cadence is met
!!
!! Writes ASCII output in columns.
!! v-slices are marked with comment lines: # v = X.XXXXX
subroutine write_output(u_val, v_val, h_N, h_S, h_E, h_W, du, dv, grid_cfg, model_cfg)
  real(dp), intent(in)               :: u_val, v_val, du, dv
  real(dp), dimension(:), intent(in) :: h_N, h_S, h_E, h_W
  type(grid_config), intent(in)      :: grid_cfg
  type(model_config), intent(in)     :: model_cfg

  real(dp), dimension(NEQ) :: h_P, dhdu_P, dhdv_P, dhduv_P
  logical  :: u_ok, v_ok
  real(dp) :: temp, u_P, v_P

  real(dp), save :: last_v_marked_val = -1.0e99_dp

  if (.not. (fields_open .and. derivs_open)) error stop 'write_output: output files not open'

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
    write(derivs_unit, '(a)')
    write(derivs_unit, '(a,f10.6)') '# v = ', v_P
    last_v_marked_val = v_val
  end if

  h_P    = 0.25_dp * (h_N + h_S + h_E - h_W)
  dhdu_P = (h_W - h_S + h_N - h_E) * 0.5_dp / du
  dhdv_P = (h_E - h_S + h_N - h_W) * 0.5_dp / dv

  write(fields_unit, '(7e16.8)') u_P, h_P(1), h_P(2), h_P(3)
  write(derivs_unit, '(7e16.8)') u_P, dhdu_P(2), dhdv_P(2)
end subroutine write_output

end module model_mod
