! Template model implementation
! Copy this file to create a new model and replace TODOs.
! Implement the stable API expected by callers:
! - module model_config_mod: defines type(model_config), load
! - module model_mod: provides NEQ and routines:
!     F, init_cond,
!     open_output_files, write_output, close_output_files

module model_config_mod
  use precision
  implicit none
  private
  public :: model_config, load

  ! TODO: customize configuration fields for your physics
  type :: model_config
    real(dp) :: q  = 0.1_dp              ! Electric charge
    real(dp) :: m0 = 1.0_dp              ! Initial mass parameter
  end type model_config

contains

!====================================================================================
! Read &physics namelist and initialize model_cfg
subroutine load(model_cfg, filename)
  type(model_config), intent(out) :: model_cfg
  character(len=*), intent(in)    :: filename

  ! Local variables for namelist reading
  ! TODO
  real(dp) :: q, m0
  namelist /physics/ q, m0

  integer :: unit, ierr

  ! Initialize with type defaults
  ! TODO
  model_cfg = model_config()
  q         = model_cfg%q
  m0        = model_cfg%m0

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
  ! TODO
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

  ! TODO: set number of equations/fields for your model
  integer, parameter :: NEQ = 1

  private
  public :: NEQ
  public :: F, init_cond
  public :: open_output_files, write_output, close_output_files

contains

!====================================================================================
! PDE right-hand side: dhduv = F(h, dhdu, dhdv)
subroutine F(dhduv, h, dhdu, dhdv, model_cfg)
  implicit none

  real(dp), dimension(NEQ), intent(out) :: dhduv
  real(dp), dimension(NEQ), intent(in)  :: h, dhdu, dhdv
  type(model_config), intent(in)        :: model_cfg

  ! TODO: implement your RHS
  dhduv = 0.0_dp
end subroutine F

!====================================================================================
! Initialize boundary conditions on u_min and v_min
subroutine init_cond(h_u0, h_v0, grid_cfg, model_cfg)
  real(dp), dimension(:,:), intent(inout) :: h_u0, h_v0
  type(grid_config),  intent(in)          :: grid_cfg
  type(model_config), intent(in)          :: model_cfg

  ! TODO: implement meaningful initial data for your physics
  h_u0 = 0.0_dp
  h_v0 = 0.0_dp
end subroutine init_cond

!====================================================================================
! Open output files and write headers
subroutine open_output_files(out_dir)
  character(len=*), intent(in) :: out_dir
  ! TODO
end subroutine open_output_files

!====================================================================================
! Close output files
subroutine close_output_files()
  ! TODO
end subroutine close_output_files

!====================================================================================
! Write outputs at midpoint P
subroutine write_output(u_val, v_val, h_N, h_S, h_E, h_W, du, dv, grid_cfg, model_cfg)
  real(dp), intent(in)               :: u_val, v_val, du, dv
  real(dp), dimension(:), intent(in) :: h_N, h_S, h_E, h_W
  type(grid_config), intent(in)      :: grid_cfg
  type(model_config), intent(in)     :: model_cfg

  real(dp), dimension(NEQ) :: h_P, dhdu_P, dhdv_P, dhduv_P
  real(dp) :: u_P, v_P
  real(dp), save :: last_v_marked_val = -1.0e99_dp

  h_P    = 0.25_dp * (h_N + h_S + h_E - h_W)
  dhdu_P = (h_W - h_S + h_N - h_E) * 0.5_dp / du
  dhdv_P = (h_E - h_S + h_N - h_W) * 0.5_dp / dv

  ! TODO
end subroutine write_output

end module model_mod
