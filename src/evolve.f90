! Wrapper module to interface with pde_stepper

module evolve_wrapper
  use precision
  use pde_stepper, only: pde_step
  use model_config_mod
  use model_mod, only: F
  implicit none
  private
  public :: step

  ! Module-level config storage for rhs_wrapper to access
  ! (Fortran doesn't support closures with explicit context passing)
  ! This is updated each call to step() to match the input model_cfg
  type(model_config), save :: cfg_module

contains

!> Solve PDE step using the modular pde_stepper
!!
!! Note: The physics_config is stored in module scope so that the
!! RHS wrapper (rhs_wrapper) can access it. This is a limitation of
!! Fortran's procedure pointer/abstract interface mechanism, which
!! doesn't support passing extra context through the call stack.
!!
!! @param[out] h_N       Solution at point (u+du, v+dv)
!! @param[in]  h_S, h_E, h_W  Boundary values
!! @param[in]  du, dv    Step sizes
!! @param[in]  cfg       Physics configuration (stored in module scope)
!! @param[in]  n_picard  Number of Picard iterations
subroutine step(h_N, h_S, h_E, h_W, du, dv, model_cfg, n_picard)
  real(dp), dimension(:), intent(out) :: h_N
  real(dp), dimension(:), intent(in)  :: h_S, h_E, h_W
  real(dp), intent(in)                :: du, dv
  type(model_config), intent(in)      :: model_cfg
  integer, intent(in)                 :: n_picard

  ! Store cfg in module scope for rhs_wrapper to access
  cfg_module = model_cfg

  call pde_step(h_N, h_S, h_E, h_W, du, dv, rhs_wrapper, n_picard)
end subroutine step

!> RHS wrapper that accesses physics config from module scope
subroutine rhs_wrapper(dhduv, h, dhdu, dhdv)
  real(dp), dimension(:), intent(out) :: dhduv
  real(dp), dimension(:), intent(in)  :: h, dhdu, dhdv

  call F(dhduv, h, dhdu, dhdv, cfg_module)
end subroutine rhs_wrapper

end module evolve_wrapper
