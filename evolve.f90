! Wrapper module to interface with pde_stepper

module evolve_wrapper
  use pde_stepper, only: pde_step
  use physics_config_mod
  implicit none
  private
  public :: step

  ! Module-level config storage for rhs_wrapper to access
  ! (Fortran doesn't support closures with explicit context passing)
  ! This is updated each call to step() to match the input cfg
  type(physics_config), save :: cfg_module

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
!! @param[in]  n_picard  Optional Picard iterations (default=1)
subroutine step(h_N, h_S, h_E, h_W, du, dv, cfg, n_picard)
  double precision, dimension(:), intent(out) :: h_N
  double precision, dimension(:), intent(in)  :: h_S, h_E, h_W
  double precision, intent(in)                :: du, dv
  type(physics_config), intent(in)            :: cfg
  integer, intent(in), optional               :: n_picard

  ! Store cfg in module scope for rhs_wrapper to access
  cfg_module = cfg

  call pde_step(h_N, h_S, h_E, h_W, du, dv, rhs_wrapper, NEQ, n_picard)
end subroutine step

!> RHS wrapper that accesses physics config from module scope
subroutine rhs_wrapper(dhduv, h, dhdu, dhdv, neq)
  use functions, only: F
  integer, intent(in) :: neq
  double precision, dimension(neq), intent(out) :: dhduv
  double precision, dimension(neq), intent(in)  :: h, dhdu, dhdv

  call F(dhduv, h, dhdu, dhdv, neq, cfg_module)
end subroutine rhs_wrapper

end module evolve_wrapper
