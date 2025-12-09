! Wrapper module to interface with pde_stepper

module evolve_wrapper
  use pde_stepper, only: pde_step
  use physics_config_mod
  implicit none
  private
  public :: step, set_cfg

  ! Module-level config (set by main program before first evolve call)
  type(physics_config) :: cfg_global

contains

  !> Set physics config for subsequent evolve calls
  subroutine set_cfg(cfg)
    type(physics_config), intent(in) :: cfg
    cfg_global = cfg
  end subroutine set_cfg

  !> Solve PDE step using the modular pde_stepper
  subroutine step(h_N, h_S, h_E, h_W, du, dv, cfg, n_picard)
    double precision, dimension(:), intent(out) :: h_N
    double precision, dimension(:), intent(in)  :: h_S, h_E, h_W
    double precision, intent(in)                :: du, dv
    type(physics_config), intent(in)            :: cfg
    integer, intent(in), optional               :: n_picard

    call pde_step(h_N, h_S, h_E, h_W, du, dv, rhs_wrapper, cfg%neq, n_picard)
  end subroutine step

  !> RHS wrapper that captures physics config
  subroutine rhs_wrapper(dhduv, h, dhdu, dhdv, neq)
    use functions, only: F
    integer, intent(in) :: neq
    double precision, dimension(neq), intent(out) :: dhduv
    double precision, dimension(neq), intent(in)  :: h, dhdu, dhdv

    ! Note: This accesses 'cfg_global' from the enclosing module scope
    ! In a full refactoring, you'd pass cfg through the call stack
    ! For now, we maintain a module-level cfg (set in main program)
    call F(dhduv, h, dhdu, dhdv, neq, cfg_global)
  end subroutine rhs_wrapper

end module evolve_wrapper

