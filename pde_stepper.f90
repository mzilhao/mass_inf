module pde_stepper
  use precision
  implicit none
  private
  public :: pde_step

  ! Abstract interface for RHS function
  ! Allows any subroutine with this signature to be passed
  abstract interface
    subroutine rhs_interface(dhduv, h, dhdu, dhdv, neq)
      import :: dp
      integer, intent(in) :: neq
      real(dp), dimension(neq), intent(out) :: dhduv
      real(dp), dimension(neq), intent(in)  :: h, dhdu, dhdv
    end subroutine rhs_interface
  end interface

contains

!> Predictor-corrector PDE time stepper
!!
!! Solves: ∂_{uv} h_j(u,v) = F_j(h_k, ∂_u h_k, ∂_v h_k)
!!
!! Uses a predictor-corrector scheme with optional Picard iterations.
!!
!! @param[out] h_N         Solution at point (u+du, v+dv) [North]
!! @param[in]  h_S         Values at point (u, v) [South]
!! @param[in]  h_E         Values at point (u, v+dv) [East]
!! @param[in]  h_W         Values at point (u+du, v) [West]
!! @param[in]  du          Step size in u direction
!! @param[in]  dv          Step size in v direction
!! @param[in]  rhs_func    Right-hand side F(dhduv, h, dhdu, dhdv, neq)
!! @param[in]  neq         Number of equations (system size)
!! @param[in]  n_picard    Number of Picard iterations for refinement (optional, default=1)
subroutine pde_step(h_N, h_S, h_E, h_W, du, dv, rhs_func, neq, n_picard)
  real(dp), dimension(:), intent(out) :: h_N
  real(dp), dimension(:), intent(in)  :: h_S, h_E, h_W
  real(dp), intent(in)                :: du, dv
  procedure(rhs_interface)                    :: rhs_func
  integer, intent(in)                         :: neq
  integer, intent(in), optional               :: n_picard

  integer :: j, n_iter
  real(dp) :: max_err
  real(dp), dimension(neq) :: h_P, dhdu_P, dhdv_P, dhduv_P, dhduv_P_new

  ! Validate input array sizes
  if (size(h_N) /= neq .or. size(h_S) /= neq .or. &
      size(h_E) /= neq .or. size(h_W) /= neq) then
    error stop "evolve: input array size mismatch with neq"
  end if

  ! Set number of Picard iterations
  n_iter = 1
  if (present(n_picard)) then
    if (n_picard < 1) error stop "evolve: n_picard must be >= 1"
    n_iter = n_picard
  end if

  ! PREDICTOR: Zero-order approximation
  h_N = h_W + h_E - h_S

  ! Evaluate at midpoint P = (u+du/2, v+dv/2)
  h_P = 0.25_dp * (h_N + h_S + h_E - h_W)

  ! Compute derivatives at P using centered differences
  dhdu_P = (h_W - h_S + h_N - h_E) * 0.5_dp / du
  dhdv_P = (h_E - h_S + h_N - h_W) * 0.5_dp / dv

  ! Evaluate RHS at P
  call rhs_func(dhduv_P, h_P, dhdu_P, dhdv_P, neq)

  ! CORRECTOR: Update with RHS contribution
  h_N = h_W + h_E - h_S + dhduv_P * du * dv

  ! PICARD ITERATIONS: Refine solution
  picard_loop: do j = 1, n_iter
    ! Recompute h_P and derivatives with refined h_N
    h_P = 0.25_dp * (h_N + h_S + h_E - h_W)
    dhdu_P = (h_W - h_S + h_N - h_E) * 0.5_dp / du
    dhdv_P = (h_E - h_S + h_N - h_W) * 0.5_dp / dv

    ! Re-evaluate RHS at P
    call rhs_func(dhduv_P_new, h_P, dhdu_P, dhdv_P, neq)

    ! Check convergence
    max_err = maxval(abs(dhduv_P_new - dhduv_P))
    ! print *, "Picard iteration ", j, ": max_err = ", max_err
    if (max_err < 1.0e-8_dp) exit picard_loop
    dhduv_P = dhduv_P_new

    ! Update solution
    h_N = h_W + h_E - h_S + dhduv_P_new * du * dv
  end do picard_loop

end subroutine pde_step

end module pde_stepper
