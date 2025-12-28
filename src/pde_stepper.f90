module pde_stepper_mod
  use precision
  implicit none
  private
  public :: pde_step

  ! Abstract interface for RHS function
  ! Allows any subroutine with this signature to be passed
  abstract interface
    subroutine rhs_interface(dhduv, h, dhdu, dhdv)
      import :: dp
      real(dp), dimension(:), intent(out) :: dhduv
      real(dp), dimension(:), intent(in)  :: h, dhdu, dhdv
    end subroutine rhs_interface
  end interface

contains

!> Predictor-corrector PDE time stepper
!!
!! Solves: ∂_{uv} h_j(u,v) = F_j(h_k, ∂_u h_k, ∂_v h_k)
!!
!! Uses a predictor-corrector scheme with Picard iterations.
!!
!! @param[out] h_N         Solution at point (u+du, v+dv) [North]
!! @param[in]  h_S         Values at point (u, v) [South]
!! @param[in]  h_E         Values at point (u, v+dv) [East]
!! @param[in]  h_W         Values at point (u+du, v) [West]
!! @param[in]  du          Step size in u direction
!! @param[in]  dv          Step size in v direction
!! @param[in]  rhs_func    Right-hand side F(dhduv, h, dhdu, dhdv)
!! @param[in]  n_iter      Number of Picard iterations for refinement
!! @param[in]  tol         Convergence tolerance (optional, default=1.0e-8)
subroutine pde_step(h_N, h_S, h_E, h_W, du, dv, rhs_func, n_iter, tol)
  real(dp), dimension(:), intent(out) :: h_N
  real(dp), dimension(:), intent(in)  :: h_S, h_E, h_W
  real(dp), intent(in)                :: du, dv
  procedure(rhs_interface)            :: rhs_func
  integer, intent(in)                 :: n_iter
  real(dp), intent(in), optional      :: tol

  integer  :: j, neq
  real(dp) :: max_err, tolerance
  real(dp), dimension(size(h_S)) :: h_P, dhdu_P, dhdv_P, dhduv_P, h_N_old

  ! Set convergence tolerance for the Picard iterations
  tolerance = 1.0e-8_dp
  if (present(tol)) tolerance = tol

  ! Validate input array sizes
  neq = size(h_S)
  if (size(h_N) /= neq .or. size(h_E) /= neq .or. size(h_W) /= neq) then
    error stop "pde_step: input array size mismatch"
  end if

  ! PREDICTOR: Zero-order approximation
  h_N = h_W + h_E - h_S

  ! Evaluate at midpoint P = (u+du/2, v+dv/2)
  h_P = 0.25_dp * (h_N + h_S + h_E - h_W)

  ! Compute derivatives at P using centered differences
  dhdu_P = (h_W - h_S + h_N - h_E) * 0.5_dp / du
  dhdv_P = (h_E - h_S + h_N - h_W) * 0.5_dp / dv

  ! Evaluate RHS at P
  call rhs_func(dhduv_P, h_P, dhdu_P, dhdv_P)

  ! CORRECTOR: Update with RHS contribution
  h_N = h_W + h_E - h_S + dhduv_P * du * dv

  ! PICARD ITERATIONS: Refine solution
  picard_loop: do j = 1, n_iter
    ! Store previous h_N for convergence check later
    h_N_old = h_N

    ! Recompute h_P and derivatives with refined h_N
    h_P = 0.25_dp * (h_N + h_S + h_E - h_W)
    dhdu_P = (h_W - h_S + h_N - h_E) * 0.5_dp / du
    dhdv_P = (h_E - h_S + h_N - h_W) * 0.5_dp / dv

    ! Re-evaluate RHS at P
    call rhs_func(dhduv_P, h_P, dhdu_P, dhdv_P)

    ! Update solution
    h_N = h_W + h_E - h_S + dhduv_P * du * dv

    ! Check convergence
    max_err = maxval(abs(h_N_old - h_N))
    ! print *, "Picard iteration ", j, ": max_err = ", max_err
    if (max_err < tolerance) exit picard_loop

  end do picard_loop

end subroutine pde_step

end module pde_stepper_mod