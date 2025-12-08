module pde_stepper
  implicit none
  private
  public :: pde_step

  ! Abstract interface for RHS function
  ! Allows any subroutine with this signature to be passed
  abstract interface
    subroutine rhs_interface(dhduv, h, dhdu, dhdv, neq)
      integer, intent(in) :: neq
      double precision, dimension(neq), intent(out) :: dhduv
      double precision, dimension(neq), intent(in)  :: h, dhdu, dhdv
    end subroutine rhs_interface
  end interface

contains

  !> Predictor-corrector PDE time stepper
  !!
  !! Solves: ∂_{uv} h_j(u,v) = F_j(h_k, ∂_u h_k, ∂_v h_k)
  !!
  !! Uses a predictor-corrector scheme with optional Picard iterations.
  !!
  !! @param[out] h_N         Solution at point (u+Du, v+Dv) [North]
  !! @param[in]  h_S         Values at point (u, v) [South]
  !! @param[in]  h_E         Values at point (u, v+Dv) [East]
  !! @param[in]  h_W         Values at point (u+Du, v) [West]
  !! @param[in]  Du          Step size in u direction
  !! @param[in]  Dv          Step size in v direction
  !! @param[in]  rhs_func    Right-hand side F(dhduv, h, dhdu, dhdv, neq)
  !! @param[in]  neq         Number of equations (system size)
  !! @param[in]  n_picard    Number of Picard iterations for refinement (optional, default=1)
  subroutine pde_step(h_N, h_S, h_E, h_W, Du, Dv, rhs_func, neq, n_picard)
    double precision, dimension(:), intent(out) :: h_N
    double precision, dimension(:), intent(in)  :: h_S, h_E, h_W
    double precision, intent(in)                :: Du, Dv
    procedure(rhs_interface)                    :: rhs_func
    integer, intent(in)                         :: neq
    integer, intent(in), optional               :: n_picard

    integer :: j, n_iter
    double precision, dimension(neq) :: h_P, dhdu_P, dhdv_P, dhduv_P, dhduv_P_new

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

    ! PREDICTOR: First-order approximation
    ! h_N ≈ h_W + h_E - h_S (from explicit scheme)
    h_N = h_W + h_E - h_S

    ! Evaluate at midpoint P = (u+Du/2, v+Dv/2)
    h_P = 0.5d0 * (h_S + h_N)

    ! Compute derivatives at P using centered differences
    dhdu_P = (h_W - h_S + h_N - h_E) * 0.5d0 / Du
    dhdv_P = (h_E - h_S + h_N - h_W) * 0.5d0 / Dv

    ! Evaluate RHS at P
    call rhs_func(dhduv_P, h_P, dhdu_P, dhdv_P, neq)

    ! CORRECTOR: Update with RHS contribution
    h_N = h_W + h_E - h_S + dhduv_P * Du * Dv

    ! PICARD ITERATIONS: Refine solution
    picard_loop: do j = 1, n_iter
      ! Recompute derivatives with refined h_N
      dhdu_P = (h_W - h_S + h_N - h_E) * 0.5d0 / Du
      dhdv_P = (h_E - h_S + h_N - h_W) * 0.5d0 / Dv

      ! Re-evaluate RHS at P
      call rhs_func(dhduv_P_new, h_P, dhdu_P, dhdv_P, neq)

      ! Update solution
      h_N = h_W + h_E - h_S + Du * Dv * dhduv_P_new
    end do picard_loop

  end subroutine pde_step

end module pde_stepper
