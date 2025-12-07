! PDE time stepper
!
! Solves: \partial_{uv} h_j(u,v) = F_j(h_k, \partial_u h_k, \partial_v h_k)
!
! Input:
!   h_S, h_E, h_W : values at S=(u,v), E=(u,v+Dv), W=(u+Du,v)
!   Du, Dv        : step sizes in u and v directions
!   F             : RHS subroutine with signature F(dhduv, h, dhdu, dhdv)
! Output:
!   h_N           : solution at N=(u+Du, v+Dv)
!
subroutine evolve(h_N, h_S, h_E, h_W, Du, Dv)
  use functions
  implicit none

  double precision, dimension(NN), intent(out) :: h_N
  double precision, dimension(NN), intent(in)  :: h_S, h_E, h_W
  double precision, intent(in)                 :: Du, Dv

  integer :: j
  double precision, dimension(NN) :: h_P, dhdu_P, dhdv_P, dhduv_P, dhduv_P_new, h_N_old

  ! First-order approximation for h_N, then compute derivatives at P = (u+Du/2, v+Dv/2)
  h_N = h_W + h_E - h_S
  h_P = 0.5d0*(h_S + h_N)
  !h_P = 0.5*(h_E + h_W)
  !h_P = 0.25*(h_S + h_N + h_E + h_W)
  dhdu_P = (h_W - h_S + h_N - h_E)*0.5d0/Du
  dhdv_P = (h_E - h_S + h_N - h_W)*0.5d0/Dv

  ! Evaluate RHS at P
  call F(dhduv_P, h_P, dhdu_P, dhdv_P)

  ! Update h_N using RHS
  h_N = h_W + h_E - h_S + dhduv_P*Du*Dv

  ! Single Picard iteration for refinement (TODO: make this configurable)
  do j = 1, 1
    h_N_old = h_N

     !  h_P = 0.25*(h_S + h_N + h_E + h_W)
    dhdu_P = (h_W - h_S + h_N - h_E)*0.5d0/Du
    dhdv_P = (h_E - h_S + h_N - h_W)*0.5d0/Dv

    ! Re-evaluate RHS at P
    call F(dhduv_P_new, h_P, dhdu_P, dhdv_P)

    ! Final update to h_N
    h_N = h_W + h_E - h_S + Du*Dv*dhduv_P_new
  end do


end subroutine evolve
