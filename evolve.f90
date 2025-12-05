! PDE time stepper
!
! Solves: \partial_{uv} h_j(u,v) = F_j(h_k, \partial_u h_k, \partial_v h_k)
!
! Input:
!   hp1, hp2, hp3 : values at p1=(u,v), p2=(u,v+Dv), p3=(u+Du,v)
!   Du, Dv        : step sizes in u and v directions
!   F              : RHS subroutine with signature F(dhduv, h, dhdu, dhdv)
! Output:
!   hp4            : solution at p4=(u+Du, v+Dv)
!
subroutine evolve(hp4, hp1, hp2, hp3, Du, Dv)
  use functions
  implicit none

  double precision, dimension(NN), intent(out) :: hp4
  double precision, dimension(NN), intent(in)  :: hp1, hp2, hp3
  double precision, intent(in)                 :: Du, Dv

  integer :: j
  double precision, dimension(NN) :: hp5, dhdup5, dhdvp5, dhduvp5, dhduvp5new, hp4old

  ! First-order approximation for hp4, then compute derivatives at p5 = (u+Du/2, v+Dv/2)
  hp4 = hp3 + hp2 - hp1
  hp5 = 0.5d0*(hp1 + hp4)
  !hp5 = 0.5*(hp2 + hp3)
  !hp5 = 0.25*(hp1 + hp4 + hp2 + hp3)
  dhdup5 = (hp3 - hp1 + hp4 - hp2)*0.5d0/Du
  dhdvp5 = (hp2 - hp1 + hp4 - hp3)*0.5d0/Dv

  ! Evaluate RHS at p5
  call F(dhduvp5, hp5, dhdup5, dhdvp5)

  ! Update hp4 using RHS
  hp4 = hp3 + hp2 - hp1 + dhduvp5*Du*Dv

  ! Single Picard iteration for refinement (TODO: make this configurable)
  do j = 1, 1
    hp4old = hp4

     !  hp5 = 0.25*(hp1 + hp4 + hp2 + hp3)
    dhdup5 = (hp3 - hp1 + hp4 - hp2)*0.5d0/Du
    dhdvp5 = (hp2 - hp1 + hp4 - hp3)*0.5d0/Dv

    ! Re-evaluate RHS at p5
    call F(dhduvp5new, hp5, dhdup5, dhdvp5)

    ! Final update to hp4
    hp4 = hp3 + hp2 - hp1 + Du*Dv*dhduvp5new
  end do


end subroutine evolve
