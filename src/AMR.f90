module amr_mod
  use precision
  use polint_mod
  use utils_mod, only: relative_difference
  implicit none
  private
  public :: refine_u_grid
contains

subroutine refine_u_grid(u, h_v0, h_W, h_S, j, next_idx, u_max, reldiff_max, plus, minus)
  real(dp),          intent(inout) :: u(:)
  real(dp),          intent(inout) :: h_v0(:,:)
  real(dp),          intent(inout) :: h_W(:)
  real(dp),          intent(in)    :: h_S(:)
  integer,           intent(inout) :: j
  integer,           intent(inout) :: next_idx
  real(dp),          intent(in)    :: u_max
  real(dp),          intent(in)    :: reldiff_max
  integer,           intent(inout) :: plus(:), minus(:)

  real(dp) :: reldiff_r
  integer  :: jm1, jm2, jm3, jp1
  integer  :: k, neq
  real(dp), dimension(5) :: interp_x, interp_y

  neq = size(h_W)

  reldiff_r = relative_difference(h_W(1), h_S(1))

  do while (reldiff_r > reldiff_max .and. j >= 4)
    jm1 = minus(j)
    jm2 = minus(jm1)
    jm3 = minus(jm2)
    jp1 = plus(j)

    do k = 1, neq
      if (u(jp1) < u_max) then
        interp_x = [ u(1), u(jm2), u(jm1), u(j), u(jp1) ]
        interp_y = [ h_v0(1,k), h_v0(jm2,k), h_v0(jm1,k), h_v0(j,k), &
                     h_v0(jp1,k) ]
      else
        interp_x = [ u(1), u(jm3), u(jm2), u(jm1), u(j) ]
        interp_y = [ h_v0(1,k), h_v0(jm3,k), h_v0(jm2,k), h_v0(jm1,k), &
                     h_v0(j,k) ]
      end if
      h_W(k) = polint((u(j) + u(jm1)) * 0.5_dp, interp_x, interp_y)
    end do

    u(next_idx)       = (u(j) + u(jm1)) * 0.5_dp
    h_v0(next_idx, :) = h_W(:)

    minus(next_idx) = jm1
    plus(next_idx)  = j
    plus(jm1)       = next_idx
    minus(j)        = next_idx
    j               = next_idx
    next_idx        = next_idx + 1

    reldiff_r = relative_difference(h_W(1), h_S(1))
  end do
end subroutine refine_u_grid

end module amr_mod
