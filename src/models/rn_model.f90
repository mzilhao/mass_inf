module model_config_mod
  use precision
  implicit none
  private
  public :: model_config, load

  !> Model configuration type - encapsulates all model parameters
  type :: model_config
    real(dp) :: A = 0.0_dp               ! Scalar field amplitude
    real(dp) :: Delta = 1.0_dp           ! Scalar field width
    real(dp) :: q  = 0.95_dp             ! Electric charge
    real(dp) :: m0 = 1.0_dp              ! Initial mass parameter
  end type model_config

contains

!====================================================================================
!> Read model configuration from a namelist file
subroutine load(model_cfg, filename)
  type(model_config), intent(out) :: model_cfg
  character(len=*), intent(in)    :: filename

  ! Local variables for namelist reading
  real(dp) :: A, Delta, q, m0
  namelist /model/ A, Delta, q, m0

  integer :: unit, ierr

  ! Initialize with type defaults
  model_cfg = model_config()
  A      = model_cfg%A
  Delta  = model_cfg%Delta
  q      = model_cfg%q
  m0     = model_cfg%m0

  ! Read namelist
  open(newunit=unit, file=filename, status='old', action='read', iostat=ierr)
  if (ierr /= 0) then
    write(*, '(a,a,a)') 'Error: could not open parameter file "', trim(filename), '"'
    call exit(1)
  end if

  read(unit, nml=model, iostat=ierr)
  if (ierr /= 0) then
    write(*, '(a,a,a)') 'Error: could not read &model namelist from "', trim(filename), '"'
    close(unit)
    call exit(1)
  end if
  close(unit)

  ! Update cfg with (possibly modified) namelist values
  model_cfg%A      = A
  model_cfg%Delta  = Delta
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

  integer, parameter  :: NEQ = 3      ! Number of equations/fields: r, phi, sigma

  real(dp), parameter :: Pi = 4.0_dp * atan(1.0_dp)

  integer, save :: diag_unit = -1, fields_unit = -1, derivs_unit = -1, constraints_unit = -1
  logical, save :: diag_open = .false., fields_open = .false., derivs_open = .false.
  logical, save :: constraints_open = .false.

  private
  public :: NEQ
  public :: F, init_cond
  public :: open_output_files, write_output, write_constraints, close_output_files

contains

!====================================================================================
!> Right-hand side of the PDE system
!! Inputs: h (field values), dhdu, dhdv (derivatives)
!! Output: dhduv (mixed derivatives)
subroutine F(dhduv, h, dhdu, dhdv, model_cfg)
  real(dp), dimension(NEQ), intent(out) :: dhduv
  real(dp), dimension(NEQ), intent(in)  :: h, dhdu, dhdv
  type(model_config), intent(in)        :: model_cfg

  real(dp) :: r, sigma, phi, q

  q = model_cfg%q

  r = h(1); sigma = h(2); phi = h(3)

  dhduv(1) = 0.5_dp * exp(2*sigma) / r * (q**2/r**2 - 1) &
           - dhdu(1) * dhdv(1) / r

  dhduv(2) = -dhdu(3) * dhdv(3)                               &
           - 0.5_dp * exp(2*sigma) / r**2 * (2*q**2/r**2 - 1) &
           + dhdu(1) * dhdv(1) / r**2

  dhduv(3) = -1.0_dp / r * (dhdv(1)*dhdu(3) + dhdu(1)*dhdv(3))

end subroutine F

!====================================================================================
!> Initialize boundary conditions at u=u_min and v=v_min
!! Returns: h_u0 (IC along u_min), h_v0 (IC along v_min)
subroutine init_cond(h_u0, h_v0, grid_cfg, model_cfg)
  real(dp), dimension(:,:), intent(inout) :: h_u0, h_v0
  type(grid_config),  intent(in)          :: grid_cfg
  type(model_config), intent(in)          :: model_cfg

  ! Local variables
  integer  :: i
  real(dp) :: u, v
  real(dp) :: r00, sigma_0, ru0, v0, v1
  real(dp) :: A, Delta, q, m0

  ! Extract config values
  A      = model_cfg%A
  Delta  = model_cfg%Delta
  m0     = model_cfg%m0
  q      = model_cfg%q
  v0     = grid_cfg%v_min

  r00     = v0
  sigma_0 = -0.5_dp * log(2.0_dp)
  ru0     = 0.25_dp * (2.0_dp / r00 * (m0 - 0.5_dp*q*q/r00) - 1.0_dp)

  v1 = v0 + Delta

  ! Boundary conditions at v = v_min
  do i = 1, grid_cfg%Nu
    u = grid_cfg%u_min + (i-1) * grid_cfg%du
    h_v0(i, 1) = r00 + u * ru0  ! r(u,v_min)
    h_v0(i, 2) = sigma_0        ! sigma(u,v_min)
    h_v0(i, 3) = 0.0_dp         ! phi(u,v_min)
  end do

  ! Boundary conditions at u = u_min
  do i = 1, grid_cfg%Nv
    v = v0 + (i-1) * grid_cfg%dv
    h_u0(i, 1) = v  ! r(u_min,v)

    if (v <= v0 + Delta) then
      ! sigma(r_umin,v)
      h_u0(i, 2) = sigma_0              &
                 + A*A/(256*Pi*Pi) * (  &
                      15*v0**2 - 24*Pi*Pi*v0**2 - 30*v0*v1  &
                    + 15*v1**2 + 24*Pi*Pi*v**2              &
                    - 16*Delta**2*cos(2*Pi*(v0-v)/Delta)    &
                    + Delta**2*cos(4*Pi*(v0-v)/Delta)       &
                    - 32*Pi*v0*v*sin(2*Pi*(v0-v)/Delta)     &
                    + 32*Pi*v1*v*sin(2*Pi*(v0-v)/Delta)     &
                    + 4*Pi*v0*v*sin(4*Pi*(v0-v)/Delta)      &
                    - 4*Pi*v1*v*sin(4*Pi*(v0-v)/Delta)      &
                    )

      ! phi(r_umin,v)
      h_u0(i, 3) = A/(4*Pi) * (2*Pi*(v-v0) - Delta*sin(2*Pi*(v - v0)/Delta))
    else
      h_u0(i, 2) = sigma_0 + 3*A*A*Delta*(v0+v1)/32
      h_u0(i, 3) = 0.5_dp * A * Delta
    end if
  end do

end subroutine init_cond

!====================================================================================
!> Open all output files
subroutine open_output_files(out_dir)
  character(len=*), intent(in) :: out_dir

  if (diag_open .or. fields_open .or. derivs_open .or. constraints_open) call close_output_files()

  open(newunit=fields_unit,      file=trim(out_dir)//'/fields.dat',      status='replace')
  open(newunit=diag_unit,        file=trim(out_dir)//'/diagnostics.dat', status='replace')
  open(newunit=derivs_unit,      file=trim(out_dir)//'/derivatives.dat', status='replace')
  open(newunit=constraints_unit, file=trim(out_dir)//'/constraints.dat', status='replace')

  write(fields_unit,      '(a)') '# Columns: u, r, sigma, phi'
  write(diag_unit,        '(a)') '# Columns: u, mass, Ricci'
  write(derivs_unit,      '(a)') '# Columns: u, drdu, drdv'
  write(constraints_unit, '(a)') '# Columns: u, Guu'

  fields_open      = .true.
  diag_open        = .true.
  derivs_open      = .true.
  constraints_open = .true.
end subroutine open_output_files

!====================================================================================
!> Close all output files if open
subroutine close_output_files()
  if (fields_open      .and. fields_unit      > 0) close(fields_unit)
  if (diag_open        .and. diag_unit        > 0) close(diag_unit)
  if (derivs_open      .and. derivs_unit      > 0) close(derivs_unit)
  if (constraints_open .and. constraints_unit > 0) close(constraints_unit)
  fields_unit      = -1; fields_open      = .false.
  diag_unit        = -1; diag_open        = .false.
  derivs_unit      = -1; derivs_open      = .false.
  constraints_unit = -1; constraints_open = .false.
end subroutine close_output_files

!====================================================================================
!> Compute diagnostics and write them if output cadence is met
!!
!! Writes ASCII output in columns.
!! v-slices are marked with comment lines: # v = X.XXXXX
subroutine write_output(u_val, v_val, h_N, h_S, h_E, h_W, du, dv, grid_cfg, model_cfg)
  real(dp), intent(in)               :: u_val, v_val, du, dv
  real(dp), dimension(:), intent(in) :: h_N, h_S, h_E, h_W
  type(grid_config), intent(in)      :: grid_cfg
  type(model_config), intent(in)     :: model_cfg

  real(dp), dimension(NEQ) :: h_P, dhdu_P, dhdv_P, dhduv_P
  logical  :: u_ok, v_ok
  real(dp) :: temp, u_P, v_P, mass, ricci

  real(dp), save :: last_v_marked_val = -1.0e99_dp

  if (.not. (fields_open .and. diag_open .and. derivs_open)) error stop 'write_output: output files not open'

  ! Check output condition first, before doing any work
  ! Output condition: write when current (u,v) aligns with sampling spacings.
  ! We check if (u - u_min)/output_du and (v - v_min)/output_dv are near integers.
  ! This is robust to floating-point drift and local AMR changes.
  if (grid_cfg%output_du > 0.0_dp) then
    temp = abs((u_val - grid_cfg%u_min) / grid_cfg%output_du)
    u_ok = abs(temp - nint(temp)) < 1.0e-7_dp
  else
    u_ok = .true.
  end if

  if (grid_cfg%output_dv > 0.0_dp) then
    temp = abs((v_val - grid_cfg%v_min) / grid_cfg%output_dv)
    v_ok = abs(temp - nint(temp)) < 1.0e-6_dp
  else
    v_ok = .true.
  end if

  if (.not. (u_ok .and. v_ok)) return

  ! The diagnostics are computed at the midpoint P = (u+du/2, v+dv/2)
  u_P = u_val + 0.5_dp * du
  v_P = v_val + 0.5_dp * dv

  ! Write a new v-slice block if v_val is more than one half output_dv away from last marked
  if (abs(v_val - last_v_marked_val) > 0.5_dp * grid_cfg%output_dv) then
    write(fields_unit, '(a)')
    write(fields_unit, '(a,f10.6)') '# v = ', v_P
    write(diag_unit,   '(a)')
    write(diag_unit,   '(a,f10.6)') '# v = ', v_P
    write(derivs_unit, '(a)')
    write(derivs_unit, '(a,f10.6)') '# v = ', v_P
    last_v_marked_val = v_val
  end if

  h_P    = 0.25_dp * (h_N + h_S + h_E + h_W)
  dhdu_P = (h_W - h_S + h_N - h_E) * 0.5_dp / du
  dhdv_P = (h_E - h_S + h_N - h_W) * 0.5_dp / dv

  call F(dhduv_P, h_P, dhdu_P, dhdv_P, model_cfg)
  call compute_diagnostics(mass, ricci, h_P, dhdu_P, dhdv_P, dhduv_P, model_cfg)

  ! Write columnar output
  write(fields_unit, '(7e16.8)') u_P, h_P(1), h_P(2), h_P(3) ! u, r, sigma, phi
  write(diag_unit,   '(7e16.8)') u_P, mass, ricci            ! u, mass, Ricci
  write(derivs_unit, '(7e16.8)') u_P, dhdu_P(1), dhdv_P(1)   ! u, drdu, drdv

end subroutine write_output

!====================================================================================
!> Compute model-dependent diagnostics (mass, Ricci)
subroutine compute_diagnostics(mass, ricci, h, dhdu, dhdv, dhduv, model_cfg)
  real(dp), intent(out)               :: mass, ricci
  real(dp), dimension(:), intent(in)  :: h, dhdu, dhdv, dhduv
  type(model_config), intent(in)      :: model_cfg

  real(dp) :: q

  q = model_cfg%q

  mass = 0.5_dp * h(1) * (                         &
         + 1                                       &
         + q*q/(h(1)*h(1))                         &
         + 2*dhdu(1) * dhdv(1) * exp(-2*h(2)) )

  ricci = 2/(h(1)*h(1)) * (                             &
              1 + 2*exp(-2*h(2)) * dhdu(1)*dhdv(1)      &
              )                                         &
          + 4*exp(-2*h(2))*dhduv(2) + 8*dhduv(1)*exp(-2*h(2)) / h(1)

end subroutine compute_diagnostics

!====================================================================================
!> Compute and write constraint violations.
!!
!! This routine needs to be called at v = const slices since it needs 2nd derivatives in u,
!! and it's therefore simpler to implement with access to all u grid points at the given v.
!! The constraint violation along u = const slices is not computed here, since it would require
!! storing the entire v-grid in memory.
!! To avoid complications with AMR, we do only the non-AMR part of the grid,
!! i.e., up to the original Nu.
subroutine write_constraints(h_v, u, v_val, grid_cfg, model_cfg)
  real(dp), intent(in)             :: h_v(:,:)
  real(dp), intent(in)             :: u(:)
  real(dp), intent(in)             :: v_val
  type(grid_config),  intent(in)   :: grid_cfg
  type(model_config), intent(in)   :: model_cfg

  real(dp) :: temp, r, r_uu, r_u, sigma_u, phi_u, du, Guu
  integer  :: j, Nu
  logical  :: v_ok
  real(dp), save :: last_v_marked_val = -1.0e99_dp

  Nu = size(u)
  ! this is the grid spacing of the original grid, before AMR
  du = grid_cfg%du

  if (.not. constraints_open) error stop 'write_output: output files not open'

  ! As above, check output condition first, before doing any work.
  ! Output condition: check if (v - v_min)/output_dv are near integers.
  if (grid_cfg%output_dv > 0.0_dp) then
    temp = abs((v_val - grid_cfg%v_min) / grid_cfg%output_dv)
    v_ok = abs(temp - nint(temp)) < 1.0e-6_dp
  else
    v_ok = .true.
  end if

  if (.not. v_ok) return

  ! For this output, we compute everything at the slice v = v_val
  if (abs(v_val - last_v_marked_val) > 0.5_dp * grid_cfg%output_dv) then
    write(constraints_unit, '(a)')
    write(constraints_unit, '(a,f10.6)') '# v = ', v_val
    last_v_marked_val = v_val
  end if


  do j = 2, Nu - 1
    r       = h_v(j, 1)
    ! first derivatives
    r_u     = (h_v(j+1, 1) - h_v(j-1, 1)) * 0.5_dp / du
    sigma_u = (h_v(j+1, 2) - h_v(j-1, 2)) * 0.5_dp / du
    phi_u   = (h_v(j+1, 3) - h_v(j-1, 3)) * 0.5_dp / du
    ! second derivative
    r_uu    = (h_v(j+1, 1) - 2.0_dp*h_v(j, 1) + h_v(j-1, 1)) / (du*du)

    Guu = r_uu - 2*r_u*sigma_u + r*phi_u*phi_u

    write(constraints_unit, '(7e16.8)') u(j), Guu
  end do

end subroutine write_constraints

end module model_mod
