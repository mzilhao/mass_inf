program mass_inflation
  use physics_config_mod
  use simulation_config_mod
  use functions
  use evolve_wrapper, only: step, set_cfg
  use polint_mod
  use imprime_mod
  implicit none

  ! Physics and simulation configuration
  type(physics_config)    :: cfg
  type(simulation_config) :: sim_cfg
  integer                 :: neq
  double precision        :: upos, v, grad_r

  double precision, allocatable, dimension(:,:) :: h_u0, h_v0
  integer,          allocatable, dimension(:)   :: plus, minus
  double precision, allocatable                 :: u(:)
  double precision, allocatable, dimension(:,:) :: h_v0_new
  double precision, allocatable, dimension(:)   :: h_S, h_E, h_W, h_N
  double precision, allocatable, dimension(:)   :: h_P, dhdu_P, dhdv_P, dhduv_P

  double precision :: mass, drdv, ricci
  double precision :: du, dv, u0, v0, uf, vf
  integer          :: Nv, Nu, Nu_max
  double precision :: tempu, tempv
  double precision, dimension(5) :: interp_x, interp_y

  integer :: i, j, k, jm1, jm2, jm3, jp1, countlast
  character(len=20), parameter :: filename = 'data.dat'

  call init_physics_config(cfg)
  call init_simulation_config(sim_cfg)
  neq = cfg%neq

  ! Create local aliases for readability
  du = sim_cfg%du
  dv = sim_cfg%dv
  u0 = sim_cfg%u0
  v0 = sim_cfg%v0
  uf = sim_cfg%uf
  vf = sim_cfg%vf
  Nu = sim_cfg%Nu
  Nv = sim_cfg%Nv
  Nu_max = sim_cfg%Nu_max

  ! Set global config for evolve_wrapper to use
  ! FIXME: This is a temporary workaround until we refactor evolve_wrapper
  call set_cfg(cfg)

  open(unit=10, file=filename)
  call print_simulation_header(10, cfg, 0.0d0)

  ! Initialize boundary conditions (returns allocated h_u0, h_v0)
  call init_cond(h_u0, h_v0, sim_cfg, cfg)

  allocate(plus(Nu_max), minus(Nu_max))
  allocate(u(Nu_max))
  allocate(h_v0_new(neq, Nu_max))
  allocate(h_S(neq), h_E(neq), h_W(neq), h_N(neq))
  allocate(h_P(neq), dhdu_P(neq), dhdv_P(neq), dhduv_P(neq))

  u = 0.0d0
  plus = 0
  minus = 0

  upos = u0
  do i = 1, Nu
    u(i) = u0 + (i - 1) * du
  end do

  do i = 1, Nu_max
    minus(i) = i - 1
    plus(i)  = i + 1
  end do

  countlast = Nu + 1
  v = sim_cfg%v0


  write(10,'(a)') '# | u | v | r | phi | sigma | mass | drdv | Ricci'

  ! inicio integracao. i sera' o passo em 'v', j o passo em 'u'.
  ! em cada passo assumimos que estamos no ponto (u,v).
  do i = 1, Nv - 1

    ! h_N <- h(u, v + dv)
    h_N(:) = h_u0(:, i + 1)

    ! preparar novos valores de h_v0
    h_v0_new(:, 1) = h_N(:)

    ! reinicializar posicao em u de cada vez que chegamos ao fim da grelha.
    upos = u0

    j = plus(1)
    v = v + dv

    tempv = abs(v - v0) * sim_cfg%resv
    if (tempv - int(tempv + dv*0.1d0) < 1.0d-6) then
      write(10,'(a)') ''
    end if

    if (mod(i, 100) == 0 .or. i == 1) then
      write(*,'(a,g10.4,a,g10.4)') 'v = ', v, '|   ', vf
    end if

    ! avancar de u para u + du
    do while (upos - sim_cfg%uf < 0.d0)

      jm1 = minus(j)
      ! dizemos que h_S tem o valor de h_v no ponto u, ie, o valor de h no ponto (u,v);
      ! h_E tem o valor que h_N teve no passo anterior (ie, o valor de h no ponto (u, v + dv) );
      ! h_W tem o valor de h no ponto (u + du, v).
      h_S(:) = h_v0(:, jm1)  ! h(u, v)
      h_E(:) = h_N(:)        ! h(u, v + dv)
      h_W(:) = h_v0(:, j)    ! h(u + du, v)

      if (sim_cfg%AMR) then
        grad_r = abs((h_W(1) - h_S(1))/(h_W(1) + h_S(1) + 1.0d-16))

        ! testamos a condicao desejada para decidir se diminuimos ou nao o
        ! passo de integracao. enquanto a condicao desejada nao for satisfeita, vamos
        ! adicionando pontos adicionais 'a grelha.
        do while (grad_r > sim_cfg%gradmax .and. j >= 4)
          jm1 = minus(j)
          jm2 = minus(jm1)
          jm3 = minus(jm2)
          jp1 = plus(j)

          do k = 1, neq
            if (u(jp1) < uf) then
              interp_x = (/ u(1), u(jm2), u(jm1), u(j), u(jp1) /)
              interp_y = (/ h_v0(k,1), h_v0(k,jm2), h_v0(k,jm1), h_v0(k,j), &
                            h_v0(k,jp1) /)
            else
              interp_x = (/ u(1), u(jm3), u(jm2), u(jm1), u(j) /)
              interp_y = (/ h_v0(k,1), h_v0(k,jm3), h_v0(k,jm2), h_v0(k,jm1), &
                            h_v0(k,j) /)
            end if
            h_W(k) = polint((u(j) + u(jm1))*0.5d0, interp_x, interp_y)
          end do

          u(countlast)       = (u(j) + u(jm1))*0.5d0
          h_v0(:, countlast) = h_W(:)

          jm1              = minus(j)
          minus(countlast) = jm1
          plus(countlast)  = j
          plus(jm1)        = countlast
          minus(j)         = countlast
          j                = countlast
          countlast        = countlast + 1

          grad_r = abs((h_W(1) - h_S(1))/(h_W(1) + h_S(1) + 1.0d-16))
        end do
      end if

      du = u(j) - u(minus(j))  ! local du may change during AMR

      ! step returns h_N = h(u + du, v + dv)
      call step(h_N, h_S, h_E, h_W, du, dv, cfg, 4)

      h_P    = 0.5d0*(h_E + h_W)
      dhdu_P = (h_W - h_S + h_N - h_E)*0.5d0 / du
      dhdv_P = (h_E - h_S + h_N - h_W)*0.5d0 / dv

      call F(dhduv_P, h_P, dhdu_P, dhdv_P, neq, cfg)
      call compute_diagnostics(h_P, dhdu_P, dhdv_P, dhduv_P, cfg, mass, drdv, ricci)

      h_v0_new(:, j) = h_N(:)
      upos = upos + du

      tempv = abs(v - v0) * sim_cfg%resv
      tempu = abs(u(j) - u0) * sim_cfg%resu
      if (abs(tempu - int(tempu + du*0.5d0)) < 1.0d-7 .and. &
          abs(tempv - int(tempv + dv*0.5d0)) < 1.0d-7) then
        call imprime(10, u(j), v, h_N, (/ mass, drdv, ricci /))
      end if

      j = plus(j)

    end do  ! ciclo while

    ! finalmente podemos escrever h_v0_new em h_v0 (e' ineficiente copiar o array todo)
    do k = 1, countlast + 1
      h_v0(:, k) = h_v0_new(:, k)
    end do

  end do

  close(10)

  deallocate(h_u0, h_v0, u, minus, plus, h_v0_new)

end program mass_inflation
