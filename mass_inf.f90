
! necessario esta interface pois a rotina 'imprime' e a funcao 'polint' aceitam
! arrays com dimensao assumida
module iface
  interface imprime
    subroutine imprime( id, u, v, h1, h2 )
     implicit none
     integer,                         intent(in) :: id
     double precision,                intent(in) :: u, v
     double precision, dimension(:),  intent(in) :: h1, h2
     integer dim1, dim2
     integer k
    end subroutine imprime
  end interface

  interface polint
    function polint( x, xa, ya )
     implicit none
     double precision polint
     double precision, dimension(:), intent(in) :: xa, ya
     double precision, intent(in) :: x
     integer i, m, ns, N
     double precision ho, hp, den, diff, difftemp, w
     double precision, allocatable, dimension(:) :: C, D
    end function polint
  end interface

end module iface
!
!=======================================================================================
!
program mass_inflation
  use iface
  use physics_config_mod
  use simulation_config_mod
  use functions
  use evolve_wrapper, only: step, set_cfg
  implicit none

  ! Physics and simulation configuration
  type(physics_config) :: cfg
  type(simulation_config) :: sim_cfg
  integer             :: neq, D
  double precision    :: lambda, q2, upos, v

  double precision, allocatable, dimension(:,:) :: h_u0, h_v0

  integer,          allocatable, dimension(:)   :: plus, minus
  double precision, allocatable                 :: u(:)
  double precision, allocatable, dimension(:,:) :: h_v0_new

  double precision, allocatable, dimension(:) :: h_S, h_E, h_W, h_N, grad
  double precision, allocatable, dimension(:) :: h_P, dhdu_P, dhdv_P, dhduv_P

  double precision mass, Ricci

  ! para fazer a interpolacao (ordem 5, neste caso)
  double precision, dimension(5) :: interp_x, interp_y

  ! ficheiro para imprimir os valores
  character(len=20), parameter :: filename = 'data.dat'

  double precision tempu, tempv
  double precision :: Du, Dv, u0, v0, uf, vf
  integer          :: i, j, k, jm1, jm2, jm3, jp1, countlast

  ! Initialize physics configuration
  ! Initialize configurations
  call init_physics_config(cfg)
  call init_simulation_config(sim_cfg)
  neq     = cfg%neq
  D       = cfg%D
  lambda  = cfg%lambda
  q2      = cfg%q2

  ! Create local aliases for readability
  Du = sim_cfg%Du
  Dv = sim_cfg%Dv
  u0 = sim_cfg%u0
  v0 = sim_cfg%v0
  uf = sim_cfg%uf
  vf = sim_cfg%vf

  ! Set global config for evolve_wrapper to use
  call set_cfg(cfg)

  open(unit=10, file=filename)
  call print_simulation_header(10, cfg, 0.0d0)

  ! Initialize boundary conditions (returns allocated h_u0, h_v0)
  call init_cond(h_u0, h_v0, sim_cfg, cfg)

  allocate( plus(sim_cfg%big_dim), minus(sim_cfg%big_dim) )
  allocate( u(sim_cfg%big_dim) )
  allocate( h_v0_new(neq, sim_cfg%big_dim) )
  allocate( h_S(neq), h_E(neq), h_W(neq), h_N(neq), grad(neq) )
  allocate( h_P(neq), dhdu_P(neq), dhdv_P(neq), dhduv_P(neq) )

  u     = 0
  plus  = 0
  minus = 0

  upos = sim_cfg%u0
  do i = 1, sim_cfg%Nu
    u(i) = sim_cfg%u0 + (i-1) * sim_cfg%Du
  end do

  ! write(10,'(a1,8a16)') '#', 'u','v','r','phi','sigma', 'mass', 'drdv', 'Ricci'
  write(10,'(a)') '# | u | v | r | phi | sigma | mass | drdv | Ricci'
  ! ! vamos imprimir as condicoes iniciais
  ! para escrever isto aqui temos de meter tb uma maneira de sacar o drdv e o Ricci inicialmente...
  !
  ! do j = 1, Nu
  ! do j = 2, Nu
  !    h_N(:) = h_v0(:,j)
  !    tempu = abs(u(j) - u0) * resu
  !    if( tempu - int(tempu + Du*0.1) < 1.0d-6 )                                        &
  !         call imprime(10, u(j), v0, h_N, (/ mass /) )
  ! end do

  do i = 1, sim_cfg%big_dim
    minus(i) = i-1
    plus(i)  = i+1
  end do

  countlast = sim_cfg%Nu + 1
  v = sim_cfg%v0

  ! inicio integracao. i sera' o passo em 'v', j o passo em 'u'.
  ! em cada passo assumimos que estamos no ponto (u,v).
  do i = 1, sim_cfg%Nv - 1      ! 'i' faz-nos avancar de v para v + Dv.

    ! comecamos por dar a h_N o valor de h_u0, ie o valor de h em (u, v + Dv).
    h_N(:) = h_u0(:,i+1)

    ! para os 'novos' valores de h_v0
    h_v0_new(:,1) = h_N(:)

    ! reinicializamos a posicao em 'u' de cada vez que chegamos ao fim da grelha.
    upos = u0

    ! se quisermos imprimir os dados aqui, temos de corrigir a expressao para a massa...
    !
    ! jm1 = minus(j)
    ! h_S(:) = h_v0(:,jm1)
    ! h_E(:) = h_N(:)
    ! h_W(:) = h_v0(:,j)
    ! h_P    = 0.5*(h_E + h_W)
    ! dhdu_P = (h_W - h_S + h_N - h_E)*0.5d0/Du
    ! dhdv_P = (h_E - h_S + h_N - h_W)*0.5d0/Dv
    ! mass = 0.5d0*h_P(1)**(D-3) * ( 1.d0 - lambda/3.d0 * h_P(1)*h_P(1)              &
    !      + q2/((h_P(1)*h_P(1))**(D-3)) + 2.d0*dhdu_P(1)*dhdv_P(1)/exp(2.0*h_P(3)) )

    j = plus(1)
    v = v + Dv

    tempv = abs(v - v0) * sim_cfg%resv
    if( tempv - int(tempv + Dv*0.1d0) < 1.0d-6 ) then
      write(10,'(a)') ''
!     call imprime(10, u0, v, h_N, mass)          ! imprimimos os valores neste ponto
    end if

    if( mod(i,100) == 0 .or. i == 1 )                                                &
      write(*,'(a,g10.4,a,g10.4)') 'v = ', v, '|   ', vf

    ! 'j' vai-nos fazer avancar de u para u + Du.
    do while ( upos - sim_cfg%uf < 0.d0 )

      jm1 = minus(j)
      ! dizemos que h_S tem o valor de h_v
      h_S(:) = h_v0(:,jm1)  ! no ponto u, ie, o valor de h no ponto (u,v);
      h_E(:) = h_N(:)       ! h_E tem o valor que h_N teve no passo
      h_W(:) = h_v0(:,j)    ! anterior (ie, o valor de h no ponto (u, v + Dv) );
      ! h_W tem o valor de h no ponto (u + Du, v).

      if (sim_cfg%AMR) then

      ! testamos a condicao desejada para decidir se diminuimos ou nao o
      ! passo de integracao.
      grad = abs( (h_W - h_S)/h_W )

        do while ( grad(1) > sim_cfg%gradmax .and. j >= 4 )
           ! enquanto a condicao desejada nao for satisfeita, vamos
           ! adicionando pontos adicionais 'a grelha.

        jm1 = minus(j)
        jm2 = minus(jm1)
        jm3 = minus(jm2)
        jp1 = plus(j)

           do k = 1, neq
              if( u(jp1) < uf ) then
                 interp_x = (/ u(1), u(jm2), u(jm1), u(j), u(jp1) /)
                 interp_y = (/ h_v0(k,1), h_v0(k,jm2), h_v0(k,jm1), h_v0(k,j),    &
                        h_v0(k,jp1) /)
              else
                 interp_x = (/ u(1), u(jm3), u(jm2), u(jm1), u(j) /)
                 interp_y = (/ h_v0(k,1), h_v0(k,jm3), h_v0(k,jm2), h_v0(k,jm1),  &
                        h_v0(k,j) /)
              end if
              h_W(k) = polint( (u(j) + u(jm1))*0.5d0 , interp_x, interp_y)
           end do

        u(countlast)      = ( u(j) + u(jm1) )*0.5d0
        h_v0(:,countlast) = h_W(:)

        jm1               = minus(j)
        minus(countlast)  = jm1
        plus(countlast)   = j
        plus(jm1)         = countlast
        minus(j)          = countlast
        j                 = countlast
        countlast         = countlast + 1

           ! Du = u(j) - u(minus(j))
           ! call step(h_N,h_S,h_E,h_W,Du,Dv)

           grad = abs( (h_W - h_S)/h_W )
        end do
      end if

      Du = u(j) - u(minus(j))     ! Update local Du (may have changed during AMR)

      ! a rotina 'step' vai-nos entao devolver h_N, que e' o valor de h
      ! no ponto (u + Du, v + Dv).
      call step(h_N, h_S, h_E, h_W, Du, Dv, cfg, 4)

      h_P    = 0.5d0*(h_E + h_W)
      dhdu_P = (h_W - h_S + h_N - h_E)*0.5d0/Du
      dhdv_P = (h_E - h_S + h_N - h_W)*0.5d0/Dv

      mass = 0.5d0*h_P(1)**(D-3) * ( 1.d0 - lambda/3.d0 * h_P(1)*h_P(1)              &
            + q2/((h_P(1)*h_P(1))**(D-3)) + 2.d0*dhdu_P(1)*dhdv_P(1)/exp(2.0d0*h_P(3)) )

      call F(dhduv_P, h_P, dhdu_P, dhdv_P, neq, cfg)

      Ricci = (D-3)*(D-2)/( h_P(1)*h_P(1) ) * (                                      &
            1 + 2*exp(-2*h_P(3)) * dhdu_P(1)*dhdv_P(1)                                &
            )                                                                         &
            + 4*exp(-2*h_P(3))*dhduv_P(3) + 4*(D-2)*dhduv_P(1)*exp(-2*h_P(3)) / h_P(1)

      h_v0_new(:,j) = h_N(:)
      upos = upos + Du

      ! output
      tempv = abs(v - v0) * sim_cfg%resv
      tempu = abs(u(j) - u0) * sim_cfg%resu
      if( abs(tempu - int(tempu + Du*0.5d0)) < 1.0d-7   .and.                          &
            abs(tempv - int(tempv + Dv*0.5d0)) < 1.0d-7 )                             &
            call imprime(10, u(j), v, h_N, (/ mass, dhdv_P(1), Ricci /) )

      j = plus(j)

    end do  ! ciclo while

    ! finalmente podemos escrever h_v0_new em h_v0 (e' ineficiente copiar o array todo)
    do k = 1, countlast + 1
      h_v0(:,k) = h_v0_new(:,k)
    end do

  end do

  close(10)

  deallocate( h_u0, h_v0, u, minus, plus, h_v0_new )

end program mass_inflation
!
!=======================================================================================
!
subroutine imprime( id, u, v, h1, h2 )
  implicit none

  integer,                         intent(in) :: id
  double precision,                intent(in) :: u, v
  double precision, dimension(:),  intent(in) :: h1, h2

  ! integer dim1, dim2
  ! integer k

  ! dim1 = size(h1)
  ! dim2 = size(h2)

  ! write(unit=id,fmt='(2g18.10)', advance='no') u, v
  ! do k = 1, dim1
  !    write(id,fmt='(g18.10)', advance='no') h1(k)
  ! end do
  ! do k = 1, dim2
  !    write(id,fmt='(g20.10)', advance='no') h2(k)
  ! end do

  ! write(id,*) ''

  write(id,*) (/ u, v, h1, h2 /)

end subroutine imprime
