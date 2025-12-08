
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
  use functions
  use evolve_wrapper, only: step, set_cfg
  implicit none

  ! Physics configuration
  type(physics_config) :: cfg
  integer             :: neq, D
  double precision    :: lambda, q2

  ! definidos pelas condicoes iniciais
  double precision :: Du, Dv, u0, v0, uf, vf, gradmax, gradmin, upos, v
  logical          :: AMR
  integer          :: Nu, Nv, resu, resv

  integer,          allocatable, dimension(:)   :: plus, minus
  double precision, allocatable                 :: u(:)
  double precision, allocatable, dimension(:,:) :: h_v0_new

  double precision, allocatable, dimension(:) :: hp1, hp2, hp3, hp4, grad
  double precision, allocatable, dimension(:) :: hp5, dhdup5, dhdvp5, dhduvp5

  double precision mass, Ricci

  ! para fazer a interpolacao (ordem 5, neste caso)
  double precision, dimension(5) :: interp_x, interp_y

  ! ficheiro para imprimir os valores
  character(len=20), parameter :: filename = 'data.dat'

  double precision tempu, tempv
  integer          i, j, k, jm1, jm2, jm3, jp1, countlast

  ! Initialize physics configuration
  call init_physics_config(cfg)
  neq     = cfg%neq
  D       = cfg%D
  lambda  = cfg%lambda
  q2      = cfg%q2

  ! Set global config for evolve_wrapper to use
  call set_cfg(cfg)

  open(unit=10, file=filename)

  call init_cond( Du, Dv, u0, v0, mass, uf, vf, Nu, Nv, gradmax, gradmin, AMR, resu, resv, 10, cfg)

  allocate( plus(big_dim), minus(big_dim) )
  allocate( u(big_dim) )
  allocate( h_v0_new(neq, big_dim) )
  allocate( hp1(neq), hp2(neq), hp3(neq), hp4(neq), grad(neq) )
  allocate( hp5(neq), dhdup5(neq), dhdvp5(neq), dhduvp5(neq) )

  u     = 0
  plus  = 0
  minus = 0

  upos = u0
  do i = 1, Nu
    u(i) = u0 + (i-1)*Du
  end do

  ! write(10,'(a1,8a16)') '#', 'u','v','r','phi','sigma', 'mass', 'drdv', 'Ricci'
  write(10,'(a)') '# | u | v | r | phi | sigma | mass | drdv | Ricci'
  ! ! vamos imprimir as condicoes iniciais
  ! para escrever isto aqui temos de meter tb uma maneira de sacar o drdv e o Ricci inicialmente...
  !
  ! do j = 1, Nu
  ! do j = 2, Nu
  !    hp4(:) = h_v0(:,j)
  !    tempu = abs(u(j) - u0) * resu
  !    if( tempu - int(tempu + Du*0.1) < 1.0d-6 )                                        &
  !         call imprime(10, u(j), v0, hp4, (/ mass /) )
  ! end do

  do i = 1, big_dim
    minus(i) = i-1
    plus(i)  = i+1
  end do

  countlast = Nu + 1
  v = v0

  ! inicio integracao. i sera' o passo em 'v', j o passo em 'u'.
  ! em cada passo assumimos que estamos no ponto (u,v).
  do i = 1, Nv - 1      ! 'i' faz-nos avancar de v para v + Dv.

    ! comecamos por dar a hp4 o valor de h_u0, ie o valor de h em (u, v + Dv).
    hp4(:) = h_u0(:,i+1)

    ! para os 'novos' valores de h_v0
    h_v0_new(:,1) = hp4(:)

    ! reinicializamos a posicao em 'u' de cada vez que chegamos ao fim da grelha.
    upos = u0

    ! se quisermos imprimir os dados aqui, temos de corrigir a expressao para a massa...
    !
    ! jm1 = minus(j)
    ! hp1(:) = h_v0(:,jm1)
    ! hp2(:) = hp4(:)
    ! hp3(:) = h_v0(:,j)
    ! hp5    = 0.5*(hp2 + hp3)
    ! dhdup5 = (hp3 - hp1 + hp4 - hp2)*0.5d0/Du
    ! dhdvp5 = (hp2 - hp1 + hp4 - hp3)*0.5d0/Dv
    ! mass = 0.5d0*hp5(1)**(D-3) * ( 1.d0 - lambda/3.d0 * hp5(1)*hp5(1)              &
    !      + q2/((hp5(1)*hp5(1))**(D-3)) + 2.d0*dhdup5(1)*dhdvp5(1)/exp(2.0*hp5(3)) )

    j = plus(1)
    v = v + Dv

    tempv = abs(v - v0) * resv
    if( tempv - int(tempv + Dv*0.1) < 1.0d-6 ) then
      write(10,'(a)') ''
!     call imprime(10, u0, v, hp4, mass)          ! imprimimos os valores neste ponto
    end if

    if( mod(i,100) == 0 .or. i == 1 )                                                &
      write(*,'(a,g10.4,a,g10.4)') 'v = ', v, '|   ', vf

    ! 'j' vai-nos fazer avancar de u para u + Du.
    do while ( upos - uf < 0.d0 )

      jm1 = minus(j)
      ! dizemos que hp1 tem o valor de h_v
      hp1(:) = h_v0(:,jm1)  ! no ponto u, ie, o valor de h no ponto (u,v);
      hp2(:) = hp4(:)       ! hp2 tem o valor que hp4 teve no passo
      hp3(:) = h_v0(:,j)    ! anterior (ie, o valor de h no ponto (u, v + Dv) );
      ! hp3 tem o valor de h no ponto (u + Du, v).

      if (AMR) then

      ! testamos a condicao desejada para decidir se diminuimos ou nao o
      ! passo de integracao.
      grad = abs( (hp3 - hp1)/hp3 )

        do while ( grad(1) > gradmax .and. j >= 4 )
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
              hp3(k) = polint( (u(j) + u(jm1))*0.5d0 , interp_x, interp_y)
           end do

        u(countlast)      = ( u(j) + u(jm1) )*0.5d0
        h_v0(:,countlast) = hp3(:)

        jm1               = minus(j)
        minus(countlast)  = jm1
        plus(countlast)   = j
        plus(jm1)         = countlast
        minus(j)          = countlast
        j                 = countlast
        countlast         = countlast + 1

           ! Du = u(j) - u(minus(j))
           ! call step(hp4,hp1,hp2,hp3,Du,Dv)

           grad = abs( (hp3 - hp1)/hp3 )
        end do
      end if

      Du = u(j) - u(minus(j))              ! usando AMR, o passo de integracao
                                                            ! sera' variavel.

      ! a rotina 'step' vai-nos entao devolver hp4, que e' o valor de h
      ! no ponto (u + Du, v + Dv).
      call step(hp4, hp1, hp2, hp3, Du, Dv, cfg, 4)

      hp5    = 0.5*(hp2 + hp3)
      dhdup5 = (hp3 - hp1 + hp4 - hp2)*0.5d0/Du
      dhdvp5 = (hp2 - hp1 + hp4 - hp3)*0.5d0/Dv

      mass = 0.5d0*hp5(1)**(D-3) * ( 1.d0 - lambda/3.d0 * hp5(1)*hp5(1)              &
            + q2/((hp5(1)*hp5(1))**(D-3)) + 2.d0*dhdup5(1)*dhdvp5(1)/exp(2.0*hp5(3)) )

      call F(dhduvp5, hp5, dhdup5, dhdvp5, neq, cfg)

      Ricci = (D-3)*(D-2)/( hp5(1)*hp5(1) ) * (                                      &
            1 + 2*exp(-2*hp5(3)) * dhdup5(1)*dhdvp5(1)                                &
            )                                                                         &
            + 4*exp(-2*hp5(3))*dhduvp5(3) + 4*(D-2)*dhduvp5(1)*exp(-2*hp5(3)) / hp5(1)

      h_v0_new(:,j) = hp4(:)
      upos = upos + Du

      ! output
      tempv = abs(v - v0) * resv
      tempu = abs(u(j) - u0) * resu
      if( abs(tempu - int(tempu + Du*0.5)) < 1.0d-7   .and.                          &
            abs(tempv - int(tempv + Dv*0.5)) < 1.0d-7 )                               &
            call imprime(10, u(j), v, hp4, (/ mass, dhdvp5(1), Ricci /) )

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
