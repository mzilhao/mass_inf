module functions
  implicit none

  ! constantes globais
  integer,          parameter :: NN     = 3                    ! numero equacoes
  ! integer,          parameter :: D      = 5                    ! dimensao espaco-tempo
  integer,          parameter :: D      = 4                    ! dimensao espaco-tempo
  double precision, parameter :: lambda = 0.d0                 ! constante cosmologica
  double precision, parameter :: q      = 0.95d0               ! carga electrica
  double precision, parameter :: q2     = q*q
  double precision, parameter :: qq2    = 0.5d0*q2*(D-3)*(D-2)
  double precision, parameter :: qq     = sqrt(qq2)
  double precision, parameter :: Pi     = 3.1415926535897932384626433d0

  integer                        big_dim

  ! variaveis globais
  double precision, allocatable, dimension(:,:) :: h_u0, h_v0

  !====================================================================================
contains
  ! h(1) -> r
  ! h(2) -> phi
  ! h(3) -> sigma
  !
  ! funcao que nos da' o lado direito da equacao que queremos resolver:
  subroutine F( dhduv, h, dhdu, dhdv )
    implicit none

    double precision, dimension(NN), intent(out) :: dhduv
    double precision, dimension(NN), intent(in)  :: h, dhdu, dhdv

    dhduv(1) =  qq2/(D-2) *  exp(2 * h(3)) / ( h(1) ** (2*D - 5) )              &
         + (D-1)/6.d0 * lambda * h(1) * exp(2 * h(3))                           &
         - (D-3)/2.d0 * exp(2 * h(3)) / h(1)                                    &
         - (D-3) * dhdu(1) * dhdv(1) / h(1)

    dhduv(2) = - (D-2)/(2.d0*h(1)) * ( dhdv(1) * dhdu(2) + dhdu(1) * dhdv(2) ) 

    dhduv(3) = - dhdu(2) * dhdv(2)                                              &
         - (3*D-8)/(2.d0*D-4) * qq2 * exp(2 * h(3)) / ( h(1) ** ( 2*(D - 2)) )  &
         - (D-4)*(D-1) * lambda/12.d0 * exp(2 * h(3))                           &
         + (D-3)*(D-2)/4.d0 * exp(2 * h(3)) / (h(1)*h(1))                       &
         + (D-3)*(D-2)/2.d0 * dhdu(1) * dhdv(1) / (h(1)*h(1)) 

  end subroutine F
  !
  !====================================================================================
  !
  subroutine init_cond( Du, Dv, u0, v0, m0, uf, vf, Nu, Nv, gradmax, gradmin, AMR, resu, &
       resv, id )

    implicit none

    double precision, intent(out) :: Du, Dv, u0, v0, uf, vf, gradmax, gradmin, m0
    integer,          intent(out) :: Nu, Nv, resu, resv
    logical,          intent(out) :: AMR
    integer,          intent(in)  :: id

    double precision :: r00, su0v0, ru0, v1
    double precision :: A, Delta
    logical          :: scalarfield

    integer i
    double precision u, v

    AMR = .true.

    ! limites de integracao
    u0 = 0.d0
    v0 = 5.d0
    ! uf = 21.d0
    uf = 30.d0

    ! vf = 20.d0
    vf = 15.d0

    resu  = 20                   ! resolucao que queremos em u ao imprimir os resultados
    resv  = 20                   ! resolucao que queremos em v ao imprimir os resultados

    gradmax = 0.0001d0           ! para controlo do AMR. maxima variacao
                                 ! permitida nas variaveis

    gradmin = 0.1d0              ! por enquanto nao estamos a usar esta variavel...
  
    ! estimativa para o tamanho dos arrays
    big_dim = int( 2*(uf - u0)/gradmax )

    allocate( h_u0(NN, big_dim), h_v0(NN, big_dim) )
    h_u0  = 0
    h_v0  = 0

    !Du = (uf - u0)/(Nu - 1)
    !Dv = (vf - v0)/(Nv - 1)

    Du = 0.01d0
    Dv = 0.0005d0

    ! numero de pontos da grelha (inicialmente)
    ! (o + 0.001 e' para garantir que nao se perde informacao na conversao para int...)
    Nu = int( (uf - u0)/Du + 1 + 0.001 )          ! grelha em u
    Nv = int( (vf - v0)/Dv + 1 + 0.001 )          ! grelha em v


    ! condicoes iniciais
    r00    = v0
    su0v0  = - 0.5d0*log(2.d0)
    m0     = 1.
    ru0    = 0.25d0*(2.d0/( r00**(D-3) ) * ( m0 - q2/(2.d0* r00**(D-3)) )           & 
         - 1.d0 + lambda * r00*r00/3.d0)

    ! A      = 0.05d0
    ! A      = 0.1d0
    A      = 0.0d0
    Delta  = 1.0d0

    scalarfield = .true.

    write(id,'(a,i2)')    '# D         = ', D
    write(id,'(a,g10.4)') '# lambda    = ', lambda
    write(id,'(a,g10.4)') '# q         = ', q
    write(id,'(a,g10.4)') '# A         = ', A
    write(id,'(a,g10.4)') '# sigma_0   = ', su0v0
    write(id,'(a,g10.4)') '# m0        = ', m0
    write(id,'(a)') '#'

    v1 = v0 + Delta

    ! em u = u0:
    do i = 1, Nv
       v = v0 + (i-1)*Dv
       ! h_u0[v]:
       h_u0(1,i) = v               ! r(u0,v)

       if(scalarfield) then
          if ( v <= v0 + Delta ) then

             ! phi(u0,v):
             h_u0(2,i) = A/(4*Pi) * ( 2*Pi*(v-v0)                           &
                  - Delta*sin(2.d0*Pi*(v-v0)/Delta) )

             ! sigma(u0,v):
             ! h_u0[3][i] = su0v0 
             !   + A*A / ( (D-2)*128.*Pi*Pi) * 
             !   ( 15.*Delta*Delta + 24*Pi*Pi*v*v - 24.*Pi*Pi*v0*v0 
             !     - 16.*Delta*Delta * cos(2.*Pi*(Delta + v - v0)/Delta) 
             !     - Delta*Delta * cos(Pi*(-3.*Delta - 4.*v + 4.*v0)/Delta) 
             !     - 32.*Delta*Pi*v * sin(2.*Pi*(Delta + v - v0)/Delta) 
             !     + 4.*Delta*Pi*v * sin(Pi*(-3.*Delta - 4.*v + 4.*v0)/Delta)
             !     ) ;
       
             ! pedro:
             h_u0(3,i) = su0v0 +                                                  &
                  2.0d0/(D-2.0) * A*A/(256*Pi*Pi) *                               &
                  ( 15*v0*v0 - 24*Pi*Pi*v0*v0 - 30*v0*v1                          &
                  + 15*v1*v1 + 24*Pi*Pi*v*v                                       &
                  - 16*Delta*Delta*cos(2*Pi*(v0-v)/Delta)                         &
                  + Delta*Delta*cos(4*Pi*(v0-v)/Delta)                            &
                  - 32*Pi*v0*v*sin(2*Pi*(v0-v)/Delta)                             &
                  + 32*Pi*v1*v*sin(2*Pi*(v0-v)/Delta)                             &
                  + 4*Pi*v0*v*sin(4*Pi*(v0-v)/Delta)                              &
                  - 4*Pi*v1*v*sin(4*Pi*(v0-v)/Delta)                              &
                  )
             
          else 
             ! phi(u0,v):
             h_u0(2,i) = 0.5d0 *A * Delta 

             ! sigma(u0,v):
             ! h_u0(3,i) = su0v0 + 3./(16.*(D-2)) * A*A *(v*v - v0*v0)

             ! pedro:
             h_u0(3,i) = su0v0 + 2.0d0/(D-2.0d0) * 3*A*A*Delta*(v0+v1)/32.0d0
 
          end if
          
       else ! sem campo escalar
          ! phi(u0,v):
          h_u0(2,i) = 0.d0
          ! sigma(u0,v):
          h_u0(3,i) = su0v0
          
       end if
  
    end do

    ! em v = v0:  
    do i = 1, Nu
       u = u0 + (i-1)*Du
       ! h_v0[u]:
       h_v0(1,i) = r00 + u * ru0    ! r(u,v0)
       h_v0(2,i) = 0.d0             ! phi(u,v0)
       h_v0(3,i) = su0v0            ! sigma(u,v0)
       
    end do

  end subroutine init_cond
  !
  !======================================================================================
  !
end module functions
