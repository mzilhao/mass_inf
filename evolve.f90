!
! rotina para resolver sistemas de equacoes do tipo
!
! \partial_{uv} h_j(u,v) = F_j ( h_k(u,v), \partial_u h_k, \partial_v h_k) 
! j = 1, ..., numero de equacoes
!
! esta rotina aceita: 
!  - tres vectores, hp1, hp2, hp3 -- os valores (conhecidos) das funcoes h
!    nos pontos p1 = (u,v); p2 = (u, v + Dv); p3 = (u + Du, v) ; 
!  - Du e Dv, os passos de integracao em u e v, resp;
!  - uma funcao F de quatro argumentos:
!     - os vectores hP, dhduP e dhdvP, sendo estes as funcoes h e as suas
!       derivadas em u e v, avaliados num ponto generico P;
!     - o primeiro argumento de F devolve um vector dhduvP, com o lado 
!       direito da equacao acima, ie, \partial_{uv} h_j(u,v), no ponto P.
!     - neste algoritmo o ponto usado sera' p5 = (u + Du/2, v + Dv/2).
!  - a rotina 'evolve' devolve (no primeiro argumento) o vector hp4, ie, 
!    as funcoes h no ponto p4 = (u + Du, v + Dv);
!
subroutine evolve( hp4, hp1, hp2, hp3, Du, Dv )
  use functions
  implicit none

  double precision, dimension(NN), intent(out)   :: hp4
  double precision, dimension(NN), intent(in)    :: hp1, hp2, hp3
  double precision, intent(in)                   :: Du, Dv

  integer j
  double precision, dimension(NN) :: hp5, dhdup5, dhdvp5, dhduvp5
  double precision, dimension(NN) :: dhduvp5new, hp4old

  ! fazemos uma primeira aproximacao para as funcoes h no ponto p4 e com 
  ! esse valor calculamos uma primeira aproximacao para h e suas derivadas
  ! no ponto p5:
  hp4 = hp3 + hp2 - hp1
  hp5 = 0.5d0*(hp1 + hp4)
  !hp5 = 0.5*(hp2 + hp3)
  !hp5 = 0.25*(hp1 + hp4 + hp2 + hp3)
  dhdup5 = (hp3 - hp1 + hp4 - hp2)*0.5d0/Du
  dhdvp5 = (hp2 - hp1 + hp4 - hp3)*0.5d0/Dv
  
  ! calculamos agora o 'lado direito' em p5: 
  call F(dhduvp5, hp5, dhdup5, dhdvp5)

  ! recalculamos agora h em p4:
  hp4 = hp3 + hp2 - hp1 + dhduvp5 * Du * Dv

  ! e agora podemos recalcular h e suas derivadas em p5.
  ! e iterar o processo???
  do j = 1,1

     ! para possivel controlo do erro, para ja nao estamos a utilizar...
     hp4old = hp4

     !  hp5 = 0.25*(hp1 + hp4 + hp2 + hp3)
     dhdup5 = (hp3 - hp1 + hp4 - hp2)*0.5d0/Du
     dhdvp5 = (hp2 - hp1 + hp4 - hp3)*0.5d0/Dv

    ! calculamos uma vez mais o lado direito em p5: 
    call F(dhduvp5new, hp5, dhdup5, dhdvp5)

    ! e por fim temos uma melhor aproximacao para h em p4
    hp4 = hp3 + hp2 - hp1 + Du * Dv * dhduvp5new 

 end do


end subroutine evolve
