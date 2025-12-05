
!  rotina 'polint', copiada do Numerical Recipes em C++, segunda edicao.  esta
!  rotina efectua uma interpolacao ou extrapolacao polinomial.  ela aceita um
!  vector 'xa', um vector 'ya' de pontos conhecidos e o ponto 'x' que se
!  pretende interpolar (ou extrapolar), e devolve um double 'y', o valor da
!  interpolacao polinomial no ponto 'x'.

function polint( x, xa, ya )
  implicit none

  double precision polint
  double precision, dimension(:), intent(in) :: xa, ya
  double precision, intent(in) :: x

  integer i, m, ns, N
  double precision ho, hp, den, diff, difftemp, w
  double precision, allocatable, dimension(:) :: C, D

  ns = 1
  N  = size(ya)

  allocate( C(N), D(N) )

  diff = abs(x - xa(1))
  do i = 1, N                   ! achamos o ponto dado mais proximo de 'x'

     difftemp = abs(x - xa(i))
     if( difftemp < diff ) then
        ns = i                  ! e guardamos a posicao em 'ns'.
        diff = difftemp 
     end if

     C(i) = ya(i)               ! inicializamos os vectores C's e D's
     D(i) = ya(i)

  end do
 
  polint = ya(ns)               ! primeira aproximacao para polint: dizemos que e'
                                ! igual ao valor do ponto mais proximo de 'x',
                                ! obtido acima
  ns = ns - 1

  do m = 1, N - 1               ! 'm' percorre as colunas do quadro (ver NR)
     do i = 1, N - m            ! para cada coluna calculamos todos os C's
        ho = xa(i) - x          ! e D's.
        hp = xa(i+m) - x
        w = C(i+1) - D(i)

        den = ho - hp
        ! este erro so' acontecera' se os dois xa's forem identicos (dentro 
        ! da aproximacao numerica)
        if( abs(den) < 1.d-12 ) then
           write(*,*) "erro na rotina polint!"
           stop
        end if

        den = w/den
        D(i) = hp * den          ! actualizamos os C's e D's. 
        C(i) = ho * den
     end do

     if( 2 * ns < N - m ) then
        polint = polint + C(ns+1)
     else
        polint  = polint + D(ns)
        ns = ns - 1
     end if

  end do

  return

end function polint
