! Output subroutine for diagnostics

module imprime_mod
  implicit none
  private
  public :: imprime

contains

  !> Output solution diagnostics to file
  !! Writes u, v, field values (h1), and derived quantities (h2)
  subroutine imprime(id, u, v, h1, h2)
    implicit none

    integer, intent(in) :: id
    double precision, intent(in) :: u, v
    double precision, dimension(:), intent(in) :: h1, h2

    write(id,*) (/ u, v, h1, h2 /)

  end subroutine imprime

end module imprime_mod
