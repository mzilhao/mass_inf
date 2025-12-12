module progress_utils
  implicit none
  private
  public :: print_banner

contains


  !> Print ASCII banner to stdout at startup
  subroutine print_banner()
    write(*,'(a)') ''
    write(*,'(a)') ''
    write(*,'(a)') '      x\                /x'
    write(*,'(a)') '      x \              / x'
    write(*,'(a)') '      x  \            /  x'
    write(*,'(a)') '      x   \          /   x'
    write(*,'(a)') '      x    \        /    x'
    write(*,'(a)') '      x     \      /     x'
    write(*,'(a)') '      x      \    /      x'
    write(*,'(a)') '      x       \  /       x'
    write(*,'(a)') '      x        \/        x'
    write(*,'(a)') '      x        /\        x'
    write(*,'(a)') '      x       /  \       x'
    write(*,'(a)') '      x      /    \      x'
    write(*,'(a)') '      x     /  /\_/\     x'
    write(*,'(a)') '      x    /  ( O_O )    x'
    write(*,'(a)') '      x   /    > ^ < \   x'
    write(*,'(a)') '      x  /            \  x'
    write(*,'(a)') '      x /              \ x'
    write(*,'(a)') '      x/                \x'
    write(*,'(a)') '       \                /\'
    write(*,'(a)') '        \              /  \'
    write(*,'(a)') '         \            /    \'
    write(*,'(a)') '          \          /      \'
    write(*,'(a)') '           \        /        \'
    write(*,'(a)') '            \      /          \'
    write(*,'(a)') '             \    /            \'
    write(*,'(a)') '              \  /              \'
    write(*,'(a)') '               \/                \'
    write(*,'(a)') '                \                /'
    write(*,'(a)') '                 \              /'
    write(*,'(a)') '                  \            /'
    write(*,'(a)') '                   \          /'
    write(*,'(a)') '                    \        /'
    write(*,'(a)') '                     \      /'
    write(*,'(a)') '                      \    /'
    write(*,'(a)') '                       \  /'
    write(*,'(a)') '                        \/'
    write(*,'(a)') ''
    write(*,'(a)') ''
  end subroutine print_banner

end module progress_utils
