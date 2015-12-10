module writemod

  implicit none

contains

  subroutine write(nom, num, coord, valeur)
    
    implicit none

    integer, intent(in) :: num
    real(8), dimension(:,:), intent(in) :: coord
    real(8), dimension(:), intent(in) :: valeur
    character(len=*), intent(in) :: nom
    integer :: ci, di, ui, i, j
    character(len=1) :: c, d, u

    ci = int(num/100)
    di = modulo(int(num/10),10)
    ui = modulo(num,10)

    write(c,'(I1)') ci
    write(d,'(I1)') di
    write(u,'(I1)') ui

    open(unit=22, file=nom//c//d//u//".txt", status="replace")

    write(22,*) '"x", "y", "level"'

    do i = 1, size(coord,1)
       write(22,*) coord(i,1)+0.1,", ", coord(i,2)+0.1,", ", valeur(i)
    end do

    close(22)
    
  end subroutine write

end module writemod
