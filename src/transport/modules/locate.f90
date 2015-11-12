module locatemod

  implicit none

contains

  subroutine locate(coord, dx, indices)

    implicit none

    real(8), dimension(2), intent(in) :: coord
    real(8), intent(in) :: dx
    integer, dimension(4,2), intent(out) :: indices
    integer :: a, b

    a = ceiling(coord(1)/dx)
    b = ceiling(coord(2)/dx)

    indices(:,1) = (/a, a+1, a+1, a/) 
    indices(:,2) = (/b+1, b+1, b, b/)

  end subroutine locate

end module locatemod
